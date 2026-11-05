"""Lemmus：read-fate 证据与磁盘守恒索引。"""
from __future__ import annotations
from .io import project_temporary_parent, atomic_output

import argparse
import csv
import gzip
import json
import sqlite3
import tempfile
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from itertools import groupby, islice, zip_longest
from pathlib import Path
from .io import (sha256_file)


class ReadFateError(ValueError):
    """The read-fate contract was violated."""


@dataclass(frozen=True, order=True)
class ReadKey:
    sample_id: str
    modality: str
    read_id: str


@dataclass(frozen=True)
class FateRow:
    key: ReadKey
    fate: str
    source: str
    detail: str = ''

    def as_dict(self) -> dict[str, str]:
        return {'sample_id': self.key.sample_id, 'modality': self.key.modality, 'read_id': self.key.read_id, 'fate': self.fate, 'source': self.source, 'detail': self.detail}


_INSERT_BATCH_SIZE = 8192


_INSERT_SQL = 'INSERT INTO fates VALUES (?, ?, ?, ?, ?, ?)'


def _validate_row(sample_id: str, modality: str, fate: str) -> None:
    _validate_sample_and_modality(sample_id, modality)
    _validate_fate(modality, fate)


class SqliteFateStore:
    """Disk-backed fate index so full ledgers do not scale with Python RAM."""

    def __init__(self, path: Path):
        self.path = path
        self.connection = sqlite3.connect(path)
        self.connection.execute('PRAGMA journal_mode=MEMORY')
        self.connection.execute('PRAGMA synchronous=OFF')
        self.connection.execute('PRAGMA temp_store=FILE')
        self.connection.execute('\n            CREATE TABLE fates (\n                sample_id TEXT NOT NULL,\n                modality TEXT NOT NULL,\n                read_id TEXT NOT NULL,\n                fate TEXT NOT NULL,\n                source TEXT NOT NULL,\n                detail TEXT NOT NULL,\n                PRIMARY KEY (sample_id, modality, read_id)\n            ) WITHOUT ROWID\n            ')

    def __enter__(self) -> SqliteFateStore:
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        if exc_type is None:
            self.connection.commit()
        else:
            self.connection.rollback()
        self.connection.close()

    def insert(self, row: FateRow) -> None:
        _validate_row(row.key.sample_id, row.key.modality, row.fate)
        try:
            self.connection.execute(_INSERT_SQL, self._row_values(row))
        except sqlite3.IntegrityError as exc:
            previous = self.get(row.key)
            if previous is not None and previous.fate != row.fate:
                raise ReadFateError(f'Conflicting fates for {row.key}: {previous.fate} vs {row.fate}') from exc
            raise ReadFateError(f'Duplicate fate evidence for {row.key}: {row.fate}') from exc

    @staticmethod
    def _row_values(row: FateRow) -> tuple[str, str, str, str, str, str]:
        return (row.key.sample_id, row.key.modality, row.key.read_id, row.fate, row.source, row.detail)

    @classmethod
    def _validated_row_values(cls, rows: Sequence[FateRow]) -> Iterator[tuple[str, str, str, str, str, str]]:
        """Yield one validated parameter tuple at a time in source order."""
        for row in rows:
            _validate_row(row.key.sample_id, row.key.modality, row.fate)
            yield cls._row_values(row)

    def _insert_batch(self, rows: Sequence[FateRow]) -> None:
        """Insert one batch while preserving the single-row failure contract."""
        savepoint = 'rna_darlin_fate_insert_batch'
        if not self.connection.in_transaction:
            self.connection.execute('BEGIN')
        self.connection.execute(f'SAVEPOINT {savepoint}')
        try:
            self.connection.executemany(_INSERT_SQL, self._validated_row_values(rows))
        except sqlite3.IntegrityError as exc:
            self.connection.execute(f'ROLLBACK TO SAVEPOINT {savepoint}')
            self.connection.execute(f'RELEASE SAVEPOINT {savepoint}')
            for row in rows:
                self.insert(row)
            raise RuntimeError('Batch insert failed but single-row replay did not reproduce it') from exc
        except BaseException:
            self.connection.execute(f'ROLLBACK TO SAVEPOINT {savepoint}')
            self.connection.execute(f'RELEASE SAVEPOINT {savepoint}')
            raise
        else:
            self.connection.execute(f'RELEASE SAVEPOINT {savepoint}')

    def insert_many(self, rows: Iterable[FateRow]) -> int:
        """Insert a stream in fixed batches and return its observed row count."""
        iterator = iter(rows)
        inserted = 0
        while (batch := list(islice(iterator, _INSERT_BATCH_SIZE))):
            self._insert_batch(batch)
            inserted += len(batch)
        return inserted

    def get(self, key: ReadKey) -> FateRow | None:
        record = self.connection.execute('\n            SELECT fate, source, detail FROM fates\n            WHERE sample_id = ? AND modality = ? AND read_id = ?\n            ', (key.sample_id, key.modality, key.read_id)).fetchone()
        if record is None:
            return None
        return FateRow(key, record[0], record[1], record[2])

    def replace(self, key: ReadKey, expected_fate: str, replacement: FateRow) -> None:
        previous = self.get(key)
        if previous is None:
            raise ReadFateError(f'Fate evidence is absent from input inventory: {key}')
        if previous.fate != expected_fate:
            if previous.fate == replacement.fate:
                raise ReadFateError(f'Duplicate fate evidence for {key}: {replacement.fate}')
            raise ReadFateError(f'Conflicting fates for {key}: {previous.fate} vs {replacement.fate}')
        _validate_row(key.sample_id, key.modality, replacement.fate)
        self.connection.execute('\n            UPDATE fates SET fate = ?, source = ?, detail = ?\n            WHERE sample_id = ? AND modality = ? AND read_id = ?\n            ', (replacement.fate, replacement.source, replacement.detail, key.sample_id, key.modality, key.read_id))

    def rows(self) -> Iterator[FateRow]:
        cursor = self.connection.execute('\n            SELECT sample_id, modality, read_id, fate, source, detail\n            FROM fates ORDER BY sample_id, modality, read_id\n            ')
        for sample_id, modality, read_id, fate, source, detail in cursor:
            yield FateRow(ReadKey(sample_id, modality, read_id), fate, source, detail)

    def count_by(self, column: str) -> dict[str, int]:
        if column not in {'sample_id', 'modality', 'fate'}:
            raise ReadFateError(f'Unsupported read-fate count column: {column}')
        return {str(key): int(count) for key, count in self.connection.execute(f'SELECT {column}, COUNT(*) FROM fates GROUP BY {column} ORDER BY {column}')}

    def total(self) -> int:
        return int(self.connection.execute('SELECT COUNT(*) FROM fates').fetchone()[0])


SCHEMA_VERSION = 1


TSV_COLUMNS = ('sample_id', 'modality', 'read_id', 'fate', 'source', 'detail')


DARLIN_SIDECAR_COLUMNS = ('read_id', 'template', 'fate', 'cell_barcode', 'clone_id', 'detail')


FATES_BY_MODALITY = {'RNA': {'RNA_ASSIGNED', 'RNA_UNASSIGNED'}, 'DARLIN': {'DARLIN_ASSIGNED', 'DARLIN_UNASSIGNED', 'DARLIN_FAILED_QC'}}


DARLIN_SIDECAR_FATES = frozenset(FATES_BY_MODALITY['DARLIN'])


def normalize_read_id(value: str | bytes) -> str:
    """规范化 FASTQ/SAM identifier，不改变其生物学名称。"""
    if isinstance(value, bytes):
        value = value.decode('utf-8', errors='strict')
    normalized = value.strip()
    if normalized.startswith('@'):
        normalized = normalized[1:]
    normalized = normalized.split(None, 1)[0]
    if normalized.endswith(('/1', '/2')):
        normalized = normalized[:-2]
    if not normalized:
        raise ReadFateError('Read identifier is empty after normalization')
    return normalized


def _validate_sample_and_modality(sample_id: str, modality: str) -> None:
    if not sample_id.strip():
        raise ReadFateError('sample_id cannot be empty')
    if modality not in FATES_BY_MODALITY:
        raise ReadFateError(f'Unsupported modality {modality!r}; expected one of {sorted(FATES_BY_MODALITY)}')


def _validate_fate(modality: str, fate: str) -> None:
    if fate not in FATES_BY_MODALITY[modality]:
        raise ReadFateError(f'Fate {fate!r} is not valid for modality {modality}; expected one of {sorted(FATES_BY_MODALITY[modality])}')


def _open_binary(path: Path):
    return gzip.open(path, 'rb') if path.suffix == '.gz' else path.open('rb')


def _open_text(path: Path):
    if path.suffix == '.gz':
        return gzip.open(path, 'rt', newline='', encoding='utf-8')
    return path.open('r', newline='', encoding='utf-8')


def _fastq_records(path: Path) -> Iterator[tuple[str, bytes]]:
    with _open_binary(path) as handle:
        record_number = 0
        while True:
            header = handle.readline()
            if not header:
                return
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            record_number += 1
            if not sequence or not plus or (not quality):
                raise ReadFateError(f'Truncated FASTQ record {record_number} in {path}')
            if not header.startswith(b'@') or not plus.startswith(b'+'):
                raise ReadFateError(f'Malformed FASTQ record {record_number} in {path}')
            if len(sequence.rstrip(b'\r\n')) != len(quality.rstrip(b'\r\n')):
                raise ReadFateError(f'Sequence/quality length mismatch at record {record_number} in {path}')
            yield (normalize_read_id(header), sequence.rstrip(b'\r\n'))


def paired_fastq_records(r1: Path, r2: Path) -> Iterator[tuple[str, bytes, bytes]]:
    """产出配对序列前，严格证明 mate 数量、顺序和 identifier 一致。"""
    missing = object()
    for record_number, (record_r1, record_r2) in enumerate(zip_longest(_fastq_records(r1), _fastq_records(r2), fillvalue=missing), start=1):
        if record_r1 is missing or record_r2 is missing:
            raise ReadFateError(f'R1/R2 FASTQ record counts differ for {r1} and {r2}')
        r1_id, r1_sequence = record_r1
        r2_id, r2_sequence = record_r2
        if r1_id != r2_id:
            raise ReadFateError(f'R1/R2 read-ID mismatch at record {record_number}: {r1_id!r} != {r2_id!r}')
        yield (r1_id, r1_sequence, r2_sequence)


def paired_fastq_ids(r1: Path, r2: Path) -> Iterator[str]:
    """产出 read ID 前，严格验证配对的完整性。"""
    for read_id, _, _ in paired_fastq_records(r1, r2):
        yield read_id


def inventory_from_fastq(sample_id: str, modality: str, r1: Path, r2: Path) -> set[ReadKey]:
    _validate_sample_and_modality(sample_id, modality)
    inventory: set[ReadKey] = set()
    for read_id in paired_fastq_ids(r1, r2):
        key = ReadKey(sample_id, modality, read_id)
        if key in inventory:
            raise ReadFateError(f'Duplicate normalized input read ID for {sample_id}/{modality}: {read_id}')
        inventory.add(key)
    return inventory


def _parse_darlin_sidecar_record(record: Mapping[str, str | None], path: Path, line_number: int) -> tuple[str, str, str, str, str, str]:
    if None in record or any((value is None for value in record.values())):
        raise ReadFateError(f'DARLIN sidecar {path}:{line_number} does not match its fixed columns')
    values = {column: str(record[column]).strip() for column in DARLIN_SIDECAR_COLUMNS}
    read_id = values['read_id']
    if normalize_read_id(read_id) != read_id:
        raise ReadFateError(f'DARLIN sidecar {path}:{line_number} read_id is not normalized: {read_id!r}')
    template = values['template']
    if not template:
        raise ReadFateError(f'DARLIN sidecar {path}:{line_number} has an empty template')
    fate = values['fate']
    if fate not in DARLIN_SIDECAR_FATES:
        raise ReadFateError(f'DARLIN sidecar {path}:{line_number} has invalid fate {fate!r}')
    cell_barcode = values['cell_barcode']
    clone_id = values['clone_id']
    detail = values['detail']
    if not detail:
        raise ReadFateError(f'DARLIN sidecar {path}:{line_number} has an empty detail')
    if fate == 'DARLIN_UNASSIGNED' and (cell_barcode or clone_id):
        raise ReadFateError(f'DARLIN sidecar {path}:{line_number} gives an unassigned read cell/clone evidence')
    if fate in {'DARLIN_ASSIGNED', 'DARLIN_FAILED_QC'} and (not (cell_barcode and clone_id)):
        raise ReadFateError(f'DARLIN sidecar {path}:{line_number} fate {fate} requires both cell_barcode and clone_id')
    return (read_id, template, fate, cell_barcode, clone_id, detail)


def _populate_darlin_sidecar_store(store: SqliteFateStore, sample_id: str, raw_r1: Path, raw_r2: Path, sidecars: Sequence[Path], templates: Sequence[str]=()) -> None:
    """Stream complete template evidence into SQLite and emit one fate per read."""
    _validate_sample_and_modality(sample_id, 'DARLIN')
    if not sidecars:
        raise ReadFateError('At least one DARLIN per-template sidecar is required')
    if templates and len(templates) != len(sidecars):
        raise ReadFateError('DARLIN template names must have one entry per sidecar')
    if templates and (any((not str(template).strip() for template in templates)) or len(set(templates)) != len(templates)):
        raise ReadFateError('DARLIN template names must be non-blank and unique')
    store.connection.executescript('\n        CREATE TABLE IF NOT EXISTS darlin_template_sources (\n            sample_id TEXT NOT NULL,\n            template TEXT NOT NULL,\n            path TEXT NOT NULL,\n            PRIMARY KEY (sample_id, template)\n        ) WITHOUT ROWID;\n        CREATE TABLE IF NOT EXISTS darlin_template_evidence (\n            sample_id TEXT NOT NULL,\n            read_id TEXT NOT NULL,\n            template TEXT NOT NULL,\n            fate TEXT NOT NULL,\n            cell_barcode TEXT NOT NULL,\n            clone_id TEXT NOT NULL,\n            detail TEXT NOT NULL,\n            PRIMARY KEY (sample_id, template, read_id)\n        ) WITHOUT ROWID;\n        ')
    for read_id in paired_fastq_ids(raw_r1, raw_r2):
        store.insert(FateRow(ReadKey(sample_id, 'DARLIN', read_id), 'DARLIN_UNASSIGNED', 'raw_fastq_inventory', ''))
    inventory_count = int(store.connection.execute("\n            SELECT COUNT(*) FROM fates\n            WHERE sample_id = ? AND modality = 'DARLIN'\n            ", (sample_id,)).fetchone()[0])
    for sidecar_index, path in enumerate(sidecars):
        expected_template = str(templates[sidecar_index]).strip() if templates else None
        template_name: str | None = expected_template
        if template_name is not None:
            try:
                store.connection.execute('INSERT INTO darlin_template_sources VALUES (?, ?, ?)', (sample_id, template_name, str(path)))
            except sqlite3.IntegrityError as exc:
                raise ReadFateError(f'Duplicate DARLIN template sidecar for sample {sample_id!r}, template {template_name!r}') from exc
        try:
            with _open_text(path) as handle:
                reader = csv.DictReader(handle, delimiter='\t')
                if tuple(reader.fieldnames or ()) != DARLIN_SIDECAR_COLUMNS:
                    raise ReadFateError(f'DARLIN sidecar {path} has columns {reader.fieldnames}; expected {DARLIN_SIDECAR_COLUMNS}')
                for line_number, record in enumerate(reader, start=2):
                    read_id, template, fate, cell_barcode, clone_id, detail = _parse_darlin_sidecar_record(record, path, line_number)
                    if expected_template is not None and template != expected_template:
                        raise ReadFateError(f'DARLIN sidecar {path} declares template {template!r}; expected {expected_template!r}')
                    if template_name is None:
                        template_name = template
                        try:
                            store.connection.execute('INSERT INTO darlin_template_sources VALUES (?, ?, ?)', (sample_id, template_name, str(path)))
                        except sqlite3.IntegrityError as exc:
                            raise ReadFateError(f'Duplicate DARLIN template sidecar for sample {sample_id!r}, template {template_name!r}') from exc
                    elif template != template_name:
                        raise ReadFateError(f'DARLIN sidecar {path} mixes templates {template_name!r} and {template!r}')
                    key = ReadKey(sample_id, 'DARLIN', read_id)
                    if store.get(key) is None:
                        raise ReadFateError(f'DARLIN sidecar evidence is absent from input inventory: {key}')
                    try:
                        store.connection.execute('\n                            INSERT INTO darlin_template_evidence\n                            VALUES (?, ?, ?, ?, ?, ?, ?)\n                            ', (sample_id, read_id, template_name, fate, cell_barcode, clone_id, detail))
                    except sqlite3.IntegrityError as exc:
                        raise ReadFateError(f'Duplicate DARLIN sidecar read_id {read_id!r} for template {template_name!r} in {path}:{line_number}') from exc
        except (OSError, UnicodeError, csv.Error) as exc:
            raise ReadFateError(f'Cannot read DARLIN sidecar {path}: {exc}') from exc
        if template_name is None:
            raise ReadFateError(f'DARLIN sidecar {path} contains no evidence rows')
        missing_count = int(store.connection.execute("\n                SELECT COUNT(*)\n                FROM fates AS f\n                LEFT JOIN darlin_template_evidence AS e\n                  ON e.sample_id = f.sample_id\n                 AND e.read_id = f.read_id\n                 AND e.template = ?\n                WHERE f.sample_id = ?\n                  AND f.modality = 'DARLIN'\n                  AND e.read_id IS NULL\n                ", (template_name, sample_id)).fetchone()[0])
        if missing_count:
            missing = [str(read_id) for read_id, in store.connection.execute("\n                    SELECT f.read_id\n                    FROM fates AS f\n                    LEFT JOIN darlin_template_evidence AS e\n                      ON e.sample_id = f.sample_id\n                     AND e.read_id = f.read_id\n                     AND e.template = ?\n                    WHERE f.sample_id = ?\n                      AND f.modality = 'DARLIN'\n                      AND e.read_id IS NULL\n                    ORDER BY f.read_id\n                    LIMIT 5\n                    ", (template_name, sample_id))]
            raise ReadFateError(f'Missing DARLIN sidecar evidence for {missing_count} input read pair(s) in template {template_name!r}: ' + ', '.join(missing))
    template_count = int(store.connection.execute('\n            SELECT COUNT(*) FROM darlin_template_sources WHERE sample_id = ?\n            ', (sample_id,)).fetchone()[0])
    aggregated_count = 0
    evidence_cursor = store.connection.execute('\n        SELECT read_id, template, fate, cell_barcode, clone_id, detail\n        FROM darlin_template_evidence\n        WHERE sample_id = ?\n        ORDER BY read_id, template\n        ', (sample_id,))
    for read_id, grouped_rows in groupby(evidence_cursor, key=lambda row: row[0]):
        template_evidence = [{'template': str(template), 'fate': str(fate), 'cell_barcode': str(cell_barcode) or None, 'clone_id': str(clone_id) or None, 'detail': str(detail) or None} for _, template, fate, cell_barcode, clone_id, detail in grouped_rows]
        if len(template_evidence) != template_count:
            raise ReadFateError(f'DARLIN template evidence count changed during aggregation for {read_id!r}')
        assigned_templates = [str(item['template']) for item in template_evidence if item['fate'] == 'DARLIN_ASSIGNED']
        if len(assigned_templates) == 1:
            terminal_fate = 'DARLIN_ASSIGNED'
            reason = 'single_assigned_template'
        elif len(assigned_templates) > 1:
            terminal_fate = 'DARLIN_FAILED_QC'
            reason = 'ambiguous_templates'
        elif any((item['fate'] == 'DARLIN_FAILED_QC' for item in template_evidence)):
            terminal_fate = 'DARLIN_FAILED_QC'
            reason = 'template_failed_qc'
        else:
            terminal_fate = 'DARLIN_UNASSIGNED'
            reason = 'all_templates_unassigned'
        key = ReadKey(sample_id, 'DARLIN', str(read_id))
        store.replace(key, 'DARLIN_UNASSIGNED', FateRow(key, terminal_fate, 'darlin_template_sidecar_consensus', json.dumps({'reason': reason, 'assigned_templates': assigned_templates, 'template_count': template_count, 'template_evidence': template_evidence}, sort_keys=True, separators=(',', ':'))))
        aggregated_count += 1
    if aggregated_count != inventory_count:
        raise ReadFateError(f'DARLIN aggregation did not conserve the raw inventory: expected {inventory_count}, aggregated {aggregated_count}')


def _load_rna_whitelist(path: Path) -> tuple[set[bytes], int]:
    """读取并严格验证原始 RNA barcode whitelist。"""
    barcodes: set[bytes] = set()
    barcode_length = None
    with path.open('rb') as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip(b'\r\n')
            if not line:
                raise ReadFateError(f'RNA whitelist {path}:{line_number} contains a blank line')
            fields = line.split(b'\t')
            if len(fields) not in {1, 2}:
                raise ReadFateError(f'RNA whitelist {path}:{line_number} must have one or two tab-separated columns')
            barcode = fields[0]
            if any((base not in b'ACGT' for base in barcode)):
                raise ReadFateError(f'RNA whitelist {path}:{line_number} first column is not raw uppercase A/C/G/T')
            if barcode_length is None:
                barcode_length = len(barcode)
            elif len(barcode) != barcode_length:
                raise ReadFateError(f'RNA whitelist {path} contains mixed barcode lengths')
            if barcode in barcodes:
                raise ReadFateError(f'RNA whitelist {path} contains duplicate barcode {barcode.decode()}')
            barcodes.add(barcode)
    if not barcodes or barcode_length not in {8, 10}:
        raise ReadFateError(f'RNA whitelist {path} must contain one or more uniform 8bp or 10bp barcodes')
    return (barcodes, barcode_length)


def _populate_rna_store(store: SqliteFateStore, sample_id: str, barcode_read: Path, insert_read: Path, whitelist_path: Path) -> None:
    """将逐 read 的 RNA whitelist 判定写入既有 store。"""
    _validate_sample_and_modality(sample_id, 'RNA')
    whitelist, barcode_length = _load_rna_whitelist(whitelist_path)
    minimum_layout_length = barcode_length + 8
    for read_id, barcode_sequence, _ in paired_fastq_records(barcode_read, insert_read):
        key = ReadKey(sample_id, 'RNA', read_id)
        barcode = barcode_sequence[:barcode_length]
        if len(barcode_sequence) < minimum_layout_length:
            store.insert(FateRow(key, 'RNA_UNASSIGNED', 'rna_barcode_gate', 'short_CxN8_layout'))
        elif barcode in whitelist:
            store.insert(FateRow(key, 'RNA_ASSIGNED', 'rna_barcode_gate', barcode.decode('ascii')))
        else:
            store.insert(FateRow(key, 'RNA_UNASSIGNED', 'rna_barcode_gate', barcode.decode('ascii', errors='replace')))


@dataclass(frozen=True)
class RnaBarcodeSourceAdapter:
    sample_id: str
    barcode_read: Path
    insert_read: Path
    whitelist: Path
    modality: str = 'RNA'
    source_name: str = 'rna_barcode_gate'

    def populate(self, store: SqliteFateStore) -> None:
        _populate_rna_store(store, self.sample_id, self.barcode_read, self.insert_read, self.whitelist)


@dataclass(frozen=True)
class DarlinTemplateSidecarSourceAdapter:
    """Combine complete per-template DARLIN evidence over one raw inventory."""
    sample_id: str
    raw_r1: Path
    raw_r2: Path
    sidecars: Sequence[Path]
    templates: Sequence[str] = ()
    modality: str = 'DARLIN'
    source_name: str = 'darlin_template_sidecar_consensus'

    def populate(self, store: SqliteFateStore) -> None:
        _populate_darlin_sidecar_store(store, self.sample_id, self.raw_r1, self.raw_r2, self.sidecars, self.templates)


def _insert_unique(rows: dict[ReadKey, FateRow], row: FateRow) -> None:
    previous = rows.get(row.key)
    if previous is not None:
        if previous.fate != row.fate:
            raise ReadFateError(f'Conflicting fates for {row.key}: {previous.fate} vs {row.fate}')
        raise ReadFateError(f'Duplicate fate evidence for {row.key}: {row.fate}')
    _validate_fate(row.key.modality, row.fate)
    rows[row.key] = row


def read_evidence_tsv(path: Path, sample_id: str, modality: str) -> list[FateRow]:
    """Read explicit per-read evidence; aggregate count files are not accepted."""
    _validate_sample_and_modality(sample_id, modality)
    with path.open(newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if reader.fieldnames is None or not {'read_id', 'fate'}.issubset(reader.fieldnames):
            raise ReadFateError(f'Evidence TSV {path} must contain read_id and fate columns')
        rows = []
        for line_number, record in enumerate(reader, start=2):
            row_sample = (record.get('sample_id') or sample_id).strip()
            row_modality = (record.get('modality') or modality).strip().upper()
            if row_sample != sample_id or row_modality != modality:
                raise ReadFateError(f'Evidence TSV {path} line {line_number} changes inventory identity')
            read_id = normalize_read_id(record.get('read_id', ''))
            fate = (record.get('fate') or '').strip()
            source = (record.get('source') or path.name).strip()
            detail = (record.get('detail') or '').strip()
            rows.append(FateRow(ReadKey(sample_id, modality, read_id), fate, source, detail))
    return rows


def collect_explicit(inventory: set[ReadKey], evidence_rows: Iterable[FateRow]) -> list[FateRow]:
    """Require exactly one explicit fate row for every inventory key."""
    evidence: dict[ReadKey, FateRow] = {}
    for row in evidence_rows:
        if row.key not in inventory:
            raise ReadFateError(f'Fate evidence is absent from input inventory: {row.key}')
        _insert_unique(evidence, row)
    missing = sorted(inventory.difference(evidence))
    if missing:
        preview = ', '.join((key.read_id for key in missing[:5]))
        raise ReadFateError(f'Missing fate evidence for {len(missing)} input read pair(s): {preview}')
    return [evidence[key] for key in sorted(evidence)]


def iter_ledger(path: Path) -> Iterator[FateRow]:
    with path.open(newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if tuple(reader.fieldnames or ()) != TSV_COLUMNS:
            raise ReadFateError(f'Read-fate TSV {path} has columns {reader.fieldnames}; expected {TSV_COLUMNS}')
        for record in reader:
            key = ReadKey(record['sample_id'], record['modality'], normalize_read_id(record['read_id']))
            yield FateRow(key, record['fate'], record['source'], record['detail'])


def validate_ledger_bundle(tsv: Path, summary_json: Path, checksum_path: Path) -> int:
    """Validate one upstream bundle and return its declared row count."""
    for label, path in (('TSV', tsv), ('summary', summary_json), ('checksum', checksum_path)):
        if not path.is_file():
            raise ReadFateError(f'Read-fate bundle {label} does not exist: {path}')
    digest = sha256_file(tsv)
    expected_checksum = f'{digest}  {tsv.name}\n'
    try:
        observed_checksum = checksum_path.read_text(encoding='utf-8')
        summary = json.loads(summary_json.read_text(encoding='utf-8'))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ReadFateError(f'Cannot read read-fate bundle for {tsv}: {exc}') from exc
    if observed_checksum != expected_checksum:
        raise ReadFateError(f'Read-fate checksum sentinel does not match {tsv.name}')
    if not isinstance(summary, Mapping):
        raise ReadFateError(f'Read-fate summary is not an object: {summary_json}')
    tsv_summary = summary.get('tsv')
    conservation = summary.get('conservation')
    if not isinstance(tsv_summary, Mapping) or not isinstance(conservation, Mapping):
        raise ReadFateError(f'Read-fate summary lacks tsv/conservation objects: {summary_json}')
    if summary.get('schema_version') != SCHEMA_VERSION or summary.get('kind') != 'rna_darlin_read_fate_summary' or summary.get('coverage_status') != 'complete' or (tsv_summary.get('path') != tsv.name) or (tsv_summary.get('columns') != list(TSV_COLUMNS)) or (tsv_summary.get('sha256') != digest) or (conservation.get('status') != 'passed') or (conservation.get('conserved') is not True):
        raise ReadFateError(f'Read-fate summary does not certify a complete conserved bundle: {summary_json}')
    sample_counts = summary.get('sample_counts')
    modality_counts = summary.get('modality_counts')
    fate_counts = summary.get('fate_counts')
    branch_conservation = summary.get('branch_conservation')
    source_status = summary.get('source_status')
    selected_modalities = summary.get('selected_modalities')
    if not all((isinstance(value, Mapping) for value in (sample_counts, modality_counts, fate_counts, branch_conservation, source_status))) or not isinstance(selected_modalities, list):
        raise ReadFateError(f'Read-fate summary lacks count/scope objects: {summary_json}')
    if any((isinstance(value, bool) or not isinstance(value, int) or value < 0 for counts in (sample_counts, modality_counts, fate_counts) for value in counts.values())):
        raise ReadFateError(f'Read-fate summary contains invalid counts: {summary_json}')
    row_count = conservation.get('output_rows')
    if isinstance(row_count, bool) or not isinstance(row_count, int) or row_count < 0:
        raise ReadFateError(f'Read-fate summary contains invalid counts: {summary_json}')
    wired = {modality for modality, status in source_status.items() if isinstance(status, Mapping) and status.get('status') == 'wired'}
    if selected_modalities != sorted(wired) or set(source_status) != set(FATES_BY_MODALITY) or any((not isinstance(source_status.get(modality), Mapping) or source_status[modality].get('status') != ('wired' if modality in wired else 'not_selected') for modality in FATES_BY_MODALITY)) or (sum(sample_counts.values()) != row_count) or (sum(modality_counts.values()) != row_count) or (sum(fate_counts.values()) != row_count) or (conservation.get('input_read_pairs') != row_count) or (conservation.get('output_rows') != row_count):
        raise ReadFateError(f'Read-fate summary counts/scope do not match {tsv.name}')
    for modality in FATES_BY_MODALITY:
        branch = branch_conservation.get(modality)
        if not isinstance(branch, Mapping):
            raise ReadFateError(f'Read-fate summary lacks {modality} branch conservation')
        count = modality_counts.get(modality, 0)
        if modality in wired:
            expected_branch = ('passed', count, count, True)
        else:
            expected_branch = ('not_assessed', None, None, None)
        observed_branch = (branch.get('status'), branch.get('input_read_pairs'), branch.get('output_rows'), branch.get('conserved'))
        if observed_branch != expected_branch:
            raise ReadFateError(f'Read-fate summary {modality} branch does not match {tsv.name}')
    return row_count


def write_outputs(rows: Sequence[FateRow], output_tsv: Path, summary_json: Path, checksum_path: Path, source_status: Mapping[str, Mapping[str, str]] | None=None) -> dict[str, object]:
    """Write an in-memory conformance result through the production disk store."""
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.read-fate-store-', dir=project_temporary_parent(output_tsv)) as temporary:
        with SqliteFateStore(Path(temporary) / 'fates.sqlite') as store:
            for row in rows:
                store.insert(row)
            return write_store_outputs(store, output_tsv, summary_json, checksum_path, source_status=source_status)


def write_store_outputs(store: SqliteFateStore, output_tsv: Path, summary_json: Path, checksum_path: Path, source_status: Mapping[str, Mapping[str, str]] | None=None) -> dict[str, object]:
    """Atomically stream stable TSV, checksum, and summary from a disk index."""
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    summary_json.parent.mkdir(parents=True, exist_ok=True)
    checksum_path.parent.mkdir(parents=True, exist_ok=True)
    with atomic_output(checksum_path) as checksum_tmp, atomic_output(summary_json) as summary_tmp, atomic_output(output_tsv) as tsv_tmp:
        with tsv_tmp.open('w', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=TSV_COLUMNS, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows((row.as_dict() for row in store.rows()))
        digest = sha256_file(tsv_tmp)
        counts = store.count_by('fate')
        sample_counts = store.count_by('sample_id')
        modality_counts = store.count_by('modality')
        total_rows = store.total()
        statuses = {modality: {'status': 'wired' if modality_counts.get(modality, 0) else 'not_selected', 'detail': 'per-read evidence present' if modality_counts.get(modality, 0) else 'modality was not selected for this ledger'} for modality in sorted(FATES_BY_MODALITY)}
        if source_status:
            for modality, status in source_status.items():
                if modality not in FATES_BY_MODALITY:
                    raise ReadFateError(f'Source status uses unsupported modality: {modality}')
                statuses[modality] = dict(status)
        wired_modalities = {modality for modality, status in statuses.items() if status.get('status') == 'wired'}
        selected_modalities = {modality for modality, status in statuses.items() if status.get('status') in {'wired', 'not_wired'}}
        invalid_statuses = {modality: status.get('status') for modality, status in statuses.items() if status.get('status') not in {'wired', 'not_wired', 'not_selected'}}
        if invalid_statuses:
            raise ReadFateError(f'Invalid source status values: {invalid_statuses}')
        for modality, count in modality_counts.items():
            if count and modality not in wired_modalities:
                raise ReadFateError(f'Read-fate rows exist for {modality}, but its source is marked not_wired')
        if selected_modalities and wired_modalities == selected_modalities:
            coverage_status = 'complete'
            global_conservation_status = 'passed'
            globally_conserved = True
        elif wired_modalities:
            coverage_status = 'partial'
            global_conservation_status = 'partial'
            globally_conserved = None
        else:
            coverage_status = 'not_assessed'
            global_conservation_status = 'not_assessed'
            globally_conserved = None
        branch_conservation = {}
        for modality in sorted(FATES_BY_MODALITY):
            count = modality_counts.get(modality, 0)
            wired = modality in wired_modalities
            branch_conservation[modality] = {'status': 'passed' if wired else 'not_assessed', 'input_read_pairs': count if wired else None, 'output_rows': count if wired else None, 'conserved': True if wired else None}
        summary: dict[str, object] = {'schema_version': SCHEMA_VERSION, 'kind': 'rna_darlin_read_fate_summary', 'coverage_status': coverage_status, 'selected_modalities': sorted(selected_modalities), 'tsv': {'path': output_tsv.name, 'columns': list(TSV_COLUMNS), 'sha256': digest}, 'conservation': {'status': global_conservation_status, 'input_read_pairs': total_rows, 'output_rows': total_rows, 'missing_fates': 0, 'duplicate_fates': 0, 'conflicting_fates': 0, 'mate_id_mismatches': 0, 'conserved': globally_conserved}, 'branch_conservation': branch_conservation, 'sample_counts': sample_counts, 'modality_counts': modality_counts, 'fate_counts': counts, 'source_status': statuses}
        summary_tmp.write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n', encoding='utf-8')
        checksum_tmp.write_text(f'{digest}  {output_tsv.name}\n', encoding='utf-8')
        checksum_path.unlink(missing_ok=True)
        summary_json.unlink(missing_ok=True)
    return summary


def _output_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument('--output-tsv', required=True, type=Path)
    parser.add_argument('--summary-json', required=True, type=Path)
    parser.add_argument('--checksum', required=True, type=Path)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest='command', required=True)
    explicit = subparsers.add_parser('collect', help='Collect explicit per-read evidence.')
    explicit.add_argument('--sample-id', required=True)
    explicit.add_argument('--modality', required=True, choices=sorted(FATES_BY_MODALITY))
    explicit.add_argument('--raw-r1', required=True, type=Path)
    explicit.add_argument('--raw-r2', required=True, type=Path)
    explicit.add_argument('--evidence', required=True, nargs='+', type=Path)
    _output_arguments(explicit)
    rna = subparsers.add_parser('collect-rna', help='Collect raw RNA barcode-gate fates.')
    rna.add_argument('--sample-id', required=True)
    rna.add_argument('--barcode-read', required=True, type=Path, help='FASTQ read containing the cell barcode and UMI (Cabernet raw R2).')
    rna.add_argument('--insert-read', required=True, type=Path, help='Paired insert FASTQ used for mate/header conservation (Cabernet raw R1).')
    rna.add_argument('--whitelist', required=True, type=Path)
    _output_arguments(rna)
    darlin = subparsers.add_parser('collect-darlin', help='Collect complete per-template DARLIN sidecars against raw paired FASTQ.')
    darlin.add_argument('--sample-id', required=True)
    darlin.add_argument('--raw-r1', required=True, type=Path)
    darlin.add_argument('--raw-r2', required=True, type=Path)
    darlin.add_argument('--sidecar', required=True, nargs='+', type=Path)
    darlin.add_argument('--template', required=True, nargs='+', help='Template names aligned one-to-one with --sidecar, including zero-read runs.')
    _output_arguments(darlin)
    merge = subparsers.add_parser('merge', help='Merge validated read-fate ledgers.')
    merge.add_argument('--inputs', nargs='*', type=Path, default=[])
    merge.add_argument('--input-summary', nargs='*', type=Path, default=[])
    merge.add_argument('--input-checksum', nargs='*', type=Path, default=[])
    merge.add_argument('--wired-modality', action='append', choices=sorted(FATES_BY_MODALITY), default=[], help='Modality whose production collector was wired, including a valid zero-read run.')
    _output_arguments(merge)
    return parser


def read_fate_main(argv: Sequence[str] | None=None) -> int:
    args = build_parser().parse_args(argv)
    source_status = None
    if args.command in {'collect', 'collect-rna', 'collect-darlin', 'merge'}:
        if args.command == 'collect':
            wired_modalities = {args.modality}
        elif args.command == 'collect-rna':
            wired_modalities = {'RNA'}
        elif args.command == 'collect-darlin':
            wired_modalities = {'DARLIN'}
        else:
            wired_modalities = set(args.wired_modality)
        source_status = {modality: {'status': 'not_selected', 'detail': f'{modality} was not selected for this ledger'} for modality in FATES_BY_MODALITY}
        if 'RNA' in wired_modalities:
            if args.command == 'collect':
                rna_status_detail = 'explicit terminal per-read RNA evidence supplied'
            else:
                rna_status_detail = 'raw paired FASTQ conserved against the configured exact C{barcode}N8 whitelist gate'
            source_status['RNA'] = {'status': 'wired', 'detail': rna_status_detail}
        if 'DARLIN' in wired_modalities:
            source_status['DARLIN'] = {'status': 'wired', 'detail': 'raw paired FASTQ conserved against complete per-template DARLIN sidecars with deterministic cross-template consensus'}
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.read-fate-store-', dir=project_temporary_parent(args.output_tsv)) as temporary:
        with SqliteFateStore(Path(temporary) / 'fates.sqlite') as store:
            if args.command == 'collect':
                inventory = inventory_from_fastq(args.sample_id, args.modality, args.raw_r1, args.raw_r2)
                evidence = []
                for path in args.evidence:
                    evidence.extend(read_evidence_tsv(path, args.sample_id, args.modality))
                for row in collect_explicit(inventory, evidence):
                    store.insert(row)
            elif args.command == 'collect-rna':
                RnaBarcodeSourceAdapter(args.sample_id, args.barcode_read, args.insert_read, args.whitelist).populate(store)
            elif args.command == 'collect-darlin':
                DarlinTemplateSidecarSourceAdapter(args.sample_id, args.raw_r1, args.raw_r2, args.sidecar, args.template).populate(store)
            else:
                if not len(args.inputs) == len(args.input_summary) == len(args.input_checksum):
                    raise ReadFateError('merge requires one --input-summary and --input-checksum for every --inputs ledger')
                for path, summary_path, checksum_path in zip(args.inputs, args.input_summary, args.input_checksum, strict=True):
                    expected_rows = validate_ledger_bundle(path, summary_path, checksum_path)
                    observed_rows = store.insert_many(iter_ledger(path))
                    if observed_rows != expected_rows:
                        raise ReadFateError(f'Read-fate summary counts/scope do not match {path.name}')
            write_store_outputs(store, args.output_tsv, args.summary_json, args.checksum, source_status=source_status)
    return 0
