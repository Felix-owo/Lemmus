"""Lemmus：输入发现、manifest 与 barcode 身份。"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
import re
import sys
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
import yaml
from .paths import PIPELINE_ROOT
from .config import (
    CONFIG_TEMPLATE_YAML,
    DROPLET_WHITELIST_ROWS,
    DROPLET_WHITELIST_SHA256,
    MANIFEST_FIELDS,
    SAMPLE_ID_PATTERN,
    SAMPLE_MANIFEST_TEMPLATE,
    load_config,
    validate_config,
)
from .io import (
    has_complete_record,
    open_binary,
    read_fastq_record,
)


REQUIRED_COLUMNS = ('RNA_Barcode', 'PlateID', 'Cell_Order')


ALLOWED_BARCODE_LENGTHS = frozenset({8, 10})


class BarcodeConfigError(ValueError):
    """Raised when Barcode_Map.csv violates the standalone RNA contract."""


@dataclass(frozen=True)
class BarcodeAudit:
    lengths: tuple[int, ...]
    minimum_hamming_distance: int | None
    single_base_ambiguities: tuple[str, ...]


@dataclass(frozen=True)
class BarcodeEntry:
    rna_barcode: str
    plate_id: str
    cell_order: str
    line_number: int


@dataclass(frozen=True)
class BarcodeConfig:
    path: Path
    rows: tuple[BarcodeEntry, ...]
    rna_audit: BarcodeAudit

    @property
    def rna_barcodes(self) -> tuple[str, ...]:
        return tuple((row.rna_barcode for row in self.rows))


def _hamming(left: str, right: str) -> int:
    return sum((a != b for a, b in zip(left, right, strict=True)))


def _audit(barcodes: tuple[str, ...]) -> BarcodeAudit:
    distances = [_hamming(left, right) for index, left in enumerate(barcodes) for right in barcodes[index + 1:] if len(left) == len(right)]
    whitelist = set(barcodes)
    candidates: dict[str, set[str]] = {}
    for barcode in barcodes:
        for index, original in enumerate(barcode):
            for base in 'ACGT':
                if base == original:
                    continue
                observed = barcode[:index] + base + barcode[index + 1:]
                if observed not in whitelist:
                    candidates.setdefault(observed, set()).add(barcode)
    ambiguous = tuple(sorted((item for item, values in candidates.items() if len(values) > 1)))
    return BarcodeAudit(lengths=tuple(sorted({len(barcode) for barcode in barcodes})), minimum_hamming_distance=min(distances) if distances else None, single_base_ambiguities=ambiguous)


def load_barcode_config(path: str | Path, *, require_rna: bool=True) -> BarcodeConfig:
    barcode_path = Path(path).expanduser().resolve()
    try:
        handle = barcode_path.open(encoding='utf-8-sig', newline='')
    except OSError as exc:
        raise BarcodeConfigError(f'Cannot read Barcode_Map.csv: {barcode_path}') from exc
    with handle:
        reader = csv.DictReader(handle)
        fields = tuple(reader.fieldnames or ())
        missing = [field for field in REQUIRED_COLUMNS if field not in fields]
        if missing:
            raise BarcodeConfigError('Barcode_Map.csv is missing column(s): ' + ', '.join(missing))
        rows: list[BarcodeEntry] = []
        seen_barcodes: set[str] = set()
        seen_cells: set[tuple[str, str]] = set()
        expected_length: int | None = None
        for line_number, raw in enumerate(reader, start=2):
            barcode = str(raw.get('RNA_Barcode', '') or '').strip()
            plate_id = str(raw.get('PlateID', '') or '').strip()
            cell_order = str(raw.get('Cell_Order', '') or '').strip()
            if not barcode and (not require_rna):
                continue
            if set(barcode) - set('ACGT') or len(barcode) not in ALLOWED_BARCODE_LENGTHS:
                raise BarcodeConfigError(f'Barcode_Map.csv line {line_number} has invalid RNA_Barcode {barcode!r}')
            expected_length = len(barcode) if expected_length is None else expected_length
            if len(barcode) != expected_length:
                raise BarcodeConfigError('RNA barcodes must use one uniform length')
            if not plate_id or not cell_order:
                raise BarcodeConfigError(f'Barcode_Map.csv line {line_number} requires PlateID and Cell_Order')
            if barcode in seen_barcodes:
                raise BarcodeConfigError(f'Duplicate RNA_Barcode: {barcode}')
            if (plate_id, cell_order) in seen_cells:
                raise BarcodeConfigError(f'Duplicate PlateID/Cell_Order: {plate_id}/{cell_order}')
            seen_barcodes.add(barcode)
            seen_cells.add((plate_id, cell_order))
            rows.append(BarcodeEntry(barcode, plate_id, cell_order, line_number))
    if require_rna and (not rows):
        raise BarcodeConfigError('Barcode_Map.csv contains no RNA barcodes')
    barcodes = tuple((row.rna_barcode for row in rows))
    audit = _audit(barcodes)
    # Hard gate mirroring the DNA pipeline contract: a minimum pairwise Hamming
    # distance <=2 makes unique Hamming-1 correction impossible, so the map is
    # rejected up front instead of mis-correcting later in alevin-fry/DARLIN.
    if barcodes and audit.minimum_hamming_distance is not None and audit.minimum_hamming_distance <= 2:
        raise BarcodeConfigError(
            'Barcode_Map.csv minimum pairwise Hamming distance is '
            f'{audit.minimum_hamming_distance}; whitelist correction requires >2'
        )
    return BarcodeConfig(barcode_path, tuple(rows), audit)


def write_rna_whitelist(config: BarcodeConfig, output_path: str | Path) -> Path:
    output = Path(output_path).expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f'.{output.name}.', dir=output.parent, text=True)
    try:
        with os.fdopen(fd, 'w', encoding='utf-8', newline='\n') as handle:
            for barcode in config.rna_barcodes:
                handle.write(barcode + '\n')
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, output)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise
    return output


REQUIRED_MANIFEST_FIELDS = MANIFEST_FIELDS[:3]


FASTQ_PATTERNS = (
    re.compile(r'(?P<sample>.+)\.R(?P<read>[12])\.raw\.f(?:ast)?q\.gz$', re.I),
    re.compile(r'(?P<sample>.+)_L\d{3}_R(?P<read>[12])_\d{3}\.f(?:ast)?q\.gz$', re.I),
    re.compile(r'(?P<sample>.+)_R(?P<read>[12])(?:_\d{3})?\.f(?:ast)?q\.gz$', re.I),
    re.compile(r'(?P<sample>.+)_(?P<read>[12])\.f(?:ast)?q\.gz$', re.I),
)


def parse_fastq_name(path: Path) -> tuple[str | None, str | None]:
    for pattern in FASTQ_PATTERNS:
        match = pattern.match(path.name)
        if match:
            return (match.group('sample'), f"R{match.group('read')}")
    return (None, None)


def discover_fastq_pairs(raw_dir: str | Path) -> tuple[dict[str, dict[str, Path]], list[str]]:
    root_dir = Path(raw_dir)
    pairs: dict[str, dict[str, Path]] = {}
    warnings: list[str] = []
    invalid: set[str] = set()
    if not root_dir.exists():
        return ({}, [f'raw directory does not exist: {root_dir}'])
    for root, _, files in os.walk(root_dir):
        for filename in files:
            if not filename.lower().endswith(('.fastq.gz', '.fq.gz')):
                continue
            path = Path(root) / filename
            sample, read = parse_fastq_name(path)
            if sample is None or read is None:
                warnings.append(f'ignored FASTQ with unsupported name: {path}')
                continue
            entry = pairs.setdefault(sample, {})
            if read in entry and entry[read] != path:
                warnings.append(f'multiple {read} FASTQs for sample {sample}')
                invalid.add(sample)
            else:
                entry[read] = path
    complete: dict[str, dict[str, Path]] = {}
    for sample, reads in sorted(pairs.items()):
        if sample in invalid:
            continue
        if set(reads) == {'R1', 'R2'}:
            complete[sample] = {'R1': reads['R1'], 'R2': reads['R2']}
        else:
            warnings.append(f'incomplete FASTQ pair for sample {sample} under {root_dir}')
    return (complete, warnings)


def require_fastq_pairs(raw_dir: str | Path, label: str, *, allow_empty: bool=False) -> dict[str, dict[str, Path]]:
    pairs, warnings = discover_fastq_pairs(raw_dir)
    if warnings:
        raise ValueError(f'Invalid {label} FASTQ inventory: ' + '; '.join(warnings))
    if not pairs and (not allow_empty):
        raise ValueError(f'No paired {label} FASTQ files were discovered under {raw_dir}')
    return pairs


def read_manifest(path: str | Path) -> tuple[list[str], list[dict[str, str]]]:
    manifest = Path(path)
    with manifest.open(encoding='utf-8-sig', newline='') as handle:
        delimiter = '\t' if manifest.suffix.lower() in {'.tsv', '.txt'} else ','
        reader = csv.DictReader(handle, delimiter=delimiter)
        fields = [str(field).strip() for field in reader.fieldnames or []]
        missing = [field for field in REQUIRED_MANIFEST_FIELDS if field not in fields]
        if missing:
            raise ValueError('sample_manifest.tsv is missing column(s): ' + ', '.join(missing))
        rows = []
        sample_lines: dict[str, int] = {}
        for line_number, raw in enumerate(reader, start=2):
            row = {field: (raw.get(field) or '').strip() for field in fields}
            if not SAMPLE_ID_PATTERN.fullmatch(row['sample_id']):
                raise ValueError(f"sample_manifest.tsv line {line_number} has invalid sample_id {row['sample_id']!r}")
            previous_line = sample_lines.get(row['sample_id'])
            if previous_line is not None:
                raise ValueError(f"sample_manifest.tsv has duplicate sample_id {row['sample_id']!r} on lines {previous_line} and {line_number}")
            sample_lines[row['sample_id']] = line_number
            rows.append(row)
    return (fields, rows)


def map_fastq_pairs_to_manifest(pairs: Mapping[str, Mapping[str, Path]], rows: list[dict[str, str]], *, column: str, label: str) -> dict[str, dict[str, object]]:
    mapped: dict[str, dict[str, object]] = {}
    declared_raw: set[str] = set()
    for row in rows:
        raw_sample = row.get(column, '')
        if not raw_sample:
            continue
        sample_id = row['sample_id']
        if raw_sample in declared_raw or sample_id in mapped:
            raise ValueError(f'Duplicate {label} mapping in sample_manifest.tsv: {raw_sample}')
        declared_raw.add(raw_sample)
        if raw_sample not in pairs:
            raise ValueError(f'{label} sample lacks paired FASTQ files: {raw_sample}')
        mapped[sample_id] = {'R1': pairs[raw_sample]['R1'], 'R2': pairs[raw_sample]['R2'], 'raw_sample': raw_sample, 'sample_id': sample_id}
    extras = sorted(set(pairs) - declared_raw)
    if extras:
        raise ValueError(f'Paired {label} FASTQs have no manifest mapping: ' + ', '.join(extras))
    return mapped


def refresh(project: Path, *, dry_run: bool=False) -> list[dict[str, str]]:
    project = project.expanduser().resolve()
    manifest = project / '00_config' / 'sample_manifest.tsv'
    existing: list[dict[str, str]] = []
    if manifest.is_file():
        _, existing = read_manifest(manifest)
    by_rna = {row['rna_sample']: row for row in existing if row.get('rna_sample')}
    by_darlin = {row['darlin_sample']: row for row in existing if row.get('darlin_sample')}
    rna, rna_warnings = discover_fastq_pairs(project / '01_raw' / 'RNA')
    darlin, darlin_warnings = discover_fastq_pairs(project / '01_raw' / 'DARLIN')
    warnings = rna_warnings + darlin_warnings
    if warnings:
        raise ValueError('; '.join(warnings))
    rows: list[dict[str, str]] = []
    used: set[int] = set()
    for raw_sample in sorted(rna):
        row = by_rna.get(raw_sample)
        if row is None and raw_sample in by_darlin:
            row = by_darlin[raw_sample]
            row['rna_sample'] = raw_sample
        if row is None:
            row = {'sample_id': raw_sample, 'rna_sample': raw_sample, 'darlin_sample': raw_sample if raw_sample in darlin else '', 'species': '', 'notes': ''}
        used.add(id(row))
        rows.append({field: row.get(field, '') for field in MANIFEST_FIELDS})
    for raw_sample in sorted(darlin):
        if any((row['darlin_sample'] == raw_sample for row in rows)):
            continue
        row = by_darlin.get(raw_sample)
        if row is None and raw_sample in by_rna:
            row = by_rna[raw_sample]
            row['darlin_sample'] = raw_sample
        if row is None:
            row = {'sample_id': raw_sample, 'rna_sample': '', 'darlin_sample': raw_sample, 'species': '', 'notes': ''}
        used.add(id(row))
        rendered = {field: row.get(field, '') for field in MANIFEST_FIELDS}
        for index, existing_row in enumerate(rows):
            if existing_row['sample_id'] == rendered['sample_id']:
                rows[index] = rendered
                break
        else:
            rows.append(rendered)
    stale = [row['sample_id'] for row in existing if id(row) not in used]
    if stale:
        raise ValueError('Manifest contains samples absent from current RNA/DARLIN FASTQs: ' + ', '.join(stale))
    if dry_run:
        writer = csv.DictWriter(sys.stdout, fieldnames=MANIFEST_FIELDS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    else:
        temporary = manifest.with_name(f'.{manifest.name}.tmp')
        with temporary.open('w', encoding='utf-8', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)
        temporary.replace(manifest)
    return rows


def manifest_refresh_main(argv: list[str] | None=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument('--project', type=Path, required=True)
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args(argv)
    refresh(args.project, dry_run=args.dry_run)
    return 0


def droplet_whitelist_path() -> Path:
    """返回固定来源的 DD-MET5 设计白名单；不依赖项目 Barcode Map。"""
    return PIPELINE_ROOT / "resources/droplet/U3CB_methylation.txt.gz"


def barcode_map_projection(path: Path) -> tuple[tuple[str, str, str], ...]:
    """读取孔板身份投影；拒绝空表、重复身份和损坏行，保留身份字符串。"""
    columns = ('plateid', 'rna_barcode', 'cell_order')
    with path.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.reader(handle)
        header = [field.strip().lower() for field in next(reader, [])]
        if not header or any(not field for field in header) or len(set(header)) != len(header):
            raise ValueError(f'{path}: Barcode_Map 表头为空或重复')
        if not set(columns).issubset(header):
            raise ValueError(f'{path}: Barcode_Map 缺少 PlateID/RNA_Barcode/Cell_Order')
        indices = [header.index(field) for field in columns]
        rows = []
        seen = [set() for _ in columns]
        for number, raw in enumerate(reader, 2):
            if len(raw) != len(header):
                raise ValueError(f'{path}:{number}: Barcode_Map 行列数与表头不一致')
            row = tuple(raw[index].strip() for index in indices)
            if any(not value for value in row):
                raise ValueError(f'{path}:{number}: Barcode_Map 身份为空')
            if any(value in values for value, values in zip(row, seen)):
                raise ValueError(f'{path}:{number}: Barcode_Map 身份重复')
            for value, values in zip(row, seen):
                values.add(value)
            rows.append(row)
    if not rows:
        raise ValueError(f'{path}: Barcode_Map 没有数据行')
    return tuple(sorted(rows))


def check_cross_pipeline_alignment(pipeline_root: Path, dna_root: Path | None = None) -> None:
    """只读核对 Alopex 的 Cabernet 孔板映射与 Droplet 固定设计；显式错误路径必须失败。"""
    explicit = dna_root is not None
    dna_root = (dna_root if explicit else pipeline_root.parent / 'Alopex').expanduser().resolve()
    if not dna_root.is_dir():
        message = f'未找到 Alopex Pipeline ({dna_root})；Cabernet/Droplet 跨 Pipeline 对齐未验证'
        if explicit:
            raise ValueError(message)
        print(f'WARN: {message}（可用 DNA_PIPELINE_ROOT 指定）')
        return
    print(f'Alopex Pipeline: {dna_root}')
    rna_map = barcode_map_projection(pipeline_root / 'resources/Barcode_Map.csv')
    dna_map = barcode_map_projection(dna_root / 'resources/Barcode_Map.csv')
    if rna_map != dna_map:
        raise ValueError('Cabernet: Lemmus 与 Alopex 的 (PlateID, RNA_Barcode, Cell_Order) 不一致')
    print(f'Cabernet Barcode_Map 跨 Pipeline 一致性: PASS ({len(rna_map)} cells)')
    rna_design = pipeline_root / 'resources/droplet/U3CB_methylation.txt.gz'
    dna_design = dna_root / 'resources/droplet_ME5_U3CB_methylation.txt.gz'
    for path in (rna_design, dna_design):
        if path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest() != DROPLET_WHITELIST_SHA256:
            raise ValueError(f'Droplet: 固定 ME5 设计白名单 SHA256 不匹配: {path}')
    missing = [str(path) for path in (rna_design, dna_design) if not path.is_file()]
    if missing:
        print('WARN: Droplet 跨 Pipeline 对齐未验证；缺少固定设计白名单: ' + ', '.join(missing)
              + '；Lemmus 白名单由 core/doctor.sh calling 准备')
    else:
        print(f'Droplet ME5 设计白名单跨 Pipeline 一致性: PASS ({DROPLET_WHITELIST_ROWS} barcodes)')


def load_droplet_whitelist(path: Path) -> set[str]:
    """核对白名单内容和唯一性，避免对大型设计表进行全量两两距离检查。"""
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != DROPLET_WHITELIST_SHA256:
        raise ValueError("DD-MET5 design whitelist SHA256 mismatch")
    values = gzip.decompress(raw).decode('ascii').splitlines()
    whitelist = set(values)
    if (len(values) != DROPLET_WHITELIST_ROWS or len(whitelist) != len(values)
            or any(len(x) != 17 or set(x) - set('AGT') for x in values)):
        raise ValueError("Invalid DD-MET5 design whitelist")
    return whitelist


def _barcode_map_digest(path):
    """Return (sha256_hex, data_row_count) stamping which Barcode_Map built a result."""
    import hashlib
    digest = hashlib.sha256()
    row_count = 0
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        next(handle, None)
        for line in handle:
            if line.strip():
                row_count += 1
    return (digest.hexdigest(), row_count)


def assert_fastq_main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fastq", type=Path, help="FASTQ or FASTQ.gz to validate.")
    parser.add_argument("label", help="Human-readable label used in error messages.")
    args = parser.parse_args(argv)

    if not args.fastq.exists():
        parser.error(f"{args.label} does not exist: {args.fastq}")
    if not has_complete_record(args.fastq):
        raise SystemExit(f"{args.label} has no complete FASTQ records: {args.fastq}")
    print(f"{args.label} contains FASTQ records: {args.fastq}")
    return 0


CARLIN_5PRIME_SEEDS = (
    b"AGCTGTACAAGTAAGCGGC",  # cCARLIN
    b"TCGGTACCTCGCGAA",      # Tigre
    b"GTACAAGTAAAGCGGCC",    # Rosa
)


def _rna_prefilter_seed_hit(sequence, seed):
    """Cheap high-recall prefilter; the full indel-aware check happens later."""
    if seed in sequence:
        return True
    # Two separated conserved sub-seeds make one sequencing error/1-bp indel
    # unlikely to erase both anchors, while avoiding DP across every RNA read.
    flank = max(7, min(10, len(seed) // 2))
    return seed[:flank] in sequence or seed[-flank:] in sequence


def prefilter(args: argparse.Namespace) -> None:
    """Stream paired FASTQ and keep reads whose R1 contains a DARLIN seed."""
    input_r1 = Path(args.input_r1)
    input_r2 = Path(args.input_r2)
    output_r1 = Path(args.output_r1)
    output_r2 = Path(args.output_r2)
    metrics = Path(args.metrics)
    output_r1.parent.mkdir(parents=True, exist_ok=True)
    output_r2.parent.mkdir(parents=True, exist_ok=True)
    metrics.parent.mkdir(parents=True, exist_ok=True)

    total = 0
    kept = 0
    seed_counts = {seed: 0 for seed in CARLIN_5PRIME_SEEDS}

    with open_binary(input_r1) as r1_in, open_binary(input_r2) as r2_in, gzip.open(
        output_r1, "wb", compresslevel=1
    ) as r1_out, gzip.open(output_r2, "wb", compresslevel=1) as r2_out:
        while True:
            r1 = read_fastq_record(r1_in)
            r2 = read_fastq_record(r2_in)
            if r1 is None and r2 is None:
                break
            if r1 is None or r2 is None:
                raise ValueError("R1/R2 FASTQ record counts differ")

            total += 1
            sequence = r1[1].strip().upper()
            matched = [seed for seed in CARLIN_5PRIME_SEEDS if _rna_prefilter_seed_hit(sequence, seed)]
            if not matched:
                continue

            kept += 1
            for seed in matched:
                seed_counts[seed] += 1
            r1_out.write(b"".join(r1))
            r2_out.write(b"".join(r2))

    with metrics.open("w") as handle:
        handle.write("sample\ttotal_reads\tkept_reads\tkept_fraction\tseed\tseed_reads\n")
        fraction = kept / total if total else 0.0
        for seed, count in seed_counts.items():
            seed_name = seed.decode("ascii")
            handle.write(f"{args.sample}\t{total}\t{kept}\t{fraction:.8f}\t{seed_name}\t{count}\n")

    print(f"sample={args.sample}")
    print(f"total_reads={total}")
    print(f"kept_reads={kept}")
    print(f"kept_fraction={kept / total if total else 0.0:.8f}")


def prefilter_rna_darlin_main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-r1", required=True, help="Input RNA R1 FASTQ, plain or gzipped.")
    parser.add_argument("--input-r2", required=True, help="Input RNA R2 FASTQ, plain or gzipped.")
    parser.add_argument("--output-r1", required=True, help="Filtered R1 FASTQ.gz output.")
    parser.add_argument("--output-r2", required=True, help="Filtered R2 FASTQ.gz output.")
    parser.add_argument("--metrics", required=True, help="TSV metrics file describing kept DARLIN-seed reads.")
    parser.add_argument("--sample", required=True, help="Sample name written to the metrics file.")
    prefilter(parser.parse_args(argv))


def write_rna_whitelist_main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--barcode-map", type=Path)
    group.add_argument("--design-whitelist", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)

    if args.design_whitelist is not None:
        values = load_droplet_whitelist(args.design_whitelist)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(mode='w', dir=args.output.parent, delete=False) as handle:
            handle.write('\n'.join(sorted(values))+'\n')
        os.replace(handle.name, args.output)
        return 0

    config = load_barcode_config(args.barcode_map, require_rna=True)
    output = write_rna_whitelist(config, args.output)
    audit = config.rna_audit
    print(
        "RNA whitelist generated from project Barcode_Map.csv: "
        f"rows={len(config.rna_barcodes)} "
        f"lengths={','.join(str(length) for length in audit.lengths)} "
        f"minimum_raw_hamming={audit.minimum_hamming_distance} "
        f"single_base_ambiguities={len(audit.single_base_ambiguities)} "
        f"output={output}"
    )
    return 0


def init_project_main(argv=None):
    """初始化项目目录、当前模板与指向真实 Pipeline 的 launcher 软连接。"""
    parser = argparse.ArgumentParser(prog="lemmus init-project")
    parser.add_argument("--project", type=Path, default=Path.cwd())
    args = parser.parse_args(argv)
    project = args.project.resolve()
    config_path = project / "00_config/config.yaml"
    template_text = CONFIG_TEMPLATE_YAML
    validate_config(yaml.safe_load(template_text))
    if config_path.exists():
        existing = load_config(config_path)
        print(f"保留已有 {existing['pipeline_mode']} 配置: {config_path}")
    for relative in ("00_config", "01_raw/RNA", "01_raw/DARLIN", "02_work", "03_results", "04_logs", "05_tmp", "06_downstream"):
        (project / relative).mkdir(parents=True, exist_ok=True)
    defaults = {
        config_path: template_text,
        project / "00_config/sample_manifest.tsv": SAMPLE_MANIFEST_TEMPLATE,
        project / "00_config/Barcode_Map.csv": (PIPELINE_ROOT / "resources/Barcode_Map.csv").read_text(encoding="utf-8"),
    }
    for path, payload in defaults.items():
        if not path.exists():
            path.write_text(payload, encoding="utf-8")
    launcher = project / "run_pipeline.sh"
    pipeline_launcher = PIPELINE_ROOT / "core/run_pipeline.sh"
    if not launcher.exists() and not launcher.is_symlink():
        launcher.symlink_to(pipeline_launcher)
    elif not launcher.is_symlink():
        print(f"WARNING: {launcher} 已存在且不是软连接，未替换。")
    print(f"项目已初始化: {project}")
    print("下一步:")
    print("  1. 编辑 00_config/config.yaml 的 pipeline_mode（默认 cabernet，可选 droplet）及对应模块。")
    print("     Cabernet 检查 Barcode_Map.csv；RNA-only 设 cabernet.DARLIN.enabled=false。")
    print("     Droplet 使用 R1=CB17+UMI12、R2=cDNA；运行前用 core/doctor.sh calling 准备环境和 ME5 白名单。")
    print("     RNA FASTQ 放入 01_raw/RNA/，targeted DARLIN FASTQ 放入 01_raw/DARLIN/。")
    print("  2. ./run_pipeline.sh --refresh-manifest")
    print("  3. ./run_pipeline.sh --dry-run")
    print("  4. ./run_pipeline.sh")
    return 0
