"""Lemmus：基础文件哈希与 FASTQ 记录原语。"""
from __future__ import annotations

import gzip
import csv
import hashlib
import os
import errno
import uuid
from contextlib import contextmanager
from pathlib import Path


def iter_tree_files(root: Path, *, exclude_names=()):
    """完整枚举非目录成员，不跟随目录链接；扫描失败不能被当作空分支。"""
    pending = [Path(root)]
    excluded = frozenset(exclude_names)
    while pending:
        directory = pending.pop()
        members = []
        with os.scandir(directory) as entries:
            for entry in entries:
                if entry.name in excluded:
                    continue
                path = Path(entry.path)
                if entry.is_dir(follow_symlinks=False):
                    pending.append(path)
                else:
                    members.append(path)
        yield from members


def validate_delimited_records(handle, *, source, sep=','):
    """在pandas改写列名或猜索引前核对表头、引号与记录宽度。"""
    reader = csv.reader(handle, delimiter=sep, strict=True)
    try:
        fields = next(reader, [])
        if (not fields or any(not field or field != field.strip() for field in fields)
                or len({field.casefold() for field in fields}) != len(fields)):
            raise ValueError(f'{source}: header has blank, duplicate or padded columns')
        for row in reader:
            if row and len(row) != len(fields):
                raise ValueError(f'{source}:{reader.line_num}: field count differs from its header')
    except csv.Error as exc:
        raise ValueError(f'{source}:{reader.line_num}: malformed delimited table: {exc}') from exc
    return fields


def read_delimited_frame(path, *, sep=',', dtype=str):
    """保留原始身份与显式空值。"""
    import pandas as pd

    path = Path(path)
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, mode='rt', encoding='utf-8-sig', newline='') as handle:
        validate_delimited_records(handle, source=path, sep=sep)
        handle.seek(0)
        return pd.read_csv(handle, sep=sep, dtype=dtype, keep_default_na=False)


def validate_file_outputs(inputs, outputs):
    """写入前拒绝输入覆盖及输出互相包含，含链接与缺失父目录的规范路径。"""
    sources = [Path(path).resolve() for path in inputs]
    destinations = [Path(path).resolve() for path in outputs]
    for index, destination in enumerate(destinations):
        for other in sources + destinations[:index]:
            if destination.is_relative_to(other) or other.is_relative_to(destination):
                raise ValueError(f"Input/output paths overlap: {destination} and {other}")


def project_temporary_parent(destination: str | Path, *, project_dir: str | Path | None = None) -> Path:
    """项目运行的工作/提交临时项写入05_tmp；无项目上下文的安装事务保留目标父目录。"""
    destination = Path(destination)
    absolute = Path(os.path.abspath(destination))
    project = Path(project_dir) if project_dir is not None else None
    if project is None:
        configured = os.environ.get("RNA_DARLIN_PROJECT_DIR")
        if configured:
            candidate = Path(os.path.abspath(configured))
            project = candidate
    if project is None:
        project = next((parent for parent in absolute.parents
                        if (parent / "00_config").is_dir()), None)
    if project is None:
        return destination.parent
    temporary = project / "05_tmp"
    if temporary.is_symlink() or (temporary.exists() and not temporary.is_dir()):
        raise ValueError(f"Project temporary directory must be a real directory: {temporary}")
    temporary.mkdir(parents=True, exist_ok=True)
    # 先拒绝跨挂载提交；不得让后续工具退化为复制覆盖已有科学结果。
    parent = absolute.parent
    while not parent.exists():
        parent = parent.parent
    if parent.stat().st_dev != temporary.stat().st_dev:
        raise OSError(errno.EXDEV, "05_tmp and destination must share a filesystem", str(destination))
    return temporary


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        while (chunk := handle.read(1024 * 1024)):
            digest.update(chunk)
    return digest.hexdigest()


def file_stat_identity(path: Path) -> tuple[int, int, int, int, int]:
    info = path.stat()
    return info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns


def assert_file_states_unchanged(states):
    """提交候选前复核本次读取窗口，不额外扫描原始FASTQ。"""
    for path, expected in states.items():
        try:
            observed = file_stat_identity(path)
        except OSError as exc:
            raise ValueError(f'Input changed or became unreadable before commit: {path}') from exc
        if observed != expected:
            raise ValueError(f'Input changed before commit: {path}')


def stable_file_identity(path: str | Path):
    """读取哈希时检测替换/改写；stat仅用于本次读取窗口，持久来源记录使用SHA256。"""
    path = Path(path)
    before = file_stat_identity(path)
    digest = sha256_file(path)
    if file_stat_identity(path) != before:
        raise ValueError(f'Input changed while recording identity: {path}')
    return digest, before


def open_text(path: Path):
    """Open plain or gzipped FASTQ text."""
    with path.open("rb") as handle:
        magic = handle.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="ascii", errors="replace")
    return path.open("rt", encoding="ascii", errors="replace")


def has_complete_record(path: Path) -> bool:
    """Return true when the file starts with one complete FASTQ record."""
    with open_text(path) as handle:
        record = [handle.readline() for _ in range(4)]
    if not all(record):
        return False
    return bool(record[1].strip()) and bool(record[3].strip())


FastqRecord = tuple[bytes, bytes, bytes, bytes]

_PHRED33_QUALITIES = bytes(range(33, 127))


def is_valid_fastq_quality(quality: str | bytes) -> bool:
    """Phred+33 grammar; an empty quality is valid for an empty sequence."""
    if isinstance(quality, bytes):
        return not quality.translate(None, _PHRED33_QUALITIES)
    return quality.isascii() and (not quality or quality.isprintable() and ' ' not in quality)


def open_binary(path: Path):
    """Open plain or gzipped FASTQ as bytes."""
    return gzip.open(path, "rb") if path.suffix == ".gz" else path.open("rb")


def read_fastq_record(handle) -> FastqRecord | None:
    """Read one FASTQ record and raise on truncation."""
    header = handle.readline()
    if not header:
        return None
    sequence = handle.readline()
    plus = handle.readline()
    quality = handle.readline()
    if not sequence or not plus or not quality:
        raise ValueError("Truncated FASTQ record")
    if not header.startswith(b'@') or not plus.startswith(b'+'):
        raise ValueError("Malformed FASTQ record")
    sequence_body = sequence.removesuffix(b'\n').removesuffix(b'\r')
    quality_body = quality.removesuffix(b'\n').removesuffix(b'\r')
    if len(sequence_body) != len(quality_body):
        raise ValueError("FASTQ sequence/quality length mismatch")
    if not is_valid_fastq_quality(quality_body):
        raise ValueError("FASTQ quality must contain ASCII 33..126")
    return header, sequence, plus, quality


@contextmanager
def atomic_output(path: str | Path, *, project_dir: str | Path | None = None):
    """独占创建本次提交文件，保留umask；失败保旧且仅清理本次候选。"""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = project_temporary_parent(path, project_dir=project_dir) / f".{path.name}.tmp.{uuid.uuid4().hex}"
    descriptor = os.open(temporary, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o666)
    try:
        os.close(descriptor)
        yield temporary
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
