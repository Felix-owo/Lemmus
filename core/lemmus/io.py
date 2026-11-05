"""Lemmus：基础文件哈希与 FASTQ 记录原语。"""
from __future__ import annotations

import gzip
import hashlib
import os
import errno
import uuid
from contextlib import contextmanager
from pathlib import Path


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
    return header, sequence, plus, quality


@contextmanager
def atomic_output(path: str | Path):
    """独占创建本次提交文件，保留umask；失败保旧且仅清理本次候选。"""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = project_temporary_parent(path) / f".{path.name}.tmp.{uuid.uuid4().hex}"
    descriptor = os.open(temporary, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o666)
    try:
        os.close(descriptor)
        yield temporary
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
