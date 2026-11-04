"""Lemmus：基础文件哈希与 FASTQ 记录原语。"""
from __future__ import annotations

import gzip
import hashlib
from pathlib import Path


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
