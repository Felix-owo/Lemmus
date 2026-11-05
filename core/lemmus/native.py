"""固定 RNA 原生定量器的安装与只读契约检查。"""

from __future__ import annotations
import argparse
import hashlib
import json
import os
import platform
import shutil
import struct
import subprocess
import sys
import tempfile
from pathlib import Path
from .config import RNA_GENE_EC_ORDER
from .paths import PIPELINE_ROOT


def _distribution():
    root = PIPELINE_ROOT / "resources/alevin-fry"
    data = json.loads((root / "distribution.json").read_text())
    if (
        data.get("schema") != 1
        or data.get("version") != "0.18.3"
        or data.get("gene_level_ec_order") != RNA_GENE_EC_ORDER
    ):
        raise ValueError("Invalid native RNA distribution contract")
    machine = platform.machine().lower()
    key = (
        "osx-arm64"
        if sys.platform == "darwin" and machine == "arm64"
        else "linux-64"
        if sys.platform.startswith("linux") and machine == "x86_64"
        else None
    )
    if key is None or key not in data["platforms"]:
        raise ValueError("No validated native RNA distribution for this platform")
    record = data["platforms"][key]
    binary = (root / record["path"]).resolve(strict=True)
    if not binary.is_relative_to(root.resolve()):
        raise ValueError("Native RNA distribution path escapes its resource root")
    return binary, record["sha256"]


def _check_binary(binary: Path):
    result = subprocess.run(
        [str(binary), "lemmus-contract"],
        capture_output=True,
        text=True,
        check=True,
        timeout=30,
    )
    if result.stdout.strip() != RNA_GENE_EC_ORDER:
        raise ValueError(
            "Native RNA quantifier lacks canonical projected-gene identities"
        )


def check_quantifier(prefix: Path):
    """只执行轻量契约查询；不访问 Conda solver、网络或科学参考。"""
    _check_binary(prefix / "bin/alevin-fry")


def _smoke(binary: Path, scratch: Path):
    """五条已知 reads 反转顺序后仍为同一 PG 歧义结果，防止仅声明新契约。"""
    pack = lambda fmt, *values: struct.pack("<" + fmt, *values)
    tag = lambda name, kind: pack("H", len(name)) + name.encode() + pack("B", kind)
    header = pack("BQ", 0, 5)
    for i in range(5):
        name = "t" + str(i)
        header += pack("H", len(name)) + name.encode()
    header += pack("Q", 1)
    header += pack("H", 2) + tag("cblen", 2) + tag("ulen", 2)
    header += pack("H", 2) + tag("b", 4) + tag("u", 4)
    header += pack("H", 1) + tag("compressed_ori_refid", 3) + pack("HH", 17, 12)
    records = []
    for umi, genes, count in [
        (11179698, [0, 1, 3], 2),
        (2791090, [0, 2, 4], 1),
        (11179634, [1, 2, 3], 2),
    ]:
        record = pack("IQQ", len(genes), 2329, umi)
        record += b"".join(pack("I", gene) for gene in genes)
        records.extend([record] * count)
    mapping = scratch / "t2g.tsv"
    mapping.write_text("".join(f"t{i}\tg{i}\n" for i in range(5)))
    for name, reads in (("forward", records), ("reverse", list(reversed(records)))):
        source = scratch / name
        source.mkdir()
        content = b"".join(reads)
        (source / "map.collated.rad").write_bytes(
            header + pack("II", len(content) + 8, len(reads)) + content
        )
        (source / "collate.json").write_text("{}")
        (source / "generate_permit_list.json").write_text("{}")
        output = scratch / (name + "_out")
        subprocess.run(
            [
                str(binary),
                "quant",
                "-i",
                str(source),
                "-o",
                str(output),
                "-m",
                str(mapping),
                "-r",
                "parsimony-gene",
                "--small-thresh",
                "0",
                "-t",
                "2",
            ],
            check=True,
            capture_output=True,
            text=True,
            timeout=60,
        )
        metadata = json.loads((output / "quant.json").read_text())
        if (
            metadata.get("gene_level_ec_order") != RNA_GENE_EC_ORDER
            or metadata.get("resolution_strategy") != "ParsimonyGene"
            or metadata.get("usa_mode") is not False
        ):
            raise ValueError(
                "Native RNA smoke produced an incorrect scientific contract"
            )
        matrix = (output / "alevin/quants_mat.mtx").read_text().splitlines()
        values = [line.split() for line in matrix if line and not line.startswith("%")]
        if values != [["1", "5", "0"]]:
            raise ValueError(
                "Native RNA smoke changed the known ambiguous molecule result"
            )
        expected = "".join("ACGT"[(2329 >> (2 * i)) & 3] for i in reversed(range(17)))
        if (output / "alevin/quants_mat_rows.txt").read_text().splitlines() != [
            expected
        ]:
            raise ValueError("Native RNA smoke changed barcode identity")
        if (output / "alevin/quants_mat_cols.txt").read_text().splitlines() != [
            f"g{i}" for i in range(5)
        ]:
            raise ValueError("Native RNA smoke changed gene identity")


def prepare_quantifier(prefix: Path):
    """核验固定发行工件和真实原生行为后原子替换；失败保留原程序。"""
    try:
        check_quantifier(prefix)
        return
    except (OSError, ValueError, subprocess.SubprocessError):
        pass
    if sys.platform.startswith("linux") and os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError(
            "Prepare the native RNA quantifier on the login node: bash core/doctor.sh conda"
        )
    source, expected = _distribution()
    target = prefix / "bin/alevin-fry"
    fd, name = tempfile.mkstemp(prefix=".alevin-fry.", dir=target.parent)
    candidate = Path(name)
    try:
        with os.fdopen(fd, "wb") as writer, source.open("rb") as reader:
            shutil.copyfileobj(reader, writer)
        if hashlib.sha256(candidate.read_bytes()).hexdigest() != expected:
            raise ValueError("Native RNA distribution SHA256 mismatch")
        candidate.chmod(0o755)
        _check_binary(candidate)
        with tempfile.TemporaryDirectory(prefix=".native-rna.", dir=prefix) as tmp:
            _smoke(candidate, Path(tmp))
        os.replace(candidate, target)
    finally:
        candidate.unlink(missing_ok=True)


def native_quantifier_main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("operation", choices=("check", "prepare"))
    parser.add_argument("--prefix", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.operation == "check":
        check_quantifier(args.prefix)
    else:
        prepare_quantifier(args.prefix)
    print("Native RNA projected-gene identity contract: PASS")
    return 0
