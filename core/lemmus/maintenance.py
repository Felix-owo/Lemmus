"""Lemmus：环境规格、参考与运行自检。"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import subprocess
import sys
import tempfile
from pathlib import Path
import yaml
from .paths import PIPELINE_ROOT
from .config import (DROPLET_WHITELIST_URL)
from .inputs import (
    check_cross_pipeline_alignment,
    droplet_whitelist_path,
    init_project_main,
    load_droplet_whitelist,
    manifest_refresh_main,
)


def cross_pipeline_check_main(argv=None):
    """Doctor 的跨 Pipeline 检查入口；默认定位同级 Alopex，不按当前工作目录猜测。"""
    argparse.ArgumentParser(prog='lemmus cross-pipeline-check').parse_args(argv)
    override = os.environ.get('DNA_PIPELINE_ROOT')
    try:
        if override is not None and not override.strip():
            raise ValueError('DNA_PIPELINE_ROOT 不能显式为空')
        check_cross_pipeline_alignment(PIPELINE_ROOT,
                                       Path(override) if override is not None else None)
    except (OSError, ValueError, csv.Error) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 1
    return 0


def prepare_droplet_whitelist_main(argv=None):
    """下载并校验固定版本的设计白名单，仅供 Doctor 构建调用。"""
    import urllib.request
    parser = argparse.ArgumentParser()
    parser.parse_args(argv)
    path = droplet_whitelist_path()
    if path.exists():
        load_droplet_whitelist(path)
        return 0
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.whitelist.', dir=path.parent) as tmp:
        candidate = Path(tmp) / path.name
        with urllib.request.urlopen(DROPLET_WHITELIST_URL, timeout=120) as response:
            candidate.write_bytes(response.read())
        load_droplet_whitelist(candidate)
        os.replace(candidate, path)
    return 0


def prepare_calling_main(argv=None):
    """在独立候选环境安装固定 Bioconductor 调用依赖，保留 RNA 定量环境。"""
    parser = argparse.ArgumentParser()
    parser.add_argument('--prefix', type=Path, required=True)
    parser.add_argument('--conda', type=Path, required=True)
    args = parser.parse_args(argv)
    code = r'''
# 外部 curl 逐包下载，避免 libcurl 批量连接全部耗尽超时；HTTP错误必须失败。
stopifnot(nzchar(Sys.which("curl")))
options(repos=c(CRAN="https://cloud.r-project.org"), timeout=1200, Ncpus=4,
        download.file.method="curl", download.file.extra=c("-fL", "--max-time", "1200"))
stopifnot(getRversion() >= "4.5", getRversion() < "4.6")
if (!requireNamespace("DropletUtils", quietly=TRUE)) {
    BiocManager::install(version="3.22", ask=FALSE, update=FALSE)
    BiocManager::install(c("SingleCellExperiment", "BiocParallel", "SparseArray", "DelayedArray",
        "DelayedMatrixStats", "HDF5Array", "rhdf5", "edgeR", "beachmat", "scuttle", "assorthead",
        "Rhdf5lib"), ask=FALSE, update=FALSE, type="source")
    install.packages("https://bioconductor.org/packages/3.22/bioc/src/contrib/DropletUtils_1.30.0.tar.gz",
        repos=NULL, type="source")
}
stopifnot(packageVersion("DropletUtils") == "1.30.0")
cat(capture.output(sessionInfo()), sep="\\n")
'''
    makevars = args.prefix / 'calling-Makevars'
    makevars.write_text('CC17 = $(CC)\nCC17STD = -std=gnu17\nC17FLAGS = $(CFLAGS)\nC17PICFLAGS = $(CPICFLAGS)\n')
    environment = dict(os.environ, R_MAKEVARS_USER=str(makevars), R_LIBS_USER=str(args.prefix / 'lib/R/library'))
    subprocess.run([str(args.conda), 'run', '--no-capture-output', '--prefix', str(args.prefix),
                    str(args.prefix / 'bin/Rscript'), '--vanilla', '-e', code], check=True, env=environment)
    return 0


def rna_reference_main(argv=None):
    """内部CLI把固定设计准备、构建后提交和只读检查连接到Doctor。"""
    parser = argparse.ArgumentParser()
    parser.add_argument('operation', choices=['prepare', 'seal', 'check'])
    parser.add_argument('--species', required=True)
    parser.add_argument('--path', type=Path, required=True)
    args = parser.parse_args(argv)
    from .rna import prepare_rrna_reference, seal_rna_reference, validate_rna_reference
    if args.operation == 'prepare': prepare_rrna_reference(args.species, args.path)
    elif args.operation == 'seal': seal_rna_reference(args.path, args.species)
    else: validate_rna_reference(args.path, args.species)
    return 0


def env_spec_main(argv=None):
    """从单一 envs.yaml 提取指定角色为 Conda 可消费 spec。"""
    parser = argparse.ArgumentParser(prog="lemmus env-spec")
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--role", choices=("control", "rules"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    import yaml
    payload = yaml.safe_load(args.source.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or args.role not in payload or not isinstance(payload[args.role], dict):
        raise SystemExit(f"invalid environment role: {args.role}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(yaml.safe_dump(payload[args.role], sort_keys=False), encoding="utf-8")
    return 0


def parser_fixture_main(argv=None):
    """生成只用于 Snakemake parser/full-DAG smoke 的最小项目 fixture。"""
    parser = argparse.ArgumentParser(prog="lemmus parser-fixture")
    parser.add_argument("--project", type=Path, required=True)
    args = parser.parse_args(argv)
    project = args.project.resolve()
    init_project_main(["--project", str(project), "--protocol", "cabernet"])
    # 只准备parser输入，不伪造可用reference或声称比对真值。
    import gzip
    raw_rna = project / "01_raw/RNA"
    raw_darlin = project / "01_raw/DARLIN"
    seq = "ACGT" * 30
    qual = "I" * len(seq)
    for base, sample in ((raw_rna, "SMOKE"), (raw_darlin, "SMOKE")):
        for mate in (1, 2):
            with gzip.open(base / f"{sample}_R{mate}.fastq.gz", "wt", encoding="ascii") as handle:
                handle.write(f"@smoke/"+str(mate)+f"\n{seq}\n+\n{qual}\n")
    manifest_refresh_main(["--project", str(project)])
    return 0


def self_check_main(argv=None):
    """发行包自检：执行原生barcode纠正与H5AD序列化，不依赖开发测试树。"""
    argparse.ArgumentParser(prog="lemmus self-check").parse_args(argv)
    import anndata as ad
    import numpy as np
    import pandas as pd
    from .darlin import correct_darlin_cell_barcodes
    from .downstream.io import atomic_write_h5ad, write_lock

    raw = pd.DataFrame({
        "cell_bc": ["AAAAAAAA", "TAAAAAAA", "AAAAAAAC", "CCCCCCCC"],
        "library": ["001"] * 4, "umi": ["ACGTACGT"] * 4,
        "cell_id": ["", "", "", ""], "umi_id": ["", "", "", ""],
    })
    corrected = correct_darlin_cell_barcodes(raw, ["AAAAAAAA", "AAAAAACC"])
    expected = ["exact", "corrected_hamming1", "ambiguous", "unmatched"]
    if (corrected["cell_bc_correction_status"].tolist() != expected
            or corrected["cell_bc"].tolist() != ["AAAAAAAA", "AAAAAAAA", "", ""]
            or corrected["cell_id"].tolist()[:2] != ["001_AAAAAAAA"] * 2):
        raise RuntimeError("DARLIN barcode self-check failed")
    counts = np.array([[0, 3], [2, 0]], dtype=np.int64)
    cells, genes = ["001_NA", "001_01"], ["G1", "mt-Co1"]
    with tempfile.TemporaryDirectory(prefix="lemmus-selfcheck.") as tmp:
        path = Path(tmp) / "RNA.h5ad"
        value = ad.AnnData(counts, obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=genes))
        with write_lock(path.parent):
            atomic_write_h5ad(value, path)
        restored = ad.read_h5ad(path)
        if (restored.obs_names.tolist() != cells or restored.var_names.tolist() != genes
                or restored.X.dtype != counts.dtype or not np.array_equal(restored.X, counts)):
            raise RuntimeError("RNA H5AD self-check failed")
    print("Runtime barcode/H5AD self-check: PASS (synthetic, not scientific acceptance)")
    return 0
