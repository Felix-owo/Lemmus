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
stopifnot(packageVersion("DropletUtils") == "1.30.0", requireNamespace("DropletUtils", quietly=TRUE))
cat(capture.output(sessionInfo()), sep="\\n")
'''
    makevars = args.prefix / 'calling-Makevars'
    makevars.write_text('CC17 = $(CC)\nCC17STD = -std=gnu17\nC17FLAGS = $(CFLAGS)\nC17PICFLAGS = $(CPICFLAGS)\n')
    environment = dict(os.environ, R_MAKEVARS_USER=str(makevars), R_LIBS_USER=str(args.prefix / 'lib/R/library'))
    subprocess.run([str(args.conda), 'run', '--no-capture-output', '--prefix', str(args.prefix),
                    str(args.prefix / 'bin/Rscript'), '--vanilla', '-e', code], check=True, env=environment)
    return 0


def check_calling_main(argv=None):
    """读取 R 包 DESCRIPTION 核对版本；重复检查不启动 R 会话。"""
    from email.parser import Parser
    import re
    parser = argparse.ArgumentParser(prog='lemmus check-calling')
    parser.add_argument('--prefix', type=Path, required=True)
    args = parser.parse_args(argv)
    versions = {}
    for package in ('base', 'DropletUtils'):
        path = args.prefix / 'lib/R/library' / package / 'DESCRIPTION'
        metadata = Parser().parsestr(path.read_text(encoding='utf-8'))
        if metadata.get_all('Package') != [package] or len(metadata.get_all('Version', [])) != 1:
            raise ValueError(f'Invalid R package identity: {path}')
        version = metadata['Version'].strip()
        if not re.fullmatch(r'[0-9]+(?:[.\-][0-9]+)*', version):
            raise ValueError(f'Invalid R package version: {path}')
        parts = tuple(map(int, re.split(r'[.\-]', version)))
        while parts and parts[-1] == 0:
            parts = parts[:-1]
        versions[package] = parts
    if not ((4, 5) <= versions['base'] < (4, 6) and versions['DropletUtils'] == (1, 30)):
        raise ValueError('Calling requires R>=4.5,<4.6 and DropletUtils1.30.0')
    print('Calling package versions: PASS (R4.5 / DropletUtils1.30.0)')
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


def census_reference_main(argv=None):
    """Doctor: login-node preparation or read-only offline validation."""
    parser = argparse.ArgumentParser()
    parser.add_argument('operation', choices=['prepare', 'check'])
    args = parser.parse_args(argv)
    from .downstream.annotation_human import prepare_census_reference, check_census_reference
    model_dir = PIPELINE_ROOT / 'resources/hg38_reference/census_2025-11-08_scvi'
    if args.operation == 'prepare':
        prepare_census_reference(model_dir)
    else:
        check_census_reference(model_dir)
    print('Census offline reference: PASS (model, gene map, balanced embedding, identity)')
    return 0


def _check_environment_packages(prefix: Path, role: str, expected_subdir: str, section: dict) -> None:
    """只读核对已安装包；不导入目标环境的科学、绘图或 Notebook 模块。"""
    import json
    import re
    from conda.models.match_spec import MatchSpec
    from conda.models.version import VersionSpec
    from importlib.metadata import Distribution
    from packaging.markers import default_environment
    from packaging.requirements import Requirement

    records: dict[str, dict] = {}
    native = False
    for path in (prefix / "conda-meta").glob("*.json"):
        record = json.loads(path.read_text(encoding="utf-8"))
        name = record["name"]
        if name in records:
            raise ValueError(f"{role}: duplicate installed package {name}")
        records[name] = record
        subdir = record.get("subdir")
        if subdir and subdir != "noarch":
            if subdir != expected_subdir:
                raise ValueError(f"{role}: installed package {name} has platform {subdir}")
            native = True
    if not native:
        raise ValueError(f"{role}: no native package records")
    python_version = records.get("python", {}).get("version", "")
    python_short = ".".join(python_version.split(".")[:2])
    site = prefix / "lib" / f"python{python_short}" / "site-packages"
    marker_environment = default_environment()
    marker_environment.update(python_version=python_short, python_full_version=python_version)
    normalize = lambda value: re.sub(r"[-_.]+", "-", value).lower()
    # Python 包可能由 pip 安装或覆盖；有效版本来自 dist-info，不由残留 Conda 记录代替。
    python_distributions = {name: name for name in (
        "anndata", "ipykernel", "jupyterlab", "numba", "numpy", "packaging", "pandas", "polars",
        "pyyaml", "scikit-learn", "seaborn", "tqdm", "jsonschema", "snakemake",
        "snakemake-executor-plugin-slurm", "pip", "pysam", "cutadapt", "multiqc", "awscli",
        "leidenalg", "cykhash", "umi_tools",
    )}
    python_distributions.update({"matplotlib-base": "matplotlib", "python-igraph": "igraph"})
    def distribution_for(name):
        candidates = list(Distribution.discover(name=name, path=[str(site)]))
        if len(candidates) != 1:
            raise ValueError(f"{role}: missing or ambiguous Python metadata for {name}")
        distribution = candidates[0]
        if Path(distribution.locate_file("")).resolve() != site.resolve():
            raise ValueError(f"{role}: Python metadata escaped the managed prefix for {name}")
        if normalize(distribution.metadata["Name"] or "") != normalize(name):
            raise ValueError(f"{role}: Python distribution identity mismatch for {name}")
        return distribution
    for dependency in section["dependencies"]:
        if isinstance(dependency, str):
            spec = MatchSpec(dependency)
            record = records.get(spec.name)
            distribution_name = python_distributions.get(spec.name)
            if distribution_name is not None:
                distribution = distribution_for(distribution_name)
                constraint = spec.get_raw_value("version")
                if constraint and not VersionSpec(constraint).match(distribution.version):
                    raise ValueError(f"{role}: effective Python package does not satisfy {dependency}")
                if (spec.get_raw_value("build") is not None or spec.get_raw_value("channel") is not None):
                    if record is None or not spec.match(record):
                        raise ValueError(f"{role}: installed artifact does not satisfy {dependency}")
            elif record is None or not spec.match(record):
                raise ValueError(f"{role}: installed package does not satisfy {dependency}")
        else:
            for item in dependency["pip"]:
                requirement = Requirement(item)
                if requirement.marker and not requirement.marker.evaluate(marker_environment):
                    continue
                distribution = distribution_for(requirement.name)
                if distribution.version not in requirement.specifier:
                    raise ValueError(f"{role}: installed pip package does not satisfy {item}")


def check_environments_main(argv=None):
    """Doctor只读读取安装记录；完整模块导入留在候选环境创建阶段。"""
    parser = argparse.ArgumentParser(prog='lemmus check-environments')
    parser.add_argument('--role', action='append', choices=('control', 'rules', 'atlas'), required=True)
    args = parser.parse_args(argv)
    import platform
    machine = platform.machine().lower()
    os_name = 'osx' if sys.platform == 'darwin' else 'linux' if sys.platform.startswith('linux') else None
    if os_name is None:
        raise ValueError(f'unsupported managed environment platform: {sys.platform}')
    expected = os_name + '-' + ('64' if machine == 'x86_64' else machine)
    sections = yaml.safe_load((PIPELINE_ROOT / 'core/envs.yaml').read_text())
    for role in args.role:
        prefix = (PIPELINE_ROOT / 'conda/current' / role).resolve(strict=True)
        _check_environment_packages(prefix, role, expected, sections[role])
    print('Managed environment package versions/platform: PASS (' + ', '.join(args.role) + ')')
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
    init_project_main(["--project", str(project)])
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
