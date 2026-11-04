#!/usr/bin/env bash
set -euo pipefail

RNA_DARLIN_PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
RNA_DARLIN_CONDA_ROOT="$RNA_DARLIN_PIPELINE_ROOT/conda"
RNA_DARLIN_CONTROL_ENV="$RNA_DARLIN_CONDA_ROOT/current/control"
RNA_DARLIN_RULE_ENV="$RNA_DARLIN_CONDA_ROOT/current/rules"
RNA_DARLIN_ATLAS_ENV="$RNA_DARLIN_CONDA_ROOT/current/atlas"
RNA_DARLIN_CALLING_ENV="$RNA_DARLIN_CONDA_ROOT/calling"
RNA_DARLIN_ATLAS_MODEL="$RNA_DARLIN_PIPELINE_ROOT/resources/hg38_reference/census_2025-11-08_scvi/model.pt"
RNA_DARLIN_ATLAS_MODEL_SIZE=931490141
RNA_DARLIN_CONTROL_EXECUTABLES=(bin/python bin/snakemake)
RNA_DARLIN_RULE_EXECUTABLES=(bin/python bin/simpleaf bin/piscem bin/alevin-fry bin/umi_tools bin/fastqc bin/multiqc lib/jvm/bin/java)
RNA_DARLIN_ATLAS_EXECUTABLES=(bin/python bin/jupyter bin/aws)
RNA_DARLIN_CALLING_EXECUTABLES=(bin/Rscript bin/cutadapt)
export RNA_DARLIN_PIPELINE_ROOT PYTHONPATH="$RNA_DARLIN_PIPELINE_ROOT/core"

usage() {
  cat <<'TXT'
用法: doctor.sh [conda [--check|--update] [--conda PATH] | calling [--check] | atlas [--check] | reference <species> ... | test | --help]

无参数                         准备全部受管环境、hg38/mm10 PG+45S、白名单、Atlas，并检查协议 DAG
conda [--check]               检查 control/rules/atlas/calling 环境和固定设计白名单
conda [--update] [--conda P]  原子重建 control/rules/atlas，并准备独立 calling 环境
atlas [--check]               登录节点准备 Census 离线参考；--check 只读核验本地资源
calling [--check]             独立构建/检查 DD-MET5 calling 环境及固定设计白名单
reference <species> ...        安装或构建到固定 Pipeline resources/<species>_reference/pg45s
  --source DIR [--link|--copy]
  --fasta FILE --gtf FILE
test                           执行安装自检和可用的 Snakemake DAG 检查
TXT
}

# 构建事务只序列化会写共享状态的环境/reference 提交，不参与普通 Pipeline 运行。
rna_darlin_transaction_lock() {
  local name="$1" lock="$RNA_DARLIN_CONDA_ROOT/.locks/$1.lock"
  mkdir -p "$(dirname "$lock")"
  if ! mkdir "$lock" 2>/dev/null; then
    echo "ERROR: $name 构建事务正在进行: $lock" >&2
    return 75
  fi
  printf '%s\n' "$$" > "$lock/pid"
  RNA_DARLIN_ACTIVE_LOCK="$lock"
  RNA_DARLIN_TRANSACTION_KIND="$name"
  trap 'if [[ -n "${RNA_DARLIN_ACTIVE_CANDIDATE:-}" ]]; then if [[ "${RNA_DARLIN_TRANSACTION_KIND:-}" == reference ]]; then echo "保留失败 reference 候选供诊断/恢复: $RNA_DARLIN_ACTIVE_CANDIDATE" >&2; else rm -rf -- "$RNA_DARLIN_ACTIVE_CANDIDATE"; fi; fi; [[ -n "${RNA_DARLIN_ACTIVE_LOCK:-}" ]] && rm -rf -- "$RNA_DARLIN_ACTIVE_LOCK"' EXIT
  trap 'exit 130' INT
  trap 'exit 143' TERM
}

rna_darlin_release_transaction_lock() {
  [[ -n "${RNA_DARLIN_ACTIVE_LOCK:-}" ]] || return 0
  rm -rf -- "$RNA_DARLIN_ACTIVE_LOCK"
  RNA_DARLIN_ACTIVE_LOCK=""
  trap - EXIT INT TERM
}

# required executables 是环境 readiness 的唯一快速清单；不比较脚本 identity/checksum/mtime。
rna_darlin_check_env() {
  local role exe prefix missing=0
  for role in control rules; do
    if [[ "$role" == control ]]; then prefix="$RNA_DARLIN_CONTROL_ENV"; list=("${RNA_DARLIN_CONTROL_EXECUTABLES[@]}")
    else prefix="$RNA_DARLIN_RULE_ENV"; list=("${RNA_DARLIN_RULE_EXECUTABLES[@]}"); fi
    for exe in "${list[@]}"; do [[ -x "$prefix/$exe" ]] || { echo "缺失: $prefix/$exe" >&2; missing=1; }; done
  done
  return "$missing"
}

rna_darlin_check_atlas_env() {
  local exe missing=0
  for exe in "${RNA_DARLIN_ATLAS_EXECUTABLES[@]}"; do [[ -x "$RNA_DARLIN_ATLAS_ENV/$exe" ]] || { echo "缺失: $RNA_DARLIN_ATLAS_ENV/$exe" >&2; missing=1; }; done
  (( missing == 0 )) && "$RNA_DARLIN_ATLAS_ENV/bin/python" -c 'import anndata, cellxgene_census, igraph, leidenalg, scanpy, scvi, sklearn, skimage, hdf5plugin, snapatac2; from ALLCools.clustering.doublets import MethylScrublet' >/dev/null 2>&1 || missing=1
  return "$missing"
}

rna_darlin_check_all_env() { rna_darlin_check_env && rna_darlin_check_atlas_env; }

# Droplet 调用环境独立发布；不重建已有 RNA 定量环境。
rna_darlin_check_calling_env() {
  local exe
  for exe in "${RNA_DARLIN_CALLING_EXECUTABLES[@]}"; do
    [[ -x "$RNA_DARLIN_CALLING_ENV/$exe" ]] || return 1
  done
  "$RNA_DARLIN_CALLING_ENV/bin/Rscript" --vanilla -e 'stopifnot(getRversion() >= "4.5", getRversion() < "4.6", packageVersion("DropletUtils") == "1.30.0")' >/dev/null 2>&1
}

# HPC 依赖安装在登录节点完成；已有环境可在 allocation 内直接复用。
rna_darlin_require_login_for_conda() {
  if [[ "$(uname -s)" == Linux && -n "${SLURM_JOB_ID:-}" ]]; then
    echo "ERROR: Conda 构建/发布须在登录节点执行: bash core/doctor.sh conda；完成后在现有 Jupyter Terminal 重跑 Doctor。" >&2
    return 2
  fi
}

# R 原生依赖仅在候选环境构建；成功核验后原子更新 calling 指针。
rna_darlin_setup_calling() {
  local conda_bin="${RNA_DARLIN_CONTROL_ENV}/bin/conda"
  if [[ "${1:-}" == --check ]]; then
    (($#==1)) || { usage >&2; return 2; }
    rna_darlin_check_calling_env
    "$RNA_DARLIN_RULE_ENV/bin/python" -c 'from lemmus.inputs import load_droplet_whitelist,droplet_whitelist_path; load_droplet_whitelist(droplet_whitelist_path())'
    return
  fi
  (($#==0)) || { usage >&2; return 2; }
  if rna_darlin_check_calling_env; then
    if [[ "$(uname -s)" == Linux && -n "${SLURM_JOB_ID:-}" ]]; then
      rna_darlin_setup_calling --check
    else
      "$RNA_DARLIN_RULE_ENV/bin/python" -m lemmus prepare-droplet-whitelist
    fi
    return
  fi
  rna_darlin_require_login_for_conda
  [[ -x "$conda_bin" ]] || { echo "ERROR: control Conda is required" >&2; return 127; }
  rna_darlin_transaction_lock calling
  local candidate spec
  candidate="$RNA_DARLIN_CONDA_ROOT/environments/calling-$(date -u +%Y%m%dT%H%M%SZ)-$$"
  RNA_DARLIN_ACTIVE_CANDIDATE="$candidate"
  mkdir -p "$candidate" "$RNA_DARLIN_CONDA_ROOT/package-cache"
  spec="$candidate/spec.yaml"
  rna_darlin_extract_env_role calling "$spec"
  if [[ "$(uname -s)" == Linux ]]; then
    printf '%s\n' '- bioconductor-dropletutils=1.30.0' >> "$spec"
  fi
  CONDA_ALWAYS_COPY=true CONDA_PKGS_DIRS="$RNA_DARLIN_CONDA_ROOT/package-cache" "$conda_bin" env create \
    --solver libmamba --prefix "$candidate/env" --file "$spec"
  "$RNA_DARLIN_RULE_ENV/bin/python" -m lemmus prepare-calling \
    --prefix "$candidate/env" --conda "$conda_bin"
  "$candidate/env/bin/Rscript" --vanilla -e 'stopifnot(getRversion() >= "4.5", getRversion() < "4.6", packageVersion("DropletUtils") == "1.30.0")'
  "$candidate/env/bin/cutadapt" --version
  "$RNA_DARLIN_RULE_ENV/bin/python" -m lemmus prepare-droplet-whitelist
  "$RNA_DARLIN_RULE_ENV/bin/python" -c 'import os,sys; p=sys.argv[1]; t=p+".tmp."+str(os.getpid()); os.symlink(sys.argv[2],t); os.replace(t,p)' \
    "$RNA_DARLIN_CALLING_ENV" "$candidate/env"
  RNA_DARLIN_ACTIVE_CANDIDATE=""
  rna_darlin_release_transaction_lock
}

# 从单一 envs.yaml 抽取一个角色，避免恢复 core/envs/ 多文件规格。
rna_darlin_extract_env_role() {
  local role="$1" output="$2"
  awk -v role="$role" '
    $0 == role ":" {inside=1; next}
    inside && /^[^[:space:]][^:]*:/ {exit}
    inside {sub(/^  /, ""); print}
  ' "$RNA_DARLIN_PIPELINE_ROOT/core/envs.yaml" > "$output"
  [[ -s "$output" ]] || { echo "ERROR: envs.yaml 缺少角色: $role" >&2; return 2; }
}

# rules 候选环境只验证当前生产工具链可执行；真实 mapping 由 Snakemake smoke/DAG 负责。
rna_darlin_rna_toolchain_smoke() {
  local rules_env="$1"
  local ALEVIN_FRY_HOME="$rules_env/.simpleaf_home"
  local PATH="$rules_env/bin:$rules_env/lib/jvm/bin:$PATH"
  export ALEVIN_FRY_HOME PATH
  "$rules_env/bin/simpleaf" --version >/dev/null
  "$rules_env/bin/piscem" --version >/dev/null
  "$rules_env/bin/alevin-fry" --version >/dev/null
  "$rules_env/bin/umi_tools" --version >/dev/null
  "$rules_env/bin/fastqc" --version >/dev/null
  "$rules_env/bin/multiqc" --version >/dev/null
}

# Pipeline 使用显式 custom geometry 和项目 whitelist；固定空 registry，避免 simpleaf 隐式联网。
rna_darlin_prepare_simpleaf_home() {
  local home="$1" rules_env="$2" registry temporary
  registry="$home/chemistries.json"
  temporary="$home/.chemistries.json.tmp.$$"
  mkdir -p "$home"
  printf '{}\n' > "$temporary"
  mv "$temporary" "$registry"
  ALEVIN_FRY_HOME="$home" PATH="$rules_env/bin:$PATH" "$rules_env/bin/simpleaf" set-paths \
    --piscem "$rules_env/bin/piscem" \
    --alevin-fry "$rules_env/bin/alevin-fry" >/dev/null
}

# simpleaf 0.28 的 splici index 只生成 USA 三列 t2g；生产定量需要标准两列 transcript->gene 映射。
rna_darlin_write_two_column_t2g() {
  local source="$1" output="$2"
  "$RNA_DARLIN_RULE_ENV/bin/python" - "$source" "$output" <<'PY'
import os
import sys
from pathlib import Path

source = Path(sys.argv[1])
output = Path(sys.argv[2])
if not source.is_file():
    raise SystemExit(f"ERROR: simpleaf index did not produce {source}")

mappings = {}
order = []
with source.open(encoding="utf-8") as handle:
    for line_number, raw in enumerate(handle, 1):
        fields = raw.rstrip("\n").split("\t")
        if len(fields) != 3 or not fields[0] or not fields[1] or fields[2] not in {"S", "U"}:
            raise SystemExit(f"ERROR: invalid simpleaf three-column t2g at {source}:{line_number}")
        transcript, gene, _state = fields
        previous = mappings.get(transcript)
        if previous is not None and previous != gene:
            raise SystemExit(f"ERROR: transcript maps to multiple genes in {source}: {transcript}")
        if previous is None:
            mappings[transcript] = gene
            order.append(transcript)
if not order:
    raise SystemExit(f"ERROR: empty simpleaf three-column t2g: {source}")

output.parent.mkdir(parents=True, exist_ok=True)
temporary = output.with_name(f".{output.name}.tmp.{os.getpid()}")
try:
    with temporary.open("x", encoding="utf-8", newline="\n") as handle:
        for transcript in order:
            handle.write(f"{transcript}\t{mappings[transcript]}\n")
    os.replace(temporary, output)
finally:
    temporary.unlink(missing_ok=True)
PY
}

# 环境构建采用 candidate -> smoke -> current 原子指针；失败候选不污染当前环境。
rna_darlin_setup_conda() (
  local update=0 conda_bin=""
  while (($#)); do
    case "$1" in
      --update) update=1; shift ;;
      --conda) (($#>=2)) || { echo "ERROR: --conda requires PATH" >&2; return 2; }; conda_bin="$2"; shift 2 ;;
      --conda=*) conda_bin="${1#--conda=}"; shift ;;
      *) echo "ERROR: unsupported conda option: $1" >&2; return 2 ;;
    esac
  done
  if (( update == 0 )) && rna_darlin_check_all_env >/dev/null 2>&1; then
    echo "受管环境可用；使用 --update 可原子重建。"
    return 0
  fi
  rna_darlin_require_login_for_conda
  [[ -n "$conda_bin" ]] || conda_bin="$(command -v conda 2>/dev/null || true)"
  [[ -n "$conda_bin" && -x "$conda_bin" ]] || { echo "ERROR: modern Conda not found; pass --conda PATH." >&2; return 127; }
  local version major
  version="$($conda_bin --version | awk '{print $2}')"; major="${version%%.*}"
  [[ "$major" == 26 ]] || { echo "ERROR: Conda >=26,<27 required; found $version" >&2; return 2; }

  rna_darlin_transaction_lock environments
  local generation env_root setup control rules atlas current_tmp
  generation="$(date -u +%Y%m%dT%H%M%SZ)-$$"
  env_root="$RNA_DARLIN_CONDA_ROOT/environments/$generation"
  control="$env_root/control"; rules="$env_root/rules"; atlas="$env_root/atlas"; setup="$env_root/setup"
  mkdir -p "$setup" "$RNA_DARLIN_CONDA_ROOT/package-cache"
  trap 's=$?; if ((s!=0)); then chmod -R u+w "'"$env_root"'" 2>/dev/null || true; rm -rf -- "'"$env_root"'"; fi; [[ -n "${RNA_DARLIN_ACTIVE_LOCK:-}" ]] && rm -rf -- "$RNA_DARLIN_ACTIVE_LOCK"; exit $s' EXIT
  rna_darlin_extract_env_role control "$setup/control.yaml"
  rna_darlin_extract_env_role rules "$setup/rules.yaml"
  rna_darlin_extract_env_role atlas "$setup/atlas.yaml"
  unset CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_EXE CONDA_PYTHON_EXE CONDA_SHLVL PYTHONHOME PYTHONPATH PYTHONSTARTUP PYTHONUSERBASE VIRTUAL_ENV
  export CONDA_PKGS_DIRS="$RNA_DARLIN_CONDA_ROOT/package-cache" PYTHONNOUSERSITE=1
  "$conda_bin" env create --solver libmamba --prefix "$control" --file "$setup/control.yaml"
  "$conda_bin" env create --solver libmamba --prefix "$rules" --file "$setup/rules.yaml"
  "$conda_bin" env create --solver libmamba --prefix "$atlas" --file "$setup/atlas.yaml"
  "$control/bin/python" -c 'import jsonschema, yaml'
  "$rules/bin/python" -c 'import anndata, hdf5plugin, numpy, pandas, pysam, yaml'
  "$atlas/bin/python" -c 'import anndata, cellxgene_census, igraph, leidenalg, scanpy, scvi, sklearn, skimage, hdf5plugin, snapatac2; from ALLCools.clustering.doublets import MethylScrublet'
  local exe
  for exe in "${RNA_DARLIN_RULE_EXECUTABLES[@]}"; do [[ -x "$rules/$exe" ]] || { echo "ERROR: candidate missing $exe" >&2; return 1; }; done
  PATH="$rules/bin:$rules/lib/jvm/bin:$PATH" "$rules/bin/fastqc" --version >/dev/null
  PATH="$rules/bin:$rules/lib/jvm/bin:$PATH" "$rules/bin/multiqc" --version >/dev/null
  rna_darlin_rna_toolchain_smoke "$rules"
  rm -rf -- "$setup"
  "$control/bin/python" - "$RNA_DARLIN_CONDA_ROOT/current" "environments/$generation" <<'PY'
import os, sys
from pathlib import Path
active=Path(sys.argv[1]); target=sys.argv[2]
temporary=active.with_name(f".{active.name}.activate.{os.getpid()}")
if os.path.lexists(temporary): temporary.unlink()
os.symlink(target, temporary)
os.replace(temporary, active)
PY
  trap - EXIT
  rna_darlin_release_transaction_lock
  echo "Activated environment generation: $generation"
)

# Atlas asset 不进入 Snakemake DAG；doctor 统一准备，run_pipeline 仅要求 control/rules。
rna_darlin_check_atlas_assets() {
  [[ -s "$RNA_DARLIN_ATLAS_MODEL" ]] || { echo "缺失: $RNA_DARLIN_ATLAS_MODEL" >&2; return 1; }
  [[ "$(wc -c < "$RNA_DARLIN_ATLAS_MODEL" | tr -d ' ')" == "$RNA_DARLIN_ATLAS_MODEL_SIZE" ]] || { echo "ERROR: Census scVI model size mismatch" >&2; return 1; }
  "$RNA_DARLIN_ATLAS_ENV/bin/python" -m lemmus census-reference check
}

rna_darlin_setup_ontology() {
  "$RNA_DARLIN_ATLAS_ENV/bin/python" - "${1:-prepare}" "$RNA_DARLIN_PIPELINE_ROOT/resources/hg38_reference/cell_ontology/cl-basic.obo" <<'PYONTOLOGY'
import os, sys, tempfile
from pathlib import Path
from urllib.request import urlretrieve
from lemmus.downstream.annotation_human import CELL_ONTOLOGY_URL, _parse_obo
mode, target = sys.argv[1], Path(sys.argv[2])
def validate(path):
    parents, names = _parse_obo(path)
    if names.get('CL:0000000') != 'cell' or not parents:
        raise ValueError(f'Invalid Cell Ontology: {path}')
if target.is_file():
    validate(target)
elif mode == 'check':
    raise FileNotFoundError(f'Missing Cell Ontology: {target}')
else:
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.ontology.', dir=target.parent) as work:
        candidate=Path(work)/target.name
        urlretrieve(CELL_ONTOLOGY_URL, candidate)
        validate(candidate)
        os.replace(candidate, target)
print('Cell Ontology: PASS')
PYONTOLOGY
}

rna_darlin_setup_atlas() {
  local check=0
  [[ "${1:-}" == --check ]] && { check=1; shift; }
  (($#==0)) || { usage >&2; return 2; }
  rna_darlin_check_atlas_env || { echo "ERROR: atlas environment incomplete; run doctor.sh conda" >&2; return 127; }
  if (( check )); then rna_darlin_setup_ontology check; rna_darlin_check_atlas_assets; return; fi
  rna_darlin_prepare_atlas_assets
}

# 完整 Doctor 已核验环境；资产准备不重复导入整套 atlas 科学依赖。
rna_darlin_prepare_atlas_assets() {
  if [[ "$(uname -s)" == Linux && -n "${SLURM_JOB_ID:-}" ]]; then
    rna_darlin_setup_ontology check
    rna_darlin_check_atlas_assets
    return
  fi
  rna_darlin_transaction_lock atlas
  rna_darlin_setup_ontology
  if [[ ! -s "$RNA_DARLIN_ATLAS_MODEL" || "$(wc -c < "$RNA_DARLIN_ATLAS_MODEL" | tr -d ' ')" != "$RNA_DARLIN_ATLAS_MODEL_SIZE" ]]; then
    mkdir -p "$(dirname "$RNA_DARLIN_ATLAS_MODEL")"
    local tmp
    tmp="$(mktemp "$RNA_DARLIN_ATLAS_MODEL.part.XXXXXXXX")"
    RNA_DARLIN_ACTIVE_CANDIDATE="$tmp"
    PATH="$RNA_DARLIN_ATLAS_ENV/bin:$PATH" "$RNA_DARLIN_ATLAS_ENV/bin/aws" s3 cp --no-sign-request --no-progress --only-show-errors \
      s3://cellxgene-contrib-public/models/scvi/2025-11-08/homo_sapiens/model.pt "$tmp"
    [[ "$(wc -c < "$tmp" | tr -d ' ')" == "$RNA_DARLIN_ATLAS_MODEL_SIZE" ]] || { rm -f -- "$tmp"; echo "ERROR: downloaded model size mismatch" >&2; return 1; }
    rm -f -- "$(dirname "$RNA_DARLIN_ATLAS_MODEL")/census_reference.json"
    mv "$tmp" "$RNA_DARLIN_ATLAS_MODEL"
    RNA_DARLIN_ACTIVE_CANDIDATE=""
  fi
  "$RNA_DARLIN_ATLAS_ENV/bin/python" -m lemmus census-reference prepare
  rna_darlin_check_atlas_assets
  rna_darlin_release_transaction_lock
}

# 固定视图只提交已验证候选；旧目录保留，提交失败立即恢复，不递归删除共享参考。
rna_darlin_publish_reference() {
  local candidate="$1" target="$2" previous="" status
  # mktemp 候选在构建期间私有；发布时恢复当前 umask 的目录权限。
  if [[ ! -L "$candidate" ]]; then chmod "$(umask -S)" "$candidate"; fi
  if [[ -e "$target" || -L "$target" ]]; then
    previous="$(mktemp -d "${target}.previous.XXXXXXXX")"
    rmdir "$previous"
    mv -- "$target" "$previous"
  fi
  if mv -- "$candidate" "$target"; then
    [[ -z "$previous" ]] || echo "保留原 reference: $previous"
  else
    status=$?
    if [[ -n "$previous" ]]; then
      mv -- "$previous" "$target" || echo "ERROR: reference 回滚失败；保留于 $previous" >&2
    fi
    return "$status"
  fi
}

# RNA reference 使用 simpleaf 构建 splici + piscem index；一次构建，多项目共享。
rna_darlin_setup_reference() {
  (($#>=1)) || { usage >&2; return 2; }
  local species="$1"; shift
  case "$species" in hg38|mm10) ;; *) echo "ERROR: supported reference species: hg38, mm10" >&2; return 2 ;; esac
  local source_dir="" fasta="" gtf="" mode=copy
  while (($#)); do
    case "$1" in
      --source) (($#>=2)) && [[ -n "$2" ]] || { usage >&2; return 2; }; source_dir="$2"; shift 2 ;;
      --fasta) (($#>=2)) && [[ -n "$2" ]] || { usage >&2; return 2; }; fasta="$2"; shift 2 ;;
      --gtf) (($#>=2)) && [[ -n "$2" ]] || { usage >&2; return 2; }; gtf="$2"; shift 2 ;;
      --link) mode=link; shift ;;
      --copy) mode=copy; shift ;;
      *) echo "ERROR: unsupported reference option: $1" >&2; return 2 ;;
    esac
  done
  [[ -z "$source_dir" || ( -z "$fasta" && -z "$gtf" ) ]] || { echo "ERROR: --source cannot be combined with --fasta/--gtf" >&2; return 2; }
  [[ "$mode" != link || -n "$source_dir" ]] || { echo "ERROR: --link requires --source" >&2; return 2; }
  local resource_root target candidate required
  resource_root="$RNA_DARLIN_PIPELINE_ROOT/resources"
  mkdir -p "$resource_root/${species}_reference"
  # 新科学参考独占视图；父目录中的既有参考与Census/ontology不参与替换。
  target="$resource_root/${species}_reference/pg45s"
  if [[ -e "$target" ]]; then
    local complete=1
    for required in genome.fa genes/genes.gtf simpleaf/ref/t2g.tsv; do [[ -s "$target/$required" ]] || complete=0; done
    compgen -G "$target/simpleaf/index/piscem_idx*" >/dev/null || complete=0
    if (( complete )) && "$RNA_DARLIN_RULE_ENV/bin/python" -m lemmus rna-reference check --species "$species" --path "$target" >/dev/null 2>&1; then
      echo "RNA PG+45S reference already ready: $target"; return 0
    fi
    [[ ! -L "$target" && -d "$target" ]] || { echo "ERROR: incomplete reference target is not a directory: $target" >&2; return 2; }
  fi
  if [[ -z "$source_dir" && -z "$fasta" && -z "$gtf" ]]; then
    fasta="$RNA_DARLIN_PIPELINE_ROOT/resources/reference_source/${species}.fa"
    gtf="$RNA_DARLIN_PIPELINE_ROOT/resources/reference_source/${species}.genes.gtf"
  fi
  rna_darlin_transaction_lock reference
  if [[ -n "$source_dir" ]]; then
    source_dir="$(cd "$source_dir" && pwd -P)"
    for required in genome.fa genes/genes.gtf simpleaf/ref/t2g.tsv; do [[ -s "$source_dir/$required" ]] || { echo "ERROR: source reference lacks $required" >&2; return 1; }; done
    compgen -G "$source_dir/simpleaf/index/piscem_idx*" >/dev/null || { echo "ERROR: source reference lacks piscem index" >&2; return 1; }
    "$RNA_DARLIN_RULE_ENV/bin/python" -m lemmus rna-reference check --species "$species" --path "$source_dir"
    candidate="$(mktemp -d "$resource_root/${species}_reference.building.XXXXXXXX")"
    RNA_DARLIN_ACTIVE_CANDIDATE="$candidate"
    if [[ "$mode" == link ]]; then rmdir "$candidate"; ln -s "$source_dir" "$candidate"
    else cp -R "$source_dir/." "$candidate/"; fi
    "$RNA_DARLIN_RULE_ENV/bin/python" -m lemmus rna-reference check --species "$species" --path "$candidate"
    rna_darlin_publish_reference "$candidate" "$target"
  else
    [[ -n "$fasta" && -n "$gtf" && -s "$fasta" && -s "$gtf" ]] || { echo "ERROR: --fasta and --gtf are required and must be non-empty" >&2; return 2; }
    [[ -x "$RNA_DARLIN_RULE_ENV/bin/simpleaf" ]] || { echo "ERROR: rules environment missing; run doctor.sh conda" >&2; return 127; }
    candidate="$(mktemp -d "$resource_root/${species}_reference.building.XXXXXXXX")"
    RNA_DARLIN_ACTIVE_CANDIDATE="$candidate"
    mkdir -p "$candidate/genes"
    cp "$fasta" "$candidate/genome.fa"; cp "$gtf" "$candidate/genes/genes.gtf"
    "$RNA_DARLIN_RULE_ENV/bin/python" -m lemmus rna-reference prepare \
      --species "$species" --path "$candidate/rrna_competitor.fa"
    export ALEVIN_FRY_HOME="$RNA_DARLIN_CONDA_ROOT/simpleaf_home"
    rna_darlin_prepare_simpleaf_home "$ALEVIN_FRY_HOME" "$RNA_DARLIN_RULE_ENV"
    mkdir -p "$candidate/.simpleaf_work" "$candidate/.simpleaf_tmp"
    ulimit -n 2048 || true
    PATH="$RNA_DARLIN_RULE_ENV/bin:$PATH" "$RNA_DARLIN_RULE_ENV/bin/simpleaf" index \
      --output "$candidate/simpleaf" --fasta "$candidate/genome.fa" --gtf "$candidate/genes/genes.gtf" \
      --unspliced "$candidate/rrna_competitor.fa" \
      --threads 4 --work-dir "$candidate/.simpleaf_work" --tmp-dir "$candidate/.simpleaf_tmp"
    rna_darlin_write_two_column_t2g "$candidate/simpleaf/ref/t2g_3col.tsv" "$candidate/simpleaf/ref/t2g.tsv"
    [[ -s "$candidate/simpleaf/ref/t2g.tsv" ]] || { echo "ERROR: simpleaf reference build incomplete: t2g.tsv" >&2; return 1; }
    compgen -G "$candidate/simpleaf/index/piscem_idx*" >/dev/null || { echo "ERROR: simpleaf reference build incomplete: piscem index" >&2; return 1; }
    "$RNA_DARLIN_RULE_ENV/bin/python" -m lemmus rna-reference seal --species "$species" --path "$candidate"
    rm -rf "$candidate/.simpleaf_work" "$candidate/.simpleaf_tmp"
    rna_darlin_publish_reference "$candidate" "$target"
  fi
  RNA_DARLIN_ACTIVE_CANDIDATE=""; rna_darlin_release_transaction_lock
  echo "RNA simpleaf/piscem reference ready: $target"
}

# test 只验证当前生产 contract；不恢复已删除 helper、compatibility 或 identity gate。
rna_darlin_test() (
  local python="${PYTHON:-$(command -v python3)}" unit_python validated=0 tmp required="${1:-0}"
  unit_python="$python"
  if rna_darlin_check_env >/dev/null 2>&1; then
    python="$RNA_DARLIN_CONTROL_ENV/bin/python"
    unit_python="$RNA_DARLIN_RULE_ENV/bin/python"
    validated=1
  fi
  tmp="$(mktemp -d "${TMPDIR:-/tmp}/lemmus-selfcheck.XXXXXXXX")"; trap 'rm -rf -- "'"$tmp"'"' EXIT
  # 只检查语法，不为整个依赖树重建pycache；共享文件系统上尤其避免数百次小文件提交。
  unset PYTHONPYCACHEPREFIX
  export PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR="$tmp/matplotlib" XDG_CACHE_HOME="$tmp/cache" PYTHONPATH="$RNA_DARLIN_PIPELINE_ROOT/core"
  echo "[1/6] Python compile (in memory)"; "$python" - "$RNA_DARLIN_PIPELINE_ROOT/core/lemmus" <<'PY'
import sys
from pathlib import Path
package = Path(sys.argv[1])
for source in [*package.glob('*.py'), *(package / 'downstream').glob('*.py')]:
    compile(source.read_bytes(), str(source), 'exec')
PY
  echo "[2/6] Shell syntax：检查入口脚本语法"; bash -n "$RNA_DARLIN_PIPELINE_ROOT/core/doctor.sh" "$RNA_DARLIN_PIPELINE_ROOT/core/run_pipeline.sh"
  echo "[3/6] envs/config/schema"; "$python" - "$RNA_DARLIN_PIPELINE_ROOT" <<'PY'
import jsonschema, yaml, sys
from pathlib import Path
root=Path(sys.argv[1]); from lemmus import config as p
envs=yaml.safe_load((root/'core/envs.yaml').read_text()); assert set(envs)=={'control','rules','atlas','calling'}
schema=yaml.safe_load(p.CONFIG_SCHEMA_YAML)
jsonschema.Draft202012Validator.check_schema(schema)
for protocol, template in p.CONFIG_TEMPLATES.items():
    cfg=yaml.safe_load(template)
    p.validate_config(cfg); jsonschema.Draft202012Validator(schema).validate(p.select_config(cfg))
    assert cfg['pipeline_mode'] == protocol
    print(f'{protocol} config: PASS')
print('envs/config/schema: PASS')
PY
  echo "[4/6] Runtime behavior：合成输入自检，不代表真实数据验收"; "$unit_python" -m lemmus self-check
  echo "[5/6] stale dependency boundary：核对当前运行入口"
  local stale=0
  if grep -RInE 'core/(scripts|internal|profiles|envs/)|(^|[^[:alnum:]_/])src/|\.\.?/src/|pyproject\.toml|config\.template\.yaml|config\.schema\.yaml' "$RNA_DARLIN_PIPELINE_ROOT/core/Snakefile" "$RNA_DARLIN_PIPELINE_ROOT"/core/lemmus/*.py "$RNA_DARLIN_PIPELINE_ROOT/README.md"; then stale=1; fi
  if grep -RInE '(^|[^[:alnum:]_])\.runtime([^[:alnum:]_]|$)' "$RNA_DARLIN_PIPELINE_ROOT/core/Snakefile" "$RNA_DARLIN_PIPELINE_ROOT"/core/lemmus/*.py "$RNA_DARLIN_PIPELINE_ROOT/README.md"; then stale=1; fi
  (( stale == 0 )) || { echo "ERROR: stale production reference found" >&2; return 1; }
  if (( validated )); then
    # 每个 case 建独立 Workflow；同一解释器复用 Snakemake/NumPy/Pandas 导入。
    "$python" - "$RNA_DARLIN_PIPELINE_ROOT" "$tmp" "$required" "$RNA_DARLIN_RULE_ENV" <<'PYDAG'
import contextlib, io, os, sys
from pathlib import Path
import yaml
from lemmus.config import load_config
from lemmus.maintenance import parser_fixture_main
from lemmus.rna import validate_rna_reference
from snakemake.cli import main
root, tmp, required, rule_env = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
original = Path.cwd()
diagnostics=io.StringIO()
with open(os.devnull, 'w') as sink:
    for species in ('hg38', 'mm10'):
        for protocol in ('cabernet', 'droplet'):
            project=tmp/f'{species}-{protocol}'
            with contextlib.redirect_stdout(sink):parser_fixture_main(['--project', str(project)])
            path=project/'00_config/config.yaml';config=yaml.safe_load(path.read_text())
            config['pipeline_mode']=protocol;config[protocol]['RNA']['reference']=species
            path.write_text(yaml.safe_dump(config, sort_keys=False))
            # 使用同一生产配置解析器；schema/输入错误不能冒充 reference 缺失。
            active=load_config(path)
            try:validate_rna_reference(root/f"resources/{active['RNA']['reference']}_reference/pg45s", species)
            except (OSError, ValueError):
                if required:raise
                print(f'Snakemake full-DAG dry-run ({species} / {protocol}): SKIPPED (shared Pipeline reference unavailable)')
                continue
            os.environ.update(RNA_DARLIN_CONFIG=str(path), RNA_DARLIN_PROJECT_DIR=str(project), RNA_DARLIN_RULE_ENV=rule_env, RNA_DARLIN_PROTOCOL=protocol)
            try:
                os.chdir(project)
                diagnostic_start=diagnostics.tell()
                with contextlib.redirect_stdout(diagnostics):
                    try:main(['--snakefile',str(root/'core/Snakefile'),'--executor','local','--cores','1','--dry-run','--runtime-source-cache-path',str(tmp/'source-cache'),'--configfile',str(path)])
                    except SystemExit as exc:
                        if exc.code not in (None,0):
                            print(diagnostics.getvalue()[diagnostic_start:], file=sys.stderr, end='')
                            raise
            finally:os.chdir(original)
            print(f'Snakemake full-DAG dry-run ({species} / {protocol}): PASS')
PYDAG
  else
    if (( required )); then echo "ERROR: full-DAG requires managed environments" >&2; return 1; fi
    echo "Snakemake full-DAG dry-run: SKIPPED (managed environment unavailable)"
  fi
  echo "[6/6] Alopex cross-pipeline alignment (cabernet / droplet)"
  rna_darlin_check_pipeline_alignment
  rm -rf -- "$tmp"; trap - EXIT
)

# 项目只选择 species；reference 始终来自 Pipeline resources/<species>_reference/pg45s。
rna_darlin_check_project_reference() {
  local config="$1" python="$RNA_DARLIN_CONTROL_ENV/bin/python"
  [[ -x "$python" ]] || python="${PYTHON:-$(command -v python3)}"
  "$python" - "$config" <<'PY'
import sys, yaml
from pathlib import Path
from lemmus.config import load_config; from lemmus.rna import validate_rna_reference
p=Path(sys.argv[1]); c=load_config(p); species=str(c['RNA']['reference'])
r=Path(__import__('os').environ['RNA_DARLIN_PIPELINE_ROOT'])/'resources'/f'{species}_reference'/'pg45s'
validate_rna_reference(r, species)
print('reference check: PASS')
PY
}

# 只读核对 Alopex 的孔板映射与固定 Droplet 设计，不调用 DNA doctor 或修改资源。
rna_darlin_check_pipeline_alignment() {
  local python="$RNA_DARLIN_RULE_ENV/bin/python"
  [[ -x "$python" ]] || python="${PYTHON:-$(command -v python3)}"
  "$python" -m lemmus cross-pipeline-check
}

# 子进程保持独立 errexit；网络失败不妨碍其它物种/reference 与安装自检完成。
rna_darlin_all() {
  local status=0 species started=$SECONDS
  if [[ "$(uname -s)" == Linux && -n "${SLURM_JOB_ID:-}" ]]; then
    bash "$RNA_DARLIN_PIPELINE_ROOT/core/doctor.sh" conda --check
  else
    bash "$RNA_DARLIN_PIPELINE_ROOT/core/doctor.sh" conda
  fi
  echo "Doctor phase environments: $((SECONDS-started))s"
  for species in hg38 mm10; do
    started=$SECONDS
    if bash "$RNA_DARLIN_PIPELINE_ROOT/core/doctor.sh" reference "$species"; then :
    else echo "ERROR: $species reference 准备失败" >&2; status=1; fi
    echo "Doctor phase reference $species: $((SECONDS-started))s"
  done
  started=$SECONDS
  if bash -c 'source "$1"; rna_darlin_test 1' -- "$RNA_DARLIN_PIPELINE_ROOT/core/doctor.sh"; then :
  else echo "ERROR: 安装自检失败" >&2; status=1; fi
  echo "Doctor phase test: $((SECONDS-started))s"
  if [[ -f "$PWD/00_config/config.yaml" ]]; then
    # 项目配置只读；失败必须阻止 READY。
    if (set -e; rna_darlin_check_project_reference "$PWD/00_config/config.yaml"); then :
    else echo "ERROR: 当前项目 reference 检查失败" >&2; status=1; fi
  fi
  started=$SECONDS
  if bash -c 'source "$1"; rna_darlin_prepare_atlas_assets' -- "$RNA_DARLIN_PIPELINE_ROOT/core/doctor.sh"; then :
  else echo "ERROR: Atlas 离线参考未就绪；请在登录节点运行 bash core/doctor.sh atlas" >&2; status=1; fi
  echo "Doctor phase atlas offline reference: $((SECONDS-started))s"
  (( status == 0 )) || return "$status"
  echo "Lemmus Doctor: READY (control/rules/atlas/calling, hg38/mm10 PG+45S, ME5, protocol DAGs)"
}

# 被 source 时只提供共享清单与函数；不解析参数、不执行维护动作。
if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
  return 0
fi

case "${1:-all}" in
  -h|--help|help) usage ;;
  conda)
    shift
    if [[ "${1:-}" == --check ]]; then shift; (($#==0)) || { usage >&2; exit 2; }; rna_darlin_check_all_env; rna_darlin_setup_calling --check
    else rna_darlin_setup_conda "$@"; rna_darlin_setup_calling; fi ;;
  atlas) shift; rna_darlin_setup_atlas "$@" ;;
  calling) shift; rna_darlin_setup_calling "$@" ;;
  reference) shift; rna_darlin_setup_reference "$@" ;;
  test) shift; (($#==0)) || { usage >&2; exit 2; }; rna_darlin_test ;;
  all)
    (($#==0)) || { usage >&2; exit 2; }
    rna_darlin_all ;;
  *) usage >&2; exit 2 ;;
esac
