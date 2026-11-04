#!/usr/bin/env bash
# Parse the complete launcher before starting a long-running controller.
# Bash must not resume reading this file at an obsolete offset after an edit.
rna_darlin_main() {
set -euo pipefail

# RNA_DARLIN_Pipeline 唯一运行入口：先解析项目软连接；macOS -> local；Linux/HPC login/compute -> Slurm executor。
rna_darlin_resolve_launcher_source() {
  local source="${BASH_SOURCE[0]}" directory target
  while [[ -L "$source" ]]; do
    directory="$(cd -P "$(dirname "$source")" && pwd)"
    target="$(readlink "$source")"
    if [[ "$target" == /* ]]; then
      source="$target"
    else
      source="$directory/$target"
    fi
  done
  directory="$(cd -P "$(dirname "$source")" && pwd)"
  printf '%s/%s\n' "$directory" "$(basename "$source")"
}
RNA_DARLIN_LAUNCHER_SOURCE="$(rna_darlin_resolve_launcher_source)"
RNA_DARLIN_PIPELINE_ROOT="$(cd "$(dirname "$RNA_DARLIN_LAUNCHER_SOURCE")/.." && pwd -P)"
export RNA_DARLIN_PIPELINE_ROOT PYTHONPATH="$RNA_DARLIN_PIPELINE_ROOT/core"
source "$RNA_DARLIN_PIPELINE_ROOT/core/doctor.sh"

usage() {
  cat <<'TXT'
用法: run_pipeline.sh [--init-project --protocol NAME | --refresh-manifest | --dry-run|-n | --help]

无参数              正式运行当前项目
--init-project       在当前目录初始化项目后立即退出
--protocol NAME      初始化时必填：cabernet 或 droplet；同一 config 展示两个模块，选择其中一个生效
--refresh-manifest   刷新 sample manifest 后立即退出
--dry-run, -n        只解析完整 DAG，不提交 rule job

平台：macOS 使用 local；Linux/HPC login 或已有 compute allocation 都由当前 shell
作为 Snakemake controller，并使用 Slurm executor 提交 rule jobs。

Cabernet 初始化: run_pipeline.sh --init-project --protocol cabernet
DD-MET5 初始化:  run_pipeline.sh --init-project --protocol droplet
TXT
}

MODE=run
INIT_PROTOCOL=""
while (($#)); do
  case "$1" in
    --init-project) [[ "$MODE" == run ]] || { echo "ERROR: maintenance flags cannot be combined" >&2; exit 2; }; MODE=init ;;
    --protocol)
      [[ $# -ge 2 && -z "$INIT_PROTOCOL" ]] || { echo "ERROR: --protocol requires one value: cabernet or droplet" >&2; exit 2; }
      case "$2" in cabernet|droplet) INIT_PROTOCOL="$2" ;; *) echo "ERROR: --protocol must be cabernet or droplet" >&2; exit 2 ;; esac
      shift ;;
    --refresh-manifest) [[ "$MODE" == run ]] || { echo "ERROR: maintenance flags cannot be combined" >&2; exit 2; }; MODE=refresh ;;
    --dry-run|-n) [[ "$MODE" == run ]] || { echo "ERROR: maintenance flags cannot be combined" >&2; exit 2; }; MODE=dryrun ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unsupported launcher argument/target: $1" >&2; usage >&2; exit 2 ;;
  esac
  shift
done

if [[ "$MODE" == init && -z "$INIT_PROTOCOL" ]]; then
  echo "ERROR: initialization requires --protocol cabernet or --protocol droplet" >&2
  exit 2
fi
if [[ "$MODE" != init && -n "$INIT_PROTOCOL" ]]; then
  echo "ERROR: --protocol is only valid with --init-project; existing projects use config.yaml" >&2
  exit 2
fi

rna_darlin_find_project() {
  local dir="$PWD" parent
  while :; do
    [[ -f "$dir/00_config/config.yaml" ]] && { printf '%s\n' "$dir"; return 0; }
    parent="$(dirname "$dir")"; [[ "$parent" != "$dir" ]] || break; dir="$parent"
  done
  return 2
}

if [[ "$MODE" == init ]]; then
  init_python="$RNA_DARLIN_CONTROL_ENV/bin/python"
  [[ -x "$init_python" ]] || init_python="$(command -v python3 || command -v python)"
  PYTHONPATH="$RNA_DARLIN_PIPELINE_ROOT/core" "$init_python" -m lemmus init-project --project "$PWD" --protocol "$INIT_PROTOCOL"
  exit $?
fi

PROJECT="$(rna_darlin_find_project)" || { echo "ERROR: cannot find <project>/00_config/config.yaml from current directory upward" >&2; exit 2; }
CONFIG="$PROJECT/00_config/config.yaml"

if [[ "$MODE" == refresh ]]; then
  refresh_python="$RNA_DARLIN_CONTROL_ENV/bin/python"
  [[ -x "$refresh_python" ]] || refresh_python="$(command -v python3 || command -v python)"
  PYTHONPATH="$RNA_DARLIN_PIPELINE_ROOT/core" "$refresh_python" -m lemmus manifest-refresh --project "$PROJECT"
  echo "manifest 已刷新: $PROJECT/00_config/sample_manifest.tsv"
  exit 0
fi

rna_darlin_check_env || { echo "ERROR: managed environment incomplete; run core/doctor.sh conda" >&2; exit 127; }
CONTROL_ENV="$RNA_DARLIN_CONTROL_ENV"; RULE_ENV="$RNA_DARLIN_RULE_ENV"
PYTHON="$CONTROL_ENV/bin/python"; SNAKEMAKE="$CONTROL_ENV/bin/snakemake"
# Check the exporter before eval: eval of empty stdout would hide its failure.
CONFIG_SCALARS="$("$PYTHON" -m lemmus config-scalars --config "$CONFIG")" || exit $?
eval "$CONFIG_SCALARS"
# Keep this array nonempty: macOS Bash 3.2 treats empty arrays as unset under -u.
LOCAL_ARGS=(--executor local --cores "$RNA_CORES")
if [[ "$PIPELINE_MODE" == droplet ]]; then
  rna_darlin_check_calling_env || { echo "ERROR: DD-MET5 requires core/doctor.sh calling" >&2; exit 127; }
  LOCAL_ARGS+=(--resources "mem_mb=$RNA_MEM_MB")
fi

export RNA_DARLIN_CONFIG="$CONFIG" RNA_DARLIN_PROJECT_DIR="$PROJECT" RNA_DARLIN_RULE_ENV="$RULE_ENV"
export PATH="$RULE_ENV/bin:$RULE_ENV/lib/jvm/bin:$CONTROL_ENV/bin:$PATH" PYTHONNOUSERSITE=1
export XDG_CACHE_HOME="$PROJECT/05_tmp/xdg_cache" MPLCONFIGDIR="$PROJECT/05_tmp/matplotlib"
unset PYTHONHOME PYTHONSTARTUP PYTHONUSERBASE VIRTUAL_ENV
mkdir -p "$XDG_CACHE_HOME" "$MPLCONFIGDIR" "$PROJECT/04_logs/controller" "$PROJECT/04_logs/slurm"

# Configure simpleaf once per project/controller run. Worker rules only read this cache,
# avoiding concurrent `simpleaf set-paths` writes from multiple samples.
export ALEVIN_FRY_HOME="$PROJECT/05_tmp/simpleaf_home"
rna_darlin_prepare_simpleaf_home "$ALEVIN_FRY_HOME" "$RULE_ENV"

COMMON=(--keep-going --keep-incomplete --retries 1 --rerun-incomplete --rerun-triggers mtime input params code software-env --show-failed-logs --latency-wait "$LATENCY" --configfile "$CONFIG")
cd "$PROJECT"
if [[ "$MODE" == dryrun ]]; then
  export RNA_DARLIN_EXECUTOR=local
  exec "$SNAKEMAKE" --snakefile "$RNA_DARLIN_PIPELINE_ROOT/core/Snakefile" "${LOCAL_ARGS[@]}" "${COMMON[@]}" --dry-run
fi

RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)-$$"; LOG="$PROJECT/04_logs/controller/$RUN_ID.log"
echo "INFO: Lemmus 开始；协议：${PIPELINE_MODE}；项目：${PROJECT}"
echo "INFO: 配置：${CONFIG}；RNA CPU：${RNA_CORES}；RNA 内存预算：${RNA_MEM_MB} MiB"
case "$(uname -s)" in
  Darwin)
    export RNA_DARLIN_EXECUTOR=local
    echo "运行模式: macOS/local; controller log=$LOG"
    "$SNAKEMAKE" --snakefile "$RNA_DARLIN_PIPELINE_ROOT/core/Snakefile" "${LOCAL_ARGS[@]}" "${COMMON[@]}" 2>&1 | tee "$LOG"
    ;;
  Linux)
    command -v sbatch >/dev/null 2>&1 || { echo "ERROR: Linux/HPC requires sbatch for the Slurm executor" >&2; exit 127; }
    export RNA_DARLIN_EXECUTOR=slurm
    if [[ -n "${SLURM_JOB_ID:-}" ]]; then echo "运行模式: Linux/HPC compute allocation controller; log=$LOG"; else echo "运行模式: Linux/HPC login controller; log=$LOG"; fi
    "$SNAKEMAKE" --snakefile "$RNA_DARLIN_PIPELINE_ROOT/core/Snakefile" \
      --executor slurm --scheduler greedy --jobs "$JOBS" --local-cores "$LOCAL_CORES" \
      --max-jobs-per-timespan 30/1s --max-status-checks-per-second 5 \
      --default-resources "slurm_partition=$PARTITION" \
      --slurm-logdir "$PROJECT/04_logs/slurm" --slurm-delete-logfiles-older-than 0 --slurm-qos "$QOS" \
      "${COMMON[@]}" 2>&1 | tee "$LOG"
    ;;
  *) echo "ERROR: unsupported platform: $(uname -s)" >&2; exit 2 ;;
esac
}

# 调用和退出属于同一解析单元，运行期间更新源码后不再读取旧脚本尾部。
rna_darlin_main "$@"; exit "$?"
