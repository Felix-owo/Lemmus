"""Pipeline安装根与项目实际结果定位；不依赖调用者当前工作目录。"""
from pathlib import Path

PIPELINE_ROOT = Path(__file__).resolve().parents[2]


def project_results_root(project_dir: Path, config: dict) -> Path:
    """与DAG消费相同的结果根；不搜索默认目录中的旧产物。"""
    value = (config.get('paths') or {}).get('results_root','03_results')
    path = Path(value).expanduser()
    return (path if path.is_absolute() else Path(project_dir)/path).resolve()
