"""Pipeline 安装根定位；不依赖当前工作目录。"""
from pathlib import Path

PIPELINE_ROOT = Path(__file__).resolve().parents[2]
