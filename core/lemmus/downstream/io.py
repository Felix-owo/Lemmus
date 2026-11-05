"""下游原子写出、序列化与完整Processor源码身份。"""
from __future__ import annotations

import fcntl
import stat as _qc_lock_stat
import threading
from contextlib import ExitStack, contextmanager

import hashlib
import json
import os
from pathlib import Path
import pandas as pd
from .config import DownstreamQCError
from ..io import atomic_output




def atomic_write_bytes(data: bytes, destination: Path) -> None:
    with atomic_output(destination) as temporary:
        temporary.write_bytes(data)


def atomic_write_csv(frame: pd.DataFrame, path: Path, *, index: bool = False) -> None:
    with atomic_output(path) as temporary:
        frame.to_csv(temporary, index=index)


def atomic_write_json(payload: dict, path: Path) -> None:
    with atomic_output(path) as temporary:
        temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _plain_string_frame(frame):
    import pandas as pd

    frame = frame.copy()
    frame.index = pd.Index(frame.index.astype(str).to_numpy(dtype=object), name=frame.index.name)
    frame.columns = pd.Index(frame.columns.astype(str).to_numpy(dtype=object))
    for column in frame.columns:
        if isinstance(frame[column].dtype, pd.StringDtype):
            frame[column] = frame[column].astype(object)
    return frame


def _sanitize_h5ad_value(value):
    import pandas as pd

    if isinstance(value, pd.DataFrame):
        return _plain_string_frame(value)
    if isinstance(value, pd.Series):
        result = value.copy()
        result.index = pd.Index(result.index.astype(str).to_numpy(dtype=object))
        if isinstance(result.dtype, pd.StringDtype):
            result = result.astype(object)
        return result
    if isinstance(value, dict):
        return {key: _sanitize_h5ad_value(item) for key, item in value.items()}
    if isinstance(value, list) and value and all(isinstance(x, dict) for x in value):
        return _plain_string_frame(pd.DataFrame(value))
    if isinstance(value, list):
        return [_sanitize_h5ad_value(item) for item in value]
    return value


def sanitize_for_h5ad(adata):
    """Last-line H5AD serialization guard for tables and mixed containers."""
    adata.obs = _plain_string_frame(adata.obs)
    adata.var = _plain_string_frame(adata.var)
    if adata.raw is not None:
        raw = adata.raw.to_adata()
        raw.obs = _plain_string_frame(raw.obs)
        raw.var = _plain_string_frame(raw.var)
        adata.raw = raw
    for key, value in list(adata.uns.items()):
        adata.uns[key] = _sanitize_h5ad_value(value)
    return adata


def atomic_write_h5ad(adata, destination: Path) -> None:
    """只提交本次调用的完整H5AD；失败保留已有结果并清理私有候选。"""
    with atomic_output(destination) as temporary:
        adata.write_h5ad(temporary)


def downstream_source_sha256() -> str:
    """绑定完整下游Processor及实际消费的核心配置/IO/路径源码，不执行被读取源码。"""
    root = Path(__file__).resolve().parent
    core = root.parent
    paths = [path for path in sorted(root.rglob("*.py"))
             if not any(part in {"__pycache__", ".ipynb_checkpoints"} for part in path.relative_to(root).parts)]
    paths += [core / name for name in ("config.py", "io.py", "paths.py")]
    digest = hashlib.sha256()
    for path in sorted(paths):
        relative = path.relative_to(core).as_posix().encode()
        data = path.read_bytes()
        digest.update(len(relative).to_bytes(4, "big"))
        digest.update(relative)
        digest.update(len(data).to_bytes(8, "big"))
        digest.update(data)
    return digest.hexdigest()


_WRITE_GUARD = threading.Lock()


_WRITE_LOCKS: dict[Path, tuple[int, int]] = {}


def _locks_after_fork() -> None:
    # fork 子进程不能延长父事务的锁寿命，也不能对共享 open description 执行 unlock。
    for _, descriptor in _WRITE_LOCKS.values():
        os.close(descriptor)
    _WRITE_LOCKS.clear()
    _WRITE_GUARD.release()


os.register_at_fork(before=_WRITE_GUARD.acquire,
                    after_in_parent=_WRITE_GUARD.release,
                    after_in_child=_locks_after_fork)


@contextmanager
def write_lock(directory: Path):
    """保护下游目录事务；跨进程/线程冲突立即拒绝，同线程嵌套复用锁。

    锁文件保留原 inode，绝不 unlink；异常或进程退出由关闭描述符释放锁。
    上游仍只使用 Snakemake 锁，sealed marker 的前后复验保持独立。
    """
    directory = Path(directory)
    if directory.is_symlink():
        raise DownstreamQCError(f"QC transaction directory must not be a symlink: {directory}")
    directory = directory.resolve()
    owner = threading.get_ident()
    process_id = os.getpid()
    busy = f"Downstream QC is busy: {directory}; retry after the current writer finishes"
    with _WRITE_GUARD:
        held = _WRITE_LOCKS.get(directory)
        nested = held is not None and held[0] == owner
        if held is not None and not nested:
            raise DownstreamQCError(busy)
        if not nested:
            descriptor = os.open(directory / ".qc-write.lock", os.O_CREAT | os.O_RDWR | os.O_NOFOLLOW, 0o600)
            try:
                if not _qc_lock_stat.S_ISREG(os.fstat(descriptor).st_mode):
                    raise DownstreamQCError(f"QC lock must be a regular file: {directory}")
                try:
                    fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
                except BlockingIOError as exc:
                    raise DownstreamQCError(busy) from exc
                _WRITE_LOCKS[directory] = (owner, descriptor)
            except BaseException:
                os.close(descriptor)
                raise
    try:
        yield
    finally:
        if not nested and os.getpid() == process_id:
            with _WRITE_GUARD:
                del _WRITE_LOCKS[directory]
                # close 足以释放 flock；避免任何清理路径删除另一个进程使用的锁文件。
                os.close(descriptor)


def invalidate_qc_completion(directory: Path) -> None:
    """在首个科学文件替换前撤销整组完成证明；失败不得留下旧成功标记。"""
    for name in ("analysis_summary.json", "qc_summary.json", "integration_summary.json"):
        (Path(directory) / name).unlink(missing_ok=True)


def clear_analysis_outputs(directory: Path, *, input_id: str) -> None:
    """清除失效的已声明产物；完整身份仍匹配的独立Droplet结果保留。"""
    branch = directory / 'droplet_integration'
    clear_branch = branch.exists()
    if clear_branch and (branch / 'integration_summary.json').is_file():
        saved = json.loads((branch / 'integration_summary.json').read_text())
        clear_branch = saved.get('rna_qc_snapshot_id') != input_id
    with ExitStack() as locks:
        # 先取齐相关写锁，再撤标记或删文件；busy不能留下半次清理。
        if clear_branch:
            locks.enter_context(write_lock(branch))
            child = branch / 'dna_doublet'
            if child.exists():
                locks.enter_context(write_lock(child))
        for output in ((directory, branch) if clear_branch else (directory,)):
            if (output/'plots').is_symlink():
                raise DownstreamQCError(f'QC plots directory must not be a symlink: {output}')
        invalidate_qc_completion(directory)
        names = (
            'RNA_Leiden_top30_markers.csv', 'RNA_expression_support_audit.csv',
            'RNA_pre_qc_by_library.csv', 'RNA_post_qc_by_library.csv',
            'RNA_mouse_cluster_marker_review.csv', 'RNA_mouse_marker_annotated.h5ad',
            'RNA_cluster_annotation_review.csv', 'RNA_fine_label_rejection_by_library.csv',
            'RNA_conditioned_reference_label_qc.csv', 'RNA_census_reference_tissue_summary.csv',
            'RNA_Census_consensus_annotated.h5ad', 'RNA_cluster_by_library.csv',
            'RNA_broad_annotation_by_library.csv', 'RNA_detail_annotation_by_library.csv',
            'RNA_detail_level_by_library.csv', 'RNA_sample_tissue_summary.csv',
            'Integrated_QC_Information.csv', 'DNA_RNA_cell_links.csv', 'RNA_cell_qc_metrics.csv',
        )
        for name in names:
            (directory / name).unlink(missing_ok=True)
        for stem in ('01_RNA_DARLIN_QC_overview', '01_RNA_QC_by_library',
                     '02_RNA_UMAP_QC_overview', '03_marker_dotplot',
                     '08_mouse_annotation_umaps', '08_census_annotation_umaps',
                     '09_dname_overlay_umap', '09_dname_qc_umaps', '09_dname_by_celltype'):
            for suffix in ('.pdf', '.png'):
                (directory/'plots'/f'{stem}{suffix}').unlink(missing_ok=True)
        if clear_branch:
            for name in ('integration_summary.json', 'RNA_doublet_scores.csv',
                         'DNA_RNA_cell_links.csv', 'Integrated_QC_Information.csv',
                         'Integrated_QC_singlets.csv', 'RNA_DNA_singlets.h5ad',
                         'SingleCpG_Adata.h5ad'):
                (branch / name).unlink(missing_ok=True)
            for stem in ('01_membership_doublets', '02_joint_qc_doublets', '03_singlet_qc'):
                for suffix in ('.pdf', '.png'):
                    (branch/'plots'/f'{stem}{suffix}').unlink(missing_ok=True)
            if child.exists():
                for name in ('DNA_doublet_summary.json', 'DNA_doublet_scores.csv', 'SingleCpG_Adata.h5ad'):
                    (child / name).unlink(missing_ok=True)
