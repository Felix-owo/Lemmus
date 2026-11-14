"""Lemmus RNA/DARLIN 科学处理与项目契约。"""

import hashlib
import importlib.machinery
import importlib.util
import os
import sys
from pathlib import Path


class _ReadOnlySourceLoader(importlib.machinery.SourceFileLoader):
    # Inspect Python's actual cache policy without regenerating installation files.
    def set_data(self, path, data, **kwargs):
        pass


def _downstream_source_identity():
    """完整科学源码的磁盘身份；只读取文件，不导入科学依赖。"""
    # Python's prefixed cache path preserves the selected import path's spelling.
    core = Path(__file__).absolute().parent
    # Bootstrap must pin source before loading any project child, including core IO.
    paths, pending = [], [core]
    while pending:
        with os.scandir(pending.pop()) as entries:
            for entry in entries:
                if entry.name in {"__pycache__", ".ipynb_checkpoints"}:
                    continue
                path = Path(entry.path)
                if entry.is_dir(follow_symlinks=False):
                    pending.append(path)
                elif path.suffix == '.py':
                    paths.append(path)
    dependencies = {core/name for name in ("__init__.py", "config.py", "io.py", "paths.py", "rna.py", "inputs.py")}
    paths = sorted(set(paths) | dependencies)
    def state(path):
        link = path.lstat()
        value = path.stat()
        return (link.st_dev, link.st_ino, link.st_size, link.st_mtime_ns, link.st_ctime_ns,
                value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns)
    states = {path.relative_to(core).as_posix(): state(path) for path in paths}
    digest = hashlib.sha256()
    for path in paths:
        relative = path.relative_to(core).as_posix().encode()
        data = path.read_bytes()
        name = __name__ + '.' + path.relative_to(core).with_suffix('').as_posix().replace('/', '.')
        if name.endswith('.__init__'):
            name = name[:-9]
        loader = _ReadOnlySourceLoader(name, str(path))
        if loader.get_code(name) != loader.source_to_code(data, str(path)):
            cache = importlib.util.cache_from_source(str(path))
            raise RuntimeError(f'Lemmus bytecode differs from scientific source: {cache}; remove this stale cache and restart the kernel or use a fresh process')
        # All package code must execute current bytes; downstream reuse hashes only its dependencies.
        if path in dependencies or path.is_relative_to(core/'downstream'):
            digest.update(len(relative).to_bytes(4, "big"))
            digest.update(relative)
            digest.update(len(data).to_bytes(8, "big"))
            digest.update(data)
    if any(state(path) != states[path.relative_to(core).as_posix()] for path in paths):
        raise RuntimeError('Lemmus source changed while recording identity; restart the kernel or use a fresh process')
    return digest.hexdigest(), states


# Reloading only the package must not approve already-loaded children from older files.
if '_IMPORTED_DOWNSTREAM_SOURCE_IDENTITY' not in globals():
    if any(name.startswith(__name__ + '.') for name in sys.modules):
        raise RuntimeError('Lemmus children were already loaded without source identity; restart the kernel or use a fresh process')
    _IMPORTED_DOWNSTREAM_SOURCE_IDENTITY = _downstream_source_identity()
