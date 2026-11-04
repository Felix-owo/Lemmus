from __future__ import annotations

from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, version as package_version
from pathlib import Path
from urllib.request import urlretrieve
import json
import random
import re
import time
import os
import subprocess

import anndata
import numpy as np
import pandas as pd
import scipy.sparse as sp

MODULE_VERSION = "v35"
CENSUS_VERSION = "2025-11-08"
SCVI_MODEL_URI = "s3://cellxgene-contrib-public/models/scvi/2025-11-08/homo_sapiens/model.pt"
CELL_ONTOLOGY_URL = "https://purl.obolibrary.org/obo/cl/cl-basic.obo"


@dataclass
class ConsensusAnnotationResult:
    cell_annotation: pd.DataFrame
    cluster_annotation: pd.DataFrame
    method_predictions: pd.DataFrame
    reference_size: int
    census_version: str
    model_uri: str
    reference_tissue_summary: pd.DataFrame | None = None


def download_scvi_model(model_dir: str | Path) -> Path:
    model_dir = Path(model_dir)
    model_dir.mkdir(parents=True, exist_ok=True)
    target = model_dir / "model.pt"
    if not target.exists():
        subprocess.run(
            ["aws", "s3", "cp", "--no-sign-request", "--no-progress",
             SCVI_MODEL_URI, str(target)],
            check=True,
        )
    return target


def download_cell_ontology(path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        urlretrieve(CELL_ONTOLOGY_URL, path)
    return path


def _parse_obo(path: str | Path):
    parents, names = {}, {}
    term = None

    def commit(current):
        if current and current["id"]:
            parents[current["id"]] = tuple(current["parents"])
            names[current["id"]] = current["name"] or current["id"]

    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line == "[Term]":
            commit(term)
            term = {"id": None, "name": None, "parents": []}
        elif line.startswith("[") and line.endswith("]"):
            commit(term)
            term = None
        elif term is not None:
            if line.startswith("id: "):
                term["id"] = line[4:]
            elif line.startswith("name: "):
                term["name"] = line[6:]
            elif line.startswith("is_a: "):
                term["parents"].append(line[6:].split(" ! ", 1)[0])
    commit(term)
    return parents, names


def _ancestors(node: str, parents: dict[str, tuple[str, ...]], cache: dict[str, set[str]]):
    if not node or node not in parents:
        return set()
    if node in cache:
        return cache[node]
    out = set()
    stack = list(parents.get(node, ()))
    while stack:
        x = stack.pop()
        if x in out:
            continue
        out.add(x)
        stack.extend(parents.get(x, ()))
    cache[node] = out
    return out


def _depth(node: str, parents: dict[str, tuple[str, ...]], cache: dict[str, int]):
    if not node or node not in parents:
        return 0
    if node in cache:
        return cache[node]
    ps = parents.get(node, ())
    d = 1 if not ps else 1 + max(_depth(p, parents, cache) for p in ps)
    cache[node] = d
    return d


def _census_var_map(census) -> pd.DataFrame:
    var = (
        census["census_data"]["homo_sapiens"].ms["RNA"].var
        .read(column_names=["feature_id", "feature_name"])
        .concat().to_pandas()
    )
    var["feature_id"] = var["feature_id"].astype(str)
    var["feature_name"] = var["feature_name"].astype(str)
    return var


def _query_for_scvi(adata: anndata.AnnData, census) -> anndata.AnnData:
    X = adata.layers["counts"] if "counts" in adata.layers else adata.X
    q = anndata.AnnData(X=sp.csr_matrix(X.copy()), obs=adata.obs.copy())
    q.var_names = adata.var_names.astype(str)
    if np.mean(q.var_names.str.startswith("ENSG")) <= 0.5:
        vm = _census_var_map(census)
        unique = vm[~vm["feature_name"].duplicated(keep=False)].set_index("feature_name")["feature_id"]
        keep = q.var_names.isin(unique.index)
        q = q[:, keep].copy()
        q.var["gene_symbol"] = q.var_names.astype(str)
        q.var_names = unique.loc[q.var_names].to_numpy()
    else:
        q.var_names = q.var_names.str.replace(r"\.\d+$", "", regex=True)
    q.var_names_make_unique()
    q.obs["batch"] = "query"
    q.obs["joinid"] = np.arange(q.n_obs)
    return q


def diagnose_census_gene_overlap(
    adata: anndata.AnnData,
    model_dir: str | Path,
    census_version: str = CENSUS_VERSION,
) -> dict[str, float | int]:
    import cellxgene_census
    import scvi

    model_dir = Path(model_dir)
    if not (model_dir / "model.pt").is_file():
        raise FileNotFoundError(f"Missing Census scVI model: {model_dir / 'model.pt'}")

    with cellxgene_census.open_soma(census_version=census_version) as census:
        q = _query_for_scvi(adata, census)

    query_before = int(adata.n_vars)
    query_after = int(q.n_vars)
    query_names = pd.Index(q.var_names.astype(str))

    reference_var_names = scvi.model.SCVI.prepare_query_anndata(
        q.copy(),
        str(model_dir),
        return_reference_var_names=True,
    )
    reference_var_names = pd.Index(np.asarray(reference_var_names).astype(str))

    overlap = query_names.intersection(reference_var_names)
    reference_n = int(len(reference_var_names))
    overlap_n = int(len(overlap))
    overlap_pct = 100.0 * overlap_n / max(reference_n, 1)

    print("=== Census gene-overlap diagnostic ===")
    print(f"Full query genes before Census mapping : {query_before:,}")
    print(f"Query genes after symbol→Ensembl map  : {query_after:,}")
    print(f"Census model reference genes          : {reference_n:,}")
    model_norm = {_normalize_ensembl_gene_id(x) for x in reference_var_names} - {""}
    query_norm = {_normalize_ensembl_gene_id(x) for x in query_names} - {""}
    norm_overlap = len(model_norm & query_norm)
    print(f"Normalized Ensembl overlap: {norm_overlap:,} / {len(model_norm):,} = {norm_overlap/max(1,len(model_norm)):.2%}")
    print(f"Model/query overlap                    : {overlap_n:,}")
    print(f"Model/query overlap                    : {overlap_pct:.2f}%")

    if overlap_pct < 20.0:
        print(
            "WARNING: Census model/query gene overlap is very low; "
            "do not trust downstream annotation until the gene universe/mapping is fixed."
        )

    return {
        "query_genes_before_mapping": query_before,
        "query_genes_after_mapping": query_after,
        "reference_genes": reference_n,
        "overlap_genes": overlap_n,
        "overlap_pct": float(overlap_pct),
    }


def _package_version(name: str) -> str:
    try:
        return package_version(name)
    except PackageNotFoundError:
        return "unknown"


def _set_deterministic_seed(seed: int) -> dict[str, object]:
    """Set every RNG used by scArches and require deterministic CPU kernels."""
    import scvi
    import torch

    seed = int(seed)
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    scvi.settings.seed = seed
    torch.use_deterministic_algorithms(True)
    return {
        "seed": seed,
        "accelerator": "cpu",
        "torch_deterministic_algorithms": True,
        "numpy_version": np.__version__,
        "torch_version": str(torch.__version__),
        "scvi_tools_version": _package_version("scvi-tools"),
    }


def project_query_scvi(
    adata: anndata.AnnData,
    model_dir: str | Path,
    census_version: str = CENSUS_VERSION,
    latent_key: str = "X_census_scvi",
    min_reference_var_overlap: float = 0.20,
    query_max_epochs: int = 200,
    seed: int = 42,
) -> np.ndarray:
    import cellxgene_census
    import scvi

    deterministic = _set_deterministic_seed(seed)

    model_dir = Path(model_dir)
    if not (model_dir / "model.pt").is_file():
        raise FileNotFoundError(model_dir / "model.pt")

    with cellxgene_census.open_soma(census_version=census_version) as census:
        q = _query_for_scvi(adata, census)

    query_before = int(adata.n_vars)
    query_after = int(q.n_vars)
    query_names = pd.Index(q.var_names.astype(str))

    reference_var_names = scvi.model.SCVI.prepare_query_anndata(
        q.copy(),
        str(model_dir),
        return_reference_var_names=True,
    )
    reference_var_names = pd.Index(np.asarray(reference_var_names).astype(str))

    overlap = query_names.intersection(reference_var_names)
    reference_n = int(len(reference_var_names))
    overlap_n = int(len(overlap))
    overlap_fraction = overlap_n / max(reference_n, 1)

    print("=== Census/scArches gene-overlap diagnostic ===", flush=True)
    print(f"Module version/path                     : {MODULE_VERSION} | {Path(__file__).resolve()}", flush=True)
    print(f"Incoming AnnData shape                  : {adata.n_obs:,} cells × {adata.n_vars:,} genes", flush=True)
    print(f"Full query genes before Census mapping : {query_before:,}", flush=True)
    print(f"Query genes after symbol→Ensembl map  : {query_after:,}", flush=True)
    print(f"Census model reference genes          : {reference_n:,}", flush=True)
    print(f"Model/query overlap                    : {overlap_n:,}", flush=True)
    print(f"Model/query overlap                    : {100.0 * overlap_fraction:.2f}%", flush=True)

    adata.uns["census_gene_overlap"] = {
        "module_version": MODULE_VERSION,
        "module_path": str(Path(__file__).resolve()),
        "query_genes_before_mapping": query_before,
        "query_genes_after_mapping": query_after,
        "reference_genes": reference_n,
        "overlap_genes": overlap_n,
        "overlap_pct": float(100.0 * overlap_fraction),
        "seed": int(seed),
    }

    if overlap_fraction < float(min_reference_var_overlap):
        raise RuntimeError(
            "scArches projection aborted: model/query reference-var overlap is "
            f"{100.0 * overlap_fraction:.2f}% ({overlap_n:,}/{reference_n:,}), "
            f"below the required {100.0 * min_reference_var_overlap:.1f}%."
        )

    # scArches query-to-reference mapping:
    # prepare → load_query_data → minimal fine tuning → latent representation.
    scvi.model.SCVI.prepare_query_anndata(q, str(model_dir))
    query_model = scvi.model.SCVI.load_query_data(q, str(model_dir))

    print(
        f"[scArches] Fine-tuning query model for {query_max_epochs} epochs "
        "with weight_decay=0.0 ...",
        flush=True,
    )
    train_t0 = time.perf_counter()
    query_model.train(
        max_epochs=int(query_max_epochs),
        plan_kwargs={"weight_decay": 0.0},
        accelerator="cpu",
        devices=1,
        enable_progress_bar=False,
    )
    print(
        f"[scArches] Query fine-tuning complete | "
        f"{_fmt_elapsed(time.perf_counter() - train_t0)}",
        flush=True,
    )

    z = query_model.get_latent_representation().astype(np.float32)
    adata.obsm[latent_key] = z
    adata.uns["census_scarches_training"] = {
        **deterministic,
        "module_version": MODULE_VERSION,
        "model_uri": SCVI_MODEL_URI,
        "census_version": str(census_version),
        "query_max_epochs": int(query_max_epochs),
        "weight_decay": 0.0,
        "devices": 1,
    }
    return z



def _progress(stage: int, message: str) -> None:
    print(f"[Census {stage}/6] {message}", flush=True)


def _fmt_elapsed(seconds: float) -> str:
    seconds = max(0, int(seconds))
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    return f"{h:02d}:{m:02d}:{s:02d}" if h else f"{m:02d}:{s:02d}"


def _reference_cache_meta(
    census_version: str,
    candidate_cells: int,
    max_cells_per_type: int,
    max_reference_cells: int,
    normal_only: bool,
    tissues: tuple[str, ...] | None,
    seed: int,
    metadata_block_size: int,
    embedding_batch_size: int,
) -> dict:
    return {
        "cache_schema": 3,
        "census_version": str(census_version),
        "candidate_cells": int(candidate_cells),
        "max_cells_per_type": int(max_cells_per_type),
        "max_reference_cells": int(max_reference_cells),
        "normal_only": bool(normal_only),
        "tissues": sorted(map(str, tissues or ())),
        "seed": int(seed),
        "embedding": "scvi",
        "organism": "Homo sapiens",
        "sampling": "distributed_contiguous_blocks_diversity_cap_v2",
        "metadata_block_size": int(metadata_block_size),
        "embedding_batch_size": int(embedding_batch_size),
    }


def _reference_cache_matches(path: str | Path, expected: dict) -> bool:
    path = Path(path)
    if not path.is_file():
        return False
    try:
        with np.load(path, allow_pickle=False) as ref:
            if "cache_meta" not in ref.files:
                return False
            raw = ref["cache_meta"]
            text = str(raw.item() if raw.ndim == 0 else raw.reshape(-1)[0])
            return json.loads(text) == expected
    except (OSError, ValueError, KeyError, json.JSONDecodeError):
        return False


def _sample_obs_blocks(obs, candidate_cells: int, block_size: int, seed: int) -> pd.DataFrame:
    """Sample Census metadata via contiguous SOMA ranges rather than millions of random coordinates."""
    rng = np.random.default_rng(seed)
    lo, hi = map(int, obs.domain[0])
    span = hi - lo + 1
    target = min(int(candidate_cells), span)
    block_size = max(1, min(int(block_size), target))
    n_blocks = max(1, int(np.ceil(target / block_size)))
    edges = np.linspace(lo, hi + 1, n_blocks + 1, dtype=np.int64)
    frames = []
    columns = [
        "soma_joinid", "cell_type", "cell_type_ontology_term_id",
        "tissue_general", "disease", "is_primary_data", "dataset_id",
    ]
    stage_t0 = time.perf_counter()
    rows_read = 0
    for i in range(n_blocks):
        block_t0 = time.perf_counter()
        seg_lo, seg_hi = int(edges[i]), int(edges[i + 1])
        width = min(block_size, max(1, seg_hi - seg_lo))
        max_start = max(seg_lo, seg_hi - width)
        start = int(rng.integers(seg_lo, max_start + 1)) if max_start > seg_lo else seg_lo
        stop = min(hi, start + width - 1)
        print(
            f"    metadata block {i + 1}/{n_blocks} | "
            f"SOMA {start:,}:{stop:,} | requesting...",
            flush=True,
        )
        frame = obs.read(coords=(slice(start, stop),), column_names=columns).concat().to_pandas()
        rows_read += len(frame)
        if not frame.empty:
            frames.append(frame)
        print(
            f"    metadata block {i + 1}/{n_blocks} | "
            f"{len(frame):,} rows | cumulative {rows_read:,} | "
            f"block {_fmt_elapsed(time.perf_counter() - block_t0)} | "
            f"stage {_fmt_elapsed(time.perf_counter() - stage_t0)}",
            flush=True,
        )
    if not frames:
        return pd.DataFrame(columns=columns)
    return pd.concat(frames, ignore_index=True).drop_duplicates("soma_joinid")


def _balanced_reference_sample(
    meta: pd.DataFrame,
    max_cells_per_type: int,
    max_reference_cells: int,
    seed: int,
) -> pd.DataFrame:
    """Per-type balanced sampling followed by a diversity-first global cap."""
    if meta.empty:
        return meta.copy()

    rng = np.random.default_rng(seed)
    sampled = meta.assign(_sample_order=rng.random(len(meta)))
    sampled = sampled.sort_values(["cell_type", "_sample_order"], kind="stable")
    sampled = (
        sampled.groupby("cell_type", observed=True, sort=True)
        .head(int(max_cells_per_type))
        .drop(columns="_sample_order")
        .reset_index(drop=True)
    )

    limit = int(max_reference_cells)
    if limit <= 0 or len(sampled) <= limit:
        return sampled

    # Round-robin by within-type rank. This preserves cell-type diversity
    # instead of letting abundant labels dominate a simple global random cap.
    sampled = sampled.assign(
        _type_rank=sampled.groupby("cell_type", observed=True, sort=True).cumcount(),
        _global_tie=rng.random(len(sampled)),
    )
    sampled = (
        sampled.sort_values(
            ["_type_rank", "cell_type", "_global_tie"],
            kind="stable",
        )
        .head(limit)
        .drop(columns=["_type_rank", "_global_tie"])
        .reset_index(drop=True)
    )
    return sampled



def _fetch_scvi_embeddings(
    meta: pd.DataFrame,
    census_version: str,
    batch_size: int,
) -> tuple[np.ndarray, pd.DataFrame]:
    from cellxgene_census.experimental import (
        get_embedding,
        get_embedding_metadata_by_name,
    )

    if meta.empty:
        raise RuntimeError("No reference cells selected before scVI embedding fetch.")

    embedding_meta = get_embedding_metadata_by_name(
        "scvi",
        "homo_sapiens",
        census_version=census_version,
    )
    if not embedding_meta:
        raise RuntimeError(
            f"No official scVI embedding metadata found for Census {census_version}."
        )

    embedding_uri = (
        "s3://cellxgene-contrib-public/contrib/cell-census/soma/"
        f"{embedding_meta['census_version']}/{embedding_meta['id']}"
    )

    print(
        f"    embedding source | census={embedding_meta['census_version']} | "
        f"id={embedding_meta['id']}",
        flush=True,
    )
    print(
        "    retrieval=experimental.get_embedding | scVI only | "
        "reuse Stage-1 metadata",
        flush=True,
    )

    ids = meta["soma_joinid"].to_numpy(np.int64)
    total = len(ids)
    total_batches = max(1, (total + batch_size - 1) // batch_size)
    stage_t0 = time.perf_counter()

    embeddings: list[np.ndarray] = []
    obs_frames: list[pd.DataFrame] = []

    for batch_no, start in enumerate(range(0, total, batch_size), start=1):
        batch_t0 = time.perf_counter()
        stop = min(start + batch_size, total)

        local_order = np.argsort(ids[start:stop], kind="stable")
        batch_ids = ids[start:stop][local_order]
        batch_meta = meta.iloc[start:stop].iloc[local_order].reset_index(drop=True)

        print(
            f"    embedding batch {batch_no}/{total_batches} | "
            f"{len(batch_ids):,} cells | requesting 50-D scVI only...",
            flush=True,
        )

        arr = get_embedding(
            embedding_meta["census_version"],
            embedding_uri,
            batch_ids,
        )
        arr = np.asarray(arr, dtype=np.float32)

        if arr.ndim != 2:
            raise RuntimeError(
                f"Unexpected scVI embedding shape for batch {batch_no}: {arr.shape}"
            )
        if arr.shape[0] != len(batch_ids):
            raise RuntimeError(
                f"scVI row mismatch in batch {batch_no}: "
                f"requested {len(batch_ids):,}, received {arr.shape[0]:,}"
            )
        if arr.shape[1] != 50:
            raise RuntimeError(
                f"Unexpected scVI dimension in batch {batch_no}: "
                f"expected 50, received {arr.shape[1]}"
            )

        embeddings.append(arr)
        obs_frames.append(batch_meta)

        print(
            f"    embedding batch {batch_no}/{total_batches} | "
            f"{stop:,}/{total:,} cells | shape={arr.shape} | "
            f"matrix={arr.nbytes / (1024 ** 2):.1f} MiB | "
            f"batch {_fmt_elapsed(time.perf_counter() - batch_t0)} | "
            f"stage {_fmt_elapsed(time.perf_counter() - stage_t0)}",
            flush=True,
        )

    z = np.vstack(embeddings)
    ref_obs = pd.concat(obs_frames, ignore_index=True)

    if len(z) != len(ref_obs):
        raise RuntimeError(
            f"Embedding/metadata row mismatch: {len(z):,} vs {len(ref_obs):,}"
        )

    print(
        f"    embedding fetch complete | {z.shape[0]:,} cells × {z.shape[1]} dims | "
        f"matrix={z.nbytes / (1024 ** 2):.1f} MiB | "
        f"stage {_fmt_elapsed(time.perf_counter() - stage_t0)}",
        flush=True,
    )
    return z, ref_obs


def build_reference_cache(
    cache_path: str | Path,
    census_version: str = CENSUS_VERSION,
    candidate_cells: int = 10_000_000,
    max_cells_per_type: int = 3_000,
    max_reference_cells: int = 500_000,
    normal_only: bool = True,
    tissues: tuple[str, ...] | None = None,
    seed: int = 42,
    metadata_block_size: int = 500_000,
    embedding_batch_size: int = 50_000,
) -> Path:
    import cellxgene_census

    build_t0 = time.perf_counter()
    cache_path = Path(cache_path)
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    print(
        f"[Census cache build] version={census_version} | target={candidate_cells:,} | "
        f"metadata_block={metadata_block_size:,} | embedding_batch={embedding_batch_size:,} | "
        f"max/type={max_cells_per_type:,} | max/reference={max_reference_cells:,}",
        flush=True,
    )
    print(
        "[Census cache build] Remote CELLxGENE/TileDB-SOMA reads are active during stages 1 and 3.",
        flush=True,
    )

    with cellxgene_census.open_soma(census_version=census_version) as census:
        _progress(1, f"Reading ~{candidate_cells:,} candidate cells from Census metadata...")
        obs = census["census_data"]["homo_sapiens"].obs
        meta = _sample_obs_blocks(obs, candidate_cells, metadata_block_size, seed)
        required_meta = {
            "soma_joinid", "cell_type", "cell_type_ontology_term_id",
            "tissue_general", "disease", "is_primary_data", "dataset_id",
        }
        missing_meta = required_meta.difference(meta.columns)
        if missing_meta:
            raise RuntimeError(
                f"Census metadata is missing required columns: {sorted(missing_meta)}"
            )
        n_raw = len(meta)
        meta = meta[meta["is_primary_data"].eq(True)].copy()
        n_primary = len(meta)
        if normal_only:
            meta = meta[meta["disease"].astype(str).eq("normal")]
        n_normal = len(meta)
        if tissues:
            meta = meta[meta["tissue_general"].astype(str).isin(tissues)]
        n_tissue = len(meta)
        bad = {"unknown", "native cell", "cell", "primary cultured cell"}
        meta = meta[~meta["cell_type"].astype(str).isin(bad)]
        print(
            f"    filtering | raw {n_raw:,} → primary {n_primary:,} → "
            f"normal {n_normal:,} → tissue {n_tissue:,} → usable {len(meta):,}",
            flush=True,
        )
        if meta.empty:
            raise RuntimeError("Census reference filtering returned zero cells")

        _progress(2, f"Balancing {len(meta):,} eligible cells; max {max_cells_per_type:,} cells/type...")
        stage2_t0 = time.perf_counter()
        pre_balance_types = int(meta["cell_type"].nunique())
        meta = _balanced_reference_sample(
            meta,
            max_cells_per_type=max_cells_per_type,
            max_reference_cells=max_reference_cells,
            seed=seed,
        )
        # Stage-3 remote I/O optimization: globally sort selected cells by SOMA join ID.
        # This preserves every metadata/label row while making embedding batches spatially localized.
        meta = meta.sort_values("soma_joinid", kind="stable").reset_index(drop=True)
        print("    Stage-3 I/O order: global soma_joinid sort before batching", flush=True)
        ids = meta["soma_joinid"].to_numpy(np.int64)
        n_batches = max(1, (len(ids) + embedding_batch_size - 1) // embedding_batch_size)
        print(
            f"    selected {len(ids):,} cells across {meta['cell_type'].nunique():,} cell types "
            f"(eligible types={pre_balance_types:,}) | "
            f"embedding batches={n_batches} | "
            f"stage {_fmt_elapsed(time.perf_counter() - stage2_t0)}",
            flush=True,
        )

        _progress(3, f"Fetching official 50-D scVI embeddings via experimental.get_embedding(); {n_batches} batches × up to {embedding_batch_size:,} cells...")
        z, ref_obs = _fetch_scvi_embeddings(meta, census_version, embedding_batch_size)

    ok = np.isfinite(z).all(axis=1)
    cache_meta = _reference_cache_meta(
        census_version, candidate_cells, max_cells_per_type, max_reference_cells,
        normal_only, tissues, seed, metadata_block_size, embedding_batch_size,
    )
    np.savez_compressed(
        cache_path,
        cache_meta=np.asarray(json.dumps(cache_meta, sort_keys=True), dtype=np.str_),
        embedding=np.asarray(z[ok], dtype=np.float32),
        cell_type=np.asarray(ref_obs["cell_type"].astype(str).to_numpy()[ok], dtype=np.str_),
        ontology=np.asarray(ref_obs["cell_type_ontology_term_id"].astype(str).to_numpy()[ok], dtype=np.str_),
        tissue=np.asarray(ref_obs["tissue_general"].astype(str).to_numpy()[ok], dtype=np.str_),
        dataset_id=np.asarray(ref_obs["dataset_id"].astype(str).to_numpy()[ok], dtype=np.str_),
    )
    size_mib = cache_path.stat().st_size / (1024 ** 2)
    print(
        f"[Census cache ready] {cache_path} | {ok.sum():,} finite embeddings | "
        f"{size_mib:.1f} MiB | total {_fmt_elapsed(time.perf_counter() - build_t0)}",
        flush=True,
    )
    return cache_path

def _centroid_predict(Xr, yr, Xq):
    from sklearn.preprocessing import normalize

    labels = np.unique(yr)
    centroids = np.vstack([Xr[yr == x].mean(axis=0) for x in labels])
    centroids = normalize(centroids)
    q = normalize(Xq)
    sim = q @ centroids.T
    p = np.exp((sim - sim.max(axis=1, keepdims=True)) / 0.08)
    p /= p.sum(axis=1, keepdims=True)
    return labels, p


def _fit_predictors(Xr, yr, Xq, seed=42):
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.preprocessing import StandardScaler

    scaler = StandardScaler()
    Xrs = scaler.fit_transform(Xr)
    Xqs = scaler.transform(Xq)

    models = {
        "knn_scvi": KNeighborsClassifier(n_neighbors=25, weights="distance", n_jobs=-1),
        "logreg_scvi": LogisticRegression(
            max_iter=600, class_weight="balanced", solver="lbfgs",
        ),
        "rf_scvi": RandomForestClassifier(
            n_estimators=400, n_jobs=-1, random_state=seed,
            class_weight="balanced_subsample", max_features="sqrt",
        ),
    }

    out = {}
    for name, model in models.items():
        train = Xrs if name != "rf_scvi" else Xr
        query = Xqs if name != "rf_scvi" else Xq
        model.fit(train, yr)
        out[name] = (model.classes_.astype(str), model.predict_proba(query))

    out["centroid_scvi"] = _centroid_predict(Xr, yr, Xq)
    return out


def _normalize_cl_id(value: object) -> str:
    text = str(value or "")
    match = re.search(r"CL:\d{7}", text)
    return match.group(0) if match else ""


def _ontology_consensus(
    method_outputs: dict[str, tuple[np.ndarray, np.ndarray]],
    ontology_map: dict[str, str],
    parents: dict[str, tuple[str, ...]],
    names: dict[str, str],
    ancestor_weight: float = 0.60,
):
    methods = list(method_outputs)
    n = next(iter(method_outputs.values()))[1].shape[0]
    pred_table = pd.DataFrame(index=np.arange(n))
    prob_table = pd.DataFrame(index=np.arange(n))

    for method, (classes, prob) in method_outputs.items():
        top = prob.argmax(axis=1)
        pred_table[method] = classes[top]
        prob_table[method] = prob[np.arange(n), top]

    ancestor_cache, depth_cache = {}, {}
    label_by_ontology = {}
    for label, oid in ontology_map.items():
        if oid:
            label_by_ontology.setdefault(oid, label)

    final_label, final_oid, agreement, mean_prob = [], [], [], []
    ontology_score, ontology_margin, exact_votes, exact_vote_margin = [], [], [], []

    for i in range(n):
        labels = pred_table.iloc[i].astype(str).tolist()
        probs = prob_table.iloc[i].to_numpy(float)
        exact = pd.Series(labels).value_counts()
        top_exact = str(exact.index[0])
        top_votes = int(exact.iloc[0])
        second_votes = int(exact.iloc[1]) if len(exact) > 1 else 0

        scores = {}
        for label in labels:
            oid = ontology_map.get(label, "")
            if not oid:
                continue
            scores[oid] = scores.get(oid, 0.0) + 1.0
            for anc in _ancestors(oid, parents, ancestor_cache):
                scores[anc] = scores.get(anc, 0.0) + ancestor_weight

        if scores:
            ranked = sorted(
                scores.items(),
                key=lambda kv: (kv[1], _depth(kv[0], parents, depth_cache), kv[0]),
                reverse=True,
            )
            oid, score = ranked[0]
            margin = score - (ranked[1][1] if len(ranked) > 1 else 0.0)
            if top_votes >= 3:
                label = top_exact
                oid = ontology_map.get(label, oid)
            else:
                label = label_by_ontology.get(oid, names.get(oid, top_exact))
        else:
            label, oid, score, margin = top_exact, ontology_map.get(top_exact, ""), float(top_votes), 0.0

        final_label.append(label)
        final_oid.append(oid)
        agreement.append(top_votes / len(methods))
        mean_prob.append(float(probs.mean()))
        ontology_score.append(float(score))
        ontology_margin.append(float(margin))
        exact_votes.append(top_votes)
        exact_vote_margin.append(top_votes - second_votes)

    result = pd.DataFrame({
        "consensus_label": final_label,
        "consensus_ontology_id": final_oid,
        "method_agreement": agreement,
        "mean_method_probability": mean_prob,
        "ontology_score": ontology_score,
        "ontology_margin": ontology_margin,
        "exact_votes": exact_votes,
        "exact_vote_margin": exact_vote_margin,
    })
    return result, pred_table, prob_table



def _normalize_ensembl_gene_id(value: object) -> str:
    x = str(value).strip().split("|", 1)[0]
    x = re.sub(r"\.\d+$", "", x)
    x = re.sub(r"_PAR_Y$", "", x)
    return x


def _broad_lineage_from_label(label: str) -> str:
    x = str(label).lower()
    rules = (
        ("Erythroid cell", (
            "erythroid", "erythrocyte", "erythroblast", "red blood cell", "normoblast",
        )),
        ("Megakaryocyte", ("megakaryocyte", "platelet")),
        ("B-lineage cell", ("b cell", "b-cell", "pre-b", "plasma cell", "plasmablast")),
        ("T/NK cell", ("t cell", "t-cell", "natural killer", "nk cell")),
        ("Myeloid cell", (
            "monocyte", "macrophage", "dendritic", "myeloid", "neutrophil",
            "granulocyte", "mast cell", "microglia", "microglial",
        )),
        ("Endothelial cell", ("endothelial", "endothelium")),
        ("Epithelial cell", ("epithelial", "epithelium", "keratinocyte", "pneumocyte")),
        ("Schwann/glial cell", (
            "schwann", "glial", "astrocyte", "oligodendrocyte",
            "oligodendroglial", "satellite glial",
        )),
        ("Neural cell", (
            "neuron", "neuronal", "neural", "neuroblast", "neuroplacodal",
            "neuroplacode", "sensory",
        )),
        ("Fibroblast/stromal cell", ("fibroblast", "stromal", "mesenchymal", "pericyte", "smooth muscle", "myofibroblast")),
        ("Hematopoietic progenitor", (
            "hematopoietic stem", "hematopoietic progenitor",
            "hematopoietic precursor", "multipotent progenitor",
        )),
    )
    for broad, terms in rules:
        if any(
            re.search(rf"(?<![a-z0-9]){re.escape(term)}(?![a-z0-9])", x)
            for term in terms
        ):
            return broad
    return "Other"


def _scarches_weighted_knn_transfer(
    Xr: np.ndarray,
    yr: np.ndarray,
    Xq: np.ndarray,
    n_neighbors: int = 50,
    uncertainty_threshold: float = 0.50,
) -> pd.DataFrame:
    from sklearn.neighbors import NearestNeighbors

    Xr = np.asarray(Xr, dtype=np.float32)
    Xq = np.asarray(Xq, dtype=np.float32)
    yr = np.asarray(yr).astype(str)

    k = max(1, min(int(n_neighbors), len(Xr)))
    knn = NearestNeighbors(
        n_neighbors=k,
        algorithm="brute",
        metric="euclidean",
        n_jobs=-1,
    )
    knn.fit(Xr)
    distances, indices = knn.kneighbors(Xq, return_distance=True)

    # scArches weighted-kNN implementation:
    # std = sd(neighbor distances)
    # weights = exp(-distance / (2/std)^2), normalized per query cell.
    std = np.std(distances, axis=1)
    std = np.maximum(std, np.finfo(np.float32).eps)
    scale = (2.0 / std) ** 2
    weights = np.exp(-distances / scale[:, None])
    weights /= np.maximum(weights.sum(axis=1, keepdims=True), np.finfo(np.float32).eps)

    pred_label = np.empty(len(Xq), dtype=object)
    pred_prob = np.zeros(len(Xq), dtype=np.float32)

    for i in range(len(Xq)):
        labels = yr[indices[i]]
        probs = {}
        for label, weight in zip(labels, weights[i], strict=False):
            probs[label] = probs.get(label, 0.0) + float(weight)
        best_label, best_prob = max(probs.items(), key=lambda kv: kv[1])
        pred_label[i] = best_label
        pred_prob[i] = best_prob

    uncertainty = 1.0 - pred_prob
    resolved = uncertainty <= float(uncertainty_threshold)
    final_label = pred_label.copy()
    final_label[~resolved] = "Unknown"

    return pd.DataFrame({
        "predicted_label": pred_label.astype(str),
        "prediction_probability": pred_prob,
        "uncertainty": uncertainty,
        "resolved": resolved,
        "cell_type": final_label.astype(str),
    })



def _weighted_knn_transfer_labels(
    Xr: np.ndarray,
    yr: np.ndarray,
    Xq: np.ndarray,
    n_neighbors: int = 50,
    uncertainty_threshold: float = 0.50,
) -> pd.DataFrame:
    """Weighted-kNN transfer for an arbitrary categorical reference label."""
    return _scarches_weighted_knn_transfer(
        Xr,
        np.asarray(yr).astype(str),
        Xq,
        n_neighbors=n_neighbors,
        uncertainty_threshold=uncertainty_threshold,
    )


def _hierarchical_weighted_knn_transfer(Xr, yr, Xq, n_neighbors=50, uncertainty_threshold=0.50):
    yr = np.asarray(yr).astype(str)
    broad_ref = np.asarray([_broad_lineage_from_label(x) for x in yr], dtype=str)
    broad = _scarches_weighted_knn_transfer(
        Xr, broad_ref, Xq, n_neighbors=n_neighbors,
        uncertainty_threshold=uncertainty_threshold,
    )
    global_sub = _scarches_weighted_knn_transfer(
        Xr, yr, Xq, n_neighbors=n_neighbors,
        uncertainty_threshold=uncertainty_threshold,
    )
    out = global_sub.copy()
    for lineage in pd.unique(broad["cell_type"]):
        if lineage in {"Unknown", "Other"}:
            continue
        qmask = broad["cell_type"].to_numpy() == lineage
        rmask = broad_ref == lineage
        if qmask.sum() == 0 or rmask.sum() < 25:
            continue
        local = _scarches_weighted_knn_transfer(
            Xr[rmask], yr[rmask], Xq[qmask],
            n_neighbors=min(n_neighbors, int(rmask.sum())),
            uncertainty_threshold=uncertainty_threshold,
        )
        out.loc[qmask, :] = local.to_numpy()
    return out.reset_index(drop=True), broad.reset_index(drop=True)


def annotate_from_census_consensus(
    adata: anndata.AnnData,
    model_dir: str | Path,
    reference_cache: str | Path,
    ontology_path: str | Path,
    cluster_key: str = "leiden",
    census_version: str = CENSUS_VERSION,
    candidate_cells: int = 10_000_000,
    max_cells_per_type: int = 3_000,
    max_reference_cells: int = 500_000,
    normal_only: bool = True,
    tissues: tuple[str, ...] | None = None,
    min_agreement: float = 0.50,
    min_mean_probability: float = 0.35,
    min_ontology_depth: int = 3,
    seed: int = 42,
    metadata_block_size: int = 50_000,
    embedding_batch_size: int = 25_000,
    weighted_knn_k: int = 50,
    uncertainty_threshold: float = 0.50,
    query_max_epochs: int = 200,
) -> ConsensusAnnotationResult:
    # ontology_path/min_* are retained in the signature for notebook compatibility,
    # but are not used as primary gates in the scArches paper workflow.
    reference_cache = Path(reference_cache)

    expected_cache_meta = _reference_cache_meta(
        census_version, candidate_cells, max_cells_per_type, max_reference_cells,
        normal_only, tissues, seed, metadata_block_size, embedding_batch_size,
    )
    cache_valid = _reference_cache_matches(reference_cache, expected_cache_meta)
    if not cache_valid:
        print(
            f"[Census cache] MISS/INVALID: {reference_cache}\n"
            "    building a new compact reference; stages 1 and 3 use network I/O.",
            flush=True,
        )
        tmp_cache = reference_cache.with_name(
            f".{reference_cache.stem}.rebuild.{os.getpid()}.npz"
        )
        tmp_cache.unlink(missing_ok=True)
        try:
            build_reference_cache(
                tmp_cache,
                census_version=census_version,
                candidate_cells=candidate_cells,
                max_cells_per_type=max_cells_per_type,
                max_reference_cells=max_reference_cells,
                normal_only=normal_only,
                tissues=tissues,
                seed=seed,
                metadata_block_size=metadata_block_size,
                embedding_batch_size=embedding_batch_size,
            )
            os.replace(tmp_cache, reference_cache)
        finally:
            tmp_cache.unlink(missing_ok=True)
    else:
        _progress(1, f"Reference cache valid: {reference_cache.name}")
        _progress(2, "Balanced reference already cached; skipping metadata sampling.")
        _progress(3, "scVI reference embeddings already cached; skipping cloud fetch.")

    _progress(
        4,
        f"scArches query-to-reference mapping for {adata.n_obs:,} cells "
        f"({query_max_epochs} fine-tuning epochs)...",
    )
    stage4_t0 = time.perf_counter()
    Xq = project_query_scvi(
        adata,
        model_dir,
        census_version,
        query_max_epochs=query_max_epochs,
        seed=seed,
    )
    print(
        f"    query latent shape: {Xq.shape[0]:,} cells × {Xq.shape[1]} dimensions | "
        f"stage {_fmt_elapsed(time.perf_counter() - stage4_t0)}",
        flush=True,
    )

    with np.load(reference_cache, allow_pickle=False) as ref:
        Xr = ref["embedding"].astype(np.float32)
        yr = ref["cell_type"].astype(str)
        ref_ontology = ref["ontology"].astype(str)
        ref_tissue = (
            ref["tissue"].astype(str)
            if "tissue" in ref.files
            else np.full(len(Xr), "Unknown", dtype=str)
        )

    reference_tissue_summary = (
        pd.Series(ref_tissue, name="census_tissue_general")
        .value_counts(dropna=False)
        .rename_axis("census_tissue_general")
        .rename("n_reference_cells")
        .to_frame()
    )
    reference_tissue_summary["fraction"] = (
        reference_tissue_summary["n_reference_cells"] / max(1, len(ref_tissue))
    )

    _progress(
        5,
        f"Hierarchical weighted-kNN label transfer: k={weighted_knn_k}, "
        f"reference={len(Xr):,} cells...",
    )
    stage5_t0 = time.perf_counter()
    transfer, broad_transfer = _hierarchical_weighted_knn_transfer(
        Xr,
        yr,
        Xq,
        n_neighbors=weighted_knn_k,
        uncertainty_threshold=uncertainty_threshold,
    )
    tissue_transfer = _weighted_knn_transfer_labels(
        Xr,
        ref_tissue,
        Xq,
        n_neighbors=weighted_knn_k,
        uncertainty_threshold=uncertainty_threshold,
    )
    print(
        f"    weighted-kNN complete | "
        f"median uncertainty={transfer['uncertainty'].median():.3f} | "
        f"stage {_fmt_elapsed(time.perf_counter() - stage5_t0)}",
        flush=True,
    )

    # Attach ontology ID only as metadata from the reference label; it is not a gate.
    ontology_map = (
        pd.DataFrame({"cell_type": yr, "ontology": ref_ontology})
        .drop_duplicates("cell_type")
        .set_index("cell_type")["ontology"]
        .to_dict()
    )

    _progress(
        6,
        f"Applying scArches uncertainty threshold <= {uncertainty_threshold:.2f} "
        "(paper uses >0.5 as unknown)...",
    )

    adata.obs["cell_type_broad_atlas"] = broad_transfer["cell_type"].to_numpy()
    adata.obs["cell_type_broad_probability"] = broad_transfer["prediction_probability"].to_numpy()
    adata.obs["cell_type_broad_uncertainty"] = broad_transfer["uncertainty"].to_numpy()
    adata.obs["cell_type_auto"] = transfer["cell_type"].to_numpy()
    adata.obs["cell_type"] = transfer["cell_type"].to_numpy()
    adata.obs["cell_type_probability"] = transfer["prediction_probability"].to_numpy()
    adata.obs["cell_type_uncertainty"] = transfer["uncertainty"].to_numpy()
    adata.obs["cell_type_confidence"] = np.where(
        transfer["uncertainty"].to_numpy() <= 0.25,
        "High",
        np.where(transfer["resolved"].to_numpy(), "Medium", "Low"),
    )
    # Explicit semantic aliases.  The legacy cell_type* columns are retained
    # for compatibility, but these names make the annotation level unambiguous.
    adata.obs["atlas_cell_type_cell"] = adata.obs["cell_type"].astype(str)
    adata.obs["atlas_probability_cell"] = adata.obs["cell_type_probability"].to_numpy()
    adata.obs["atlas_uncertainty_cell"] = adata.obs["cell_type_uncertainty"].to_numpy()
    adata.obs["atlas_confidence_cell"] = adata.obs["cell_type_confidence"].astype(str)
    adata.obs["atlas_tissue_cell"] = tissue_transfer["cell_type"].to_numpy()
    adata.obs["atlas_tissue_probability"] = tissue_transfer["prediction_probability"].to_numpy()
    adata.obs["atlas_tissue_uncertainty"] = tissue_transfer["uncertainty"].to_numpy()
    adata.obs["atlas_tissue_confidence"] = np.where(
        tissue_transfer["uncertainty"].to_numpy() <= 0.25,
        "High",
        np.where(tissue_transfer["resolved"].to_numpy(), "Medium", "Low"),
    )
    adata.obs["cell_type_ontology_id"] = [
        ontology_map.get(label, "") if label != "Unknown" else ""
        for label in adata.obs["cell_type"].astype(str)
    ]

    # Compatibility columns retained with paper-native meanings.
    adata.obs["cell_type_agreement"] = transfer["prediction_probability"].to_numpy()
    adata.obs["pred_weighted_knn"] = transfer["predicted_label"].to_numpy()
    adata.obs["prob_weighted_knn"] = transfer["prediction_probability"].to_numpy()

    cell_cols = [
        "cell_type_auto", "cell_type", "cell_type_ontology_id",
        "cell_type_probability", "cell_type_uncertainty",
        "cell_type_confidence", "cell_type_agreement",
        "pred_weighted_knn", "prob_weighted_knn",
        "cell_type_broad_atlas", "cell_type_broad_probability",
        "cell_type_broad_uncertainty", "atlas_cell_type_cell",
        "atlas_probability_cell", "atlas_uncertainty_cell",
    ]
    cell = adata.obs[cell_cols].copy()

    rows = []
    if cluster_key in adata.obs:
        for cluster, idx in adata.obs.groupby(cluster_key, observed=True).groups.items():
            x = adata.obs.loc[idx]
            resolved = x[x["cell_type_broad_atlas"].astype(str).isin(BROAD_LINEAGE_MARKERS)]
            resolved_fraction = float(len(resolved) / len(x))

            if resolved.empty:
                cluster_label = "Unknown"
                cluster_fraction = 0.0
            else:
                vc = resolved["cell_type_broad_atlas"].astype(str).value_counts()
                cluster_label = str(vc.index[0])
                cluster_fraction = float(vc.iloc[0] / len(resolved))

            rows.append({
                cluster_key: str(cluster),
                "atlas_broad_majority": cluster_label,
                "atlas_broad_fraction_among_resolved": cluster_fraction,
                "atlas_broad_resolved_fraction": resolved_fraction,
                "atlas_broad_support": resolved_fraction * cluster_fraction,
                "atlas_median_probability": float(x["cell_type_probability"].median()),
                "atlas_median_uncertainty": float(x["cell_type_uncertainty"].median()),
                "atlas_high_conf_fraction": float((x["cell_type_confidence"] == "High").mean()),
                "atlas_broad_unknown_fraction": float(
                    (~x["cell_type_broad_atlas"].astype(str).isin(BROAD_LINEAGE_MARKERS)).mean()
                ),
                "n_cells": int(len(x)),
            })
    cluster = pd.DataFrame(rows).set_index(cluster_key) if rows else pd.DataFrame()

    adata.uns["census_consensus_annotation"] = {
        "module_version": MODULE_VERSION,
        "method": "scArches_hierarchical_reference_mapping_weighted_kNN",
        "paper": "Lotfollahi et al., Nature Biotechnology 2022",
        "census_version": census_version,
        "model_uri": SCVI_MODEL_URI,
        "reference_cache": str(reference_cache),
        "reference_size": int(len(Xr)),
        "candidate_cells": int(candidate_cells),
        "max_cells_per_type": int(max_cells_per_type),
        "max_reference_cells": int(max_reference_cells),
        "weighted_knn_k": int(weighted_knn_k),
        "uncertainty_threshold": float(uncertainty_threshold),
        "query_max_epochs": int(query_max_epochs),
        "seed": int(seed),
        "accelerator": "cpu",
        "deterministic_training": True,
        "numpy_version": np.__version__,
        "scvi_tools_version": _package_version("scvi-tools"),
        "torch_version": _package_version("torch"),
        "normal_only": bool(normal_only),
        "tissues": list(tissues or ()),
    }

    resolved_n = int((adata.obs["cell_type"] != "Unknown").sum())
    print(
        f"[Census/scArches raw evidence] transfer complete before broad-conditioned refinement: "
        f"{resolved_n:,}/{adata.n_obs:,} cells resolved | "
        f"{adata.n_obs - resolved_n:,} unknown",
        flush=True,
    )

    return ConsensusAnnotationResult(
        cell_annotation=cell,
        cluster_annotation=cluster,
        method_predictions=transfer,
        reference_size=int(len(Xr)),
        census_version=census_version,
        model_uri=SCVI_MODEL_URI,
        reference_tissue_summary=reference_tissue_summary,
    )



BROAD_LINEAGE_MARKERS = {
    "Epithelial cell": {
        "EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "KRT14",
        "CLDN3", "CLDN4", "CLDN7", "MUC1", "KRT17",
    },
    "Endothelial cell": {
        "PECAM1", "VWF", "KDR", "ESAM", "EMCN", "CLDN5",
        "RAMP2", "ENG", "PLVAP", "CAVIN2",
    },
    "B-lineage cell": {
        "CD79A", "CD79B", "MS4A1", "CD74", "CD37", "HLA-DRA",
        "CD22", "VPREB1", "VPREB3", "CD19",
    },
    "Myeloid cell": {
        "AIF1", "TYROBP", "FCER1G", "CTSS", "LST1", "SPI1",
        "C1QA", "C1QB", "C1QC", "CSF1R", "GMFG",
    },
    "T/NK cell": {
        "CD3D", "CD3E", "TRBC1", "TRBC2", "IL7R",
        "NKG7", "GNLY", "KLRD1", "CCL5",
    },
    "Neural cell": {
        "STMN2", "GAP43", "UCHL1", "TUBB3", "RBFOX3",
        "MAP2", "SNAP25", "SYT1", "SYP", "NEFL",
    },
    "Schwann/glial cell": {
        "S100B", "SOX10", "S100A1", "PLP1", "MPZ",
        "PMP22", "SLC1A3", "GFAP", "FABP7",
    },
    "Fibroblast/stromal cell": {
        "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM",
        "COL6A1", "COL6A2", "PDGFRA", "C7", "COL5A1",
    },
    "Hematopoietic progenitor": {
        "CD34", "SPINK2", "HLF", "MECOM", "GATA2", "PROM1",
        "KIT", "AVP", "CRHBP", "HOXA9", "MEIS1",
    },
    "Erythroid cell": {
        "HBB", "HBA1", "HBA2", "ALAS2", "GYPA", "KLF1",
        "AHSP", "EPOR", "SLC4A1",
    },
    "Megakaryocyte": {
        "PPBP", "PF4", "GP9", "ITGA2B", "GP1BA", "GP6",
        "NRGN", "TREML1",
    },
}


def _atlas_broad_lineage(label: object) -> str:
    x = str(label or "")
    if x.lower() in ("", "unknown", "unresolved", "nan"):
        return "Unknown"
    return _broad_lineage_from_label(x)


def _score_marker_evidence(marker_rows: pd.DataFrame) -> dict[str, object]:
    """Apply the preregistered marker gates and score broad lineage panels."""
    required = {
        "names", "pvals_adj", "logfoldchanges",
        "pct_nz_group", "pct_nz_reference",
    }
    missing = sorted(required.difference(marker_rows.columns))
    if missing:
        raise KeyError(f"marker rows missing evidence columns: {missing}")

    eligible = marker_rows.loc[
        pd.to_numeric(marker_rows["pvals_adj"], errors="coerce").le(0.05)
        & pd.to_numeric(marker_rows["logfoldchanges"], errors="coerce").ge(1.0)
        & (
            pd.to_numeric(marker_rows["pct_nz_group"], errors="coerce")
            - pd.to_numeric(marker_rows["pct_nz_reference"], errors="coerce")
        ).ge(0.10)
    ].copy()
    strict_upper = list(dict.fromkeys(eligible["names"].astype(str).str.upper()))
    strict_scores = {
        lineage: sum(g in markers for g in strict_upper)
        for lineage, markers in BROAD_LINEAGE_MARKERS.items()
    }
    strict_ranked = sorted(
        strict_scores.items(), key=lambda kv: (kv[1], kv[0]), reverse=True,
    )
    strict_lineage, strict_score = strict_ranked[0]
    strict_second_score = strict_ranked[1][1] if len(strict_ranked) > 1 else 0
    strict_margin = int(strict_score - strict_second_score)

    # A small cluster tested against the full transcriptome can have a coherent
    # lineage program but no individual genes surviving global multiple-testing
    # correction. Rescue only a large, unambiguous multi-gene program; the final
    # resolver additionally requires Atlas broad-majority agreement.
    program_eligible = marker_rows.loc[
        pd.to_numeric(marker_rows["logfoldchanges"], errors="coerce").ge(1.0)
        & (
            pd.to_numeric(marker_rows["pct_nz_group"], errors="coerce")
            - pd.to_numeric(marker_rows["pct_nz_reference"], errors="coerce")
        ).ge(0.10)
    ].copy()
    program_upper = list(dict.fromkeys(
        program_eligible["names"].astype(str).str.upper()
    ))
    program_scores = {
        lineage: sum(g in markers for g in program_upper)
        for lineage, markers in BROAD_LINEAGE_MARKERS.items()
    }
    program_ranked = sorted(
        program_scores.items(), key=lambda kv: (kv[1], kv[0]), reverse=True,
    )
    program_lineage, program_score = program_ranked[0]
    program_second_score = program_ranked[1][1] if len(program_ranked) > 1 else 0
    program_margin = int(program_score - program_second_score)

    scoring_lineage = strict_lineage
    best_score = strict_score
    second_score = strict_second_score
    margin = strict_margin
    deciding_upper = strict_upper
    if strict_score >= 3 and strict_margin >= 2:
        strength = "strong"
        best_lineage = scoring_lineage
    elif strict_score >= 2 and strict_margin >= 1:
        strength = "moderate"
        best_lineage = scoring_lineage
    elif program_score >= 5 and program_margin >= 3:
        strength = "coherent_program"
        best_lineage = program_lineage
        scoring_lineage = program_lineage
        best_score = program_score
        second_score = program_second_score
        margin = program_margin
        deciding_upper = program_upper
    else:
        strength = "weak"
        best_lineage = "Uncertain"
    hits = [g for g in deciding_upper if g in BROAD_LINEAGE_MARKERS[scoring_lineage]]
    return {
        "marker_lineage": best_lineage,
        "marker_score": int(best_score),
        "marker_second_score": int(second_score),
        "marker_score_margin": margin,
        "marker_strength": strength,
        "marker_hits": ", ".join(hits),
        "marker_eligible_genes": ", ".join(strict_upper),
        "marker_strict_score": int(strict_score),
        "marker_strict_hits": ", ".join(
            g for g in strict_upper if g in BROAD_LINEAGE_MARKERS[strict_lineage]
        ),
        "marker_program_score": int(program_score),
        "marker_program_margin": int(program_margin),
        "marker_program_hits": ", ".join(
            g for g in program_upper if g in BROAD_LINEAGE_MARKERS[program_lineage]
        ),
        "marker_program_eligible_genes": ", ".join(program_upper),
    }


def _resolve_cluster_final_label(
    atlas_lineage: str,
    atlas_broad_support: float,
    marker_lineage: str,
    marker_strength: str,
    strong_atlas_support: float = 0.50,
    weak_atlas_support: float = 0.20,
) -> tuple[str, str, str, bool]:
    """Resolve broad annotation conservatively from independent Atlas/marker evidence."""
    atlas_lineage = str(atlas_lineage)
    marker_lineage = str(marker_lineage)
    atlas_valid = atlas_lineage not in {"Unknown", "Other", "Uncertain"}
    marker_valid = marker_lineage not in {"Unknown", "Other", "Uncertain"}
    atlas_strong = atlas_valid and atlas_broad_support >= float(strong_atlas_support)
    atlas_weak = atlas_valid and atlas_broad_support >= float(weak_atlas_support)
    marker_conflict = (
        marker_valid and atlas_valid and marker_lineage != atlas_lineage
    )

    if atlas_strong and marker_strength in {"strong", "coherent_program"} and marker_conflict:
        return "Uncertain", "Low", "high_support_conflict", True
    if marker_strength == "strong" and marker_valid:
        if atlas_valid and not marker_conflict:
            confidence = "High" if atlas_strong else "Medium"
            return marker_lineage, confidence, "strong_marker_atlas_agree", False
        return marker_lineage, "Medium", "strong_marker_low_atlas", marker_conflict
    if marker_strength == "coherent_program" and marker_valid:
        if atlas_valid and not marker_conflict:
            return (
                marker_lineage, "Medium",
                "coherent_marker_program_atlas_agree", False,
            )
        return (
            "Uncertain", "Low",
            "coherent_marker_program_requires_atlas_agreement", marker_conflict,
        )
    if marker_strength == "moderate" and marker_valid and atlas_weak and not marker_conflict:
        return marker_lineage, "Medium", "moderate_marker_atlas_agree", False
    if atlas_strong:
        return atlas_lineage, "Medium", "atlas_broad_supported", marker_conflict
    return "Uncertain", "Low", "insufficient_consensus", marker_conflict


def make_cluster_review(
    adata: anndata.AnnData,
    markers_top30: pd.DataFrame,
    cluster_key: str = "leiden",
    marker_cluster_col: str = "leiden",
    marker_gene_col: str = "names",
    n_markers: int = 30,
    strong_atlas_support: float = 0.50,
    weak_atlas_support: float = 0.20,
    fine_probability_threshold: float = 0.50,
    annotation_result: ConsensusAnnotationResult | None = None,
) -> pd.DataFrame:
    evidence_columns = {
        marker_cluster_col, marker_gene_col, "pvals_adj", "logfoldchanges",
        "pct_nz_group", "pct_nz_reference",
    }
    missing = sorted(evidence_columns.difference(markers_top30.columns))
    if missing:
        raise KeyError(
            f"markers_top30 missing required columns {missing}; "
            f"available columns: {list(markers_top30.columns)}"
        )

    marker_groups = {
        str(cluster): frame.head(int(n_markers)).rename(columns={marker_gene_col: "names"})
        for cluster, frame in markers_top30.groupby(marker_cluster_col, observed=True, sort=False)
    }

    rows = []
    final_by_cluster: dict[str, str] = {}

    for cluster, idx in adata.obs.groupby(cluster_key, observed=True).groups.items():
        x = adata.obs.loc[idx]
        resolved = x[x["cell_type_broad_atlas"].astype(str).isin(BROAD_LINEAGE_MARKERS)]
        resolved_fraction = float(len(resolved) / len(x))

        if resolved.empty:
            atlas_lineage = "Unknown"
            cluster_fraction = 0.0
        else:
            vc = resolved["cell_type_broad_atlas"].astype(str).value_counts()
            atlas_lineage = str(vc.index[0])
            cluster_fraction = float(vc.iloc[0] / len(resolved))

        atlas_broad_support = resolved_fraction * cluster_fraction
        marker_frame = marker_groups.get(str(cluster), pd.DataFrame(columns=list(evidence_columns)))
        marker_evidence = _score_marker_evidence(marker_frame)

        final_label, final_confidence, final_reason, marker_conflict = _resolve_cluster_final_label(
            atlas_lineage=atlas_lineage,
            atlas_broad_support=atlas_broad_support,
            marker_lineage=str(marker_evidence["marker_lineage"]),
            marker_strength=str(marker_evidence["marker_strength"]),
            strong_atlas_support=strong_atlas_support,
            weak_atlas_support=weak_atlas_support,
        )
        final_by_cluster[str(cluster)] = final_label

        rows.append({
            cluster_key: str(cluster),
            "atlas_broad_majority": atlas_lineage,
            "atlas_broad_fraction_among_resolved": cluster_fraction,
            "atlas_broad_resolved_fraction": resolved_fraction,
            "atlas_broad_support": atlas_broad_support,
            **marker_evidence,
            "marker_conflict": bool(marker_conflict),
            "cell_type_final": final_label,
            "cell_type_final_confidence": final_confidence,
            "cell_type_final_reason": final_reason,
            "atlas_median_probability": float(x["cell_type_probability"].median()),
            "atlas_median_uncertainty": float(x["cell_type_uncertainty"].median()),
            "atlas_unknown_fraction": float((x["cell_type"] == "Unknown").mean()),
            "top_markers": ", ".join(marker_frame["names"].astype(str)),
            "n_cells": int(len(x)),
        })

    review = pd.DataFrame(rows).set_index(cluster_key).sort_index()

    cluster_as_str = adata.obs[cluster_key].astype(str)
    adata.obs["cell_type_final"] = cluster_as_str.map(final_by_cluster).astype(str)
    adata.obs["cell_type_final_confidence"] = cluster_as_str.map(
        review["cell_type_final_confidence"].to_dict()
    ).astype(str)
    adata.obs["cell_type_final_reason"] = cluster_as_str.map(
        review["cell_type_final_reason"].to_dict()
    ).astype(str)
    adata.obs["cell_type_final_marker_hits"] = cluster_as_str.map(
        review["marker_hits"].to_dict()
    ).astype(str)
    adata.obs["cell_type_final_atlas_broad_support"] = cluster_as_str.map(
        review["atlas_broad_support"].to_dict()
    ).astype(float)
    # Compatibility aliases now have the same broad semantics as cell_type_final.
    adata.obs["final_cluster_annotation"] = adata.obs["cell_type_final"].astype(str)
    adata.obs["final_cluster_annotation_source"] = adata.obs["cell_type_final_reason"].astype(str)
    adata.obs["cell_type_final_source"] = adata.obs["cell_type_final_reason"].astype(str)

    raw_fine = adata.obs["pred_weighted_knn"].astype(str)
    fine_broad = raw_fine.map(_atlas_broad_lineage).astype(str)
    fine_probability = pd.to_numeric(adata.obs["atlas_probability_cell"], errors="coerce")
    final_broad = adata.obs["cell_type_final"].astype(str)
    rejection = np.select(
        [
            adata.obs["atlas_cell_type_cell"].astype(str).eq("Unknown"),
            fine_probability.lt(float(fine_probability_threshold)),
            final_broad.eq("Uncertain"),
            fine_broad.ne(final_broad),
        ],
        ["atlas_unknown", "low_probability", "final_uncertain", "broad_mismatch"],
        default="accepted",
    )
    accepted = rejection == "accepted"
    adata.obs["atlas_fine_broad_lineage"] = fine_broad
    adata.obs["cell_type_fine_confident"] = np.where(accepted, raw_fine, "Unknown")
    adata.obs["cell_type_fine_rejection_reason"] = rejection

    # Store as an AnnData-native DataFrame rather than list[dict].
    # list[dict] becomes an object array during H5AD serialization and fails
    # when rows contain mixed strings/bools/floats/ints.
    review_for_uns = review.reset_index().copy()
    for col in review_for_uns.select_dtypes(include=["object"]).columns:
        review_for_uns[col] = review_for_uns[col].fillna("").astype(str)
    adata.uns["cell_type_cluster_review"] = review_for_uns
    adata.uns["cell_type_annotation_decision_contract"] = {
        "module_version": MODULE_VERSION,
        "primary_label": "cell_type_final",
        "atlas_fine_evidence": "atlas_cell_type_cell",
        "fine_confident_label": "cell_type_fine_confident",
        "marker_top_n": int(n_markers),
        "marker_adjusted_p_max": 0.05,
        "marker_logfc_min": 1.0,
        "marker_detection_difference_min": 0.10,
        "strong_marker_min_hits": 3,
        "strong_marker_min_margin": 2,
        "moderate_marker_min_hits": 2,
        "moderate_marker_min_margin": 1,
        "coherent_program_logfc_min": 1.0,
        "coherent_program_detection_difference_min": 0.10,
        "coherent_program_min_hits": 5,
        "coherent_program_min_margin": 3,
        "coherent_program_requires_atlas_majority_agreement": True,
        "strong_atlas_support": float(strong_atlas_support),
        "weak_atlas_support": float(weak_atlas_support),
        "fine_probability_threshold": float(fine_probability_threshold),
        "sample_or_tissue_prior_used": False,
    }
    if annotation_result is not None:
        annotation_result.cell_annotation = adata.obs[
            [
                "cell_type_final", "cell_type_final_confidence",
                "cell_type_final_reason", "cell_type_final_marker_hits",
                "cell_type_final_atlas_broad_support", "atlas_cell_type_cell",
                "atlas_probability_cell", "atlas_uncertainty_cell",
                "cell_type_fine_confident", "cell_type_fine_rejection_reason",
                "final_cluster_annotation",
            ]
        ].copy()
        annotation_result.cluster_annotation = review.copy()
    return review


REFERENCE_FINE_LABEL_EXCLUDE_EXACT = frozenset({
    "", "unknown", "unresolved", "cell", "native cell", "primary cultured cell",
    "malignant cell", "pluripotent stem cell", "supporting cell",
})
REFERENCE_FINE_LABEL_EXCLUDE_SUBSTRINGS = (
    "sensu nematoda", "sensu protostomia",
)


def _fine_reference_label_is_usable(label: object, broad_lineage: str) -> bool:
    text = str(label or "").strip()
    lower = text.lower()
    if lower == str(broad_lineage).strip().lower():
        return False
    if lower in REFERENCE_FINE_LABEL_EXCLUDE_EXACT:
        return False
    if any(term in lower for term in REFERENCE_FINE_LABEL_EXCLUDE_SUBSTRINGS):
        return False
    return _broad_lineage_from_label(text) == str(broad_lineage)


def _conditioned_detail_transfer(
    Xr: np.ndarray,
    yr: np.ndarray,
    ontology_ids: np.ndarray,
    Xq: np.ndarray,
    broad_lineage: str,
    parents: dict[str, tuple[str, ...]],
    ontology_names: dict[str, str],
    n_neighbors: int,
    fine_probability_threshold: float,
    ontology_probability_threshold: float,
    min_ontology_depth: int,
) -> pd.DataFrame:
    from sklearn.neighbors import NearestNeighbors

    k = max(1, min(int(n_neighbors), len(Xr)))
    knn = NearestNeighbors(
        n_neighbors=k, algorithm="brute", metric="euclidean", n_jobs=-1,
    ).fit(np.asarray(Xr, dtype=np.float32))
    distances, indices = knn.kneighbors(np.asarray(Xq, dtype=np.float32))
    std = np.maximum(np.std(distances, axis=1), np.finfo(np.float32).eps)
    scale = (2.0 / std) ** 2
    weights = np.exp(-distances / scale[:, None])
    weights /= np.maximum(weights.sum(axis=1, keepdims=True), np.finfo(np.float32).eps)

    ancestor_cache: dict[str, set[str]] = {}
    depth_cache: dict[str, int] = {}
    rows: list[dict[str, object]] = []
    for cell_no in range(len(Xq)):
        local_labels = yr[indices[cell_no]].astype(str)
        local_ontology = ontology_ids[indices[cell_no]].astype(str)
        fine_votes: dict[str, float] = {}
        ontology_votes: dict[str, float] = {}
        for label, ontology_id, weight in zip(
            local_labels, local_ontology, weights[cell_no], strict=True,
        ):
            fine_votes[label] = fine_votes.get(label, 0.0) + float(weight)
            oid = _normalize_cl_id(ontology_id)
            if not oid:
                continue
            for node in {oid, *_ancestors(oid, parents, ancestor_cache)}:
                ontology_votes[node] = ontology_votes.get(node, 0.0) + float(weight)

        fine_label, fine_probability = max(fine_votes.items(), key=lambda item: item[1])
        detail_label = str(broad_lineage)
        detail_level = "broad"
        detail_probability = 1.0
        detail_reason = "broad_backoff"
        if fine_probability >= float(fine_probability_threshold):
            detail_label = fine_label
            detail_level = "fine"
            detail_probability = float(fine_probability)
            detail_reason = "conditioned_fine_supported"
        else:
            candidates = []
            for ontology_id, probability in ontology_votes.items():
                name = ontology_names.get(ontology_id, "")
                if probability < float(ontology_probability_threshold):
                    continue
                if _depth(ontology_id, parents, depth_cache) < int(min_ontology_depth):
                    continue
                if not _fine_reference_label_is_usable(name, broad_lineage):
                    continue
                candidates.append((
                    _depth(ontology_id, parents, depth_cache),
                    float(probability), ontology_id, name,
                ))
            if candidates:
                _, detail_probability, _, detail_label = max(candidates)
                detail_level = "ontology"
                detail_reason = "conditioned_ontology_backoff"

        rows.append({
            "conditioned_fine_raw": fine_label,
            "conditioned_fine_probability": float(fine_probability),
            "conditioned_fine_uncertainty": float(1.0 - fine_probability),
            "cell_type_fine_confident": (
                fine_label
                if fine_probability >= float(fine_probability_threshold)
                else "Unknown"
            ),
            "cell_type_detail_final": detail_label,
            "cell_type_detail_level": detail_level,
            "cell_type_detail_probability": float(detail_probability),
            "cell_type_detail_reason": detail_reason,
        })
    return pd.DataFrame(rows)


def refine_fine_annotation_with_broad_conditioning(
    adata: anndata.AnnData,
    reference_cache: str | Path,
    ontology_path: str | Path,
    *,
    final_broad_key: str = "cell_type_final",
    latent_key: str = "X_census_scvi",
    n_neighbors: int = 50,
    fine_probability_threshold: float = 0.50,
    ontology_probability_threshold: float = 0.50,
    min_ontology_depth: int = 3,
    min_reference_cells_per_label: int = 100,
    min_reference_datasets_per_label: int = 2,
    min_reference_tissues_per_label: int = 2,
    annotation_result: ConsensusAnnotationResult | None = None,
) -> pd.DataFrame:
    """Refine only within a validated broad lineage, then back off hierarchically."""
    if final_broad_key not in adata.obs:
        raise KeyError(f"Missing formal broad annotation: {final_broad_key}")
    if latent_key not in adata.obsm:
        raise KeyError(f"Missing Census query latent representation: {latent_key}")

    with np.load(reference_cache, allow_pickle=False) as ref:
        Xr = ref["embedding"].astype(np.float32)
        yr = ref["cell_type"].astype(str)
        ontology_ids = ref["ontology"].astype(str)
        reference_tissues = ref["tissue"].astype(str)
        reference_datasets = ref["dataset_id"].astype(str)
    parents, ontology_names = _parse_obo(ontology_path)
    reference_broad = np.asarray([_broad_lineage_from_label(label) for label in yr])
    reference_label_stats = (
        pd.DataFrame({
            "label": yr,
            "dataset": reference_datasets,
            "tissue": reference_tissues,
        })
        .groupby("label", observed=True)
        .agg(
            reference_cells=("label", "size"),
            reference_datasets=("dataset", "nunique"),
            reference_tissues=("tissue", "nunique"),
        )
    )
    reproducible_labels = set(reference_label_stats.index[
        reference_label_stats["reference_cells"].ge(int(min_reference_cells_per_label))
        & reference_label_stats["reference_datasets"].ge(int(min_reference_datasets_per_label))
        & reference_label_stats["reference_tissues"].ge(int(min_reference_tissues_per_label))
    ])

    output = pd.DataFrame(index=adata.obs_names)
    output["conditioned_fine_raw"] = "Unknown"
    output["conditioned_fine_probability"] = 0.0
    output["conditioned_fine_uncertainty"] = 1.0
    output["cell_type_fine_confident"] = "Unknown"
    output["cell_type_detail_final"] = adata.obs[final_broad_key].astype(str)
    output["cell_type_detail_level"] = "broad"
    output["cell_type_detail_probability"] = 1.0
    output["cell_type_detail_reason"] = "broad_backoff"

    reference_qc: dict[str, dict[str, int]] = {}
    final_broad = adata.obs[final_broad_key].astype(str)
    for broad_lineage in sorted(final_broad.unique()):
        if broad_lineage == "Uncertain":
            mask = final_broad.eq(broad_lineage)
            output.loc[mask, "cell_type_detail_final"] = "Uncertain"
            output.loc[mask, "cell_type_detail_level"] = "uncertain"
            output.loc[mask, "cell_type_detail_probability"] = 0.0
            output.loc[mask, "cell_type_detail_reason"] = "final_broad_uncertain"
            continue
        query_mask = final_broad.eq(broad_lineage).to_numpy()
        broad_mask = reference_broad == broad_lineage
        structurally_usable_mask = broad_mask & np.asarray([
            _fine_reference_label_is_usable(label, broad_lineage) for label in yr
        ])
        usable_mask = structurally_usable_mask & np.isin(yr, list(reproducible_labels))
        reference_qc[broad_lineage] = {
            "broad_reference_cells": int(broad_mask.sum()),
            "structurally_usable_reference_cells": int(structurally_usable_mask.sum()),
            "usable_fine_reference_cells": int(usable_mask.sum()),
            "excluded_reference_cells": int(broad_mask.sum() - usable_mask.sum()),
            "excluded_low_diversity_reference_cells": int(
                structurally_usable_mask.sum() - usable_mask.sum()
            ),
            "usable_fine_labels": int(pd.Series(yr[usable_mask]).nunique()),
        }
        if not query_mask.any() or usable_mask.sum() < max(25, int(n_neighbors)):
            continue
        local = _conditioned_detail_transfer(
            Xr[usable_mask], yr[usable_mask], ontology_ids[usable_mask],
            adata.obsm[latent_key][query_mask], broad_lineage,
            parents, ontology_names, n_neighbors,
            fine_probability_threshold, ontology_probability_threshold,
            min_ontology_depth,
        )
        output.loc[query_mask, local.columns] = local.to_numpy()

    for column in output.columns:
        adata.obs[column] = output[column].to_numpy()
    adata.obs["cell_type_fine_rejection_reason"] = np.where(
        adata.obs["cell_type_fine_confident"].astype(str).ne("Unknown"),
        "accepted_conditioned",
        np.where(
            final_broad.eq("Uncertain"),
            "final_broad_uncertain",
            "conditioned_low_probability",
        ),
    )
    adata.obs["cell_type_detail_confidence"] = np.select(
        [
            adata.obs["cell_type_detail_level"].eq("fine")
            & adata.obs["cell_type_detail_probability"].ge(0.75),
            adata.obs["cell_type_detail_level"].isin(["fine", "ontology"]),
        ],
        ["High", "Medium"],
        default="Low",
    )
    reference_qc_table = pd.DataFrame.from_dict(reference_qc, orient="index")
    reference_qc_table.index.name = "broad_lineage"
    adata.uns["conditioned_fine_annotation"] = {
        "method": "broad_conditioned_weighted_kNN_with_ontology_backoff",
        "reference_cache": str(reference_cache),
        # Store as a DataFrame: HDF5 treats '/' in dict keys (for example,
        # Schwann/glial cell) as a group separator and corrupts round-trips.
        "reference_label_qc": reference_qc_table,
        "fine_probability_threshold": float(fine_probability_threshold),
        "ontology_probability_threshold": float(ontology_probability_threshold),
        "min_ontology_depth": int(min_ontology_depth),
        "min_reference_cells_per_label": int(min_reference_cells_per_label),
        "min_reference_datasets_per_label": int(min_reference_datasets_per_label),
        "min_reference_tissues_per_label": int(min_reference_tissues_per_label),
        "sample_or_tissue_prior_used": False,
        "excluded_exact_labels": sorted(REFERENCE_FINE_LABEL_EXCLUDE_EXACT),
        "excluded_label_substrings": list(REFERENCE_FINE_LABEL_EXCLUDE_SUBSTRINGS),
    }
    if annotation_result is not None:
        annotation_result.cell_annotation = adata.obs[
            [
                "cell_type_final", "cell_type_final_confidence",
                "cell_type_final_reason", "cell_type_final_marker_hits",
                "cell_type_final_atlas_broad_support", "atlas_cell_type_cell",
                "atlas_probability_cell", "atlas_uncertainty_cell",
                "conditioned_fine_raw", "conditioned_fine_probability",
                "conditioned_fine_uncertainty", "cell_type_fine_confident",
                "cell_type_fine_rejection_reason", "cell_type_detail_final",
                "cell_type_detail_level", "cell_type_detail_probability",
                "cell_type_detail_confidence", "cell_type_detail_reason",
                "final_cluster_annotation",
            ]
        ].copy()
    return output



def annotate_immune_with_celltypist(
    adata: anndata.AnnData,
    model: str = "Immune_All_High.pkl",
    majority_voting: bool = True,
):
    import celltypist
    from celltypist import models
    models.download_models(model=model)
    source = adata.raw.to_adata() if adata.raw is not None else adata.copy()
    prediction = celltypist.annotate(source, model=model, majority_voting=majority_voting)
    obs = prediction.predicted_labels.copy()
    obs.columns = [f"immune_{c}" for c in obs.columns]
    for c in obs:
        adata.obs[c] = obs[c].reindex(adata.obs_names)
    return prediction


def export_consensus_annotation(
    output_dir: str | Path,
    result: ConsensusAnnotationResult,
    cluster_review: pd.DataFrame | None = None,
) -> dict[str, Path]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "cell_annotation": output_dir / "RNA_Census_consensus_cell_annotation.csv",
        "cluster_annotation": output_dir / "RNA_Census_consensus_by_Leiden_cluster.csv",
        "method_predictions": output_dir / "RNA_Census_method_predictions.csv",
    }
    result.cell_annotation.to_csv(paths["cell_annotation"])
    result.cluster_annotation.to_csv(paths["cluster_annotation"])
    result.method_predictions.to_csv(paths["method_predictions"])
    if result.reference_tissue_summary is not None and not result.reference_tissue_summary.empty:
        paths["reference_tissue_summary"] = (
            output_dir / "RNA_Census_reference_tissue_summary.csv"
        )
        result.reference_tissue_summary.to_csv(paths["reference_tissue_summary"])
    if cluster_review is not None:
        paths["cluster_review"] = output_dir / "RNA_Census_cluster_marker_review.csv"
        cluster_review.to_csv(paths["cluster_review"])
    return paths


__all__ = [
    "MODULE_VERSION", "CENSUS_VERSION", "SCVI_MODEL_URI", "CELL_ONTOLOGY_URL",
    "ConsensusAnnotationResult", "BROAD_LINEAGE_MARKERS", "download_scvi_model", "download_cell_ontology",
    "diagnose_census_gene_overlap", "project_query_scvi", "build_reference_cache",
    "annotate_from_census_consensus", "make_cluster_review",
    "refine_fine_annotation_with_broad_conditioning",
    "annotate_immune_with_celltypist", "export_consensus_annotation",
]
