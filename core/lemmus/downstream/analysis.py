"""QC后聚类、嵌入与物种注释编排；过滤细胞保留原始计数和缺失语义。"""
from __future__ import annotations

import hashlib
from pathlib import Path
from .config import (
    ANALYSIS_STAGE_VERSION,
    AnalysisError,
    DEFAULT_ANALYSIS_PARAMS,
    DEFAULT_CENSUS_PARAMS,
    RNA_QC_ANNOTATION_COLUMNS,
    validate_analysis_params,
)
from .integration import (
    integrate_dna_qc,
    load_dna_qc,
)
from .io import (
    write_lock,
    atomic_write_csv,
    atomic_write_bytes,
    invalidate_qc_completion,
    clear_analysis_outputs,
    atomic_write_h5ad,
    atomic_write_json,
    sanitize_for_h5ad,
)
from .qc import (
    _qc_snapshot_id,
    build_rna_qc_table,
    cell_qc_mask,
    cell_qc_metrics_table,
    compute_rna_qc,
    detect_species,
    load_rna_matrix,
    publish_rna_qc_table,
)
from ..paths import PIPELINE_ROOT


def sample_tissue_from_library(library: str) -> str:
    """Tissue token: text between the first and second '_' of the library name."""
    parts = str(library).split("_")
    return parts[1] if len(parts) > 2 else "NA"


def publish_full_annotated(adata_full, adata_analyzed, mask):
    """Complete-cell published H5AD: no data is lost at publication time.

    ``adata_full`` is the pre-filter object (raw counts + per-cell QC metrics
    for every cell); ``adata_analyzed`` holds the scanpy results computed on
    the QC-passing subset only (clustering/annotation must not be polluted by
    low-quality cells). The published object keeps all cells and raw counts,
    carries ``qc_pass`` for every cell, and attaches leiden/annotation
    columns plus PCA/UMAP coordinates only where they exist: filtered-out
    cells get empty annotations and all-NaN embedding rows (invisible in
    scatter plots, but recoverable for any re-analysis).
    """
    import numpy as np
    import pandas as pd

    published = adata_full.copy()
    published.obs["qc_pass"] = (
        mask.to_numpy().astype(bool) if isinstance(mask, pd.Series) else list(mask)
    )
    analyzed_names = adata_analyzed.obs_names.astype(str)
    extra = [c for c in adata_analyzed.obs.columns if c not in published.obs.columns]
    if extra:
        carried = adata_analyzed.obs[extra].reindex(published.obs_names.astype(str))
        for column in carried.columns:
            published.obs[column] = carried[column]
    position = published.obs_names.get_indexer(analyzed_names)
    if bool((position < 0).any()):
        raise AnalysisError(
            "analyzed cells are not a subset of the full cell set; cannot "
            "publish the complete-cell H5AD"
        )
    for key in list(adata_analyzed.obsm.keys()):
        embedding = np.asarray(adata_analyzed.obsm[key], dtype=float)
        full = np.full((published.n_obs, embedding.shape[1]), np.nan, dtype=float)
        full[position] = embedding
        published.obsm[key] = full
    for column in adata_analyzed.var.columns:
        if column not in published.var.columns:
            published.var[column] = adata_analyzed.var[column]
    published.uns.update(adata_analyzed.uns)
    return published


def normalize_log_hvg(adata, params: dict) -> None:
    """Full-gene invariant: no gene subsetting before projection; HVGs only
    mask the PCA computed later."""
    import scanpy as sc

    if "counts" not in adata.layers:
        raise AnalysisError("Missing raw counts layer before normalization.")
    if adata.layers["counts"].shape != adata.shape:
        raise AnalysisError(
            f"counts layer shape {adata.layers['counts'].shape} != adata shape {adata.shape}"
        )
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata.copy()
    sc.pp.highly_variable_genes(adata, n_top_genes=params["n_top_hvg"], flavor="seurat")
    gene_symbols = adata.var["gene_symbol"].astype(str).str.upper()
    adata.var["is_ig"] = gene_symbols.str.match(r"^IG[HKL]")
    adata.var["is_mt"] = gene_symbols.str.startswith("MT-")
    pca_mask = adata.var["highly_variable"].copy()
    pca_mask &= adata.var["n_cells_supported"].ge(params["min_cells_per_gene"])
    if params["exclude_ig_from_pca"]:
        pca_mask &= ~adata.var["is_ig"]
    if params["exclude_mt_from_pca"]:
        pca_mask &= ~adata.var["is_mt"]
    adata.var["highly_variable_for_pca"] = pca_mask
    if int(adata.var["highly_variable_for_pca"].sum()) < 800:
        raise AnalysisError(
            "排除 IG/MT 后用于 PCA 的 HVG <800；请检查基因命名或关闭 exclude_ig_from_pca 做诊断"
        )


def cluster_and_embed(adata, params: dict) -> None:
    import scanpy as sc

    seed = params["random_state"]
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(
        adata, mask_var="highly_variable_for_pca",
        n_comps=params["n_pcs"], svd_solver="arpack", random_state=seed,
    )
    sc.pp.neighbors(adata, n_neighbors=params["n_neighbors"], n_pcs=params["n_pcs"], random_state=seed)
    sc.tl.leiden(
        adata, resolution=params["leiden_resolution"], key_added="leiden",
        flavor="igraph", n_iterations=2, random_state=seed,
    )
    sc.tl.umap(adata, min_dist=params["umap_min_dist"], random_state=seed)


def compute_cluster_markers(adata):
    """Wilcoxon top-30 markers per Leiden cluster, symbols + stable gene IDs."""
    import pandas as pd
    import scanpy as sc

    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", use_raw=True, pts=True)
    marker_tables = []
    for cluster in adata.obs["leiden"].cat.categories:
        tmp = sc.get.rank_genes_groups_df(adata, group=cluster).head(30).copy()
        # rank_genes_groups returns var_names (versioned Ensembl IDs). Keep
        # those for traceability, and expose symbols for marker-panel matching.
        tmp.insert(0, "gene_id", tmp["names"].astype(str))
        symbol_by_id = adata.var["gene_symbol"].astype(str)
        tmp["names"] = tmp["gene_id"].map(symbol_by_id).fillna(tmp["gene_id"])
        tmp.insert(0, "leiden", str(cluster))
        marker_tables.append(tmp)
    return pd.concat(marker_tables, ignore_index=True)


def annotate_mouse(adata, markers_top30):
    """Deterministic mouse marker-program annotation (see lemmus.downstream.annotation_mouse)."""
    import pandas as pd

    from . import annotation_mouse as mouse_anno

    if getattr(mouse_anno, "MODULE_VERSION", None) != "v1":
        raise AnalysisError(
            f"Wrong lemmus.downstream.annotation_mouse module loaded: {getattr(mouse_anno, 'MODULE_VERSION', None)!r}"
        )
    result = mouse_anno.annotate_mouse_marker_program(adata, markers_top30=markers_top30, cluster_key="leiden")
    for column in result.cell_assignments.columns:
        adata.obs[column] = result.cell_assignments[column].astype(str)
    # Downstream DNAme plotting colors/groups on 'cell_type' for both species.
    adata.obs["cell_type"] = adata.obs["cell_type_final"].astype(str).astype("category")
    adata.uns["mouse_marker_annotation"] = result.uns
    adata.obs["sample_tissue"] = adata.obs["library"].astype(str).map(sample_tissue_from_library)
    return result


def annotate_human(adata, markers_top30, census_params: dict, pipeline_root: Path):
    """CELLxGENE Census scVI consensus (human only; pretrained models are human)."""
    from . import annotation_human as cta

    EXPECTED_ANNOTATION_MODULE_VERSION = "v35"
    if getattr(cta, "MODULE_VERSION", None) != EXPECTED_ANNOTATION_MODULE_VERSION:
        raise AnalysisError(
            f"Wrong lemmus.downstream.annotation_human module loaded: version={getattr(cta, 'MODULE_VERSION', None)!r}; "
            f"expected {EXPECTED_ANNOTATION_MODULE_VERSION}"
        )
    model_dir = pipeline_root / "resources" / "hg38_reference" / f"census_{census_params['census_version']}_scvi"
    if not (model_dir / "model.pt").is_file():
        raise AnalysisError(
            f"缺少 {model_dir / 'model.pt'}\n"
            f"aws s3 cp --no_sign-request {cta.SCVI_MODEL_URI} {model_dir / 'model.pt'}"
        )
    if cta.CENSUS_VERSION != census_params["census_version"]:
        raise AnalysisError(f"Module Census version {cta.CENSUS_VERSION} != {census_params['census_version']}")

    census_result = cta.annotate_from_census_consensus(
        adata,
        model_dir=model_dir,
        reference_cache=model_dir / f"census_{census_params['census_version']}_scvi_reference.npz",
        ontology_path=pipeline_root / "resources" / "hg38_reference" / "cell_ontology" / "cl-basic.obo",
        cluster_key="leiden",
        census_version=census_params["census_version"],
        candidate_cells=census_params["candidate_cells"],
        max_cells_per_type=census_params["max_cells_per_type"],
        max_reference_cells=census_params["max_reference_cells"],
        metadata_block_size=census_params["metadata_block_size"],
        embedding_batch_size=census_params["embedding_batch_size"],
        normal_only=census_params["normal_only"],
        tissues=tuple(census_params["tissues"]) if census_params["tissues"] else None,
        weighted_knn_k=census_params["weighted_knn_k"],
        uncertainty_threshold=census_params["uncertainty_threshold"],
        query_max_epochs=census_params["query_max_epochs"],
        seed=DEFAULT_ANALYSIS_PARAMS["random_state"],
    )
    cluster_review = cta.make_cluster_review(
        adata, markers_top30, annotation_result=census_result,
        cluster_key="leiden", marker_cluster_col="leiden", marker_gene_col="names", n_markers=30,
    )
    cta.refine_fine_annotation_with_broad_conditioning(
        adata,
        reference_cache=model_dir / f"census_{census_params['census_version']}_scvi_reference.npz",
        ontology_path=pipeline_root / "resources" / "hg38_reference" / "cell_ontology" / "cl-basic.obo",
        final_broad_key="cell_type_final",
        n_neighbors=census_params["weighted_knn_k"],
        fine_probability_threshold=0.50,
        ontology_probability_threshold=0.50,
        min_ontology_depth=3,
        min_reference_cells_per_label=100,
        min_reference_datasets_per_label=2,
        min_reference_tissues_per_label=2,
        annotation_result=census_result,
    )
    adata.obs["sample_tissue"] = adata.obs["library"].astype(str).map(sample_tissue_from_library)
    adata.obs["cell_type"] = adata.obs["cell_type_final"].astype(str).astype("category")
    return census_result, cluster_review


def _analysis_qc_tables(adata, stage: str):
    return adata.obs.groupby("library", observed=True).agg(
        n_cells=("library", "size"),
        median_UMI=("total_counts", "median"),
        median_supported_genes=("n_genes_supported", "median"),
        median_mito_pct=("pct_counts_mt", "median"),
    )


def _expression_support_audit(adata):
    return adata.obs.groupby("library", observed=True).agg(
        n_cells=("library", "size"),
        median_expected_umi=("total_counts", "median"),
        median_supported_expected_umi=("total_supported_expected_counts", "median"),
        median_subunit_expected_umi=("subunit_expected_counts", "median"),
        median_positive_genes=("n_genes_positive_expected", "median"),
        median_supported_genes=("n_genes_supported", "median"),
        median_fractional_only_genes=("n_genes_fractional_only", "median"),
        median_fractional_nnz_fraction=("fractional_nnz_fraction", "median"),
        median_subunit_count_mass_fraction=("subunit_count_mass_fraction", "median"),
    )


def run_analysis(
    project_dir: Path,
    *,
    params: dict | None = None,
    dna_project_root: Path | None = None,
    enable_immune_optional: bool = False,
    output_root: Path | None = None,
) -> Path:
    """从原始RNA计数计算QC、聚类和注释；固定路径覆盖，完成摘要最后提交。"""
    import importlib.metadata
    from datetime import datetime, timezone
    from time import perf_counter

    started = perf_counter()

    import pandas as pd

    project_dir = Path(project_dir).expanduser().resolve()
    pipeline_root = PIPELINE_ROOT
    params = validate_analysis_params(params)

    species = detect_species(project_dir)
    print(f"INFO: RNA 分析开始；物种：{species}；项目：{project_dir}", flush=True)

    # MT比例从原始生物gene计数计算；先校验所有输入，再替换正式QC。
    full, inputs = build_rna_qc_table(project_dir)
    snapshot_id = _qc_snapshot_id(project_dir)
    out_dir = (output_root or project_dir / "06_downstream").expanduser()
    if out_dir.is_symlink():
        raise AnalysisError(f"QC output directory must not be a symlink: {out_dir}")
    out_dir = out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    with write_lock(out_dir):
        adata = load_rna_matrix(project_dir)
        compute_rna_qc(adata, params)
        _mt = adata.obs["pct_counts_mt"].reindex(full["Sample_ID"].astype(str))
        if bool(_mt.isna().any()):
            missing = _mt[_mt.isna()].index.tolist()[:10]
            raise AnalysisError(
                "RNA QC table cells absent from the expression matrix; cannot "
                "attach pct_counts_mt: " + ", ".join(missing)
            )
        full["pct_counts_mt"] = _mt.to_numpy()
        integrated = None
        if dna_project_root is not None:
            dna_qc, dna_manifest = load_dna_qc(Path(dna_project_root))
            integrated = integrate_dna_qc(
                publish_rna_qc_table(full), dna_qc, dna_manifest, rna_identity=full
            )
            identity_columns = ["Sample_ID", "sample_id", "rna_sample", "PlateID", "Cell_Order", "RNA_Barcode"]
            identity_text = full[identity_columns].sort_values("Sample_ID").to_csv(index=False)
            inputs["dna"] = {**dna_qc.attrs["dna_provenance"], "library_links": integrated.attrs["library_links"],
                             "rna_identity_sha256": hashlib.sha256(identity_text.encode()).hexdigest(),
                             "library_map": integrated.attrs.get("library_map", {})}
        clear_analysis_outputs(out_dir, input_id=snapshot_id)
        atomic_write_csv(publish_rna_qc_table(full), out_dir / "RNA_QC_Information.csv")
        manifest_src = project_dir / "00_config" / "sample_manifest.tsv"
        atomic_write_bytes(manifest_src.read_bytes(), out_dir / "sample_manifest.tsv")

        # 聚类只消费通过QC的子集；正式H5AD仍保留全部细胞和原始计数。
        stage_started = perf_counter()
        print("INFO: 聚类开始；QC过滤、标准化、降维和marker计算", flush=True)
        pre_qc = _analysis_qc_tables(adata, "pre")
        mask = cell_qc_mask(adata.obs, params)
        atomic_write_csv(cell_qc_metrics_table(adata.obs, mask), out_dir / "RNA_cell_qc_metrics.csv")
        adata_full = adata.copy()  # complete-cell snapshot (raw counts) for publication
        adata = adata[mask].copy()
        post_qc = _analysis_qc_tables(adata, "post")
        print(f"INFO: QC 通过：{adata.n_obs}/{int(pre_qc['n_cells'].sum())} 个细胞；原始基因数：{adata.n_vars:,}", flush=True)
        normalize_log_hvg(adata, params)
        cluster_and_embed(adata, params)
        markers_top30 = compute_cluster_markers(adata)
        atomic_write_csv(markers_top30, out_dir / "RNA_Leiden_top30_markers.csv")

        expression_support_audit = _expression_support_audit(adata)
        atomic_write_csv(expression_support_audit, out_dir / "RNA_expression_support_audit.csv", index=True)
        atomic_write_csv(pre_qc, out_dir / "RNA_pre_qc_by_library.csv", index=True)
        atomic_write_csv(post_qc, out_dir / "RNA_post_qc_by_library.csv", index=True)

        cluster_by_library = pd.crosstab(adata.obs["leiden"], adata.obs["library"])

        print(f"INFO: 聚类完成；耗时：{perf_counter() - stage_started:.2f}s", flush=True)
        stage_started = perf_counter()
        print(f"INFO: 细胞注释开始；物种分支：{species}", flush=True)
        if species == "mouse":
            mouse_result = annotate_mouse(adata, markers_top30)
            atomic_write_csv(mouse_result.cluster_review, out_dir / "RNA_mouse_cluster_marker_review.csv", index=True)
            broad_by_library = pd.crosstab(adata.obs["cell_type_final"], adata.obs["library"])
            detail_by_library = pd.crosstab(adata.obs["cell_type_detail_final"], adata.obs["library"])
            detail_level_by_library = pd.crosstab(adata.obs["cell_type_detail_level"], adata.obs["library"])
            annotated_h5ad = out_dir / "RNA_mouse_marker_annotated.h5ad"
            annotation_summary = {
                "method": f"mouse_marker_program_v1",
                "module_version": mouse_anno_version(),
            }
        else:
            census_params = {**DEFAULT_CENSUS_PARAMS}
            census_result, cluster_review = annotate_human(adata, markers_top30, census_params, pipeline_root)
            atomic_write_csv(cluster_review, out_dir / "RNA_cluster_annotation_review.csv", index=True)
            broad_by_library = pd.crosstab(adata.obs["cell_type_final"], adata.obs["library"])
            detail_by_library = pd.crosstab(adata.obs["cell_type_detail_final"], adata.obs["library"])
            detail_level_by_library = pd.crosstab(adata.obs["cell_type_detail_level"], adata.obs["library"])
            fine_rejection_by_library = pd.crosstab(adata.obs["cell_type_fine_rejection_reason"], adata.obs["library"])
            reference_label_qc = adata.uns["conditioned_fine_annotation"]["reference_label_qc"].copy()
            atomic_write_csv(fine_rejection_by_library, out_dir / "RNA_fine_label_rejection_by_library.csv", index=True)
            atomic_write_csv(reference_label_qc, out_dir / "RNA_conditioned_reference_label_qc.csv", index=True)
            atomic_write_csv(census_result.reference_tissue_summary, out_dir / "RNA_census_reference_tissue_summary.csv", index=True)
            if enable_immune_optional:
                _annotate_immune_optional(adata)
            annotated_h5ad = out_dir / "RNA_Census_consensus_annotated.h5ad"
            annotation_summary = {"method": "census_scvi_consensus", "census_version": census_params["census_version"]}

        print(f"INFO: 细胞注释完成；耗时：{perf_counter() - stage_started:.2f}s", flush=True)
        # 未通过QC的细胞保持空注释，不把缺失解释为阴性。
        annotations = adata.obs[RNA_QC_ANNOTATION_COLUMNS].copy()
        annotations.index = adata.obs_names.astype(str)
        atomic_write_csv(
            publish_rna_qc_table(full, annotations=annotations),
            out_dir / "RNA_QC_Information.csv",
        )

        atomic_write_csv(cluster_by_library, out_dir / "RNA_cluster_by_library.csv", index=True)
        atomic_write_csv(broad_by_library, out_dir / "RNA_broad_annotation_by_library.csv", index=True)
        atomic_write_csv(detail_by_library, out_dir / "RNA_detail_annotation_by_library.csv", index=True)
        atomic_write_csv(detail_level_by_library, out_dir / "RNA_detail_level_by_library.csv", index=True)
        sample_tissue_summary = (
            adata.obs.groupby(["sample_tissue", "library"], observed=True)
            .size().rename("n_cells").reset_index().sort_values(["sample_tissue", "library"])
        )
        atomic_write_csv(sample_tissue_summary, out_dir / "RNA_sample_tissue_summary.csv")

        # Complete-cell publication: raw counts and qc_pass for every cell;
        # analysis results (leiden/annotations/embeddings) only for qc_pass cells.
        published = publish_full_annotated(adata_full, adata, mask)
        sanitize_for_h5ad(published)
        atomic_write_h5ad(published, annotated_h5ad)
        print(
            f"INFO: 注释矩阵已写入：{published.n_obs} 个细胞；"
            f"QC 通过：{int(published.obs['qc_pass'].sum())}；完成摘要待提交"
        )

        integrated_path = None
        if integrated is not None:
            updated = publish_rna_qc_table(full, annotations=annotations).set_index("Sample_ID")
            for column in RNA_QC_ANNOTATION_COLUMNS:
                integrated[column] = integrated["Sample_ID"].map(updated[column])
            integrated_path = out_dir / "Integrated_QC_Information.csv"
            if 'cell_links' in integrated.attrs:
                ledger = integrated.attrs['cell_links'].copy()
                ledger['rna_qc_pass'] = ledger['rna_cell_id'].map(published.obs['qc_pass'])
                atomic_write_csv(ledger, out_dir / 'DNA_RNA_cell_links.csv')
            atomic_write_csv(integrated, integrated_path)
            print(f"INFO: RNA×DNA QC 已写入：{integrated.shape[0]} 个细胞；{integrated_path}", flush=True)

        # 最后提交实际参数、输入身份及版本，前面的单文件成功不代表整组完成。
        versions = {}
        for package in ("scanpy", "anndata", "numpy", "pandas", "scipy"):
            try:
                versions[package] = importlib.metadata.version(package)
            except importlib.metadata.PackageNotFoundError:
                versions[package] = None
        if _qc_snapshot_id(project_dir) != snapshot_id:
            raise AnalysisError('RNA输入在分析期间改变；本次输出未完成，请重试')
        atomic_write_json(
            {
                "analysis_stage_version": ANALYSIS_STAGE_VERSION,
                "generated_at_utc": datetime.now(timezone.utc).isoformat(),
                "qc_snapshot_id": snapshot_id,
                "project_dir": str(project_dir),
                "species": species,
                "params": params,
                "random_state": params["random_state"],
                "annotation": annotation_summary,
                "package_versions": versions,
                "n_cells": int(adata.n_obs),
                "n_cells_total": int(published.n_obs),
                "n_genes": int(adata.n_vars),
                "n_leiden_clusters": int(adata.obs["leiden"].nunique()),
                "annotated_h5ad": annotated_h5ad.name,
                "integrated_table": integrated_path.name if integrated_path else None,
                "inputs": inputs,
            },
            out_dir / "analysis_summary.json",
        )
        print(f"INFO: RNA 分析完成；QC 通过：{adata.n_obs}/{published.n_obs} 个细胞；总耗时：{perf_counter() - started:.2f}s；输出：{out_dir}", flush=True)
        return out_dir


def mouse_anno_version() -> str:
    from . import annotation_mouse as mouse_anno

    return str(getattr(mouse_anno, "MODULE_VERSION", "unknown"))


def _annotate_immune_optional(adata) -> None:
    """Optional CellTypist immune refinement (human-trained models only)."""
    from . import annotation_human as cta

    cta.annotate_immune_with_celltypist(
        adata, model="Immune_All_High.pkl", majority_voting=True, download_model=True, prefix="immune_",
    )
