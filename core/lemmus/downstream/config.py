"""下游参数、发布schema和错误类型；无写事务。"""
from __future__ import annotations

import math


SCHEMA_VERSION = 1


IDENTITY_COLUMNS = [
    "Sample_ID",
    "sample_id",
    "PlateID",
    "Cell_Order",
    "RNA_Barcode",
    "rna_sample",
    "library",
]


RNA_OBS_EXCLUDE = set(IDENTITY_COLUMNS) | {"cell_id"}


RNA_QC_ANNOTATION_COLUMNS = [
    "cell_type_final",
    "cell_type_final_confidence",
    "cell_type_detail_final",
    "cell_type_detail_level",
]


RNA_QC_COLUMNS = [
    "Sample_ID",
    "CloneID",
    "PlateID",
    "mapping_rate",
    "dedup_rate",
    "num_genes_expressed",
    "num_genes_over_mean",
    "pct_counts_mt",
    "DARLIN_Clone_ID",
    *RNA_QC_ANNOTATION_COLUMNS,
]


class DownstreamQCError(RuntimeError):
    """Raised for any fail-closed downstream input/output problem."""


ANALYSIS_STAGE_VERSION = "v1"


HUMAN_REFERENCES = {"hg38", "grch38"}


MOUSE_REFERENCES = {"mm10", "grcm38", "grcm39", "mouse"}


DEFAULT_ANALYSIS_PARAMS = {
    "min_genes_per_cell": 300,
    "gene_support_min_expected_umi": 1.0,
    "min_counts_per_cell": 1000,
    "min_cells_per_gene": 3,
    "max_mito_pct": 10.0,
    "n_top_hvg": 2000,
    "n_pcs": 30,
    "n_neighbors": 15,
    "leiden_resolution": 0.7,
    "umap_min_dist": 0.45,
    "exclude_ig_from_pca": True,
    "exclude_mt_from_pca": True,
    "random_state": 42,
}


DEFAULT_CENSUS_PARAMS = {
    "census_version": "2025-11-08",
    "candidate_cells": 10_000_000,
    "max_cells_per_type": 3_000,
    "max_reference_cells": 500_000,
    "metadata_block_size": 500_000,
    "embedding_batch_size": 50_000,
    "normal_only": True,
    "tissues": (),
    "weighted_knn_k": 50,
    "uncertainty_threshold": 0.50,
    "query_max_epochs": 200,
}


class AnalysisError(RuntimeError):
    """Raised for fail-closed analysis-stage problems."""


def validate_analysis_params(params: dict) -> dict:
    merged = {**DEFAULT_ANALYSIS_PARAMS, **(params or {})}
    positive_ints = (
        "min_genes_per_cell", "min_counts_per_cell", "min_cells_per_gene",
        "n_top_hvg", "n_pcs", "n_neighbors", "random_state",
    )
    for key in positive_ints:
        value = merged[key]
        if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
            raise AnalysisError(f"analysis param {key} must be a positive int, got {value!r}")
    for key in ("max_mito_pct", "gene_support_min_expected_umi", "leiden_resolution", "umap_min_dist"):
        value = merged[key]
        if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(value) or value <= 0:
            raise AnalysisError(f"analysis param {key} must be a positive number, got {value!r}")
    if merged['max_mito_pct'] > 100:
        raise AnalysisError('analysis param max_mito_pct must not exceed 100')
    for key in ("exclude_ig_from_pca", "exclude_mt_from_pca"):
        if not isinstance(merged[key], bool):
            raise AnalysisError(f"analysis param {key} must be a bool, got {merged[key]!r}")
    return merged
