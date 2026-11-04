#!/usr/bin/env python3
"""Mouse marker-program annotation for RNA_DARLIN_Pipeline downstream.

Human projects are annotated by the CELLxGENE Census scVI consensus path in
``lemmus.downstream.annotation_human``. CELLxGENE publishes pretrained scVI embedding
models for human only, so mouse projects (mm10/GRCm38/GRCm39) instead get the
deterministic cluster-level marker-program annotation implemented here:

- Broad ``cell_type_final`` is decided once per Leiden cluster from curated
  mouse lineage marker programs (mean within-cluster detection fraction), with
  a margin/confidence rule and a cluster-marker gate that mirrors the human
  contract (adjusted P <= 0.05, logFC >= 1, detection difference >= 0.10).
- B-lineage and plasma-lineage clusters additionally receive a state-level
  ``cell_type_detail_final`` (Naive/GC/Memory B; Plasmablast/Plasma Cell).
- The annotation is fully deterministic (no randomness, no seeds); thresholds,
  marker programs, and per-cluster coverage are recorded for ``uns`` export.

Only ``numpy``/``pandas`` are imported at module level; ``anndata`` is imported
lazily inside :func:`annotate_mouse_marker_program` so unit tests can drive
the pure-DataFrame core without the atlas stack.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping, Sequence

import numpy as np
import pandas as pd

MODULE_VERSION = "v1"

# Broad lineage programs (mouse MGI symbols; matched case-insensitively so
# both MGI-style "Cd79a" and uppercase exports resolve to the same gene).
MOUSE_LINEAGE_PROGRAMS: Mapping[str, tuple[str, ...]] = {
    "B cell": ("Cd79a", "Cd79b", "Ms4a1", "Cd19", "Bank1", "Pax5", "Cd22"),
    "T cell": ("Cd3d", "Cd3e", "Cd3g", "Trac", "Cd2", "Cd28"),
    "NK cell": ("Nkg7", "Gzma", "Klrb1c", "Ncr1", "Prf1"),
    "Monocyte/Macrophage": ("Lyz2", "Itgam", "Csf1r", "Fcgr1", "C1qa", "C1qb", "Adgre1"),
    "Dendritic cell": ("Flt3", "H2-Ab1", "Itgax", "Xcr1", "Fcer1a"),
    "Neutrophil": ("S100a8", "S100a9", "Lst1", "Csf3r"),
    "Plasma cell": ("Jchain", "Xbp1", "Prdm1", "Sdc1", "Mzb1", "Derl3"),
    "Erythroid": ("Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Alas2"),
    "Epithelial cell": ("Epcam", "Krt8", "Krt18", "Cdh1"),
    "Fibroblast": ("Col1a1", "Col1a2", "Dcn", "Pdgfra"),
    "Endothelial cell": ("Pecam1", "Cdh5", "Vwf", "Kdr"),
}

# State programs resolved inside the B-cell and plasma-cell broad lineages.
MOUSE_B_STATE_PROGRAMS: Mapping[str, tuple[str, ...]] = {
    "Naive B": ("Cd79a", "Ms4a1", "Bank1", "Pax5", "Fcer2a"),
    "GC B": ("Aicda", "Bcl6", "Myc", "Chek1"),
    "Memory B": ("Cd27", "Tnfrsf13b", "Cd86"),
}
MOUSE_PLASMA_STATE_PROGRAMS: Mapping[str, tuple[str, ...]] = {
    "Plasmablast": ("Xbp1", "Prdm1", "Jchain", "Mki67", "Top2a"),
    "Plasma Cell": ("Sdc1", "Mzb1", "Derl3", "Prdm1", "Xbp1", "Jchain"),
}
STATE_PROGRAMS_BY_BROAD: Mapping[str, Mapping[str, tuple[str, ...]]] = {
    "B cell": MOUSE_B_STATE_PROGRAMS,
    "Plasma cell": MOUSE_PLASMA_STATE_PROGRAMS,
}

# Decision thresholds. These mirror the human annotation contract where a
# cluster-level decision is required before any cell-level label is written.
# Detection is support-aware: parsimony-gene-em emits fractional expected
# counts, so a gene counts as detected only at >= 1.0 expected UMI
# (log1p >= log(2)) — the same support rule as the pipeline's
# n_genes_supported metric. Plain expr > 0 would let fractional allocations
# inflate every program's detection.
DETECTION_LOG1P_MIN = float(np.log1p(1.0))
MIN_PROGRAM_GENE_COVERAGE = 0.50  # program is unscored when fewer genes map
BROAD_MIN_DETECTION = 0.35        # top lineage mean detection fraction (high confidence)
BROAD_MEDIUM_DETECTION = 0.25     # top lineage mean detection fraction (medium confidence)
BROAD_MIN_MARGIN = 0.10           # top - second lineage detection margin (high confidence)
BROAD_MEDIUM_MARGIN = 0.05        # medium confidence margin
STATE_MIN_MARGIN = 0.08           # B/plasma state programs
MARKER_PADJ_MAX = 0.05            # cluster marker gate: adjusted P
MARKER_LOGFC_MIN = 1.0            # cluster marker gate: logFC
MARKER_DETECTION_DIFF_MIN = 0.10  # cluster marker gate: detection difference

UNCERTAIN_LABEL = "Uncertain"
UNRESOLVED_LABEL = "Unresolved"


class MouseAnnotationError(RuntimeError):
    """Raised for fail-closed mouse annotation input problems."""


@dataclass
class MouseMarkerAnnotation:
    """Cluster-level annotation result mapped back to every cell."""

    cluster_review: pd.DataFrame
    cell_assignments: pd.DataFrame
    uns: dict = field(default_factory=dict)


def _uppercase_expression(expression: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """Uppercase gene columns, collapsing duplicate symbols deterministically.

    Duplicate uppercase symbols are collapsed onto the column with the highest
    mean expression; ties keep the leftmost column so the result is stable.
    """
    frame = expression.copy()
    frame.columns = pd.Index(frame.columns.astype(str).str.strip().str.upper())
    positions: dict[str, int] = {}
    means = frame.to_numpy(dtype=float).mean(axis=0)
    for index, label in enumerate(frame.columns):
        previous = positions.get(label)
        if previous is None or means[index] > means[previous]:
            positions[label] = index
    keep = sorted(positions.values())
    duplicates = int(frame.shape[1] - len(keep))
    if duplicates:
        frame = frame.iloc[:, keep]
    return frame, duplicates


def _program_detection_scores(
    cluster_block: pd.DataFrame,
    programs: Mapping[str, tuple[str, ...]],
) -> tuple[dict[str, float], dict[str, float], dict[str, list[str]]]:
    """Mean within-cluster support-aware detection fraction per program.

    A gene is detected at log1p expression >= log(2), i.e. >= 1.0 expected
    UMI on the linear scale; fractional EM allocations do not count.
    """
    scores: dict[str, float] = {}
    coverage: dict[str, float] = {}
    present_genes: dict[str, list[str]] = {}
    columns = set(cluster_block.columns)
    for label, genes in programs.items():
        mapped = [gene.upper() for gene in genes if gene.upper() in columns]
        coverage[label] = len(mapped) / len(genes)
        present_genes[label] = mapped
        if len(mapped) / len(genes) < MIN_PROGRAM_GENE_COVERAGE or not mapped:
            scores[label] = float("nan")
            continue
        detection = cluster_block[mapped].ge(DETECTION_LOG1P_MIN).mean(axis=0)
        scores[label] = float(detection.mean())
    return scores, coverage, present_genes


def _gated_marker_hits(
    program_genes: Sequence[str],
    cluster_markers: pd.DataFrame,
) -> list[str]:
    """Program genes that pass the cluster-level marker gate."""
    if cluster_markers is None or len(cluster_markers) == 0:
        return []
    columns = set(cluster_markers.columns)
    if not {"pvals_adj", "logfoldchanges"}.issubset(columns):
        return []
    gated = cluster_markers[
        (cluster_markers["pvals_adj"] <= MARKER_PADJ_MAX)
        & (cluster_markers["logfoldchanges"] >= MARKER_LOGFC_MIN)
    ].copy()
    if {"pct_nz_group", "pct_nz_reference"}.issubset(columns):
        gated = gated[
            (gated["pct_nz_group"] - gated["pct_nz_reference"])
            >= MARKER_DETECTION_DIFF_MIN
        ]
    gated_symbols = set(gated["names"].astype(str).str.strip().str.upper())
    return [gene for gene in program_genes if gene.upper() in gated_symbols]


def _assign_confidence(
    scores: Mapping[str, float], gated_counts: Mapping[str, int]
) -> tuple[str, str, float]:
    """Return (lineage, confidence, margin); alphabetical order breaks ties."""
    ranked = sorted(
        ((label, score) for label, score in scores.items() if np.isfinite(score)),
        key=lambda item: (-item[1], item[0]),
    )
    if not ranked:
        return UNCERTAIN_LABEL, "low", float("nan")
    top_label, top_score = ranked[0]
    margin = top_score - ranked[1][1] if len(ranked) > 1 else top_score
    if (
        top_score >= BROAD_MIN_DETECTION
        and margin >= BROAD_MIN_MARGIN
        and gated_counts.get(top_label, 0) >= 1
    ):
        return top_label, "high", float(margin)
    if top_score >= BROAD_MEDIUM_DETECTION and margin >= BROAD_MEDIUM_MARGIN:
        return top_label, "medium", float(margin)
    # Neither gate passes: fail closed to Uncertain instead of keeping the
    # weakly-supported argmax lineage.
    return UNCERTAIN_LABEL, "low", float(margin)


def annotate_mouse_marker_program_core(
    expression: pd.DataFrame,
    clusters: pd.Series,
    markers_top30: pd.DataFrame | None = None,
    marker_cluster_column: str = "leiden",
) -> MouseMarkerAnnotation:
    """Deterministic cluster-level marker-program annotation (pure pandas).

    ``expression`` must be cells x genes log-normalized values (detection only
    needs expr > 0). ``clusters`` is the per-cell cluster label (string-cast).
    ``markers_top30`` optionally carries per-cluster top markers with
    ``names`` (gene symbols), ``pvals_adj`` and ``logfoldchanges``.
    """
    if expression.shape[0] != len(clusters):
        raise MouseAnnotationError(
            f"expression rows ({expression.shape[0]}) != cluster labels ({len(clusters)})"
        )
    if expression.index.astype(str).tolist() != clusters.index.astype(str).tolist():
        raise MouseAnnotationError("expression and cluster labels are not aligned")

    frame, duplicates = _uppercase_expression(expression)
    cluster_labels = clusters.astype(str)

    marker_tables: dict[str, pd.DataFrame] = {}
    if markers_top30 is not None and len(markers_top30):
        required = {marker_cluster_column, "names", "pvals_adj", "logfoldchanges"}
        missing = required.difference(markers_top30.columns)
        if missing:
            raise MouseAnnotationError(
                f"markers_top30 lacks required columns: {sorted(missing)}"
            )
        marker_tables = {
            str(key): table
            for key, table in markers_top30.groupby(marker_cluster_column, observed=True)
        }

    rows: list[dict] = []
    assignments: list[pd.DataFrame] = []
    for cluster in sorted(cluster_labels.unique()):
        in_cluster = cluster_labels == cluster
        block = frame.loc[in_cluster]
        scores, coverage, present_genes = _program_detection_scores(
            block, MOUSE_LINEAGE_PROGRAMS
        )
        cluster_markers = marker_tables.get(cluster, pd.DataFrame())
        gated_counts = {
            label: len(_gated_marker_hits(genes, cluster_markers))
            for label, genes in MOUSE_LINEAGE_PROGRAMS.items()
        }
        lineage, confidence, margin = _assign_confidence(scores, gated_counts)

        detail, detail_level, detail_reason = lineage, "broad", ""
        state_programs = STATE_PROGRAMS_BY_BROAD.get(lineage)
        if state_programs is not None:
            state_scores, state_coverage, state_genes = _program_detection_scores(
                block, state_programs
            )
            ranked_states = sorted(
                (
                    (label, score)
                    for label, score in state_scores.items()
                    if np.isfinite(score)
                ),
                key=lambda item: (-item[1], item[0]),
            )
            if ranked_states:
                top_state, top_score = ranked_states[0]
                state_margin = (
                    top_score - ranked_states[1][1] if len(ranked_states) > 1 else top_score
                )
                if state_margin >= STATE_MIN_MARGIN:
                    detail, detail_level = top_state, "state"
                    detail_reason = (
                        f"{top_state} detection={top_score:.3f} margin={state_margin:.3f}"
                    )
                else:
                    detail_reason = (
                        f"state margin {state_margin:.3f} < {STATE_MIN_MARGIN}; "
                        "kept at broad level"
                    )
        if lineage == UNCERTAIN_LABEL:
            detail, detail_level, detail_reason = (
                UNRESOLVED_LABEL,
                "unresolved",
                "no lineage program reached the detection/margin gate",
            )

        reason = (
            f"{lineage} detection={scores.get(lineage, float('nan')):.3f} "
            f"margin={margin:.3f} gated={gated_counts.get(lineage, 0)}"
            if lineage != UNCERTAIN_LABEL
            else f"no lineage program passed the detection/margin gates (margin={margin:.3f})"
        )
        rows.append(
            {
                "cluster": cluster,
                "n_cells": int(in_cluster.sum()),
                "cell_type_final": lineage,
                "cell_type_final_confidence": confidence,
                "cell_type_final_reason": reason,
                "cell_type_detail_final": detail,
                "cell_type_detail_level": detail_level,
                "cell_type_detail_reason": detail_reason,
                **{f"detection_{label}": round(score, 4) for label, score in scores.items()},
                **{f"gated_{label}": count for label, count in gated_counts.items()},
                "program_gene_coverage": ";".join(
                    f"{label}:{coverage[label]:.2f}" for label in sorted(coverage)
                ),
            }
        )
        assignments.append(
            pd.DataFrame(
                {
                    "cell_type_final": [lineage] * int(in_cluster.sum()),
                    "cell_type_final_confidence": [confidence] * int(in_cluster.sum()),
                    "cell_type_detail_final": [detail] * int(in_cluster.sum()),
                    "cell_type_detail_level": [detail_level] * int(in_cluster.sum()),
                    "cell_type_detail_reason": [detail_reason] * int(in_cluster.sum()),
                    "annotation_method": [f"mouse_marker_program_{MODULE_VERSION}"] * int(in_cluster.sum()),
                },
                index=block.index,
            )
        )

    cluster_review = pd.DataFrame(rows).set_index("cluster", drop=True)
    cell_assignments = pd.concat(assignments).loc[expression.index]
    uns = {
        "annotation_method": f"mouse_marker_program_{MODULE_VERSION}",
        "module_version": MODULE_VERSION,
        "deterministic": True,
        "thresholds": {
            "detection_log1p_min": DETECTION_LOG1P_MIN,
            "min_program_gene_coverage": MIN_PROGRAM_GENE_COVERAGE,
            "broad_min_detection": BROAD_MIN_DETECTION,
            "broad_medium_detection": BROAD_MEDIUM_DETECTION,
            "broad_min_margin": BROAD_MIN_MARGIN,
            "broad_medium_margin": BROAD_MEDIUM_MARGIN,
            "state_min_margin": STATE_MIN_MARGIN,
            "marker_padj_max": MARKER_PADJ_MAX,
            "marker_logfc_min": MARKER_LOGFC_MIN,
            "marker_detection_diff_min": MARKER_DETECTION_DIFF_MIN,
        },
        "lineage_programs": {
            label: list(genes) for label, genes in MOUSE_LINEAGE_PROGRAMS.items()
        },
        "state_programs": {
            broad: {label: list(genes) for label, genes in programs.items()}
            for broad, programs in STATE_PROGRAMS_BY_BROAD.items()
        },
        "duplicate_uppercase_symbols_collapsed": duplicates,
    }
    return MouseMarkerAnnotation(
        cluster_review=cluster_review, cell_assignments=cell_assignments, uns=uns
    )


def annotate_mouse_marker_program(
    adata,
    markers_top30: pd.DataFrame | None = None,
    cluster_key: str = "leiden",
    gene_symbol_column: str = "gene_symbol",
) -> MouseMarkerAnnotation:
    """anndata wrapper: expression comes from ``adata.raw`` (log1p, full genes).

    The notebook scales ``adata.X`` for PCA, so the marker programs must be
    scored on the pre-scaling log-normalized matrix retained in ``raw``.
    """
    if cluster_key not in adata.obs.columns:
        raise MouseAnnotationError(f"adata.obs lacks cluster column: {cluster_key}")
    if gene_symbol_column not in adata.var.columns:
        raise MouseAnnotationError(f"adata.var lacks gene symbol column: {gene_symbol_column}")
    if adata.raw is not None:
        source, var_symbols = adata.raw.to_adata(), adata.raw.var[gene_symbol_column]
    else:
        source, var_symbols = adata, adata.var[gene_symbol_column]
    import scipy.sparse as sp

    matrix = source.X
    if sp.issparse(matrix):
        matrix = matrix.toarray()
    expression = pd.DataFrame(
        np.asarray(matrix, dtype=float),
        index=source.obs_names.astype(str),
        columns=var_symbols.astype(str),
    )
    return annotate_mouse_marker_program_core(
        expression,
        clusters=adata.obs[cluster_key],
        markers_top30=markers_top30,
        marker_cluster_column=cluster_key,
    )
