"""finalized RNA读取、细胞QC及命令分派；不改变原始分子计数。"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from datetime import datetime, timezone
from time import perf_counter
from pathlib import Path
import pandas as pd
from lemmus.config import select_config
from .config import (
    AnalysisError,
    DownstreamQCError,
    HUMAN_REFERENCES,
    IDENTITY_COLUMNS,
    MOUSE_REFERENCES,
    RNA_OBS_EXCLUDE,
    RNA_QC_ANNOTATION_COLUMNS,
    RNA_QC_COLUMNS,
    SCHEMA_VERSION,
    validate_analysis_params,
)
from ..io import (
    sha256_file,
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
    atomic_write_json,
)

DARLIN_CLONE_TABLE = Path("03_results/DARLIN/DARLIN_clone_id/DARLIN_Clone_ID.csv")


def _load_rna_h5ad_obs(h5ad_path: Path) -> pd.DataFrame:
    import anndata as ad

    adata = ad.read_h5ad(h5ad_path, backed="r")
    obs = adata.obs.copy()
    obs.insert(0, "Sample_ID", obs.index.astype(str))
    adata.file.close()
    return obs.reset_index(drop=True)


def _load_clone_table(clone_csv: Path) -> pd.DataFrame:
    clones = pd.read_csv(clone_csv, dtype=str, keep_default_na=False)
    if "Sample_ID" not in clones.columns:
        raise DownstreamQCError(f"DARLIN clone table has no Sample_ID column: {clone_csv}")
    duplicates = clones.loc[clones["Sample_ID"].duplicated(), "Sample_ID"].unique().tolist()
    if duplicates:
        raise DownstreamQCError(
            "DARLIN clone table Sample_ID is not unique: "
            + ", ".join(duplicates[:10])
        )
    return clones


def build_rna_qc_table(project_dir: Path) -> tuple[pd.DataFrame, dict]:
    """Merge per-sample RNA.h5ad obs with the DARLIN clone table (one row/cell)."""
    results_dir = project_dir / "03_results"
    count_dir = results_dir / "count"
    if not count_dir.is_dir():
        raise DownstreamQCError(f"Project has no 03_results/count directory: {project_dir}")

    h5ad_paths = sorted(count_dir.glob("*/RNA.h5ad"))
    if not h5ad_paths:
        raise DownstreamQCError(f"No 03_results/count/<sample_id>/RNA.h5ad found in {project_dir}")

    frames = []
    for h5ad_path in h5ad_paths:
        frame = _load_rna_h5ad_obs(h5ad_path)
        frame["Sample_ID"] = frame["Sample_ID"].astype(str)
        duplicates = frame.loc[frame["Sample_ID"].duplicated(), "Sample_ID"].unique().tolist()
        if duplicates:
            raise DownstreamQCError(
                f"RNA.h5ad obs_names are not unique in {h5ad_path}: "
                + ", ".join(duplicates[:10])
            )
        frames.append(frame)

    rna = pd.concat(frames, ignore_index=True) if len(frames) > 1 else frames[0]

    ordered_identity = [c for c in IDENTITY_COLUMNS if c in rna.columns]
    technical = [c for c in rna.columns if c not in RNA_OBS_EXCLUDE]
    rna = rna[ordered_identity + technical]

    clone_csv = project_dir / DARLIN_CLONE_TABLE
    if clone_csv.is_file():
        clones = _load_clone_table(clone_csv)
        rna_ids = set(rna["Sample_ID"])
        clone_ids = set(clones["Sample_ID"])
        darlin_only = sorted(clone_ids - rna_ids)
        if darlin_only:
            alignment = results_dir / "QC_Results" / "rna_darlin_cell_alignment.tsv"
            raise DownstreamQCError(
                "DARLIN clone table contains cells absent from RNA.h5ad: "
                f"{darlin_only[:10]}. See {alignment} for the pipeline-level audit."
            )
        rna = rna.merge(clones, on="Sample_ID", how="left", validate="one_to_one")
        cells_without_darlin = int(rna["DARLIN_Clone_ID"].isna().sum()) if "DARLIN_Clone_ID" in rna.columns else 0
        if cells_without_darlin:
            print(
                f"INFO: {cells_without_darlin}/{len(rna)} RNA cells have no DARLIN clone-table "
                "row (no DARLIN evidence); DARLIN columns are left empty for them."
            )
        darlin_inputs = {"clone_table": sha256_file(clone_csv)}
    else:
        print(
            f"WARN: no DARLIN clone table at {clone_csv}; "
            "publishing RNA-only QC table."
        )
        darlin_inputs = {}

    inputs = {"rna_h5ad": [str(p) for p in h5ad_paths], **darlin_inputs}
    return rna, inputs


def publish_rna_qc_table(
    rna: pd.DataFrame, annotations: pd.DataFrame | None = None
) -> pd.DataFrame:
    """Project the merged per-cell table onto the published CSV schema.

    Takes the full frame from :func:`build_rna_qc_table` and returns exactly
    ``RNA_QC_COLUMNS`` in order. ``CloneID`` is derived from ``sample_id``;
    ``DARLIN_Clone_ID`` is published as an empty column when the project has
    no DARLIN clone table, so the CSV schema is stable across projects. The
    four trailing annotation columns are filled from ``annotations`` (index =
    Sample_ID) when the analysis stage provides them, and stay empty
    otherwise; annotation rows pointing at cells absent from the QC table
    fail closed.
    """
    if "CloneID" not in rna.columns:
        if "sample_id" not in rna.columns:
            raise DownstreamQCError(
                "Cannot derive CloneID: RNA QC table has no sample_id column"
            )
        rna = rna.copy()
        rna["CloneID"] = rna["sample_id"].astype(str)
    if "DARLIN_Clone_ID" not in rna.columns:
        rna = rna.copy()
        rna["DARLIN_Clone_ID"] = pd.NA
    rna = rna.copy()
    if 'n_genes_supported' not in rna:
        raise DownstreamQCError('Finalized RNA results lack n_genes_supported')
    rna['num_genes_expressed'] = rna['n_genes_supported']
    if annotations is not None:
        missing_annotation = [
            column for column in RNA_QC_ANNOTATION_COLUMNS
            if column not in annotations.columns
        ]
        if missing_annotation:
            raise DownstreamQCError(
                "annotation frame lacks columns: " + ", ".join(missing_annotation)
            )
        unknown = sorted(
            set(annotations.index.astype(str)) - set(rna["Sample_ID"].astype(str))
        )
        if unknown:
            raise DownstreamQCError(
                "annotation rows reference cells absent from the RNA QC table: "
                + ", ".join(unknown[:10])
            )
        for column in RNA_QC_ANNOTATION_COLUMNS:
            values = pd.Series(
                annotations[column].to_numpy(),
                index=annotations.index.astype(str),
            )
            rna[column] = rna["Sample_ID"].astype(str).map(values)
    else:
        for column in RNA_QC_ANNOTATION_COLUMNS:
            rna[column] = pd.NA
    missing = [column for column in RNA_QC_COLUMNS if column not in rna.columns]
    if missing:
        raise DownstreamQCError(
            "RNA QC table lacks published-schema columns: " + ", ".join(missing)
        )
    return rna[RNA_QC_COLUMNS].reset_index(drop=True)


def _qc_snapshot_id(project_dir: Path) -> str:
    """Stable id binding the downstream table to the project inputs it was built from."""
    digest = hashlib.sha256()
    import yaml
    cfg = yaml.safe_load((project_dir / '00_config/config.yaml').read_text())
    files = ['00_config/config.yaml', '00_config/sample_manifest.tsv']
    files += [str(p.relative_to(project_dir)) for p in sorted((project_dir / '03_results/count').glob('*/RNA.h5ad'))]
    if cfg.get('pipeline_mode') != 'droplet':
        files += ['00_config/Barcode_Map.csv']
    if (project_dir / DARLIN_CLONE_TABLE).is_file():
        files.append(str(DARLIN_CLONE_TABLE))
    for rel in files:
        path = project_dir / rel
        if not path.is_file():
            raise DownstreamQCError(f"Missing project config file: {path}")
        digest.update(rel.encode("utf-8"))
        digest.update(b"\0")
        digest.update(sha256_file(path).encode("ascii"))
        digest.update(b"\0")
    return digest.hexdigest()


def _load_cached_qc(path: Path, rna: pd.DataFrame) -> pd.DataFrame:
    published = pd.read_csv(path, dtype={'Sample_ID':str, 'CloneID':str, 'PlateID':str},
                            keep_default_na=False, float_precision='round_trip')
    if (list(published.columns) != RNA_QC_COLUMNS
            or published['Sample_ID'].duplicated().any()
            or set(published['Sample_ID']) != set(rna['Sample_ID'])):
        raise DownstreamQCError('Existing RNA QC disagrees with finalized RNA results')
    current = rna.set_index('Sample_ID')
    aligned = published.set_index('Sample_ID').loc[current.index]
    expected_clones = current['CloneID'] if 'CloneID' in current else current['sample_id']
    if (aligned['CloneID'].tolist() != expected_clones.astype(str).tolist()
            or aligned['PlateID'].tolist() != current['PlateID'].astype(str).tolist()):
        raise DownstreamQCError('Existing RNA QC cell identities differ from finalized RNA results')
    return aligned.reset_index()


def run_processor(project_dir: Path, output_root: Path | None = None, *,
                  dna_project_root: Path | None = None,
                  params: dict | None = None) -> Path:
    started = perf_counter()
    project_dir = project_dir.expanduser().resolve()
    print(f"INFO: RNA QC 开始；项目：{project_dir}", flush=True)
    qc_results = (output_root or project_dir / "06_downstream").expanduser()
    if qc_results.is_symlink():
        raise DownstreamQCError(f"QC output directory must not be a symlink: {qc_results}")
    qc_results = qc_results.resolve()

    rna, inputs = build_rna_qc_table(project_dir)
    snapshot_id = _qc_snapshot_id(project_dir)
    out_dir = qc_results
    out_dir.mkdir(parents=True, exist_ok=True)

    with write_lock(out_dir):
        qc_path = out_dir / 'RNA_QC_Information.csv'
        is_droplet = 'protocol' in rna and rna['protocol'].eq('droplet').all()
        qc_keys = {'min_genes_per_cell', 'min_counts_per_cell', 'max_mito_pct'}
        if params and (not is_droplet or set(params) - qc_keys):
            raise DownstreamQCError('QC-only overrides require Droplet and cell QC thresholds; use --analyze for analysis parameters')
        metrics_path = out_dir/'RNA_cell_qc_metrics.csv'
        prior_params = None
        analysis_path = out_dir/'analysis_summary.json'
        summary_path = out_dir/'qc_summary.json'
        # 目录固定后，必须从完成摘要核对输入；相同细胞名不代表相同计数或clone表。
        prior = None
        if analysis_path.is_file():
            prior = json.loads(analysis_path.read_text())
        elif summary_path.is_file():
            prior = json.loads(summary_path.read_text())
        completed_qc = (qc_path.is_file() and isinstance(prior, dict)
                        and prior.get('qc_snapshot_id') == snapshot_id)
        if is_droplet and completed_qc:
            if analysis_path.is_file():
                prior_params = json.loads(analysis_path.read_text())['params']
            elif summary_path.is_file():
                prior_params = json.loads(summary_path.read_text()).get('inputs', {}).get('rna_quality_filter')
        effective_params = validate_analysis_params({**(prior_params or {}), **(params or {})})
        if completed_qc:
            published = _load_cached_qc(qc_path, rna)
            if is_droplet:
                metrics = pd.read_csv(metrics_path, dtype={'Sample_ID':str}, keep_default_na=False, float_precision='round_trip')
                if metrics['Sample_ID'].duplicated().any() or set(metrics['Sample_ID']) != set(rna['Sample_ID']):
                    raise DownstreamQCError('Existing RNA quality-filter table disagrees with finalized RNA results')
                mask = cell_qc_mask(metrics, effective_params)
                analysis_changed = (analysis_path.is_file()
                                    and any(effective_params[key] != prior_params[key] for key in qc_keys))
                if analysis_changed or not metrics['qc_pass'].equals(mask):
                    clear_analysis_outputs(out_dir, input_id=snapshot_id)
                    for column in RNA_QC_ANNOTATION_COLUMNS:
                        published[column] = pd.NA
                    metrics['qc_pass'] = mask
                    atomic_write_csv(metrics, metrics_path)
                    atomic_write_csv(published, qc_path)
                rna['qc_pass'] = rna['Sample_ID'].map(metrics.set_index('Sample_ID')['qc_pass'])
        elif is_droplet:
            adata = load_rna_matrix(project_dir)
            compute_rna_qc(adata, effective_params)
            mask = cell_qc_mask(adata.obs, effective_params)
            rna['pct_counts_mt'] = rna['Sample_ID'].map(adata.obs['pct_counts_mt'])
            rna['qc_pass'] = rna['Sample_ID'].map(pd.Series(mask, index=adata.obs_names))
            published = publish_rna_qc_table(rna)
            clear_analysis_outputs(out_dir, input_id=snapshot_id)
            atomic_write_csv(cell_qc_metrics_table(adata.obs, mask), metrics_path)
            atomic_write_csv(published, qc_path)
        else:
            # finalized obs保留原生计数，线粒体比例由下游从生物gene矩阵计算。
            adata = load_rna_matrix(project_dir)
            compute_rna_qc(adata, effective_params)
            rna['pct_counts_mt'] = rna['Sample_ID'].map(adata.obs['pct_counts_mt'])
            published = publish_rna_qc_table(rna)
            clear_analysis_outputs(out_dir, input_id=snapshot_id)
            atomic_write_csv(published, qc_path)
        if is_droplet:
            inputs['rna_quality_filter'] = effective_params

        if dna_project_root is not None:
            dna, manifest = load_dna_qc(dna_project_root)
            integrated = integrate_dna_qc(published, dna, manifest, rna_identity=rna)
            invalidate_qc_completion(out_dir)
            if 'cell_links' in integrated.attrs:
                atomic_write_csv(integrated.attrs['cell_links'], out_dir/'DNA_RNA_cell_links.csv')
            atomic_write_csv(integrated, out_dir/'Integrated_QC_Information.csv')
            provenance = {**dna.attrs['dna_provenance'],
                          'library_links':integrated.attrs['library_links'],
                          'library_map':integrated.attrs.get('library_map',{}),
                          'rna_qc_snapshot_id':snapshot_id, 'integrated_cells':len(integrated)}
            atomic_write_json(provenance, out_dir/'integration_summary.json')
            inputs['dna'] = provenance

        elif any((out_dir/name).exists() for name in (
                'integration_summary.json', 'Integrated_QC_Information.csv', 'DNA_RNA_cell_links.csv')):
            # 关闭可选整合时撤下旧整合，不触碰独立Droplet分支或用户文件。
            if (out_dir/'plots').is_symlink():
                raise DownstreamQCError(f'QC plots directory must not be a symlink: {out_dir}')
            invalidate_qc_completion(out_dir)
            for name in ('Integrated_QC_Information.csv', 'DNA_RNA_cell_links.csv'):
                (out_dir/name).unlink(missing_ok=True)
            for stem in ('09_dname_overlay_umap', '09_dname_qc_umaps', '09_dname_by_celltype'):
                for suffix in ('.pdf', '.png'):
                    (out_dir/'plots'/f'{stem}{suffix}').unlink(missing_ok=True)

        manifest_src = project_dir / "00_config" / "sample_manifest.tsv"
        manifest_copy = out_dir / "sample_manifest.tsv"
        atomic_write_bytes(manifest_src.read_bytes(), manifest_copy)

        if _qc_snapshot_id(project_dir) != snapshot_id:
            invalidate_qc_completion(out_dir)
            raise DownstreamQCError('RNA输入在QC期间改变；本次输出未完成，请重试')
        atomic_write_json(
            {
                "schema_version": SCHEMA_VERSION,
                "generated_at_utc": datetime.now(timezone.utc).isoformat(),
                "qc_snapshot_id": snapshot_id,
                "project_dir": str(project_dir),
                "n_cells": int(len(rna)),
                "n_samples": int(rna["sample_id"].nunique()) if "sample_id" in rna.columns else None,
                "columns": list(published.columns),
                "inputs": inputs,
            },
            out_dir / "qc_summary.json",
        )

        print(f"INFO: RNA QC 完成；细胞：{len(rna)}；耗时：{perf_counter() - started:.2f}s；输出：{qc_path}", flush=True)
        return out_dir


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m lemmus.downstream.qc", description=__doc__)
    parser.add_argument("--project-dir", required=True, type=Path, help="RNA_DARLIN project root")
    parser.add_argument(
        "--output-root", type=Path, default=None,
        help="Alternate output root (default: <project>/06_downstream)",
    )
    parser.add_argument(
        "--analyze", action="store_true",
        help="Run the full analysis stage (QC -> scanpy pipeline -> species "
             "annotation -> optional DNAme integration) in addition to the "
             "per-cell QC table.",
    )
    parser.add_argument(
        "--dna-project-root", type=Path, default=None,
        help="Optional Alopex project root; enables RNA x DNAme integration.",
    )
    parser.add_argument("--enable-immune-optional", action="store_true",
                        help="Human only: CellTypist immune refinement after the Census consensus.")
    parser.add_argument("--min-genes-per-cell", type=int, default=None)
    parser.add_argument("--min-counts-per-cell", type=int, default=None)
    parser.add_argument("--min-cells-per-gene", type=int, default=None)
    parser.add_argument("--max-mito-pct", type=float, default=None)
    parser.add_argument("--n-top-hvg", type=int, default=None)
    parser.add_argument("--n-pcs", type=int, default=None)
    parser.add_argument("--n-neighbors", type=int, default=None)
    parser.add_argument("--leiden-resolution", type=float, default=None)
    parser.add_argument("--umap-min-dist", type=float, default=None)
    parser.add_argument("--random-state", type=int, default=None)
    parser.add_argument("--keep-ig-in-pca", action="store_true",
                        help="Do not exclude IG genes from the PCA mask.")
    parser.add_argument("--keep-mt-in-pca", action="store_true",
                        help="Do not exclude mitochondrial genes from the PCA mask.")
    args = parser.parse_args(argv)
    try:
        if args.analyze:
            from .analysis import run_analysis

            overrides = {
                key: getattr(args, key)
                for key in (
                    "min_genes_per_cell", "min_counts_per_cell", "min_cells_per_gene",
                    "max_mito_pct", "n_top_hvg", "n_pcs", "n_neighbors",
                    "leiden_resolution", "umap_min_dist", "random_state",
                )
            }
            params = {key: value for key, value in overrides.items() if value is not None}
            if args.keep_ig_in_pca:
                params["exclude_ig_from_pca"] = False
            if args.keep_mt_in_pca:
                params["exclude_mt_from_pca"] = False
            run_analysis(
                args.project_dir,
                params=params or None,
                dna_project_root=args.dna_project_root,
                enable_immune_optional=args.enable_immune_optional,
                output_root=args.output_root,
            )
        else:
            analysis_only = ('min_cells_per_gene', 'n_top_hvg', 'n_pcs', 'n_neighbors',
                             'leiden_resolution', 'umap_min_dist', 'random_state')
            if any(getattr(args, key) is not None for key in analysis_only) or args.keep_ig_in_pca or args.keep_mt_in_pca or args.enable_immune_optional:
                raise DownstreamQCError('Analysis parameters require --analyze')
            qc_params = {key: getattr(args, key) for key in
                         ('min_genes_per_cell', 'min_counts_per_cell', 'max_mito_pct')
                         if getattr(args, key) is not None}
            run_processor(args.project_dir, args.output_root, dna_project_root=args.dna_project_root,
                          params=qc_params or None)
    except (DownstreamQCError, AnalysisError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    return 0


def detect_species(project_dir: Path) -> str:
    """Species branch from the project config `RNA.reference`."""
    import yaml

    config_path = Path(project_dir) / "00_config" / "config.yaml"
    if not config_path.is_file():
        raise AnalysisError(f"缺少项目配置: {config_path}")
    config = select_config(yaml.safe_load(config_path.read_text(encoding="utf-8")) or {})
    reference = str((config or {}).get("RNA", {}).get("reference", "")).strip().lower()
    if reference in HUMAN_REFERENCES:
        return "human"
    if reference in MOUSE_REFERENCES:
        return "mouse"
    raise AnalysisError(
        f"无法识别的 RNA reference: {reference!r}; "
        f"支持 human {sorted(HUMAN_REFERENCES)} 或 mouse {sorted(MOUSE_REFERENCES)}"
    )


def load_rna_matrix(project_dir: Path):
    """Load and concat every finalized per-sample RNA.h5ad (full-gene counts)."""
    import anndata
    import scipy.sparse as sp

    count_dir = Path(project_dir) / "03_results" / "count"
    paths = sorted(count_dir.glob("*/RNA.h5ad"))
    if not paths:
        raise AnalysisError(f"在 {count_dir} 中未找到 */RNA.h5ad")

    adatas = []
    seen: set[str] = set()
    required_obs = {"RNA_Barcode", "PlateID", "Cell_Order", "sample_id", "rna_sample"}
    for path in paths:
        sample_id = path.parent.name
        x = anndata.read_h5ad(path)
        x.obs_names = x.obs_names.astype(str)
        missing = required_obs.difference(x.obs.columns)
        if missing:
            raise AnalysisError(f"{path} 缺少 Pipeline finalization metadata: {sorted(missing)}")
        if not x.obs_names.is_unique:
            raise AnalysisError(f"{path} 中 project cell ID 不唯一")
        expected_prefix = f"{sample_id}_"
        if not all(str(cell).startswith(expected_prefix) for cell in x.obs_names):
            raise AnalysisError(f"{path} 的 cell ID 不符合 <sample_id>_<PlateID> contract")
        overlap = set(x.obs_names) & seen
        if overlap:
            raise AnalysisError(f"跨 RNA sample 出现重复 cell ID，例如: {sorted(overlap)[:5]}")
        seen.update(x.obs_names)
        observed_sample_ids = set(x.obs["sample_id"].astype(str))
        if observed_sample_ids != {sample_id}:
            raise AnalysisError(
                f"{path} 的 sample_id metadata 与结果目录不一致: {sorted(observed_sample_ids)}"
            )
        # library is the user-facing manifest sample_id; rna_sample stays as the
        # original FASTQ/library provenance identifier.
        x.obs["library"] = x.obs["sample_id"].astype(str)
        adatas.append(x)

    # Preserve gene_id/gene_symbol across samples; anndata.concat defaults
    # to merge=None and would silently discard all var metadata here.
    adata = anndata.concat(adatas, join="outer", merge="same", fill_value=0)
    required_var = {"gene_id", "gene_symbol"}
    missing_var = required_var.difference(adata.var.columns)
    if missing_var:
        raise AnalysisError(f"合并后缺少 gene metadata: {sorted(missing_var)}")
    gene_symbols = adata.var["gene_symbol"].astype("string")
    invalid_symbols = gene_symbols.isna() | gene_symbols.str.strip().eq("")
    if bool(invalid_symbols.any()):
        bad = adata.var_names[invalid_symbols][:5].tolist()
        raise AnalysisError(f"gene_symbol 为空，无法安全进行 QC/marker annotation，例如: {bad}")
    adata.X = sp.csr_matrix(adata.X)
    adata.layers["counts"] = adata.X.copy()
    return adata


def compute_rna_qc(adata, params: dict) -> None:
    """Support-aware QC metrics; var_names are versioned Ensembl IDs, so all
    biological gene matching goes through the retained gene_symbol column."""
    import numpy as np
    import scanpy as sc
    import scipy.sparse as sp

    gene_symbols = adata.var["gene_symbol"].astype(str).str.upper()
    adata.var["mt"] = gene_symbols.str.startswith("MT-")
    if int(adata.var["mt"].sum()) == 0:
        raise AnalysisError("未从 gene_symbol 识别到线粒体基因；请检查 reference annotation")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=[50], layer="counts", inplace=True)
    counts = sp.csr_matrix(adata.layers["counts"])
    support = float(params["gene_support_min_expected_umi"])
    supported = counts >= support
    adata.obs["n_genes_positive_expected"] = np.diff(counts.indptr)
    adata.obs["n_genes_supported"] = np.asarray(supported.sum(axis=1)).ravel().astype(int)
    adata.obs["n_genes_fractional_only"] = (
        adata.obs["n_genes_positive_expected"] - adata.obs["n_genes_supported"]
    )
    rows = np.repeat(np.arange(adata.n_obs), np.diff(counts.indptr))
    fractional = ~np.isclose(counts.data, np.rint(counts.data), atol=1e-6)
    subunit = (counts.data > 0) & (counts.data < support)
    fractional_nnz = np.bincount(rows, weights=fractional, minlength=adata.n_obs)
    subunit_mass = np.bincount(rows, weights=np.where(subunit, counts.data, 0.0), minlength=adata.n_obs)
    positive = adata.obs["n_genes_positive_expected"].to_numpy()
    total = adata.obs["total_counts"].to_numpy()
    adata.obs["fractional_nnz_fraction"] = np.divide(
        fractional_nnz, positive, out=np.zeros(adata.n_obs), where=positive > 0
    )
    adata.obs["subunit_count_mass_fraction"] = np.divide(
        subunit_mass, total, out=np.zeros(adata.n_obs), where=total > 0
    )
    adata.obs["subunit_expected_counts"] = subunit_mass
    adata.obs["total_supported_expected_counts"] = total - subunit_mass
    adata.var["n_cells_supported"] = np.asarray(supported.sum(axis=0)).ravel().astype(int)


def cell_qc_mask(obs, params: dict):
    """Pure-pandas cell filter so the threshold contract stays unit-testable."""
    import pandas as pd

    return (
        (obs["n_genes_supported"] >= params["min_genes_per_cell"])
        & (obs["total_counts"] >= params["min_counts_per_cell"])
        & (obs["pct_counts_mt"] < params["max_mito_pct"])
    ) if isinstance(obs, pd.DataFrame) else None


def cell_qc_metrics_table(obs, mask) -> "pd.DataFrame":
    """Per-cell gate metrics and verdict for ALL cells (pre-filter).

    Backs ``RNA_cell_qc_metrics.csv`` so the notebook violin plots can show
    the complete cell population (including cells the gate removed) together
    with the exact thresholds that removed them.
    """
    import pandas as pd

    verdict = mask.to_numpy().astype(bool) if isinstance(mask, pd.Series) else list(mask)
    return pd.DataFrame(
        {
            "Sample_ID": obs.index.astype(str),
            "library": obs["library"].astype(str).to_numpy(),
            "n_genes_supported": obs["n_genes_supported"].to_numpy(),
            "total_counts": obs["total_counts"].to_numpy(),
            "pct_counts_mt": obs["pct_counts_mt"].to_numpy(),
            "qc_pass": verdict,
        }
    )


if __name__ == "__main__":
    raise SystemExit(main())
