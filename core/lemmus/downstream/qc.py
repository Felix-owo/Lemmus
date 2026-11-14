"""finalized RNA读取、细胞QC及命令分派；不改变原始分子计数。"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import stat
import sys
from datetime import datetime, timezone
from time import perf_counter
from pathlib import Path
import pandas as pd
from lemmus.config import ConfigError, load_config
from lemmus.rna import validate_prepared_rna, validate_rna_count_matrix
from lemmus.paths import PIPELINE_ROOT, project_results_root
from lemmus.inputs import read_manifest, validate_barcode_map_provenance
from lemmus.io import read_delimited_frame
from .config import (
    AnalysisError,
    DownstreamQCError,
    DEFAULT_CENSUS_PARAMS,
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
    file_stat_identity,
    sha256_file,
    stable_file_identity,
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
    downstream_source_sha256,
    human_annotation_reference_identity,
)

DARLIN_CLONE_TABLE = Path("DARLIN/DARLIN_clone_id/DARLIN_Clone_ID.csv")


def _rna_h5ad_paths(project_dir: Path, config: dict) -> list[Path]:
    """枚举全部已完成文库；仅缺席的RNA.h5ad可跳过，权限与I/O错误必须传播。"""
    count_dir = project_results_root(project_dir,config)/"count"
    paths = []
    with os.scandir(count_dir) as entries:
        for entry in entries:
            if not stat.S_ISDIR(entry.stat().st_mode):
                continue
            path = Path(entry.path) / "RNA.h5ad"
            try:
                mode = path.stat().st_mode
            except FileNotFoundError:
                continue
            if not stat.S_ISREG(mode):
                raise DownstreamQCError(f"Finalized RNA result must be a regular file: {path}")
            paths.append(path)
    return sorted(paths)


def _darlin_clone_path(project_dir: Path, config: dict) -> Path | None:
    if not config.get('DARLIN',{}).get('enabled',True):
        return None
    path = project_results_root(project_dir,config)/DARLIN_CLONE_TABLE
    try:
        mode = path.stat().st_mode
    except FileNotFoundError:
        return None
    if not stat.S_ISREG(mode):
        raise DownstreamQCError(f'DARLIN clone table must be a regular file: {path}')
    return path


def _validate_rna_libraries(project_dir: Path, obs: pd.DataFrame) -> Path:
    """固定输出目录可能保留旧文库；原始与项目文库身份必须共同匹配当前manifest。"""
    config = load_config(Path(project_dir)/'00_config/config.yaml')
    manifest_path = Path(config.get('sample_manifest','00_config/sample_manifest.tsv')).expanduser()
    if not manifest_path.is_absolute():
        manifest_path = Path(project_dir)/manifest_path
    try:
        fields, rows = read_manifest(manifest_path)
    except ValueError as exc:
        raise DownstreamQCError(f'Invalid RNA manifest: {exc}') from exc
    manifest = pd.DataFrame(rows, columns=fields)
    columns = ['sample_id', 'rna_sample']
    if any(column not in frame for column in columns for frame in (manifest, obs)):
        raise DownstreamQCError('RNA manifest/results lack sample_id/rna_sample identity')
    selected = manifest.loc[manifest.rna_sample.ne(''), columns]
    observed = obs[columns].drop_duplicates()
    if (any(selected[column].duplicated().any() or observed[column].duplicated().any() for column in columns)
            or set(map(tuple, selected.to_numpy())) != set(map(tuple, observed.to_numpy()))):
        raise DownstreamQCError('Finalized RNA libraries disagree with current manifest sample_id/rna_sample; '
                                'rerun the Pipeline and resolve stale result directories')

    required = {'RNA_Barcode', 'PlateID', 'Cell_Order', 'library'}
    if not required.issubset(obs):
        raise DownstreamQCError(f'Finalized RNA cell identity lacks metadata: {sorted(required.difference(obs))}')
    if not obs['library'].astype(str).equals(obs['sample_id'].astype(str)):
        raise DownstreamQCError('Finalized RNA library identity differs from canonical sample_id')
    sample = obs['sample_id'].astype(str)
    barcode = obs['RNA_Barcode'].astype(str)
    actual_ids = obs['Sample_ID'].astype(str) if 'Sample_ID' in obs else obs.index.astype(str)
    if config.get('pipeline_mode') == 'droplet':
        if 'protocol' not in obs or not obs['protocol'].eq('droplet').all():
            raise DownstreamQCError('Finalized RNA protocol identity disagrees with Droplet config')
        if not obs['PlateID'].eq('').all() or not barcode.str.fullmatch('[ACGT]{17}').all():
            raise DownstreamQCError('Finalized Droplet cell identity requires an explicit CB17 barcode and empty PlateID')
        for _, group in obs.groupby('sample_id', observed=True, sort=False):
            orders = group['Cell_Order'].astype(str)
            if (not orders.str.fullmatch('[1-9][0-9]*').all()
                    or {int(value) for value in orders} != set(range(1,len(group)+1))):
                raise DownstreamQCError('Finalized Droplet Cell_Order identity must rank each called cell once')
            if 'total_expected_counts' not in group:
                raise DownstreamQCError('Finalized Droplet Cell_Order identity lacks native count totals')
            totals = pd.to_numeric(group['total_expected_counts'],errors='coerce')
            if totals.dtype.kind not in 'iuf' or not (totals.ge(0) & totals.mod(1).eq(0)).all():
                raise DownstreamQCError('Finalized Droplet Cell_Order identity requires finite integer count totals')
            ranked = sorted(zip((-int(value) for value in totals),group['RNA_Barcode'].astype(str),strict=True))
            ranks = {value:str(index+1) for index,(_,value) in enumerate(ranked)}
            if orders.tolist() != group['RNA_Barcode'].astype(str).map(ranks).tolist():
                raise DownstreamQCError('Finalized Droplet Cell_Order identity disagrees with native count/barcode ranking')
        expected_ids = sample + '_' + barcode
    else:
        if 'protocol' in obs and not obs['protocol'].eq('cabernet').all():
            raise DownstreamQCError('Finalized RNA protocol identity disagrees with Cabernet config')
        barcode_map = Path(config.get('barcode_map','00_config/Barcode_Map.csv')).expanduser()
        if not barcode_map.is_absolute():
            barcode_map = Path(project_dir)/barcode_map
        cells = read_delimited_frame(barcode_map)
        columns = ['RNA_Barcode','PlateID','Cell_Order']
        if not set(columns).issubset(cells):
            raise DownstreamQCError('Barcode Map lacks required cell identity columns')
        cells = cells[columns].apply(lambda values: values.str.strip())
        if any(cells[column].eq('').any() or cells[column].duplicated().any() for column in columns):
            raise DownstreamQCError('Barcode Map cell identities must be nonempty and one-to-one')
        cells = cells.set_index('RNA_Barcode')
        for column in ('PlateID','Cell_Order'):
            expected = barcode.map(cells[column])
            if expected.isna().any() or not obs[column].astype(str).equals(expected):
                raise DownstreamQCError(f'Finalized RNA {column} identity disagrees with current Barcode Map')
        expected_ids = sample + '_' + obs['PlateID'].astype(str)
    if list(actual_ids) != expected_ids.tolist():
        raise DownstreamQCError('Finalized RNA cell identities disagree with explicit sample/barcode/plate metadata')
    return manifest_path


def _validate_rna_reference(adata, reference: str) -> None:
    if adata.uns['rna_quantification'].get('species') != reference:
        raise ValueError(f'Finalized RNA reference disagrees with project reference {reference!r}')


def _load_rna_h5ad_obs(h5ad_path: Path, *, reference: str | None = None) -> pd.DataFrame:
    import anndata as ad

    adata = ad.read_h5ad(h5ad_path, backed="r")
    try:
        validate_prepared_rna(adata)
        if reference is not None:
            _validate_rna_reference(adata, reference)
        obs = adata.obs.copy()
        obs.insert(0, "Sample_ID", obs.index.astype(str))
    except ValueError as exc:
        raise DownstreamQCError(f'{h5ad_path}: {exc}') from exc
    finally:
        adata.file.close()
    return obs.reset_index(drop=True)


def _load_clone_table(clone_csv: Path, barcode_map: Path) -> pd.DataFrame:
    clones = read_delimited_frame(clone_csv)
    if 'DARLIN_Clone_ID' not in clones:
        raise DownstreamQCError('DARLIN clone table lacks DARLIN_Clone_ID')
    try:
        validate_barcode_map_provenance(clones,barcode_map)
    except ValueError as exc:
        raise DownstreamQCError(str(exc)) from exc
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
    config = load_config(project_dir/'00_config/config.yaml')
    reference = config['RNA']['reference']
    results_dir = project_results_root(project_dir,config)
    count_dir = results_dir / "count"
    if not count_dir.is_dir():
        raise DownstreamQCError(f"Project has no selected RNA count directory: {count_dir}")

    h5ad_paths = _rna_h5ad_paths(project_dir,config)
    if not h5ad_paths:
        raise DownstreamQCError(f"No finalized <sample_id>/RNA.h5ad found in {count_dir}")

    frames = []
    for h5ad_path in h5ad_paths:
        frame = _load_rna_h5ad_obs(h5ad_path, reference=reference)
        frame["Sample_ID"] = frame["Sample_ID"].astype(str)
        duplicates = frame.loc[frame["Sample_ID"].duplicated(), "Sample_ID"].unique().tolist()
        if duplicates:
            raise DownstreamQCError(
                f"RNA.h5ad obs_names are not unique in {h5ad_path}: "
                + ", ".join(duplicates[:10])
            )
        frames.append(frame)

    rna = pd.concat(frames, ignore_index=True) if len(frames) > 1 else frames[0]
    manifest_path = _validate_rna_libraries(project_dir, rna)

    ordered_identity = [c for c in IDENTITY_COLUMNS if c in rna.columns]
    technical = [c for c in rna.columns if c not in RNA_OBS_EXCLUDE]
    rna = rna[ordered_identity + technical]

    clone_csv = _darlin_clone_path(project_dir,config)
    if clone_csv is not None:
        barcode_map = Path(config.get('barcode_map','00_config/Barcode_Map.csv')).expanduser()
        if not barcode_map.is_absolute():
            barcode_map = project_dir/barcode_map
        clones = _load_clone_table(clone_csv,barcode_map)
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
            f"INFO: DARLIN disabled or clone table absent under {results_dir}; "
            "publishing RNA-only QC table."
        )
        darlin_inputs = {}

    inputs = {"rna_h5ad": [str(p) for p in h5ad_paths], "rna_reference": reference,
              "sample_manifest": str(manifest_path), **darlin_inputs}
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
    """Bind QC to its inputs and executed scientific logic before reusing parameters."""
    project_dir = Path(project_dir).expanduser().resolve()
    digest = hashlib.sha256()
    digest.update(downstream_source_sha256().encode("ascii"))
    digest.update(b"\0")
    config_rel = '00_config/config.yaml'
    config_path = project_dir / config_rel
    config_digest, config_state = stable_file_identity(config_path)
    cfg = load_config(config_path)
    files = [config_rel, str(cfg.get('sample_manifest','00_config/sample_manifest.tsv'))]
    rna_paths = _rna_h5ad_paths(project_dir, cfg)
    files += [os.path.relpath(p,project_dir) for p in rna_paths]
    if cfg.get('pipeline_mode') != 'droplet':
        files += [str(cfg.get('barcode_map','00_config/Barcode_Map.csv'))]
    clone_csv = _darlin_clone_path(project_dir,cfg)
    if clone_csv is not None:
        files.append(os.path.relpath(clone_csv,project_dir))
    states = {}
    for rel in files:
        path = project_dir / Path(rel).expanduser()
        if not path.is_file():
            raise DownstreamQCError(f"Missing project config file: {path}")
        digest.update(rel.encode("utf-8"))
        digest.update(b"\0")
        file_digest, state = (config_digest, config_state) if rel == config_rel else stable_file_identity(path)
        states[path] = state
        digest.update(file_digest.encode("ascii"))
        digest.update(b"\0")
    if (_rna_h5ad_paths(project_dir, cfg) != rna_paths
            or _darlin_clone_path(project_dir, cfg) != clone_csv):
        raise DownstreamQCError("RNA QC input files changed while recording snapshot")
    for path, state in states.items():
        if file_stat_identity(path) != state:
            raise DownstreamQCError(f"RNA QC input changed while recording snapshot: {path}")
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


def run_processor(project_dir: Path, *,
                  dna_project_root: Path | None = None,
                  params: dict | None = None) -> Path:
    started = perf_counter()
    project_dir = project_dir.expanduser().resolve()
    print(f"INFO: RNA QC 开始；项目：{project_dir}", flush=True)
    qc_results = project_dir / "06_downstream"
    if qc_results.is_symlink():
        raise DownstreamQCError(f"QC output directory must not be a symlink: {qc_results}")
    qc_results = qc_results.resolve()

    snapshot_id = _qc_snapshot_id(project_dir)
    rna, inputs = build_rna_qc_table(project_dir)
    if dna_project_root is not None:
        dna, manifest = load_dna_qc(dna_project_root, expected_reference=inputs['rna_reference'])
    out_dir = qc_results
    out_dir.mkdir(parents=True, exist_ok=True)

    with write_lock(out_dir):
        qc_path = out_dir / 'RNA_QC_Information.csv'
        is_droplet = 'protocol' in rna and rna['protocol'].eq('droplet').all()
        qc_keys = {'min_genes_per_cell', 'min_counts_per_cell'}
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
        same_input = isinstance(prior, dict) and prior.get('qc_snapshot_id') == snapshot_id
        completed_qc = same_input
        reused_artifacts = {}
        reference_states = {}
        reused_annotation = None
        annotation_stale = False
        if same_input:
            required = [qc_path.name]
            if is_droplet or analysis_path.is_file():
                required.append(metrics_path.name)
            if analysis_path.is_file():
                annotated = prior.get('annotated_h5ad')
                if annotated not in ('RNA_mouse_marker_annotated.h5ad', 'RNA_Census_consensus_annotated.h5ad'):
                    completed_qc = False
                else:
                    required.append(annotated)
            recorded = prior.get('qc_artifacts', {})
            if not isinstance(recorded, dict):
                completed_qc = False
            if completed_qc:
                for name in required:
                    if name not in recorded:
                        completed_qc = False
                        break
                    try:
                        identity = stable_file_identity(out_dir/name)
                    except FileNotFoundError:
                        completed_qc = False
                        break
                    if identity[0] != recorded[name]:
                        completed_qc = False
                        break
                    reused_artifacts[name] = identity[1]
            if not completed_qc:
                reused_artifacts.clear()
        if completed_qc:
            annotation = prior.get('annotation') or {}
            if (annotation.get('method') == 'census_scvi_consensus'
                    or prior.get('annotated_h5ad') == 'RNA_Census_consensus_annotated.h5ad'):
                recorded_reference = annotation.get('reference')
                current_reference = None
                if recorded_reference:
                    try:
                        current_reference, reference_states = human_annotation_reference_identity(
                            PIPELINE_ROOT, DEFAULT_CENSUS_PARAMS['census_version'])
                    except FileNotFoundError:
                        pass
                annotation_stale = not recorded_reference or current_reference != recorded_reference
                if annotation_stale:
                    reference_states.clear()
                else:
                    reused_annotation = annotation
        if is_droplet and same_input:
            if analysis_path.is_file():
                prior_params = prior['params']
            elif summary_path.is_file():
                prior_params = prior.get('inputs', {}).get('rna_quality_filter')
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
                if annotation_stale or analysis_changed or not metrics['qc_pass'].equals(mask):
                    clear_analysis_outputs(out_dir, input_id=snapshot_id)
                    reused_artifacts.clear()
                    reference_states.clear()
                    reused_annotation = None
                    for column in RNA_QC_ANNOTATION_COLUMNS:
                        published[column] = pd.NA
                    metrics['qc_pass'] = mask
                    atomic_write_csv(metrics, metrics_path)
                    atomic_write_csv(published, qc_path)
                rna['qc_pass'] = rna['Sample_ID'].map(metrics.set_index('Sample_ID')['qc_pass'])
            elif annotation_stale:
                clear_analysis_outputs(out_dir, input_id=snapshot_id)
                reused_artifacts.clear()
                for column in RNA_QC_ANNOTATION_COLUMNS:
                    published[column] = pd.NA
                atomic_write_csv(published, qc_path)
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
            integrated = integrate_dna_qc(published, dna, manifest, rna_identity=rna)
            invalidate_qc_completion(out_dir)
            ledger_path = out_dir / 'DNA_RNA_cell_links.csv'
            if 'cell_links' in integrated.attrs:
                atomic_write_csv(integrated.attrs['cell_links'], ledger_path)
            else:
                ledger_path.unlink(missing_ok=True)
            atomic_write_csv(integrated, out_dir/'Integrated_QC_Information.csv')
            provenance = {**dna.attrs['dna_provenance'],
                          'library_links':integrated.attrs['library_links'],
                          'library_map':integrated.attrs.get('library_map',{}),
                          'rna_qc_snapshot_id':snapshot_id, 'integrated_cells':len(integrated)}
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

        manifest_src = Path(inputs['sample_manifest'])
        manifest_copy = out_dir / "sample_manifest.tsv"
        atomic_write_bytes(manifest_src.read_bytes(), manifest_copy)

        try:
            qc_summary = {
                "schema_version": SCHEMA_VERSION,
                "generated_at_utc": datetime.now(timezone.utc).isoformat(),
                "qc_snapshot_id": snapshot_id,
                "project_dir": str(project_dir),
                "n_cells": int(len(rna)),
                "n_samples": int(rna["sample_id"].nunique()) if "sample_id" in rna.columns else None,
                "columns": list(published.columns),
                "qc_artifacts": {path.name: sha256_file(path) for path in (
                    (qc_path, metrics_path) if is_droplet else (qc_path,))},
                "inputs": inputs,
                **({'annotation': reused_annotation} if reused_annotation is not None else {}),
            }
            for name, state in reused_artifacts.items():
                info = (out_dir/name).stat()
                if tuple(getattr(info, key) for key in (
                        'st_dev', 'st_ino', 'st_size', 'st_mtime_ns', 'st_ctime_ns')) != state:
                    raise DownstreamQCError('Cached RNA QC artifacts changed during QC; retry from current inputs')
            if _qc_snapshot_id(project_dir) != snapshot_id:
                raise DownstreamQCError('RNA输入在QC期间改变；本次输出未完成，请重试')
            if any(file_stat_identity(path) != state for path, state in reference_states.items()):
                raise DownstreamQCError('Human annotation reference changed during cached QC; retry from current inputs')
            if dna_project_root is not None:
                current, _ = load_dna_qc(dna_project_root, expected_reference=inputs['rna_reference'])
                if current.attrs['dna_provenance'] != dna.attrs['dna_provenance']:
                    raise DownstreamQCError('DNA输入在QC期间改变；本次整合未完成，请重试')
                atomic_write_json(provenance, out_dir/'integration_summary.json')
            atomic_write_json(qc_summary, summary_path)
        except BaseException:
            invalidate_qc_completion(out_dir)
            raise

        print(f"INFO: RNA QC 完成；细胞：{len(rna)}；耗时：{perf_counter() - started:.2f}s；输出：{qc_path}", flush=True)
        return out_dir


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m lemmus.downstream.qc", description=__doc__)
    parser.add_argument("--project-dir", required=True, type=Path, help="RNA_DARLIN project root")
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
    parser.add_argument("--min-genes-per-cell", type=int, default=None)
    parser.add_argument("--min-counts-per-cell", type=int, default=None)
    parser.add_argument("--min-cells-per-gene", type=int, default=None)
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
    if sys.platform == 'darwin' and os.uname().machine == 'arm64':
        # This single-operation CLI avoids the Torch/Numba duplicate-libomp crash.
        os.environ['NUMBA_THREADING_LAYER'] = 'workqueue'
        if 'numba' in sys.modules:
            sys.modules['numba'].config.THREADING_LAYER = 'workqueue'
    try:
        if args.analyze:
            from .analysis import run_analysis

            overrides = {
                key: getattr(args, key)
                for key in (
                    "min_genes_per_cell", "min_counts_per_cell", "min_cells_per_gene",
                    "n_top_hvg", "n_pcs", "n_neighbors",
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
            )
        else:
            analysis_only = ('min_cells_per_gene', 'n_top_hvg', 'n_pcs', 'n_neighbors',
                             'leiden_resolution', 'umap_min_dist', 'random_state')
            if any(getattr(args, key) is not None for key in analysis_only) or args.keep_ig_in_pca or args.keep_mt_in_pca:
                raise DownstreamQCError('Analysis parameters require --analyze')
            qc_params = {key: getattr(args, key) for key in
                         ('min_genes_per_cell', 'min_counts_per_cell')
                         if getattr(args, key) is not None}
            run_processor(args.project_dir, dna_project_root=args.dna_project_root,
                          params=qc_params or None)
    except (ConfigError, DownstreamQCError, AnalysisError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    return 0


def detect_species(project_dir: Path) -> str:
    """Species branch from the project config `RNA.reference`."""
    config_path = Path(project_dir) / "00_config" / "config.yaml"
    if not config_path.is_file():
        raise AnalysisError(f"缺少项目配置: {config_path}")
    config = load_config(config_path)
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
    import numpy as np
    import scipy.sparse as sp

    config = load_config(Path(project_dir)/'00_config/config.yaml')
    count_dir = project_results_root(project_dir,config)/'count'
    paths = _rna_h5ad_paths(project_dir,config)
    if not paths:
        raise AnalysisError(f"在 {count_dir} 中未找到 */RNA.h5ad")

    adatas = []
    seen: set[str] = set()
    required_obs = {"RNA_Barcode", "PlateID", "Cell_Order", "sample_id", "rna_sample"}
    for path in paths:
        sample_id = path.parent.name
        x = anndata.read_h5ad(path)
        try:
            validate_prepared_rna(x)
            _validate_rna_reference(x, config['RNA']['reference'])
            validate_rna_count_matrix(x.X)
            x.X = sp.csr_matrix(x.X)
        except ValueError as exc:
            raise AnalysisError(f'{path}: {exc}') from exc
        if (not x.var_names.is_unique or 'gene_id' not in x.var
                or x.var['gene_id'].astype(str).tolist() != x.var_names.tolist()):
            raise AnalysisError(f'{path}: RNA gene identities must be unique and match gene_id')
        if adatas and set(x.var_names) != set(adatas[0].var_names):
            raise AnalysisError(f'{path}: RNA libraries must contain the same full gene set; never fill missing features')
        if 'gene_symbol' not in x.var:
            raise AnalysisError(f'{path}: RNA matrix lacks gene_symbol metadata')
        symbols = x.var['gene_symbol'].astype('string')
        invalid_symbols = symbols.isna() | symbols.str.strip().eq('')
        if invalid_symbols.any():
            bad = x.var_names[invalid_symbols][:5].tolist()
            raise AnalysisError(f'{path}: empty gene_symbol for {bad}')
        if adatas:
            reference_symbols = adatas[0].var['gene_symbol'].reindex(x.var_names)
            conflicts = symbols.ne(reference_symbols)
            if conflicts.any():
                bad = x.var_names[conflicts][:5].tolist()
                raise AnalysisError(f'{path}: gene_symbol disagrees across RNA libraries for {bad}')
        # Compare text values independently of H5AD string or categorical storage.
        x.var['gene_id'] = x.var_names.array
        x.var['gene_symbol'] = symbols.astype(object)
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
        adatas.append(x)

    dtypes = {x.X.dtype for x in adatas}
    if np.result_type(*dtypes).kind == 'f' and any(dtype.kind in 'iu' for dtype in dtypes):
        # Mixed storage can promote exact integer counts to a lossy floating dtype.
        for x in adatas:
            values = x.X.data
            if values.dtype.kind in 'uf' and (values >= 2**63).any():
                raise AnalysisError('RNA counts exceed the signed 64-bit range for combining mixed storage types')
            x.X = x.X.astype(np.int64, copy=False)

    # Preserve gene_id/gene_symbol across samples; anndata.concat defaults
    # to merge=None and would silently discard all var metadata here.
    adata = anndata.concat(adatas, join="inner", merge="same")
    _validate_rna_libraries(project_dir, adata.obs.assign(Sample_ID=adata.obs_names.astype(str)))
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
    counts = sp.csr_matrix(adata.layers["counts"])
    low_precision = counts.dtype.kind == 'f' and counts.dtype.itemsize < 8
    if (counts.data == 0).any() or low_precision:
        from anndata import AnnData

        # Protect raw buffers and accumulate lower precision storage in float64.
        counts = counts.astype(np.float64) if low_precision else counts.copy()
        counts.eliminate_zeros()
        qc_input = AnnData(counts, obs=pd.DataFrame(index=adata.obs_names), var=adata.var[["mt"]].copy())
        obs_metrics, var_metrics = sc.pp.calculate_qc_metrics(qc_input, qc_vars=["mt"], percent_top=[50])
        adata.obs[obs_metrics.columns] = obs_metrics
        adata.var[var_metrics.columns] = var_metrics
    else:
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=[50], layer="counts", inplace=True)
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
    # Weighted bincount forces integer molecules through floating-point sums.
    subunit_counts = sp.csr_array(
        (np.where(subunit, counts.data, 0), counts.indices, counts.indptr),
        shape=counts.shape, copy=False,
    )
    subunit_mass = np.asarray(subunit_counts.sum(
        axis=1, dtype=np.float64 if counts.dtype.kind == 'f' else None,
    )).ravel()
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
    """Filter by supported genes and UMI counts; MT% remains diagnostic only."""
    import pandas as pd

    return (
        (obs["n_genes_supported"] >= params["min_genes_per_cell"])
        & (obs["total_counts"] >= params["min_counts_per_cell"])
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
