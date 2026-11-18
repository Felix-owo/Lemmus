"""独立RNA/DNA doublet判定及同代SingleCpG导出。"""
from __future__ import annotations

import hashlib
import json
import math
import os
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from zipfile import BadZipFile
import pandas as pd
from .config import (
    DownstreamQCError,
    RNA_QC_ANNOTATION_COLUMNS,
    RNA_QC_COLUMNS,
    validate_analysis_params,
)
from ..io import (
    sha256_file,
)
from .integration import (
    integrate_droplet_qc,
    load_dna_qc,
    validate_droplet_integration_links,
)
from ..io import project_temporary_parent
from ..inputs import read_manifest
from .io import (
    write_lock,
    read_verified_qc_csv,
    atomic_write_csv,
    atomic_write_h5ad,
    atomic_write_json,
    downstream_source_sha256,
    sanitize_for_h5ad,
)
from .qc import (
    _qc_snapshot_id,
    build_rna_qc_table,
    compute_rna_qc,
    load_rna_matrix,
    publish_rna_qc_table,
)


def score_droplet_doublets(adata, *, expected_doublet_rate=0.05, random_state=0, doublet_threshold_scale=1.0):
    """逐文库运行 Scrublet，以自动阈值乘倍率判定；降低倍率会增加剔除，失败不得伪装成 singlet。"""
    import numpy as np
    import scanpy as sc
    from ..rna import validate_rna_count_matrix
    if not math.isfinite(expected_doublet_rate) or not 0 < expected_doublet_rate < 1:
        raise DownstreamQCError('expected_doublet_rate must be finite and between 0 and 1')
    if not math.isfinite(doublet_threshold_scale) or not 0 < doublet_threshold_scale <= 1:
        raise DownstreamQCError('doublet_threshold_scale must be finite and in (0, 1]')
    if not adata.obs_names.is_unique:
        raise DownstreamQCError('Duplicate RNA cell IDs')
    try:
        validate_rna_count_matrix(adata.X)
    except ValueError as error:
        raise DownstreamQCError(f'Scrublet requires finite nonnegative integer raw counts: {error}') from error
    records, libraries = [], {}
    for library in sorted(adata.obs['sample_id'].unique()):
        subset = adata[adata.obs['sample_id'].eq(library)].copy()
        if subset.n_obs <= 30 or subset.n_vars <= 30:
            raise DownstreamQCError(f'{library}: 数据不足以运行 30-PC Scrublet；不自动判为 singlet')
        subset.obs.drop(columns=['doublet_score', 'predicted_doublet'], errors='ignore', inplace=True)
        subset.uns.pop('scrublet', None)
        try:
            sc.pp.scrublet(subset, expected_doublet_rate=expected_doublet_rate,
                           random_state=random_state, n_prin_comps=30,
                           use_approx_neighbors=False, verbose=False)
        except (ValueError, RuntimeError, ImportError, KeyError) as error:
            raise DownstreamQCError(f'{library}: Scrublet 失败，未发布 singlet 结果：{error}') from error
        if not {'doublet_score', 'predicted_doublet'}.issubset(subset.obs):
            raise DownstreamQCError(f'{library}: Scrublet 未返回完整判定；自动阈值可能失败')
        scores = subset.obs['doublet_score']
        labels = subset.obs['predicted_doublet']
        threshold = subset.uns.get('scrublet', {}).get('threshold')
        if (threshold is None or not np.isfinite(threshold) or scores.dtype.kind not in 'iuf'
                or not scores.between(0, 1).fillna(False).all() or labels.isna().any()):
            raise DownstreamQCError(f'{library}: Scrublet 分数或自动阈值无效；不把未判定细胞当 singlet')
        if not pd.api.types.is_bool_dtype(labels.dtype):
            raise DownstreamQCError(f'{library}: Scrublet 返回非布尔判定')
        auto_threshold = float(threshold)
        if auto_threshold <= 0:
            raise DownstreamQCError(f'{library}: Scrublet 自动阈值无效')
        threshold = auto_threshold * doublet_threshold_scale
        n_doublets_auto = int(scores.gt(auto_threshold).sum())
        labels = scores.gt(threshold)
        frame = pd.DataFrame({'Sample_ID': subset.obs_names, 'CloneID': library,
                              'doublet_score': scores.to_numpy(), 'predicted_doublet': labels.to_numpy(dtype=bool)})
        records.append(frame)
        libraries[str(library)] = {'n_cells': subset.n_obs, 'auto_threshold': auto_threshold,
                                  'threshold': float(threshold), 'threshold_scale': doublet_threshold_scale,
                                  'n_doublets_auto': n_doublets_auto, 'n_doublets': int(labels.sum())}
    result = pd.concat(records, ignore_index=True).set_index('Sample_ID').loc[adata.obs_names].reset_index()
    result = result.rename(columns={result.columns[0]: 'Sample_ID'})
    return result, libraries


def _droplet_rna_libraries(rna_root: Path) -> pd.DataFrame:
    from lemmus.config import load_config

    config = load_config(rna_root / '00_config/config.yaml')
    if config.get('pipeline_mode') != 'droplet':
        raise DownstreamQCError('此入口只接受 Droplet RNA 项目')
    manifest_path = Path(config.get('sample_manifest', '00_config/sample_manifest.tsv')).expanduser()
    fields, rows = read_manifest(rna_root / manifest_path)
    manifest = pd.DataFrame(rows, columns=fields)
    return manifest.loc[manifest.rna_sample.ne(''), ['sample_id', 'rna_sample']]


def run_droplet_qc_integration(rna_project_root, dna_project_root, *, expected_doublet_rate=0.05, random_state=0,
                              doublet_threshold_scale=1.0):
    """仅整合 Droplet QC、预测并去除 RNA doublet，发布 raw-count singlet AnnData。"""
    import importlib.metadata
    import numpy as np

    roots = []
    for value in (rna_project_root, dna_project_root):
        if value is None or not str(value).strip():
            raise DownstreamQCError('请填写 RNA_PROJECT_ROOT 和 DNA_PROJECT_ROOT 两个项目目录')
        roots.append(Path(value).expanduser().resolve())
    rna_root, dna_root = roots
    snapshot = _qc_snapshot_id(rna_root)
    selected = _droplet_rna_libraries(rna_root)
    full, rna_inputs = build_rna_qc_table(rna_root)
    if 'protocol' not in full or not full.protocol.eq('droplet').all():
        raise DownstreamQCError('RNA.h5ad 必须包含 finalized Droplet 身份')
    observed = full[['sample_id', 'rna_sample']].drop_duplicates()
    if (any(selected[c].duplicated().any() for c in selected)
            or set(map(tuple, selected.to_numpy())) != set(map(tuple, observed.to_numpy()))):
        raise DownstreamQCError(
            f'RNA finalized 文库与项目 manifest 不一致：manifest={selected.to_dict("records")}；'
            f'RNA.h5ad={observed.to_dict("records")}。'
            'sample_id 是项目样本名，rna_sample 是原始 FASTQ 文库名；须与已完成结果保持一致。')
    if not full.Sample_ID.eq(full.sample_id.astype(str) + '_' + full.RNA_Barcode.astype(str)).all():
        raise DownstreamQCError('RNA cell ID 与 finalized sample_id/barcode 不一致')
    dna, dm = load_dna_qc(dna_root, expected_reference=rna_inputs['rna_reference'])
    adata = load_rna_matrix(rna_root)
    compute_rna_qc(adata, validate_analysis_params(None))
    full['pct_counts_mt'] = full.Sample_ID.map(adata.obs['pct_counts_mt'])
    rna = publish_rna_qc_table(full)
    integrated = integrate_droplet_qc(rna, dna, dm, full)
    doublets, libraries = score_droplet_doublets(
        adata, expected_doublet_rate=expected_doublet_rate, random_state=random_state,
        doublet_threshold_scale=doublet_threshold_scale)
    by_cell = doublets.set_index('Sample_ID')
    for column in ('doublet_score', 'predicted_doublet'):
        integrated[column] = integrated.Sample_ID.map(by_cell[column])
    ledger = integrated.attrs['cell_links'].copy()
    ledger['doublet_score'] = ledger.rna_cell_id.map(by_cell.doublet_score)
    ledger['predicted_doublet'] = ledger.rna_cell_id.map(by_cell.predicted_doublet).astype('boolean')
    ledger['doublet_status'] = np.where(ledger.rna_cell_id.eq(''), 'not_assessed_no_rna',
                                      np.where(ledger.predicted_doublet.fillna(False), 'doublet', 'singlet'))
    keep = integrated.loc[~integrated.predicted_doublet.astype(bool)].copy()
    singlets = adata[keep.Sample_ID.tolist()].copy()
    attached = keep.set_index('Sample_ID').reindex(singlets.obs_names)
    for column in attached:
        if column in singlets.obs and column in ('CloneID', 'PlateID'):
            if not singlets.obs[column].astype(str).equals(attached[column].astype(str)):
                raise DownstreamQCError(f'AnnData/QC identity conflict: {column}')
        singlets.obs[column] = attached[column].to_numpy()
    singlets.obs['Sample_ID'] = singlets.obs_names.astype(str)
    for column in (*RNA_QC_ANNOTATION_COLUMNS, 'DARLIN_Clone_ID'):
        singlets.obs[column] = singlets.obs[column].fillna('').astype(str)
    versions = {name: importlib.metadata.version(name) for name in ('scanpy', 'anndata', 'numpy', 'pandas', 'scipy', 'scikit-image')}
    provenance = {'rna_project': str(rna_root), 'rna_qc_snapshot_id': snapshot,
                  'dna': dna.attrs['dna_provenance'], 'library_map': integrated.attrs['library_map'],
                  'library_links': integrated.attrs['library_links'],
                  'doublet_method': 'scanpy.pp.scrublet', 'expected_doublet_rate': expected_doublet_rate,
                  'doublet_threshold_scale': doublet_threshold_scale,
                  'doublet_decision': 'doublet_score > auto_threshold * doublet_threshold_scale',
                  'random_state': random_state, 'n_prin_comps': 30, 'use_approx_neighbors': False,
                  'input_counts': 'unrounded native molecule counts',
                  'selection': 'both modalities QC available and RNA predicted_doublet == False; no additional cell QC filter',
                  'versions': versions, 'processor_sha256': downstream_source_sha256()}
    out = rna_root / '06_downstream' / 'droplet_integration'
    singlets.uns['droplet_qc_integration'] = json.dumps(provenance, sort_keys=True)
    def check_inputs(message):
        if _qc_snapshot_id(rna_root) != snapshot:
            raise DownstreamQCError(message)
        current, _ = load_dna_qc(dna_root, expected_reference=rna_inputs['rna_reference'])
        if current.attrs['dna_provenance'] != provenance['dna']:
            raise DownstreamQCError(message)

    out.mkdir(parents=True, exist_ok=True)
    with write_lock(out):
        # 输出锁不保护DNA项目；覆盖旧产物前与发布摘要前均复核完整身份。
        check_inputs('DNA/RNA 输入在分析期间改变，请重试')
        child = out / 'dna_doublet'
        plots = out / 'plots'
        for directory, label in ((child, 'DNA doublet'), (plots, 'QC plots')):
            if directory.is_symlink():
                raise DownstreamQCError(f'{label} directory must not be a symlink: {directory}')
            if directory.exists() and not directory.is_dir():
                raise DownstreamQCError(f'{label} path is not a directory: {directory}')
        marker = out / 'integration_summary.json'
        marker.unlink(missing_ok=True)
        for directory, names in (
            (out, ('SingleCpG_Adata.h5ad',)),
            (child, ('DNA_doublet_summary.json', 'DNA_doublet_scores.csv', 'SingleCpG_Adata.h5ad')),
        ):
            for name in names:
                (directory / name).unlink(missing_ok=True)
        for stem in ('01_membership_doublets', '02_joint_qc_doublets', '03_singlet_qc'):
            for suffix in ('.pdf', '.png'):
                (plots/f'{stem}{suffix}').unlink(missing_ok=True)
        atomic_write_csv(doublets, out / 'RNA_doublet_scores.csv')
        atomic_write_csv(ledger, out / 'DNA_RNA_cell_links.csv')
        atomic_write_csv(integrated, out / 'Integrated_QC_Information.csv')
        atomic_write_csv(keep, out / 'Integrated_QC_singlets.csv')
        atomic_write_h5ad(sanitize_for_h5ad(singlets), out / 'RNA_DNA_singlets.h5ad')
        summary = {**provenance, 'libraries': libraries, 'rna_cells': adata.n_obs,
                   'integrated_cells': len(integrated), 'integrated_doublets': int(integrated.predicted_doublet.sum()),
                   'singlet_cells': singlets.n_obs, 'generated_at_utc': datetime.now(timezone.utc).isoformat(),
                   'qc_artifacts': {name: sha256_file(out / name) for name in (
                       'Integrated_QC_Information.csv', 'DNA_RNA_cell_links.csv', 'RNA_doublet_scores.csv',
                       'Integrated_QC_singlets.csv', 'RNA_DNA_singlets.h5ad')}}
        check_inputs('DNA/RNA 输入在写入期间改变；本次整合未完成，请重试')
        atomic_write_json(summary, marker)
        print(f'QC 整合 {len(integrated)} cells → {singlets.n_obs} singlets；输出：{out}')
        return out


def _dna_doublet_input(dna):
    """绑定小型整数计数产物。状态复验无需重复解码矩阵。"""
    provenance = dna.attrs["dna_provenance"]
    source = Path(provenance["qc_path"]).parent / "DNA_Doublet_1Mb_Counts.npz"
    try:
        state = source.lstat()
        signature = dict(name=source.name, device=state.st_dev, inode=state.st_ino,
                         mtime_ns=state.st_mtime_ns, size=state.st_size)
        if (source.is_symlink() or not source.is_file() or not state.st_size
                or provenance["doublet_counts"] != signature):
            raise ValueError("stale integer-count artifact")
    except (KeyError, OSError, ValueError) as error:
        raise DownstreamQCError("同代 1 Mb DNA 整数计数缺失或已改变，请重跑 Alopex QC") from error
    return dict(source=str(source), source_identity=signature, ctime_ns=state.st_ctime_ns,
                cpg_generation=provenance["cpg_generation"], delivery_id=provenance["delivery_id"],
                cpg_coordinate_base=provenance["cpg_coordinate_base"])


def _load_dna_doublet_counts(dna, identity):
    import numpy as np

    try:
        with np.load(identity["source"], allow_pickle=False) as archive:
            if set(archive.files) != {"ids", "regions", "mc", "cov"}:
                raise ValueError("invalid integer-count fields")
            ids = archive["ids"]
            regions = archive["regions"]
            mc, cov = archive["mc"], archive["cov"]
        if (ids.ndim != 1 or ids.dtype.kind != "U" or len(ids) != len(dna)
                or set(ids) != set(dna.Sample_ID) or regions.ndim != 1 or regions.dtype.kind != "U"
                or len(regions) < 2 or mc.shape != (len(ids), len(regions)) or mc.shape != cov.shape
                or mc.dtype != np.int64 or cov.dtype != np.int64 or (mc < 0).any() or (mc > cov).any()):
            raise ValueError("invalid integer-count shape, type, cell IDs or counts")
        base = identity["cpg_coordinate_base"]
        last_chrom, last_end, seen = None, None, set()
        for label in regions:
            match = re.fullmatch(r"(chr(?:[1-9]|1[0-9]|2[0-2])):([0-9]+)-([0-9]+)", label)
            if match is None:
                raise ValueError("invalid autosomal 1 Mb window name")
            chrom, start, end = match.groups()
            start, end = int(start), int(end)
            if (not 0 < end-start <= 1_000_000
                    or (chrom == last_chrom and (start != last_end or (start-base)%1_000_000 != 0))
                    or (chrom != last_chrom and (start != base or chrom in seen))):
                raise ValueError("invalid autosomal 1 Mb window coordinates")
            seen.add(chrom)
            last_chrom, last_end = chrom, end
    except (BadZipFile, EOFError, KeyError, ValueError, OSError) as error:
        raise DownstreamQCError("同代 DNA 1 Mb 计数无法读取或内容非法，请重跑 Alopex QC") from error
    if _dna_doublet_input(dna) != identity:
        raise DownstreamQCError("DNA 1 Mb 计数在读取期间改变，请重试")
    return mc, cov, ids


def score_dna_doublets(mc, cov, identities, *, expected_doublet_rate=0.06, random_state=0,
                       doublet_threshold_scale=1.0):
    """逐 DNA 文库以真实区域 mc/cov 调用 MethylScrublet；低信息细胞保持未判定，降低阈值倍率增加标记。"""
    import numpy as np
    import numba
    from threadpoolctl import threadpool_limits
    from ALLCools.clustering.doublets import MethylScrublet
    from ALLCools.mcds.utilities import calculate_posterior_mc_frac

    if (not math.isfinite(expected_doublet_rate) or not 0 < expected_doublet_rate < 1
            or not math.isfinite(doublet_threshold_scale) or not 0 < doublet_threshold_scale <= 1):
        raise DownstreamQCError('DNA expected_doublet_rate 须在 (0,1)，doublet_threshold_scale 须在 (0,1]')
    mc, cov = np.asarray(mc), np.asarray(cov)
    if (mc.ndim != 2 or mc.shape != cov.shape or mc.shape[0] != len(identities)
            or not identities.Sample_ID.is_unique or identities.CloneID.eq('').any()
            or not np.isfinite(mc).all() or not np.isfinite(cov).all()
            or (mc < 0).any() or (cov < mc).any()
            or not np.equal(mc, np.floor(mc)).all() or not np.equal(cov, np.floor(cov)).all()):
        raise DownstreamQCError('DNA doublet 输入必须为身份唯一、有限且满足 0≤mc≤cov 的整数计数')
    records, libraries = [], {}
    for library in sorted(identities.CloneID.unique()):
        rows = np.flatnonzero(identities.CloneID.eq(library))
        m, c = mc[rows].astype(float), cov[rows].astype(float)
        covered = (c > 0).sum(axis=1)
        with np.errstate(divide='ignore', invalid='ignore'):
            posterior = calculate_posterior_mc_frac(m, c)
        assessed = (covered >= 2) & np.isfinite(posterior).all(axis=1)
        informative = c[assessed].sum(axis=0) > 0
        n_cells, n_regions = int(assessed.sum()), int(informative.sum())
        if n_regions < 2 or n_cells < max(12, min(200, n_regions)):
            raise DownstreamQCError(f'{library}: DNA 数据不足以运行 MethylScrublet；不自动判为 singlet')
        model = MethylScrublet(sim_doublet_ratio=2.0, n_neighbors=10,
                               expected_doublet_rate=expected_doublet_rate, stdev_doublet_rate=0.02,
                               metric='euclidean', random_state=random_state, n_jobs=1)
        previous_threads, previous_rng = numba.get_num_threads(), np.random.get_state()
        try:
            numba.set_num_threads(1)
            with threadpool_limits(limits=1):
                model.fit(m[assessed][:, informative], c[assessed][:, informative])
        except (ValueError, RuntimeError, ImportError, TypeError) as error:
            raise DownstreamQCError(f'{library}: MethylScrublet 失败，未发布 DNA 判定：{error}') from error
        finally:
            numba.set_num_threads(previous_threads)
            np.random.set_state(previous_rng)
        scores = np.asarray(model.doublet_scores_obs_, dtype=float)
        auto = model.threshold_
        if (scores.shape != (n_cells,) or not np.isfinite(scores).all() or (scores < 0).any()
                or (scores > 1).any() or auto is None or not np.isfinite(auto) or auto <= 0):
            raise DownstreamQCError(f'{library}: MethylScrublet 分数或自动阈值无效')
        threshold = float(auto) * doublet_threshold_scale
        frame = identities.iloc[rows][['Sample_ID', 'CloneID']].copy()
        frame['DNA_doublet_score'] = np.nan
        frame['DNA_predicted_doublet'] = pd.array([pd.NA] * len(rows), dtype='boolean')
        frame['DNA_doublet_status'] = 'not_assessed_low_information'
        frame.loc[assessed, 'DNA_doublet_score'] = scores
        frame.loc[assessed, 'DNA_predicted_doublet'] = scores > threshold
        frame.loc[assessed, 'DNA_doublet_status'] = np.where(scores > threshold, 'doublet', 'singlet')
        records.append(frame)
        libraries[str(library)] = dict(n_cells=len(rows), n_assessed=n_cells, n_unassessed=int((~assessed).sum()),
                                      n_regions=n_regions, auto_threshold=float(auto), threshold=threshold,
                                      threshold_scale=doublet_threshold_scale,
                                      n_doublets_auto=int((scores > auto).sum()), n_doublets=int((scores > threshold).sum()))
    result = pd.concat(records).set_index('Sample_ID').loc[identities.Sample_ID].reset_index()
    return result, libraries


def run_dna_doublet_analysis(result_dir, dna_project_root, *, expected_doublet_rate=0.06,
                            random_state=0, doublet_threshold_scale=1.0):
    """从同代常染色体 1 Mb mc/cov 分析全部 DNA，独立发布判定，不改变 RNA singlet 筛选。"""
    import importlib.metadata

    parent = Path(result_dir).expanduser().resolve()
    dna_root = Path(dna_project_root).expanduser().resolve()
    marker = parent / 'integration_summary.json'
    marker_bytes = marker.read_bytes()
    integration = json.loads(marker_bytes)
    dna, dna_manifest = load_dna_qc(dna_root)
    rna_root = Path(integration['rna_project'])
    validate_droplet_integration_links(integration, dna_manifest, _droplet_rna_libraries(rna_root))
    if integration['dna'] != dna.attrs['dna_provenance']:
        raise DownstreamQCError('RNA 整合与当前 DNA 交付不一致')
    if _qc_snapshot_id(rna_root) != integration['rna_qc_snapshot_id']:
        raise DownstreamQCError('RNA 结果已改变，请重新运行整合')
    identity = _dna_doublet_input(dna)
    ledger_path = parent / 'DNA_RNA_cell_links.csv'
    if sha256_file(ledger_path) != integration['qc_artifacts']['DNA_RNA_cell_links.csv']:
        raise DownstreamQCError('RNA doublet 对应表已改变')
    versions = {name: importlib.metadata.version(name) for name in
                ('ALLCools', 'numpy', 'pandas', 'scikit-learn', 'pynndescent', 'numba')}
    provenance = dict(dna=integration['dna'], counts=identity,
                      integration_summary_sha256=hashlib.sha256(marker_bytes).hexdigest(),
                      integration_dir=str(parent), method='ALLCools.clustering.doublets.MethylScrublet',
                      bin_size=1000000, chromosomes='chr1-chr22 numeric autosomes present in reference',
                      input_counts='Alopex CpG: summed methylation and coverage integers per 1 Mb bin',
                      expected_doublet_rate=expected_doublet_rate, random_state=random_state,
                      sim_doublet_ratio=2.0, n_neighbors=10, stdev_doublet_rate=0.02, n_jobs=1,
                      doublet_threshold_scale=doublet_threshold_scale,
                      decision='DNA_doublet_score > auto_threshold * doublet_threshold_scale',
                      selection='all DNA cells; marks only; RNA singlet filter unchanged',
                      versions=versions, processor_sha256=downstream_source_sha256())

    def check_inputs():
        if (_dna_doublet_input(dna) != identity
                or marker.read_bytes() != marker_bytes
                or _qc_snapshot_id(rna_root) != integration['rna_qc_snapshot_id']
                or sha256_file(ledger_path) != integration['qc_artifacts']['DNA_RNA_cell_links.csv']):
            raise DownstreamQCError('DNA doublet 输入在分析期间改变，拒绝发布或复用')
        current, _ = load_dna_qc(dna_root)
        if current.attrs['dna_provenance'] != integration['dna']:
            raise DownstreamQCError('DNA doublet 当前交付/QC已改变，拒绝发布或复用')

    out = parent / 'dna_doublet'
    out.mkdir(parents=True, exist_ok=True)
    with write_lock(parent), write_lock(out):
        completed = out / 'DNA_doublet_summary.json'
        if completed.is_file():
            saved = json.loads(completed.read_bytes())
            artifacts = saved.get('artifacts', {})
            if (all(saved.get(key) == value for key, value in provenance.items())
                    and isinstance(artifacts, dict) and set(artifacts) == {'DNA_doublet_scores.csv'}):
                score_path = out/'DNA_doublet_scores.csv'
                try:
                    state = score_path.stat()
                    matches = sha256_file(score_path) == artifacts[score_path.name]
                except FileNotFoundError:
                    matches = False
                if matches:
                    check_inputs()
                    after = score_path.stat()
                    if any(getattr(state, key) != getattr(after, key) for key in (
                            'st_dev', 'st_ino', 'st_size', 'st_mtime_ns', 'st_ctime_ns')):
                        raise DownstreamQCError('Cached DNA doublet scores changed during reuse; retry')
                    return out
        mc, cov, ids = _load_dna_doublet_counts(dna, identity)
        identities = dna.set_index('Sample_ID').loc[ids].rename_axis('Sample_ID').reset_index()
        scores, libraries = score_dna_doublets(mc, cov, identities, expected_doublet_rate=expected_doublet_rate,
                                             random_state=random_state, doublet_threshold_scale=doublet_threshold_scale)
        links = read_verified_qc_csv(
            ledger_path, integration['qc_artifacts'][ledger_path.name], dtype=str,
        ).set_index('dna_cell_id')
        links = links.loc[links.index != '']
        if not links.index.is_unique or not set(ids).issubset(links.index):
            raise DownstreamQCError('DNA/RNA 对应表缺少 DNA 细胞或身份重复')
        rna_labels = links.reindex(scores.Sample_ID).predicted_doublet
        if not rna_labels.isin(['True', 'False', '']).all():
            raise DownstreamQCError('RNA doublet 判定无效')
        scores['RNA_predicted_doublet'] = pd.array(rna_labels.map({'True': True, 'False': False}), dtype='boolean')
        scores['joint_predicted_doublet'] = scores.DNA_predicted_doublet | scores.RNA_predicted_doublet
        check_inputs()
        out.mkdir(parents=True, exist_ok=True)
        completed.unlink(missing_ok=True)
        (out / 'SingleCpG_Adata.h5ad').unlink(missing_ok=True)
        atomic_write_csv(scores, out / 'DNA_doublet_scores.csv')
        summary = dict(provenance, libraries=libraries, dna_cells=len(scores),
                       dna_doublets=int(scores.DNA_predicted_doublet.sum()),
                       dna_unassessed=int(scores.DNA_predicted_doublet.isna().sum()),
                       joint_doublets=int(scores.joint_predicted_doublet.sum()),
                       generated_at_utc=datetime.now(timezone.utc).isoformat(),
                       artifacts={'DNA_doublet_scores.csv': sha256_file(out / 'DNA_doublet_scores.csv')})
        check_inputs()
        atomic_write_json(summary, completed)
        return out


def export_droplet_single_cpg(result_dir, dna_project_root, *, dna_doublet_result_dir=None):
    """保留同代 SingleCpG 全部 DNA 细胞及矩阵，加入 DNA/RNA QC 和可空 doublet 判定后原子导出。"""
    import shutil
    import tempfile
    import h5py
    import hdf5plugin
    import numpy as np
    from anndata.io import read_elem, write_elem

    out = Path(result_dir).expanduser().resolve()
    dna_root = Path(dna_project_root).expanduser().resolve()
    marker = out / 'integration_summary.json'
    marker_bytes = marker.read_bytes()
    summary = json.loads(marker_bytes)
    dna, dna_manifest = load_dna_qc(dna_root)
    rna_root = Path(summary['rna_project'])
    validate_droplet_integration_links(summary, dna_manifest, _droplet_rna_libraries(rna_root))
    provenance = dna.attrs['dna_provenance']
    if summary['dna'] != provenance:
        raise DownstreamQCError('整合结果与当前 DNA 交付/QC 不一致，请重新运行整合')
    dna_doublet = None
    destination = out
    if dna_doublet_result_dir is not None:
        destination = Path(dna_doublet_result_dir).expanduser().resolve()
        dna_marker = destination / 'DNA_doublet_summary.json'
        dna_marker_bytes = dna_marker.read_bytes()
        dna_doublet = json.loads(dna_marker_bytes)
        if (dna_doublet['dna'] != provenance or dna_doublet['integration_dir'] != str(out)
                or dna_doublet['integration_summary_sha256'] != hashlib.sha256(marker_bytes).hexdigest()):
            raise DownstreamQCError('DNA doublet 分析与当前 DNA/RNA 整合不一致')

    def check_inputs():
        if (marker.read_bytes() != marker_bytes
                or _qc_snapshot_id(rna_root) != summary['rna_qc_snapshot_id']):
            raise DownstreamQCError('SingleCpG 导出输入已改变，请重新运行整合')
        for name in ('Integrated_QC_Information.csv', 'DNA_RNA_cell_links.csv'):
            expected = summary.get('qc_artifacts', {}).get(name)
            if not expected or sha256_file(out / name) != expected:
                raise DownstreamQCError('整合 QC 文件缺少有效哈希，请重新运行整合')
        if dna_doublet is not None:
            if (dna_marker.read_bytes() != dna_marker_bytes
                    or _dna_doublet_input(dna) != dna_doublet['counts']
                    or sha256_file(destination / 'DNA_doublet_scores.csv')
                    != dna_doublet['artifacts']['DNA_doublet_scores.csv']):
                raise DownstreamQCError('DNA doublet 输入/判定已改变，请重新运行 DNA 分析')
        current, _ = load_dna_qc(dna_root)
        if current.attrs['dna_provenance'] != provenance:
            raise DownstreamQCError('SingleCpG 当前交付/QC已改变，请重新运行整合')

    check_inputs()
    rna_counts = {'num_genes_expressed', 'num_genes_over_mean'}
    rna_rates = {'mapping_rate', 'dedup_rate', 'pct_counts_mt'}
    rna_dtypes = dict.fromkeys(RNA_QC_COLUMNS, str)
    rna_dtypes.update(dict.fromkeys(rna_counts, 'Int64'))
    rna_dtypes.update(dict.fromkeys(rna_rates, float))
    artifacts = summary['qc_artifacts']
    common = read_verified_qc_csv(
        out / 'Integrated_QC_Information.csv', artifacts['Integrated_QC_Information.csv'],
        dtype=rna_dtypes, na_values=dict.fromkeys(rna_counts | rna_rates, ['']), float_precision='round_trip',
    )
    ledger = read_verified_qc_csv(
        out / 'DNA_RNA_cell_links.csv', artifacts['DNA_RNA_cell_links.csv'],
        dtype=defaultdict(lambda: str, {'doublet_score': float}),
        na_values={'doublet_score': ['']}, float_precision='round_trip',
    )
    if not set(RNA_QC_COLUMNS).issubset(common) or common.Sample_ID.duplicated().any():
        raise DownstreamQCError('整合表缺少完整 RNA QC 或包含重复细胞')
    dna_links = ledger.loc[ledger.dna_cell_id.ne('')].set_index('dna_cell_id')
    if not dna_links.index.is_unique or not set(dna.Sample_ID).issubset(dna_links.index):
        raise DownstreamQCError('DNA/RNA 对应表缺少 DNA 细胞或包含重复身份')
    links = dna_links.reindex(dna.Sample_ID)
    matched = links.rna_cell_id.ne('') & links.rna_qc_available.eq('True')
    if (links.loc[matched, 'rna_cell_id'].duplicated().any()
            or set(links.loc[matched, 'rna_cell_id']) != set(common.Sample_ID)):
        raise DownstreamQCError('RNA QC 与 DNA/RNA 对应表不一致')
    labels = links.predicted_doublet
    if (not labels[matched].isin(['True', 'False']).all()
            or labels[~matched].ne('').any()):
        raise DownstreamQCError('RNA doublet 判定缺失或 DNA-only 被错误判定')

    source = Path(provenance['qc_path']).parent / 'SingleCpG_Adata.h5ad'
    schema_path = source.with_name('SingleCpG_Adata.input_schema.json')
    schema_bytes = schema_path.read_bytes()
    schema = json.loads(schema_bytes)
    before = source.stat()
    signature = dict(device=before.st_dev, inode=before.st_ino, mtime_ns=before.st_mtime_ns,
                     name=source.name, size=before.st_size)
    if (schema.get('input', {}).get('delivery_id') != provenance['delivery_id']
            or schema.get('h5ad') != signature):
        raise DownstreamQCError('SingleCpG 矩阵未绑定当前 DNA 交付或文件身份已改变')
    generation = provenance["cpg_generation"]
    if (schema["input"].get("cpg_generation") != generation
            or schema["input"].get("cpg_coordinate_base") != provenance["cpg_coordinate_base"]):
        raise DownstreamQCError("SingleCpG 与 DNA QC 的 CpG 代际或坐标不一致")
    if dna_doublet is not None and dna_doublet["counts"]["cpg_generation"] != generation:
        raise DownstreamQCError("SingleCpG 与 DNA doublet 的 CpG 代际不一致")
    with h5py.File(source, 'r') as handle:
        obs = read_elem(handle['obs'])
    ids = obs.index.astype(str)
    if not ids.is_unique or set(ids) != set(dna.Sample_ID):
        raise DownstreamQCError('SingleCpG 与当前 DNA QC 的细胞集合不一致')
    aligned_dna = dna.set_index('Sample_ID').reindex(ids)
    for column in ('Sample_ID', 'CloneID', 'PlateID'):
        expected = pd.Series(ids, index=ids) if column == 'Sample_ID' else aligned_dna[column]
        if column in obs and obs[column].astype(object).fillna('').astype(str).tolist() != expected.tolist():
            raise DownstreamQCError(f'SingleCpG 身份与 DNA QC 冲突：{column}')
    obs['Sample_ID'] = ids
    for column in aligned_dna:
        obs[column] = aligned_dna[column].to_numpy()
    links = dna_links.reindex(ids)
    rna = common.set_index('Sample_ID').reindex(links.rna_cell_id)
    for column in RNA_QC_COLUMNS:
        values = pd.Series(links.rna_cell_id.to_numpy()) if column == 'Sample_ID' else rna[column]
        if column in rna_counts:
            projected = values.array
        elif column in rna_rates:
            projected = values.to_numpy()
        else:
            projected = pd.Categorical(values.mask(values.eq('')))
        obs['RNA_' + column] = projected
    for column in links:
        values = links[column].reset_index(drop=True)
        if column == 'doublet_score':
            obs[column] = values.to_numpy()
        elif column in ('predicted_doublet', 'dna_qc_available', 'rna_qc_available', 'rna_qc_pass'):
            if not values.isin(['True', 'False', '']).all():
                raise DownstreamQCError(f'对应表包含非法布尔值：{column}')
            obs[column] = pd.array(values.map({'True': True, 'False': False}), dtype='boolean')
        else:
            obs[column] = pd.Categorical(values.mask(values.eq('')))
    if dna_doublet is not None:
        scores = read_verified_qc_csv(
            destination / 'DNA_doublet_scores.csv', dna_doublet['artifacts']['DNA_doublet_scores.csv'],
            dtype=defaultdict(lambda: str, {'DNA_doublet_score': float}),
            na_values={'DNA_doublet_score': ['']}, float_precision='round_trip',
        )
        if (not scores.Sample_ID.is_unique or set(scores.Sample_ID) != set(ids)
                or not {'DNA_doublet_score', 'DNA_predicted_doublet', 'DNA_doublet_status',
                        'RNA_predicted_doublet', 'joint_predicted_doublet', 'CloneID'}.issubset(scores)):
            raise DownstreamQCError('DNA doublet 表缺少完整细胞集合或判定列')
        scores = scores.set_index('Sample_ID').reindex(ids)
        if scores.CloneID.tolist() != aligned_dna.CloneID.astype(str).tolist():
            raise DownstreamQCError('DNA doublet 文库身份与 DNA QC 冲突')
        for column in ('DNA_predicted_doublet', 'RNA_predicted_doublet', 'joint_predicted_doublet'):
            values = scores[column]
            if not values.isin(['True', 'False', '']).all():
                raise DownstreamQCError(f'DNA doublet 包含非法布尔值：{column}')
            obs[column] = pd.array(values.map({'True': True, 'False': False}), dtype='boolean')
        if (not obs.RNA_predicted_doublet.equals(obs.predicted_doublet)
                or not obs.joint_predicted_doublet.equals(obs.DNA_predicted_doublet | obs.RNA_predicted_doublet)):
            raise DownstreamQCError('DNA/RNA 联合 doublet 与独立判定不一致')
        obs['DNA_doublet_score'] = scores.DNA_doublet_score.to_numpy()
        obs['DNA_doublet_status'] = pd.Categorical(scores.DNA_doublet_status)
    export_provenance = dict(source=str(source), source_identity=signature,
                             dna=provenance, integration_summary_sha256=hashlib.sha256(marker_bytes).hexdigest(),
                             qc_artifacts=summary['qc_artifacts'], selection='all DNA cells including doublets and DNA-only',
                             rna_qc_prefix='RNA_', processor_sha256=downstream_source_sha256())
    if dna_doublet is not None:
        export_provenance['dna_doublet'] = dna_doublet
    target = destination / 'SingleCpG_Adata.h5ad'
    if target == source:
        raise DownstreamQCError('SingleCpG 输出不得覆盖输入矩阵')
    with write_lock(out), write_lock(destination):
        fd, temporary = tempfile.mkstemp(prefix='.SingleCpG-', suffix='.h5ad', dir=project_temporary_parent(target))
        os.close(fd)
        tmp = Path(temporary)
        try:
            shutil.copyfile(source, tmp)
            with h5py.File(tmp, 'r+') as handle:
                del handle['obs']
                write_elem(handle, 'obs', obs)
                write_elem(handle.require_group('uns'), 'droplet_single_cpg_export',
                           json.dumps(export_provenance, sort_keys=True))
            check_inputs()
            after = source.stat()
            if ((before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
                    != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
                    or schema_path.read_bytes() != schema_bytes):
                raise DownstreamQCError('SingleCpG 输入在导出期间改变，拒绝发布')
            os.replace(tmp, target)
        finally:
            tmp.unlink(missing_ok=True)
        print(f'SingleCpG 全部 DNA cells：{len(obs)}；RNA QC 已匹配：{int(matched.sum())}；输出：{target}')
        return target
