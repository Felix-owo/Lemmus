"""Lemmus：RNA calling、原生计数与项目细胞身份。"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import tempfile
import time
from collections.abc import Sequence
from pathlib import Path
from typing import Any
import numpy as np
import pandas as pd
from .paths import PIPELINE_ROOT
from .config import (
    DROPLET_R2_ADAPTERS,
    DROPLET_SOURCE_COMMIT,
    DROPLET_WHITELIST_ROWS,
    DROPLET_WHITELIST_SHA256,
    DROPLET_WHITELIST_URL,
    PROTOCOLS,
    RNA_RESOLUTION,
    RRNA_COMPETITOR_ID,
    RRNA_REFERENCES,
    SAMPLE_ID_PATTERN,
)
from .inputs import (
    _barcode_map_digest,
    droplet_whitelist_path,
    load_barcode_config,
    load_droplet_whitelist,
)


def droplet_geometry() -> str:
    """返回 DD-MET5 的物理 R1 条码和 UMI、物理 R2 表达序列布局。"""
    return '1{b[17]u[12]x:}2{r:}'


def preprocess_droplet_rna(r1: Path, r2: Path, output_r1: Path, output_r2: Path, qc: Path) -> dict:
    """按官方次序调用 Cutadapt 匹配器，配对流式清理 RNA 并记录互斥的 read-pair 去向。"""
    import cutadapt
    import dnaio
    from cutadapt.adapters import BackAdapter, RightmostFrontAdapter
    adapters = [BackAdapter(seq, max_errors=0.1, min_overlap=3) if end == '3'
                else RightmostFrontAdapter(seq, max_errors=0.1, min_overlap=3)
                for seq, end in DROPLET_R2_ADAPTERS]
    destinations = (output_r1, output_r2, qc)
    for path in destinations:
        path.parent.mkdir(parents=True, exist_ok=True)
    metrics = dict(input_pairs=0, short_cb_umi=0, invalid_cb=0, invalid_umi=0,
                   short_r2=0, passed_pairs=0, adapter_trimmed_pairs=0)
    hits = [0] * len(adapters)
    started = time.monotonic()
    with tempfile.TemporaryDirectory(prefix='.droplet-trim.', dir=output_r1.parent) as tmp:
        a, b = Path(tmp) / 'R1.fq.gz', Path(tmp) / 'R2.fq.gz'
        with dnaio.open(r1, r2) as reads, dnaio.open(a, b, mode='w', compression_level=1) as writer:
            for first, second in reads:
                metrics['input_pairs'] += 1
                if len(first) < 29:
                    metrics['short_cb_umi'] += 1
                    continue
                if set(first.sequence[:17]) - set('ACGT'):
                    metrics['invalid_cb'] += 1
                    continue
                if set(first.sequence[17:29]) - set('ACGT'):
                    metrics['invalid_umi'] += 1
                    continue
                trimmed = False
                for i, adapter in enumerate(adapters):
                    match = adapter.match_to(second.sequence)
                    if match is not None:
                        second = match.trimmed(second)
                        hits[i] += 1
                        trimmed = True
                metrics['adapter_trimmed_pairs'] += int(trimmed)
                if len(second) < 60:
                    metrics['short_r2'] += 1
                    continue
                writer.write(first[:29], second)
                metrics['passed_pairs'] += 1
        if metrics['input_pairs'] != sum(metrics[k] for k in ('short_cb_umi','invalid_cb','invalid_umi','short_r2','passed_pairs')):
            raise ValueError('DD-MET5 preprocessing pair accounting failed')
        metrics.update(protocol='droplet', chemistry='DD-MET5', cutadapt_version=cutadapt.__version__,
                       source_commit=DROPLET_SOURCE_COMMIT, adapters=list(DROPLET_R2_ADAPTERS),
                       adapter_hits=hits, max_error_rate=0.1, min_overlap=3, min_r2_length=60,
                       elapsed_seconds=time.monotonic()-started)
        candidate_qc = Path(tmp) / 'qc.json'
        candidate_qc.write_text(json.dumps(metrics, indent=2)+'\n')
        qc.unlink(missing_ok=True)
        os.replace(a, output_r1)
        os.replace(b, output_r2)
        os.replace(candidate_qc, qc)
    return metrics


def preprocess_droplet_main(argv=None):
    """执行文库级 RNA 预处理；QC 最后提交，失败不得留下完成标记。"""
    parser = argparse.ArgumentParser()
    for arg in ('r1','r2','output-r1','output-r2','qc'):
        parser.add_argument('--'+arg, type=Path, required=True)
    args = parser.parse_args(argv)
    preprocess_droplet_rna(args.r1, args.r2, args.output_r1, args.output_r2, args.qc)
    return 0


EMPTY_DROPS_R = r'''
args <- commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages(library(DropletUtils))
stopifnot(packageVersion("DropletUtils") == "1.30.0")
x <- as(Matrix::readMM(args[1]), "CsparseMatrix")
rownames(x) <- paste0("g", seq_len(nrow(x)))
colnames(x) <- readLines(args[2])
stopifnot(ncol(x) == length(colnames(x)))
set.seed(1)
e <- emptyDrops(x, lower=100, niters=10000, round=TRUE, BPPARAM=BiocParallel::SerialParam())
called <- !is.na(e$FDR) & e$FDR <= 0.001
write.table(data.frame(barcode=colnames(x), as.data.frame(e), called=called),
    args[3], sep="\t", quote=FALSE, row.names=FALSE, na="NA")
m <- S4Vectors::metadata(e)
jsonlite::write_json(list(method="DropletUtils::emptyDrops", version="1.30.0",
    lower=100L, niters=10000L, seed=1L, fdr=0.001, round=TRUE,
    auto_retain_threshold=if(is.null(m$retain)) NULL else m$retain,
    automatically_retained_cells=sum(e$Total >= m$retain, na.rm=TRUE),
    tested=sum(!is.na(e$PValue)), limited=sum(e$Limited, na.rm=TRUE),
    matrix_barcodes=ncol(x), called_cells=sum(called),
    metadata=capture.output(str(m)), session=capture.output(sessionInfo())),
    args[4], auto_unbox=TRUE, pretty=TRUE, null="null")
'''


def call_droplet_rna(input_h5ad: Path, calls: Path, qc: Path, rscript: Path) -> None:
    """仅对临时稀疏矩阵运行标准取整空液滴检验；原始 EM 表达矩阵保持只读。"""
    import anndata as ad
    import hdf5plugin
    import scipy.sparse as sp
    from scipy.io import mmwrite
    for path in (calls, qc):
        path.parent.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(input_h5ad)
    if RRNA_COMPETITOR_ID in adata.var_names:
        raise ValueError('Use prepare-rna-counts to separate the 45S technical feature before cell calling')
    if not adata.obs_names.is_unique:
        raise ValueError('Duplicate RNA barcode in raw quantification')
    whitelist = load_droplet_whitelist(droplet_whitelist_path())
    if set(map(str, adata.obs_names)) - whitelist:
        raise ValueError('Raw RNA quantification contains a barcode outside the design whitelist')
    matrix = sp.csc_matrix(adata.X.T)
    if matrix.data.size == 0 or not np.isfinite(matrix.data).all() or (matrix.data < 0).any():
        raise ValueError('emptyDrops requires a finite nonnegative, nonempty background matrix')
    started = time.monotonic()
    with tempfile.TemporaryDirectory(prefix='.emptydrops.', dir=calls.parent) as tmp:
        folder = Path(tmp)
        mmwrite(folder/'counts.mtx', matrix, precision=17)
        (folder/'barcodes.txt').write_text('\n'.join(map(str, adata.obs_names))+'\n')
        (folder/'call.R').write_text(EMPTY_DROPS_R)
        try:
            subprocess.run([str(rscript), '--vanilla', str(folder/'call.R'), str(folder/'counts.mtx'),
                            str(folder/'barcodes.txt'), str(folder/'calls.tsv'), str(folder/'qc.json')], check=True)
            frame = pd.read_csv(folder/'calls.tsv', sep='\t', dtype={'barcode': str})
            if frame['barcode'].tolist() != list(map(str, adata.obs_names)):
                raise ValueError('emptyDrops barcode order or set changed')
            payload = json.loads((folder/'qc.json').read_text())
            payload['elapsed_seconds'] = time.monotonic()-started
            payload['expected_count_sum'] = float(matrix.sum())
            payload['status'] = 'PASS' if payload['called_cells'] else 'NO_CALLED_CELLS'
            (folder/'qc.json').write_text(json.dumps(payload, indent=2)+'\n')
            qc.unlink(missing_ok=True)
            os.replace(folder/'calls.tsv', calls)
            os.replace(folder/'qc.json', qc)
            if not payload['called_cells']:
                raise ValueError('emptyDrops called no cells; no threshold fallback is permitted')
        except BaseException as exc:
            failure = {'status':'FAILED', 'error':str(exc), 'lower':100, 'niters':10000, 'seed':1,
                       'elapsed_seconds':time.monotonic()-started}
            qc.with_suffix('.failed.json').write_text(json.dumps(failure, indent=2)+'\n')
            raise


def call_droplet_rna_main(argv=None):
    """执行独立 RNA 细胞调用，不读取 DNA 结果或文库对应表。"""
    parser = argparse.ArgumentParser()
    for arg in ('input', 'calls', 'qc', 'rscript'):
        parser.add_argument('--'+arg, type=Path, required=True)
    args = parser.parse_args(argv)
    call_droplet_rna(args.input, args.calls, args.qc, args.rscript)
    return 0


def build_rna_cell_metadata(
    barcode_map: str | Path,
    canonical_barcodes: Sequence[str],
    *,
    sample_id: str,
    rna_sample: str,
) -> list[dict[str, str]]:
    """Map canonical RNA barcodes to stable project cell identities.

    The scientific identity contract is ``<sample_id>_<PlateID>``.  Canonical
    RNA barcodes remain explicit metadata so the simpleaf/alevin-fry output is
    always traceable back to the project Barcode_Map.csv.
    """
    if not SAMPLE_ID_PATTERN.fullmatch(str(sample_id)):
        raise ValueError(f"invalid sample_id for RNA finalization: {sample_id!r}")
    config = load_barcode_config(barcode_map, require_rna=True)
    by_barcode = {row.rna_barcode: row for row in config.rows}
    metadata: list[dict[str, str]] = []
    missing: list[str] = []
    seen_cell_ids: dict[str, str] = {}

    for raw_barcode in canonical_barcodes:
        barcode = str(raw_barcode)
        row = by_barcode.get(barcode)
        if row is None:
            missing.append(barcode)
            continue
        cell_id = f"{sample_id}_{row.plate_id}"
        previous = seen_cell_ids.get(cell_id)
        if previous is not None and previous != barcode:
            raise ValueError(
                "RNA project cell ID collision: "
                f"{previous!r} and {barcode!r} both map to {cell_id!r}. "
                "<sample_id>_<PlateID> requires PlateID to be unique among retained RNA cells."
            )
        seen_cell_ids[cell_id] = barcode
        metadata.append(
            {
                "cell_id": cell_id,
                "RNA_Barcode": barcode,
                "PlateID": row.plate_id,
                "Cell_Order": row.cell_order,
                "sample_id": str(sample_id),
                "rna_sample": str(rna_sample),
                # ``sample_id`` is the project-facing biological sample name.
                # Keep the raw FASTQ/library identifier separately in
                # ``rna_sample`` for provenance.
                "library": str(sample_id),
            }
        )

    if missing:
        preview = ", ".join(repr(x) for x in missing[:5])
        raise ValueError(
            "alevin-fry returned canonical RNA barcode(s) absent from Barcode_Map.csv: "
            f"{preview}" + (" ..." if len(missing) > 5 else "")
        )
    if len(metadata) != len(canonical_barcodes):
        raise ValueError("RNA finalization metadata length mismatch")
    return metadata


def rrna_reference_spec(species: str) -> dict[str, str]:
    """返回固定物种的45S竞争设计；不把人类序列静默应用于其它物种。"""
    if species not in RRNA_REFERENCES:
        raise ValueError(f'No pinned 45S competitor for RNA reference {species!r}; supported: hg38, mm10')
    accession, digest = RRNA_REFERENCES[species]
    return dict(species=species, accession=accession, feature_id=RRNA_COMPETITOR_ID,
                fasta_sha256=digest, source=f'https://www.ncbi.nlm.nih.gov/nuccore/{accession}')


def prepare_rrna_reference(species: str, output: Path) -> None:
    """把源码内固定的RefSeq序列复制至reference候选，不改原始FASTA/GTF。"""
    spec = rrna_reference_spec(species)
    source = PIPELINE_ROOT / 'resources/rrna' / f'{species}_45s.fa'
    payload = source.read_bytes()
    if hashlib.sha256(payload).hexdigest() != spec['fasta_sha256']:
        raise ValueError(f'Pinned 45S FASTA checksum mismatch: {source}')
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=output.parent, delete=False) as handle:
        temporary = Path(handle.name)
        handle.write(payload)
    try:
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)


def _validate_rrna_reference_files(reference: Path, species: str) -> dict[str, str]:
    spec = rrna_reference_spec(species)
    required = ['genome.fa', 'genes/genes.gtf', 'rrna_competitor.fa', 'simpleaf/ref/t2g.tsv',
                'simpleaf/index/piscem_idx.ctab', 'simpleaf/index/piscem_idx.ectab',
                'simpleaf/index/piscem_idx.refinfo', 'simpleaf/index/gene_id_to_name.tsv']
    missing = [name for name in required if not (reference/name).is_file() or not (reference/name).stat().st_size]
    if missing:
        raise ValueError('PG+45S RNA reference incomplete: '+', '.join(missing))
    if hashlib.sha256((reference/'rrna_competitor.fa').read_bytes()).hexdigest() != spec['fasta_sha256']:
        raise ValueError('RNA reference 45S competitor differs from the pinned species sequence')
    found = 0
    with (reference/'simpleaf/ref/t2g.tsv').open() as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 2 or not all(fields):
                raise ValueError('PG+45S requires a nonempty two-column transcript-to-gene map')
            if RRNA_COMPETITOR_ID in fields:
                if fields != [RRNA_COMPETITOR_ID, RRNA_COMPETITOR_ID]:
                    raise ValueError('45S technical transcript/gene identity conflict')
                found += 1
    if found != 1:
        raise ValueError('RNA reference must contain exactly one explicit 45S competitor mapping')
    return spec


def seal_rna_reference(reference: Path, species: str) -> None:
    """原生index成功后声明当前竞争设计；旧index不能仅靠现有文件被误认就绪。"""
    spec = _validate_rrna_reference_files(reference, species)
    # 完整扫描只发生在构建边界，普通项目启动不反复扫描数GB expanded FASTA。
    fasta = reference/'simpleaf/ref/roers_ref.fa'
    sequence = []; active = False; records = 0
    with fasta.open() as handle:
        for line in handle:
            if line.startswith('>'):
                active = line[1:].split()[0] == RRNA_COMPETITOR_ID
                records += int(active)
            elif active:
                sequence.append(line.strip())
    pinned = ''.join((reference/'rrna_competitor.fa').read_text().splitlines()[1:])
    if records != 1 or ''.join(sequence) != pinned:
        raise ValueError('Native expanded reference did not retain the exact pinned 45S sequence')
    metadata = {**spec, 'schema_version': 1, 'resolution': RNA_RESOLUTION,
                'technical_feature_excluded_from_expression': True}
    output = reference/'simpleaf/ref/rrna_competitor.json'
    with tempfile.NamedTemporaryFile(mode='w', dir=output.parent, delete=False) as handle:
        temporary = Path(handle.name); json.dump(metadata, handle, indent=2)
    try:
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)


def validate_rna_reference(reference: Path, species: str) -> dict[str, Any]:
    """Doctor和工作流共用的只读reference检查，拒绝没有竞争设计的旧参考。"""
    spec = _validate_rrna_reference_files(reference, species)
    marker = reference/'simpleaf/ref/rrna_competitor.json'
    if not marker.is_file():
        raise ValueError('PG+45S RNA reference is not sealed; rebuild with doctor.sh reference')
    metadata = json.loads(marker.read_text())
    expected = {**spec, 'schema_version': 1, 'resolution': RNA_RESOLUTION,
                'technical_feature_excluded_from_expression': True}
    if metadata != expected:
        raise ValueError('RNA reference competitor metadata disagrees with the current PG+45S policy')
    return metadata


def export_rna_counts(quant_dir: Path, gene_map: Path, output_h5ad: Path) -> None:
    """原生 MTX 只做格式导出；直接读取参考注释，避免消费失败的文件副本。"""
    import anndata as ad
    import scipy.sparse as sp

    quant = json.loads((quant_dir / 'quant.json').read_text())
    if quant.get('resolution_strategy') != 'ParsimonyGene' or quant.get('usa_mode') is not False:
        raise ValueError('RNA export requires two-column parsimony-gene quantification')
    folder = quant_dir / 'alevin'

    def axis(name):
        values = (folder / name).read_text().splitlines()
        if not values or any(not v or '\t' in v for v in values) or len(set(values)) != len(values):
            raise ValueError('RNA cell/gene identities must be nonempty unique single-column strings')
        return values

    barcodes = axis('quants_mat_rows.txt')
    genes = axis('quants_mat_cols.txt')
    names = pd.read_csv(gene_map, sep='\t', header=None, dtype=str, keep_default_na=False)
    if (names.shape[1] != 2 or (names == '').any().any() or names[0].duplicated().any()
            or set(names[0]) != set(genes)):
        raise ValueError('Reference gene map must match every native gene exactly once')
    names.columns = ['gene_id', 'gene_symbol']
    var = names.set_index('gene_id', drop=False).loc[genes].copy()
    var.index.name = 'index'

    with (folder / 'quants_mat.mtx').open() as handle:
        header = handle.readline().strip().lower().split()
        if (len(header) != 5 or header[:3] != ['%%matrixmarket', 'matrix', 'coordinate']
                or header[3] not in {'integer', 'real'} or header[4] != 'general'):
            raise ValueError('Expected general coordinate integer/real native Matrix Market counts')
        line = handle.readline()
        while line.startswith('%'):
            line = handle.readline()
        rows, cols, nnz = map(int, line.split())
        if ((rows, cols) != (len(barcodes), len(genes)) or nnz < 0
                or quant.get('num_genes') != cols or quant.get('num_quantified_cells') != rows):
            raise ValueError('Native quantification metadata, axes and matrix shape disagree')
        dtype = [('row', np.int64), ('col', np.int64),
                 ('value', np.int64 if header[3] == 'integer' else np.float64)]
        entries = np.loadtxt(handle, dtype=dtype, ndmin=1) if nnz else np.empty(0, dtype=dtype)
        if nnz == 0 and any(line.strip() and not line.startswith('%') for line in handle):
            raise ValueError('Native matrix has undeclared entries')
    if (len(entries) != nnz or (entries['row'] < 1).any() or (entries['row'] > rows).any()
            or (entries['col'] < 1).any() or (entries['col'] > cols).any()
            or not np.isfinite(entries['value']).all() or (entries['value'] < 0).any()):
        raise ValueError('Invalid native matrix coordinates, entry count or molecule values')
    matrix = sp.coo_matrix((entries['value'], (entries['row'] - 1, entries['col'] - 1)),
                           shape=(rows, cols)).tocsr()
    if matrix.nnz != nnz:
        raise ValueError('Duplicate native matrix coordinates must not be silently summed')

    rename = dict(CorrectedReads='corrected_reads', MappedReads='mapped_reads',
                  DeduplicatedReads='deduplicated_reads', MappingRate='mapping_rate',
                  DedupRate='dedup_rate', MeanByMax='mean_by_max',
                  NumGenesExpressed='num_genes_expressed', NumGenesOverMean='num_genes_over_mean')
    feature_types = {key: (np.int64 if key in {'CorrectedReads', 'MappedReads',
                     'NumGenesExpressed', 'NumGenesOverMean'} else np.float64) for key in rename}
    features = pd.read_csv(quant_dir / 'featureDump.txt', sep='\t',
                           dtype=str, keep_default_na=False)
    if list(features.columns) != ['CB', *rename] or features['CB'].tolist() != barcodes:
        raise ValueError('Native featureDump barcode identity/order or columns disagree')
    obs = features.drop(columns='CB').rename(columns=rename)
    for column in obs:
        # 原生空背景条码的比率可为 NaN；保留未定义值，不补零、不改分子计数。
        source_name = next(key for key, value in rename.items() if value == column)
        obs[column] = np.asarray(obs[column], dtype=feature_types[source_name])
    obs.index = pd.Index(barcodes, name='index')
    obs.insert(0, 'barcodes', barcodes)
    uns = {key: json.dumps(json.loads((quant_dir / file).read_text()), separators=(',', ':'))
           for key, file in [('gpl_info', 'generate_permit_list.json'), ('collate_info', 'collate.json'),
                             ('quant_info', 'quant.json'), ('simpleaf_map_info', 'simpleaf_map_info.json')]}
    data = ad.AnnData(matrix, obs=obs, var=var, uns=uns)
    output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix='.' + output_h5ad.name + '.', suffix='.h5ad', dir=output_h5ad.parent)
    os.close(fd)
    temporary = Path(name)
    try:
        data.write_h5ad(temporary)
        os.replace(temporary, output_h5ad)
    finally:
        temporary.unlink(missing_ok=True)


def export_rna_counts_main(argv=None):
    """内部格式导出入口；不执行 barcode correction 或分子重分配。"""
    parser = argparse.ArgumentParser()
    for key in ('quant-dir', 'gene-map', 'output'):
        parser.add_argument('--' + key, type=Path, required=True)
    args = parser.parse_args(argv)
    export_rna_counts(args.quant_dir, args.gene_map, args.output)
    return 0


def prepare_rna_counts(input_h5ad: Path, quant_json: Path, reference: Path,
                       species: str, output_h5ad: Path) -> None:
    """保存45S技术计数再投影生物gene；不取整，不更改任何细胞/生物gene计数。"""
    import anndata as ad
    import hdf5plugin
    import scipy.sparse as sp
    policy = validate_rna_reference(reference, species)
    quant = json.loads(quant_json.read_text())
    if quant.get('resolution_strategy') != 'ParsimonyGene' or quant.get('usa_mode') is not False:
        raise ValueError('Production RNA counts require two-column parsimony-gene quantification')
    data = ad.read_h5ad(input_h5ad)
    if not data.obs_names.is_unique or not data.var_names.is_unique:
        raise ValueError('RNA count cell/gene identities must be unique')
    if 'gene_id' in data.var and list(map(str, data.var['gene_id'])) != list(map(str, data.var_names)):
        raise ValueError('RNA var index and gene_id disagree')
    if RRNA_COMPETITOR_ID not in data.var_names:
        raise ValueError('Native RNA matrix lacks the explicit 45S competitor')
    matrix = sp.csr_matrix(data.X)
    if (not np.isfinite(matrix.data).all() or (matrix.data < 0).any()
            or (matrix.data != np.floor(matrix.data)).any()):
        raise ValueError('Native parsimony-gene molecule counts must be finite nonnegative integers; never round')
    matrix.eliminate_zeros()
    if quant.get('num_genes') != data.n_vars or quant.get('num_quantified_cells') != data.n_obs:
        raise ValueError('Native quantification metadata and matrix shape disagree')
    j = data.var_names.get_loc(RRNA_COMPETITOR_ID)
    sink = np.asarray(matrix[:,j].sum(axis=1)).ravel()
    keep = np.arange(data.n_vars) != j
    biological = data[:,keep].copy()
    biological.X = matrix[:,keep].copy()
    biological.obs['rna45s_competitor_counts'] = sink
    biological.uns['rna_quantification'] = {**policy, 'count_semantics': 'native unrounded non-EM molecule counts',
        'quantifier_version': str(quant.get('version_str', 'not_recorded')),
        'background_competitor_count_sum': float(sink.sum()),
        'background_competitor_cells_positive': int((sink > 0).sum())}
    output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix='.'+output_h5ad.name+'.', suffix='.h5ad', dir=output_h5ad.parent)
    os.close(fd); temporary = Path(name)
    try:
        biological.write_h5ad(temporary)
        os.replace(temporary, output_h5ad)
    finally:
        temporary.unlink(missing_ok=True)


def prepare_rna_counts_main(argv=None):
    """内部矩阵投影入口；原始定量与技术计数的审计含义保持明确。"""
    parser = argparse.ArgumentParser()
    for key in ['input', 'quant-json', 'reference', 'output']:
        parser.add_argument('--'+key, type=Path, required=True)
    parser.add_argument('--species', required=True)
    a = parser.parse_args(argv)
    prepare_rna_counts(a.input, a.quant_json, a.reference, a.species, a.output)
    return 0


def finalize_rna_h5ad(
    input_h5ad: str | Path,
    output_h5ad: str | Path,
    barcode_map: str | Path | None,
    *,
    sample_id: str,
    rna_sample: str,
    protocol: str = 'cabernet',
    calls: Path | None = None,
    calling_qc: Path | None = None,
    preprocess_qc: Path | None = None,
    permit_qc: Path | None = None,
    map_info: Path | None = None,
) -> Path:
    """Create the project-level RNA.h5ad from raw simpleaf quantification."""
    try:
        import anndata as ad
        import hdf5plugin  # noqa: F401 - registers simpleaf's Blosc HDF5 filter
        import pandas as pd
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "RNA finalization requires anndata, hdf5plugin, and pandas in the rules environment"
        ) from exc

    source = Path(input_h5ad)
    if not source.is_file():
        raise FileNotFoundError(f"simpleaf quants.h5ad does not exist: {source}")
    adata = ad.read_h5ad(source)
    if RRNA_COMPETITOR_ID in adata.var_names:
        raise ValueError('Use prepare-rna-counts to separate the 45S technical feature before finalization')
    canonical = [str(value) for value in adata.obs_names]
    if protocol == 'droplet':
        if barcode_map is not None or any(p is None for p in (calls, calling_qc, preprocess_qc, permit_qc, map_info)):
            raise ValueError('Droplet finalization requires independent RNA calls and QC, without Barcode Map')
        if not SAMPLE_ID_PATTERN.fullmatch(sample_id) or sample_id in {'.', '..'}:
            raise ValueError('Invalid Droplet RNA sample ID')
        frame = pd.read_csv(calls, sep='\t', dtype={'barcode': str})
        if frame['barcode'].tolist() != canonical or frame['barcode'].duplicated().any():
            raise ValueError('RNA calls disagree with quantification barcodes')
        if frame['called'].dtype != bool:
            raise ValueError('RNA calling mask must contain explicit booleans')
        selected = frame['called'].to_numpy()
        if not selected.any():
            raise ValueError('Cannot finalize a failed or empty RNA cell call')
        trim = json.loads(preprocess_qc.read_text())
        mapping = json.loads(map_info.read_text())
        permit = json.loads(permit_qc.read_text())
        correction = permit['correction_stats']
        barcode_reads = [correction[key+'_reads'] for key in ('exact', 'corrected', 'ambiguous', 'not_found')]
        if (any(type(n) is not int or n < 0 for n in barcode_reads)
                or mapping['num_reads'] != trim['passed_pairs']
                or not 0 <= sum(barcode_reads) <= mapping['num_mapped'] <= mapping['num_reads']):
            raise ValueError('Droplet trim/mapping/barcode read accounting failed')
        call_metrics = json.loads(calling_qc.read_text())
        if call_metrics.get('status') != 'PASS' or call_metrics['called_cells'] != int(selected.sum()):
            raise ValueError('Droplet calling QC disagrees with the selected cells')
        adata = adata[selected].copy()
        frame = frame.loc[selected].copy()
        canonical = list(map(str, adata.obs_names))
        counts = np.asarray(adata.X.sum(axis=1)).ravel()
        rank = {barcode: str(i+1) for i, (_, barcode) in enumerate(sorted(zip(-counts, canonical)))}
        metadata = [dict(cell_id=f'{sample_id}_{barcode}', RNA_Barcode=barcode, PlateID='',
                         Cell_Order=rank[barcode], sample_id=sample_id, rna_sample=rna_sample,
                         library=sample_id) for barcode in canonical]
        for key in ('Total', 'LogProb', 'PValue', 'FDR'):
            adata.obs['emptyDrops_'+key] = frame[key].to_numpy()
        adata.obs['emptyDrops_Limited'] = pd.array(frame['Limited'], dtype='boolean')
        adata.obs['rna_called'] = True
        adata.obs['protocol'] = 'droplet'
    else:
        metadata = build_rna_cell_metadata(
            barcode_map, canonical, sample_id=sample_id, rna_sample=rna_sample
        )

    # Preserve all simpleaf/alevin-fry obs columns and add the project identity
    # columns in a deterministic order.
    for column in ("RNA_Barcode", "PlateID", "Cell_Order", "sample_id", "rna_sample", "library"):
        adata.obs[column] = [row[column] for row in metadata]
    adata.obs_names = pd.Index([row["cell_id"] for row in metadata], name="cell_id")
    if not adata.obs_names.is_unique:
        raise ValueError("RNA.h5ad project cell IDs are not unique")

    # 保留原生分子计数和≥1支持度审计；任何消费端都不通过取整制造观测。
    import scipy.sparse as sp

    matrix = sp.csr_matrix(adata.X)
    if matrix.data.size and (
        not np.isfinite(matrix.data).all() or np.any(matrix.data < 0)
    ):
        raise ValueError("RNA count matrix contains non-finite or negative values")
    total_expected = np.asarray(matrix.sum(axis=1)).ravel()
    positive_genes = np.diff(matrix.indptr).astype(np.int64)
    supported_genes = np.asarray((matrix >= 1.0).sum(axis=1)).ravel().astype(np.int64)
    rows = np.repeat(np.arange(matrix.shape[0]), np.diff(matrix.indptr))
    fractional = ~np.isclose(matrix.data, np.rint(matrix.data), atol=1e-6)
    subunit = (matrix.data > 0) & (matrix.data < 1.0)
    fractional_nnz = np.bincount(rows, weights=fractional, minlength=matrix.shape[0])
    subunit_mass = np.bincount(
        rows, weights=np.where(subunit, matrix.data, 0.0), minlength=matrix.shape[0]
    )
    adata.obs["total_expected_counts"] = total_expected
    adata.obs["subunit_expected_counts"] = subunit_mass
    adata.obs["total_supported_expected_counts"] = total_expected - subunit_mass
    adata.obs["n_genes_positive_expected"] = positive_genes
    adata.obs["n_genes_supported"] = supported_genes
    adata.obs["n_genes_fractional_only"] = positive_genes - supported_genes
    adata.obs["fractional_nnz_fraction"] = np.divide(
        fractional_nnz, positive_genes,
        out=np.zeros(matrix.shape[0], dtype=float), where=positive_genes > 0,
    )
    adata.obs["subunit_count_mass_fraction"] = np.divide(
        subunit_mass, total_expected,
        out=np.zeros(matrix.shape[0], dtype=float), where=total_expected > 0,
    )

    provenance = dict(adata.uns.get("rna_darlin_pipeline", {}))
    if protocol == 'droplet':
        provenance.update(protocol='droplet', chemistry='DD-MET5', geometry=droplet_geometry(),
                          expected_orientation='rc', whitelist_sha256=DROPLET_WHITELIST_SHA256,
                          whitelist_source=DROPLET_WHITELIST_URL, whitelist_rows=DROPLET_WHITELIST_ROWS,
                          calling=call_metrics, preprocessing=trim, barcode_correction=correction,
                          barcode_correction_version=permit['version_str'],
                          mapping={k:mapping[k] for k in ('num_reads','num_mapped')})
        barcode_provenance = {'cell_id_format': '<sample_id>_<barcode>'}
    else:
        barcode_map_sha256, barcode_map_rows = _barcode_map_digest(barcode_map)
        barcode_provenance = dict(cell_id_format='<sample_id>_<PlateID>',
            barcode_map_path=str(Path(barcode_map).expanduser().resolve()),
            barcode_map_sha256=barcode_map_sha256, barcode_map_rows=int(barcode_map_rows))
    provenance.update(
        {
            "sample_id": str(sample_id),
            "rna_sample": str(rna_sample),
            **barcode_provenance,
            "canonical_barcode_column": "RNA_Barcode",
            "source_quants_h5ad": str(source),
            "gene_detection_primary_metric": "n_genes_supported",
            "gene_detection_expected_count_threshold": 1.0,
            "gene_detection_raw_positive_metric": "n_genes_positive_expected",
        }
    )
    adata.uns["rna_darlin_pipeline"] = provenance
    if 'rna45s_competitor_counts' in adata.obs:
        adata.uns['rna_quantification']['finalized_competitor_count_sum'] = float(adata.obs['rna45s_competitor_counts'].sum())
        adata.uns['rna_quantification']['finalized_competitor_cells_positive'] = int((adata.obs['rna45s_competitor_counts'] > 0).sum())

    destination = Path(output_h5ad)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temp_path = destination.with_name(destination.name + ".tmp.h5ad")
    try:
        adata.write_h5ad(temp_path)
        os.replace(temp_path, destination)
    finally:
        if temp_path.exists():
            temp_path.unlink()
    return destination


def finalize_rna_h5ad_main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Finalize simpleaf quants.h5ad into project-level RNA.h5ad."
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--barcode-map", type=Path)
    parser.add_argument("--protocol", choices=sorted(PROTOCOLS), required=True)
    parser.add_argument("--calls", type=Path)
    parser.add_argument("--calling-qc", type=Path)
    parser.add_argument("--preprocess-qc", type=Path)
    parser.add_argument("--permit-qc", type=Path)
    parser.add_argument("--map-info", type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--rna-sample", required=True)
    args = parser.parse_args(argv)
    output = finalize_rna_h5ad(
        args.input,
        args.output,
        args.barcode_map,
        sample_id=args.sample_id,
        rna_sample=args.rna_sample,
        protocol=args.protocol, calls=args.calls, calling_qc=args.calling_qc,
        preprocess_qc=args.preprocess_qc,
        permit_qc=args.permit_qc, map_info=args.map_info,
    )
    print(f"RNA project matrix finalized: {output}")
    return 0
