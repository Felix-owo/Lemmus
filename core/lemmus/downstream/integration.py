"""当前sealed DNA交付与RNA身份关联；不猜测文库或细胞映射。"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import pandas as pd
from .config import (
    DownstreamQCError,
)


def _identity_frame(frame: pd.DataFrame, columns: tuple[str, ...], label: str) -> pd.DataFrame:
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise DownstreamQCError(f"{label} lacks required identity columns: {missing}")
    result = frame.copy()
    for column in columns:
        if result[column].isna().any():
            raise DownstreamQCError(f"{label} has missing {column}")
        result[column] = result[column].astype(str).str.strip()
        if result[column].eq("").any():
            raise DownstreamQCError(f"{label} has empty {column}")
    return result


def _match_cells_by_clone_plate(keys: pd.DataFrame, rna: pd.DataFrame) -> dict[str, str]:
    pairs = list(zip(rna["CloneID"], rna["PlateID"]))
    if len(pairs) != len(set(pairs)):
        raise DownstreamQCError("RNA (CloneID, PlateID) is not unique")
    index = dict(zip(pairs, rna["Sample_ID"]))
    mapping = {}
    for row in keys.itertuples(index=False):
        match = index.get((row.CloneID, row.PlateID))
        if match is None:
            raise DownstreamQCError(
                f"Cannot uniquely align DNA cell {row.Sample_ID!r} to RNA "
                f"(CloneID={row.CloneID!r}, PlateID={row.PlateID!r}); "
                "check the sealed DNA rna_sample link and RNA.h5ad identities."
            )
        mapping[row.Sample_ID] = match
    if len(set(mapping.values())) != len(mapping):
        raise DownstreamQCError("Multiple DNA cells map to the same RNA cell")
    if set(mapping.values()) != set(rna["Sample_ID"]):
        raise DownstreamQCError("RNA cells have no DNAme counterpart; both QC tables must cover the same cells")
    return mapping


def resolve_droplet_library_map(dna_manifest, rna_identity):
    """只从 sealed DNA manifest 解析可选 RNA raw 关联；留空的文库保持未关联。"""
    if dna_manifest is None or rna_identity is None:
        raise DownstreamQCError('Droplet linking requires sealed DNA manifest and RNA.h5ad identity')
    dna = _identity_frame(dna_manifest, ('project_sample_id',), 'DNA sealed manifest')
    if 'rna_sample' not in dna:
        raise DownstreamQCError('DNA sealed manifest lacks rna_sample')
    dna['rna_sample'] = dna['rna_sample'].fillna('').astype(str).str.strip()
    pairs = dna[['project_sample_id', 'rna_sample']].drop_duplicates()
    if pairs['project_sample_id'].duplicated().any():
        raise DownstreamQCError('DNA sealed manifest has inconsistent RNA links within a library')
    pairs = pairs.loc[pairs.rna_sample.ne('')].rename(columns={'project_sample_id': 'dna_project_sample_id'})
    rna = _identity_frame(rna_identity, ('sample_id', 'rna_sample'), 'RNA identity')
    rna = rna[['sample_id', 'rna_sample']].drop_duplicates()
    if any(rna[c].duplicated().any() for c in rna):
        raise DownstreamQCError('RNA canonical/raw library identities must be one-to-one')
    if pairs.rna_sample.duplicated().any():
        raise DownstreamQCError('DNA/RNA library links must be one-to-one')
    if set(pairs.rna_sample) - set(rna.rna_sample):
        raise DownstreamQCError('DNA rna_sample references a library absent from finalized RNA results')
    result = pairs.sort_values('dna_project_sample_id').reset_index(drop=True)
    result.attrs['source'] = {'mode': 'sealed_dna_manifest',
        'sha256': hashlib.sha256(result.to_csv(index=False).encode()).hexdigest()}
    return result


def validate_droplet_integration_links(summary, dna_manifest):
    """续跑读取 Raw/Single 前，核对整合证明确实绑定当前 sealed 文库关联。"""
    if not isinstance(summary.get('library_links'), list):
        raise DownstreamQCError('整合缺少 sealed 文库关联证明，请重新运行整合')
    columns = ('dna_project_sample_id', 'rna_sample', 'rna_project_sample_id', 'protocol')
    links = _identity_frame(pd.DataFrame(summary['library_links'], columns=columns), columns, 'Integration library links')
    if not links.protocol.eq('droplet').all():
        raise DownstreamQCError('整合文库协议不一致，请重新运行整合')
    rna = links.rename(columns={'rna_project_sample_id': 'sample_id'})
    expected = resolve_droplet_library_map(dna_manifest, rna)
    pairs = links[['dna_project_sample_id', 'rna_sample']].sort_values('dna_project_sample_id')
    if (summary.get('library_map') != expected.attrs['source']
            or pairs.to_dict('records') != expected.to_dict('records')):
        raise DownstreamQCError('整合文库关联与当前 sealed manifest 不一致，请重新运行整合')


def integrate_droplet_qc(rna_qc, dna_qc, dna_manifest, rna_identity):
    """按明确文库和完整条码关联独立细胞名单，保留缺失和单模态结果的对应清单。"""
    import re
    if rna_identity is None or dna_manifest is None:
        raise DownstreamQCError('Droplet linking requires RNA.h5ad identity, sealed DNA manifest with optional rna_sample links')
    rna = _identity_frame(rna_qc, ('Sample_ID', 'CloneID'), 'RNA QC')
    dna = _identity_frame(dna_qc, ('Sample_ID', 'CloneID'), 'DNA QC')
    rm = _identity_frame(rna_identity, ('Sample_ID', 'sample_id', 'rna_sample', 'RNA_Barcode'), 'RNA.h5ad identity')
    dm = _identity_frame(dna_manifest, ('sample_id', 'project_sample_id', 'protocol', 'dna_barcode'), 'DNA sealed manifest')
    links = resolve_droplet_library_map(dm, rm)
    if set(dm['protocol']) != {'droplet'}:
        raise DownstreamQCError('Explicit Droplet library mapping requires a Droplet DNA delivery')
    for frame, col, label in ((rna,'Sample_ID','RNA QC'), (dna,'Sample_ID','DNA QC'),
                              (rm,'Sample_ID','RNA identity'), (dm,'sample_id','DNA manifest')):
        if frame[col].duplicated().any():
            raise DownstreamQCError(f'{label} has duplicate cell IDs')
    for frame, col in ((rna,'PlateID'), (dna,'PlateID'), (rm,'PlateID'), (dm,'plate_id')):
        if col not in frame or frame[col].fillna('').astype(str).ne('').any():
            raise DownstreamQCError('Droplet cells require empty PlateID')
        frame[col] = ''
    for frame, lib, bc in ((rm,'rna_sample','RNA_Barcode'), (dm,'project_sample_id','dna_barcode')):
        if frame[[lib,bc]].duplicated().any() or not frame[bc].map(lambda x: bool(re.fullmatch('[AGT]{17}',x))).all():
            raise DownstreamQCError('Droplet canonical barcodes must be unique AGT17 within each library')
    raw_libraries = rm[['rna_sample','sample_id']].drop_duplicates()
    if any(raw_libraries[c].duplicated().any() for c in raw_libraries):
        raise DownstreamQCError('RNA raw library and project sample mapping is not one-to-one')
    raw_to_sample = dict(zip(raw_libraries['rna_sample'],raw_libraries['sample_id']))
    for qc, manifest, idcol, clonecol, label in ((rna,rm,'Sample_ID','sample_id','RNA'),(dna,dm,'sample_id','project_sample_id','DNA')):
        known = manifest.set_index(idcol)
        if not set(qc['Sample_ID']).issubset(known.index):
            raise DownstreamQCError(f'{label} QC contains cells absent from finalized results')
        if qc['CloneID'].tolist() != known.loc[qc['Sample_ID'],clonecol].tolist():
            raise DownstreamQCError(f'{label} QC library identity disagrees with finalized results')
    by_dna = dict(zip(links['dna_project_sample_id'],links['rna_sample']))
    by_rna = {v:k for k,v in by_dna.items()}
    dna_ids, rna_ids = set(dna['Sample_ID']), set(rna['Sample_ID'])
    rna_index = {(r.rna_sample,r.RNA_Barcode): r for r in rm.itertuples(index=False)}
    records, consumed = [], set()
    for d in dm.itertuples(index=False):
        raw = by_dna.get(d.project_sample_id, '')
        r = rna_index.get((raw,d.dna_barcode))
        if r is not None:
            consumed.add(r.Sample_ID)
        dq, rq = d.sample_id in dna_ids, r is not None and r.Sample_ID in rna_ids
        missing = ([] if dq else ['dna_result_missing']) + ([] if r is None or rq else ['rna_qc_missing'])
        records.append(dict(dna_project_sample_id=d.project_sample_id, rna_sample=raw,
            rna_project_sample_id=raw_to_sample.get(raw,''), barcode=d.dna_barcode,
            dna_cell_id=d.sample_id, rna_cell_id='' if r is None else r.Sample_ID,
            membership='unmapped_dna_library' if not raw else ('both' if r is not None else 'dna_only'),
            result_status=';'.join(missing) or 'available', dna_qc_available=dq, rna_qc_available=rq,
            dna_analysis_status=getattr(d,'analysis_status',''),
            rna_qc_pass=getattr(r,'qc_pass',None) if r is not None else None))
    for r in rm.itertuples(index=False):
        if r.Sample_ID in consumed:
            continue
        dq = by_rna.get(r.rna_sample,'')
        rq = r.Sample_ID in rna_ids
        records.append(dict(dna_project_sample_id=dq, rna_sample=r.rna_sample,
            rna_project_sample_id=r.sample_id, barcode=r.RNA_Barcode,
            dna_cell_id='', rna_cell_id=r.Sample_ID,
            membership='rna_only' if dq else 'unmapped_rna_library',
            result_status='available' if rq else 'rna_qc_missing',
            dna_qc_available=False, rna_qc_available=rq, dna_analysis_status='',
            rna_qc_pass=getattr(r,'qc_pass',None)))
    ledger = pd.DataFrame(records).sort_values(['dna_project_sample_id','rna_sample','barcode']).reset_index(drop=True)
    common = ledger.loc[ledger['membership'].eq('both') & ledger['dna_qc_available'] & ledger['rna_qc_available']]
    rpart = rna.set_index('Sample_ID').loc[common['rna_cell_id']].reset_index()
    dpart = dna.set_index('Sample_ID').loc[common['dna_cell_id']].reset_index()
    dpart['Sample_ID'] = common['rna_cell_id'].to_numpy()
    dpart['CloneID'] = rpart['CloneID'].to_numpy()
    result = _merge_on_ids(rpart,dpart)
    result.attrs['cell_links'] = ledger
    result.attrs['library_links'] = [dict(dna_project_sample_id=k, rna_sample=v,
        rna_project_sample_id=raw_to_sample[v], protocol='droplet') for k,v in sorted(by_dna.items())]
    result.attrs['library_map'] = links.attrs['source']
    return result


def integrate_dna_qc(
    rna_qc: pd.DataFrame,
    dna_qc: pd.DataFrame,
    dna_manifest: pd.DataFrame | None = None,
    *,
    rna_identity: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """按已发布 DNA 关联、RNA.h5ad 原始文库身份和 PlateID 一对一整合，保留 RNA 侧身份。"""
    is_droplet = dna_manifest is not None and 'protocol' in dna_manifest and dna_manifest['protocol'].eq('droplet').any()
    if is_droplet:
        return integrate_droplet_qc(rna_qc, dna_qc, dna_manifest, rna_identity)
    identity = ("Sample_ID", "CloneID", "PlateID")
    rna = _identity_frame(rna_qc, identity, "RNA QC table")
    dna = _identity_frame(dna_qc, identity, "DNAme QC table")
    for frame, label in ((rna, "RNA"), (dna, "DNAme")):
        if frame["Sample_ID"].duplicated().any():
            raise DownstreamQCError(f"{label} Sample_ID is not unique")
        if frame[["CloneID", "PlateID"]].duplicated().any():
            raise DownstreamQCError(f"{label} (CloneID, PlateID) is not unique")
    if dna_manifest is None:
        if set(rna["Sample_ID"]) != set(dna["Sample_ID"]):
            raise DownstreamQCError("Sample_ID sets differ; provide the sealed DNA sample_manifest.tsv")
        return _merge_on_ids(rna, dna)

    manifest = _identity_frame(
        dna_manifest, ("sample_id", "project_sample_id", "plate_id", "protocol"),
        "DNA sealed sample_manifest.tsv",
    )
    for column in ("dna_raw_sample", "rna_sample"):
        if column not in manifest.columns:
            raise DownstreamQCError(f"DNA sealed sample_manifest.tsv lacks required column: {column}")
        manifest[column] = manifest[column].fillna("").astype(str).str.strip()
    if manifest["sample_id"].duplicated().any():
        raise DownstreamQCError("DNA sealed sample_manifest.tsv sample_id is not unique")
    if not set(manifest["protocol"]).issubset({"cabernet", "taps", "srd"}):
        raise DownstreamQCError("DNA sealed sample_manifest.tsv has an unsupported protocol")
    libraries = manifest[["project_sample_id", "dna_raw_sample", "rna_sample", "protocol"]].drop_duplicates()
    if libraries["project_sample_id"].duplicated().any():
        raise DownstreamQCError("DNA manifest maps one project_sample_id to multiple raw-sample identities")
    nonempty_dna = libraries.loc[libraries["dna_raw_sample"].ne("")]
    if nonempty_dna["dna_raw_sample"].duplicated().any():
        raise DownstreamQCError("DNA manifest maps one dna_raw_sample to multiple project samples")
    linked = libraries.loc[libraries["protocol"].isin(["cabernet", "taps"]) & libraries["rna_sample"].ne("")]
    if linked["rna_sample"].duplicated().any():
        raise DownstreamQCError("DNA manifest maps multiple DNA libraries to one external rna_sample")
    selected = manifest.set_index("sample_id")
    if not set(dna["Sample_ID"]).issubset(selected.index):
        raise DownstreamQCError("DNAme cells are absent from the sealed DNA sample_manifest.tsv")
    selected = selected.loc[dna["Sample_ID"]]
    for qc_column, manifest_column in (("CloneID", "project_sample_id"), ("PlateID", "plate_id")):
        if dna[qc_column].tolist() != selected[manifest_column].tolist():
            raise DownstreamQCError(f"DNAme {qc_column} disagrees with the sealed DNA manifest")

    libraries = libraries.loc[libraries["project_sample_id"].isin(selected["project_sample_id"])]
    linked = linked.loc[linked["project_sample_id"].isin(selected["project_sample_id"])]
    targets = dict(zip(libraries["project_sample_id"], libraries["project_sample_id"]))
    rna_cells = None
    if rna_identity is not None or len(linked):
        source = rna_identity if rna_identity is not None else rna_qc
        rna_cells = _identity_frame(
            source, ("Sample_ID", "sample_id", "rna_sample", "PlateID"), "RNA.h5ad identity",
        )
        if rna_cells["Sample_ID"].duplicated().any() or set(rna_cells["Sample_ID"]) != set(rna["Sample_ID"]):
            raise DownstreamQCError("RNA.h5ad identity and RNA QC table must contain the same unique cells")
        rna_cells = rna_cells.set_index("Sample_ID").loc[rna["Sample_ID"]]
        if rna_cells["sample_id"].tolist() != rna["CloneID"].tolist() or rna_cells["PlateID"].tolist() != rna["PlateID"].tolist():
            raise DownstreamQCError("RNA QC CloneID/PlateID disagrees with RNA.h5ad identity")
        raw_libraries = rna_cells[["rna_sample", "sample_id"]].drop_duplicates()
        if raw_libraries["rna_sample"].duplicated().any() or raw_libraries["sample_id"].duplicated().any():
            raise DownstreamQCError("RNA.h5ad raw library and sample_id mapping is not one-to-one")
        raw_to_sample = dict(zip(raw_libraries["rna_sample"], raw_libraries["sample_id"]))
        for row in linked.itertuples(index=False):
            if row.rna_sample not in raw_to_sample:
                raise DownstreamQCError(f"External RNA raw library {row.rna_sample!r} is absent from RNA.h5ad")
            targets[row.project_sample_id] = raw_to_sample[row.rna_sample]
    if len(set(targets.values())) != len(targets):
        raise DownstreamQCError("Multiple DNA project samples map to one RNA sample_id")
    keys = pd.DataFrame({
        "Sample_ID": dna["Sample_ID"].to_numpy(),
        "CloneID": selected["project_sample_id"].map(targets).to_numpy(),
        "PlateID": selected["plate_id"].to_numpy(),
    })
    mapping = _match_cells_by_clone_plate(keys, rna)
    if rna_cells is not None:
        for source_column, target_column in (("rna_barcode", "RNA_Barcode"), ("cell_order", "Cell_Order")):
            if source_column not in selected.columns or target_column not in rna_cells.columns:
                continue
            for sample, value in selected[source_column].items():
                value = str(value).strip()
                if value and value != str(rna_cells.loc[mapping[sample], target_column]).strip():
                    raise DownstreamQCError(f"DNA/RNA {target_column} disagrees for DNA cell {sample!r}")
    dna["Sample_ID"] = dna["Sample_ID"].map(mapping)
    dna["CloneID"] = keys["CloneID"].to_numpy()
    integrated = _merge_on_ids(rna, dna)
    integrated.attrs["library_links"] = [
        {"dna_project_sample_id": row.project_sample_id, "dna_raw_sample": row.dna_raw_sample,
         "external_rna_sample": row.rna_sample if row.protocol in {"cabernet", "taps"} else "",
         "rna_project_sample_id": targets[row.project_sample_id], "protocol": row.protocol}
        for row in libraries.itertuples(index=False)
    ]
    return integrated


def _merge_on_ids(rna: pd.DataFrame, dna: pd.DataFrame) -> pd.DataFrame:
    shared_identity = [column for column in ("CloneID", "PlateID", "RNA_Barcode") if column in rna.columns and column in dna.columns]
    for column in shared_identity:
        by_id = dna.set_index("Sample_ID")[column].reindex(rna["Sample_ID"])
        if rna[column].tolist() != by_id.tolist():
            raise DownstreamQCError(f"RNA and DNAme {column} disagree for Sample_ID")
    dna_part = dna.drop(columns=[c for c in ("CloneID", "PlateID", "RNA_Barcode") if c in dna.columns])
    collisions = sorted(set(dna_part.columns) & set(rna.columns) - {"Sample_ID"})
    if collisions:
        raise DownstreamQCError("RNA and DNAme QC tables share non-identity columns: " + ", ".join(collisions))
    integrated = rna.merge(dna_part, on="Sample_ID", how="inner", validate="one_to_one")
    if len(integrated) != len(rna) or len(integrated) != len(dna_part):
        raise DownstreamQCError("Integrated table row count does not match its inputs")
    return integrated


def load_dna_qc(dna_project_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """读取当前 complete 交付对应的 DNA QC 与已发布关联，不读取可变 config 或按 mtime 选旧代。"""
    import csv
    import io
    import re

    root = dna_project_dir.expanduser().resolve()
    marker_path = root / "03_results/run_manifest.json"
    try:
        marker_bytes = marker_path.read_bytes()
        marker = json.loads(marker_bytes)
    except (OSError, ValueError) as error:
        raise DownstreamQCError(f"DNA project has no readable complete run_manifest.json: {root}") from error
    if not isinstance(marker, dict) or not all(isinstance(marker.get(key), dict) for key in ("delivery", "run", "pipeline")):
        raise DownstreamQCError("DNA run_manifest.json has invalid delivery metadata")
    delivery = marker["delivery"].get("id", "")
    if (marker.get("schema_version") != 8 or marker.get("status") != "complete"
            or not isinstance(delivery, str) or re.fullmatch(r"[0-9a-f]{64}", delivery) is None
            or marker.get("run", {}).get("snapshot_id") != delivery
            or marker.get("pipeline", {}).get("name") != "Alopex"
            or Path(str(marker.get("project_dir", ""))).resolve() != root):
        raise DownstreamQCError("DNA run_manifest.json is not a current complete delivery")
    manifest_path = root / "03_results/QC_Results/sample_manifest.tsv"
    outputs = marker.get("outputs")
    if not isinstance(outputs, list) or not all(isinstance(entry, dict) for entry in outputs):
        raise DownstreamQCError("DNA run_manifest.json has invalid output inventory")
    entries = [x for x in outputs if x.get("path") == "QC_Results/sample_manifest.tsv"]
    if len(entries) != 1:
        raise DownstreamQCError("DNA delivery does not declare one sealed sample_manifest.tsv")
    try:
        manifest_bytes = manifest_path.read_bytes()
    except OSError as error:
        raise DownstreamQCError(f"Cannot read DNA sealed sample_manifest.tsv: {manifest_path}") from error
    manifest_hash = hashlib.sha256(manifest_bytes).hexdigest()
    if entries[0].get("size_bytes") != len(manifest_bytes) or entries[0].get("sha256") != manifest_hash:
        raise DownstreamQCError("DNA sealed sample_manifest.tsv disagrees with the published inventory")
    manifest = pd.read_csv(io.BytesIO(manifest_bytes), sep="\t", dtype=str, keep_default_na=False)
    required = {"sample_id", "project_sample_id", "plate_id", "dna_raw_sample", "rna_sample", "protocol",
                "dna_status", "analysis_status", "backend_qc_status", "cpg_path"}
    missing = sorted(required - set(manifest.columns))
    if missing:
        raise DownstreamQCError(f"DNA sealed sample_manifest.tsv lacks required columns: {missing}")
    if manifest["sample_id"].duplicated().any():
        raise DownstreamQCError("DNA sealed sample_manifest.tsv sample_id is not unique")
    downstream_root = root / "06_downstream"
    qc_path = downstream_root / "DNAme_QC_Information.csv"
    completion = downstream_root / ".cache/DNAme_QC_Input_Identity.json"
    if any(path.is_symlink() for path in (downstream_root, completion.parent, completion, qc_path)):
        raise DownstreamQCError("DNA downstream QC paths must not be symlinks")
    try:
        completion_state = completion.stat()
        completion_bytes = completion.read_bytes()
        identity = json.loads(completion_bytes)
    except (OSError, ValueError) as error:
        raise DownstreamQCError("Current DNA delivery has no completed downstream QC; run DNA downstream") from error
    if (not isinstance(identity, dict) or identity.get('delivery_id') != delivery
            or identity.get('snapshot_id') != delivery
            or identity.get('run_manifest_sha256') != hashlib.sha256(marker_bytes).hexdigest()):
        raise DownstreamQCError("Current DNA delivery disagrees with downstream QC completion; run DNA downstream")
    try:
        qc_bytes = qc_path.read_bytes()
    except OSError as error:
        raise DownstreamQCError(f"Current DNA delivery has no downstream QC table; run DNA downstream: {qc_path}") from error
    fields = next(csv.reader(io.StringIO(qc_bytes.decode("utf-8-sig"))))
    identity_columns = ("Sample_ID", "CloneID", "PlateID")
    na_values = {column: ["", "nan", "NaN", "NA", "N/A", "None", "null"] for column in fields if column not in identity_columns}
    dna_qc = pd.read_csv(io.BytesIO(qc_bytes), dtype={column: str for column in identity_columns},
                         keep_default_na=False, na_values=na_values, float_precision="round_trip")
    dna_qc = _identity_frame(dna_qc, ("Sample_ID", "CloneID"), "DNAme QC table")
    dna_qc["PlateID"] = dna_qc["PlateID"].fillna("").astype(str)
    active = manifest.loc[(manifest["dna_status"] == "Pass") & (manifest["analysis_status"] == "PASS")
                          & (manifest["backend_qc_status"] == "PASS") & manifest["cpg_path"].ne("")]
    if dna_qc["Sample_ID"].duplicated().any() or set(active["sample_id"]) != set(dna_qc["Sample_ID"]):
        raise DownstreamQCError("Current DNA QC cells disagree with the active sealed manifest cells")
    try:
        after = completion.stat()
        state_keys = ('st_dev', 'st_ino', 'st_size', 'st_mtime_ns', 'st_ctime_ns')
        unchanged = all(getattr(after, key) == getattr(completion_state, key) for key in state_keys)
        unchanged = unchanged and completion.read_bytes() == completion_bytes
    except OSError:
        unchanged = False
    if not unchanged or marker_path.read_bytes() != marker_bytes:
        raise DownstreamQCError("DNA delivery or QC completion changed while reading QC inputs; retry integration")
    dna_qc.attrs["dna_provenance"] = {
        "project_dir": str(root), "delivery_id": delivery, "run_manifest_sha256": hashlib.sha256(marker_bytes).hexdigest(),
        "sealed_manifest_sha256": manifest_hash, "qc_path": str(qc_path), "qc_sha256": hashlib.sha256(qc_bytes).hexdigest(),
    }
    return dna_qc, manifest
