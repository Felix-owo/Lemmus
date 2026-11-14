"""Lemmus：持久科学指标与 MultiQC 展示。"""
from __future__ import annotations
from .config import DARLIN_SOURCE_CONTRACT
from .io import atomic_output, read_delimited_frame

import argparse
import csv
import glob
import json
import math
import re
from collections.abc import Sequence
from pathlib import Path
import numpy as np
import pandas as pd


def _flatten_json_scalars(payload: object, prefix: str = "") -> dict[str, object]:
    """Flatten scalar JSON values while preserving their original key path."""
    flat: dict[str, object] = {}
    if isinstance(payload, dict):
        for key, value in payload.items():
            child = f"{prefix}.{key}" if prefix else str(key)
            flat.update(_flatten_json_scalars(value, child))
    elif isinstance(payload, list):
        # Lists are intentionally skipped: MultiQC summary tables should contain
        # unambiguous scalar metrics, not positional arrays with unclear semantics.
        return flat
    elif isinstance(payload, (str, int, float, bool)) or payload is None:
        flat[prefix] = payload
    return flat


def _read_json_mapping(path: Path) -> dict[str, object]:
    import json
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return payload


def _report_count(value: object, label: str) -> int:
    """报告保留producer的非负整数计数；不能用int把小数或非法值改成观测。"""
    if not re.fullmatch(r'[0-9]+', str(value)):
        raise ValueError(f'{label}: expected a nonnegative integer count, got {value!r}')
    return int(value)


def _report_fraction(value: object, label: str) -> float:
    """barcode纠错比例的分母在producer定义；报告只接受有限的0–1原始比例。"""
    fraction = float(value)
    if not math.isfinite(fraction) or not 0 <= fraction <= 1:
        raise ValueError(f'{label}: expected a finite fraction in [0, 1]')
    return fraction


def build_multiqc_custom_content(
    *,
    samples: list[str],
    rna_h5ads: list[Path],
    map_infos: list[Path],
    quant_jsons: list[Path],
    output_dir: Path,
    darlin_clone_table: Path | None = None,
    darlin_barcode_metrics: list[Path] | None = None,
    alignment_output: Path | None = None,
    sample_manifest: Path | None = None,
    barcode_map: Path | None = None,
) -> list[Path]:
    """Create deterministic MultiQC custom-content JSON files.

    FastQC remains a native MultiQC module. RNA expression summaries are derived
    from finalized project RNA.h5ad files. simpleaf/piscem/alevin-fry JSON is
    represented using exact scalar key paths so no third-party metric is silently
    reinterpreted. DARLIN summaries come from the final clone table when enabled.
    """
    import json
    import numpy as np
    import anndata as ad
    import hdf5plugin  # noqa: F401 - register compression filters before H5AD reads
    import pandas as pd

    if not (len(samples) == len(rna_h5ads) == len(map_infos) == len(quant_jsons)):
        raise ValueError("MultiQC RNA inputs must align one-to-one with sample IDs")
    if len(set(samples)) != len(samples):
        raise ValueError("MultiQC sample IDs must be unique")
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []

    expression_data: dict[str, dict[str, object]] = {}
    simpleaf_data: dict[str, dict[str, object]] = {}
    rna_cell_ids: dict[str, set[str]] = {}
    for sample, h5ad_path, map_path, quant_path in zip(samples, rna_h5ads, map_infos, quant_jsons, strict=True):
        adata = ad.read_h5ad(h5ad_path, backed="r")
        try:
            if not adata.obs_names.is_unique:
                raise ValueError(f'Finalized RNA cell identities are not unique: {h5ad_path}')
            n_cells = int(adata.n_obs)
            n_genes = int(adata.n_vars)
            rna_cell_ids[sample] = set(str(cell_id) for cell_id in adata.obs_names)
            # backed sparse matrices are not uniformly summable without materialising;
            # per-cell detected-gene metrics are therefore only calculated when present
            # in obs. Core dimensions are always reliable and cheap.
            row: dict[str, object] = {"Cells": n_cells, "Genes": n_genes}
            metric_columns = (
                ("total_expected_counts", "Median expected UMIs/cell", 1.0),
                ("total_supported_expected_counts", "Median supported expected UMIs/cell", 1.0),
                ("subunit_expected_counts", "Median subunit expected UMIs/cell", 1.0),
                ("n_genes_supported", "Median supported genes/cell (>=1 expected UMI)", 1.0),
                ("n_genes_positive_expected", "Median positive genes/cell (>0 expected UMI)", 1.0),
                ("n_genes_fractional_only", "Median fractional-only genes/cell", 1.0),
                ("fractional_nnz_fraction", "Median fractional nonzero (%)", 100.0),
                ("subunit_count_mass_fraction", "Median subunit count mass (%)", 100.0),
            )
            for candidate, label, multiplier in metric_columns:
                if candidate in adata.obs.columns and n_cells:
                    values = pd.to_numeric(adata.obs[candidate], errors="coerce").dropna()
                    if len(values):
                        row[label] = float(np.median(values.to_numpy())) * multiplier
            provenance = adata.uns.get('rna_darlin_pipeline', {})
            quantification = adata.uns.get('rna_quantification', {})
            for key, value in _flatten_json_scalars(quantification, 'rna_quantification').items():
                if isinstance(value, (int, float, bool, str, np.number)):
                    row[key] = value.item() if isinstance(value, np.generic) else value
            if provenance.get('protocol') == 'droplet':
                for section in ('calling', 'preprocessing', 'barcode_correction'):
                    for key, value in _flatten_json_scalars(provenance.get(section, {}), section).items():
                        if isinstance(value, (int, float, bool, np.number)):
                            row[key] = value.item() if isinstance(value, np.generic) else value
            expression_data[sample] = row
        finally:
            if getattr(adata, "file", None) is not None:
                adata.file.close()

        map_flat = _flatten_json_scalars(_read_json_mapping(map_path), "map")
        quant_flat = _flatten_json_scalars(_read_json_mapping(quant_path), "quant")
        # Preserve exact tool key paths. Numeric/bool scalars are most useful for
        # comparison tables; short strings are kept for modes / versions / status.
        merged: dict[str, object] = {}
        for key, value in {**map_flat, **quant_flat}.items():
            if isinstance(value, (int, float, bool)) and not isinstance(value, complex):
                merged[key] = value
            elif isinstance(value, str) and len(value) <= 80:
                merged[key] = value
        simpleaf_data[sample] = merged

    if samples:
        expression_payload = {
            "id": "rna_expression_summary",
            "section_name": "RNA expression summary",
            "description": "Finalized RNA.h5ad support-aware QC. For EM expected counts, >=1 expected UMI is the primary detected-gene definition; raw >0 values are retained for audit.",
            "plot_type": "table",
            "pconfig": {"id": "rna_expression_summary_table", "title": "RNA expression summary", "col1_header": "Sample"},
            "data": expression_data,
        }
        expression_path = output_dir / "rna_expression_summary_mqc.json"
        expression_path.write_text(json.dumps(expression_payload, indent=2, sort_keys=True), encoding="utf-8")
        outputs.append(expression_path)

        simpleaf_payload = {
            "id": "simpleaf_raw_metrics",
            "section_name": "simpleaf / piscem / alevin-fry metrics",
            "description": "Scalar values copied from simpleaf map_info.json and quant.json. Column names retain their original JSON key paths to avoid semantic reinterpretation.",
            "plot_type": "table",
            "pconfig": {"id": "simpleaf_raw_metrics_table", "title": "simpleaf / piscem / alevin-fry raw metrics", "col1_header": "Sample"},
            "data": simpleaf_data,
        }
        simpleaf_path = output_dir / "simpleaf_raw_metrics_mqc.json"
        simpleaf_path.write_text(json.dumps(simpleaf_payload, indent=2, sort_keys=True), encoding="utf-8")
        outputs.append(simpleaf_path)

    if darlin_clone_table is not None:
        clone_df = read_delimited_frame(darlin_clone_table)
        expected = {"Sample_ID", "DARLIN_QC_Pass", "DARLIN_Total_Read_Count", "DARLIN_Total_UMI_Count", "DARLIN_Barcodes_Used", "DARLIN_Double_Cell_Suspect"}
        missing = expected.difference(clone_df.columns)
        if missing:
            raise ValueError(f"DARLIN clone table lacks MultiQC columns: {sorted(missing)}")
        if clone_df['Sample_ID'].eq('').any() or clone_df['Sample_ID'].duplicated().any():
            raise ValueError('DARLIN clone table requires nonempty unique Sample_ID')
        def _table_bool(value: object) -> bool:
            if isinstance(value, bool):
                return value
            text = str(value).strip().lower()
            if text in {"true", "1", "yes"}:
                return True
            if text in {"false", "0", "no"}:
                return False
            raise ValueError(f"Unexpected boolean value in DARLIN clone table: {value!r}")

        darlin_data: dict[str, dict[str, object]] = {}
        for row in clone_df.itertuples(index=False):
            record = row._asdict()
            sid = str(record["Sample_ID"])
            darlin_data[sid] = {
                "QC pass": _table_bool(record["DARLIN_QC_Pass"]),
                "Total reads": _report_count(record["DARLIN_Total_Read_Count"], f'{sid} total reads'),
                "Total UMIs": _report_count(record["DARLIN_Total_UMI_Count"], f'{sid} total UMIs'),
                "Barcodes used": _report_count(record["DARLIN_Barcodes_Used"], f'{sid} barcodes used'),
                "Double-cell suspect": _table_bool(record["DARLIN_Double_Cell_Suspect"]),
            }
        darlin_payload = {
            "id": "darlin_clone_qc",
            "section_name": "DARLIN clone QC",
            "description": "Per-cell/project identity DARLIN QC metrics from the final DARLIN_Clone_ID.csv table.",
            "plot_type": "table",
            "pconfig": {"id": "darlin_clone_qc_table", "title": "DARLIN clone QC", "col1_header": "Sample / cell ID"},
            "data": darlin_data,
        }
        darlin_path = output_dir / "darlin_clone_qc_mqc.json"
        darlin_path.write_text(json.dumps(darlin_payload, indent=2, sort_keys=True), encoding="utf-8")
        outputs.append(darlin_path)

        if alignment_output is not None:
            from .inputs import read_manifest, load_barcode_config, validate_barcode_map_provenance

            if sample_manifest is None or barcode_map is None:
                raise ValueError('RNA/DARLIN alignment requires the project manifest and Barcode_Map')
            validate_barcode_map_provenance(clone_df,barcode_map)
            _, libraries = read_manifest(sample_manifest)
            rna_libraries = {row['sample_id'] for row in libraries if row['rna_sample']}
            if set(samples) != rna_libraries:
                raise ValueError('RNA report libraries disagree with the project manifest')
            cells = load_barcode_config(barcode_map).rows
            owners = {}
            for library in libraries:
                for cell in cells:
                    cell_id = f"{library['sample_id']}_{cell.plate_id}"
                    if cell_id in owners:
                        raise ValueError(f'Ambiguous project cell identity in manifest/Barcode_Map: {cell_id!r}')
                    owners[cell_id] = library['sample_id']
            for sample, identities in rna_cell_ids.items():
                if any(owners.get(cell_id) != sample for cell_id in identities):
                    raise ValueError(f'RNA cell identities disagree with manifest/Barcode_Map: {sample!r}')
            darlin_ids_by_sample: dict[str, set[str]] = {row['sample_id']: set() for row in libraries}
            for sid in clone_df["Sample_ID"].astype(str):
                owner = owners.get(sid)
                if owner is None:
                    raise ValueError(
                        f"DARLIN Sample_ID {sid!r} is absent from the manifest/Barcode_Map cell identities"
                    )
                darlin_ids_by_sample[owner].add(sid)
            alignment_data: dict[str, dict[str, object]] = {}
            tsv_rows = []
            for sample in darlin_ids_by_sample:
                rna_ids = rna_cell_ids.get(sample, set())
                darlin_ids = darlin_ids_by_sample.get(sample, set())
                matched = rna_ids & darlin_ids
                rna_only = sorted(rna_ids - darlin_ids)
                darlin_only = sorted(darlin_ids - rna_ids)
                alignment_data[sample] = {
                    "RNA cells": len(rna_ids),
                    "DARLIN cells": len(darlin_ids),
                    "Matched cells": len(matched),
                    "RNA-only cells": len(rna_only),
                    "DARLIN-only cells": len(darlin_only),
                }
                tsv_rows.append({
                    "sample_id": sample,
                    "rna_cells": len(rna_ids),
                    "darlin_cells": len(darlin_ids),
                    "matched_cells": len(matched),
                    "rna_only_cells": len(rna_only),
                    "darlin_only_cells": len(darlin_only),
                    "rna_only_examples": ";".join(rna_only[:5]),
                    "darlin_only_examples": ";".join(darlin_only[:5]),
                })
            alignment_payload = {
                "id": "rna_darlin_cell_alignment",
                "section_name": "RNA/DARLIN cell alignment",
                "description": "Per-sample cell identity parity between RNA.h5ad obs_names and DARLIN_Clone_ID.csv Sample_ID (both <sample_id>_<PlateID>). Differences are expected when RNA min_reads drops a cell that DARLIN still calls, or vice versa; examples live in rna_darlin_cell_alignment.tsv.",
                "plot_type": "table",
                "pconfig": {"id": "rna_darlin_cell_alignment_table", "title": "RNA/DARLIN cell alignment", "col1_header": "Sample"},
                "data": alignment_data,
            }
            alignment_mqc = output_dir / "rna_darlin_cell_alignment_mqc.json"
            alignment_mqc.write_text(json.dumps(alignment_payload, indent=2, sort_keys=True), encoding="utf-8")
            outputs.append(alignment_mqc)
            alignment_output = Path(alignment_output)
            alignment_output.parent.mkdir(parents=True, exist_ok=True)
            with atomic_output(alignment_output) as alignment_tmp:
                with alignment_tmp.open("w", encoding="utf-8", newline="") as handle:
                    writer = csv.DictWriter(
                        handle,
                        fieldnames=["sample_id", "rna_cells", "darlin_cells", "matched_cells", "rna_only_cells", "darlin_only_cells", "rna_only_examples", "darlin_only_examples"],
                        delimiter="\t",
                        lineterminator="\n",
                    )
                    writer.writeheader()
                    writer.writerows(tsv_rows)

    darlin_barcode_metrics = list(darlin_barcode_metrics or [])
    if darlin_barcode_metrics:
        barcode_data: dict[str, dict[str, object]] = {}
        expected_columns = {
            "sample", "total_reads", "exact_barcode_reads", "corrected_hamming1_reads",
            "ambiguous_barcode_reads", "unmatched_barcode_reads", "accepted_barcode_reads",
            "accepted_fraction", "correction_fraction_of_accepted",
            "source_type", "source_identity_contract",
            "darlin_candidate_reads_any_template", "dual_anchor_reads_any_template",
        }
        for metrics_path in darlin_barcode_metrics:
            metrics = read_delimited_frame(metrics_path, sep="\t")
            missing = expected_columns.difference(metrics.columns)
            if missing or len(metrics) != 1:
                raise ValueError(f"Invalid DARLIN barcode metrics {metrics_path}: missing={sorted(missing)}, rows={len(metrics)}")
            row = metrics.iloc[0]
            if (row['source_type'] not in ('targeted','rna_rescue')
                    or row['source_identity_contract'] != DARLIN_SOURCE_CONTRACT):
                raise ValueError(f'Invalid DARLIN source metadata; rerun extract-barcodes-multi: {metrics_path}')
            sample = str(row["sample"])
            if not sample or sample in barcode_data:
                raise ValueError(f'DARLIN barcode metrics require a nonempty unique sample: {sample!r}')
            counts = {key: _report_count(row[key], f'{metrics_path}: {key}') for key in (
                'total_reads', 'exact_barcode_reads', 'corrected_hamming1_reads',
                'ambiguous_barcode_reads', 'unmatched_barcode_reads', 'accepted_barcode_reads')}
            if (counts['total_reads'] != sum(counts[key] for key in (
                    'exact_barcode_reads', 'corrected_hamming1_reads', 'ambiguous_barcode_reads', 'unmatched_barcode_reads'))
                    or counts['accepted_barcode_reads'] != counts['exact_barcode_reads'] + counts['corrected_hamming1_reads']):
                raise ValueError(f'DARLIN barcode read accounting failed: {metrics_path}')
            candidates = _report_count(row['darlin_candidate_reads_any_template'], f'{sample} candidate reads')
            dual_anchor = _report_count(row['dual_anchor_reads_any_template'], f'{sample} dual-anchor reads')
            if not dual_anchor <= candidates <= counts['accepted_barcode_reads']:
                raise ValueError(f'DARLIN candidate read accounting failed: {metrics_path}')
            barcode_data[sample] = {
                "Total reads": counts["total_reads"],
                "Exact barcode reads": counts["exact_barcode_reads"],
                "Hamming-1 corrected reads": counts["corrected_hamming1_reads"],
                "Ambiguous barcode reads": counts["ambiguous_barcode_reads"],
                "Unmatched barcode reads": counts["unmatched_barcode_reads"],
                "Accepted barcode reads": counts["accepted_barcode_reads"],
                "Accepted barcode (%)": _report_fraction(row["accepted_fraction"], f'{sample} accepted fraction') * 100.0,
                "Correction among accepted (%)": _report_fraction(row["correction_fraction_of_accepted"], f'{sample} correction fraction') * 100.0,
                "Source type": row['source_type'],
                "DARLIN candidate reads": candidates,
                "Dual-anchor reads": dual_anchor,
            }
        barcode_payload = {
            "id": "darlin_barcode_qc",
            "section_name": "DARLIN barcode correction QC",
            "description": "DARLIN cell-barcode recovery using exact matching followed by fail-closed unique Hamming-1 correction, aligned with the RNA permit-list policy.",
            "plot_type": "table",
            "pconfig": {"id": "darlin_barcode_qc_table", "title": "DARLIN barcode correction QC", "col1_header": "Sample"},
            "data": barcode_data,
        }
        barcode_path = output_dir / "darlin_barcode_qc_mqc.json"
        barcode_path.write_text(json.dumps(barcode_payload, indent=2, sort_keys=True), encoding="utf-8")
        outputs.append(barcode_path)
    return outputs


def multiqc_custom_main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Generate MultiQC custom-content files for RNA/DARLIN outputs.")
    parser.add_argument("--sample", action="append", default=[])
    parser.add_argument("--rna-h5ad", action="append", type=Path, default=[])
    parser.add_argument("--map-info", action="append", type=Path, default=[])
    parser.add_argument("--quant-json", action="append", type=Path, default=[])
    parser.add_argument("--darlin-clone-table", type=Path, default=None)
    parser.add_argument("--darlin-barcode-metrics", action="append", type=Path, default=[])
    parser.add_argument("--alignment-output", type=Path, default=None, help="Per-sample RNA/DARLIN cell parity audit TSV.")
    parser.add_argument("--sample-manifest", type=Path, default=None)
    parser.add_argument("--barcode-map", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    build_multiqc_custom_content(
        samples=args.sample,
        rna_h5ads=args.rna_h5ad,
        map_infos=args.map_info,
        quant_jsons=args.quant_json,
        output_dir=args.output_dir,
        darlin_clone_table=args.darlin_clone_table,
        darlin_barcode_metrics=args.darlin_barcode_metrics,
        alignment_output=args.alignment_output,
        sample_manifest=args.sample_manifest,
        barcode_map=args.barcode_map,
    )
    return 0


def multiqc_report_main(argv: Sequence[str] | None = None) -> int:
    """独立报告进程只读原生默认值与Pipeline配置，禁止隐式配置改写样本或删模块。"""
    parser = argparse.ArgumentParser(description="Generate the isolated RNA/DARLIN MultiQC report.")
    parser.add_argument("--source", action="append", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    import multiqc
    from multiqc.core.update_config import ClConfig

    find_user_files = multiqc.config.find_user_files
    multiqc.config.find_user_files = lambda: None
    try:
        result = multiqc.run(
            *(glob.escape(str(path.resolve(strict=True))) for path in args.source),
            cfg=ClConfig(output_dir=args.output_dir, filename="multiqc_report.html", force=True,
                         cl_config=[json.dumps({'data_dir_name': 'multiqc_data', 'fn_clean_trim': [],
                                               'fn_clean_exts': [{'type': 'regex', 'module': 'fastqc',
                                                                  'pattern': r'\.fastq\.gz$'}]})]),
            interactive=False,
        )
        return result.sys_exit_code
    finally:
        multiqc.config.find_user_files = find_user_files
