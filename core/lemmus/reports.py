"""Lemmus：持久科学指标与 MultiQC 展示。"""
from __future__ import annotations

import argparse
import csv
import json
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
        clone_df = pd.read_csv(darlin_clone_table)
        expected = {"Sample_ID", "DARLIN_QC_Pass", "DARLIN_Total_Read_Count", "DARLIN_Total_UMI_Count", "DARLIN_Barcodes_Used", "DARLIN_Double_Cell_Suspect"}
        missing = expected.difference(clone_df.columns)
        if missing:
            raise ValueError(f"DARLIN clone table lacks MultiQC columns: {sorted(missing)}")
        def _table_bool(value: object) -> bool:
            if isinstance(value, bool):
                return value
            text = str(value).strip().lower()
            if text in {"true", "1", "yes"}:
                return True
            if text in {"false", "0", "no", "", "nan"}:
                return False
            raise ValueError(f"Unexpected boolean value in DARLIN clone table: {value!r}")

        darlin_data: dict[str, dict[str, object]] = {}
        for row in clone_df.itertuples(index=False):
            record = row._asdict()
            sid = str(record["Sample_ID"])
            darlin_data[sid] = {
                "QC pass": _table_bool(record["DARLIN_QC_Pass"]),
                "Total reads": int(record["DARLIN_Total_Read_Count"]),
                "Total UMIs": int(record["DARLIN_Total_UMI_Count"]),
                "Barcodes used": int(record["DARLIN_Barcodes_Used"]),
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
            darlin_ids_by_sample: dict[str, set[str]] = {sample: set() for sample in samples}
            for sid in clone_df["Sample_ID"].astype(str):
                owners = [sample for sample in samples if sid.startswith(f"{sample}_")]
                if len(owners) != 1:
                    raise ValueError(
                        f"DARLIN Sample_ID {sid!r} does not map to exactly one RNA sample_id; "
                        "refusing to guess the per-sample cell alignment"
                    )
                darlin_ids_by_sample[owners[0]].add(sid)
            alignment_data: dict[str, dict[str, object]] = {}
            tsv_rows = []
            for sample in samples:
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
            alignment_tmp = alignment_output.with_name(alignment_output.name + ".tmp")
            with alignment_tmp.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["sample_id", "rna_cells", "darlin_cells", "matched_cells", "rna_only_cells", "darlin_only_cells", "rna_only_examples", "darlin_only_examples"],
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(tsv_rows)
            alignment_tmp.replace(alignment_output)

    darlin_barcode_metrics = list(darlin_barcode_metrics or [])
    if darlin_barcode_metrics:
        barcode_data: dict[str, dict[str, object]] = {}
        expected_columns = {
            "sample", "total_reads", "exact_barcode_reads", "corrected_hamming1_reads",
            "ambiguous_barcode_reads", "unmatched_barcode_reads", "accepted_barcode_reads",
            "accepted_fraction", "correction_fraction_of_accepted",
        }
        for metrics_path in darlin_barcode_metrics:
            metrics = pd.read_csv(metrics_path, sep="\t")
            missing = expected_columns.difference(metrics.columns)
            if missing or len(metrics) != 1:
                raise ValueError(f"Invalid DARLIN barcode metrics {metrics_path}: missing={sorted(missing)}, rows={len(metrics)}")
            row = metrics.iloc[0]
            sample = str(row["sample"])
            barcode_data[sample] = {
                "Total reads": int(row["total_reads"]),
                "Exact barcode reads": int(row["exact_barcode_reads"]),
                "Hamming-1 corrected reads": int(row["corrected_hamming1_reads"]),
                "Ambiguous barcode reads": int(row["ambiguous_barcode_reads"]),
                "Unmatched barcode reads": int(row["unmatched_barcode_reads"]),
                "Accepted barcode reads": int(row["accepted_barcode_reads"]),
                "Accepted barcode (%)": float(row["accepted_fraction"]) * 100.0,
                "Correction among accepted (%)": float(row["correction_fraction_of_accepted"]) * 100.0,
                "Source type": str(row["source_type"]) if "source_type" in metrics.columns else ("RNA rescue" if sample.endswith("__RNA") else "targeted"),
                "DARLIN candidate reads": int(row["darlin_candidate_reads_any_template"]) if "darlin_candidate_reads_any_template" in metrics.columns else 0,
                "Dual-anchor reads": int(row["dual_anchor_reads_any_template"]) if "dual_anchor_reads_any_template" in metrics.columns else 0,
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
    )
    return 0
