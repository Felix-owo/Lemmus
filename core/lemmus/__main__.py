"""Lemmus：内部命令分派。"""
from __future__ import annotations

import importlib
import sys


def main(argv=None):
    """统一内部 CLI；用户入口仍只有 doctor.sh 与 run_pipeline.sh。"""
    import sys
    args = list(sys.argv[1:] if argv is None else argv)
    if not args:
        raise SystemExit("usage: python -m lemmus <command> ...")
    dispatch = {
        "cross-pipeline-check": ("maintenance", "cross_pipeline_check_main"),
        "init-project": ("inputs", "init_project_main"),
        "prepare-droplet-whitelist": ("maintenance", "prepare_droplet_whitelist_main"),
        "prepare-calling": ("maintenance", "prepare_calling_main"),
        "preprocess-droplet-rna": ("rna", "preprocess_droplet_main"),
        "call-droplet-rna": ("rna", "call_droplet_rna_main"),
        "manifest-refresh": ("inputs", "manifest_refresh_main"),
        "read-fate": ("read_fate", "read_fate_main"),
        "extract-barcodes-multi": ("darlin", "extract_barcodes_multi_main"),
        "merge-barcodes": ("darlin", "barcode_table_main"),
        "build-clone-table": ("darlin", "barcode_table_main"),
        "assert-fastq": ("inputs", "assert_fastq_main"),
        "prefilter-rna-darlin": ("inputs", "prefilter_rna_darlin_main"),
        "write-rna-whitelist": ("inputs", "write_rna_whitelist_main"),
        "finalize-rna-h5ad": ("rna", "finalize_rna_h5ad_main"),
        "rna-reference": ("maintenance", "rna_reference_main"),
        "prepare-rna-counts": ("rna", "prepare_rna_counts_main"),
        "multiqc-custom": ("reports", "multiqc_custom_main"),
        "env-spec": ("maintenance", "env_spec_main"),
        "parser-fixture": ("maintenance", "parser_fixture_main"),
        "self-check": ("maintenance", "self_check_main"),
        "config-scalars": ("config", "config_scalars_main"),
    }
    if args[0] in ("-h", "--help"):
        import argparse
        parser = argparse.ArgumentParser(prog="lemmus", description="Lemmus 内部阶段命令；项目操作使用 Shell 入口。")
        parser.add_argument("command", choices=dispatch)
        parser.print_help()
        return 0
    command = args.pop(0)
    try:
        module_name, function_name = dispatch[command]
    except KeyError as exc:
        raise SystemExit(f"unknown command: {command}") from exc
    handler = getattr(importlib.import_module(f".{module_name}", __package__), function_name)
    if command in ("merge-barcodes", "build-clone-table"):
        args.insert(0, command)
    result = handler(args)
    return 0 if result is None else result


if __name__ == "__main__":
    raise SystemExit(main())
