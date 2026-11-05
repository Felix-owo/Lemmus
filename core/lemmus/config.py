"""Lemmus：配置、默认值与协议能力。"""
from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any
import yaml


CONFIG_SCHEMA_VERSION = 2


RNA_RESOLUTION = 'parsimony-gene'


RRNA_COMPETITOR_ID = 'RRNA45S_COMPETITOR'


RRNA_REFERENCES = {
    'hg38': ('NR_146144.1', 'd5a3d20f8f84fe5770096d55644d6b7fc8d5116e64343b1259805df90c17ebcd'),
    'mm10': ('NR_046233.2', '1bea30ddff317d9a76f5e3827082dc4446d8f365115d9dedd1d6fdbb274875ee'),
}


class ConfigError(ValueError):
    """Raised when a project config violates the current release contract."""


def load_config(path: str | Path) -> dict[str, Any]:
    config_path = Path(path).expanduser().resolve()
    try:
        payload = yaml.safe_load(config_path.read_text(encoding='utf-8')) or {}
    except (OSError, yaml.YAMLError) as exc:
        raise ConfigError(f'Cannot read project config {config_path}: {exc}') from exc
    if not isinstance(payload, dict):
        raise ConfigError('Project config root must be a mapping')
    validate_config(payload)
    return select_config(payload)


def select_config(payload: dict[str, Any]) -> dict[str, Any]:
    """Return the active protocol's settings without modifying the document.

    Both modules may be displayed in one config. Existing flat schema-2
    documents keep their meaning; mixing flat and modular settings is rejected.
    """
    import copy
    if not isinstance(payload, dict):
        raise ConfigError('Project config root must be a mapping')
    modules = {'cabernet', 'droplet'}
    if not modules.intersection(payload):
        return copy.deepcopy(payload)
    mode = payload.get('pipeline_mode')
    if mode not in modules:
        raise ConfigError('pipeline_mode must be cabernet or droplet')
    if not modules.issubset(payload) or any(not isinstance(payload[m], dict) for m in modules):
        raise ConfigError('Modular config requires both cabernet and droplet mappings')
    settings = {'RNA', 'DARLIN', 'QC', 'barcode_map'}
    if settings.intersection(payload):
        raise ConfigError('Do not mix top-level RNA/DARLIN/QC/barcode_map with protocol modules')
    common = {key: value for key, value in payload.items() if key not in modules}
    if set(payload[mode]) - settings:
        raise ConfigError(f'Unsupported setting in {mode} module: {sorted(set(payload[mode]) - settings)}')
    active = copy.deepcopy({**common, **payload[mode]})
    if not isinstance(active.get('RNA'), dict):
        raise ConfigError(f'{mode}.RNA must be a mapping')
    fixed_rna = {'resolution'} | ({'min_reads', 'expected_orientation'} if mode == 'droplet' else set())
    if fixed_rna.intersection(active['RNA']):
        raise ConfigError(f'{mode}.RNA contains fixed settings: {sorted(fixed_rna.intersection(active["RNA"]))}')
    if mode == 'droplet' and {'DARLIN', 'QC', 'barcode_map'}.intersection(payload[mode]):
        raise ConfigError('droplet only accepts tunable RNA settings; DARLIN/QC/barcode_map are not configurable')
    # 固定科学策略只进入有效配置，不作为用户可调叶子展示。
    active['RNA']['resolution'] = RNA_RESOLUTION
    if mode == 'droplet':
        active['RNA'].update(min_reads=1, expected_orientation='rc')
        active['DARLIN'] = {'enabled': False}
        active['QC'] = {'read_fate': {'enabled': False}}
    return active


def config_document_schema() -> dict[str, Any]:
    """Schema for a single document displaying both protocol modules."""
    import copy
    active = yaml.safe_load(CONFIG_SCHEMA_YAML)
    module_keys = ('RNA', 'DARLIN', 'QC', 'barcode_map')
    common_keys = ('schema_version', 'pipeline_mode', 'raw_data_dir', 'sample_manifest', 'paths', 'slurm')
    properties = {key: copy.deepcopy(active['properties'][key]) for key in common_keys}
    for mode in ('cabernet', 'droplet'):
        keys = module_keys if mode == 'cabernet' else ('RNA',)
        properties[mode] = {
            'type': 'object', 'additionalProperties': False,
            'required': list(keys),
            'properties': {key: copy.deepcopy(active['properties'][key]) for key in keys},
        }
        rna_properties = properties[mode]['properties']['RNA']['properties']
        del rna_properties['resolution']
        if mode == 'droplet':
            del rna_properties['min_reads']
            del rna_properties['expected_orientation']
    return {'type': 'object', 'additionalProperties': False,
            'required': [*common_keys, 'cabernet', 'droplet'], 'properties': properties}


def _positive_int(value: object, label: str) -> int:
    if isinstance(value, bool):
        raise ConfigError(f'{label} must be a positive integer')
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise ConfigError(f'{label} must be a positive integer') from exc
    if parsed < 1:
        raise ConfigError(f'{label} must be a positive integer')
    return parsed


def validate_config(payload: dict[str, Any]) -> None:
    if {'cabernet', 'droplet'}.intersection(payload):
        import jsonschema
        # Reject mixed layouts before resolving either module.
        select_config(payload)
        try:
            jsonschema.Draft202012Validator(config_document_schema()).validate(payload)
        except jsonschema.ValidationError as exc:
            location = '.'.join(str(part) for part in exc.absolute_path) or '<root>'
            raise ConfigError(f'config schema validation failed at {location}: {exc.message}') from exc
        for mode in ('cabernet', 'droplet'):
            try:
                validate_config(select_config({**payload, 'pipeline_mode': mode}))
            except ConfigError as exc:
                raise ConfigError(f'{mode} module: {exc}') from exc
        return
    if payload.get('schema_version') != CONFIG_SCHEMA_VERSION:
        raise ConfigError(f"schema_version must be {CONFIG_SCHEMA_VERSION}; found {payload.get('schema_version')!r}")
    protocol = str(payload.get('pipeline_mode', '')).strip().lower()
    if protocol not in PROTOCOLS:
        raise ConfigError(f'pipeline_mode must be one of {sorted(PROTOCOLS)}')
    rna = payload.get('RNA')
    darlin = payload.get('DARLIN')
    slurm = payload.get('slurm')
    if not isinstance(rna, dict) or not isinstance(darlin, dict) or (not isinstance(slurm, dict)):
        raise ConfigError('RNA, DARLIN and slurm must be mappings')
    _positive_int(rna.get('cores'), 'RNA.cores')
    if rna.get('resolution', RNA_RESOLUTION) != RNA_RESOLUTION:
        raise ConfigError(f'RNA.resolution must be {RNA_RESOLUTION}; multi-gene EM is not the production policy')
    if 'quant_threads' in rna:
        _positive_int(rna.get('quant_threads'), 'RNA.quant_threads')
        if int(rna['quant_threads']) > int(rna['cores']):
            raise ConfigError('RNA.quant_threads must not exceed RNA.cores')
    _positive_int(slurm.get('jobs'), 'slurm.jobs')
    _positive_int(slurm.get('local_cores'), 'slurm.local_cores')
    if not isinstance(darlin.get('enabled'), bool):
        raise ConfigError('DARLIN.enabled must be a boolean')
    if protocol == 'droplet':
        if darlin['enabled'] or darlin.get('include_rna_fastq', False):
            raise ConfigError('DD-MET5 RNA requires DARLIN.enabled=false and no RNA DARLIN rescue')
        for key, expected in (('min_reads', 1), ('expected_orientation', 'rc')):
            if key in rna and rna[key] != expected:
                raise ConfigError(f'DD-MET5 RNA requires RNA.{key}={expected}')
        if payload.get('barcode_map') is not None:
            raise ConfigError('Droplet uses the pinned ME5 design whitelist; remove barcode_map')
        if payload.get('QC', {}).get('read_fate', {}).get('enabled', False):
            raise ConfigError('Droplet read accounting is provided by preprocessing and native quantification QC')
    # JSON Schema 是项目配置唯一完整结构 gate；与上述跨字段语义检查共享同一加载链。
    try:
        import jsonschema
        schema = yaml.safe_load(CONFIG_SCHEMA_YAML)
        jsonschema.Draft202012Validator(schema).validate(payload)
    except jsonschema.ValidationError as exc:
        location = '.'.join(str(part) for part in exc.absolute_path) or '<root>'
        raise ConfigError(f'config schema validation failed at {location}: {exc.message}') from exc


def resolve_project_config(path: str | Path) -> tuple[Path, Path, dict[str, Any]]:
    config_path = Path(path).expanduser().resolve()
    if config_path.name != 'config.yaml' or config_path.parent.name != '00_config':
        raise ConfigError(f'RNA_DARLIN_Pipeline only accepts <project>/00_config/config.yaml; found {config_path}')
    project_dir = config_path.parent.parent.resolve()
    return (config_path, project_dir, load_config(config_path))


MANIFEST_FIELDS = ('sample_id', 'rna_sample', 'darlin_sample', 'species', 'notes')


SAMPLE_ID_PATTERN = re.compile('^[A-Za-z0-9_.-]+$')


@dataclass(frozen=True)
class Protocol:
    name: str
    barcode_mate: str
    insert_mate: str
    expression_enabled: bool
    darlin_enabled: bool


PROTOCOLS = {'cabernet': Protocol('cabernet', 'R2', 'R1', True, True),
             'droplet': Protocol('droplet', 'R1', 'R2', True, False)}


def get_protocol(name: str) -> Protocol:
    normalized = str(name).strip().lower()
    try:
        return PROTOCOLS[normalized]
    except KeyError as exc:
        raise ValueError(f'Unsupported RNA/DARLIN protocol: {name!r}') from exc


DROPLET_SOURCE_COMMIT = "69ba000a7df409c88d813eae2b3457643763ec28"


DROPLET_WHITELIST_SHA256 = "67981bee6acda8db275257a54e727d1e0a0bf4b884c10d8e761d04eaa9b53588"


DROPLET_WHITELIST_ROWS = 829440


DROPLET_WHITELIST_URL = (
    f"https://raw.githubusercontent.com/seekgene/SeekSoulMethyl/{DROPLET_SOURCE_COMMIT}/"
    "dependence/seeksoultools/utils/barcode/ME5/U3CB_methylation.txt.gz"
)


DROPLET_R2_ADAPTERS = (
    ("GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC", "5"),
    ("AAAAAAAAAAAA", "3"),
    ("AAGCAGTGGTATCAACGCAGAGTACATGG", "5"),
    ("CCCATATAAGAAA", "3"),
)


DEFAULT_EXTRACTION = {'min_reads_per_allele_group': 1, 'seq_5prime_max_edits': 1, 'seq_3prime_upper_n': 15, 'seq_3prime_max_mismatches': 1, 'min_cell_read_fraction': 0.8, 'read_cutoff': 3, 'rna_rescue_read_cutoff': 1, 'rna_novel_min_reads': 2, 'rna_novel_min_umis': 2, 'rna_catalog_max_edits': 0, 'rna_catalog_min_partial_length': 40, 'chunk_size': 100000}


DEFAULT_QC = {'min_read_count': 3, 'min_umi_count': 1, 'require_edited_barcode': True, 'require_3prime_detected': True, 'max_3prime_mismatches': 1, 'allow_read_end_rescue': False, 'min_read_end_read_count': 20, 'min_read_end_umi_count': 3, 'min_read_end_cell_read_fraction': 1.0, 'min_read_end_barcode_length': 250, 'max_ambiguous_bases': 0, 'max_barcode_complexity': 6, 'max_reference_sample_count': 4, 'require_reference_evidence': True}


CONFIG_TEMPLATE_YAML = r"""# Lemmus project configuration: two protocol modules in one file.
# Only the module selected by pipeline_mode supplies RNA/DARLIN/QC settings.
# rc is a mapping orientation, not a protocol selector.
schema_version: 2
pipeline_mode: cabernet

# ===== Shared project paths and scheduler =====
raw_data_dir: 01_raw/RNA
sample_manifest: 00_config/sample_manifest.tsv
paths:
  work_root: 02_work
  results_root: 03_results
  logs_root: 04_logs
slurm:
  jobs: 30
  local_cores: 2
  latency_wait: 120
  partition: amd-ep2
  qos: huge

# ===== Cabernet: R1=cDNA; R2=CB8/10+UMI8 =====
# Only active when pipeline_mode: cabernet. RNA-only: set cabernet.DARLIN.enabled=false.
cabernet:
  barcode_map: 00_config/Barcode_Map.csv  # 与 Alopex Cabernet 共享 PlateID/RNA_Barcode/Cell_Order；DNA 外部关联填写 RNA raw basename
  RNA:
    reference: mm10
    cores: 16
    # quant_threads: 13   # 可选：simpleaf quant 线程预算，缺省 = cores。设低于 cores
                           # 可让单线程 DARLIN 提取与 RNA 定量并行（piscem 实测
                           # map_threads=12 时 13 即可），不得大于 cores。
    mem_mb: 16384
    runtime: 720
    # Standard 3' scRNA-seq: R1=cDNA; R2=CB+UMI. CB length is inferred from Barcode_Map.csv; UMI=8 bp.
    min_reads: 10
    expected_orientation: both

  DARLIN:
    enabled: true
    mem_mb: 13312
    runtime: 240
    fastq_dir: 01_raw/DARLIN
    # Recover paired DARLIN transcripts from the corresponding RNA FASTQ as secondary evidence.
    include_rna_fastq: false
    require_all_rna_samples: false
    template:
      - Tigre_2022_v2
      - Rosa_v2
      - cCARLIN
    extraction:
      # Performance-only streaming bound; biological thresholds below retain their defaults.
      chunk_size: 100000
      seq_5prime_max_edits: 1
      rna_rescue_read_cutoff: 1
      rna_novel_min_reads: 2
      rna_novel_min_umis: 2
      rna_catalog_max_edits: 0
      rna_catalog_min_partial_length: 40

  QC:
    double_cell:
      enabled: true
    # Full per-read fate ledger re-reads the complete FASTQ set; keep it out of
    # the default production DAG.
    read_fate:
      enabled: false


# ===== DD-MET5 Droplet: R1=CB17+UMI12; R2=cDNA =====
# Only active when pipeline_mode: droplet. Official ME5 whitelist; no plate map.
droplet:  # 对应 Alopex analysis.protocol: droplet；共享 ME5 设计，DNA/RNA 独立调用细胞，下游用 DNA sealed manifest 的可选 rna_sample 关联
  RNA:
    reference: hg38
    cores: 16
    # quant_threads: 16  # Optional explicit quant budget; default = cores, 1..cores.
    mem_mb: 49152
    runtime: 720
"""


CONFIG_SCHEMA_YAML = r"""$schema: https://json-schema.org/draft/2020-12/schema
type: object
additionalProperties: false
required: [schema_version, pipeline_mode, raw_data_dir, sample_manifest, paths, RNA, DARLIN, QC, slurm]
allOf:
  - if:
      properties:
        pipeline_mode: {const: cabernet}
    then:
      required: [barcode_map]
properties:
  schema_version: {const: 2}
  pipeline_mode: {enum: [cabernet, droplet]}
  raw_data_dir: {type: string, minLength: 1}
  sample_manifest: {type: string, minLength: 1}
  barcode_map: {type: string, minLength: 1}
  paths:
    type: object
    additionalProperties: false
    required: [work_root, results_root, logs_root]
    properties:
      work_root: {type: string, minLength: 1}
      results_root: {type: string, minLength: 1}
      logs_root: {type: string, minLength: 1}
  RNA:
    type: object
    additionalProperties: false
    required: [reference, cores]
    properties:
      reference: {type: string, minLength: 1}
      cores: {type: integer, minimum: 1}
      quant_threads: {type: integer, minimum: 1}
      mem_mb: {type: integer, minimum: 1}
      runtime: {type: integer, minimum: 1}
      min_reads: {type: integer, minimum: 1}
      expected_orientation: {enum: [fw, rc, both]}
      resolution: {const: parsimony-gene}
  DARLIN:
    type: object
    additionalProperties: false
    required: [enabled]
    properties:
      enabled: {type: boolean}
      mem_mb: {type: integer, minimum: 1}
      runtime: {type: integer, minimum: 1}
      fastq_dir: {type: string, minLength: 1}
      include_rna_fastq: {type: boolean}
      require_all_rna_samples: {type: boolean}
      template:
        type: array
        minItems: 1
        uniqueItems: true
        items: {enum: [cCARLIN, Rosa, Rosa_v2, Tigre, Tigre_2022, Tigre_2022_v2]}
      template_locus:
        type: object
        additionalProperties: {type: string, minLength: 1}
      extraction:
        type: object
        additionalProperties: false
        properties:
          min_reads_per_allele_group: {type: integer, minimum: 1}
          seq_5prime_max_edits: {type: integer, minimum: 0, maximum: 2}
          seq_3prime_upper_n: {type: integer, minimum: 1}
          seq_3prime_max_mismatches: {type: integer, minimum: 0}
          min_cell_read_fraction: {type: number, minimum: 0, maximum: 1}
          chunk_size: {type: integer, minimum: 1}
          rna_rescue_read_cutoff: {type: integer, minimum: 1}
          rna_novel_min_reads: {type: integer, minimum: 1}
          rna_novel_min_umis: {type: integer, minimum: 1}
          rna_catalog_max_edits: {type: integer, minimum: 0, maximum: 2}
          rna_catalog_min_partial_length: {type: integer, minimum: 20}
          read_cutoff:
            anyOf:
              - {type: integer, minimum: 1}
              - {type: string, minLength: 1}
      qc:
        type: object
        additionalProperties: false
        properties:
          min_read_count: {type: integer, minimum: 0}
          min_umi_count: {type: integer, minimum: 0}
          require_edited_barcode: {type: boolean}
          require_3prime_detected: {type: boolean}
          max_3prime_mismatches: {type: integer, minimum: 0}
          allow_read_end_rescue: {type: boolean}
          min_read_end_read_count: {type: integer, minimum: 0}
          min_read_end_umi_count: {type: integer, minimum: 0}
          min_read_end_cell_read_fraction: {type: number, minimum: 0, maximum: 1}
          min_read_end_barcode_length: {type: integer, minimum: 0}
          max_ambiguous_bases: {type: integer, minimum: 0}
          max_barcode_complexity: {type: integer, minimum: 1}
          max_reference_sample_count: {type: integer, minimum: 0}
          require_reference_evidence: {type: boolean}
  QC:
    type: object
    additionalProperties: false
    properties:
      double_cell:
        type: object
        additionalProperties: false
        required: [enabled]
        properties:
          enabled: {type: boolean}
      read_fate:
        type: object
        additionalProperties: false
        required: [enabled]
        properties:
          enabled: {type: boolean}
  slurm:
    type: object
    additionalProperties: false
    required: [jobs, local_cores, latency_wait, partition, qos]
    properties:
      jobs: {type: integer, minimum: 1}
      local_cores: {type: integer, minimum: 1}
      latency_wait: {type: integer, minimum: 0}
      partition: {type: string, minLength: 1}
      qos: {type: string, minLength: 1}
"""


SAMPLE_MANIFEST_TEMPLATE = r"""sample_id	rna_sample	darlin_sample	species	notes
"""


def config_scalars_main(argv=None):
    """导出 launcher 所需的 implementation/runtime 标量，不复制 config 默认。"""
    parser = argparse.ArgumentParser(prog="lemmus config-scalars")
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args(argv)
    cfg = load_config(args.config)
    slurm = cfg["slurm"]
    values = {
        "JOBS": slurm["jobs"],
        "PIPELINE_MODE": cfg['pipeline_mode'],
        "RNA_MEM_MB": cfg['RNA'].get('mem_mb', 16384),
        "RNA_CORES": cfg["RNA"]["cores"],
        "LOCAL_CORES": slurm["local_cores"],
        "LATENCY": slurm["latency_wait"],
        "PARTITION": slurm["partition"],
        "QOS": slurm["qos"],
    }
    import shlex
    for key, value in values.items():
        print(f"{key}={shlex.quote(str(value))}")
    return 0
