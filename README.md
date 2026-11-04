# Lemmus

当前版本为Git tag `v5.10`，正式发行包以包内README为准。源码计算位于`core/lemmus/`，下游实现归其`downstream/`子包，顶层`downstream/`只承载Notebook。

Cabernet和Droplet都固定使用 **PG（parsimony-gene）＋物种对应45S竞争参考**。PG拒绝残余多gene歧义，不运行EM，也不将原始分子数取整。支持hg38（人NR_146144.1）与mm10（小鼠NR_046233.2），固定序列及来源见 [45S参考](resources/rrna/README.md)。技术feature `RRNA45S_COMPETITOR` 的计数单独保存在obs/uns；`biological.h5ad`、细胞调用及下游表达不含该feature。它不是生物学单拷贝gene；零技术计数不能证明没有rRNA污染。

本版本采用上述策略；已有其它resolution配置需在Pipeline外明确调整，旧index须用Doctor重建，历史结果不自动重算或覆盖。

Lemmus v5.10（原 RNA_DARLIN_Pipeline）是从 Camellia 拆出的独立 RNA + targeted DARLIN workflow，基于 Snakemake。它可以完成：

- RNA：Cabernet 3′ 与 DD-MET5 Droplet RNA 独立调用和定量（simpleaf 单一 backend）；
- DARLIN：targeted lineage barcode 的 streaming 提取、UMI 级等位基因共识与 directional clone denoising；
- FastQC / MultiQC 项目级 QC 汇总与 RNA/DARLIN cell alignment 审计；
- 下游 QC 可视化、Census 细胞类型注释与跨 Pipeline 整合 QC 表。

普通用户只需要使用两个入口：

```text
core/doctor.sh        安装、构建 reference 并验证 Pipeline
core/run_pipeline.sh  创建项目、检查输入并运行分析
```

Cabernet RNA 路径按标准 TSO-based 3′ UMI scRNA-seq 处理，不再保留 STAR / featureCounts / RNA UMI-tools legacy backend；Git history 负责版本回退。

公开仓库发布可安装的正式版本；开发规范、完整回归和实验记录由私有开发仓库维护。运行代码、Doctor、固定资源与下游 Notebook 随本仓库发布。

本项目沿用 Alopex 的 GNU GPL v3 许可，见 [LICENSE](LICENSE)。

## 全流程

```text
从 GitHub 下载 Pipeline
  → 准备物种 FASTA + GTF
  → 运行 Doctor（conda 构建 / reference 构建 / 回归测试）
  → 显式选择协议初始化项目（--init-project --protocol cabernet 或 droplet）
  → 检查 config.yaml；Cabernet 检查 Barcode_Map.csv
  → 将 RNA / DARLIN paired FASTQ 放入 01_raw/
  → 刷新并检查 sample_manifest.tsv
  → dry-run
  → 正式运行
  → 检查 03_results/ 与 QC 输出
  → 下游 Notebook：RNA QC 可视化 + DNA/RNA 整合 QC 表
```

## 1. 从 GitHub 下载 Pipeline

从 [v5.10 Release](https://github.com/Felix-owo/Lemmus/releases/tag/v5.10) 下载 `Lemmus-v5.10.tar.gz` 和 `SHA256SUMS`，校验后解压：

```bash
shasum -a 256 -c SHA256SUMS
tar -xzf Lemmus-v5.10.tar.gz
cd Lemmus-v5.10
```

也可使用正式 tag 获取代码：

```bash
git clone --branch v5.10 --single-branch https://github.com/Felix-owo/Lemmus.git
cd Lemmus
```

开始前需要：

- macOS 或 Linux；Linux/HPC 正式运行需要 Slurm；
- 一个可以从终端调用的 Conda，例如 Miniforge；
- `git` 与 `curl`；
- HPC 上可以使用 `sbatch` 和 `squeue`。

不要跨机器复制已经构建的 `conda/` 或 reference index；在目标机器上让 Doctor 构建。

## 2. 准备标准 reference

RNA reference 由 genome FASTA + GTF 构建。Pipeline 不会自动下载原始 reference；首次构建前准备好目标物种的 FASTA 与 GTF，然后交给 Doctor：

```bash
core/doctor.sh reference mm10 --fasta /path/mm10.fa --gtf /path/mm10.gtf
```

构建使用 `simpleaf index`，产出结构固定为：

```text
resources/<species>_reference/pg45s/
├── genome.fa
├── genes/genes.gtf
├── rrna_competitor.fa
└── simpleaf/
    ├── ref/t2g.tsv
    └── index/piscem_idx*
```

reference readiness 检查 FASTA、GTF、固定45S序列、两列t2g、piscem index和竞争设计标记；不检查、也不构建 STAR index。`simpleaf/ref/t2g.tsv` 固定为两列 transcript→gene map，输出普通 gene×cell 矩阵；不默认使用三列 USA t2g（USA 输出语义不同且 resolution 约束不同）。

项目 config 只通过 `RNA.reference` 选择 species；生产 reference 固定解析为 Pipeline 根目录的 `resources/<species>_reference/pg45s/`，项目配置不能覆盖 `reference_dir`、`simpleaf_reference`、`genome_fasta` 或 `gtf` 路径，`--init-project` 也不会创建项目级 `resources/`。

## 3. 构建并验证 Pipeline

HPC 标准顺序与 Alopex 一致：

1. 在本机校验正确的 reference、索引、固定45S及封存标记，再上传至 HPC 对应路径并核对版本/校验值。只补传缺失或不匹配的组件，保留 PG＋物种45S 的科学约定。
2. 在登录节点提前下载和准备 Conda 等依赖。Lemmus 需要安装环境时使用 `bash core/doctor.sh conda`；已有正确环境直接复用，不因源码更新重复安装。
3. 在**现有 JupyterLab 的 Terminal** 中进入 Pipeline 根目录，运行 `bash core/doctor.sh`。完整 Doctor 与科学调用在 Jupyter 的原 allocation 中执行，不在登录节点运行，不另起独立 Doctor allocation。

普通部署只要求更新 Pipeline 后跑通标准 Doctor，不默认追加全仓开发回归或额外 benchmark。

分段维护命令如下；环境准备在登录节点，验证在 Jupyter 中执行：

```bash
core/doctor.sh conda          # 构建/更新受管 Conda 环境
core/doctor.sh test           # 编译、配置/schema、运行自检与可用的full-DAG dry-run
```

rules 环境至少提供：Python、simpleaf、piscem、alevin-fry、umi_tools（仅 DARLIN 使用）、FastQC、MultiQC、Java。当前生产版本策略以可从 Bioconda 一致安装为优先（simpleaf 0.28 系列、piscem 0.22.1 系列、alevin-fry 0.18 系列）；版本升级必须同时验证 index-format compatibility 和 doctor smoke。

Doctor 下载 Census model 时使用 atlas 环境自己的 AWS CLI 与 Python，不依赖登录 shell 的 Python。

`core/doctor.sh test` 不依赖Git或开发测试目录：执行Python compile、Shell语法、配置/schema、barcode纠正与整数H5AD写读自检；受管环境和reference存在时还执行Snakemake full-DAG dry-run。自检使用临时合成输入，不能替代原始数据的科学验收。完整开发回归由私有开发仓库维护。

`doctor.sh test` 的末步只读核对同级 `../Alopex`（可用 `DNA_PIPELINE_ROOT` 指定实际安装根目录）：Cabernet 比较两个 `resources/Barcode_Map.csv` 的 `(PlateID, RNA_Barcode, Cell_Order)` 字符串投影；使用 CSV 解析，表头去 BOM/外围空白并忽略大小写，拒绝空表、重复身份和损坏行。Droplet 比较 Lemmus `resources/droplet/U3CB_methylation.txt.gz` 与 Alopex `resources/droplet_ME5_U3CB_methylation.txt.gz`，两者必须匹配固定 ME5 SHA256。默认 Alopex 未安装时告警并标记未验证；显式根目录为空/不存在、已安装 Alopex 缺少孔板表、映射冲突或已有设计文件损坏均失败。Droplet 设计尚未准备时仅该项告警未验证（Lemmus 使用 `doctor.sh calling` 准备）；检查不构建环境、不运行 DNA doctor、不修改两侧资源。

Doctor 的 PASS 与告警未验证应分别核对；跨 Pipeline 检查通过仅证明共享身份资源一致，正式项目仍需通过输入、reference 与 DAG 检查。

## 4. 创建分析项目

Pipeline 安装目录与数据项目目录应分开。初始化必须显式选择协议，不能只靠 `expected_orientation: rc` 判断文库类型。同一份[内置配置模板](core/lemmus/config.py) 同时展示 `cabernet:` 和 `droplet:` 两个模块；顶部 `pipeline_mode` 只选择其中一组，参数不会跨模块继承。初始化写入这份完整配置，并设置所选模式。

### DD-MET5 Droplet 项目

```bash
mkdir -p /path/to/DropletProject
cd /path/to/DropletProject
/path/to/Lemmus/core/run_pipeline.sh --init-project --protocol droplet
```

生成含两个模块、默认选择 Droplet 的 `00_config/config.yaml`、`sample_manifest.tsv`；不生成 `Barcode_Map.csv`。运行前使用 `core/doctor.sh calling` 准备独立 calling 环境与官方 ME5 白名单。

### Cabernet 项目

```bash
mkdir -p /path/to/CabernetProject
cd /path/to/CabernetProject
/path/to/Lemmus/core/run_pipeline.sh --init-project --protocol cabernet
```

生成含两个模块、默认选择 Cabernet 的 `00_config/config.yaml`、`sample_manifest.tsv` 和 `Barcode_Map.csv`。RNA-only Cabernet 项目在 `cabernet.DARLIN.enabled` 处设为 false。

两种项目均生成指向 Pipeline 的 `run_pipeline.sh` 软连接以及 `01_raw/`、`02_work/`、`03_results/`、`04_logs/`、`05_tmp/` 和 `06_downstream/` 工作目录。后续命令从各自项目目录执行。重复初始化会保留已有配置；显式协议与已有配置不一致时拒绝操作，不覆盖结果或配置。`--protocol` 仅用于初始化，运行阶段以项目 config 为准。

## 5. 同一份 config.yaml 的两个模块

```yaml
pipeline_mode: droplet  # 选择生效模块；可选 cabernet 或 droplet
# schema_version、paths、raw_data_dir、sample_manifest、slurm 是共用项。
cabernet:
  # 孔板文库的 barcode_map、RNA、DARLIN、QC 参数
  ...
droplet:
  # DD-MET5 的 RNA、DARLIN、QC 参数
  ...
```

上面仅展示结构；可直接使用的完整配置见 [内置配置模板](core/lemmus/config.py)。只调整所选模块内的参数。下文的 `RNA.cores` 等简写均指所选模块中的对应字段，例如 `droplet.RNA.cores`。已有顶层 RNA/DARLIN/QC 的 schema-2 配置保持原有含义；同一文件不能混用顶层参数和模块参数，以免出现覆盖歧义。

### DD-MET5 Droplet 配置

完整模板的 `droplet:` 模块：[内置配置模板](core/lemmus/config.py)。R1 为 CB17+UMI12，R2 为 cDNA；不使用孔板映射。主要参数如下，其余工作目录和调度字段由模板提供：

```yaml
pipeline_mode: droplet
droplet:
  RNA:
    reference: hg38
    cores: 16
    mem_mb: 49152
    runtime: 720
    min_reads: 1
    expected_orientation: rc
    resolution: parsimony-gene
  DARLIN:
    enabled: false
  QC:
    read_fate:
      enabled: false
```

| 配置 | Droplet 含义 |
| --- | --- |
| `RNA.reference` | 按真实物种填写；模板为 `hg38`，reference 固定来自 Pipeline `resources/<species>_reference/pg45s/`。 |
| `RNA.cores` | local 总 CPU 预算，模板 16。 |
| `RNA.quant_threads` | 可选定量线程预算，缺省等于 cores，范围 `1..cores`。 |
| `RNA.mem_mb` / `RNA.runtime` | 定量资源，模板为 49152 MiB / 720 分钟；local 并发同时受 cores 和 mem_mb 约束。 |
| `RNA.min_reads` | **固定 1**，保留 emptyDrops 背景。 |
| `RNA.expected_orientation` | **固定 rc**，geometry 自动为 `1{b[17]u[12]x:}2{r:}`。 |
| `RNA.resolution` | **固定 parsimony-gene**，不运行多gene EM。 |
| `DARLIN.enabled` / `QC.read_fate.enabled` | **固定 false**；由修剪与原生定量提供配对/条码审计。 |

`droplet:` 模块不包含 `barcode_map`、DARLIN primer/template/extraction 参数或孔板 double-cell 参数。它们仅出现在 `cabernet:` 模块，选择 Droplet 时不会读取它们。

### Cabernet 配置

完整模板的 `cabernet:` 模块：[内置配置模板](core/lemmus/config.py)。R1 为 cDNA，R2 为 CB8/10+UMI8；运行前检查 `00_config/Barcode_Map.csv`。主要参数如下：

```yaml
pipeline_mode: cabernet
cabernet:
  barcode_map: 00_config/Barcode_Map.csv
  RNA:
    reference: mm10
    cores: 16
    mem_mb: 16384
    runtime: 720
    min_reads: 10
    expected_orientation: both
    resolution: parsimony-gene
  DARLIN:
    enabled: true
```

| 配置 | Cabernet 含义 |
| --- | --- |
| `RNA.reference` | 按真实物种填写；模板为 `mm10`，reference 固定来自 Pipeline `resources/<species>_reference/pg45s/`。 |
| `barcode_map` | canonical whitelist 与 PlateID 映射来源，模板为 `00_config/Barcode_Map.csv`。 |
| `RNA.cores` | local 总 CPU 预算，模板 16，不被 Slurm controller 的 local core 设置截断。 |
| `RNA.quant_threads` | 可选定量线程预算，缺省等于 cores，范围 `1..cores`；显式调低可留出并发 DARLIN 所需的核。 |
| `RNA.mem_mb` / `RNA.runtime` | 定量资源，模板为 16384 MiB / 720 分钟。 |
| `RNA.min_reads` | canonical whitelist 的原生 permit-list read 门槛，默认 10。 |
| `RNA.expected_orientation` | 接受 `fw` / `rc` / `both`，默认 `both`；改成 rc 不会改变 Cabernet 读端结构。 |
| `RNA.resolution` | 固定 `parsimony-gene`，其它值拒绝。 |
| `DARLIN.enabled` | 默认 true；没有 targeted DARLIN 文库时设 false。 |

DARLIN 提取参数（`DARLIN.extraction.chunk_size` 等）的默认值集中定义在 `lemmus.config.DEFAULT_EXTRACTION`，由 Snakefile 导入使用；`cabernet:` 模块列出常用项，一般项目无需调整。

整份配置使用 `schema_version: 2`，不允许项目级 `reference_dir`、`simpleaf_reference`、`genome_fasta` 或 `gtf` 覆盖。FASTQ 导入后刷新并检查 `00_config/sample_manifest.tsv`。

### 项目契约

生产 Pipeline 只接受当前 schema 的配置与 manifest，不在运行时自动迁移旧 schema、旧列名或旧路径。旧项目应先在项目外完成显式升级，再进入生产运行。

### Barcode_Map.csv

Cabernet 的 `Barcode_Map.csv` 是 canonical cell-barcode whitelist 与 PlateID 映射来源：

- 列固定为 `DNA_Barcode`、`RNA_Barcode`、`PlateID`、`Cell_Order`；
- RNA Cell Barcode 长度由该表自动推断，只允许 8 bp 或 10 bp；UMI 固定 8 bp；
- barcode 只能包含 `A/C/G/T`；
- 同长度 barcode 的最小两两 Hamming distance 必须 > 2，否则加载即拒绝（保证 unique Hamming-1 纠错无歧义）；
- `PlateID` 和 `Cell_Order` 必须唯一。

对应 simpleaf custom geometry 由 Pipeline 自动生成，不在规则里维护第二套 barcode 定义：

```text
8-bp CB : 1{r:}2{b[8]u[8]x:}
10-bp CB: 1{r:}2{b[10]u[8]x:}
```

## 6. 准备 raw FASTQ 和 sample manifest

将 paired FASTQ 放入项目的 `01_raw/`：

```text
01_raw/RNA/<sample>_R1.fastq.gz      01_raw/RNA/<sample>_R2.fastq.gz
01_raw/DARLIN/<sample>_R1.fastq.gz   01_raw/DARLIN/<sample>_R2.fastq.gz
```

RNA 和 DARLIN 均递归扫描子目录，支持 `<sample>_1.fq.gz` / `<sample>_2.fq.gz`、`<sample>_R1[_001].fastq.gz` / R2、`<sample>_L001_R1_001.fastq.gz` / R2 和 `<sample>.R1.raw.fastq.gz` / R2；扩展名均可为 `.fq.gz` 或 `.fastq.gz`。无需重命名供应商文件；移除上述读端后缀后得到原始文库名。每个文库必须恰好一对文件，缺 mate、重复读端（包括多个 lane）或不支持的 FASTQ 名称均报错。

Cabernet 文库为 R1 = cDNA、R2 = Cell Barcode + UMI；DD-MET5 Droplet 为 R1 = CB17 + UMI12、R2 = 表达序列。

导入 FASTQ 后运行：

```bash
/path/to/Lemmus/core/run_pipeline.sh --refresh-manifest
/path/to/Lemmus/core/run_pipeline.sh --dry-run
```

`--refresh-manifest` 根据 `01_raw/` 更新 `00_config/sample_manifest.tsv`；`--dry-run` 检查配置、FASTQ、当前协议的白名单、reference 和待执行任务，不会正式分析数据。

## 7. 运行 Pipeline

macOS 和 Linux/HPC 使用相同的用户命令：

```bash
/path/to/Lemmus/core/run_pipeline.sh --dry-run
/path/to/Lemmus/core/run_pipeline.sh
```

- macOS 自动在本机以 local executor 运行，总计算宽度直接来自 `RNA.cores`；
- Linux/HPC 使用 Snakemake Slurm executor 调度 rule 作业；`slurm.local_cores` 只用于 Slurm controller；
- 运行失败后修正问题重新执行同一命令即可；`02_work/` 与 `.snakemake/` 保留并用于增量恢复，由 Snakemake 根据当前输入、参数、代码和输出状态决定最小重算集合。

RNA 生产链只有一个 backend：

```text
raw R1/R2
   │
   ├─ FastQC（原始 reads QC）
   │
   └─ simpleaf quant
        ├─ piscem map-sc
        ├─ canonical whitelist
        ├─ unique + Hamming-1 cell-barcode correction
        ├─ alevin-fry collate / molecule resolution
        └─ PG AnnData → 分离45S技术计数 → biological.h5ad → RNA.h5ad
```

不再执行 RNA 的 `umi_tools extract → Python trim → Cutadapt → STAR → featureCounts → samtools → umi_tools count` 多阶段链。RNA barcode correction 使用 alevin-fry 的 deterministic `unique + hamming-1` policy；孔板路线负责 Barcode Map 的长度和 Hamming ambiguity 审计；Droplet 核验固定设计白名单的哈希与唯一性，不做全表两两距离检查。两者均不维护第二套 correction table。

DARLIN targeted library 保持独立处理：按 sample 对 paired FASTQ 单次 streaming 解析，以 `DARLIN.extraction.chunk_size`（默认 100,000 reads）分块同时评估全部配置的 templates，不为每个 template 重复读取 FASTQ。流式内存有界（read-ID inventory 流、8 字节重复检测哈希、sidecar 按 read ordinal 对接），常驻内存只随分子多样性而非 read 总数增长。提取 rule 前置 `assert-fastq` 门禁，损坏/截断的 FASTQ 秒级失败。primer-variant 决策（`use_short_3prime`）每 sample/template 只做一次并在后续 chunk 不一致时 fail closed。Cell barcode 与 RNA 路径统一：exact match 优先，其后仅允许 unique Hamming-1 rescue，多个 whitelist 邻居的 barcode fail-closed。**UMI 级等位基因共识**在 read cutoff 之前把同一 cell+UMI 的近同源序列变体合并到同一分子（read 支持相加）并同步重映射 read-level 证据；Clone-sequence directional clustering 使用 `umi_tools.UMIClusterer`，但严格限定在单个 biological cell 内。`umi_tools` 依赖仅用于 DARLIN，不参与 RNA quantification。

## 8. 数据产出与完成确认

Cabernet 的主要结果位于下表。Droplet 的正式矩阵同样为 `03_results/count/<sample_id>/RNA.h5ad`，原始定量和 mapping 指标位于 `02_work/rna_core/droplet_quant/<sample_id>/`，不产生 DARLIN 输出：

```text
03_results/
├── count/<sample_id>/
│   ├── RNA.h5ad                      正式项目级 RNA 矩阵
│   ├── af_quant/alevin/biological.h5ad  持久生物gene矩阵（45S技术计数在obs）
│   ├── af_quant/quant.json
│   └── af_map/map_info.json
├── DARLIN/DARLIN_clone_id/DARLIN_Clone_ID.csv   DARLIN 最终 clone table
└── QC_Results/
    ├── multiqc_report.html
    ├── multiqc_data/
    └── rna_darlin_cell_alignment.tsv            RNA/DARLIN 细胞对齐审计
```

Cabernet 的 `RNA.h5ad.obs_names` 固定为 `<sample_id>_<PlateID>`；Droplet 固定为 `<sample_id>_<barcode>`，并保留 `RNA_Barcode`、`PlateID`、`Cell_Order`、`sample_id`、`rna_sample` 元数据；`library` 固定等于 manifest 中的 `sample_id`，用于下游分组和图例，`rna_sample` 仅保存原始 FASTQ/library 名称。输出 barcode 不在对应协议的 whitelist 时 Pipeline fail closed。Cabernet 另须拒绝两个 barcode 映射到同一 `<sample_id>_<PlateID>`，其 `RNA.h5ad.uns` 与 `DARLIN_Clone_ID.csv` 携带 `barcode_map_sha256`；Droplet 在 `RNA.h5ad.uns` 保留固定设计白名单身份及独立调用证据。

`parsimony-gene` 输出原生整数分子数，不做多gene EM或手动取整。`n_genes_supported`（至少1分子）是正式RNA QC主复杂度指标；positive/fractional/subunit审计字段保留其数值定义，正常PG输出的fractional/subunit应为0。技术45S计数不会进入表达分母。`quants.h5ad` 仅作为投影输入的临时工件，成功写出 `biological.h5ad` 后由Snakemake回收；后者保留全部背景barcode供独立calling和恢复。

MultiQC 汇总五类证据：FastQC 由 MultiQC 原生模块解析；`RNA.h5ad` 生成 RNA expression summary；simpleaf 的 `map_info.json` / `quant.json` 以保留原始 JSON key path 的 custom-content table 展示；启用 DARLIN 时加入最终 clone table 的 DARLIN QC section、barcode correction 指标（exact / Hamming-1 corrected / ambiguous / unmatched reads 与 correction fraction）和 RNA/DARLIN cell alignment 审计。

下游合并必须保留 `RNA.h5ad.var` 中的 `gene_id` 与 `gene_symbol`：版本化 Ensembl ID 作为稳定索引，MT/IG QC 和 marker 解释使用 `gene_symbol`。主注释 H5AD 通过同目录临时文件原子发布，写入失败不会破坏旧结果。

最终必需的非空文件：

```text
03_results/QC_Results/multiqc_report.html
03_results/count/<sample_id>/RNA.h5ad
```

启用 DARLIN 时，另须生成 `03_results/DARLIN/DARLIN_clone_id/DARLIN_Clone_ID.csv`。

## 下游分析

下游 Processor 负责计算，Notebook 负责调参与可视化。使用 Doctor 管理的 `conda/current/atlas` 环境，输入为 finalized `03_results/count/<sample_id>/RNA.h5ad`；细胞身份已在 finalize 阶段生成，Notebook 不重复重命名。

### QC 与完整注释

```bash
PIPELINE=/path/to/Lemmus
# 生成 RNA/DARLIN 每细胞 QC 表
PYTHONPATH="$PIPELINE/core" "$PIPELINE/conda/current/atlas/bin/python" \
  -m lemmus.downstream.qc --project-dir /path/to/MyProject

# 完整分析：QC、聚类、marker 与物种分支注释
PYTHONPATH="$PIPELINE/core" "$PIPELINE/conda/current/atlas/bin/python" \
  -m lemmus.downstream.qc --project-dir /path/to/MyProject --analyze
```

需要 DNAme 整合时，在命令中加入 `--dna-project-root /path/to/DNA_project`。完整注释也可运行 [RNA_QC_Visualization.ipynb](downstream/RNA_QC_Visualization.ipynb)，通过 `QC_THRESHOLDS` / `ANALYSIS_OVERRIDES` 配置 Processor。

| 输出（`06_downstream/`） | 使用边界 |
| --- | --- |
| `RNA_QC_Information.csv` | 固定13列；`CloneID` 为 RNA library 身份，未运行分析或未通过 QC 时注释列为空 |
| `RNA_cell_qc_metrics.csv` | 全细胞 QC 与 `qc_pass`，供展示过滤前分布与阈值 |
| 注释 H5AD | mouse 为 `RNA_mouse_marker_annotated.h5ad`，human 为 `RNA_Census_consensus_annotated.h5ad`；保留全细胞 raw counts，聚类/注释/embedding 只填 QC 通过细胞，被滤细胞注释为空、embedding 为 NaN |
| marker / 审计 CSV | 支持注释解释与检查；完整 DARLIN 证据仍在 `03_results/DARLIN/DARLIN_clone_id/DARLIN_Clone_ID.csv` |
| `analysis_summary.json` | 实际参数、seed、依赖版本与输入哈希 |

human 使用 Census scVI consensus；小鼠使用 deterministic marker-program 注释，低支持保留 `Uncertain` / `Unresolved`。正式 broad 标签为 `cell_type_final`；样品名与来源用于审计，不作为注释先验。

Droplet QC-only CLI 支持 `--min-genes-per-cell`、`--min-counts-per-cell`、`--max-mito-pct`。阈值只更新下游 `qc_pass`，不修改分子计数矩阵或重新调用细胞；未传参数时沿用已保存的 QC 设置。已有聚类/注释结果时，修改阈值会撤销失效注释及其完成标记；需要注释时再用 `--analyze` 在原位置重建。实际参数以 Notebook 配置和 Processor 的 `--help` 为准。

### Droplet QC 与 doublet Notebook

打开 [RNA_Droplet_QC_Visualization.ipynb](downstream/RNA_Droplet_QC_Visualization.ipynb)，填写 `RNA_PROJECT_ROOT`、`DNA_PROJECT_ROOT`。RNA 须已有 finalized RNA.h5ad，DNA 须已有当前 complete delivery、同代 QC 与所需矩阵。

| 操作 | 行为 |
| --- | --- |
| RNA doublet | 逐文库运行 Scrublet，Notebook 的 `DOUBLET_THRESHOLD_SCALE` 默认0.5；整合 QC 并生成 RNA singlet AnnData，不运行聚类、Census、marker 或注释 |
| DNA doublet | `DNA_DOUBLET_THRESHOLD_SCALE` 默认0.5，使用同代 RawAdata 的原始 mc/cov 独立分析全部 DNA；只添加判定，不改变 RNA singlet 筛选 |
| 调整倍率 | 实际阈值为自动阈值×倍率，严格大于才判阳性；减小倍率会增加剔除，也可能误标正常细胞。失败或低信息不当作 singlet |
| 最后一格 | 导出全部 DNA 的 SingleCpG 副本，保留原矩阵/身份/注释，加入 DNA20、`RNA_` 前缀 RNA13 和独立/联合 doublet 判定；DNA-only 的 RNA QC 留空 |

整合结果位于 `06_downstream/droplet_integration/`，DNA 判定与 SingleCpG 位于其 `dna_doublet/`。summary 最后提交；调整输入或参数后覆盖固定位置的对应结果。续跑 DNA doublet/SingleCpG 导出前会复核当前 sealed 文库关联，失配时须重新整合。配置入口见上述 Droplet Notebook。

### DNA/RNA 文库关联

在 DNA manifest 的可选 `rna_sample` 填 RNA_DARLIN **原始 RNA basename**。例如 DNA `sample_id=DNA_alias, rna_sample=RAW_RNA`，RNA `sample_id=RNA_alias, rna_sample=RAW_RNA`；两侧项目样本名可独立命名。RNA manifest 的 `sample_id` 必须与已完成 RNA.h5ad 一致，只修改 manifest 不能重命名结果。

| 协议 | 连接方式 |
| --- | --- |
| Cabernet / TAPS | sealed RNA raw → RNA canonical sample_id，再按 PlateID 一对一对齐；两侧参与 QC 的细胞集合须一致，身份/条码/Cell_Order 或列名冲突失败 |
| SRD | DNA `rna_sample` 是本地 RNA tube，不用于外部桥接；按共享项目样本名与 PlateID 连接 |
| Droplet | 通过 sealed `rna_sample` 连接文库，再按完整 CB17 精确连接；留空保持未关联，双方独有细胞保留在 `DNA_RNA_cell_links.csv`，不补零表达 |

关联只读当前 complete DNA 的 sealed manifest 和同代 QC，不读取可变 config 或按 mtime 挑目录。修改关联后须重新发布 DNA、生成同代 downstream QC，再重跑 RNA 整合。`Integrated_QC_Information.csv` 保留 RNA 身份与实际映射溯源；RNA13/DNA20 的独立 QC schema 不变。

## 常见问题

- Doctor 报 reference 缺失：确认 FASTA/GTF 路径后重新运行 `core/doctor.sh reference <species> --fasta ... --gtf ...`。
- 找不到 paired FASTQ：确认 RNA / DARLIN FASTQ 分别位于 `01_raw/RNA/` 与 `01_raw/DARLIN/`，且每个 R1 都有对应 R2。
- 配置校验失败：根据错误修改 `00_config/config.yaml`；不要添加模板中不存在的字段（尤其不得添加 reference 路径覆盖）。
- Cabernet Barcode_Map 校验失败：检查 barcode 长度（8/10 bp）、字符集与最小 Hamming distance > 2。
- `assert-fastq` 失败：DARLIN FASTQ 损坏或截断，替换文件后重跑。
- 分析失败：先查看 `04_logs/` 中对应 rule 的日志，修正问题后重新运行 `core/run_pipeline.sh`。

## 验证

```bash
core/doctor.sh test
```

会执行内存语法检查、Shell syntax、配置/schema 和安装运行自检，不运行开发全仓或性能回归；受管环境存在时还会执行 Snakemake full-DAG dry-run。

## DD-MET5 RNA 快速配置

依次完成[初始化 Droplet 项目](#dd-met5-droplet-项目)、[检查 Droplet 配置](#dd-met5-droplet-配置)、[导入 FASTQ](#6-准备-raw-fastq-和-sample-manifest)及 dry-run，再正式运行。RNA 独立 calling 与定量完成后，按[下游关联步骤](#dnarna-文库关联)准备 DNA sealed 关联和同代 QC，运行 CLI 或专用 Notebook。
