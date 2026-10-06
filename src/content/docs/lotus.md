---
title: 'Lotus T2T 分析流程'
sidebar:
  label: 'Lotus T2T 分析流程'
  order: 0
---


## 中文

Gifu 与 MG20 两个 *Lotus japonicus* 材料的端粒到端粒（T2T）基因组组装、注释与下游分析代码。流程按 **Gifu 实际运行版本**呈现，MG20 用对应输入执行相同步骤；每个步骤目录自带 README，说明该步的做法、输入与参数。本仓库所有 README 均为中英双语，中文在前、英文在后。

### 仓库结构与设计规则

顶层：`scripts/`（分析流程，按执行顺序编号）、`environment/`（软件版本记录）、`third_party/`（他人代码，保留署名与许可）、`demo/`（可运行小示例与期望输出）、`LICENSE`（MIT，仅自研代码）、`.gitignore`。

**设计规则：目录顺序 = 执行顺序 = 阅读顺序**，`scripts/` 从 `01_assembly` 到 `08_single_nucleus` 依次执行。

### `scripts/` 的组织

分两级：**阶段**（主要分析块）与**步骤**（一次工具调用），数字前缀标明层级内的顺序。

```text
scripts/
├── 01_assembly/          # 01_hifiasm, 02_purge_haplotigs, 03_remove_organelle, 04_telomere_check,
│                         # 05_contig_ordering, 06_gapclose, 07_polish, 08_assembly_qc, 09_hic_scaffold
├── 02_annotation/        # 01_repeat_edta, 02_rnaseq_align, 03_isoseq, 04_braker, 05_geta, 06_gene_model_filter,
│                         # 07_gene_classification, 08_functional_annotation, 09_annotation_export, 10_gene_correspondence
├── 03_telomere/          # 端粒基序识别、覆盖度、跨端粒读长
├── 04_rDNA/              # rDNA 位点鉴定与阵列结构
├── 05_centromere/        # 串联重复聚类、超级簇合并
├── 06_genome_comparison/ # 全基因组比对、SV 鉴定、重排验证
├── 07_transcriptome/     # bulk RNA-seq、共表达网络、候选优先级
└── 08_single_nucleus/    # snRNA-seq 质控、聚类、注释、虚拟敲除
```

`03_telomere`、`04_rDNA`、`05_centromere` 各自独立成阶段；`06_genome_comparison` 使用完成的组装，需两个材料都组装完成后运行；`08_single_nucleus` 与 bulk RNA-seq 的数据模态和软件栈不同。步骤目录约定：Gifu 为示例，脚本保留具体文件名与设置；无配置层，输入输出路径直接写在脚本中；方法与参数记录在步骤 README；保留原始文件名（部分脚本互相按名调用）；只保留最终版本。

### 两个源码环境

| 阶段 | 项目目录 |
|---|---|
| `01_assembly` … `05_centromere` | `{PROJ_04_LOTUS_GENOME}`，部分注释脚本仍指向更早的 ``{PROJ_LOTUS_ZC}`` |
| `06_genome_comparison` … `08_single_nucleus` | ``{PROJ_LOTUS_HX}`` |

### 其他顶层目录

- `environment/` — `versions.tsv` 每个工具一行：版本、用途与来源；`#CSUB` / `#SBATCH` 提交指令保留在脚本顶部。
- `third_party/` — 他人代码，保留原作者署名与许可，列于 `third_party/THIRD_PARTY_NOTICES.md`。
- `demo/` — 按阶段组织为 `demo/<stage>/{input,expected_output,run_demo.sh}`；现有 `demo/03_telomere/` 与 `demo/05_centromere/`，仅需 Python 3，不到一秒完成，无需第三方包、参考基因组或测序数据。

### 许可说明

- 仓库根目录的 MIT 许可只覆盖**自研脚本**。
- `third_party/` 中的代码**不重新许可**，保留原作者与条款，见 `third_party/THIRD_PARTY_NOTICES.md`。
- 仅被调用的第三方工具（hifiasm、EDTA、BRAKER、GETA、HiC-Pro、BUSCO、Merqury、StainedGlass、SyRI 等）不随仓库分发，版本与来源见 `environment/versions.tsv`。

## English

Assembly, annotation and downstream analysis code for the telomere-to-telomere (T2T) genomes of the *Lotus japonicus* accessions Gifu and MG20. The pipeline is presented **as it was run on the Gifu accession**; the same steps were applied to MG20 with the corresponding inputs. Each step directory carries its own README describing what the step does, its inputs and its parameters. All READMEs in this repository are bilingual, Chinese first and English second.

### Repository layout and design rule

Top level: `scripts/` (the pipeline, numbered in execution order), `environment/` (software version records), `third_party/` (code by others, with attribution and licence), `demo/` (small runnable examples with expected output), `LICENSE` (MIT, our own code only), `.gitignore`.

**Design rule: directory order = execution order = reading order**; `scripts/` runs from `01_assembly` to `08_single_nucleus`.

### How `scripts/` is organised

Two levels: **stage** (a major analysis block) and **step** (one tool invocation); numeric prefixes make the order explicit inside each level.

```text
scripts/
├── 01_assembly/          # 01_hifiasm, 02_purge_haplotigs, 03_remove_organelle, 04_telomere_check,
│                         # 05_contig_ordering, 06_gapclose, 07_polish, 08_assembly_qc, 09_hic_scaffold
├── 02_annotation/        # 01_repeat_edta, 02_rnaseq_align, 03_isoseq, 04_braker, 05_geta, 06_gene_model_filter,
│                         # 07_gene_classification, 08_functional_annotation, 09_annotation_export, 10_gene_correspondence
├── 03_telomere/          # telomere motif calling, coverage, spanning reads
├── 04_rDNA/              # rDNA locus identification and array structure
├── 05_centromere/        # tandem-repeat clustering, supercluster merging
├── 06_genome_comparison/ # whole-genome alignment, SV calling, rearrangement validation
├── 07_transcriptome/     # bulk RNA-seq, co-expression networks, candidate prioritisation
└── 08_single_nucleus/    # snRNA-seq QC, clustering, annotation, virtual knockout
```

`03_telomere`, `04_rDNA` and `05_centromere` are separate stages; `06_genome_comparison` consumes finished assemblies and runs after both accessions are assembled; `08_single_nucleus` differs from bulk RNA-seq in data modality and software stack. Conventions inside a step directory: Gifu is the worked example and scripts keep the concrete filenames and settings used; there is no configuration layer, so paths appear directly in the script; the method and parameters live in the step README; original filenames are kept because several scripts invoke one another by name; only final script versions are kept.

### Two source environments

| Stages | Project tree |
|---|---|
| `01_assembly` … `05_centromere` | `{PROJ_04_LOTUS_GENOME}`, with some annotation scripts still pointing at the earlier ``{PROJ_LOTUS_ZC}`` |
| `06_genome_comparison` … `08_single_nucleus` | ``{PROJ_LOTUS_HX}`` |

### Other top-level directories

- `environment/` — `versions.tsv`, one row per tool: version, purpose and source; `#CSUB` / `#SBATCH` submission directives are left at the top of the scripts where they were used.
- `third_party/` — code written by others, keeping original author attribution and licence, listed in `third_party/THIRD_PARTY_NOTICES.md`.
- `demo/` — organised per stage as `demo/<stage>/{input,expected_output,run_demo.sh}`; the demos need Python 3 only, finish in under a second, and require no third-party packages, reference genome or sequencing data.

### Notes on licensing

- The repository-root MIT licence covers **our own scripts only**.
- Code in `third_party/` is **not relicensed**; it keeps its original authors and terms, see `third_party/THIRD_PARTY_NOTICES.md`.
- Third-party tools that are merely *called* (hifiasm, EDTA, BRAKER, GETA, HiC-Pro, BUSCO, Merqury, StainedGlass, SyRI, …) are not redistributed here; versions and sources are listed in `environment/versions.tsv`.

### 本阶段步骤

| 步骤 | 内容 |
| --- | --- |
| [Demo](/lotus/demo/) | — |
| [Scripts](/lotus/scripts/) | — |


