---
title: '分析流程总览'
sidebar:
  label: '分析流程总览'
  order: -1
---


## 中文

按执行顺序编号的分析流程。每个步骤目录保存该步 **Gifu 实际运行的脚本**，以及说明方法、输入、输出、参数与所用工具的 `README.md`；脚本保留原始文件名（部分脚本互相按名调用），步骤 README 指明该步的入口脚本。

### 阶段索引

| 阶段 | 步骤数 | 内容 |
|---|---|---|
| [`01_assembly/`](/_projects/lotus/scripts/01-assembly/) | 9 | 基因组组装、打磨、质控 |
| [`02_annotation/`](/_projects/lotus/scripts/02-annotation/) | 10 | 重复屏蔽、基因预测、过滤、基因分类、功能注释 |
| [`03_telomere/`](/_projects/lotus/scripts/03-telomere/) | — | 端粒基序识别与末端验证 |
| [`04_rDNA/`](/_projects/lotus/scripts/04-rdna/) | 2 | rDNA 位点鉴定、阵列结构、拷贝数 |
| [`05_centromere/`](/_projects/lotus/scripts/05-centromere/) | 4 | 串联重复聚类与着丝粒界定 |
| [`06_genome_comparison/`](/_projects/lotus/scripts/06-genome-comparison/) | 5 | 全基因组比对、结构变异、PAV、重排验证 |
| [`07_transcriptome/`](/_projects/lotus/scripts/07-transcriptome/) | 7 | bulk RNA-seq、组织图谱、WGCNA、Mfuzz、共表达 |
| [`08_single_nucleus/`](/_projects/lotus/scripts/08-single-nucleus/) | — | snRNA-seq 质控、整合、聚类、注释、虚拟敲除 |

### 执行顺序
阶段按编号顺序执行，`01_assembly` 内部顺序为 `01_hifiasm → 02_purge_haplotigs → 03_remove_organelle → 04_telomere_check → 05_contig_ordering → 06_gapclose → 07_polish → 08_assembly_qc → 09_hic_scaffold`；`02_annotation` 使用 `01_assembly` 产出的完成打磨的染色体级组装，`03_telomere`、`04_rDNA`、`05_centromere` 同样作用于完成组装的基因组且彼此独立，`06_genome_comparison` 需要两个材料都完成组装，`07_transcriptome` 与 `08_single_nucleus` 以注释为参考。

### 两个源码环境
脚本来自两位作者各自搭建的工作树，内部路径体现了这一点；工具版本见 `environment/versions.tsv`。

| 阶段 | 项目目录 | Conda 环境 |
|---|---|---|
| `01_assembly`、`02_annotation`、`03_telomere`、`04_rDNA`、`05_centromere` | `{PROJ_04_LOTUS_GENOME}`（注释脚本同时引用 ``{PROJ_LOTUS_ZC}``） | `py3`、`geta-2.7.1`、singularity 镜像 |
| `06_genome_comparison`、`07_transcriptome`、`08_single_nucleus` | ``{PROJ_LOTUS_HX}`` | `lotus`、`r_wgcna`、`singlecell` |

### 基因组命名

| 脚本中的名称 | 指代 |
|---|---|
| `Gifu_v1.0.fasta` | **Gifu T2T 组装**（早期名称） |
| `Lotus_GifuT2T_v1.0.fasta` | 同一套 Gifu T2T 组装（后期名称） |
| `Gifu_v1.0.fasta.mod.*`、`...mod.MAKER.masked` | 由早期基名派生的 EDTA 输出 |
| `MG20_v1.0.fasta`、`Lotus_MG20T2T_v1.0.fasta` | **MG20 T2T 组装**，同样情况 |

同一套组装在项目不同时期被改名，**并非不同的基因组**，因此重复屏蔽、比对与注释脚本会混用这些名称，脚本可能用 `Gifu_v1.0.fasta` 建索引、同时读取 `Lotus_GifuT2T` 注释。`ChiDou` 与 `ChiHei` 是**大豆品种名**，不是 Lotus 材料：引用它们的脚本或属大豆分析（`07_transcriptome` 中的 `Gm_*` 步骤），或为 `02_annotation/01_repeat_edta/LAI.sh` 这类由大豆流程改写的脚本，用于 Gifu 前需替换基因组文件名。

### 路径占位符
**脚本中没有机器相关的绝对路径**，所有生产路径已替换为 ``{TOKEN}`` 占位符；全树做一次查找替换成自己的位置即可运行。占位符是纯文本，不是 shell 或 Python 语法：``{TOKEN}`` 刻意不写作 `$`{TOKEN}``，以免被误认为已定义的变量。

| 占位符 | 含义 |
|---|---|
| `{PROJ_04_LOTUS_GENOME}` | 基因组组装项目树（`04_Lotus_genome`） |
| ``{PROJ_LOTUS_ZC}`` | 该作者更早的工作树（`Project/Lotus`） |
| ``{PROJ_LOTUS_HX}`` | RNA-seq / 单核项目树（`Project/Lotus`，第二位作者） |
| ``{CONDA_PY3}`` | conda 环境 `py3` |
| ``{CONDA_SINGLECELL}``、``{CONDA_LOTUS}``、``{CONDA_R_WGCNA}`` | 转录组与单核工作使用的三个 conda 环境 |
| ``{CONDA_ZC}``、``{CONDA_HX}`` | 存放上述两个环境的两套 conda 安装 |
| ``{SOFTWARE_ZC}``、``{SOFTWARE_HX}``、``{SOFTWARE_LOCAL}`` | 工具安装目录（BRAKER、HapHiC、StainedGlass 等） |
| ``{HOME_ZC}``、``{HOME_HX}``、``{HOME_LOCAL}`` | 账号家目录 |

- **只有上表 15 个 token 是占位符。** 脚本本身还含有 `{sample}`、`{chrom}`、`{SP_NAME}`、`{WD}`、`{THREADS}`、`{PREFIX}`、`{label}`、`{n}` 等 `{...}` 表达式，属于作者自己的格式字符串变量或 shell heredoc 模板，在文件内定义与使用，**不要替换**；查找替换只替换上表列出的 token。
- `#CSUB`、`#SBATCH` 调度器指令、shebang 中的 `/usr/bin/env`、计时封装中的 `/usr/bin/time`、以及 `~` 形式的家目录路径（`~/anaconda3/etc/profile.d/conda.sh`、`~/software/StainedGlass-0.6/`）均**保持原样**，它们不是机器相关路径；tilde 形式可移植，同时假定 `~/software` 布局。

### 相关目录
相关目录：`environment/versions.tsv`（工具版本、用途与来源）、`third_party/`（随仓库分发的第三方代码及其许可）、`demo/`（可运行小示例与期望输出）。

<details>
<summary>English</summary>
The analysis pipeline, numbered in execution order. Each step directory holds the scripts for that step **as they were run on the Gifu accession**, together with a `README.md` describing the method, inputs, outputs, parameters and the tools used; scripts keep their original filenames because several of them invoke one another by name, and the step README identifies the entry point for that step.

### Stage index

| Stage | Steps | What it covers |
|---|---|---|
| [`01_assembly/`](/_projects/lotus/scripts/01-assembly/) | 9 | genome assembly, polishing, quality control |
| [`02_annotation/`](/_projects/lotus/scripts/02-annotation/) | 10 | repeat masking, gene prediction, filtering, gene classification, functional annotation |
| [`03_telomere/`](/_projects/lotus/scripts/03-telomere/) | — | telomere motif calling and terminal validation |
| [`04_rDNA/`](/_projects/lotus/scripts/04-rdna/) | 2 | rDNA locus identification, array structure, copy number |
| [`05_centromere/`](/_projects/lotus/scripts/05-centromere/) | 4 | tandem-repeat clustering and centromere delimitation |
| [`06_genome_comparison/`](/_projects/lotus/scripts/06-genome-comparison/) | 5 | whole-genome alignment, structural variation, PAV, rearrangement validation |
| [`07_transcriptome/`](/_projects/lotus/scripts/07-transcriptome/) | 7 | bulk RNA-seq, tissue atlas, WGCNA, Mfuzz, co-expression |
| [`08_single_nucleus/`](/_projects/lotus/scripts/08-single-nucleus/) | — | snRNA-seq QC, integration, clustering, annotation, virtual knockout |

### Execution order
Stages run in numeric order; inside `01_assembly` the order is `01_hifiasm → 02_purge_haplotigs → 03_remove_organelle → 04_telomere_check → 05_contig_ordering → 06_gapclose → 07_polish → 08_assembly_qc → 09_hic_scaffold`. `02_annotation` consumes the polished, chromosome-level assembly produced by `01_assembly`, `03_telomere`, `04_rDNA` and `05_centromere` also operate on the finished assembly and can run independently of one another, `06_genome_comparison` requires both accessions to be assembled, and `07_transcriptome` and `08_single_nucleus` use the annotation as their reference.

### Two source environments
The scripts come from the working trees of two co-authors who set up their analyses separately, and their internal paths reflect that; tool versions are listed in `environment/versions.tsv`.

| Stages | Project tree | Conda environments |
|---|---|---|
| `01_assembly`, `02_annotation`, `03_telomere`, `04_rDNA`, `05_centromere` | `{PROJ_04_LOTUS_GENOME}` (annotation scripts also reference ``{PROJ_LOTUS_ZC}``) | `py3`, `geta-2.7.1`, singularity images |
| `06_genome_comparison`, `07_transcriptome`, `08_single_nucleus` | ``{PROJ_LOTUS_HX}`` | `lotus`, `r_wgcna`, `singlecell` |

### Genome naming

| Name seen in the scripts | Refers to |
|---|---|
| `Gifu_v1.0.fasta` | the **Gifu T2T assembly** (earlier name) |
| `Lotus_GifuT2T_v1.0.fasta` | the same Gifu T2T assembly (later name) |
| `Gifu_v1.0.fasta.mod.*`, `...mod.MAKER.masked` | EDTA outputs derived from the earlier basename |
| `MG20_v1.0.fasta`, `Lotus_MG20T2T_v1.0.fasta` | the **MG20 T2T assembly**, same situation |

The same assemblies were renamed at different periods of the project and **are not different genomes**, so repeat-masking, alignment and annotation scripts mix these names freely, and a script may build an index from `Gifu_v1.0.fasta` while reading a `Lotus_GifuT2T` annotation. `ChiDou` and `ChiHei` are **soybean cultivar names**, not Lotus accessions: scripts that reference them are either soybean analyses (the `Gm_*` steps in `07_transcriptome`) or, as with `02_annotation/01_repeat_edta/LAI.sh`, adapted from a soybean run with the genome filename substituted before use on Gifu.

### Path placeholders
**No script contains a machine-specific absolute path**: every production path has been replaced by a ``{TOKEN}`` placeholder, and one search-and-replace over the whole tree substitutes your own locations. Placeholders are plain text, not shell or Python syntax: ``{TOKEN}`` is deliberately not `$`{TOKEN}``, so a bare token cannot be mistaken for a variable that is already defined.

| Token | Stands for |
|---|---|
| `{PROJ_04_LOTUS_GENOME}` | the genome-assembly project tree (`04_Lotus_genome`) |
| ``{PROJ_LOTUS_ZC}`` | the earlier working tree of the same author (`Project/Lotus`) |
| ``{PROJ_LOTUS_HX}`` | the RNA-seq / single-nucleus project tree (`Project/Lotus`, second author) |
| ``{CONDA_PY3}`` | conda environment `py3` |
| ``{CONDA_SINGLECELL}``, ``{CONDA_LOTUS}``, ``{CONDA_R_WGCNA}`` | the three conda environments used for the transcriptome and single-nucleus work |
| ``{CONDA_ZC}``, ``{CONDA_HX}`` | the two conda installations that hold those environments |
| ``{SOFTWARE_ZC}``, ``{SOFTWARE_HX}``, ``{SOFTWARE_LOCAL}`` | tool installation directories (BRAKER, HapHiC, StainedGlass, …) |
| ``{HOME_ZC}``, ``{HOME_HX}``, ``{HOME_LOCAL}`` | the account home directories |

- **Only the 15 tokens in the table above are placeholders.** Several scripts already contained `{...}` expressions of their own — `{sample}`, `{chrom}`, `{SP_NAME}`, `{WD}`, `{THREADS}`, `{PREFIX}`, `{label}`, `{n}` and others. Those are the authors' format-string variables or shell heredoc templates, defined and used inside the file; **do not substitute them** — a search-and-replace replaces only the tokens listed above.
- `#CSUB` / `#SBATCH` scheduler directives, `/usr/bin/env` in shebangs, `/usr/bin/time` in the timing wrappers, and home-relative paths written with `~` (`~/anaconda3/etc/profile.d/conda.sh`, `~/software/StainedGlass-0.6/`) were all **left untouched**: they are not machine-specific; the tilde form is portable and assumes a `~/software` layout.

### Related directories
Related directories: `environment/versions.tsv` (tool versions, purposes and sources), `third_party/` (third-party code redistributed with this repository, and its licences), `demo/` (small runnable examples with expected output).

</details>

### 本阶段步骤

| 步骤 | 内容 |
| --- | --- |
| [01 组装](/_projects/lotus/scripts/01-assembly/) | — |
| [02 注释](/_projects/lotus/scripts/02-annotation/) | — |
| [03 端粒](/_projects/lotus/scripts/03-telomere/) | — |
| [04 rDNA](/_projects/lotus/scripts/04-rdna/) | — |
| [05 着丝粒](/_projects/lotus/scripts/05-centromere/) | — |
| [06 比较基因组](/_projects/lotus/scripts/06-genome-comparison/) | — |
| [07 转录组](/_projects/lotus/scripts/07-transcriptome/) | — |
| [08 单核转录组](/_projects/lotus/scripts/08-single-nucleus/) | — |


