---
title: '07 转录组'
sidebar:
  label: '07 转录组'
  order: -1
---


## 中文
**作用**：*Lotus japonicus* bulk RNA-seq——Gifu 与 MG20 组织图谱文库的预处理与定量、组织图谱与跨材料相关性、共表达分析（WGCNA、Mfuzz、Pearson 网络），以及突变体与公共数据定量。
**入口**：各步骤脚本见下表；阶段级 `run.sh` 调用 `Jaccard.py`。
**输入**：`00_prepare/align.conf` 与 Gifu、MG20 组织文库原始 FASTQ；`01_assembly`/`02_annotation` 的 GifuT2T FASTA 与 GTF；`07_coexpression_network/` 读取的基因、SNF 与同源列表；`01_tissue_atlas/` 与 `08_module_consistency/` 另需条件均值矩阵 `Gifu_all_samples_mean.tsv`、`MG20_all_samples_mean.tsv`。
**输出**：`00_prepare/all_samples.tpm.tsv` 供后续各步使用；各子目录的图表与表格见其 README。
**运行**：按 `00_prepare/` → `01_tissue_atlas/` → `05_wgcna/` → `06_mfuzz/` → `07_coexpression_network/` → `08_module_consistency/` 的顺序执行各目录脚本（各脚本按自身工作目录解析输入，不能写成一行）；路径占位符 ``{TOKEN}`` 见 `scripts/README.md`。
**工具**：`00_prepare/` 用 fastp、HISAT2、samtools、StringTie、GNU Awk；`01`、`02` 用 Python 3（pandas/numpy/scipy）与 R（tidyverse）；`05`–`08` 用 R 4.4.3（WGCNA、Mfuzz、Biobase、psych、data.table、igraph、ggraph、pheatmap、ggplot2）。
| 步骤 | 作用 |
|---|---|
| `00_prepare/` | fastp 质控、HISAT2 比对、StringTie 定量、TPM 矩阵 |
| `01_tissue_atlas/` | 组织图谱、SNF × 新基因相关、`spearman.R` 跨材料相关性 |
| `02_mutant_and_public_data/` | *Lotus* 结瘤突变体、大豆与 *Medicago truncatula* 定量 |
| `05_wgcna/` | 模块检测与合并、模块–性状相关、hub 基因 |
| `06_mfuzz/` | brown 与 lightcyan1 模块的时间序列软聚类 |
| `07_coexpression_network/` | Pearson 共表达网络 |
| `08_module_consistency/` | 模块基因的跨材料一致性 |
| `run.sh` + `Jaccard.py` | τ 组织特异性指数、Gifu/MG20 top-5000 基因集的 Jaccard 重叠 |

<details>
<summary>English</summary>
**Purpose**: Bulk RNA-seq of *Lotus japonicus* — preparation and quantification of the Gifu and MG20 tissue-atlas libraries, the atlas and its cross-accession correlation, co-expression analysis (WGCNA, Mfuzz, Pearson networks), and mutant and public-dataset quantification.
**Entry point**: the scripts listed in the table; the stage-level `run.sh` invokes `Jaccard.py`.
**Inputs**: `00_prepare/align.conf` and the raw FASTQ of the Gifu and MG20 tissue libraries; the GifuT2T FASTA/GTF from `01_assembly`/`02_annotation`; the gene, SNF and orthology lists read by `07_coexpression_network/`; `01_tissue_atlas/` and `08_module_consistency/` additionally need the condition-mean matrices `Gifu_all_samples_mean.tsv` and `MG20_all_samples_mean.tsv`.
**Outputs**: `00_prepare/all_samples.tpm.tsv`, consumed by every later step; the plots and tables of each sub-directory are listed in its README.
**Run**: run the scripts of `00_prepare/` → `01_tissue_atlas/` → `05_wgcna/` → `06_mfuzz/` → `07_coexpression_network/` → `08_module_consistency/` in that order (each script resolves its inputs relative to its own working directory, so the steps are not a one-liner); for the ``{TOKEN}`` path placeholders see `scripts/README.md`.
**Tools**: `00_prepare/` uses fastp, HISAT2, samtools, StringTie and GNU Awk; `01` and `02` use Python 3 (pandas/numpy/scipy) and R (tidyverse); `05`–`08` use R 4.4.3 (WGCNA, Mfuzz, Biobase, psych, data.table, igraph, ggraph, pheatmap, ggplot2).
| Step | Purpose |
|---|---|
| `00_prepare/` | fastp QC, HISAT2 alignment, StringTie quantification, TPM matrix |
| `01_tissue_atlas/` | tissue atlas, SNF × new-gene correlation, cross-accession correlation (`spearman.R`) |
| `02_mutant_and_public_data/` | *Lotus* nodulation mutants, soybean and *Medicago truncatula* quantification |
| `05_wgcna/` | module detection and merging, module–trait correlation, hub genes |
| `06_mfuzz/` | time-series soft clustering of the brown and lightcyan1 modules |
| `07_coexpression_network/` | Pearson co-expression network |
| `08_module_consistency/` | module gene consistency across accessions |
| `run.sh` + `Jaccard.py` | τ tissue-specificity index, Jaccard overlap of the Gifu/MG20 top-5000 gene sets |

</details>

### 本阶段步骤

| 步骤 | 内容 |
| --- | --- |
| [数据准备](/_projects/lotus/scripts/07-transcriptome/00-prepare/) | — |
| [组织表达图谱](/_projects/lotus/scripts/07-transcriptome/01-tissue-atlas/) | — |
| [突变体与公共数据](/_projects/lotus/scripts/07-transcriptome/02-mutant-and-public-data/) | — |
| [WGCNA 共表达模块](/_projects/lotus/scripts/07-transcriptome/05-wgcna/) | — |
| [Mfuzz 时序聚类](/_projects/lotus/scripts/07-transcriptome/06-mfuzz/) | — |
| [共表达网络](/_projects/lotus/scripts/07-transcriptome/07-coexpression-network/) | — |
| [模块一致性](/_projects/lotus/scripts/07-transcriptome/08-module-consistency/) | — |


