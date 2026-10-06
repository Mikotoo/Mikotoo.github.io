---
title: '论文分析脚本（Rscript）'
sidebar:
  label: '论文分析脚本（Rscript）'
  order: -1
---


## 1. 项目概述

本目录保存 Lotus T2T 基因组论文中使用的主要分析脚本，覆盖以下四类分析：

1. Gifu、MG20 新旧基因组间的结构变异和 PAV 分析；
2. GifuT2T 与 GifuOld、MG20Old、MG20T2T 之间的同源基因对应；
3. 新增基因分类，包括 gap-filling、annotation-corrected 和 newly annotated genes；
4. bulk RNA-seq 与 single-nucleus RNA-seq 分析，包括 WGCNA、Mfuzz、共表达、细胞聚类、细胞类型注释和 scTenifoldKnk。

脚本主要对应论文 Fig. 2–4 及相关 Supplementary Data。当前目录不是整篇论文的完整复现包：Fig. 1 的基因组组装与质量评估流程、Fig. 5 的 LotusDB 构建流程未包含在内。

截至 2026-09-15，本目录共包含 51 个原有文件：23 个 Shell 脚本、16 个 Python 脚本、10 个 R 脚本、1 个配置文件和 1 个子目录说明文件。

## 2. 目录结构

```text
Rscript/
├── Structural variation/       # MUMmer4、SyRI 和 PAV 分析，对应 Fig. 2
├── Homologous genes/           # 新旧版本及品种间基因对应，对应 Fig. 2a
├── New genes finder/           # 新增基因识别和三类新基因划分，对应 Fig. 2b–d
├── Bulk-RNAseq analysis/       # 表达相关性、WGCNA、Mfuzz 和共表达，对应 Fig. 3
└── Singlecell-RNAseq analysis/ # 单核 RNA-seq 和 scTenifoldKnk，对应 Fig. 4
```

所有脚本均按服务器目录设计，默认项目根目录为：

```text
/share/org/YZWL/yzwl_hanxuu/Project/Lotus
```

脚本中的绝对路径、线程数、队列名称和 Conda 环境需要根据实际运行环境检查后再执行。

## 3. 脚本状态说明

| 状态 | 含义 |
|---|---|
| 主流程 | 直接产生论文图表或补充数据所依赖的结果，应长期保留 |
| 上游流程 | 生成主流程所需的比对、表达矩阵或对象，应保留 |
| 辅助脚本 | 用于格式转换、质控、注释或绘图，通常不单独运行 |
| 运行前检查 | 脚本有硬编码路径、非幂等输出或依赖已有中间文件，不能直接无检查重跑 |

## 4. Structural variation

### 4.1 主要脚本

| 脚本 | 比较组合 | 主要用途 |
|---|---|---|
| `GifuT2T_vs_GifuOld_variation_finder.sh` | GifuOld vs GifuT2T | MUMmer4 比对、SyRI 变异识别及 plotsr 绘图 |
| `MG20T2T_MG20Old_variation_finder.sh` | MG20Old vs MG20T2T | MUMmer4 比对、SyRI 变异识别及 plotsr 绘图 |
| `GifuT2T_vs_MG20T2T_variation_finder.sh` | MG20T2T vs GifuT2T | 两个 T2T 品种间 SNP、SV 和重排识别 |
| `*_get_PAV_from_coords.sh` | 对应上述三个组合 | 根据比对区间的补集计算未比对区段和 PAV |

### 4.2 推荐运行顺序

每个比较组合应在独立工作目录运行，避免不同组合产生的 `syri.out`、BED 和临时文件互相覆盖。

```bash
bash <comparison>_variation_finder.sh
bash <comparison>_get_PAV_from_coords.sh
```

### 4.3 运行前检查

- 依赖：MUMmer4、SyRI、plotsr、samtools、bedtools、seqkit、awk。
- `variation_finder.sh` 内部分类型 BED 使用 `>>` 追加。重跑前必须删除旧 BED，或把脚本改成先清空输出，否则会重复计数。
- `*.PAV.summary` 同样采用追加方式，重跑前应清空。
- BED 是 0-based、half-open 格式时，区间长度应使用 `$3-$2`。当前 GifuOld 和 GifuT2T-vs-MG20T2T 的 PAV 汇总仍含 `$3-$2+1`，正式重跑前应统一确认。
- FASTA 中染色体名称必须与 coords、SyRI 输出和 genome 文件一致。

## 5. Homologous genes

### 5.1 文件说明

- `blast.sh`：依次完成 GifuT2T 对 GifuOld、MG20Old 和 MG20T2T 的 BLASTP，并合并为 `All_species_merged_by_GifuT2T.tsv`。
- `blast_merge.py`：聚合 BLAST HSP，结合蛋白长度计算 identity 和覆盖率。

### 5.2 运行方法

在同一工作目录准备以下文件：

```text
Lotus_GifuOld.pep.fa
Lotus_GifuT2T.pep.fa
Lotus_MG20Old.pep.fa
Lotus_MG20T2T.pep.fa
blast.sh
blast_merge.py
```

随后运行：

```bash
bash blast.sh
```

### 5.3 运行前检查

- 依赖：BLAST+、seqkit、Python 3、pandas、awk。
- `ALL_GENE_LIST` 当前为服务器绝对路径，应确认文件内容是完整的 GifuT2T 基因列表。
- `blast_merge.py` 当前按 HSP 长度求和；若同一 query–subject 的 HSP 重叠，覆盖率可能被高估。严格复现时应先合并 query 和 subject 上的重叠区间。
- 论文方法中的 JCVI 共线性验证不在当前脚本中；如该步骤用于最终基因对应结果，应另行补充脚本和参数。

## 6. New genes finder

### 6.1 文件说明

- `new_genes_finder.sh`：单个基因组版本比较的主体流程。
- `new_genes_find.sh`：依次启动 Gifu 和 MG20 分析的示例入口。

主体流程首先通过 BLASTP 找出 T2T 注释中没有旧版本蛋白同源关系的基因，然后将这些基因划分为：

1. `new_genes_Gap`：位于旧基因组未覆盖区域；
2. `genes_anno_correct`：旧基因组中存在对应序列，但旧注释需要修正；
3. `new_genes_anno`：排除前两类后的新注释基因。

### 6.2 主要输出

```text
all_new_genes.list
new_genes_Gap/new_genes_Gap.list
genes_anno_correct/genes_anno_correct.list
new_genes_anno/new_genes_anno.list
```

### 6.3 推荐运行方式

Gifu 和 MG20 必须在不同工作目录运行，因为主体脚本会生成同名的 `blastp/`、`blastn/`、`all_new_genes.list` 等文件。

```bash
mkdir -p Gifu && cd Gifu
bash ../new_genes_finder.sh <T2T.pep> <T2T.fa> <T2T.gff3> \
  <Old.pep> <Old.fa> <Old.gff3> <Old_region.uncovered.bed>

mkdir -p ../MG20 && cd ../MG20
bash ../new_genes_finder.sh <T2T.pep> <T2T.fa> <T2T.gff3> \
  <Old.pep> <Old.fa> <Old.gff3> <Old_region.uncovered.bed>
```

## 7. Bulk-RNAseq analysis

### 7.1 主要脚本

| 脚本 | 状态 | 用途 |
|---|---|---|
| `fastp.sh` | 上游流程 | 双端 RNA-seq 数据质控和过滤 |
| `align.conf` | 配置文件 | 样本名及双端 FASTQ 路径 |
| `align.sh` | 上游流程 | HISAT2 比对、唯一比对筛选、StringTie 定量和 TPM 矩阵合并 |
| `spearman.R` | 主流程 | 样本或组织表达 Spearman 相关性，对应 Fig. 3 |
| `WGCNA.R` | 主流程 | bulk RNA-seq 共表达模块分析 |
| `Mfuzz.R` | 主流程 | 时间序列表达模式聚类，最终论文采用 4 个 clusters |
| `co-exp_psych.R` | 主流程 | Pearson 相关和显著性过滤，生成共表达网络 |
| `module_gene_expression_correlation.R` | 主流程 | Gifu/MG20 模块基因表达一致性及审稿补充分析 |

### 7.2 推荐运行顺序

```bash
bash fastp.sh
bash align.sh
Rscript spearman.R
Rscript WGCNA.R
Rscript Mfuzz.R
Rscript co-exp_psych.R
Rscript module_gene_expression_correlation.R
```

后五个 R 脚本并非完全线性依赖，应根据各脚本中设定的输入矩阵分别运行。

## 8. Singlecell-RNAseq analysis

### 8.1 主流程顺序

| 顺序 | 入口脚本 | 作用 |
|---:|---|---|
| 00 | `00_mkref.sh` | 使用 GifuT2T FASTA/GTF 构建 Cell Ranger reference |
| 00 | `00_cellranger.sh` | 对各样本执行 `cellranger count` |
| 01 | `01_annaDataQC.sh` | 读取 Cell Ranger 矩阵，执行 QC、doublet 识别和过滤 |
| 02 | `02_integrate.sh` | 合并样本并使用 scVI 进行整合 |
| 03 | `03_dimenRedu_adjusted.sh` | UMAP、Leiden 聚类及不同高变基因数比较 |
| 04 | `04_clusterSplit.sh` | 聚类拆分、样本组成和 UMAP 绘图 |
| 05 | `05_cluster_DEGs_finder.sh` | Seurat/RDS 转换、cluster marker 和富集分析 |
| 06 | `06_1.marker_gene_finder.sh` | 与其他植物 marker 蛋白比对，辅助细胞类型判断 |
| 06 | `06_1.Individual_clusters_marker_gene_finder.sh` | 针对指定 cluster 提取 DEG 与注释 |
| 06 | `06_2.markerAnno.sh` | 根据 marker 结果绘制和整理细胞类型注释 |
| 06 | `06_3.cell_type_exp.sh` | 按细胞类型进行 DEG 分析 |
| 07 | `07_batch_multiple_scTenifoldKnk.sh` | 为候选基因生成批量 scTenifoldKnk 作业 |

### 8.2 Script 子目录

- `01_annaDataQC.py`：样本读取、质量控制、doublet 识别和合并。
- `02_integrate.py`：scVI 模型训练和样本整合。
- `03_1-dimenRedu_adjusted.py`：UMAP、Leiden、多参数聚类与聚类组成统计。
- `03_1-cluster_statitcs.py`：独立的 cluster/样本细胞数统计图。
- `03_2-fliter_to_all_gene_h5ad.py`：把聚类结果映射回全基因 AnnData 对象。
- `04_1-*`：聚类拆分、样本分布、UMAP 及细胞组成绘图。
- `05_1-Seurat_Exp_h5ad_to_rds.py` 和 `05_2-Exp_h5ad_to_rds_cluster.R`：Python/Seurat 对象转换及 marker 分析。
- `05_3-marker_cluster_DEGs_own_order.py`：按指定顺序绘制 cluster marker。
- `05_3-marker_knowngene_lotus_own_order.py`：绘制已知 Lotus marker。
- `05_3-markerAnno_Ps_v5_own.py`：整理和绘制细胞类型注释。
- `05_4-cluster8_DEG_new_UMAP.py`：early signaling cluster 候选基因 UMAP。
- `06_1-Seurat_Exp_h5ad_to_rds.py`、`06_2-Exp_h5ad_to_rds_cell_type.R`、`06_3-agg_to_fpkm.R`、`06_4-DEGs_cell_type.R`：细胞类型层面的表达聚合和 DEG 分析。
- `07_batch_multiple_scTenifoldKnk.R`：执行单个候选基因的多次 scTenifoldKnk 扰动分析。
- `makesure_02dpi_samples_all_have_Early_signal_cell.py`：验证两个 2 dpi 重复中均存在 early signaling cells。
- `plot_02dpi1_02dpi2_gene_expression_umap.py`：分别绘制两个 2 dpi 重复的候选基因表达 UMAP。

### 8.3 scTenifoldKnk 批量运行

输入基因列表至少包含两列：gene 和 group。第二个命令行参数为要分析的细胞类型。

```bash
bash 07_batch_multiple_scTenifoldKnk.sh candidate_gene.list "Early signaling cell"
```

脚本生成：

```text
06_batch_multiple_scTenifoldKnk/<group>/<gene>/<gene>_scTenifoldKnk.sh
scTenifoldKnk_run.list
```

随后按服务器集群提交规则提交 `scTenifoldKnk_run.list` 中的作业。

## 9. 软件环境与实测版本

以下版本于 2026-09-15 在服务器 `ln01n001` 上实测。Conda 版本为 25.7.0。检查时分别激活 `lotus`、`singlecell` 和 `r_wgcna` 环境，并同时使用命令行 `--version`、Python package metadata 或 R `packageVersion()` 获取版本。

### 9.1 `lotus` 环境

环境路径：

```text
/share/org/YZWL/yzwl_hanxuu/anaconda3/envs/lotus
```

该环境主要用于结构变异、同源基因、新基因和 bulk RNA-seq 上游分析。

| 软件或包 | 实测版本 | 主要用途 |
|---|---:|---|
| Bash | 4.4.20 | Shell 主流程 |
| GNU Awk | 5.1.0 | 表格筛选和统计 |
| Python | 3.10.14 | `blast_merge.py` 等辅助脚本 |
| BLAST+ (`blastp`/`makeblastdb`) | 2.16.0+ | 蛋白和核酸同源搜索 |
| bedtools | 2.31.1 | 区间交集、补集和 PAV 分析 |
| fastp | 1.1.0 | bulk RNA-seq FASTQ 质控 |
| HISAT2 | 2.2.1 | bulk RNA-seq 比对 |
| MUMmer4 (`nucmer`/`dnadiff`) | 4.0.1 | 基因组比对和差异统计 |
| samtools | 1.21 | BAM/FASTA 索引和处理 |
| seqkit | 2.10.1 | FASTA/FASTQ 筛选和长度统计 |
| StringTie | 2.1.7 | 转录本定量和 TPM 输出 |
| SyRI | 1.7.1 | 共线区和结构变异识别 |
| plotsr | 1.1.1 | SyRI 结果绘图 |
| NumPy | 1.26.4 | Python 数值计算 |
| pandas | 2.3.3 | BLAST 表格处理 |
| SciPy | 1.14.1 | Python 科学计算依赖 |

补充说明：该环境同时安装了旧版 `mummer 3.23` 和 `mummer4 4.0.1`，但当前 `nucmer --version` 实际解析为 4.0.1。运行记录中应保存 `command -v nucmer` 和 `nucmer --version`，避免 PATH 变化后误用旧版。

### 9.2 `singlecell` 环境

环境路径：

```text
/share/org/YZWL/yzwl_hanxuu/anaconda3/envs/singlecell
```

#### Python 软件包

| 软件包 | 实测版本 |
|---|---:|
| Python | 3.9.19 |
| anndata | 0.10.8 |
| scanpy | 1.10.1 |
| scvi-tools | 1.1.2 |
| rpy2 | 3.5.11 |
| NumPy | 1.24.1 |
| pandas | 2.2.2 |
| SciPy | 1.13.1 |
| matplotlib | 3.9.1 |
| seaborn | 0.13.2 |
| leidenalg | 0.10.2 |
| scikit-learn | 1.1.3 |
| PyTorch | 2.4.0 |
| scikit-misc | 未检测到 |

#### R 及主要 R 包

| 软件包 | 实测版本 |
|---|---:|
| R | 4.4.1 |
| Seurat | 4.3.0 |
| scTenifoldKnk | 1.0.2 |
| clustree | 0.5.1 |
| DESeq2 | 1.46.0 |
| Biobase | 2.66.0 |
| GenomicFeatures | 1.58.0 |
| ggplot2 | 4.0.1 |
| data.table | 1.17.8 |
| dplyr | 1.1.4 |
| gghalves | 0.1.4 |
| ggraph | 2.2.2 |
| ggrepel | 0.9.6 |
| glue | 1.8.0 |
| igraph | 2.2.1 |
| patchwork | 1.3.2 |
| pheatmap | 1.0.13 |
| RColorBrewer | 1.1.3 |
| remotes | 2.5.0 |
| scales | 1.4.0 |
| tidygraph | 1.3.0 |
| tidyverse | 2.0.0 |

### 9.3 `r_wgcna` 环境

环境路径：

```text
/share/org/YZWL/yzwl_hanxuu/anaconda3/envs/r_wgcna
```

该环境主要用于 bulk RNA-seq 的 WGCNA、Mfuzz、相关性和网络绘图。

| 软件包 | 实测版本 |
|---|---:|
| R | 4.4.3 |
| WGCNA | 1.73 |
| Mfuzz | 2.66.0 |
| psych | 2.5.6 |
| ggplot2 | 4.0.1 |
| tidyverse | 2.0.0 |
| dplyr | 1.1.4 |
| data.table | 1.17.8 |
| igraph | 2.1.4 |
| ggraph | 2.2.1 |
| tidygraph | 1.3.0 |
| pheatmap | 1.0.13 |
| RColorBrewer | 1.1.3 |
| Biobase | 2.66.0 |
| patchwork | 1.3.2 |
| ggrepel | 0.9.6 |
| glue | 1.8.0 |
| scales | 1.4.0 |

### 9.4 `tools` 目录中的独立软件和源码

部分软件不属于 `lotus`、`singlecell` 或 `r_wgcna` 的 Conda 包，而是独立保存在：

```text
/share/org/YZWL/yzwl_hanxuu/tools
```

| 软件 | 实测版本或状态 | 实际路径 | 使用说明 |
|---|---:|---|---|
| Cell Ranger | 9.0.1 | `/share/org/YZWL/yzwl_hanxuu/tools/cellranger-9.0.1` | 已安装，但默认未加入 `lotus` 和 `singlecell` 的 PATH；用于 `00_mkref.sh` 和 `00_cellranger.sh` |
| JCVI | 1.5.7 | 源码：`/share/org/YZWL/yzwl_hanxuu/tools/jcvi-main`；运行环境：`~/anaconda3/envs/jcvi_env` | `jcvi_env` 使用 Python 3.12.12；用于论文中基因对应关系的共线性验证 |
| scTenifoldKnk 源码 | 1.0.2 | `/share/org/YZWL/yzwl_hanxuu/tools/scTenifoldKnk-master` | 源码 `DESCRIPTION` 版本与 `singlecell` 环境中安装的 R 包 1.0.2 一致 |
| `Scripts.scDblFinder` | 未声明独立版本号 | 由 `singlecell` 环境中的自定义 Python `Scripts` 模块提供 | 建议保存模块源码、实际文件路径和 Git commit；不要将它与 Bioconductor 的 R 包 `scDblFinder` 混为一谈 |

### 本阶段步骤

| 步骤 | 内容 |
| --- | --- |
| [Singlecell-RNAseq analysis](/_projects/rscript/singlecell-rnaseq-analysis/) | — |


