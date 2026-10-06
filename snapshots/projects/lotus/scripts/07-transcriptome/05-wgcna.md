---
title: 'WGCNA 共表达模块'
sidebar:
  label: 'WGCNA 共表达模块'
  order: 5
---


## 中文
**作用**：由 Gifu 与 MG20 组织图谱 TPM 矩阵构建加权基因共表达网络，检测并合并模块，把模块特征基因与组织性状关联，并提取 hub 基因。
**入口**：`WGCNA.R`
**输入**：上一级目录的 `Gifu_all_samples.tpm.tsv`、`MG20_all_samples.tpm.tsv`（`gene_id` 为行名、每样本一列）；样本分组从列名解析。
**输出**：`result/` 下的样本聚类、软阈值、基因聚类、模块检测与合并 PDF；`06.modules_expr/module_<color>.csv`（每模块的基因 × 样本矩阵）；`07.module_trait_correlation.csv`、`08.module_trait_pvalue.csv`、`07.trait_*.csv`、`09.module_trait_heatmap.pdf`；`10.hub_gene_results/`、`11.cytoscape/`，以及每模块一张 `12.gene_trait_combined_network_<module>.pdf`。
**运行**：`Rscript WGCNA.R`（无参数；在打算生成 `result/` 的目录内运行，两张 TPM 矩阵放在该目录的上一级）。
**工具**：R 4.4.3 + WGCNA（`adjacency`、`TOMsimilarity`、`cutreeDynamic`、`mergeCloseModules`、`exportNetworkToCytoscape`）、tidyverse/dplyr、tidygraph/igraph/ggraph、pheatmap；参数：软阈值 power 14、`minModuleSize` 30、`deepSplit` 2、合并 `cutHeight` 0.3、hub 阈值 MM ≥ 0.7 且 GS ≥ 0.5、基因过滤为至少一个条件组均值 TPM > 0.5 后取 MAD 前 75%。

<details>
<summary>English</summary>
**Purpose**: Build a weighted gene co-expression network from the Gifu and MG20 tissue-atlas TPM matrices, detect and merge modules, relate module eigengenes to tissue traits and extract hub genes.
**Entry point**: `WGCNA.R`
**Inputs**: `Gifu_all_samples.tpm.tsv` and `MG20_all_samples.tpm.tsv` one directory above the working directory (`gene_id` as row names, one column per sample); sample groups are parsed from the column names.
**Outputs**: under `result/`, the sample-clustering, soft-threshold, gene-clustering, module-detection and module-merge PDFs; `06.modules_expr/module_<color>.csv` (per-module gene × sample matrix); `07.module_trait_correlation.csv`, `08.module_trait_pvalue.csv`, `07.trait_*.csv`, `09.module_trait_heatmap.pdf`; `10.hub_gene_results/`, `11.cytoscape/`, and one `12.gene_trait_combined_network_<module>.pdf` per module.
**Run**: `Rscript WGCNA.R` (no arguments; run it from the directory intended to hold `result/`, with the two TPM matrices in that directory's parent).
**Tools**: R 4.4.3 with WGCNA (`adjacency`, `TOMsimilarity`, `cutreeDynamic`, `mergeCloseModules`, `exportNetworkToCytoscape`), tidyverse/dplyr, tidygraph/igraph/ggraph, pheatmap; settings as coded: soft power 14, `minModuleSize` 30, `deepSplit` 2, merge `cutHeight` 0.3, hub thresholds MM ≥ 0.7 and GS ≥ 0.5, gene filter = group mean TPM > 0.5 in at least one condition, then the top 75% by MAD.

</details>


