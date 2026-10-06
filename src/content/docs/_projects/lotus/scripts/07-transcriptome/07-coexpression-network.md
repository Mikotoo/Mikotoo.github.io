---
title: '共表达网络'
sidebar:
  label: '共表达网络'
  order: 7
---


## 中文
**作用**：把 SNF 家族基因与候选基因列表在两种材料的根、接种根与结瘤样本中做相关，保留通过 Bonferroni 校正 p 值与相关系数阈值的基因对并写出网络。
**入口**：`co-exp_psych.R`
**输入**：参数 `go_id_list.txt`（首列为候选基因 ID）与 `out_dir`（输出目录）；`../../../Gifu_all_samples.tpm.tsv`、`../../../MG20_all_samples.tpm.tsv`；工作目录中的 `all_new_genes.list`、`SNF_gene.list` 与同源表 `All_species_merged_by_GifuT2T.tsv`。
**输出**：`<out_dir>/SNF_vs_ALL/` 与 `<out_dir>/SNF_vs_New/`，各含 `edges_bonferroni.txt`、`nodes.txt`、`network_igraph.rds`；`<out_dir>/network_gene_composition.txt` 给出各节点类别的计数。
**运行**：`Rscript co-exp_psych.R go_id_list.txt output_dir`（须从两张 TPM 矩阵之上三层的目录运行）。
**工具**：R 4.4.3 + psych（`corr.test(method = "pearson", adjust = "none")`）、igraph、data.table 与 tidyverse；参数：样本列匹配 `root_hpi48`/`root_uni`/`nodule`，最小 TPM > 0.5，`alpha = 0.05`、`p_cutoff = alpha / nrow(edges_all)`，保留 `p < p_cutoff` 且 `cor > 0.7` 的边。

<details>
<summary>English</summary>
**Purpose**: Correlate SNF-family genes with a candidate gene list across the root, inoculated-root and nodule samples of both accessions, keep pairs passing a Bonferroni-corrected p-value and a correlation cutoff, and write the resulting networks.
**Entry point**: `co-exp_psych.R`
**Inputs**: the arguments `go_id_list.txt` (first column = candidate gene IDs) and `out_dir` (output directory); `../../../Gifu_all_samples.tpm.tsv` and `../../../MG20_all_samples.tpm.tsv`; `all_new_genes.list`, `SNF_gene.list` and the orthology table `All_species_merged_by_GifuT2T.tsv` in the working directory.
**Outputs**: `<out_dir>/SNF_vs_ALL/` and `<out_dir>/SNF_vs_New/`, each with `edges_bonferroni.txt`, `nodes.txt` and `network_igraph.rds`; `<out_dir>/network_gene_composition.txt` with the counts per node class.
**Run**: `Rscript co-exp_psych.R go_id_list.txt output_dir` (run from a directory three levels below the two TPM matrices).
**Tools**: R 4.4.3 with psych (`corr.test(method = "pearson", adjust = "none")`), igraph, data.table and tidyverse; settings as coded: sample columns matching `root_hpi48`/`root_uni`/`nodule`, minimum TPM > 0.5, `alpha = 0.05` with `p_cutoff = alpha / nrow(edges_all)`, edges kept when `p < p_cutoff` and `cor > 0.7`.

</details>


