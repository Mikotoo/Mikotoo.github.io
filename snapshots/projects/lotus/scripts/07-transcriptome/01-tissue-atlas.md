---
title: '组织表达图谱'
sidebar:
  label: '组织表达图谱'
  order: 1
---


## 中文
**作用**：Gifu 组织图谱分析——把 SNF 家族基因与全图谱中相关性最高的 top-5“新”基因配对，并计算 Gifu 与 MG20 的跨材料表达相关性。
**入口**：`snf_new_corr.sh`、`spearman.R`；`getCount.sh` 只有计数任务的投递头，没有命令体。
**输入**：`snf_new_corr.sh` 需要 `expr.tsv`（TPM 矩阵，首列为基因 ID）与 `genes.tsv`（GeneID、Category，Category 含 `SNF` 或 `New`）；`spearman.R` 需要 `Gifu_all_samples_mean.tsv` 与 `MG20_all_samples_mean.tsv`。
**输出**：`<outprefix>.snf_new_corr/` 下的 `*.SNF_vs_New.all_pairs.tsv`、`*.SNF_top5_newgenes.tsv`、`*.SNF_top5_matrix.tsv`、`*.summary.txt`、`*.SNF_top5_heatmap.pdf|.png`、`*.SNF_top5_dotplot.pdf|.png`；`spearman.R` 写 `Gifu_MG20.sample_expression_spearman.tsv` 与 `Gifu_MG20.organ_expression_spearman.tsv` 及各自的 PDF。
**运行**：`bash snf_new_corr.sh expr.tsv genes.tsv outprefix`（须以本目录为工作目录）；`Rscript spearman.R`
**工具**：Python 3 + pandas/numpy/scipy（Pearson r 与 t 检验 p）；R + ggplot2/readr/dplyr/tidyr 与 tidyverse；集群投递头 `#CSUB -q c01 -n 64`。

<details>
<summary>English</summary>
**Purpose**: Tissue-atlas analysis for the Gifu accession — pairing SNF-family genes with the top-5 most correlated "new" genes across the atlas, and computing the cross-accession expression correlation between Gifu and MG20.
**Entry point**: `snf_new_corr.sh`, `spearman.R`; `getCount.sh` holds only the job header of the counting step and has no command body.
**Inputs**: `snf_new_corr.sh` needs `expr.tsv` (TPM matrix, first column = gene ID) and `genes.tsv` (GeneID and Category; Category containing `SNF` or `New`); `spearman.R` needs `Gifu_all_samples_mean.tsv` and `MG20_all_samples_mean.tsv`.
**Outputs**: in `<outprefix>.snf_new_corr/` — `*.SNF_vs_New.all_pairs.tsv`, `*.SNF_top5_newgenes.tsv`, `*.SNF_top5_matrix.tsv`, `*.summary.txt`, `*.SNF_top5_heatmap.pdf|.png` and `*.SNF_top5_dotplot.pdf|.png`; `spearman.R` writes `Gifu_MG20.sample_expression_spearman.tsv` and `Gifu_MG20.organ_expression_spearman.tsv`, each with a PDF.
**Run**: `bash snf_new_corr.sh expr.tsv genes.tsv outprefix` (run it with this directory as the working directory); `Rscript spearman.R`
**Tools**: Python 3 with pandas/numpy/scipy (Pearson r and t-test p); R with ggplot2/readr/dplyr/tidyr and tidyverse; cluster job header `#CSUB -q c01 -n 64`.

</details>


