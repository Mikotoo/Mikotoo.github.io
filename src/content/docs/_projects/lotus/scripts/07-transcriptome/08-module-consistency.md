---
title: '模块一致性'
sidebar:
  label: '模块一致性'
  order: 8
---


## 中文
**作用**：在 12 个配对条件上比较 WGCNA brown 与 lightcyan1 模块基因在 Gifu 与 MG20 之间的逐基因表达模式，以 Spearman 为主、Pearson 为辅，写出逐基因表与相关分布图。
**入口**：`module_gene_expression_correlation.R`
**输入**：`Gifu_all_samples_mean.tsv`、`MG20_all_samples_mean.tsv`（`gene_id` 列加 12 个条件列，脚本剥离 `Gifu_`/`MG20_` 前缀）；`result/06.modules_expr/module_brown.csv` 与 `module_lightcyan1.csv`。
**输出**：`result/13.module_gene_expression_correlation/` 下的 `brown_lightcyan1_gene_expression_correlations.tsv` 与 `.filtered.tsv`、`correlation_threshold_summary.tsv`、`spearman_correlation_distribution.pdf|.png`、`pearson_correlation_distribution.pdf|.png` 与 `README.txt`。
**运行**：`Rscript module_gene_expression_correlation.R`（无参数；脚本按自身目录解析路径，须与两张均值矩阵及 `result/` 树并列并从该目录运行）。
**工具**：R 4.4.3 + ggplot2/scales（`ggsave(device = cairo_pdf)`）与 base R（`cor.test(method = "spearman"/"pearson")`）；参数：相关基于 12 个条件下的 log2(均值 TPM + 1)，高相关阈值 R ≥ 0.8，基因需在至少 2 个条件下均值 TPM ≥ 0.5 且两向量 SD > 0。

<details>
<summary>English</summary>
**Purpose**: Compare the per-gene expression pattern of the WGCNA brown and lightcyan1 module genes between Gifu and MG20 across the 12 matched conditions, using Spearman as the primary statistic and Pearson as an auxiliary one, and write the per-gene tables and correlation distribution plots.
**Entry point**: `module_gene_expression_correlation.R`
**Inputs**: `Gifu_all_samples_mean.tsv` and `MG20_all_samples_mean.tsv` (a `gene_id` column plus the 12 condition columns; the `Gifu_`/`MG20_` prefix is stripped in the script); `result/06.modules_expr/module_brown.csv` and `module_lightcyan1.csv`.
**Outputs**: under `result/13.module_gene_expression_correlation/` — `brown_lightcyan1_gene_expression_correlations.tsv` and `.filtered.tsv`, `correlation_threshold_summary.tsv`, `spearman_correlation_distribution.pdf|.png`, `pearson_correlation_distribution.pdf|.png` and `README.txt`.
**Run**: `Rscript module_gene_expression_correlation.R` (no arguments; all paths are resolved relative to the script's own directory, so the script must sit beside the two mean matrices and the `result/` tree and be run from that directory).
**Tools**: R 4.4.3 with ggplot2/scales (`ggsave(device = cairo_pdf)`) and base R (`cor.test(method = "spearman"/"pearson")`); settings as coded: correlations on log2(mean TPM + 1) over the 12 conditions, high-correlation threshold R ≥ 0.8, and a gene is analysed when it has mean TPM ≥ 0.5 in at least 2 conditions and both vectors have SD > 0.

</details>


