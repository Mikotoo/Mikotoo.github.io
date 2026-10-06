---
title: 'Mfuzz 时序聚类'
sidebar:
  label: 'Mfuzz 时序聚类'
  order: 6
---


## 中文
**作用**：把 WGCNA brown 与 lightcyan1 模块的基因沿结瘤发育序列（root_uni → root_hpi48 → nodule_dpi10 → dpi21 → dpi45）做软聚类，再把同一批簇投影到根序列上。
**入口**：`Mfuzz.R`
**输入**：工作目录中的 `Gifu_all_samples.tpm.tsv`、`MG20_all_samples.tpm.tsv`，以及 `./result/06.modules_expr/module_brown.csv` 与 `module_lightcyan1.csv`。
**输出**：`result/12.Mfuzz_out/` 下的 `Dmin_plot.pdf`、`Mfuzz_stability.pdf`、`Mfuzz_test.pdf`、`cluster1.genes.list`…`cluster4.genes.list`、`Mfuzz_membership.tsv`、`Mfuzz_plot.pdf`、`Mfuzz_plot2.pdf`、`Mfuzz_root_projection_plot.pdf`、`Mfuzz_root_projection_plot2.pdf`。
**运行**：`Rscript Mfuzz.R`（无参数；工作目录须同时持有两张 TPM 矩阵与 `WGCNA.R` 生成的 `result/` 树）。
**工具**：R 4.4.3 + Mfuzz（`mestimate`、`Dmin`、`mfuzz`、`mfuzz.plot`、`mfuzz.plot2`）、Biobase、tidyverse；参数：阶段均值 TPM 取 `log1p`、`filter.NA(thres = 0.25)`、`fill.NA(mode = "mean")`、`filter.std(min.std = 0)`、`Dmin(crange = seq(4, 20), repeats = 3)`、最终 `c = 4`（`set.seed(1234)`）。

<details>
<summary>English</summary>
**Purpose**: Soft-cluster the genes of the WGCNA brown and lightcyan1 modules along the nodule developmental series (root_uni → root_hpi48 → nodule_dpi10 → dpi21 → dpi45) and project the same clusters onto the root series.
**Entry point**: `Mfuzz.R`
**Inputs**: `Gifu_all_samples.tpm.tsv` and `MG20_all_samples.tpm.tsv` in the working directory, plus `./result/06.modules_expr/module_brown.csv` and `module_lightcyan1.csv`.
**Outputs**: under `result/12.Mfuzz_out/` — `Dmin_plot.pdf`, `Mfuzz_stability.pdf`, `Mfuzz_test.pdf`, `cluster1.genes.list` … `cluster4.genes.list`, `Mfuzz_membership.tsv`, `Mfuzz_plot.pdf`, `Mfuzz_plot2.pdf`, `Mfuzz_root_projection_plot.pdf` and `Mfuzz_root_projection_plot2.pdf`.
**Run**: `Rscript Mfuzz.R` (no arguments; the working directory must hold both TPM matrices and the `result/` tree produced by `WGCNA.R`).
**Tools**: R 4.4.3 with Mfuzz (`mestimate`, `Dmin`, `mfuzz`, `mfuzz.plot`, `mfuzz.plot2`), Biobase and tidyverse; settings as coded: `log1p` on stage-mean TPM, `filter.NA(thres = 0.25)`, `fill.NA(mode = "mean")`, `filter.std(min.std = 0)`, `Dmin(crange = seq(4, 20), repeats = 3)`, and the final clustering `c = 4` with `set.seed(1234)`.

</details>


