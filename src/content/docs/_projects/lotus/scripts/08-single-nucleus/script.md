---
title: '单核分析脚本'
sidebar:
  label: '单核分析脚本'
  order: 999
---


## 中文
**作用**：`08_single_nucleus` 流程的 Python 与 R 实现脚本，由上一级目录的 7 个 `*.sh` 按相对路径调用。
**入口**：无独立入口，由上一级入口脚本驱动（`python ./Script/<name>.py`、`Rscript ./Script/<name>.R`）。
**输入**：`../02_processData/` 下的 `*_filtered.h5ad`、`data_clusterd_all_gene_{n}.h5ad`、`zz_00_pea.rds`、`All_species_merged_by_GifuT2T.SNF_Symbol.tsv`；`../04_Anno/` 下的 marker CSV 与 Gifu T2T GTF。
**输出**：`../02_processData/` 的 h5ad/rds/CSV、`../03_figures/` 的图、`../04_Anno/` 的 Seurat 对象与 marker 表，以及每个基因目录下的 `01_runs/`、`02_summaries/`、`03_plots/`。
**运行**：`bash ../run_singlecell.sh`、`bash ../05_cluster_DEGs_finder.sh` 等；无入口脚本的成员手动运行，目录名 `Script` 不可更改。
**工具**：scanpy、anndata、scvi-tools、leidenalg、rpy2（Python 3.9.19）；Seurat、scTenifoldKnk、DESeq2、clustree（R 4.4.1）。

<details>
<summary>English</summary>
**Purpose**: The Python and R implementation scripts of the `08_single_nucleus` workflow, called by relative path from the seven `*.sh` entry scripts in the parent directory.
**Entry point**: none of their own; they are driven by the parent entry scripts (`python ./Script/<name>.py`, `Rscript ./Script/<name>.R`).
**Inputs**: `*_filtered.h5ad`, `data_clusterd_all_gene_{n}.h5ad`, `zz_00_pea.rds` and `All_species_merged_by_GifuT2T.SNF_Symbol.tsv` under `../02_processData/`; the marker CSVs and the Gifu T2T GTF under `../04_Anno/`.
**Outputs**: h5ad/rds/CSV files in `../02_processData/`, plot files in `../03_figures/`, Seurat objects and marker tables in `../04_Anno/`, and `01_runs/`, `02_summaries/`, `03_plots/` inside each gene directory.
**Run**: `bash ../run_singlecell.sh`, `bash ../05_cluster_DEGs_finder.sh`, and so on; the members without an entry script are run by hand, and the `Script` directory name must not change.
**Tools**: scanpy, anndata, scvi-tools, leidenalg, rpy2 (Python 3.9.19); Seurat, scTenifoldKnk, DESeq2, clustree (R 4.4.1).

</details>


