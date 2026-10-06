---
title: '08 单核转录组'
sidebar:
  label: '08 单核转录组'
  order: -1
---


## 中文
**作用**：把 Gifu T2T 基因组上的单核 RNA-seq（根的两个阶段 Dpi00 与 Dpi02，各两个重复）从 Cell Ranger 定量做到聚类与细胞类型注释、簇/细胞类型 DEG 与 scTenifoldKnk 虚拟敲除。
**入口**：下表 7 个 shell 脚本；`Script/*.py` 中需手动运行的为 `03_1-cluster_statitcs.py`、`03_2-fliter_to_all_gene_h5ad.py`、`05_3-marker_cluster_DEGs_own_order.py` 与 `makesure_02dpi_samples_all_have_Early_signal_cell.py`。
**输入**：`sample.list`（`ProjectID CNName RunID SampleName`）、四个 `*_filtered.h5ad`、Gifu T2T FASTA/GTF 与 `.pep.fa`、`{Gmax,Gm2.1,GmZH,At}_marker_genes.pep.fa`、eggNOG 注释；后续步骤按当前目录读 `02_processData/` 与 `04_Anno/`，scTenifoldKnk 另需 `zz_00_pea.rds` 与 `All_species_merged_by_GifuT2T.SNF_Symbol.tsv`。
**输出**：`00_ref/`、`01_cellranger_out/<SampleName>/`、`02_processData/`（h5ad、`zz_*`、`zz_00_pea.rds`、`00_02_pseudobulk_fpkm.csv`）、`03_figures/`（降维图、`clustree_plot_*.pdf`、火山图）、`04_Anno/`（Seurat rds、marker CSV、`blast_*/`、富集结果）、`06_batch_multiple_scTenifoldKnk/<group>/<gene>/`。
**运行**：在本阶段根目录（含 `02_processData/`、`03_figures/`、`04_Anno/`）按表中顺序运行，编号脚本携带 `#CSUB` 作业指令、交由调度器提交；``{TOKEN}`` 路径占位符见 `scripts/README.md`；注意 `03_1-dimenRedu_adjusted.py` 写出的文件名与 `03_1-cluster_statitcs.py`、`03_2-fliter_to_all_gene_h5ad.py` 读取的文件名不一致，需先统一。
**工具**：Cell Ranger 9.0.1；Python 3.9.19 + scanpy/anndata/scvi-tools/leidenalg/rpy2；R 4.4.1 + Seurat/scTenifoldKnk/DESeq2/clustree；BLAST+ 与 seqkit；参数：`n_top_genes` 4000（扫描 1000–5000）、Leiden 0.4（扫描 0.2–1.2）、每基因 100 次 scTenifoldKnk 重复（K ∈ {10,30,50,70,100}）。
| 脚本 | 用途 |
|---|---|
| `00_mkref.sh` | `cellranger mkref --genome=Gifu`（Gifu T2T FASTA + GTF）→ `./00_ref` |
| `00_cellranger.sh` | 对 `sample.list` 中每个样本执行一次 `cellranger count` → `01_cellranger_out/<SampleName>` |
| `run_singlecell.sh` | 一个驱动依次跑 `Script/01_annaDataQC.py`、`02_integrate.py`、`03_1-dimenRedu_adjusted.py`、`05_3-markerAnno_Ps_v5_own.py`、`06_4-DEGs_cell_type.R` |
| `05_cluster_DEGs_finder.sh` | 先 `05_1-…py` 再 `05_2-…R`：Seurat 对象与簇 marker，随后按簇做富集 |
| `06_1.marker_gene_finder.sh` | DEG → 肽 FASTA → 对 Gmax/Gm2.1/GmZH/At marker 做 `blastp` 并合并结果 |
| `06_1.Individual_clusters_marker_gene_finder.sh` | 取各簇 top-20 DEG 序列并抽 eggNOG 注释（启用簇 `8`） |
| `07_batch_multiple_scTenifoldKnk.sh` | 每个候选基因生成一个 scTenifoldKnk 作业脚本与 `scTenifoldKnk_run.list` |

<details>
<summary>English</summary>
**Purpose**: Take the single-nucleus RNA-seq data on the Gifu T2T genome (the two root stages Dpi00 and Dpi02, two replicates each) from Cell Ranger quantification through clustering and cell-type annotation, cluster/cell-type DEGs and scTenifoldKnk virtual knockout.
**Entry point**: the seven shell scripts in the table below; the `Script/*.py` files with no entry script, run by hand, are `03_1-cluster_statitcs.py`, `03_2-fliter_to_all_gene_h5ad.py`, `05_3-marker_cluster_DEGs_own_order.py` and `makesure_02dpi_samples_all_have_Early_signal_cell.py`.
**Inputs**: `sample.list` (`ProjectID CNName RunID SampleName`), the four `*_filtered.h5ad` objects, Gifu T2T FASTA/GTF and `.pep.fa`, `{Gmax,Gm2.1,GmZH,At}_marker_genes.pep.fa`, eggNOG annotation; later steps read `02_processData/` and `04_Anno/` from the working directory, and scTenifoldKnk additionally needs `zz_00_pea.rds` and `All_species_merged_by_GifuT2T.SNF_Symbol.tsv`.
**Outputs**: `00_ref/`, `01_cellranger_out/<SampleName>/`, `02_processData/` (h5ad files, `zz_*`, `zz_00_pea.rds`, `00_02_pseudobulk_fpkm.csv`), `03_figures/` (dimension-reduction plots, `clustree_plot_*.pdf`, volcano plots), `04_Anno/` (Seurat rds, marker CSVs, `blast_*/`, enrichment results) and `06_batch_multiple_scTenifoldKnk/<group>/<gene>/`.
**Run**: from the stage root (holding `02_processData/`, `03_figures/`, `04_Anno/`) run the scripts in table order; the numbered scripts carry `#CSUB` directives and are submitted to the scheduler; for the ``{TOKEN}`` path placeholders see `scripts/README.md`; note that `03_1-dimenRedu_adjusted.py` writes a file name that `03_1-cluster_statitcs.py` and `03_2-fliter_to_all_gene_h5ad.py` do not read, so it must be unified first.
**Tools**: Cell Ranger 9.0.1; Python 3.9.19 with scanpy/anndata/scvi-tools/leidenalg/rpy2; R 4.4.1 with Seurat/scTenifoldKnk/DESeq2/clustree; BLAST+ and seqkit; settings as coded: `n_top_genes` 4000 (scanned 1000–5000), Leiden 0.4 (scanned 0.2–1.2), and 100 scTenifoldKnk replicates per gene (K ∈ {10,30,50,70,100}).
| Script | Purpose |
|---|---|
| `00_mkref.sh` | `cellranger mkref --genome=Gifu` from the Gifu T2T FASTA + GTF → `./00_ref` |
| `00_cellranger.sh` | one `cellranger count` per sample in `sample.list` → `01_cellranger_out/<SampleName>` |
| `run_singlecell.sh` | one driver running `Script/01_annaDataQC.py`, `02_integrate.py`, `03_1-dimenRedu_adjusted.py`, `05_3-markerAnno_Ps_v5_own.py`, `06_4-DEGs_cell_type.R` in order |
| `05_cluster_DEGs_finder.sh` | `05_1-…py` then `05_2-…R`: Seurat object and cluster markers, then per-cluster enrichment |
| `06_1.marker_gene_finder.sh` | DEGs → peptide FASTA → `blastp` against the Gmax/Gm2.1/GmZH/At markers, results merged |
| `06_1.Individual_clusters_marker_gene_finder.sh` | top-20 DEG sequences per cluster plus eggNOG annotation (active cluster `8`) |
| `07_batch_multiple_scTenifoldKnk.sh` | one scTenifoldKnk job script per candidate gene plus `scTenifoldKnk_run.list` |

</details>

### 本阶段步骤

| 步骤 | 内容 |
| --- | --- |
| [单核分析脚本](/_projects/lotus/scripts/08-single-nucleus/script/) | — |


