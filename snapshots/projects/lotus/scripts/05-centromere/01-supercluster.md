---
title: '超级簇合并'
sidebar:
  label: '超级簇合并'
  order: 1
---


## 中文

**作用**：TRF 全基因组找串联重复，每个重复家族取一个代表单元用 CD-HIT-EST 聚类，代表序列用 BLASTN 回贴基因组，命中区间 Jaccard ≥ 0.7 的家族用并查集并入超级簇，只在 100 kb bin 中足够致密时才保留。
**入口**：`trf.sh` → `blast/blast.sh` → `blast/cluster_density.sh`（驱动）
**输入**：`Lotus_GifuT2T_v1.0.fasta`；`Gifu_cluster_reps.fa`；`Gifu_cluster_reps.blast`（第 1–4 列 = cluster、chr、start、end）；`Gifu_cluster_representatives.txt`（第 4 列 = cluster、第 5 列 = 单元序列）。
**输出**：驱动脚本第三个参数指定目录（默认 `supercluster`）下的 `windows.100000.bed`、`clusters.list`、`merged/*.bed`、`pairs.jaccard_ge0.7.tsv`、`cluster2super.tsv`、`super_summary.tsv`、`superbeds/*.bed`、`density/density.long.tsv`、`kept_supers.tsv`、`superunits/*.units.bed`。
**运行**：`bash trf.sh`；`bash blast/blast.sh`；`bash blast/cluster_density.sh Gifu_cluster_reps.blast {PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Lotus_GifuT2T_v1.0.fasta supercluster`（`trf.sh` 与 `blast.sh` 多为注释记录，仅末行有效；`cluster_density.sh` 须在 `blast/` 下运行）。
**工具**：TRF `2 7 7 80 10 50 500 -d -h`；CD-HIT-EST `-c 0.85 -n 3 -d 0 -T 8 -M 16000`；BLAST+ `makeblastdb -dbtype nucl`、`blastn -task blastn -evalue 1e-10 -perc_identity 80 -num_threads 64`；bedtools `makewindows -w 100000`、`merge -d 50`、`jaccard`；阈值 `WIN=100000`、`JAC=0.7`、`DENS_CUT=0.5`、`BIN_KEEP=5`。

<details>
<summary>English</summary>
**Purpose**: TRF calls tandem repeats genome-wide, one representative unit per repeat family is clustered with CD-HIT-EST, the representatives are mapped back with BLASTN, and families whose merged genomic hits reach a Jaccard ≥ 0.7 are merged with Union-Find into superclusters, kept only when dense enough in 100-kb bins.
**Entry point**: `trf.sh` → `blast/blast.sh` → `blast/cluster_density.sh` (driver)
**Inputs**: `Lotus_GifuT2T_v1.0.fasta`; `Gifu_cluster_reps.fa`; `Gifu_cluster_reps.blast` (columns 1–4 = cluster, chr, start, end); `Gifu_cluster_representatives.txt` (column 4 = cluster, column 5 = unit sequence).
**Outputs**: in the directory given as the third argument of the driver (default `supercluster`): `windows.100000.bed`, `clusters.list`, `merged/*.bed`, `pairs.jaccard_ge0.7.tsv`, `cluster2super.tsv`, `super_summary.tsv`, `superbeds/*.bed`, `density/density.long.tsv`, `kept_supers.tsv`, `superunits/*.units.bed`.
**Run**: `bash trf.sh`; `bash blast/blast.sh`; `bash blast/cluster_density.sh Gifu_cluster_reps.blast {PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Lotus_GifuT2T_v1.0.fasta supercluster` (`trf.sh` and `blast.sh` are mostly commented records with only their last line live; `cluster_density.sh` must run from `blast/`).
**Tools**: TRF `2 7 7 80 10 50 500 -d -h`; CD-HIT-EST `-c 0.85 -n 3 -d 0 -T 8 -M 16000`; BLAST+ `makeblastdb -dbtype nucl`, `blastn -task blastn -evalue 1e-10 -perc_identity 80 -num_threads 64`; bedtools `makewindows -w 100000`, `merge -d 50`, `jaccard`; thresholds `WIN=100000`, `JAC=0.7`, `DENS_CUT=0.5`, `BIN_KEEP=5`.

</details>


