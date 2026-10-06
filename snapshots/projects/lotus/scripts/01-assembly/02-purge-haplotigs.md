---
title: '去冗余单倍型'
sidebar:
  label: '去冗余单倍型'
  order: 2
---


## 中文
**作用**：把 HiFi reads 回贴到 hifiasm primary contig 统计读深与覆盖度，并去除单倍型冗余（haplotig）contig。
**入口**：`02.1.minimap.sh`（比对）→ `02.2.purge.sh`（去冗余）
**输入**：`Gifu.hifiasm.fasta`（步骤 01 输出）；HiFi reads ``{PROJ_LOTUS_ZC}`/data/Gifu/hifi/Gifu_hifi.fastq.gz`。
**输出**：`aligned.bam`、`aligned.bam.200.gencov`、`coverage_stats.csv`；去冗余后的 contig FASTA 与 haplotig FASTA（`purge_haplotigs purge` 的默认名）。
**运行**：`bash 02.1.minimap.sh`，然后 `bash 02.2.purge.sh`
**工具**：minimap2 `-ax map-pb`；samtools `view -hF 256`、`sort -@ 20 -m 1G -T tmp.ali`；purge_haplotigs `readhist -t 88`、`contigcov -l 15 -m 100 -h 185`、`purge -t 88 -a 50`。
**注意**：`02.2.purge.sh` 中 `readhist` 与 `contigcov` 处于注释状态、仅 `purge` 生效，其读取的 `aligned.bam.200.gencov` 须按原名存在。

<details>
<summary>English</summary>
**Purpose**: Align HiFi reads back to the hifiasm primary contigs to measure read depth and coverage, then purge haplotig (duplicated) contigs.
**Entry point**: `02.1.minimap.sh` (alignment) → `02.2.purge.sh` (purging)
**Inputs**: `Gifu.hifiasm.fasta` (step 01 output); HiFi reads ``{PROJ_LOTUS_ZC}`/data/Gifu/hifi/Gifu_hifi.fastq.gz`.
**Outputs**: `aligned.bam`, `aligned.bam.200.gencov`, `coverage_stats.csv`; the purged contig FASTA and the haplotig FASTA (the `purge_haplotigs purge` default names).
**Run**: `bash 02.1.minimap.sh`, then `bash 02.2.purge.sh`
**Tools**: minimap2 `-ax map-pb`; samtools `view -hF 256`, `sort -@ 20 -m 1G -T tmp.ali`; purge_haplotigs `readhist -t 88`, `contigcov -l 15 -m 100 -h 185`, `purge -t 88 -a 50`.
**Notes**: In `02.2.purge.sh`, `readhist` and `contigcov` are commented out and only `purge` is active, so the `aligned.bam.200.gencov` it reads must exist under that exact name.

</details>


