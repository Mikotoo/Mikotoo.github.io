---
title: '易位核查'
sidebar:
  label: '易位核查'
  order: 3
---


## 中文

**作用**：用 read-pair 证据验证区分 Gifu 与 MG20 的 Chr1↔Chr2 易位是真实结构差异而非组装错误，并对多个样本批量分型：同一套 reads 比对到两个基因组，统计跨染色体 pair。
**入口**：`run.sh`
**输入**：`trans.conf`（3 列制表符分隔、无表头：`sample`、`bam_gifu`、`bam_mg`）及其中列出的 BAM；逐对计数另需坐标排序并建索引的 `gifu2gifu.bam`、`gifu2mg.bam`、`mg2gifu.bam`、`mg2mg.bam`；区域 `Chr1:22829936-22849936`（Gifu）与 `Chr2:99260160-99280160`（MG20）。
**输出**：`bam.list`；每样本 `{sample}.trans.gifu_anchored.merged_pairs.tsv`、`{sample}.trans.mg_anchored.merged_pairs.tsv` 与 4 个子集 BAM；汇总 `lotus_translocation_typing.tsv`（`sample`、`n_cross_in_gifu_bam`、`n_cross_in_mg_bam`、`call`、`ratio_MG20`、`note`）；逐对计数在 stdout 出计数表并写 `<bam>_chr1_chr2_pairs.txt` 等文件。
**运行**：`bash run.sh`（须在存放 `trans.conf`、BAM 与 TSV 的目录中运行；文件末尾的逐对计数命令以注释保留，需单独运行）。
**工具**：pysam（BAM 访问）、pandas、numpy；`samtools index`；Python 3（finder 经 `subprocess` 调用，`--min-support` 默认 5，两侧都达标为 `Mixed`，都低于为 `Uncertain`）；CSUB（`-q c01`、`-n 88`）。

<details>
<summary>English</summary>
**Purpose**: Read-pair evidence that the Chr1↔Chr2 translocation distinguishing Gifu from MG20 is a real structural difference rather than a mis-assembly, plus batch typing across samples: the same reads are mapped to both genomes and cross-chromosome pairs are counted.
**Entry point**: `run.sh`
**Inputs**: `trans.conf` (3 tab-separated columns, no header: `sample`, `bam_gifu`, `bam_mg`) and the BAMs it lists; per-pair counting also needs the coordinate-sorted, indexed `gifu2gifu.bam`, `gifu2mg.bam`, `mg2gifu.bam`, `mg2mg.bam`; regions `Chr1:22829936-22849936` (Gifu) and `Chr2:99260160-99280160` (MG20).
**Outputs**: `bam.list`; per sample `{sample}.trans.gifu_anchored.merged_pairs.tsv`, `{sample}.trans.mg_anchored.merged_pairs.tsv` and four subset BAMs; the summary `lotus_translocation_typing.tsv` (`sample`, `n_cross_in_gifu_bam`, `n_cross_in_mg_bam`, `call`, `ratio_MG20`, `note`); per-pair counting prints count tables to stdout and writes files such as `<bam>_chr1_chr2_pairs.txt`.
**Run**: `bash run.sh` (must run in the directory holding `trans.conf`, the BAMs and the TSVs; the per-pair counting commands are kept commented at the end of the file and must be run separately).
**Tools**: `pysam` (BAM access), pandas, numpy; `samtools index`; Python 3 (the finder is called through `subprocess`; `--min-support` defaults to 5, both sides at or above it give `Mixed`, below it on both sides gives `Uncertain`); CSUB (`-q c01`, `-n 88`).

</details>


