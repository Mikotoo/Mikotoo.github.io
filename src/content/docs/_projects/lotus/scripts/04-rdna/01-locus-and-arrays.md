---
title: 'rDNA 位点与阵列'
sidebar:
  label: 'rDNA 位点与阵列'
  order: 1
---


## 中文

**作用**：用 Barrnap 注释 Gifu T2T 的 rRNA 基因以确定 rDNA 位点，切出各阵列做 StainedGlass 自比对，并统计每个阵列的 18S/5.8S/28S/5S 特征数以得到完整 45S 重复单元数。
**入口**：`barrnap.sh`
**输入**：`Lotus_GifuT2T_v1.0.fasta`（含随附 `.fai`，Chr1–Chr6）；区段 BED `Chr2_0_3Mbp.bed`、`Chr2_24_26Mbp.bed`、`Chr5_14_16Mbp.bed`、`Chr6_8_10Mbp.bed`、`Chr6_25_28Mbp.bed`；Barrnap 的 `rRNA.gff`/`rRNA.fa`；生产 GFF `Gifu_rRNA.gff`、`MG20_rRNA.gff`。
**输出**：`rRNA.gff`/`rRNA.fa`；`result/*_output.tbl.gz`、`*_output.bed.gz` 及 `.fa`/`.fai`/`.bed`/`.fasta`/`.mmi`/`.bam` 中间文件；`results/` 下各阵列的 PDF/PNG 图；阵列组成表打印到标准输出并写入 `03_annotation/rDNA_array_summary.tsv`（45S 完整单元数 = 阵列内 `min(18S, 5.8S, 28S)`，阵列坐标硬编码在 `count_rDNA.py` 中）。
**运行**：`bash barrnap.sh` → `python3 count_rDNA.py` → `bash result/process.sh`（`process.sh` 从自身目录读 `../Lotus_GifuT2T_v1.0.fasta` 与区段 BED，结果写到当前目录，须在 `result/` 下运行）。
**工具**：Barrnap `--kingdom euk --threads 88`；minimap2 建索引 `-f 1000 -s 400 -ax ava-ont`、比对 `-t 88 -f 10000 -s 400 -ax ava-ont --dual=yes --eqx`；StainedGlass 0.6 `samIdentity.py --threads 88 --matches 400 --header`、`refmt.py --window 2000`、`aln_plot.R --threads 88`；bedtools `getfasta`、`makewindows -w 2000`；samtools `faidx`、`sort -m 4G`；`bgzip`、`Rscript`、`python3.13`（经 ``{SOFTWARE_LOCAL}``）；`count_rDNA.py` 仅用标准库。

<details>
<summary>English</summary>
**Purpose**: Annotate the Gifu T2T rRNA genes with Barrnap to define the rDNA loci, extract each array for StainedGlass self-identity analysis, and count the 18S/5.8S/28S/5S features per array to derive the number of complete 45S repeat units.
**Entry point**: `barrnap.sh`
**Inputs**: `Lotus_GifuT2T_v1.0.fasta` (with the shipped `.fai`, Chr1–Chr6); the region BEDs `Chr2_0_3Mbp.bed`, `Chr2_24_26Mbp.bed`, `Chr5_14_16Mbp.bed`, `Chr6_8_10Mbp.bed`, `Chr6_25_28Mbp.bed`; Barrnap's `rRNA.gff`/`rRNA.fa`; the production GFFs `Gifu_rRNA.gff` and `MG20_rRNA.gff`.
**Outputs**: `rRNA.gff`/`rRNA.fa`; `result/*_output.tbl.gz` and `*_output.bed.gz` plus the `.fa`/`.fai`/`.bed`/`.fasta`/`.mmi`/`.bam` intermediates; PDF/PNG plots per array under `results/`; the array composition table on stdout and in `03_annotation/rDNA_array_summary.tsv` (complete 45S units = `min(18S, 5.8S, 28S)` inside each array, with the array coordinates hardcoded in `count_rDNA.py`).
**Run**: `bash barrnap.sh` → `python3 count_rDNA.py` → `bash result/process.sh` (`process.sh` reads `../Lotus_GifuT2T_v1.0.fasta` and the region BEDs from its own directory and writes into the current directory, so it must run inside `result/`).
**Tools**: Barrnap `--kingdom euk --threads 88`; minimap2 index `-f 1000 -s 400 -ax ava-ont`, alignment `-t 88 -f 10000 -s 400 -ax ava-ont --dual=yes --eqx`; StainedGlass 0.6 `samIdentity.py --threads 88 --matches 400 --header`, `refmt.py --window 2000`, `aln_plot.R --threads 88`; bedtools `getfasta`, `makewindows -w 2000`; samtools `faidx`, `sort -m 4G`; `bgzip`, `Rscript`, `python3.13` (through ``{SOFTWARE_LOCAL}``); `count_rDNA.py` is standard library only.

</details>


