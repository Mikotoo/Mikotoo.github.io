---
title: '端粒检查'
sidebar:
  label: '端粒检查'
  order: 4
---


## 中文
**作用**：扫描端粒重复单元并列出富集它们的 1 kb 窗口，另独立运行 BUSCO 评估基因区完整度。
**入口**：`find_tel.sh`（调用 `get_tel.py`）；`busco.sh` 独立运行
**输入**：`contigs.fa`（端粒扫描与 `busco.sh` 共用）。
**输出**：`Gifu.tel_unit.bed`、`contigs.fa.fai`、`Gifu.length`、`Gifu_1k.bed`、`Gifu_tel.result`（每窗口 >50 命中）；BUSCO 输出目录 `Gifu_lotus`。
**运行**：`bash find_tel.sh`；`bash busco.sh`
**工具**：Python + Biopython `SeqIO`（基序 `CCCTAAA`、`TTTAGGG` 及各位置的单碱基替换）；samtools `faidx`；bedtools `makewindows -w 1000`、`intersect -wa`；BUSCO `-l embryophyta_odb10 -m genome --cpu 64`。

<details>
<summary>English</summary>
**Purpose**: Scan the contig set for telomeric repeat units and list the 1 kb windows dense in them, and run BUSCO separately for gene-space completeness.
**Entry point**: `find_tel.sh` (calls `get_tel.py`); `busco.sh` runs independently
**Inputs**: `contigs.fa` (shared by the telomere scan and `busco.sh`).
**Outputs**: `Gifu.tel_unit.bed`, `contigs.fa.fai`, `Gifu.length`, `Gifu_1k.bed`, `Gifu_tel.result` (>50 hits per window); BUSCO output directory `Gifu_lotus`.
**Run**: `bash find_tel.sh`; `bash busco.sh`
**Tools**: Python with Biopython `SeqIO` (motifs `CCCTAAA` and `TTTAGGG` plus all single-base substitutions); samtools `faidx`; bedtools `makewindows -w 1000`, `intersect -wa`; BUSCO `-l embryophyta_odb10 -m genome --cpu 64`.

</details>


