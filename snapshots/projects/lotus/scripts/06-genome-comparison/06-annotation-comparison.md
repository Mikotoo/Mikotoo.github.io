---
title: '注释比较'
sidebar:
  label: '注释比较'
  order: 6
---


## 中文

**作用**：比较 accession 之间的重复/转座子注释；`run.sh` 中启用的一步筛出落在特殊区域（着丝粒、rDNA、端粒）之外的 LTR 插入。
**入口**：`run.sh`
**输入**：`LTR_insertion.bed`、`LTR_in_special_regions.bed`；这两个文件由脚本中被注释的两个上游步骤生成。
**输出**：`LTR_outside.bed` — 在 `LTR_in_special_regions.bed` 中无重叠的 LTR 插入。
**运行**：`bash run.sh`（`bedtools intersect -v -a LTR_insertion.bed -b LTR_in_special_regions.bed > LTR_outside.bed`；两个输入须先存在，且须在存放它们的目录中运行）。
**工具**：bedtools 2.31.1（`intersect -v`；被注释的步骤用 `intersect -wa`）、awk（把 `MG20_LTR.insertionTime.txt` 转成 `LTR_insertion.bed`，并修正 start > end 的记录）。

<details>
<summary>English</summary>
**Purpose**: Compare the repeat/transposon annotation between the accessions; the single active step in `run.sh` selects the LTR insertions that fall outside the special regions (centromere, rDNA, telomere).
**Entry point**: `run.sh`
**Inputs**: `LTR_insertion.bed` and `LTR_in_special_regions.bed`; both are produced by the two upstream steps that are commented out in the script.
**Outputs**: `LTR_outside.bed` — LTR insertions with no overlap in `LTR_in_special_regions.bed`.
**Run**: `bash run.sh` (`bedtools intersect -v -a LTR_insertion.bed -b LTR_in_special_regions.bed > LTR_outside.bed`; both inputs must exist first and the script must run from the directory holding them).
**Tools**: bedtools 2.31.1 (`intersect -v`; the commented step uses `intersect -wa`), awk (converting `MG20_LTR.insertionTime.txt` into `LTR_insertion.bed`, including a swap fix for records where start > end).

</details>


