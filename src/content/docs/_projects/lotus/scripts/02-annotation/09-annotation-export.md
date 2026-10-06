---
title: '注释导出'
sidebar:
  label: '注释导出'
  order: 9
---


## 中文

**作用**：导出最终注释，并用导出的基因 BED 与自身相交，得到用于检查导出基因集的非冗余重叠对。
**入口**：`bedtools.sh`。
**输入**：`Lotus_Gifu_gene.bed`，即从最终 Gifu T2T 注释导出的基因区间，第 4 列必须是基因标识符。
**输出**：`overlaps.tsv`，每个重叠基因对一行，格式为 `bedtools intersect -wa -wb`（第 1–4 列来自第一个区间，第 5–9 列来自第二个）。
**运行**：`bash bedtools.sh`，即 `bedtools intersect -a Lotus_Gifu_gene.bed -b Lotus_Gifu_gene.bed -wa -wb | awk -F'\t' '$4 < $9' > overlaps.tsv`；脚本无 shebang 也无 `cd`，须在存放 `Lotus_Gifu_gene.bed` 的目录中运行。
**工具**：bedtools `intersect -a Lotus_Gifu_gene.bed -b Lotus_Gifu_gene.bed -wa -wb`；awk `-F'\t' '$4 < $9'`（每对只保留一次）。

<details>
<summary>English</summary>
**Purpose**: Export the final annotation and intersect the exported gene BED with itself to obtain the non-redundant overlapping pairs used to check the exported gene set.
**Entry point**: `bedtools.sh`.
**Inputs**: `Lotus_Gifu_gene.bed`, the gene intervals exported from the final Gifu T2T annotation, with the gene identifier in column 4.
**Outputs**: `overlaps.tsv`, one line per overlapping gene pair, in `bedtools intersect -wa -wb` format (columns 1–4 from the first interval, columns 5–9 from the second).
**Run**: `bash bedtools.sh`, i.e. `bedtools intersect -a Lotus_Gifu_gene.bed -b Lotus_Gifu_gene.bed -wa -wb | awk -F'\t' '$4 < $9' > overlaps.tsv`; the script has no shebang and no `cd`, so it must run in the directory holding `Lotus_Gifu_gene.bed`.
**Tools**: bedtools `intersect -a Lotus_Gifu_gene.bed -b Lotus_Gifu_gene.bed -wa -wb`; awk `-F'\t' '$4 < $9'` (keeps each pair once).

</details>


