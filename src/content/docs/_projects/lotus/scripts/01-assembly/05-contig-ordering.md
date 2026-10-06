---
title: 'contig 排序'
sidebar:
  label: 'contig 排序'
  order: 5
---


## 中文
**作用**：用 nucmer 把组装 contig 比对到参考组装，推导每个 query contig 的主导参考 contig、比对长度、比例与链向，用于排序与定向。
**入口**：`05.1.nucmer.sh` → `05.2.contigOrder.py`
**输入**：`Gifu_ref.fa`（参考）与 `Gifu_contig.fa`（组装 contig），文件名硬编码并在工作目录解析。
**输出**：`Gifu_scaffold.delta`、`Gifu_scaffold.filter.delta`、`Gifu_scaffold.filter.delta.coords`、`Gifu_scaffold.delta.coords`、一张 `mummerplot` PNG；`mummer_analysis_results/dominant_mapping.csv` 与 `.../query_contig_on_dominant_ref_colored.pdf`（`+` 蓝、`-` 红）。
**运行**：`bash 05.1.nucmer.sh`，然后 `python 05.2.contigOrder.py`
**工具**：MUMmer `nucmer --prefix Gifu_scaffold`、`delta-filter -i 89 -l 1000 -1`、`show-coords -THrd`、`mummerplot --layout --png --large`（需 gnuplot）；Python 3 + pandas/matplotlib（`read_csv(sep='\t', skiprows=5)`）。

<details>
<summary>English</summary>
**Purpose**: Align the assembled contigs to a reference with nucmer and derive each query contig's dominant reference contig, aligned length, ratio and strand, for ordering and orientation.
**Entry point**: `05.1.nucmer.sh` → `05.2.contigOrder.py`
**Inputs**: `Gifu_ref.fa` (reference) and `Gifu_contig.fa` (contigs) — hardcoded filenames resolved in the working directory.
**Outputs**: `Gifu_scaffold.delta`, `Gifu_scaffold.filter.delta`, `Gifu_scaffold.filter.delta.coords`, `Gifu_scaffold.delta.coords`, a `mummerplot` PNG; `mummer_analysis_results/dominant_mapping.csv` and `.../query_contig_on_dominant_ref_colored.pdf` (`+` blue, `-` red).
**Run**: `bash 05.1.nucmer.sh`, then `python 05.2.contigOrder.py`
**Tools**: MUMmer `nucmer --prefix Gifu_scaffold`, `delta-filter -i 89 -l 1000 -1`, `show-coords -THrd`, `mummerplot --layout --png --large` (needs gnuplot); Python 3 with pandas and matplotlib (`read_csv(sep='\t', skiprows=5)`).

</details>


