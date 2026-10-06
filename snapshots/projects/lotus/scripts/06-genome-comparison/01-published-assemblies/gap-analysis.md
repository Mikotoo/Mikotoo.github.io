---
title: 'gap 区域分析'
sidebar:
  label: 'gap 区域分析'
  order: 999
---


## 中文

**作用**：比较已发表的 Gifu (v1.2) 与 MG20 (gnm3) 组装的 gap 数目与长度，并写出统计表与对比图。
**入口**：`gap_analysis.py`
**输入**：与 `gap_analysis/` 同级的 `Gifu/Gifu_gap.gff`、`Gifu/Lotusjaponicus_Gifu_v1.2_genome.fa.fai`、`MG20/Lotus_MG20_gap.gff`、`MG20/Lotus_japonicus.fasta.fai`（GFF3 中 type 为 `gap`，1-based 闭区间）；染色体长度取自 `.fai`。
**输出**：`results/`（`{Gifu,MG20}_gap_summary_per_chromosome.csv`、`{key}_gap_positions.csv`、`{key}_gap_size_distribution.csv`、`gap_genome_comparison.csv`）与 `figures/`（`{key}_gap_per_chromosome.png`、`{key}_gap_size_distribution.png`、`{key}_gap_map.png`、`comparison_genome_level.png`、`comparison_per_chromosome.png`、`comparison_size_ecdf.png`）。
**运行**：`python gap_analysis.py`（无绝对路径：数据、`results/`、`figures/` 均由脚本位置推导，`Gifu/`、`MG20/` 必须与 `gap_analysis/` 同级）。
**工具**：Python 3 + pandas、numpy、matplotlib（`Agg`）；字体优先 `WenQuanYi Micro Hei`，回退 `DejaVu Sans`；染色体排序正则 `(?:chr|Lj)(\d+|[A-Za-z]+)$`。

<details>
<summary>English</summary>
**Purpose**: Compare gap number and length between the published Gifu (v1.2) and MG20 (gnm3) assemblies and write the statistics tables and comparison plots.
**Entry point**: `gap_analysis.py`
**Inputs**: in the directory next to `gap_analysis/`: `Gifu/Gifu_gap.gff`, `Gifu/Lotusjaponicus_Gifu_v1.2_genome.fa.fai`, `MG20/Lotus_MG20_gap.gff`, `MG20/Lotus_japonicus.fasta.fai` (the GFF3 records with type `gap`, 1-based closed intervals); chromosome lengths come from the `.fai`.
**Outputs**: `results/` (`{Gifu,MG20}_gap_summary_per_chromosome.csv`, `{key}_gap_positions.csv`, `{key}_gap_size_distribution.csv`, `gap_genome_comparison.csv`) and `figures/` (`{key}_gap_per_chromosome.png`, `{key}_gap_size_distribution.png`, `{key}_gap_map.png`, `comparison_genome_level.png`, `comparison_per_chromosome.png`, `comparison_size_ecdf.png`).
**Run**: `python gap_analysis.py` (no absolute paths: the data, `results/` and `figures/` are derived from the script location, so `Gifu/` and `MG20/` must sit next to `gap_analysis/`).
**Tools**: Python 3 with pandas, numpy and matplotlib (`Agg`); font preference `WenQuanYi Micro Hei`, falling back to `DejaVu Sans`; chromosome order regex `(?:chr|Lj)(\d+|[A-Za-z]+)$`.

</details>


