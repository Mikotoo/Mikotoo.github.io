---
title: '已发表组装'
sidebar:
  label: '已发表组装'
  order: -1
---


## 中文

**作用**：用 MUMmer 点阵图比较已发表的 Gifu 与 MG20 组装；gap 内容比较见 `gap_analysis/`。
**入口**：`mummer.sh`
**输入**：`ref.filter.delta`；`-R `{PROJ_LOTUS_ZC}`/Gifu/05_3.mummer/Gifu_ref.fa`、`-Q `{PROJ_LOTUS_ZC}`/MG20/05_3.mummer/MG20_ref.fa`。
**输出**：前缀 `ref.filter` 的 `mummerplot` PNG；被注释的步骤产出 `ref.filter.delta` 与 `ref.filter.delta.coords`。
**运行**：`bash mummer.sh`（仅 `mummerplot` 一行启用；缺 `ref.filter.delta` 时先运行注释中的 `nucmer`、`delta-filter`、`show-coords`）。
**工具**：MUMmer4 4.0.1 — `nucmer --prefix ref`、`delta-filter -i 89 -l 1000 -1`、`show-coords -THrd`、`mummerplot --layout --png --large`；CSUB（`-q c01`、`-n 64`、`span[hosts=1]`）。

<details>
<summary>English</summary>
**Purpose**: Dot-plot comparison of the published Gifu and MG20 assemblies; the gap content comparison is in `gap_analysis/`.
**Entry point**: `mummer.sh`
**Inputs**: `ref.filter.delta`; `-R `{PROJ_LOTUS_ZC}`/Gifu/05_3.mummer/Gifu_ref.fa`, `-Q `{PROJ_LOTUS_ZC}`/MG20/05_3.mummer/MG20_ref.fa`.
**Outputs**: the `mummerplot` PNG for prefix `ref.filter`; the commented steps produce `ref.filter.delta` and `ref.filter.delta.coords`.
**Run**: `bash mummer.sh` (only the `mummerplot` line is active; when `ref.filter.delta` is absent, run the commented `nucmer`, `delta-filter` and `show-coords` steps first).
**Tools**: MUMmer4 4.0.1 — `nucmer --prefix ref`, `delta-filter -i 89 -l 1000 -1`, `show-coords -THrd`, `mummerplot --layout --png --large`; CSUB (`-q c01`, `-n 64`, `span[hosts=1]`).

</details>

### 本阶段步骤

| 步骤 | 内容 |
| --- | --- |
| [gap 区域分析](/_projects/lotus/scripts/06-genome-comparison/01-published-assemblies/gap-analysis/) | — |


