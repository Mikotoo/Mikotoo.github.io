---
title: 'StainedGlass 分析'
sidebar:
  label: 'StainedGlass 分析'
  order: 4
---


## 中文

**作用**：用 StainedGlass 自比对展示 Gifu 着丝粒重复阵列的高阶结构，每条染色体一次：切出单条染色体、切成 2 kb 窗口、minimap2 全比全自比对，再用 StainedGlass 脚本转成 identity 表并绘图。
**入口**：`run_stainglass.sh`
**输入**：各 `Gifu_Chr<N>/` 内的 `Gifu_Chr<N>_2k.fasta` 与其 `.fai`；`get.sh` 另需 `mg20_cent.fa` 与 `.fai`。
**输出**：每个 `Gifu_Chr<N>/` 内的 `Gifu_Chr<N>_2k.mmi`、`Gifu_Chr<N>.bam`、`.tbl.gz`、`.full.tbl.gz`、`.bed.gz`，图件在 `results/Gifu_Chr<N>_figures/{pdfs,pngs}/`；`get.sh` 每条染色体写 `<id>.fa`、`<id>.fa.fai`、`<id>_2k.bed`、`<id>_2k.fasta`。
**运行**：`bash run_stainglass.sh`（在 `04_stainglass/` 下依次进入 `Gifu_Chr1`…`Gifu_Chr6`，每目录内建索引与自比对 → `samIdentity.py` → `bgzip` → `refmt.py` → `aln_plot.R`）。
**工具**：StainedGlass 0.6（`samIdentity.py --threads 88 --matches 400 --header`、`refmt.py --window 2000`、`aln_plot.R --threads 88`，工作流脚本来自 `SG=~/software/StainedGlass-0.6/workflow/scripts`）；minimap2 `-f 1000 -s 400 -ax ava-ont` 建索引、`-t 88 -f 10000 -s 400 -ax ava-ont --dual=yes --eqx` 比对；samtools `faidx`、`sort -m 4G`；bedtools `makewindows -w 2000`、`getfasta`；`bgzip`、`Rscript`。

<details>
<summary>English</summary>
**Purpose**: StainedGlass self-alignments showing the higher-order organisation of the Gifu centromeric repeat arrays, one chromosome at a time: each chromosome is cut out, windowed at 2 kb and aligned to itself all-versus-all with minimap2, then converted to identity tables with the StainedGlass scripts and plotted.
**Entry point**: `run_stainglass.sh`
**Inputs**: the 2-kb window FASTA and its `.fai` inside each `Gifu_Chr<N>/` directory (`Gifu_Chr<N>_2k.fasta`, `Gifu_Chr<N>.fa.fai`); `get.sh` additionally needs `mg20_cent.fa` and its `.fai`.
**Outputs**: inside each `Gifu_Chr<N>/` directory `Gifu_Chr<N>_2k.mmi`, `Gifu_Chr<N>.bam`, `.tbl.gz`, `.full.tbl.gz`, `.bed.gz`, with plots under `results/Gifu_Chr<N>_figures/{pdfs,pngs}/`; `get.sh` writes `<id>.fa`, `<id>.fa.fai`, `<id>_2k.bed` and `<id>_2k.fasta` per chromosome.
**Run**: `bash run_stainglass.sh` (from `04_stainglass/`, entering `Gifu_Chr1`…`Gifu_Chr6` in turn; each directory runs index and self alignment → `samIdentity.py` → `bgzip` → `refmt.py` → `aln_plot.R`).
**Tools**: StainedGlass 0.6 (`samIdentity.py --threads 88 --matches 400 --header`, `refmt.py --window 2000`, `aln_plot.R --threads 88`, with the workflow scripts taken from `SG=~/software/StainedGlass-0.6/workflow/scripts`); minimap2 index `-f 1000 -s 400 -ax ava-ont`, alignment `-t 88 -f 10000 -s 400 -ax ava-ont --dual=yes --eqx`; samtools `faidx`, `sort -m 4G`; bedtools `makewindows -w 2000`, `getfasta`; `bgzip`, `Rscript`.

</details>


