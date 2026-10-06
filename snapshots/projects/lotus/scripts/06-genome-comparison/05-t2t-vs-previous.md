---
title: 'T2T 与前版本比较'
sidebar:
  label: 'T2T 与前版本比较'
  order: 5
---


## 中文

**作用**：把每个新 T2T 组装与同一 accession 的已发表旧组装比较（Gifu T2T vs Gifu old；MG20 T2T vs MG20 old）：MUMmer4 + SyRI 变异、由比对区间补集得到的 PAV，以及 T2T 特有新序列按注释区域的分区。
**入口**：`GifuT2T_GifuOld/GifuT2T_vs_GifuOld_variation_finder.sh`
**输入**：Gifu 用 `Lotus_Gifu_Old.fa`（seqkit 清理为 `Lotus_Gifu_Old1.fa`）与 `Lotus_GifuT2T_v1.0.fasta`，MG20 用 `Lotus_MG20_Old.fa` 与 `Lotus_MG20T2T_v1.0.fasta`，以及 `genomes.txt`；PAV 需对应的 `<prefix>.filtered.coords`；`calc_bin_coverage.sh` 用 `gifuT2T_newseq.bed` 与 `{PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Gifu_100k.bed`；`run.sh`/`run.py` 用 `gifuT2T_newseq.bed`（或 `mg20T2T_newseq.bed`）与 `cent.bed`、`rDNA.bed`、`tel.bed`、`repeat2.bed`。
**输出**：`<p>.delta`/`.filtered.delta`、dnadiff 报告、`syri.out`、`Old_<type>.bed` 与 `T2T_<type>.bed`（SNP、INS、DEL、INV、TRANS、INVTR、DUP、INVDP、PAV）、stdout 计数表与 `*_findSV.pdf`；PAV 出 `.genome`、`<asm>.aligned.bed`、`<asm>.PAV.bed`、`PAV.summary`；`calc_bin_coverage.sh` 出 `gifu_100kb_coverage.tsv`；`run.sh`/`run.py` 出各区域 bp 与占比及 `new_sequences_region_proportions.pdf`。
**运行**：每个比较各自一个工作目录，先跑 variation finder 再跑 PAV：`cd GifuT2T_GifuOld && bash GifuT2T_vs_GifuOld_variation_finder.sh && bash GifuT2T_vs_GifuOld_get_PAV_from_coords.sh`（MG20 同理）；随后 `bash calc_bin_coverage.sh`，并在 `Gifu_new/`、`MG20_new/` 中各跑 `bash run.sh` 或 `python3 run.py`。
**工具**：MUMmer4 4.0.1（`nucmer -c 1000 --mum --maxgap=1000`、`delta-filter -i 90 -l 1000 -q`、`show-coords -THrd`、`dnadiff -d -p`）、SyRI 1.7.1、plotsr 1.1.1、samtools 1.21、bedtools 2.31.1（`merge -d 50`、`complement -g`、`intersect`、`subtract`）、seqkit 2.10.1（`grep -v -r -p`）、GNU Awk 5.1.0、Python 3 + pandas/matplotlib、CSUB（`-q c01`、`-n 64`）。

<details>
<summary>English</summary>
**Purpose**: Compare each new T2T assembly with the published earlier assembly of the same accession (Gifu T2T vs Gifu old; MG20 T2T vs MG20 old): MUMmer4 + SyRI variants, PAV from the complement of the aligned intervals, and a partition of the T2T-specific new sequences by the annotated regions they fall in.
**Entry point**: `GifuT2T_GifuOld/GifuT2T_vs_GifuOld_variation_finder.sh`
**Inputs**: for Gifu, `Lotus_Gifu_Old.fa` (cleaned by seqkit to `Lotus_Gifu_Old1.fa`) and `Lotus_GifuT2T_v1.0.fasta`; for MG20, `Lotus_MG20_Old.fa` and `Lotus_MG20T2T_v1.0.fasta`; plus `genomes.txt`; the PAV step needs the matching `<prefix>.filtered.coords`; `calc_bin_coverage.sh` uses `gifuT2T_newseq.bed` and `{PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Gifu_100k.bed`; `run.sh`/`run.py` use `gifuT2T_newseq.bed` (or `mg20T2T_newseq.bed`) with `cent.bed`, `rDNA.bed`, `tel.bed` and `repeat2.bed`.
**Outputs**: `<p>.delta`/`.filtered.delta`, dnadiff reports, `syri.out`, `Old_<type>.bed` and `T2T_<type>.bed` (`type` = SNP, INS, DEL, INV, TRANS, INVTR, DUP, INVDP, PAV), a count table on stdout and `*_findSV.pdf`; the PAV step writes `.genome`, `<asm>.aligned.bed`, `<asm>.PAV.bed` and `PAV.summary`; `calc_bin_coverage.sh` writes `gifu_100kb_coverage.tsv`; `run.sh`/`run.py` print per-region bp and ratios and write `new_sequences_region_proportions.pdf`.
**Run**: give each comparison its own working directory and run the variation finder before the PAV step: `cd GifuT2T_GifuOld && bash GifuT2T_vs_GifuOld_variation_finder.sh && bash GifuT2T_vs_GifuOld_get_PAV_from_coords.sh` (the same for MG20); then `bash calc_bin_coverage.sh`, and `bash run.sh` or `python3 run.py` in `Gifu_new/` and `MG20_new/`.
**Tools**: MUMmer4 4.0.1 (`nucmer -c 1000 --mum --maxgap=1000`, `delta-filter -i 90 -l 1000 -q`, `show-coords -THrd`, `dnadiff -d -p`), SyRI 1.7.1, plotsr 1.1.1, samtools 1.21, bedtools 2.31.1 (`merge -d 50`, `complement -g`, `intersect`, `subtract`), seqkit 2.10.1 (`grep -v -r -p`), GNU Awk 5.1.0, Python 3 with pandas and matplotlib, CSUB (`-q c01`, `-n 64`).

</details>


