---
title: '品种间比对'
sidebar:
  label: '品种间比对'
  order: 2
---


## 中文

**作用**：MUMmer4 + SyRI 比较 MG20T2T（参考）与 GifuT2T（查询），输出各类变异 BED、事件与 bp 计数以及由比对区间补集得到的 PAV；`syri.sh` 是另一条更早的 minimap2 路线，覆盖三对基因组并合成一张 plotsr 图。
**入口**：`GifuT2T_vs_MG20T2T_variation_finder.sh`
**输入**：`Lotus_MG20T2T_v1.0.fasta`、`Lotus_GifuT2T_v1.0.fasta`（裸文件名）；PAV 另需上一步产出的 `GifuT2T_MG20T2T.filtered.coords`；`syri.sh` 用 ``{PROJ_LOTUS_ZC}`/comp/reference/` 下的 `Gifu_v1.fa`、`Gifu_T2T.fa`、`MG20_T2T.fa`、`MG20_v1.fa` 与 `genomes.txt`。
**输出**：`GifuT2T_MG20T2T.delta`/`.filtered.delta`、dnadiff 报告、`.filtered.coords`、`syri.out`、`MG20_<type>.bed` 与 `Gifu_<type>.bed`（SNP、INS、DEL、INV、TRANS、INVTR、DUP、INVDP、PAV）、计数表与 `variation.out`；PAV 出 `<asm>.aligned.bed`、`<asm>.PAV.bed`、`PAV.summary`；`syri.sh` 出三个 BAM 及 `.bai` 与 `syri_samechromosome_plot.pdf`。
**运行**：`bash GifuT2T_vs_MG20T2T_variation_finder.sh`，随后 `bash GifuT2T_vs_MG20T2T_get_PAV_from_coords.sh`；`bash syri.sh` 为独立路线（每个比较须在各自工作目录运行，否则 `syri.out`、BED 与临时文件互相覆盖）。
**工具**：MUMmer4 4.0.1（`nucmer -c 1000 --mum --maxgap=1000`、`delta-filter -i 90 -l 1000 -q`、`show-coords -THrd`、`dnadiff -d -p`）、SyRI 1.7.1、plotsr 1.1.1、samtools 1.21、bedtools 2.31.1（`merge -d 50`、`complement -g`）、minimap2（`-ax asm20`，仅 `syri.sh`）、GNU Awk 5.1.0、CSUB（`-q c01`、`-n 64`）。

<details>
<summary>English</summary>
**Purpose**: MUMmer4 + SyRI comparison of MG20T2T (reference) and GifuT2T (query) — per-type variant BEDs, event and base-pair counts, and PAV from the complement of the aligned intervals; `syri.sh` is a separate, earlier minimap2 route over three genome pairs ending in one combined plotsr plot.
**Entry point**: `GifuT2T_vs_MG20T2T_variation_finder.sh`
**Inputs**: `Lotus_MG20T2T_v1.0.fasta` and `Lotus_GifuT2T_v1.0.fasta` (bare filenames); the PAV step also needs `GifuT2T_MG20T2T.filtered.coords` from the step before; `syri.sh` uses `Gifu_v1.fa`, `Gifu_T2T.fa`, `MG20_T2T.fa` and `MG20_v1.fa` plus `genomes.txt` under ``{PROJ_LOTUS_ZC}`/comp/reference/`.
**Outputs**: `GifuT2T_MG20T2T.delta`/`.filtered.delta`, dnadiff reports, `.filtered.coords`, `syri.out`, `MG20_<type>.bed` and `Gifu_<type>.bed` (`type` = SNP, INS, DEL, INV, TRANS, INVTR, DUP, INVDP, PAV), the count table and `variation.out`; the PAV step writes `<asm>.aligned.bed`, `<asm>.PAV.bed` and `PAV.summary`; `syri.sh` writes three BAMs with `.bai` and `syri_samechromosome_plot.pdf`.
**Run**: `bash GifuT2T_vs_MG20T2T_variation_finder.sh`, then `bash GifuT2T_vs_MG20T2T_get_PAV_from_coords.sh`; `bash syri.sh` is a separate route (each comparison must run in its own working directory, otherwise `syri.out`, the BEDs and temporary files overwrite one another).
**Tools**: MUMmer4 4.0.1 (`nucmer -c 1000 --mum --maxgap=1000`, `delta-filter -i 90 -l 1000 -q`, `show-coords -THrd`, `dnadiff -d -p`), SyRI 1.7.1, plotsr 1.1.1, samtools 1.21, bedtools 2.31.1 (`merge -d 50`, `complement -g`), minimap2 (`-ax asm20`, `syri.sh` only), GNU Awk 5.1.0, CSUB (`-q c01`, `-n 64`).

</details>


