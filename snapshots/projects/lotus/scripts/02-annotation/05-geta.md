---
title: '基因模型整合（GETA）'
sidebar:
  label: '基因模型整合（GETA）'
  order: 5
---


## 中文

**作用**：在 Gifu T2T 组装上用 GETA 结合 RNA-seq、蛋白同源、从头预测、重复模型、Pfam 结构域和 BUSCO 谱系做第二套独立基因预测。
**入口**：`geta.sh`；`bedtools.sh` 查 `gene.bed` 自身重叠。
**输入**：基因组 ``{PROJ_LOTUS_ZC}`/Gifu/00.assemble_result/Lotus_GifuT2T_v1.0.fasta`；``{PROJ_LOTUS_ZC}`/data/Gifu/RNA/clean/` 下的 `-R1`/`-R2` `fastq.gz`（`--pe1`/`--pe2`，完整列表在 `geta.sh` 第 11–12 行）；蛋白证据 `.../06.braker/douke_pep.fasta`；配置 ``{SOFTWARE_ZC}`/geta-2.7.1/conf_for_big_genome.txt`；HMM 库 ``{SOFTWARE_ZC}`/Pfam/Pfam-A`、`Pfam-B`；BUSCO 谱系 ``{SOFTWARE_ZC}`/busco_downloads/lineages/embryophyta_odb10`；`bedtools.sh` 读 `gene.bed`。
**输出**：GETA 预测以 `--out_prefix Gifu_T2T` 为前缀、基因 ID 用 `--gene_prefix Gifu_T2T`（脚本未逐个列出产出的文件名）；`bedtools.sh` 的 `overlaps.tsv`（`gene.bed` 自身相交，保留 `$4 < $8`）。
**运行**：`bash geta.sh`；`bash bedtools.sh`（用裸相对名 `gene.bed`/`overlaps.tsv`，须在存放它们的目录中运行）。
**工具**：GETA `geta.pl --cpu 176 --max_used_read_num 2000000000 --RM_species_Dfam Viridiplantae --RM_species_RepBase Viridiplantae --augustus_species Gifu_T2T --HMM_db Pfam-A,Pfam-B --config conf_for_big_genome.txt`；bedtools `intersect -a gene.bed -b gene.bed -wa -wb`。

<details>
<summary>English</summary>
**Purpose**: Run GETA on the Gifu T2T assembly with RNA-seq, protein homology, ab initio prediction, repeat models, Pfam domains and a BUSCO lineage, giving a second, independent gene prediction.
**Entry point**: `geta.sh`; `bedtools.sh` checks `gene.bed` self-overlaps.
**Inputs**: genome ``{PROJ_LOTUS_ZC}`/Gifu/00.assemble_result/Lotus_GifuT2T_v1.0.fasta`; paired `-R1`/`-R2` `fastq.gz` files under ``{PROJ_LOTUS_ZC}`/data/Gifu/RNA/clean/` (`--pe1`/`--pe2`, full lists on lines 11–12 of `geta.sh`); protein evidence `.../06.braker/douke_pep.fasta`; config ``{SOFTWARE_ZC}`/geta-2.7.1/conf_for_big_genome.txt`; HMM libraries ``{SOFTWARE_ZC}`/Pfam/Pfam-A` and `Pfam-B`; BUSCO lineage ``{SOFTWARE_ZC}`/busco_downloads/lineages/embryophyta_odb10`; `bedtools.sh` reads `gene.bed`.
**Outputs**: the GETA prediction prefixed `--out_prefix Gifu_T2T` with gene IDs under `--gene_prefix Gifu_T2T` (the script does not enumerate the output files); `overlaps.tsv` from `bedtools.sh` (self-intersect of `gene.bed`, keeping `$4 < $8`).
**Run**: `bash geta.sh`; `bash bedtools.sh` (it uses the bare relative names `gene.bed` and `overlaps.tsv`, so it must run in the directory holding them).
**Tools**: GETA `geta.pl --cpu 176 --max_used_read_num 2000000000 --RM_species_Dfam Viridiplantae --RM_species_RepBase Viridiplantae --augustus_species Gifu_T2T --HMM_db Pfam-A,Pfam-B --config conf_for_big_genome.txt`; bedtools `intersect -a gene.bed -b gene.bed -wa -wb`.

</details>


