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

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>bedtools.sh</code> · Shell · 3 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/02_annotation/05_geta/bedtools.sh)

```bash
bedtools intersect -a gene.bed -b gene.bed -wa -wb \
| awk -F'\t' '$4 < $8 {print $0}' \
> overlaps.tsv
```

</details>

<details open>
<summary><code>geta.sh</code> · Shell · 21 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/02_annotation/05_geta/geta.sh)

```bash
#!/bin/bash
#CSUB -J geta
#CSUB -q c01
#CSUB -o geta.out
#CSUB -e geta.error
#CSUB -n 176
#CSUB -R span[ptile=88]

geta.pl --genome {PROJ_LOTUS_ZC}/Gifu/00.assemble_result/Lotus_GifuT2T_v1.0.fasta \
--RM_species_Dfam Viridiplantae --RM_species_RepBase Viridiplantae \
--pe1 {PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_flower_rep1-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_flower_rep2-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_flower_rep3-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_leaf_rep1-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_leaf_rep2-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_leaf_rep3-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_nodule_dpi10_rep1-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_nodule_dpi10_rep2-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_nodule_dpi10_rep3-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_nodule_dpi21_rep1-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_nodule_dpi21_rep2-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_nodule_dpi21_rep3-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_nodule_dpi45_rep1-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_nodule_dpi45_rep2-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_nodule_dpi45_rep3-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_pod_rep1-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_pod_rep2-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_pod_rep3-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_dpi02_rep1-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_dpi02_rep2-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_dpi02_rep3-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_dpi10_rep1-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_dpi10_rep2-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_dpi10_rep3-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_dpi21_rep1-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_dpi21_rep2-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_dpi21_rep3-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_dpi45_rep1-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_dpi45_rep2-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_dpi45_rep3-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_hpi48_rep1-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_hpi48_rep2-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_hpi48_rep3-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_uni_rep1-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_uni_rep2-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_uni_rep3-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_stem_rep1-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_stem_rep2-R1.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_stem_rep3-R1.fastq.gz \
--pe2 {PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_flower_rep1-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_flower_rep2-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_flower_rep3-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_leaf_rep1-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_leaf_rep2-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_leaf_rep3-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_nodule_dpi10_rep1-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_nodule_dpi10_rep2-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_nodule_dpi10_rep3-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_nodule_dpi21_rep1-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_nodule_dpi21_rep2-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_nodule_dpi21_rep3-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_nodule_dpi45_rep1-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_nodule_dpi45_rep2-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_nodule_dpi45_rep3-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_pod_rep1-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_pod_rep2-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_pod_rep3-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_dpi02_rep1-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_dpi02_rep2-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_dpi02_rep3-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_dpi10_rep1-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_dpi10_rep2-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_dpi10_rep3-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_dpi21_rep1-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_dpi21_rep2-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_dpi21_rep3-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_dpi45_rep1-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_dpi45_rep2-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_dpi45_rep3-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_hpi48_rep1-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_hpi48_rep2-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_hpi48_rep3-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_uni_rep1-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_uni_rep2-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_root_uni_rep3-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_stem_rep1-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_stem_rep2-R2.fastq.gz,{PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/Gifu_stem_rep3-R2.fastq.gz \
--protein {PROJ_LOTUS_ZC}/Gifu/09.annotation/06.braker/douke_pep.fasta \
--config {SOFTWARE_ZC}/geta-2.7.1/conf_for_big_genome.txt \
--augustus_species Gifu_T2T \
--HMM_db {SOFTWARE_ZC}/Pfam/Pfam-A,{SOFTWARE_ZC}/Pfam/Pfam-B \
--BUSCO_lineage_dataset {SOFTWARE_ZC}/busco_downloads/lineages/embryophyta_odb10 \
--out_prefix Gifu_T2T \
--gene_prefix Gifu_T2T \
--cpu 176 \
--max_used_read_num 2000000000
```

</details>
