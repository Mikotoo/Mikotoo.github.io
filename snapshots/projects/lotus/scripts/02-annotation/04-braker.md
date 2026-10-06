---
title: '基因预测（BRAKER）'
sidebar:
  label: '基因预测（BRAKER）'
  order: 4
---


## 中文

**作用**：在 EDTA 屏蔽后的 Gifu 基因组上用蛋白、Iso-Seq 与 RNA-seq 证据跑 BRAKER3，再按最小长度、N 端甲硫氨酸、无内部终止密码子过滤蛋白及其对应的 CDS 与 GFF3。
**入口**：`Gifu_braker.sh`（BRAKER3）；`makeblastdb.sh`、`filter.sh`、`bedtools.sh`、`filter_gff_overlaps.py` 建库、过滤与去冗余。
**输入**：屏蔽基因组 ``{PROJ_LOTUS_ZC}`/Gifu/09.annotation/03.edta/Gifu_v1.0.fasta.mod.MAKER.masked`（来自 `01_repeat_edta`）；蛋白证据 `.../06.braker/douke_pep.fasta`；Iso-Seq BAM `.../05.isoseq/Gifu_isoforms.sorted.bam`；RNA-seq BAM 在 `.../04.rnasq/01.align/`（逗号分隔的 `--bam` 列表在 `Gifu_braker.sh` 第 11 行）；`filter.sh` 另读 `braker.aa`、`braker.codingseq`、`braker.gff3`，`bedtools.sh` 读 `gene_filter.bed`。
**输出**：`.../06.braker` 内的 `braker.aa`、`braker.codingseq`、`braker.gff3`；`Gifu_filtered.pep.fa`、`Gifu_filtered.cds.fa`、`Gifu_filtered.gff3`（`02_rnaseq_align` 用的是 `Gifu_filtered.gtf`，需自行转换）；`overlaps.tsv`；`filter_gff_overlaps.py` 写其第二个位置参数。
**运行**：`bash makeblastdb.sh` → `bash Gifu_braker.sh` → `bash filter.sh` → `python3 filter_gff_overlaps.py <in.gff3> <out.gff3>`；`filter.sh` 调 `python AnnotationFilter.py -l 20 -i braker.aa -o Gifu_filtered.pep.fa -c braker.codingseq -C Gifu_filtered.cds.fa -g braker.gff3 -G Gifu_filtered.gff3`；`bedtools.sh` 跑 `bedtools intersect -a gene_filter.bed -b gene_filter.bed -wa -wb | awk -F'\t' '$4 < $9' > overlaps.tsv`。
**工具**：BRAKER3 `braker.pl --species=Gifu_1 --threads 88 --gff3`（`singularity exec braker3_sandbox`，镜像是本地镜像名而非路径）；AUGUSTUS `--AUGUSTUS_CONFIG_PATH=`{SOFTWARE_ZC}`/Augustus-master/config/`；BLAST+ `makeblastdb -in douke_pep.fasta -dbtype prot -out douke_pep.fasta`；bedtools `intersect -wa -wb`；`AnnotationFilter.py -l 20`；`filter_gff_overlaps.py --strand-aware`。

<details>
<summary>English</summary>
**Purpose**: Run BRAKER3 on the EDTA-masked Gifu genome with protein, Iso-Seq and RNA-seq evidence, then filter the predicted proteins (minimum length, N-terminal methionine, no internal stop codon) together with their matching CDS and GFF3 features.
**Entry point**: `Gifu_braker.sh` (BRAKER3); `makeblastdb.sh`, `filter.sh`, `bedtools.sh` and `filter_gff_overlaps.py` build the BLAST database, filter and reduce redundancy.
**Inputs**: masked genome ``{PROJ_LOTUS_ZC}`/Gifu/09.annotation/03.edta/Gifu_v1.0.fasta.mod.MAKER.masked` (from `01_repeat_edta`); protein evidence `.../06.braker/douke_pep.fasta`; Iso-Seq BAM `.../05.isoseq/Gifu_isoforms.sorted.bam`; RNA-seq BAMs in `.../04.rnasq/01.align/` (the comma-separated `--bam` list is on line 11 of `Gifu_braker.sh`); `filter.sh` also reads `braker.aa`, `braker.codingseq`, `braker.gff3`, and `bedtools.sh` reads `gene_filter.bed`.
**Outputs**: `braker.aa`, `braker.codingseq`, `braker.gff3` in `.../06.braker`; `Gifu_filtered.pep.fa`, `Gifu_filtered.cds.fa`, `Gifu_filtered.gff3` (note that `02_rnaseq_align` uses `Gifu_filtered.gtf`); `overlaps.tsv`; `filter_gff_overlaps.py` writes its second positional argument.
**Run**: `bash makeblastdb.sh` → `bash Gifu_braker.sh` → `bash filter.sh` → `python3 filter_gff_overlaps.py <in.gff3> <out.gff3>`; `filter.sh` calls `python AnnotationFilter.py -l 20 -i braker.aa -o Gifu_filtered.pep.fa -c braker.codingseq -C Gifu_filtered.cds.fa -g braker.gff3 -G Gifu_filtered.gff3`; `bedtools.sh` runs `bedtools intersect -a gene_filter.bed -b gene_filter.bed -wa -wb | awk -F'\t' '$4 < $9' > overlaps.tsv`.
**Tools**: BRAKER3 `braker.pl --species=Gifu_1 --threads 88 --gff3` (via `singularity exec braker3_sandbox`, a local image name rather than a path); AUGUSTUS `--AUGUSTUS_CONFIG_PATH=`{SOFTWARE_ZC}`/Augustus-master/config/`; BLAST+ `makeblastdb -in douke_pep.fasta -dbtype prot -out douke_pep.fasta`; bedtools `intersect -wa -wb`; `AnnotationFilter.py -l 20`; `filter_gff_overlaps.py --strand-aware`.

</details>


