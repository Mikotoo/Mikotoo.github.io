---
title: 'Assembly qc'
sidebar:
  label: 'Assembly qc'
  order: 9
---


## 中文
**作用**：用回贴深度、BUSCO 与 Merqury 三块独立评估最终组装的完整度与碱基准确性。
**入口**：`01_mapping/{01_ont,02_hifi,03_wgs,04_pan_wgs}/`、`02_busco/`、`03_merqury/{01_hifi,02_wgs,03_pan_wgs}/` 下的脚本
**输入**：组装 ``{PROJ_LOTUS_ZC}`/Gifu/00.assemble_result/Gifu_v1.0.fasta`；reads ONT `Gifu-jing-ye.pass.ul.fa`、HiFi `Gifu_hifi.fastq.gz`、WGS `Gifu_1-R{1,2}.fastq.gz`；窗口 `Gifu_{100k,10k}.bed`；HiFi k-mer 库 `ccs.merylDB`；BUSCO lineage `embryophyta_odb10`。
**输出**：`gifu_{ont,hifi,wgs}_{100k,10k}_coverage.txt`；BUSCO `Gifu_lotus/`；Merqury `Gifu_hifimer`、`Gifu_wgsmer`、`Gifu_<sample>_wgsmer`；pan-WGS `{Gifu,MG20}_Clean_R{1,2}.fq` 与 `MG20_subsample_{summary,scientific_name}.txt`。
**运行**：进入相应子目录后 `bash <script>`
**工具**：minimap2 `-ax map-ont`/`-ax map-pb`、bwa `mem`；samtools `view -b -q 30`、`sort`、`flagstat`；bedtools `coverage -mean`（100 kb 与 10 kb 窗口）；BBMap `filterbyname.sh`；BLAST+ `blastn -outfmt 5`；BUSCO `-m genome --offline --cpu 64`；Merqury 与 meryl `count k=21`。
**注意**：`ont.sh`、`ccs.sh`、`bwa.sh` 与 `04_pan_wgs/samtools.sh` 的比对命令均被注释，脚本只对现成的 `ont_q30.bam`、`ccs_q30.bam`、`wgs_q30.bam` 与 `*_unmapped.literal.ids` 重跑下游覆盖度与 reads 过滤，须先自备这些文件。

<details>
<summary>English</summary>

**Purpose**: Assess the finished assembly in three independent blocks: read-remapping depth, BUSCO gene-space completeness and Merqury k-mer QV.
**Entry point**: scripts under `01_mapping/{01_ont,02_hifi,03_wgs,04_pan_wgs}/`, `02_busco/` and `03_merqury/{01_hifi,02_wgs,03_pan_wgs}/`
**Inputs**: assembly ``{PROJ_LOTUS_ZC}`/Gifu/00.assemble_result/Gifu_v1.0.fasta`; reads ONT `Gifu-jing-ye.pass.ul.fa`, HiFi `Gifu_hifi.fastq.gz`, WGS `Gifu_1-R{1,2}.fastq.gz`; windows `Gifu_{100k,10k}.bed`; HiFi k-mer database `ccs.merylDB`; BUSCO lineage `embryophyta_odb10`.
**Outputs**: `gifu_{ont,hifi,wgs}_{100k,10k}_coverage.txt`; BUSCO `Gifu_lotus/`; Merqury `Gifu_hifimer`, `Gifu_wgsmer`, `Gifu_<sample>_wgsmer`; pan-WGS `{Gifu,MG20}_Clean_R{1,2}.fq` and `MG20_subsample_{summary,scientific_name}.txt`.
**Run**: `cd` into the relevant sub-directory and run `bash <script>`
**Tools**: minimap2 `-ax map-ont`/`-ax map-pb`, bwa `mem`; samtools `view -b -q 30`, `sort`, `flagstat`; bedtools `coverage -mean` (100 kb and 10 kb windows); BBMap `filterbyname.sh`; BLAST+ `blastn -outfmt 5`; BUSCO `-m genome --offline --cpu 64`; Merqury with meryl `count k=21`.
**Notes**: The alignment commands in `ont.sh`, `ccs.sh`, `bwa.sh` and `04_pan_wgs/samtools.sh` are all commented out; the scripts only re-run the downstream coverage and read filtering on supplied `ont_q30.bam`, `ccs_q30.bam`, `wgs_q30.bam` and `*_unmapped.literal.ids`, so those files must be provided.

</details>

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>01_mapping/01_ont/ont.sh</code> · Shell · 10 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/01_assembly/09_assembly_qc/01_mapping/01_ont/ont.sh)

```bash
minimap2 -ax map-ont {PROJ_LOTUS_ZC}/Gifu/00.assemble_result/Gifu_v1.0.fasta {PROJ_LOTUS_ZC}/data/Gifu/ont/Gifu-jing-ye.pass.ul.fa -o ont.sam
samtools sort -@ 88 -o ont.bam -T ont_tmp.ali ont.sam
samtools flagstat ont.sam > ont.flagstat

samtools view -h -F 2304 ont.bam -o primary_only.bam

samtools view -b -q 30 -o ont_q30.bam ont.bam

bedtools coverage -a {PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Gifu_100k.bed -b ont_q30.bam -mean > gifu_ont_100k_coverage.txt
bedtools coverage -a {PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Gifu_10k.bed -b ont_q30.bam -mean > gifu_ont_10k_coverage.txt
```

</details>

<details open>
<summary><code>01_mapping/02_hifi/ccs.sh</code> · Shell · 19 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/01_assembly/09_assembly_qc/01_mapping/02_hifi/ccs.sh)

```bash
#!/bin/bash
#CSUB -J ccs
#CSUB -q c01
#CSUB -o ccs.out
#CSUB -e ccs.error
#CSUB -n 64
#CSUB -R span[hosts=1]

#minimap2 -ax map-pb {PROJ_LOTUS_ZC}/Gifu/00.assemble_result/Gifu_v1.0.fasta {PROJ_LOTUS_ZC}/data/Gifu/hifi/Gifu_hifi.fastq.gz -o ccs.sam
#samtools flagstat ccs.sam > ccs.flagstat
#samtools sort -@ 64 -m 10G -o ccs.bam -T ccs_tmp.ali ccs.sam
#samtools flagstat ccs.bam > ccs.flagstat
#samtools index ccs.bam

#samtools view -h -F 2304 ccs.bam -o primary_only.bam
#samtools view -b -q 30 -o ccs_q30.bam ccs.bam

#bedtools coverage -a {PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Gifu_100k.bed -b ccs_q30.bam -mean > gifu_hifi_100k_coverage.txt
#bedtools coverage -a {PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Gifu_10k.bed -b ccs_q30.bam -mean > gifu_hifi_10k_coverage.txt
```

</details>

<details open>
<summary><code>01_mapping/03_wgs/bwa.sh</code> · Shell · 21 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/01_assembly/09_assembly_qc/01_mapping/03_wgs/bwa.sh)

```bash
#!/bin/bash
#CSUB -J bwa
#CSUB -q c01
#CSUB -o bwa.out
#CSUB -e bwa.error
#CSUB -n 88
#CSUB -R span[hosts=1]

#bwa index -p Gifu_v1.0.fasta -a bwtsw Gifu_v1.0.fasta
#bwa mem Gifu_v1.0.fasta {PROJ_LOTUS_ZC}/data/001-NC_raw_data/clean/Gifu_1-R1.fastq.gz {PROJ_LOTUS_ZC}/data/001-NC_raw_data/clean/Gifu_1-R2.fastq.gz > bwa_1.sam
#samtools flagstat bwa_1.sam > bwa_1.flagstat

#bwa mem Gifu_v1.0.fasta {PROJ_LOTUS_ZC}/data/001-NC_raw_data/clean/Gifu_2-R1.fastq.gz {PROJ_LOTUS_ZC}/data/001-NC_raw_data/clean/Gifu_2-R2.fastq.gz > bwa_2.sam
#samtools flagstat bwa_2.sam > bwa_2.flagstat
#samtools sort -@ 88 -o bwa.bam -T bwa_tmp.ali bwa.sam

samtools view -b -q 30 -o wgs_q30.bam {PROJ_04_LOTUS_GENOME}/06_SNP/02.toCram/picard/GifuT2T_markdup.bam


bedtools coverage -a {PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Gifu_100k.bed -b wgs_q30.bam -mean > gifu_wgs_100k_coverage.txt
bedtools coverage -a {PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Gifu_10k.bed -b wgs_q30.bam -mean > gifu_wgs_10k_coverage.txt
```

</details>

<details open>
<summary><code>02_busco/busco.sh</code> · Shell · 13 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/01_assembly/09_assembly_qc/02_busco/busco.sh)

```bash
#!/bin/bash
#CSUB -J busco
#CSUB -q c01
#CSUB -o busco.out
#CSUB -e busco.error
#CSUB -n 64
#CSUB -R span[hosts=1]

busco -i {PROJ_LOTUS_ZC}/Gifu/00.assemble_result/Gifu_v1.0.fasta \
      -l {SOFTWARE_ZC}/busco_downloads/lineages/embryophyta_odb10 \
      -o Gifu_lotus \
      -m genome \
      --offline --cpu 64
```

</details>

<details open>
<summary><code>03_merqury/01_hifi/merqury.sh</code> · Shell · 9 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/01_assembly/09_assembly_qc/03_merqury/01_hifi/merqury.sh)

```bash
#!/bin/bash
#CSUB -J merqury
#CSUB -q c01
#CSUB -o merqury.out
#CSUB -e merqury.error
#CSUB -n 64
#CSUB -R span[hosts=1]

merqury.sh {PROJ_LOTUS_ZC}/Gifu/07.polish/ccs.merylDB {PROJ_LOTUS_ZC}/Gifu/00.assemble_result/Gifu_v1.0.fasta Gifu_hifimer
```

</details>

<details open>
<summary><code>03_merqury/02_wgs/build_merylDB.sh</code> · Shell · 13 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/01_assembly/09_assembly_qc/03_merqury/02_wgs/build_merylDB.sh)

```bash
#!/bin/bash
#CSUB -J meryl
#CSUB -q c01
#CSUB -o meryl.out
#CSUB -e meryl.error
#CSUB -n 88
#CSUB -R span[hosts=1]

meryl count k=21 {PROJ_LOTUS_ZC}/data/001-NC_raw_data/clean/Gifu_1-R1.fastq.gz output wgs1.merylDB
meryl count k=21 {PROJ_LOTUS_ZC}/data/001-NC_raw_data/clean/Gifu_1-R2.fastq.gz output wgs2.merylDB

meryl  union-sum output wgs.meryl wgs1.merylDB wgs2.merylDB
merqury.sh wgs.meryl {PROJ_LOTUS_ZC}/Gifu/00.assemble_result/Gifu_v1.0.fasta Gifu_wgsmer
```

</details>

<details open>
<summary><code>03_merqury/02_wgs/merqury.sh</code> · Shell · 9 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/01_assembly/09_assembly_qc/03_merqury/02_wgs/merqury.sh)

```bash
#!/bin/bash
#CSUB -J merqury
#CSUB -q c01
#CSUB -o merqury.out
#CSUB -e merqury.error
#CSUB -n 64
#CSUB -R span[hosts=1]

merqury.sh wgs.meryl {PROJ_LOTUS_ZC}/Gifu/00.assemble_result/Gifu_v1.0.fasta Gifu_wgsmer
```

</details>

<details open>
<summary><code>03_merqury/03_pan_wgs/merqury.sh</code> · Shell · 23 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/01_assembly/09_assembly_qc/03_merqury/03_pan_wgs/merqury.sh)

```bash
#!/bin/bash
#CSUB -J merqury
#CSUB -q c01
#CSUB -o merqury.out
#CSUB -e merqury.error
#CSUB -n 64
#CSUB -R span[hosts=1]

cat {PROJ_LOTUS_ZC}/data/001-NC_raw_data/clean/align.conf|while read id;
do
arr=($id)
sample=${arr[0]}
fq1=${arr[1]}
fq2=${arr[2]}

cd {PROJ_LOTUS_ZC}/Gifu/08.quality/03.merqury/03.pan_wgsmer
mkdir $sample && cd $sample
meryl count k=21 $fq1 output ${sample}_R1.merylDB
meryl count k=21 $fq2 output ${sample}_R2.merylDB
meryl union-sum output ${sample}.meryl ${sample}_R1.merylDB ${sample}_R2.merylDB

merqury.sh ${sample}.meryl {PROJ_LOTUS_ZC}/Gifu/00.assemble_result/Gifu_v1.0.fasta Gifu_${sample}_wgsmer
done
```

</details>
