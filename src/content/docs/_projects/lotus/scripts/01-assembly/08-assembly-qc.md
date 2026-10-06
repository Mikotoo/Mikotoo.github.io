---
title: '组装质量评估'
sidebar:
  label: '组装质量评估'
  order: 8
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


