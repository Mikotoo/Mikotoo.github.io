---
title: 'rDNA 拷贝数'
sidebar:
  label: 'rDNA 拷贝数'
  order: 3
---


## 中文

**作用**：从测序深度估计 Gifu 与 MG20 的 45S 与 5S rDNA 拷贝数；把 T2T 装配的 rDNA 阵列屏蔽后追加一个 45S 与一个 5S 代表重复构成折叠参考，WGS/HiFi reads 比对后用 mosdepth 量深度，组件深度除以单拷贝 BUSCO 位点深度中位数即得拷贝数，并与装配来源计数并列。
**入口**：`run_rDNA_CN.sh`
**输入**：装配 `{PROJ_04_LOTUS_GENOME}/00_assembly_result/{GifuT2T/Lotus_GifuT2T_v1.0.fasta,MG20T2T/Lotus_MG20T2T_v1.0.fasta}`；WGS `01_data/<acc>/WGS/<acc>_Clean_R1.fq` 与 `_R2.fq`；HiFi `01_data/<acc>/hifi/<acc>_hifi.fastq.gz`；BUSCO `02_assembly/<acc>/08.quality/02.busco/<acc>_lotus/run_embryophyta_odb10/full_table.tsv`；rDNA 阵列 BED `03_annotation/<acc>/rDNA/rRNA.bed`。
**输出**：`${OUTDIR}` = `{PROJ_04_LOTUS_GENOME}/03_annotation/rDNA_CN`，每材料一个子目录：`02_rdna/`（`rDNA_arrays.bed`、代表重复 `.fa`/`.gff`）、`03_busco/`（`*_BUSCO_singlecopy.bed`）、`04_ref/`（`*_CN_reference.fa` 及索引）、`05_map/`（`*_{WGS,HiFi}_CN.bam`）、`06_depth/`（mosdepth `*.regions.bed.gz`）、`07_cn/`（`*_{WGS,HiFi}_CN_summary.tsv`），另有合并表 `final_rDNA_CN.tsv`。
**运行**：`bash run_rDNA_CN.sh`（`#CSUB -J rDNA_CN`、队列 c01、`-n 64`、`span[hosts=1]`；一次提交内先 Gifu 后 MG20）。
**工具**：bwa `mem -t 50`；minimap2 `-ax map-hifi -t 50`；samtools `faidx`、`sort -@ 16`、`index`、`flagstat`；bedtools `maskfasta`、`intersect -v`；mosdepth `-t 8 --by <bed>`；Barrnap `--kingdom euk` 校验代表重复；Python 3（pandas、numpy）。

<details>
<summary>English</summary>
**Purpose**: Estimate the 45S and 5S rDNA copy number of Gifu and MG20 from read depth; the rDNA arrays of the T2T assembly are masked and one representative 45S and one representative 5S repeat are appended to form a collapsed reference, WGS and HiFi reads are mapped and depth is measured with mosdepth, and each rDNA component's depth divided by the median depth over single-copy BUSCO loci gives the copy number, set beside the assembly-derived counts.
**Entry point**: `run_rDNA_CN.sh`
**Inputs**: the assemblies `{PROJ_04_LOTUS_GENOME}/00_assembly_result/{GifuT2T/Lotus_GifuT2T_v1.0.fasta,MG20T2T/Lotus_MG20T2T_v1.0.fasta}`; WGS `01_data/<acc>/WGS/<acc>_Clean_R1.fq` and `_R2.fq`; HiFi `01_data/<acc>/hifi/<acc>_hifi.fastq.gz`; BUSCO `02_assembly/<acc>/08.quality/02.busco/<acc>_lotus/run_embryophyta_odb10/full_table.tsv`; the rDNA array BED `03_annotation/<acc>/rDNA/rRNA.bed`.
**Outputs**: `${OUTDIR}` = `{PROJ_04_LOTUS_GENOME}/03_annotation/rDNA_CN`, one subdirectory per accession: `02_rdna/` (`rDNA_arrays.bed`, the representative repeats `.fa`/`.gff`), `03_busco/` (`*_BUSCO_singlecopy.bed`), `04_ref/` (`*_CN_reference.fa` and indexes), `05_map/` (`*_{WGS,HiFi}_CN.bam`), `06_depth/` (mosdepth `*.regions.bed.gz`), `07_cn/` (`*_{WGS,HiFi}_CN_summary.tsv`), plus the merged table `final_rDNA_CN.tsv`.
**Run**: `bash run_rDNA_CN.sh` (`#CSUB -J rDNA_CN`, queue c01, `-n 64`, `span[hosts=1]`; Gifu then MG20 in one submission).
**Tools**: bwa `mem -t 50`; minimap2 `-ax map-hifi -t 50`; samtools `faidx`, `sort -@ 16`, `index`, `flagstat`; bedtools `maskfasta`, `intersect -v`; mosdepth `-t 8 --by <bed>`; Barrnap `--kingdom euk` to verify the representative repeats; Python 3 (pandas, numpy).

</details>


