---
title: 'hifiasm 组装'
sidebar:
  label: 'hifiasm 组装'
  order: 1
---


## 中文
**作用**：用 hifiasm 以 ONT ultra-long 与 PacBio HiFi reads 组装 Gifu，并把 primary contig GFA 转成 FASTA 并建索引。
**入口**：`01.hifiasm.sh`
**输入**：ONT UL reads ``{PROJ_LOTUS_ZC}`/data/Gifu/ont/Gifu-jing-ye.pass.ul.fq.gz`；HiFi reads ``{PROJ_LOTUS_ZC}`/data/Gifu/hifi/Gifu_hifi.fastq.gz`。
**输出**：`Gifu.sam.bp.p_ctg.gfa`、`Gifu.hifiasm.fasta`、`Gifu.hifiasm.fasta.fai`。
**运行**：`bash 01.hifiasm.sh`
**工具**：hifiasm `-t 88 -o Gifu.sam --telo-m CCCTAAA --ul`（UL reads 在前、HiFi FASTQ 在后）；gfatools `gfa2fa`；samtools `faidx`。
**注意**：`01.hifiasm.sh` 中 hifiasm 命令处于注释状态，生效的只有 `gfatools gfa2fa` 与 `samtools faidx`，运行前须已存在 GFA。

<details>
<summary>English</summary>
**Purpose**: Assemble the Gifu accession from ONT ultra-long and PacBio HiFi reads with hifiasm, then convert the primary contig GFA to FASTA and index it.
**Entry point**: `01.hifiasm.sh`
**Inputs**: ONT UL reads ``{PROJ_LOTUS_ZC}`/data/Gifu/ont/Gifu-jing-ye.pass.ul.fq.gz`; HiFi reads ``{PROJ_LOTUS_ZC}`/data/Gifu/hifi/Gifu_hifi.fastq.gz`.
**Outputs**: `Gifu.sam.bp.p_ctg.gfa`, `Gifu.hifiasm.fasta`, `Gifu.hifiasm.fasta.fai`.
**Run**: `bash 01.hifiasm.sh`
**Tools**: hifiasm `-t 88 -o Gifu.sam --telo-m CCCTAAA --ul` (UL reads first, then the HiFi FASTQ); gfatools `gfa2fa`; samtools `faidx`.
**Notes**: The hifiasm command in `01.hifiasm.sh` is commented out; only `gfatools gfa2fa` and `samtools faidx` are active, so a GFA must already exist.

</details>


