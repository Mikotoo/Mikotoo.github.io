---
title: 'hifiasm 组装'
sidebar:
  label: 'hifiasm 组装'
  order: 1
---


## 中文
**作用**：用 hifiasm 对 ONT ultra-long reads 与 PacBio HiFi reads 进行初步 contig 组装，并把 primary contig GFA 转成 FASTA 并建索引。
**入口**：`01.hifiasm.sh`
**输入**：ONT UL reads; HiFi reads
**输出**：`Gifu.sam.bp.p_ctg.gfa`、`Gifu.hifiasm.fasta`、`Gifu.hifiasm.fasta.fai`
**运行**：`bash 01.hifiasm.sh`
**工具**：hifiasm `-t 64 -o Gifu.sam --telo-m CCCTAAA --ul`（UL reads 在前、HiFi FASTQ 在后）；gfatools `gfa2fa`；samtools `faidx`

<details>
<summary>English</summary>

**Purpose**: Assemble the Gifu accession from ONT ultra-long and PacBio HiFi reads with hifiasm, then convert the primary contig GFA to FASTA and index it.
**Entry point**: `01.hifiasm.sh`
**Inputs**: ONT UL reads; HiFi reads
**Outputs**: `Gifu.sam.bp.p_ctg.gfa`, `Gifu.hifiasm.fasta`, `Gifu.hifiasm.fasta.fai`.
**Run**: `bash 01.hifiasm.sh`
**Tools**: hifiasm `-t 64 -o Gifu.sam --telo-m CCCTAAA --ul` (UL reads first, then the HiFi FASTQ); gfatools `gfa2fa`; samtools `faidx`.

</details>

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>01.hifiasm.sh</code> · Shell · 7 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/01_assembly/01_hifiasm/01.hifiasm.sh)

```bash
#1.1 contig组装

hifiasm -t 64 -o Gifu.sam --telo-m CCCTAAA --ul {ontReads}.fq.gz {hifiReads}.fastq.gz

gfatools gfa2fa Gifu.sam.bp.p_ctg.gfa > Gifu.hifiasm.fasta

samtools faidx Gifu.hifiasm.fasta
```

</details>
