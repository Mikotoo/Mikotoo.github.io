---
title: '基因模型过滤'
sidebar:
  label: '基因模型过滤'
  order: 6
---


## 中文

**作用**：对 GETA 基因模型做跨工具过滤——合并 RNA-seq 文库，用 BLASTP 与 Pfam 筛查预测蛋白，并估计转录丰度，据此保留有表达支持的基因模型。
**入口**：`cat.sh`、`run_BLASTP.sh`、`run_PfamScan.sh`、`run_RSEM.sh`。
**输入**：`.../data/Gifu/RNA/clean/1` 与 `/2` 中列出的各文库 clean reads（由 `cat.sh` 合并）；Gifu T2T 预测蛋白 `Gifu_T2T.protein.fasta`（来自步骤 05）；BLASTP 命令列表 `cmd.list`；RSEM 所需的 `reference.fasta`、`reference.fasta.fai`、`gene.map`。
**输出**：`Gifu_RNA_all_R1.fastq.gz`、`Gifu_RNA_all_R2.fastq.gz`；`blast*.out.tmp` → `BLASTP.OUT.TMP`；`pfam.out`；`rsem_outdir/`（含 `bowtie.bam`）、`sorted.bam`、`gene.bed`、`transcriptome_chromSizes.txt`、`cov.bed`。
**运行**：`bash cat.sh` → `bash run_BLASTP.sh` → `bash run_PfamScan.sh` → `bash run_RSEM.sh`（`cat.sh` 必须先跑，`run_RSEM.sh` 消费它生成的合并 fastq）。
**工具**：ParaFly `-c cmd.list -CPU 80`；`pfam_scan.pl -fasta Gifu_T2T.protein.fasta -dir <PfamScan/data> -cpu 88 -outfile pfam.out`；`align_and_estimate_abundance.pl --est_method RSEM --aln_method bowtie --prep_reference --thread_count 64` 加 bowtie；`samtools sort -@ 64`；`bedtools coverage -a gene.bed -b sorted.bam -sorted -g transcriptome_chromSizes.txt`；BLASTP 参数写在 `cmd.list` 里。

<details>
<summary>English</summary>

**Purpose**: Cross-tool filtering of the GETA gene models — pool the RNA-seq libraries, screen the predicted proteins by BLASTP and against Pfam, and estimate transcript abundance so that gene models are kept on expression support.
**Entry point**: `cat.sh`, `run_BLASTP.sh`, `run_PfamScan.sh`, `run_RSEM.sh`.
**Inputs**: the clean reads of the libraries listed in `.../data/Gifu/RNA/clean/1` and `/2` (pooled by `cat.sh`); Gifu T2T predicted proteins `Gifu_T2T.protein.fasta` (from step 05); `cmd.list`, the BLASTP command list; `reference.fasta`, `reference.fasta.fai` and `gene.map`, required by RSEM.
**Outputs**: `Gifu_RNA_all_R1.fastq.gz` and `Gifu_RNA_all_R2.fastq.gz`; `blast*.out.tmp` → `BLASTP.OUT.TMP`; `pfam.out`; `rsem_outdir/` (containing `bowtie.bam`), `sorted.bam`, `gene.bed`, `transcriptome_chromSizes.txt`, `cov.bed`.
**Run**: `bash cat.sh` → `bash run_BLASTP.sh` → `bash run_PfamScan.sh` → `bash run_RSEM.sh` (`cat.sh` must run first, since `run_RSEM.sh` consumes the pooled fastq files it writes).
**Tools**: ParaFly `-c cmd.list -CPU 80`; `pfam_scan.pl -fasta Gifu_T2T.protein.fasta -dir <PfamScan/data> -cpu 88 -outfile pfam.out`; `align_and_estimate_abundance.pl --est_method RSEM --aln_method bowtie --prep_reference --thread_count 64` with bowtie; `samtools sort -@ 64`; `bedtools coverage -a gene.bed -b sorted.bam -sorted -g transcriptome_chromSizes.txt`; BLASTP parameters live in `cmd.list`.

</details>

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>cat.sh</code> · Shell · 10 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/02_annotation/06_gene_model_filter/cat.sh)

```bash
#!/bin/bash
#CSUB -J cat
#CSUB -q c01
#CSUB -o cat.out
#CSUB -e cat.error
#CSUB -n 64
#CSUB -R span[host=1]

cat {PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/2|while read id;do cat ${id}>> Gifu_RNA_all_R2.fastq.gz; done
cat {PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/1|while read id;do cat ${id}>> Gifu_RNA_all_R1.fastq.gz; done
```

</details>

<details open>
<summary><code>run_BLASTP.sh</code> · Shell · 11 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/02_annotation/06_gene_model_filter/run_BLASTP.sh)

```bash
#!/bin/bash
#CSUB -J cat
#CSUB -q c01
#CSUB -o cat.out
#CSUB -e cat.error
#CSUB -n 88
#CSUB -R span[host=1]

ParaFly -c cmd.list -CPU 80
cat blast*.out.tmp > BLASTP.OUT.TMP
rm blast*.out.tmp
```

</details>

<details open>
<summary><code>run_PfamScan.sh</code> · Shell · 9 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/02_annotation/06_gene_model_filter/run_PfamScan.sh)

```bash
#!/bin/bash
#CSUB -J run_pfam
#CSUB -q c01
#CSUB -o run_pfam.out
#CSUB -e run_pfam.error
#CSUB -n 88
#CSUB -R span[host=1]

perl {SOFTWARE_ZC}/PfamScan/pfam_scan.pl -fasta  {PROJ_LOTUS_ZC}/Gifu/09.annotation/07.geta/Gifu_T2T.protein.fasta -dir {SOFTWARE_ZC}/PfamScan/data -cpu 88 -outfile pfam.out
```

</details>

<details open>
<summary><code>run_RSEM.sh</code> · Shell · 15 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/02_annotation/06_gene_model_filter/run_RSEM.sh)

```bash
#!/bin/bash
#CSUB -J rsem
#CSUB -q c01
#CSUB -o rsem.out
#CSUB -e rsem.error
#CSUB -n 64
#CSUB -R span[host=1]


mkdir rsem_outdir
{CONDA_ZC}/envs/annotation/bin/align_and_estimate_abundance.pl --transcripts reference.fasta --seqType fq --left {PROJ_LOTUS_ZC}/Gifu/09.annotation/08.geta_filter/Gifu_RNA_all_R1.fastq.gz --right {PROJ_LOTUS_ZC}/Gifu/09.annotation/08.geta_filter/Gifu_RNA_all_R2.fastq.gz --est_method RSEM --aln_method bowtie --gene_trans_map gene.map --prep_reference --output_dir rsem_outdir --thread_count 64
samtools sort -@ 64 -o sorted.bam rsem_outdir/bowtie.bam
awk '{print $1"\t1\t"$2}' reference.fasta.fai > gene.bed
awk '{print $1"\t"$2}' reference.fasta.fai > transcriptome_chromSizes.txt
bedtools coverage -a gene.bed -b sorted.bam -sorted -g transcriptome_chromSizes.txt > cov.bed
```

</details>
