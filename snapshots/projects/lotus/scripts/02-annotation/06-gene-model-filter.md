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


