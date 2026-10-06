---
title: '功能注释'
sidebar:
  label: '功能注释'
  order: 8
---


## 中文

**作用**：为最终 Gifu T2T 基因模型赋予功能描述：`eggnog.sh` 用 eggNOG-mapper 做主功能注释，同目录的 `run_BLASTP.sh`、`run_PfamScan.sh`、`run_RSEM.sh` 提供同源、结构域与表达证据。
**入口**：`eggnog.sh`；`run_BLASTP.sh`、`run_PfamScan.sh`、`run_RSEM.sh` 各自独立运行。
**输入**：Gifu T2T 预测蛋白 `00_assembly_result/GifuT2T/Lotus_GifuT2T.pep.fa`（`eggnog.sh` 使用）；GETA 输出树下的 `Gifu_T2T.protein.fasta`（`run_PfamScan.sh` 使用）；`cmd.list`（ParaFly 的 BLASTP 命令列表）；`reference.fasta`、`reference.fasta.fai`、`gene.map` 与合并的 `Gifu_RNA_all_R1/R2.fastq.gz`（`run_RSEM.sh` 使用）。
**输出**：`Lotus_Gifu.emapper.annotations` 及其余 `Lotus_Gifu.*` eggNOG-mapper 结果；`blast*.out.tmp` → `BLASTP.OUT.TMP`；`pfam.out`；`rsem_outdir/`（含 `bowtie.bam`）、`sorted.bam`、`gene.bed`、`transcriptome_chromSizes.txt`、`cov.bed`。
**运行**：`bash eggnog.sh`（`python `{SOFTWARE_ZC}`/eggnog-mapper-2.1.13/emapper.py -m hmmer -d Eukaryota -i {PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Lotus_GifuT2T.pep.fa --output Lotus_Gifu --cpu 64`）；`bash run_BLASTP.sh`；`bash run_PfamScan.sh`；`bash run_RSEM.sh`（后三者没有作业头，需要 PBS 环境变量 `$PBS_O_WORKDIR`）。
**工具**：eggNOG-mapper 2.1.13 `emapper.py -m hmmer -d Eukaryota --output Lotus_Gifu --cpu 64`；ParaFly `-c cmd.list -CPU 0`；`pfam_scan.pl -fasta Gifu_T2T.protein.fasta -dir `{SOFTWARE_ZC}`/Pfam -cpu 1 -outfile pfam.out`；`align_and_estimate_abundance.pl --est_method RSEM --aln_method bowtie --prep_reference --thread_count 1` 加 bowtie；`samtools sort -@ 1`；`bedtools coverage -a gene.bed -b sorted.bam -sorted -g transcriptome_chromSizes.txt`。

<details>
<summary>English</summary>
**Purpose**: Assign functional descriptions to the final Gifu T2T gene models: `eggnog.sh` runs the main functional annotation with eggNOG-mapper, and the copies of `run_BLASTP.sh`, `run_PfamScan.sh` and `run_RSEM.sh` in this directory provide the homology, domain and expression evidence.
**Entry point**: `eggnog.sh`; `run_BLASTP.sh`, `run_PfamScan.sh` and `run_RSEM.sh` each run on their own.
**Inputs**: Gifu T2T predicted proteins `00_assembly_result/GifuT2T/Lotus_GifuT2T.pep.fa` (used by `eggnog.sh`); `Gifu_T2T.protein.fasta` under the GETA output tree (used by `run_PfamScan.sh`); `cmd.list`, the BLASTP command list for ParaFly; `reference.fasta`, `reference.fasta.fai`, `gene.map` and the pooled `Gifu_RNA_all_R1/R2.fastq.gz` (used by `run_RSEM.sh`).
**Outputs**: `Lotus_Gifu.emapper.annotations` and the other `Lotus_Gifu.*` eggNOG-mapper outputs; `blast*.out.tmp` → `BLASTP.OUT.TMP`; `pfam.out`; `rsem_outdir/` (containing `bowtie.bam`), `sorted.bam`, `gene.bed`, `transcriptome_chromSizes.txt`, `cov.bed`.
**Run**: `bash eggnog.sh` (`python `{SOFTWARE_ZC}`/eggnog-mapper-2.1.13/emapper.py -m hmmer -d Eukaryota -i {PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Lotus_GifuT2T.pep.fa --output Lotus_Gifu --cpu 64`); `bash run_BLASTP.sh`; `bash run_PfamScan.sh`; `bash run_RSEM.sh` (the last three carry no job header and expect the PBS variable `$PBS_O_WORKDIR`).
**Tools**: eggNOG-mapper 2.1.13 `emapper.py -m hmmer -d Eukaryota --output Lotus_Gifu --cpu 64`; ParaFly `-c cmd.list -CPU 0`; `pfam_scan.pl -fasta Gifu_T2T.protein.fasta -dir `{SOFTWARE_ZC}`/Pfam -cpu 1 -outfile pfam.out`; `align_and_estimate_abundance.pl --est_method RSEM --aln_method bowtie --prep_reference --thread_count 1` with bowtie; `samtools sort -@ 1`; `bedtools coverage -a gene.bed -b sorted.bam -sorted -g transcriptome_chromSizes.txt`.

</details>


