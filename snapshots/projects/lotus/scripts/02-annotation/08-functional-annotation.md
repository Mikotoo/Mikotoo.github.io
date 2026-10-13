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

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>eggnog.sh</code> · Shell · 9 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/02_annotation/08_functional_annotation/eggnog.sh)

```bash
#!/bin/bash
#CSUB -J eggnog
#CSUB -q c01
#CSUB -o eggnog.out
#CSUB -e eggnog.error
#CSUB -n 64
#CSUB -R span[hosts=1]

python {SOFTWARE_ZC}/eggnog-mapper-2.1.13/emapper.py -m hmmer -d Eukaryota -i {PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Lotus_GifuT2T.pep.fa --output Lotus_Gifu --cpu 64
```

</details>

<details open>
<summary><code>run_BLASTP.sh</code> · Shell · 5 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/02_annotation/08_functional_annotation/run_BLASTP.sh)

```bash
#!/bin/bash
cd $PBS_O_WORKDIR
ParaFly -c cmd.list -CPU 0
cat blast*.out.tmp > BLASTP.OUT.TMP
rm blast*.out.tmp
```

</details>

<details open>
<summary><code>run_PfamScan.sh</code> · Shell · 3 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/02_annotation/08_functional_annotation/run_PfamScan.sh)

```bash
#!/bin/bash
cd $PBS_O_WORKDIR
perl {SOFTWARE_ZC}/PfamScan/pfam_scan.pl -fasta  {PROJ_LOTUS_ZC}/Gifu/09.annotation/07.geta/Gifu_T2T.protein.fasta -dir {SOFTWARE_ZC}/Pfam -cpu 1 -outfile pfam.out
```

</details>

<details open>
<summary><code>run_RSEM.sh</code> · Shell · 7 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/02_annotation/08_functional_annotation/run_RSEM.sh)

```bash
#!/bin/bash
cd $PBS_O_WORKDIR
{CONDA_ZC}/envs/annotation/bin/align_and_estimate_abundance.pl --transcripts reference.fasta --seqType fq --left {PROJ_LOTUS_ZC}/Gifu/09.annotation/08.geta_filter/Gifu_RNA_all_R1.fastq.gz --right {PROJ_LOTUS_ZC}/Gifu/09.annotation/08.geta_filter/Gifu_RNA_all_R1.fastq.gz --est_method RSEM --aln_method bowtie --gene_trans_map gene.map --prep_reference --output_dir rsem_outdir --thread_count 1
samtools sort -@ 1 -o sorted.bam rsem_outdir/bowtie.bam
awk '{print $1"\t1\t"$2}' reference.fasta.fai > gene.bed
awk '{print $1"\t"$2}' reference.fasta.fai > transcriptome_chromSizes.txt
bedtools coverage -a gene.bed -b sorted.bam -sorted -g transcriptome_chromSizes.txt > cov.bed
```

</details>
