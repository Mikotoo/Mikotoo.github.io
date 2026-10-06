---
title: '数据准备'
sidebar:
  label: '数据准备'
  order: 0
---


## 中文
**作用**：bulk RNA-seq 上游处理——fastp 质控、HISAT2 比对到 Gifu T2T 参考、StringTie 逐样本定量，并把各样本 TPM 列合并为一张矩阵。
**入口**：`fastp.sh`、`align.sh`
**输入**：`align.conf`（样本名、R1 FASTQ、R2 FASTQ 三列）及其指向的 `raw/<sample>-R1.fq.gz`、`-R2.fq.gz`；建索引用 `Gifu_v1.0.fasta`，StringTie 用 `Lotus_GifuT2T.gtf`。
**输出**：`clean/<sample>-R{1,2}.fastq.gz`、`report/<sample>.json|.html`、排序并建索引的 `<sample>.bam(.bai)`、`02.stringtie/<sample>/<sample>.gtf|.tsv`、`02.tpm/<sample>.tpm` 与合并矩阵 `02.tpm/all_samples.tpm.tsv`。
**运行**：`bash fastp.sh`、`bash align.sh`（`fastp.sh` 须在本目录内运行；`align.sh` 的 `$index` 从未赋值、`-x` 拿不到索引路径，且须以 `01.align` 为当前目录，运行前需修正脚本）。
**工具**：fastp（`-w 88`，JSON/HTML 报告）、HISAT2（`hisat2-build -p 64`；`-p 64 --new-summary`）、samtools（`view -H`、`sort -@ 64`、`index`）、StringTie（`-p 64 -G -e -B -o -A`）、Awk 与 `sort`/`join` 合并 TPM。

<details>
<summary>English</summary>
**Purpose**: Upstream bulk RNA-seq processing — fastp read QC, HISAT2 alignment to the Gifu T2T reference and StringTie per-sample quantification, joined into one TPM matrix.
**Entry point**: `fastp.sh`, `align.sh`
**Inputs**: `align.conf` (sample name, R1 FASTQ, R2 FASTQ) and the `raw/<sample>-R1.fq.gz` / `-R2.fq.gz` files it points to; `Gifu_v1.0.fasta` for the index and `Lotus_GifuT2T.gtf` for StringTie.
**Outputs**: `clean/<sample>-R{1,2}.fastq.gz`, `report/<sample>.json|.html`, the sorted and indexed `<sample>.bam(.bai)`, `02.stringtie/<sample>/<sample>.gtf|.tsv`, `02.tpm/<sample>.tpm` and the joined matrix `02.tpm/all_samples.tpm.tsv`.
**Run**: `bash fastp.sh`, `bash align.sh` (`fastp.sh` must run in this directory; in `align.sh` `$index` is never set, so `-x` receives no index path, and `01.align` must be the current directory — fix the script before running).
**Tools**: fastp (`-w 88`, JSON/HTML reports), HISAT2 (`hisat2-build -p 64`; `-p 64 --new-summary`), samtools (`view -H`, `sort -@ 64`, `index`), StringTie (`-p 64 -G -e -B -o -A`), Awk with `sort`/`join` for the TPM merge.

</details>


