---
title: 'IsoSeq 全长转录本'
sidebar:
  label: 'IsoSeq 全长转录本'
  order: 3
---


## 中文

**作用**：处理 Gifu PacBio Iso-Seq 文库：`isoseq3 refine` 去引物、`isoseq3 cluster` 得到高质量全长一致序列，再比对到 T2T 组装得到长读转录证据 BAM。
**入口**：`isoseq.sh`。
**输入**：脚本内硬编码的 `sample1` = ``{PROJ_LOTUS_ZC}`/data/Gifu/isoseq/hifi_reads/P22TR251406648-1-bc12_r84069_20250617_081306_1_A01.hifi_reads.bam` 与同前缀的 `primer1`（`.primer.fasta`）；参考基因组 ``{PROJ_LOTUS_ZC}`/Gifu/00.assemble_result/Gifu_v1.0.fasta`。
**输出**：`Gifu.refined.bam`、`Gifu.clustered.bam`、`Gifu.clustered.fastq`、`Gifu_isoforms.sam`、`Gifu_isoforms.sorted.bam`（含索引），后者被 `04_braker/Gifu_braker.sh` 使用。
**运行**：`bash isoseq.sh`；仅 bc12 文库生效，覆盖样本 1–4 的 `bamtools merge` 循环与 TAMA/`gffread`/TransDecoder 段均被注释。
**工具**：isoseq3 `refine "$sample1" "$primer1" Gifu.refined.bam`、`cluster Gifu.refined.bam Gifu.clustered.bam --verbose --use-qvs`；minimap2 `-ax splice -uf -k14 <genome> Gifu.clustered.fastq`；samtools `fastq`、`sort -o Gifu_isoforms.sorted.bam`、`index`。

<details>
<summary>English</summary>
**Purpose**: Process the Gifu PacBio Iso-Seq library — `isoseq3 refine` removes primers, `isoseq3 cluster` builds high-quality full-length consensus isoforms, and the isoforms are aligned to the T2T assembly to give a long-read transcript-evidence BAM.
**Entry point**: `isoseq.sh`.
**Inputs**: the hardcoded `sample1` = ``{PROJ_LOTUS_ZC}`/data/Gifu/isoseq/hifi_reads/P22TR251406648-1-bc12_r84069_20250617_081306_1_A01.hifi_reads.bam` and `primer1` = the `.primer.fasta` with the same prefix; reference genome ``{PROJ_LOTUS_ZC}`/Gifu/00.assemble_result/Gifu_v1.0.fasta`.
**Outputs**: `Gifu.refined.bam`, `Gifu.clustered.bam`, `Gifu.clustered.fastq`, `Gifu_isoforms.sam`, `Gifu_isoforms.sorted.bam` (plus index), the last consumed by `04_braker/Gifu_braker.sh`.
**Run**: `bash isoseq.sh`; only the bc12 library is used, while the loop over samples 1–4 with `bamtools merge` and the TAMA/`gffread`/TransDecoder blocks are commented out.
**Tools**: isoseq3 `refine "$sample1" "$primer1" Gifu.refined.bam`, `cluster Gifu.refined.bam Gifu.clustered.bam --verbose --use-qvs`; minimap2 `-ax splice -uf -k14 <genome> Gifu.clustered.fastq`; samtools `fastq`, `sort -o Gifu_isoforms.sorted.bam`, `index`.

</details>


