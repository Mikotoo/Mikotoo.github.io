---
title: '去冗余单倍型'
sidebar:
  label: '去冗余单倍型'
  order: 2
---


## 中文
**作用**：把 HiFi reads 回贴到 hifiasm primary contig 统计读深与覆盖度，并去除单倍型冗余（haplotig）contig。
**入口**：`02.1.minimap.sh`（比对）→ `02.2.purge.sh`（去冗余）
**输入**：`Gifu.hifiasm.fasta`（步骤 01 输出）；HiFi reads。
**输出**：`aligned.bam`、`aligned.bam.200.gencov`、`coverage_stats.csv`；去冗余后的 contig FASTA 与 haplotig FASTA（`purge_haplotigs purge` 的默认名）。
**运行**：`bash 02.1.minimap.sh`，然后 `bash 02.2.purge.sh`
**工具**：minimap2 `-ax map-pb`；samtools `view -hF 256`、`sort -@ 20 -m 1G -T tmp.ali`；purge_haplotigs `readhist -t 88`、`contigcov -l 15 -m 100 -h 185`、`purge -t 88 -a 50`。

<details>
<summary>English</summary>

**Purpose**: Align HiFi reads back to the hifiasm primary contigs to measure read depth and coverage, then purge haplotig (duplicated) contigs.
**Entry point**: `02.1.minimap.sh` (alignment) → `02.2.purge.sh` (purging)
**Inputs**: `Gifu.hifiasm.fasta` (step 01 output); HiFi reads.
**Outputs**: `aligned.bam`, `aligned.bam.200.gencov`, `coverage_stats.csv`; the purged contig FASTA and the haplotig FASTA (the `purge_haplotigs purge` default names).
**Run**: `bash 02.1.minimap.sh`, then `bash 02.2.purge.sh`
**Tools**: minimap2 `-ax map-pb`; samtools `view -hF 256`, `sort -@ 20 -m 1G -T tmp.ali`; purge_haplotigs `readhist -t 88`, `contigcov -l 15 -m 100 -h 185`, `purge -t 88 -a 50`.

</details>

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>02.1.minimap.sh</code> · Shell · 5 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/01_assembly/02_purge_haplotigs/02.1.minimap.sh)

```bash
#2.1 把hifi reads回帖到1.1的初步组装结果上

minimap2 -ax map-pb Gifu.hifiasm.fasta {hifiReads}.fastq.gz \
|samtools view -hF 256 - \
|samtools sort -@ 20 -m 1G -o aligned.bam -T tmp.ali
```

</details>

<details open>
<summary><code>02.2.purge.sh</code> · Shell · 5 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/01_assembly/02_purge_haplotigs/02.2.purge.sh)

```bash
#2.2 根据比对结果，去除单倍型冗余contig

purge_haplotigs readhist -b aligned.bam -g Gifu.hifiasm.fasta -t 88
purge_haplotigs contigcov -i aligned.bam.200.gencov -o coverage_stats.csv -l 15 -m 100 -h 185
purge_haplotigs purge -g Gifu.hifiasm.fasta -c coverage_stats.csv -b aligned.bam -t 88 -a 50
```

</details>
