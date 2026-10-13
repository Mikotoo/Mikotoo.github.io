---
title: 'Haphic'
sidebar:
  label: 'Haphic'
  order: 999
---


## 中文
**作用**：用 HapHiC 依据 Hi-C 邻近信息把 contig 级组装挂载到染色体级
**入口**：`05.1.alignGetBam.sh`（比对、去 PCR 重复、过滤）→ `05.2.haphic.sh`（挂载）→ `05.build/juicebox.sh`（构建 Juicebox 输入供人工检查与校正）
**输入**：`asm.fa`（contig 级组装）；Hi-C 双端 FASTQ；`05.1` 产生的 `HiC.filtered.bam`。
**输出**：`HiC.bam`、`HiC.filtered.bam`；HapHiC 挂载输出（染色体级组装及附属文件）；Juicebox 构建目录。
**运行**：`bash 05.1.alignGetBam.sh`；`bash 05.2.haphic.sh`；`bash 05.build/juicebox.sh`
**工具**：bwa `index`、`mem -5SP -t 48`；samblaster；samtools `view -@ 48 -S -h -b -F 3340`；`filter_bam HiC.bam 1 --nm 3 --threads 48`；`haphic pipeline asm.fa HiC.filtered.bam 6 --threads 48 --processes 48 --quick_view`（6 为染色体数）；Juicebox `juicer pre` 加 `juicer_tools.1.9.9_jcuda.0.8.jar`（`-Xmx32G`，需 JVM）。

<details>
<summary>English</summary>

**Purpose**: Scaffold the contig-level assembly into chromosomes with HapHiC from Hi-C proximity data.
**Entry point**: `05.1.alignGetBam.sh` (align, remove PCR duplicates, filter) → `05.2.haphic.sh` (scaffolding) → `05.build/juicebox.sh` (build the Juicebox input for manual inspection and correction)
**Inputs**: `asm.fa` (contig-level assembly); Hi-C paired FASTQ files; `HiC.filtered.bam` produced by `05.1`.
**Outputs**: `HiC.bam`, `HiC.filtered.bam`; the HapHiC scaffolding output (chromosome-level assembly and supporting files); the Juicebox build directory.
**Run**: `bash 05.1.alignGetBam.sh`; `bash 05.2.haphic.sh`; `bash 05.build/juicebox.sh`
**Tools**: bwa `index`, `mem -5SP -t 48`; samblaster; samtools `view -@ 48 -S -h -b -F 3340`; `filter_bam HiC.bam 1 --nm 3 --threads 48`; `haphic pipeline asm.fa HiC.filtered.bam 6 --threads 48 --processes 48 --quick_view` (6 = number of chromosomes); Juicebox `juicer pre` plus `juicer_tools.1.9.9_jcuda.0.8.jar` (`-Xmx32G`, needs a JVM).

</details>

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>05.1.alignGetBam.sh</code> · Shell · 8 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/01_assembly/05_hic_scaffold/haphic/05.1.alignGetBam.sh)

```bash


 (1) Align Hi-C data to the assembly, remove PCR duplicates and filter out secondary and supplementary alignments
bwa index asm.fa
bwa mem -5SP -t 48 asm.fa {PROJ_LOTUS_ZC}/data/Gifu/hic/Gifu_hic_1.fastq {PROJ_LOTUS_ZC}/data/Gifu/hic/Gifu_hic_2.fastq | samblaster | samtools view - -@ 48 -S -h -b -F 3340 -o HiC.bam

 (2) Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
{SOFTWARE_ZC}/HapHiC-main/utils/filter_bam HiC.bam 1 --nm 3 --threads 48 | samtools view - -b -@ 48 -o HiC.filtered.bam
```

</details>

<details open>
<summary><code>05.2.haphic.sh</code> · Shell · 4 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/01_assembly/05_hic_scaffold/haphic/05.2.haphic.sh)

```bash


haphic pipeline asm.fa HiC.filtered.bam 6 --threads 48 --processes 48 --quick_view
haphic pipeline asm.fa HiC.filtered.bam 6 --threads 48 --processes 48 --inflation 3
```

</details>

<details open>
<summary><code>05.build/juicebox.sh</code> · Shell · 12 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/01_assembly/05_hic_scaffold/haphic/05.build/juicebox.sh)

```bash
#!/bin/bash
#CSUB -J haphic
#CSUB -q c01
#CSUB -o haphic.out
#CSUB -e haphic.error
#CSUB -n 88
#CSUB -R span[hosts=1]

ln -s {PROJ_LOTUS_ZC}/Gifu/05.haphic/asm.fa .
samtools faidx asm.fa
{SOFTWARE_ZC}/HapHiC-main/scripts/../utils/juicer pre -a -q 1 -o out_JBAT {PROJ_LOTUS_ZC}/Gifu/05.haphic/HiC.filtered.bam scaffolds.raw.agp asm.fa.fai >out_JBAT.log 2>&1
(java -Djava.awt.headless=true -jar -Xmx32G {SOFTWARE_ZC}/HapHiC-main/scripts/../utils/juicer_tools.1.9.9_jcuda.0.8.jar pre out_JBAT.txt out_JBAT.hic.part <(cat out_JBAT.log | grep PRE_C_SIZE | awk '{print $2" "$3}')) && (mv out_JBAT.hic.part out_JBAT.hic)
```

</details>
