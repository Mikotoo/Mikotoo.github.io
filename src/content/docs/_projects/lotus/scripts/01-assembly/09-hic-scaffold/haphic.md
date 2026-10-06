---
title: 'HapHiC 分型挂载'
sidebar:
  label: 'HapHiC 分型挂载'
  order: 999
---


## 中文
**作用**：用 HapHiC 依据 Hi-C 邻近信息把 contig 级组装挂载到染色体级；父目录的 HiC-Pro 3.1.0 流程只产出接触矩阵。
**入口**：`05.1.alignGetBam.sh`（比对、去 PCR 重复、过滤）→ `05.2.haphic.sh`（挂载）→ `04.build/juicebox.sh`（构建 Juicebox 输入供人工检查与校正）
**输入**：`asm.fa`（contig 级组装）；Hi-C 双端 FASTQ；`05.1` 产生的 `HiC.filtered.bam`。
**输出**：`HiC.bam`、`HiC.filtered.bam`；HapHiC 挂载输出（染色体级组装及附属文件）；Juicebox 构建目录。
**运行**：`bash 05.1.alignGetBam.sh`；`bash 05.2.haphic.sh`；`bash 04.build/juicebox.sh`
**工具**：bwa `index`、`mem -5SP -t 48`；samblaster；samtools `view -@ 48 -S -h -b -F 3340`；`filter_bam HiC.bam 1 --nm 3 --threads 48`；`haphic pipeline asm.fa HiC.filtered.bam 6 --threads 48 --processes 48 --quick_view`（6 为染色体数）；Juicebox `juicer pre` 加 `juicer_tools.1.9.9_jcuda.0.8.jar`（`-Xmx32G`，需 JVM）。
**注意**：`05.1.alignGetBam.sh` 中 `bwa mem` 被注释、生效的是 `filter_bam` 管道；`04.build/juicebox.sh` 硬编码 `Project/Lotus/Gifu/05.haphic` 且分配 `-Xmx32G`，不能在笔记本上运行。

<details>
<summary>English</summary>
**Purpose**: Scaffold the contig-level assembly into chromosomes with HapHiC from Hi-C proximity data; the HiC-Pro 3.1.0 workflow in the parent directory only produces contact matrices.
**Entry point**: `05.1.alignGetBam.sh` (align, remove PCR duplicates, filter) → `05.2.haphic.sh` (scaffolding) → `04.build/juicebox.sh` (build the Juicebox input for manual inspection and correction)
**Inputs**: `asm.fa` (contig-level assembly); Hi-C paired FASTQ files; `HiC.filtered.bam` produced by `05.1`.
**Outputs**: `HiC.bam`, `HiC.filtered.bam`; the HapHiC scaffolding output (chromosome-level assembly and supporting files); the Juicebox build directory.
**Run**: `bash 05.1.alignGetBam.sh`; `bash 05.2.haphic.sh`; `bash 04.build/juicebox.sh`
**Tools**: bwa `index`, `mem -5SP -t 48`; samblaster; samtools `view -@ 48 -S -h -b -F 3340`; `filter_bam HiC.bam 1 --nm 3 --threads 48`; `haphic pipeline asm.fa HiC.filtered.bam 6 --threads 48 --processes 48 --quick_view` (6 = number of chromosomes); Juicebox `juicer pre` plus `juicer_tools.1.9.9_jcuda.0.8.jar` (`-Xmx32G`, needs a JVM).
**Notes**: In `05.1.alignGetBam.sh` the `bwa mem` command is commented out and the active line is the `filter_bam` pipe; `04.build/juicebox.sh` hardcodes `Project/Lotus/Gifu/05.haphic` and allocates `-Xmx32G`, so it cannot run on a laptop.

</details>


