---
title: '打磨校正'
sidebar:
  label: '打磨校正'
  order: 7
---


## 中文
**作用**：对补洞后的组装做长读长与 k-mer 迭代抛光：Winnowmap2 比对、falconc 过滤、Racon 共识、Merfin/bcftools 校正。
**入口**：`bash.sh`（调用第三方 `polish.sh`；`build_merylDB.sh` 构建 k-mer 库）
**输入**：`ontPolish.iter_2.consensus.fasta`（生效调用所用的草图）；HiFi reads ``{PROJ_LOTUS_ZC}`/data/Gifu/hifi/Gifu_hifi.fastq.gz`；读 k-mer 库 `ccs.merylDB`。
**输出**：`ccs.merylDB`；每轮 `<prefix>.iter_N.*`（`winnowmap.bam`、`falconc.sam`、`racon.fasta`、merfin 校正 VCF、`consensus.fasta`）与 `.memtime` 日志。
**运行**：`bash bash.sh`
**工具**：`polish.sh 10 2 ontPolish.iter_2.consensus.fasta <hifi.fastq.gz> ccs.merylDB hifiPolish`；meryl `count k=21`；winnowmap `-k 21 -ax map-pb`；falconc `bam-filter-clipped -F 0x104`；racon `-L -S`；merfin `-polish -peak 106.7`；bcftools `consensus -H 1`。
**注意**：`bash.sh` 中前缀 `ontPolish` 的 ONT 运行被注释，生效的 HiFi 运行以它的输出为草图，须先备好该文件；`polish.sh` 为 `third_party/` 下未修改的第三方代码，以绝对路径调用。

<details>
<summary>English</summary>
**Purpose**: Iteratively polish the gap-closed assembly with long reads and k-mers: Winnowmap2 mapping, falconc filtering, Racon consensus and Merfin/bcftools correction.
**Entry point**: `bash.sh` (calls the third-party `polish.sh`; `build_merylDB.sh` builds the k-mer database)
**Inputs**: `ontPolish.iter_2.consensus.fasta` (the draft used by the active call); HiFi reads ``{PROJ_LOTUS_ZC}`/data/Gifu/hifi/Gifu_hifi.fastq.gz`; read k-mer database `ccs.merylDB`.
**Outputs**: `ccs.merylDB`; per iteration `<prefix>.iter_N.*` (`winnowmap.bam`, `falconc.sam`, `racon.fasta`, the Merfin correction VCF, `consensus.fasta`) and `.memtime` logs.
**Run**: `bash bash.sh`
**Tools**: `polish.sh 10 2 ontPolish.iter_2.consensus.fasta <hifi.fastq.gz> ccs.merylDB hifiPolish`; meryl `count k=21`; winnowmap `-k 21 -ax map-pb`; falconc `bam-filter-clipped -F 0x104`; racon `-L -S`; merfin `-polish -peak 106.7`; bcftools `consensus -H 1`.
**Notes**: The ONT run with prefix `ontPolish` in `bash.sh` is commented out while the active HiFi run takes its output as the draft, so that file must exist first; `polish.sh` is unmodified third-party code under `third_party/`, called by an absolute path.

</details>


