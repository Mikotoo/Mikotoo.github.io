---
title: 'Hi-C 挂载'
sidebar:
  label: 'Hi-C 挂载'
  order: -1
---


## 中文
**作用**：用 HiC-Pro 3.1.0 处理 Hi-C reads 生成接触矩阵，再转成 cooler `.cool`/`.mcool` 并做平衡化。
**入口**：`reference/index.sh` → `hicpro.sh` → `cooler/matrix2cool.sh`（`haphic/` 为另一条挂载路线，见其 README）
**输入**：`rawreads/` 与 `config-hicpro.txt`（`hicpro.sh` 使用）；`Lotus_GifuT2T_v1.0.fasta`（`reference/index.sh` 建索引）；HiC-Pro 矩阵树 `process/hic_results/matrix/Gifu/raw/{200000,100000,40000}/Gifu_<res>.matrix` 与 `..._abs.bed`（在 `cooler/matrix2cool.sh` 中硬编码）。
**输出**：`process/` 下的 HiC-Pro 运行树（矩阵、`*_abs.bed`、日志）与 `Gifu_T2T.*` bowtie2 索引；`Gifu_200kb.cool`、`Gifu_100kb.cool`、`Gifu_40kb.cool` 及平衡化后的多分辨率 `Gifu.mcool`。
**运行**：`bash reference/index.sh`；`bash hicpro.sh`；`bash cooler/matrix2cool.sh`
**工具**：`bowtie2-build -f --thread 88`（索引前缀 `Gifu_T2T`）；HiC-Pro 3.1.0 `--input rawreads -c config-hicpro.txt --output process`；cooler `load -f coo --one-based --count-as-float`、`create --mode mcool`、`balance --ignore-diags 2 --force`。

<details>
<summary>English</summary>
**Purpose**: Process Hi-C reads with HiC-Pro 3.1.0 into contact matrices, then convert them to cooler `.cool`/`.mcool` files and balance them.
**Entry point**: `reference/index.sh` → `hicpro.sh` → `cooler/matrix2cool.sh` (`haphic/` is a separate scaffolding route, see its README)
**Inputs**: `rawreads/` and `config-hicpro.txt` (used by `hicpro.sh`); `Lotus_GifuT2T_v1.0.fasta` (indexed by `reference/index.sh`); the HiC-Pro matrix tree `process/hic_results/matrix/Gifu/raw/{200000,100000,40000}/Gifu_<res>.matrix` with `..._abs.bed`, hardcoded in `cooler/matrix2cool.sh`.
**Outputs**: the HiC-Pro run tree under `process/` (matrices, `*_abs.bed`, logs) and the `Gifu_T2T.*` bowtie2 index; `Gifu_200kb.cool`, `Gifu_100kb.cool`, `Gifu_40kb.cool` and the balanced multi-resolution `Gifu.mcool`.
**Run**: `bash reference/index.sh`; `bash hicpro.sh`; `bash cooler/matrix2cool.sh`
**Tools**: `bowtie2-build -f --thread 88` (index basename `Gifu_T2T`); HiC-Pro 3.1.0 `--input rawreads -c config-hicpro.txt --output process`; cooler `load -f coo --one-based --count-as-float`, `create --mode mcool`, `balance --ignore-diags 2 --force`.

</details>

### 本阶段步骤

| 步骤 | 内容 |
| --- | --- |
| [HapHiC 分型挂载](/_projects/lotus/scripts/01-assembly/09-hic-scaffold/haphic/) | — |


