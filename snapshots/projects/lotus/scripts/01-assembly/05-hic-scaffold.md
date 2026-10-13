---
title: 'Hic scaffold'
sidebar:
  label: 'Hic scaffold'
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
| [Haphic](/_projects/lotus/scripts/01-assembly/05-hic-scaffold/haphic/) | — |

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>cooler/matrix2cool.sh</code> · Shell · 48 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/01_assembly/05_hic_scaffold/cooler/matrix2cool.sh)

```bash
#!/bin/bash
#CSUB -J cooler
#CSUB -q c01
#CSUB -o cooler.out
#CSUB -e cooler.error
#CSUB -n 64
#CSUB -R span[hosts=1]

cat {PROJ_04_LOTUS_GENOME}/02_assembly/Gifu/09.hic/process/hic_results/matrix/Gifu/raw/200000/Gifu_200000.matrix \
| awk '$1<=$2{print $0}' > {PROJ_04_LOTUS_GENOME}/02_assembly/Gifu/09.hic/process/hic_results/matrix/Gifu/raw/200000/Gifu_200000_up.matrix

cooler load -f coo --one-based --count-as-float \
{PROJ_04_LOTUS_GENOME}/02_assembly/Gifu/09.hic/process/hic_results/matrix/Gifu/raw/200000/Gifu_200000_abs.bed \
{PROJ_04_LOTUS_GENOME}/02_assembly/Gifu/09.hic/process/hic_results/matrix/Gifu/raw/200000/Gifu_200000_up.matrix \
Gifu_200kb.cool


cat {PROJ_04_LOTUS_GENOME}/02_assembly/Gifu/09.hic/process/hic_results/matrix/Gifu/raw/100000/Gifu_100000.matrix \
| awk '$1<=$2{print $0}' > {PROJ_04_LOTUS_GENOME}/02_assembly/Gifu/09.hic/process/hic_results/matrix/Gifu/raw/100000/Gifu_100000_up.matrix

cooler load -f coo --one-based --count-as-float \
{PROJ_04_LOTUS_GENOME}/02_assembly/Gifu/09.hic/process/hic_results/matrix/Gifu/raw/100000/Gifu_100000_abs.bed \
{PROJ_04_LOTUS_GENOME}/02_assembly/Gifu/09.hic/process/hic_results/matrix/Gifu/raw/100000/Gifu_100000_up.matrix \
Gifu_100kb.cool


cat {PROJ_04_LOTUS_GENOME}/02_assembly/Gifu/09.hic/process/hic_results/matrix/Gifu/raw/40000/Gifu_40000.matrix \
| awk '$1<=$2{print $0}' > {PROJ_04_LOTUS_GENOME}/02_assembly/Gifu/09.hic/process/hic_results/matrix/Gifu/raw/40000/Gifu_40000_up.matrix

cooler load -f coo --one-based --count-as-float \
{PROJ_04_LOTUS_GENOME}/02_assembly/Gifu/09.hic/process/hic_results/matrix/Gifu/raw/40000/Gifu_40000_abs.bed \
{PROJ_04_LOTUS_GENOME}/02_assembly/Gifu/09.hic/process/hic_results/matrix/Gifu/raw/40000/Gifu_40000_up.matrix \
Gifu_40kb.cool



cooler create --mode mcool -o Gifu.mcool
cooler cp Gifu_200kb.cool Gifu.mcool::resolutions/200000
cooler cp Gifu_100kb.cool Gifu.mcool::resolutions/100000
cooler cp Gifu_40kb.cool Gifu.mcool::resolutions/40000


for resolution_uri in $(cooler ls Gifu.mcool); do
    cooler balance \
        --ignore-diags 2 \
        --force \
        "$resolution_uri"
done
```

</details>

<details open>
<summary><code>hicpro.sh</code> · Shell · 3 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/01_assembly/05_hic_scaffold/hicpro.sh)

```bash


HiC-Pro --input rawreads -c config-hicpro.txt --output process
```

</details>

<details open>
<summary><code>reference/index.sh</code> · Shell · 9 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/01_assembly/05_hic_scaffold/reference/index.sh)

```bash
#!/bin/bash
#CSUB -J bwa
#CSUB -q c01
#CSUB -o index.out
#CSUB -e index.error
#CSUB -n 88
#CSUB -R span[hosts=1]

bowtie2-build -f Lotus_GifuT2T_v1.0.fasta --thread 88 Gifu_T2T
```

</details>
