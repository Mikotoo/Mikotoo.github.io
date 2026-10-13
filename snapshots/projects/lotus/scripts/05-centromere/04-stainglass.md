---
title: 'StainedGlass 分析'
sidebar:
  label: 'StainedGlass 分析'
  order: 4
---


## 中文

**作用**：用 StainedGlass 自比对展示 Gifu 着丝粒重复阵列的高阶结构，每条染色体一次：切出单条染色体、切成 2 kb 窗口、minimap2 全比全自比对，再用 StainedGlass 脚本转成 identity 表并绘图。
**入口**：`run_stainglass.sh`
**输入**：各 `Gifu_Chr<N>/` 内的 `Gifu_Chr<N>_2k.fasta` 与其 `.fai`；`get.sh` 另需 `mg20_cent.fa` 与 `.fai`。
**输出**：每个 `Gifu_Chr<N>/` 内的 `Gifu_Chr<N>_2k.mmi`、`Gifu_Chr<N>.bam`、`.tbl.gz`、`.full.tbl.gz`、`.bed.gz`，图件在 `results/Gifu_Chr<N>_figures/{pdfs,pngs}/`；`get.sh` 每条染色体写 `<id>.fa`、`<id>.fa.fai`、`<id>_2k.bed`、`<id>_2k.fasta`。
**运行**：`bash run_stainglass.sh`（在 `04_stainglass/` 下依次进入 `Gifu_Chr1`…`Gifu_Chr6`，每目录内建索引与自比对 → `samIdentity.py` → `bgzip` → `refmt.py` → `aln_plot.R`）。
**工具**：StainedGlass 0.6（`samIdentity.py --threads 88 --matches 400 --header`、`refmt.py --window 2000`、`aln_plot.R --threads 88`，工作流脚本来自 `SG=~/software/StainedGlass-0.6/workflow/scripts`）；minimap2 `-f 1000 -s 400 -ax ava-ont` 建索引、`-t 88 -f 10000 -s 400 -ax ava-ont --dual=yes --eqx` 比对；samtools `faidx`、`sort -m 4G`；bedtools `makewindows -w 2000`、`getfasta`；`bgzip`、`Rscript`。

<details>
<summary>English</summary>

**Purpose**: StainedGlass self-alignments showing the higher-order organisation of the Gifu centromeric repeat arrays, one chromosome at a time: each chromosome is cut out, windowed at 2 kb and aligned to itself all-versus-all with minimap2, then converted to identity tables with the StainedGlass scripts and plotted.
**Entry point**: `run_stainglass.sh`
**Inputs**: the 2-kb window FASTA and its `.fai` inside each `Gifu_Chr<N>/` directory (`Gifu_Chr<N>_2k.fasta`, `Gifu_Chr<N>.fa.fai`); `get.sh` additionally needs `mg20_cent.fa` and its `.fai`.
**Outputs**: inside each `Gifu_Chr<N>/` directory `Gifu_Chr<N>_2k.mmi`, `Gifu_Chr<N>.bam`, `.tbl.gz`, `.full.tbl.gz`, `.bed.gz`, with plots under `results/Gifu_Chr<N>_figures/{pdfs,pngs}/`; `get.sh` writes `<id>.fa`, `<id>.fa.fai`, `<id>_2k.bed` and `<id>_2k.fasta` per chromosome.
**Run**: `bash run_stainglass.sh` (from `04_stainglass/`, entering `Gifu_Chr1`…`Gifu_Chr6` in turn; each directory runs index and self alignment → `samIdentity.py` → `bgzip` → `refmt.py` → `aln_plot.R`).
**Tools**: StainedGlass 0.6 (`samIdentity.py --threads 88 --matches 400 --header`, `refmt.py --window 2000`, `aln_plot.R --threads 88`, with the workflow scripts taken from `SG=~/software/StainedGlass-0.6/workflow/scripts`); minimap2 index `-f 1000 -s 400 -ax ava-ont`, alignment `-t 88 -f 10000 -s 400 -ax ava-ont --dual=yes --eqx`; samtools `faidx`, `sort -m 4G`; bedtools `makewindows -w 2000`, `getfasta`; `bgzip`, `Rscript`.

</details>

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>get.sh</code> · Shell · 10 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/05_centromere/04_stainglass/get.sh)

```bash
cat mg20_cent.fa.fai|awk '{print $1}'|while read id;
do 
cd {PROJ_04_LOTUS_GENOME}/08.cent/stainGlass
mkdir ${id} && cd ${id}
echo ${id} > ${id}
perl ~/script/getseq.pl ${id} ../mg20_cent.fa > ${id}.fa
samtools faidx ${id}.fa
bedtools makewindows -g ${id}.fa.fai -w 2000 > ${id}_2k.bed
bedtools getfasta -fi ${id}.fa -bed ${id}_2k.bed > ${id}_2k.fasta
done
```

</details>

<details open>
<summary><code>process.sh</code> · Shell · 37 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/05_centromere/04_stainglass/process.sh)

```bash
#!/bin/bash
#CSUB -J process
#CSUB -q c01
#CSUB -o process.out
#CSUB -e process.error
#CSUB -n 88
#CSUB -R span[hosts=1]

bedtools getfasta -fi {HOME_ZC}/Project/Ref/ChiHei/ChiHei_v1.0.fasta -bed Chr13.bed > Chr13.fa
samtools faidx Chr13.fa
bedtools makewindows -g Chr13.fa.fai -w 2000 > Chr13_output.bed
bedtools getfasta -fi Chr13.fa -bed Chr13_output.bed > Chr13_output.2000.fasta

minimap2 -f 1000 -s 400 -ax ava-ont -d Chr13_output.fasta.mmi Chr13_output.2000.fasta
minimap2 -t 88 -f 10000 -s 400 -ax ava-ont --dual=yes --eqx Chr13_output.fasta.mmi Chr13_output.2000.fasta | samtools sort -m 4G -o Chr13_output.bam
python ~/software/StainedGlass-0.6/workflow/scripts/samIdentity.py --threads 88 --matches 400 --header Chr13_output.bam > Chr13_output.tbl
bgzip -c Chr13_output.tbl > Chr13_output.tbl.gz
python ~/software/StainedGlass-0.6/workflow/scripts/refmt.py --window 2000 --fai Chr13.fa.fai --full Chr13_output.full.tbl.gz Chr13_output.tbl.gz Chr13_output.bed.gz
mkdir -p results/Chr13_figures/pdfs
mkdir -p results/Chr13_figures/pngs
Rscript ~/software/StainedGlass-0.6/workflow/scripts/aln_plot.R -b Chr13_output.bed.gz --threads 88 --prefix Chr13



bedtools getfasta -fi {HOME_ZC}/Project/Ref/ChiHei/ChiHei_v1.0.fasta -bed Chr19.bed > Chr19.fa
samtools faidx Chr19.fa
bedtools makewindows -g Chr19.fa.fai -w 2000 > Chr19_output.bed
bedtools getfasta -fi Chr19.fa -bed Chr19_output.bed > Chr19_output.2000.fasta

minimap2 -f 1000 -s 400 -ax ava-ont -d Chr19_output.fasta.mmi Chr19_output.2000.fasta
minimap2 -t 88 -f 10000 -s 400 -ax ava-ont --dual=yes --eqx Chr19_output.fasta.mmi Chr19_output.2000.fasta | samtools sort -m 4G -o Chr19_output.bam
python ~/software/StainedGlass-0.6/workflow/scripts/samIdentity.py --threads 88 --matches 400 --header Chr19_output.bam > Chr19_output.tbl
bgzip -c Chr19_output.tbl > Chr19_output.tbl.gz
python ~/software/StainedGlass-0.6/workflow/scripts/refmt.py --window 2000 --fai Chr19.fa.fai --full Chr19_output.full.tbl.gz Chr19_output.tbl.gz Chr19_output.bed.gz
mkdir -p results/Chr19_figures/pdfs
mkdir -p results/Chr19_figures/pngs
Rscript ~/software/StainedGlass-0.6/workflow/scripts/aln_plot.R -b Chr19_output.bed.gz --threads 88 --prefix Chr19
```

</details>

<details open>
<summary><code>run_stainglass.sh</code> · Shell · 36 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/05_centromere/04_stainglass/run_stainglass.sh)

```bash
#!/bin/bash
#CSUB -J process
#CSUB -q c01
#CSUB -o process.out
#CSUB -e process.error
#CSUB -n 88
#CSUB -R span[hosts=1]

# StainedGlass 0.6 self-identity run for the six Gifu centromere regions.
#
# This is the loop form of the six per-region scripts that were run one
# directory at a time (Gifu_Chr1 … Gifu_Chr6). Each directory holds the region
# FASTA produced by get.sh: <region>.fa, <region>.fa.fai and <region>_2k.fasta.
#
# Run from 04_stainglass/, after get.sh:
#   bash run_stainglass.sh

SG=~/software/StainedGlass-0.6/workflow/scripts

for dir in Gifu_Chr*; do
    [ -d "$dir" ] || continue
    region="$dir"
    echo ">>> ${region}"
    cd "$dir"

    minimap2 -f 1000 -s 400 -ax ava-ont -d ${region}_2k.mmi ${region}_2k.fasta
    minimap2 -t 88 -f 10000 -s 400 -ax ava-ont --dual=yes --eqx ${region}_2k.mmi ${region}_2k.fasta | samtools sort -m 4G -o ${region}.bam
    python ${SG}/samIdentity.py --threads 88 --matches 400 --header ${region}.bam > ${region}.tbl
    bgzip -c ${region}.tbl > ${region}.tbl.gz
    python ${SG}/refmt.py --window 2000 --fai ${region}.fa.fai --full ${region}.full.tbl.gz ${region}.tbl.gz ${region}.bed.gz
    mkdir -p results/${region}_figures/pdfs
    mkdir -p results/${region}_figures/pngs
    Rscript ${SG}/aln_plot.R -b ${region}.bed.gz --threads 88 --prefix ${region}

    cd ..
done
```

</details>
