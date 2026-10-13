---
title: '品种间比对'
sidebar:
  label: '品种间比对'
  order: 2
---


## 中文

**作用**：MUMmer4 + SyRI 比较 MG20T2T（参考）与 GifuT2T（查询），输出各类变异 BED、事件与 bp 计数以及由比对区间补集得到的 PAV；`syri.sh` 是另一条更早的 minimap2 路线，覆盖三对基因组并合成一张 plotsr 图。
**入口**：`GifuT2T_vs_MG20T2T_variation_finder.sh`
**输入**：`Lotus_MG20T2T_v1.0.fasta`、`Lotus_GifuT2T_v1.0.fasta`（裸文件名）；PAV 另需上一步产出的 `GifuT2T_MG20T2T.filtered.coords`；`syri.sh` 用 ``{PROJ_LOTUS_ZC}`/comp/reference/` 下的 `Gifu_v1.fa`、`Gifu_T2T.fa`、`MG20_T2T.fa`、`MG20_v1.fa` 与 `genomes.txt`。
**输出**：`GifuT2T_MG20T2T.delta`/`.filtered.delta`、dnadiff 报告、`.filtered.coords`、`syri.out`、`MG20_<type>.bed` 与 `Gifu_<type>.bed`（SNP、INS、DEL、INV、TRANS、INVTR、DUP、INVDP、PAV）、计数表与 `variation.out`；PAV 出 `<asm>.aligned.bed`、`<asm>.PAV.bed`、`PAV.summary`；`syri.sh` 出三个 BAM 及 `.bai` 与 `syri_samechromosome_plot.pdf`。
**运行**：`bash GifuT2T_vs_MG20T2T_variation_finder.sh`，随后 `bash GifuT2T_vs_MG20T2T_get_PAV_from_coords.sh`；`bash syri.sh` 为独立路线（每个比较须在各自工作目录运行，否则 `syri.out`、BED 与临时文件互相覆盖）。
**工具**：MUMmer4 4.0.1（`nucmer -c 1000 --mum --maxgap=1000`、`delta-filter -i 90 -l 1000 -q`、`show-coords -THrd`、`dnadiff -d -p`）、SyRI 1.7.1、plotsr 1.1.1、samtools 1.21、bedtools 2.31.1（`merge -d 50`、`complement -g`）、minimap2（`-ax asm20`，仅 `syri.sh`）、GNU Awk 5.1.0、CSUB（`-q c01`、`-n 64`）。

<details>
<summary>English</summary>

**Purpose**: MUMmer4 + SyRI comparison of MG20T2T (reference) and GifuT2T (query) — per-type variant BEDs, event and base-pair counts, and PAV from the complement of the aligned intervals; `syri.sh` is a separate, earlier minimap2 route over three genome pairs ending in one combined plotsr plot.
**Entry point**: `GifuT2T_vs_MG20T2T_variation_finder.sh`
**Inputs**: `Lotus_MG20T2T_v1.0.fasta` and `Lotus_GifuT2T_v1.0.fasta` (bare filenames); the PAV step also needs `GifuT2T_MG20T2T.filtered.coords` from the step before; `syri.sh` uses `Gifu_v1.fa`, `Gifu_T2T.fa`, `MG20_T2T.fa` and `MG20_v1.fa` plus `genomes.txt` under ``{PROJ_LOTUS_ZC}`/comp/reference/`.
**Outputs**: `GifuT2T_MG20T2T.delta`/`.filtered.delta`, dnadiff reports, `.filtered.coords`, `syri.out`, `MG20_<type>.bed` and `Gifu_<type>.bed` (`type` = SNP, INS, DEL, INV, TRANS, INVTR, DUP, INVDP, PAV), the count table and `variation.out`; the PAV step writes `<asm>.aligned.bed`, `<asm>.PAV.bed` and `PAV.summary`; `syri.sh` writes three BAMs with `.bai` and `syri_samechromosome_plot.pdf`.
**Run**: `bash GifuT2T_vs_MG20T2T_variation_finder.sh`, then `bash GifuT2T_vs_MG20T2T_get_PAV_from_coords.sh`; `bash syri.sh` is a separate route (each comparison must run in its own working directory, otherwise `syri.out`, the BEDs and temporary files overwrite one another).
**Tools**: MUMmer4 4.0.1 (`nucmer -c 1000 --mum --maxgap=1000`, `delta-filter -i 90 -l 1000 -q`, `show-coords -THrd`, `dnadiff -d -p`), SyRI 1.7.1, plotsr 1.1.1, samtools 1.21, bedtools 2.31.1 (`merge -d 50`, `complement -g`), minimap2 (`-ax asm20`, `syri.sh` only), GNU Awk 5.1.0, CSUB (`-q c01`, `-n 64`).

</details>

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>GifuT2T_vs_MG20T2T_get_PAV_from_coords.sh</code> · Shell · 90 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/06_genome_comparison/02_accession_alignment/GifuT2T_vs_MG20T2T_get_PAV_from_coords.sh)

```bash
#!/usr/bin/env bash
set -euo pipefail

############################################
# 输入文件（请按需修改）
############################################
COORDS=GifuT2T_MG20T2T.filtered.coords

Gifu_T2T=Lotus_GifuT2T_v1.0.fasta
MG20_T2T=Lotus_MG20T2T_v1.0.fasta

PREFIX=GifuT2T_vs_MG20T2T

############################################
# 0. 准备 genome 文件
############################################
samtools faidx "$Gifu_T2T"
samtools faidx "$MG20_T2T"

cut -f1,2 ${Gifu_T2T}.fai > GifuT2T.genome
cut -f1,2 ${MG20_T2T}.fai > MG20T2T.genome

############################################
# 1. 提取 GifuT2T（query）比对区段
############################################
awk -F'\t' '
{
    chr = $11
    s = $3
    e = $4
    if (s > e) {tmp=s; s=e; e=tmp}
    if (e > s)
        print chr "\t" s-1 "\t" e
}
' "$COORDS" \
| sort -k1,1 -k2,2n \
| bedtools merge -d 50 \
> ${PREFIX}.GifuT2T.aligned.bed

############################################
# 2. 提取 MG20T2T（ref）比对区段
############################################
awk -F'\t' '
{
    chr = $10
    s = $1
    e = $2
    if (s > e) {tmp=s; s=e; e=tmp}
    if (e > s)
        print chr "\t" s-1 "\t" e
}
' "$COORDS" \
| sort -k1,1 -k2,2n \
| bedtools merge -d 50 \
> ${PREFIX}.MG20T2T.aligned.bed

############################################
# 3. 反选（complement）得到 PAV
############################################
bedtools complement \
  -i ${PREFIX}.GifuT2T.aligned.bed \
  -g GifuT2T.genome \
> ${PREFIX}.GifuT2T.PAV.bed

bedtools complement \
  -i ${PREFIX}.MG20T2T.aligned.bed \
  -g MG20T2T.genome \
> ${PREFIX}.MG20T2T.PAV.bed

############################################
# 4. 统计区间数量 & 总长度
############################################

for f in ${PREFIX}.GifuT2T.PAV.bed ${PREFIX}.MG20T2T.PAV.bed
do

    echo "File: $f" >> ${PREFIX}.PAV.summary

    awk '
    {
        len = $3 - $2 + 1
        if (len > 0) {
            n++
            sum += len
        }
    }
    END {
        printf("Intervals: %d\nTotal length (bp): %d\n", n, sum)
    }' "$f" >> ${PREFIX}.PAV.summary
done
```

</details>

<details open>
<summary><code>GifuT2T_vs_MG20T2T_variation_finder.sh</code> · Shell · 52 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/06_genome_comparison/02_accession_alignment/GifuT2T_vs_MG20T2T_variation_finder.sh)

```bash
#!/bin/bash
#CSUB -J variation
#CSUB -q c01
#CSUB -o variation.out
#CSUB -e variation.error
#CSUB -n 64
#CSUB -R span[hosts=1]

# 1. 全基因组比对（根据需要调整参数）
# 使用nucmer进行比对，参数可根据前述讨论选择 --mum 或 --maxmatch
# 这里示例采用更全面的--maxmatch并设置较高阈值避免噪音
nucmer -c 1000 --mum --maxgap=1000 --prefix=GifuT2T_MG20T2T Lotus_MG20T2T_v1.0.fasta Lotus_GifuT2T_v1.0.fasta

# 2. 过滤比对结果，提高比对质量
delta-filter -i 90 -l 1000 -q GifuT2T_MG20T2T.delta > GifuT2T_MG20T2T.filtered.delta

# 可选：输出比对统计和绘图（非必要步骤，仅检查比对情况）
dnadiff -d GifuT2T_MG20T2T.filtered.delta -p GifuT2T_MG20T2T.filtered
# show-coords 可用于生成比对坐标表格
show-coords -THrd GifuT2T_MG20T2T.filtered.delta > GifuT2T_MG20T2T.filtered.coords

# 3. 使用SyRI识别变异
syri -c GifuT2T_MG20T2T.filtered.coords -d GifuT2T_MG20T2T.filtered.delta \
    -r Lotus_MG20T2T_v1.0.fasta -q Lotus_GifuT2T_v1.0.fasta

# 4. 按类别提取变异坐标并保存为BED文件
# 小变异类别：SNP、插入、缺失
awk '$11=="SNP" {print $1"\t"$2"\t"$3 >> "MG20_SNP.bed"; print $6"\t"$7"\t"$8 >> "Gifu_SNP.bed"}' syri.out
awk '$11=="INS" {print $1"\t"$2"\t"$3 >> "MG20_INS.bed"; print $6"\t"$7"\t"$8 >> "Gifu_INS.bed"}' syri.out
awk '$11=="DEL" {print $1"\t"$2"\t"$3 >> "MG20_DEL.bed"; print $6"\t"$7"\t"$8 >> "Gifu_DEL.bed"}' syri.out

# 结构变异类别：倒位、易位、倒置易位、重复、倒置重复
awk '$11=="INV" {print $1"\t"$2"\t"$3 >> "MG20_INV.bed"; print $6"\t"$7"\t"$8 >> "Gifu_INV.bed"}' syri.out
awk '$11=="TRANS" {print $1"\t"$2"\t"$3 >> "MG20_TRANS.bed"; print $6"\t"$7"\t"$8 >> "Gifu_TRANS.bed"}' syri.out
awk '$11=="INVTR" {print $1"\t"$2"\t"$3 >> "MG20_INVTR.bed"; print $6"\t"$7"\t"$8 >> "Gifu_INVTR.bed"}' syri.out
awk '$11=="DUP" {print $1"\t"$2"\t"$3 >> "MG20_DUP.bed"; print $6"\t"$7"\t"$8 >> "Gifu_DUP.bed"}' syri.out
awk '$11=="INVDP" {print $1"\t"$2"\t"$3 >> "MG20_INVDP.bed"; print $6"\t"$7"\t"$8 >> "Gifu_INVDP.bed"}' syri.out

# 存在/缺失变异（PAV）：参考或查询中特有的大片段（INS/DEL/NOTAL >50bp）
awk '($11=="INS" || $11=="DEL" || $11=="NOTAL") && (($3-$2+1>50) || ($8-$7+1>50)) \
    {print $1"\t"$2"\t"$3 >> "MG20_PAV.bed"; print $6"\t"$7"\t"$8 >> "Gifu_PAV.bed"}' syri.out

# 5. 统计各类别变异事件数量和总长度
echo -e "Category\tCount\tTotal_length_in_MG20\tTotal_length_in_Gifu"
for type in SNP INS DEL INV TRANS INVTR DUP INVDP PAV; do
    # 统计事件数（行数）
    count=$(wc -l < "MG20_${type}.bed")
    # 计算参考和查询总长度（bp）
    total_ref=$(awk '{sum += $3-$2+1} END{print sum}' "MG20_${type}.bed")
    total_qry=$(awk '{sum += $3-$2+1} END{print sum}' "Gifu_${type}.bed")
    echo -e "${type}\t${count}\t${total_ref}\t${total_qry}"
done
```

</details>

<details open>
<summary><code>syri.sh</code> · Shell · 40 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/06_genome_comparison/02_accession_alignment/syri.sh)

```bash
#!/bin/bash
#CSUB -J process
#CSUB -q c01
#CSUB -o process.out
#CSUB -e process.error
#CSUB -n 64
#CSUB -R span[hosts=1]

minimap2 -ax asm20 -t 64 --eqx {PROJ_LOTUS_ZC}/comp/reference/Gifu_v1.fa \
{PROJ_LOTUS_ZC}/comp/reference/Gifu_T2T.fa \
| samtools sort -O BAM - > Gifu_v1_Gifu_T2T.bam
samtools index Gifu_v1_Gifu_T2T.bam

syri -c Gifu_v1_Gifu_T2T.bam -r {PROJ_LOTUS_ZC}/comp/reference/Gifu_v1.fa \
-q {PROJ_LOTUS_ZC}/comp/reference/Gifu_T2T.fa -F B --prefix Gifu_v1_Gifu_T2T


minimap2 -ax asm20 -t 64 --eqx {PROJ_LOTUS_ZC}/comp/reference/Gifu_T2T.fa \
{PROJ_LOTUS_ZC}/comp/reference/MG20_T2T.fa \
| samtools sort -O BAM - > Gifu_T2T_MG20_T2T.bam
samtools index Gifu_T2T_MG20_T2T.bam

syri -c Gifu_T2T_MG20_T2T.bam -r {PROJ_LOTUS_ZC}/comp/reference/Gifu_T2T.fa \
-q {PROJ_LOTUS_ZC}/comp/reference/MG20_T2T.fa -F B --prefix Gifu_T2T_MG20_T2T


minimap2 -ax asm20 -t 64 --eqx {PROJ_LOTUS_ZC}/comp/reference/MG20_T2T.fa \
{PROJ_LOTUS_ZC}/comp/reference/MG20_v1.fa \
| samtools sort -O BAM - > MG20_T2T_MG20_v1.bam
samtools index MG20_T2T_MG20_v1.bam

syri -c MG20_T2T_MG20_v1.bam -r {PROJ_LOTUS_ZC}/comp/reference/MG20_T2T.fa \
-q {PROJ_LOTUS_ZC}/comp/reference/MG20_v1.fa -F B --prefix MG20_T2T_MG20_v1


plotsr --sr Gifu_v1_Gifu_T2Tsyri.out \
--sr Gifu_T2T_MG20_T2Tsyri.out \
--sr MG20_T2T_MG20_v1syri.out \
--genomes genomes.txt \
-o syri_samechromosome_plot.pdf
```

</details>
