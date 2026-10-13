---
title: 'T2T 与前版本比较'
sidebar:
  label: 'T2T 与前版本比较'
  order: 5
---


## 中文

**作用**：把每个新 T2T 组装与同一 accession 的已发表旧组装比较（Gifu T2T vs Gifu old；MG20 T2T vs MG20 old）：MUMmer4 + SyRI 变异、由比对区间补集得到的 PAV，以及 T2T 特有新序列按注释区域的分区。
**入口**：`GifuT2T_GifuOld/GifuT2T_vs_GifuOld_variation_finder.sh`
**输入**：Gifu 用 `Lotus_Gifu_Old.fa`（seqkit 清理为 `Lotus_Gifu_Old1.fa`）与 `Lotus_GifuT2T_v1.0.fasta`，MG20 用 `Lotus_MG20_Old.fa` 与 `Lotus_MG20T2T_v1.0.fasta`，以及 `genomes.txt`；PAV 需对应的 `<prefix>.filtered.coords`；`calc_bin_coverage.sh` 用 `gifuT2T_newseq.bed` 与 `{PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Gifu_100k.bed`；`run.sh`/`run.py` 用 `gifuT2T_newseq.bed`（或 `mg20T2T_newseq.bed`）与 `cent.bed`、`rDNA.bed`、`tel.bed`、`repeat2.bed`。
**输出**：`<p>.delta`/`.filtered.delta`、dnadiff 报告、`syri.out`、`Old_<type>.bed` 与 `T2T_<type>.bed`（SNP、INS、DEL、INV、TRANS、INVTR、DUP、INVDP、PAV）、stdout 计数表与 `*_findSV.pdf`；PAV 出 `.genome`、`<asm>.aligned.bed`、`<asm>.PAV.bed`、`PAV.summary`；`calc_bin_coverage.sh` 出 `gifu_100kb_coverage.tsv`；`run.sh`/`run.py` 出各区域 bp 与占比及 `new_sequences_region_proportions.pdf`。
**运行**：每个比较各自一个工作目录，先跑 variation finder 再跑 PAV：`cd GifuT2T_GifuOld && bash GifuT2T_vs_GifuOld_variation_finder.sh && bash GifuT2T_vs_GifuOld_get_PAV_from_coords.sh`（MG20 同理）；随后 `bash calc_bin_coverage.sh`，并在 `Gifu_new/`、`MG20_new/` 中各跑 `bash run.sh` 或 `python3 run.py`。
**工具**：MUMmer4 4.0.1（`nucmer -c 1000 --mum --maxgap=1000`、`delta-filter -i 90 -l 1000 -q`、`show-coords -THrd`、`dnadiff -d -p`）、SyRI 1.7.1、plotsr 1.1.1、samtools 1.21、bedtools 2.31.1（`merge -d 50`、`complement -g`、`intersect`、`subtract`）、seqkit 2.10.1（`grep -v -r -p`）、GNU Awk 5.1.0、Python 3 + pandas/matplotlib、CSUB（`-q c01`、`-n 64`）。

<details>
<summary>English</summary>

**Purpose**: Compare each new T2T assembly with the published earlier assembly of the same accession (Gifu T2T vs Gifu old; MG20 T2T vs MG20 old): MUMmer4 + SyRI variants, PAV from the complement of the aligned intervals, and a partition of the T2T-specific new sequences by the annotated regions they fall in.
**Entry point**: `GifuT2T_GifuOld/GifuT2T_vs_GifuOld_variation_finder.sh`
**Inputs**: for Gifu, `Lotus_Gifu_Old.fa` (cleaned by seqkit to `Lotus_Gifu_Old1.fa`) and `Lotus_GifuT2T_v1.0.fasta`; for MG20, `Lotus_MG20_Old.fa` and `Lotus_MG20T2T_v1.0.fasta`; plus `genomes.txt`; the PAV step needs the matching `<prefix>.filtered.coords`; `calc_bin_coverage.sh` uses `gifuT2T_newseq.bed` and `{PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Gifu_100k.bed`; `run.sh`/`run.py` use `gifuT2T_newseq.bed` (or `mg20T2T_newseq.bed`) with `cent.bed`, `rDNA.bed`, `tel.bed` and `repeat2.bed`.
**Outputs**: `<p>.delta`/`.filtered.delta`, dnadiff reports, `syri.out`, `Old_<type>.bed` and `T2T_<type>.bed` (`type` = SNP, INS, DEL, INV, TRANS, INVTR, DUP, INVDP, PAV), a count table on stdout and `*_findSV.pdf`; the PAV step writes `.genome`, `<asm>.aligned.bed`, `<asm>.PAV.bed` and `PAV.summary`; `calc_bin_coverage.sh` writes `gifu_100kb_coverage.tsv`; `run.sh`/`run.py` print per-region bp and ratios and write `new_sequences_region_proportions.pdf`.
**Run**: give each comparison its own working directory and run the variation finder before the PAV step: `cd GifuT2T_GifuOld && bash GifuT2T_vs_GifuOld_variation_finder.sh && bash GifuT2T_vs_GifuOld_get_PAV_from_coords.sh` (the same for MG20); then `bash calc_bin_coverage.sh`, and `bash run.sh` or `python3 run.py` in `Gifu_new/` and `MG20_new/`.
**Tools**: MUMmer4 4.0.1 (`nucmer -c 1000 --mum --maxgap=1000`, `delta-filter -i 90 -l 1000 -q`, `show-coords -THrd`, `dnadiff -d -p`), SyRI 1.7.1, plotsr 1.1.1, samtools 1.21, bedtools 2.31.1 (`merge -d 50`, `complement -g`, `intersect`, `subtract`), seqkit 2.10.1 (`grep -v -r -p`), GNU Awk 5.1.0, Python 3 with pandas and matplotlib, CSUB (`-q c01`, `-n 64`).

</details>

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>GifuT2T_GifuOld/GifuT2T_vs_GifuOld_get_PAV_from_coords.sh</code> · Shell · 90 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/06_genome_comparison/05_t2t_vs_previous/GifuT2T_GifuOld/GifuT2T_vs_GifuOld_get_PAV_from_coords.sh)

```bash
#!/usr/bin/env bash
set -euo pipefail

############################################
# 输入文件（请按需修改）
############################################
COORDS=GifuT2T_GifuOld.filtered.coords

FA_T2T=Lotus_GifuT2T_v1.0.fasta
FA_OLD=Lotus_Gifu_Old1.fa

PREFIX=GifuT2T_vs_GifuOld

############################################
# 0. 准备 genome 文件
############################################
samtools faidx "$FA_T2T"
samtools faidx "$FA_OLD"

cut -f1,2 ${FA_T2T}.fai > GifuT2T.genome
cut -f1,2 ${FA_OLD}.fai > GifuOld.genome

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
# 2. 提取 GifuOld（ref）比对区段
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
> ${PREFIX}.GifuOld.aligned.bed

############################################
# 3. 反选（complement）得到 PAV
############################################
bedtools complement \
  -i ${PREFIX}.GifuT2T.aligned.bed \
  -g GifuT2T.genome \
> ${PREFIX}.GifuT2T.PAV.bed

bedtools complement \
  -i ${PREFIX}.GifuOld.aligned.bed \
  -g GifuOld.genome \
> ${PREFIX}.GifuOld.PAV.bed

############################################
# 4. 统计区间数量 & 总长度
############################################

for f in ${PREFIX}.GifuT2T.PAV.bed ${PREFIX}.GifuOld.PAV.bed
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
<summary><code>GifuT2T_GifuOld/GifuT2T_vs_GifuOld_variation_finder.sh</code> · Shell · 71 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/06_genome_comparison/05_t2t_vs_previous/GifuT2T_GifuOld/GifuT2T_vs_GifuOld_variation_finder.sh)

```bash
#!/bin/bash
#CSUB -J variation
#CSUB -q c01
#CSUB -o variation.out
#CSUB -e variation.error
#CSUB -n 64
#CSUB -R span[hosts=1]

# 0. 去除原本基因组中的多余的染色体
seqkit grep -v -r -p "LjG1.1_chloro" -p "LjG1.1_chr0" -p "LjG1.1_mito" \
 Lotus_Gifu_Old.fa > Lotus_Gifu_Old1.fa

# 1. 全基因组比对（根据需要调整参数）
# 使用nucmer进行比对，参数可根据前述讨论选择 --mum 或 --maxmatch
# 这里示例采用更全面的--maxmatch并设置较高阈值避免噪音
nucmer -c 1000 --mum --maxgap=1000 --prefix=GifuT2T_GifuOld Lotus_Gifu_Old1.fa Lotus_GifuT2T_v1.0.fasta

# 2. 过滤比对结果，提高比对质量
delta-filter -i 90 -l 1000 -q GifuT2T_GifuOld.delta > GifuT2T_GifuOld.filtered.delta

# 可选：输出比对统计和绘图（非必要步骤，仅检查比对情况）
dnadiff -d GifuT2T_GifuOld.filtered.delta -p GifuT2T_GifuOld.filtered
# show-coords 可用于生成比对坐标表格
show-coords -THrd GifuT2T_GifuOld.filtered.delta > GifuT2T_GifuOld.filtered.coords

# 3. 使用SyRI识别变异
syri -c GifuT2T_GifuOld.filtered.coords -d GifuT2T_GifuOld.filtered.delta \
    -r Lotus_Gifu_Old1.fa -q Lotus_GifuT2T_v1.0.fasta

# 4. 按类别提取变异坐标并保存为BED文件
# 小变异类别：SNP、插入、缺失
awk '$11=="SNP" {print $1"\t"$2"\t"$3 >> "Old_SNP.bed"; print $6"\t"$7"\t"$8 >> "T2T_SNP.bed"}' syri.out
awk '$11=="INS" {print $1"\t"$2"\t"$3 >> "Old_INS.bed"; print $6"\t"$7"\t"$8 >> "T2T_INS.bed"}' syri.out
awk '$11=="DEL" {print $1"\t"$2"\t"$3 >> "Old_DEL.bed"; print $6"\t"$7"\t"$8 >> "T2T_DEL.bed"}' syri.out

# 结构变异类别：倒位、易位、倒置易位、重复、倒置重复
awk '$11=="INV" {print $1"\t"$2"\t"$3 >> "Old_INV.bed"; print $6"\t"$7"\t"$8 >> "T2T_INV.bed"}' syri.out
awk '$11=="TRANS" {print $1"\t"$2"\t"$3 >> "Old_TRANS.bed"; print $6"\t"$7"\t"$8 >> "T2T_TRANS.bed"}' syri.out
awk '$11=="INVTR" {print $1"\t"$2"\t"$3 >> "Old_INVTR.bed"; print $6"\t"$7"\t"$8 >> "T2T_INVTR.bed"}' syri.out
awk '$11=="DUP" {print $1"\t"$2"\t"$3 >> "Old_DUP.bed"; print $6"\t"$7"\t"$8 >> "T2T_DUP.bed"}' syri.out
awk '$11=="INVDP" {print $1"\t"$2"\t"$3 >> "Old_INVDP.bed"; print $6"\t"$7"\t"$8 >> "T2T_INVDP.bed"}' syri.out

# 存在/缺失变异（PAV）：参考或查询中特有的大片段（INS/DEL/NOTAL >50bp）
awk '($11=="INS" || $11=="DEL" || $11=="NOTAL") && (($3-$2+1>50) || ($8-$7+1>50)) \
    {print $1"\t"$2"\t"$3 >> "Old_PAV.bed"; print $6"\t"$7"\t"$8 >> "T2T_PAV.bed"}' syri.out

# 5. 统计各类别变异事件数量和总长度
echo -e "Category\tCount\tTotal_length_in_Old\tTotal_length_in_T2T"
for type in SNP INS DEL INV TRANS INVTR DUP INVDP PAV; do
    # 统计事件数（行数）
    count=$(wc -l < "Old_${type}.bed")
    # 计算参考和查询总长度（bp）
    total_ref=$(awk '{sum += $3-$2+1} END{print sum}' "Old_${type}.bed")
    total_qry=$(awk '{sum += $3-$2+1} END{print sum}' "T2T_${type}.bed")
    echo -e "${type}\t${count}\t${total_ref}\t${total_qry}"
done

#基于syri的结构变异结果
#plotsr --sr syri.out --genomes genomes.txt -o GifuT2T_GifuOld_findSV.pdf

#-s 250000 过滤掉 <250 kb 的小型 rearrangements
#-R 会把没有被结构重排打断的相邻 syntenic blocks 合并起来
plotsr \
    --sr syri.out \
    --genomes genomes.txt \
    -s 250000 \
    -R \
    -W 5.3 \
    -H 5 \
    -f 8 \
    -o GifuT2T_GifuOld_findSV.pdf
```

</details>

<details open>
<summary><code>Gifu_new/run.py</code> · Python · 120 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/06_genome_comparison/05_t2t_vs_previous/Gifu_new/run.py)

```python
import matplotlib.pyplot as plt

def read_bed(file_path):
    """读取BED文件，返回列表[(chrom, start, end)]"""
    regions = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == "":
                continue
            chrom, start, end = line.strip().split()[:3]
            regions.append((chrom, int(start), int(end)))
    return regions

def calculate_overlap(region1, region2):
    """计算两个区域的重叠部分长度"""
    chrom1, start1, end1 = region1
    chrom2, start2, end2 = region2
    # 确保是同一个染色体
    if chrom1 != chrom2:
        return 0
    # 计算重叠区间
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    if overlap_start < overlap_end:
        return overlap_end - overlap_start
    return 0

def subtract_region(region1, region2):
    """从region1中去除region2的部分"""
    chrom1, start1, end1 = region1
    chrom2, start2, end2 = region2
    if chrom1 != chrom2:
        return [region1]
    result = []
    # 如果region1完全在region2之外，则直接返回
    if end1 <= start2 or start1 >= end2:
        return [region1]
    # 计算去除后的区域
    if start1 < start2:
        result.append((chrom1, start1, start2))
    if end1 > end2:
        result.append((chrom1, end2, end1))
    return result

def calculate_total_length(regions):
    """计算区域总长度"""
    return sum([end - start for _, start, end in regions])

# 读取文件
new_seq_bed = "gifuT2T_newseq.bed"
cent_bed = "cent.bed"
rDNA_bed = "rDNA.bed"
tel_bed = "gifu_tel.bed"
TE_bed = "repeat2.bed"

new_seq_regions = read_bed(new_seq_bed)
cent_regions = read_bed(cent_bed)
rDNA_regions = read_bed(rDNA_bed)
tel_regions = read_bed(tel_bed)
TE_regions = read_bed(TE_bed)

# 计算新序列与各个区域的重叠部分的长度
cent_overlap = sum([calculate_overlap(seq, cent) for seq in new_seq_regions for cent in cent_regions])
rDNA_overlap = sum([calculate_overlap(seq, rDNA) for seq in new_seq_regions for rDNA in rDNA_regions])
tel_overlap = sum([calculate_overlap(seq, tel) for seq in new_seq_regions for tel in tel_regions])
TE_overlap = sum([calculate_overlap(seq, TE) for seq in new_seq_regions for TE in TE_regions])

# 计算新序列的总长度
new_seq_length = calculate_total_length(new_seq_regions)

# 逐步去除每个区域的交集部分，并计算剩余的部分
remaining_after_cent = new_seq_regions
for cent in cent_regions:
    remaining_after_cent = [seq for seq in remaining_after_cent if calculate_overlap(seq, cent) == 0]

remaining_after_rDNA = remaining_after_cent
for rDNA in rDNA_regions:
    remaining_after_rDNA = [seq for seq in remaining_after_rDNA if calculate_overlap(seq, rDNA) == 0]

remaining_after_tel = remaining_after_rDNA
for tel in tel_regions:
    remaining_after_tel = [seq for seq in remaining_after_tel if calculate_overlap(seq, tel) == 0]

remaining_after_TE = remaining_after_tel
for TE in TE_regions:
    remaining_after_TE = [seq for seq in remaining_after_TE if calculate_overlap(seq, TE) == 0]

# 计算其他区域的长度
other_overlap = calculate_total_length(remaining_after_TE)

# 计算占比
cent_ratio = cent_overlap / new_seq_length
rDNA_ratio = rDNA_overlap / new_seq_length
tel_ratio = tel_overlap / new_seq_length
TE_ratio = TE_overlap / new_seq_length
other_ratio = other_overlap / new_seq_length

# 输出各个区域的交集和占比
print(f"Centromere overlap: {cent_overlap} bp, Centromere ratio: {cent_ratio:.4f}")
print(f"rDNA overlap: {rDNA_overlap} bp, rDNA ratio: {rDNA_ratio:.4f}")
print(f"Telomere overlap: {tel_overlap} bp, Telomere ratio: {tel_ratio:.4f}")
print(f"TE overlap (unique): {TE_overlap} bp, TE ratio: {TE_ratio:.4f}")
print(f"Other regions overlap: {other_overlap} bp, Other regions ratio: {other_ratio:.4f}")

# 使用Python绘制饼图并保存为PDF文件
labels = ['Centromere', 'rDNA', 'Telomere', 'Transposable Elements', 'Other Regions']
sizes = [cent_overlap, rDNA_overlap, tel_overlap, TE_overlap, other_overlap]
colors = ['#ff9999','#66b3ff','#99ff99','#ffcc99','#c2c2f0']

# 绘制饼图
plt.figure(figsize=(7, 7))
plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140, colors=colors)
plt.axis('equal')  # 保证饼图是圆形的
plt.title('Proportion of New Sequences in Different Regions')

# 保存为PDF文件
plt.savefig('new_sequences_region_proportions.pdf', format='pdf')

# 显示图形
plt.show()
```

</details>

<details open>
<summary><code>Gifu_new/run.sh</code> · Shell · 52 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/06_genome_comparison/05_t2t_vs_previous/Gifu_new/run.sh)

```bash
#!/bin/bash

# 设置文件路径
new_seq_bed="gifuT2T_newseq.bed"
cent_bed="cent.bed"
rDNA_bed="rDNA.bed"
tel_bed="gifu_tel.bed"
TE_bed="repeat2.bed"

# 计算新序列与每个区域的重叠部分（去重后的TE区域）
cent_overlap=$(bedtools intersect -a $new_seq_bed -b $cent_bed -wa | awk '{sum+=$3-$2} END {print sum}')
rDNA_overlap=$(bedtools intersect -a $new_seq_bed -b $rDNA_bed -wa | awk '{sum+=$3-$2} END {print sum}')
tel_overlap=$(bedtools intersect -a $new_seq_bed -b $tel_bed -wa | awk '{sum+=$3-$2} END {print sum}')

# 计算新序列与TE的重叠区域并去重
TE_overlap_raw=$(bedtools intersect -a $new_seq_bed -b $TE_bed -wa)
TE_overlap_unique=$(echo "$TE_overlap_raw" | bedtools sort | bedtools merge -i - | awk '{sum+=$3-$2} END {print sum}')

# 计算新序列的总长度
new_seq_length=$(awk '{sum+=$3-$2} END {print sum}' $new_seq_bed)

# 计算其他区域的长度（不在cent, rDNA, tel, TE内的区域）
other_overlap=$(($new_seq_length - $cent_overlap - $rDNA_overlap - $tel_overlap - $TE_overlap_unique))

# 输出占比
echo "Centromere overlap: $cent_overlap bp"
echo "rDNA overlap: $rDNA_overlap bp"
echo "Telomere overlap: $tel_overlap bp"
echo "TE overlap (unique): $TE_overlap_unique bp"
echo "Other regions overlap: $other_overlap bp"

# 使用Python绘制饼图并保存为PDF文件
python3 - <<EOF
import matplotlib.pyplot as plt

# 占比数据
labels = ['Centromere', 'rDNA', 'Telomere', 'Transposable Elements', 'Other Regions']
sizes = [$cent_overlap, $rDNA_overlap, $tel_overlap, $TE_overlap_unique, $other_overlap]
colors = ['#ff9999','#66b3ff','#99ff99','#ffcc99','#c2c2f0']

# 绘制饼图
plt.figure(figsize=(7, 7))
plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140, colors=colors)
plt.axis('equal')  # 保证饼图是圆形的
plt.title('Proportion of New Sequences in Different Regions')

# 保存为PDF文件
plt.savefig('new_sequences_region_proportions.pdf', format='pdf')

# 显示图形
plt.show()
EOF
```

</details>

<details open>
<summary><code>MG20T2T_MG20Old/MG20T2T_MG20Old_get_PAV_from_coords.sh</code> · Shell · 90 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/06_genome_comparison/05_t2t_vs_previous/MG20T2T_MG20Old/MG20T2T_MG20Old_get_PAV_from_coords.sh)

```bash
#!/usr/bin/env bash
set -euo pipefail

############################################
# 输入文件（请按需修改）
############################################
COORDS=MG20T2T_MG20Old.filtered.coords

FA_T2T=Lotus_MG20T2T_v1.0.fasta
FA_OLD=Lotus_MG20_Old1.fa

PREFIX=MG20T2T_vs_MG20Old

############################################
# 0. 准备 genome 文件
############################################
samtools faidx "$FA_T2T"
samtools faidx "$FA_OLD"

cut -f1,2 ${FA_T2T}.fai > MG20T2T.genome
cut -f1,2 ${FA_OLD}.fai > MG20Old.genome

############################################
# 1. 提取 MG20T2T（query）比对区段
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
> ${PREFIX}.MG20T2T.aligned.bed

############################################
# 2. 提取 MG20Old（ref）比对区段
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
> ${PREFIX}.MG20Old.aligned.bed

############################################
# 3. 反选（complement）得到 PAV
############################################
bedtools complement \
  -i ${PREFIX}.MG20T2T.aligned.bed \
  -g MG20T2T.genome \
> ${PREFIX}.MG20T2T.PAV.bed

bedtools complement \
  -i ${PREFIX}.MG20Old.aligned.bed \
  -g MG20Old.genome \
> ${PREFIX}.MG20Old.PAV.bed

############################################
# 4. 统计区间数量 & 总长度
############################################

for f in ${PREFIX}.MG20T2T.PAV.bed ${PREFIX}.MG20Old.PAV.bed
do

    echo "File: $f" >> ${PREFIX}.PAV.summary

    awk '
    {
        len = $3 - $2
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
<summary><code>MG20T2T_MG20Old/MG20T2T_MG20Old_variation_finder.sh</code> · Shell · 69 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/06_genome_comparison/05_t2t_vs_previous/MG20T2T_MG20Old/MG20T2T_MG20Old_variation_finder.sh)

```bash
#!/bin/bash
#CSUB -J variation
#CSUB -q c01
#CSUB -o variation.out
#CSUB -e variation.error
#CSUB -n 64
#CSUB -R span[hosts=1]

# 0.去除原本基因组中的多余的染色体
seqkit grep -v -r -p "lotja.MG20.gnm3.Lj0" -p "lotja.MG20.gnm3.Ljmito" -p "lotja.MG20.gnm3.Ljchloro" \
    -p "lotja.MG20.gnm3.LjM" -p "lotja.MG20.gnm3.LjC" \
    Lotus_MG20_Old.fa > Lotus_MG20_Old1.fa

# 1. 全基因组比对（根据需要调整参数）
# 使用nucmer进行比对，参数可根据前述讨论选择 --mum 或 --maxmatch
# 这里示例采用更全面的--maxmatch并设置较高阈值避免噪音
nucmer -c 1000 --mum --maxgap=1000 --prefix=MG20T2T_MG20Old Lotus_MG20_Old1.fa Lotus_MG20T2T_v1.0.fasta

# 2. 过滤比对结果，提高比对质量
delta-filter -i 90 -l 1000 -q MG20T2T_MG20Old.delta > MG20T2T_MG20Old.filtered.delta

# 可选：输出比对统计和绘图（非必要步骤，仅检查比对情况）
dnadiff -d MG20T2T_MG20Old.filtered.delta -p MG20T2T_MG20Old.filtered
# show-coords 可用于生成比对坐标表格
show-coords -THrd MG20T2T_MG20Old.filtered.delta > MG20T2T_MG20Old.filtered.coords

# 3. 使用SyRI识别变异
syri -c MG20T2T_MG20Old.filtered.coords -d MG20T2T_MG20Old.filtered.delta \
    -r Lotus_MG20_Old1.fa -q Lotus_MG20T2T_v1.0.fasta

# 4. 按类别提取变异坐标并保存为BED文件
# 小变异类别：SNP、插入、缺失
awk '$11=="SNP" {print $1"\t"$2"\t"$3 >> "Old_SNP.bed"; print $6"\t"$7"\t"$8 >> "T2T_SNP.bed"}' syri.out
awk '$11=="INS" {print $1"\t"$2"\t"$3 >> "Old_INS.bed"; print $6"\t"$7"\t"$8 >> "T2T_INS.bed"}' syri.out
awk '$11=="DEL" {print $1"\t"$2"\t"$3 >> "Old_DEL.bed"; print $6"\t"$7"\t"$8 >> "T2T_DEL.bed"}' syri.out

# 结构变异类别：倒位、易位、倒置易位、重复、倒置重复
awk '$11=="INV" {print $1"\t"$2"\t"$3 >> "Old_INV.bed"; print $6"\t"$7"\t"$8 >> "T2T_INV.bed"}' syri.out
awk '$11=="TRANS" {print $1"\t"$2"\t"$3 >> "Old_TRANS.bed"; print $6"\t"$7"\t"$8 >> "T2T_TRANS.bed"}' syri.out
awk '$11=="INVTR" {print $1"\t"$2"\t"$3 >> "Old_INVTR.bed"; print $6"\t"$7"\t"$8 >> "T2T_INVTR.bed"}' syri.out
awk '$11=="DUP" {print $1"\t"$2"\t"$3 >> "Old_DUP.bed"; print $6"\t"$7"\t"$8 >> "T2T_DUP.bed"}' syri.out
awk '$11=="INVDP" {print $1"\t"$2"\t"$3 >> "Old_INVDP.bed"; print $6"\t"$7"\t"$8 >> "T2T_INVDP.bed"}' syri.out

# 存在/缺失变异（PAV）：参考或查询中特有的大片段（INS/DEL/NOTAL >50bp）
awk '($11=="INS" || $11=="DEL" || $11=="NOTAL") && (($3-$2+1>50) || ($8-$7+1>50)) \
    {print $1"\t"$2"\t"$3 >> "Old_PAV.bed"; print $6"\t"$7"\t"$8 >> "T2T_PAV.bed"}' syri.out

# 5. 统计各类别变异事件数量和总长度
echo -e "Category\tCount\tTotal_length_in_Old\tTotal_length_in_T2T"
for type in SNP INS DEL INV TRANS INVTR DUP INVDP PAV; do
    # 统计事件数（行数）
    count=$(wc -l < "Old_${type}.bed")
    # 计算参考和查询总长度（bp）
    total_ref=$(awk '{sum += $3-$2+1} END{print sum}' "Old_${type}.bed")
    total_qry=$(awk '{sum += $3-$2+1} END{print sum}' "T2T_${type}.bed")
    echo -e "${type}\t${count}\t${total_ref}\t${total_qry}"
done

#-s 250000 过滤掉 <250 kb 的小型 rearrangements
#-R 会把没有被结构重排打断的相邻 syntenic blocks 合并起来
plotsr \
    --sr syri.out \
    --genomes genomes.txt \
    -s 250000 \
    -R \
    -W 5.3 \
    -H 5 \
    -f 8 \
    -o MG20T2T_MG20Old_findSV.pdf
```

</details>

<details open>
<summary><code>MG20_new/run.py</code> · Python · 120 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/06_genome_comparison/05_t2t_vs_previous/MG20_new/run.py)

```python
import matplotlib.pyplot as plt

def read_bed(file_path):
    """读取BED文件，返回列表[(chrom, start, end)]"""
    regions = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == "":
                continue
            chrom, start, end = line.strip().split()[:3]
            regions.append((chrom, int(start), int(end)))
    return regions

def calculate_overlap(region1, region2):
    """计算两个区域的重叠部分长度"""
    chrom1, start1, end1 = region1
    chrom2, start2, end2 = region2
    # 确保是同一个染色体
    if chrom1 != chrom2:
        return 0
    # 计算重叠区间
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    if overlap_start < overlap_end:
        return overlap_end - overlap_start
    return 0

def subtract_region(region1, region2):
    """从region1中去除region2的部分"""
    chrom1, start1, end1 = region1
    chrom2, start2, end2 = region2
    if chrom1 != chrom2:
        return [region1]
    result = []
    # 如果region1完全在region2之外，则直接返回
    if end1 <= start2 or start1 >= end2:
        return [region1]
    # 计算去除后的区域
    if start1 < start2:
        result.append((chrom1, start1, start2))
    if end1 > end2:
        result.append((chrom1, end2, end1))
    return result

def calculate_total_length(regions):
    """计算区域总长度"""
    return sum([end - start for _, start, end in regions])

# 读取文件
new_seq_bed = "mg20T2T_newseq.bed"
cent_bed = "cent.bed"
rDNA_bed = "rDNA.bed"
tel_bed = "tel.bed"
TE_bed = "repeat2.bed"

new_seq_regions = read_bed(new_seq_bed)
cent_regions = read_bed(cent_bed)
rDNA_regions = read_bed(rDNA_bed)
tel_regions = read_bed(tel_bed)
TE_regions = read_bed(TE_bed)

# 计算新序列与各个区域的重叠部分的长度
cent_overlap = sum([calculate_overlap(seq, cent) for seq in new_seq_regions for cent in cent_regions])
rDNA_overlap = sum([calculate_overlap(seq, rDNA) for seq in new_seq_regions for rDNA in rDNA_regions])
tel_overlap = sum([calculate_overlap(seq, tel) for seq in new_seq_regions for tel in tel_regions])
TE_overlap = sum([calculate_overlap(seq, TE) for seq in new_seq_regions for TE in TE_regions])

# 计算新序列的总长度
new_seq_length = calculate_total_length(new_seq_regions)

# 逐步去除每个区域的交集部分，并计算剩余的部分
remaining_after_cent = new_seq_regions
for cent in cent_regions:
    remaining_after_cent = [seq for seq in remaining_after_cent if calculate_overlap(seq, cent) == 0]

remaining_after_rDNA = remaining_after_cent
for rDNA in rDNA_regions:
    remaining_after_rDNA = [seq for seq in remaining_after_rDNA if calculate_overlap(seq, rDNA) == 0]

remaining_after_tel = remaining_after_rDNA
for tel in tel_regions:
    remaining_after_tel = [seq for seq in remaining_after_tel if calculate_overlap(seq, tel) == 0]

remaining_after_TE = remaining_after_tel
for TE in TE_regions:
    remaining_after_TE = [seq for seq in remaining_after_TE if calculate_overlap(seq, TE) == 0]

# 计算其他区域的长度
other_overlap = calculate_total_length(remaining_after_TE)

# 计算占比
cent_ratio = cent_overlap / new_seq_length
rDNA_ratio = rDNA_overlap / new_seq_length
tel_ratio = tel_overlap / new_seq_length
TE_ratio = TE_overlap / new_seq_length
other_ratio = other_overlap / new_seq_length

# 输出各个区域的交集和占比
print(f"Centromere overlap: {cent_overlap} bp, Centromere ratio: {cent_ratio:.4f}")
print(f"rDNA overlap: {rDNA_overlap} bp, rDNA ratio: {rDNA_ratio:.4f}")
print(f"Telomere overlap: {tel_overlap} bp, Telomere ratio: {tel_ratio:.4f}")
print(f"TE overlap (unique): {TE_overlap} bp, TE ratio: {TE_ratio:.4f}")
print(f"Other regions overlap: {other_overlap} bp, Other regions ratio: {other_ratio:.4f}")

# 使用Python绘制饼图并保存为PDF文件
labels = ['Centromere', 'rDNA', 'Telomere', 'Transposable Elements', 'Other Regions']
sizes = [cent_overlap, rDNA_overlap, tel_overlap, TE_overlap, other_overlap]
colors = ['#ff9999','#66b3ff','#99ff99','#ffcc99','#c2c2f0']

# 绘制饼图
plt.figure(figsize=(7, 7))
plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140, colors=colors)
plt.axis('equal')  # 保证饼图是圆形的
plt.title('Proportion of New Sequences in Different Regions')

# 保存为PDF文件
plt.savefig('new_sequences_region_proportions.pdf', format='pdf')

# 显示图形
plt.show()
```

</details>

<details open>
<summary><code>MG20_new/run.sh</code> · Shell · 85 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/06_genome_comparison/05_t2t_vs_previous/MG20_new/run.sh)

```bash
#!/bin/bash

# 设置文件路径
new_seq_bed="mg20T2T_newseq.bed"
cent_bed="cent.bed"
rDNA_bed="rDNA.bed"
tel_bed="tel.bed"
TE_bed="repeat2.bed"

# 检查文件路径是否正确
if [[ ! -f "$new_seq_bed" || ! -f "$cent_bed" || ! -f "$rDNA_bed" || ! -f "$tel_bed" || ! -f "$TE_bed" ]]; then
  echo "Error: One or more input files do not exist."
  exit 1
fi

# 计算新序列与cent的交集的真实长度
cent_overlap=$(bedtools intersect -a $new_seq_bed -b $cent_bed -wb | awk '{sum+=$3-$2} END {print sum}')
new_seq_length=$(awk '{sum+=$3-$2} END {print sum}' $new_seq_bed)
cent_ratio=$(echo "scale=4; $cent_overlap / $new_seq_length" | bc)

# 去除cent部分后，计算剩余新序列与rDNA的交集的真实长度
remaining_after_cent=$(bedtools subtract -a $new_seq_bed -b $cent_bed)
rDNA_overlap=$(bedtools intersect -a $remaining_after_cent -b $rDNA_bed -wb | awk '{sum+=$3-$2} END {print sum}')
# 如果rDNA没有交集，则设置为0
if [ -z "$rDNA_overlap" ]; then
    rDNA_overlap=0
fi
rDNA_ratio=$(echo "scale=4; $rDNA_overlap / $new_seq_length" | bc)

# 去除rDNA部分后，计算剩余新序列与telomere的交集的真实长度
remaining_after_rDNA=$(bedtools subtract -a $remaining_after_cent -b $rDNA_bed)
tel_overlap=$(bedtools intersect -a $remaining_after_rDNA -b $tel_bed -wb | awk '{sum+=$3-$2} END {print sum}')
# 如果tel没有交集，则设置为0
if [ -z "$tel_overlap" ]; then
    tel_overlap=0
fi
tel_ratio=$(echo "scale=4; $tel_overlap / $new_seq_length" | bc)

# 去除tel部分后，计算剩余新序列与TE的交集的真实长度
remaining_after_tel=$(bedtools subtract -a $remaining_after_rDNA -b $tel_bed)
TE_overlap_raw=$(bedtools intersect -a $remaining_after_tel -b $TE_bed -wb)
TE_overlap_unique=$(echo "$TE_overlap_raw" | bedtools sort | bedtools merge -i - | awk '{sum+=$3-$2} END {print sum}')
# 如果TE没有交集，则设置为0
if [ -z "$TE_overlap_unique" ]; then
    TE_overlap_unique=0
fi
TE_ratio=$(echo "scale=4; $TE_overlap_unique / $new_seq_length" | bc)

# 最终去除TE后的区域即为Other regions
remaining_after_TE=$(bedtools subtract -a $remaining_after_tel -b $TE_bed)
other_overlap=$(echo "$remaining_after_TE" | awk '{sum+=$3-$2} END {print sum}')
# 如果Other regions没有交集，则设置为0
if [ -z "$other_overlap" ]; then
    other_overlap=0
fi

# 输出各个区域的交集和占比
echo "Centromere overlap: $cent_overlap bp, Centromere ratio: $cent_ratio"
echo "rDNA overlap: $rDNA_overlap bp, rDNA ratio: $rDNA_ratio"
echo "Telomere overlap: $tel_overlap bp, Telomere ratio: $tel_ratio"
echo "TE overlap (unique): $TE_overlap_unique bp, TE ratio: $TE_ratio"
echo "Other regions overlap: $other_overlap bp"

# 使用Python绘制饼图并保存为PDF文件
python3 - <<EOF
import matplotlib.pyplot as plt

# 占比数据
# 如果有空值，将其设置为0
sizes = [$cent_overlap, $rDNA_overlap, $tel_overlap, $TE_overlap_unique, $other_overlap]
labels = ['Centromere', 'rDNA', 'Telomere', 'Transposable Elements', 'Other Regions']
colors = ['#ff9999','#66b3ff','#99ff99','#ffcc99','#c2c2f0']

# 绘制饼图
plt.figure(figsize=(7, 7))
plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140, colors=colors)
plt.axis('equal')  # 保证饼图是圆形的
plt.title('Proportion of New Sequences in Different Regions')

# 保存为PDF文件
plt.savefig('new_sequences_region_proportions.pdf', format='pdf')

# 显示图形
plt.show()
EOF
```

</details>

<details open>
<summary><code>calc_bin_coverage.sh</code> · Shell · 48 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/06_genome_comparison/05_t2t_vs_previous/calc_bin_coverage.sh)

```bash
#!/usr/bin/env bash
set -euo pipefail

FEATURE_BED="gifuT2T_newseq.bed"
BIN_BED="{PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Gifu_100k.bed"
OUT_TSV="gifu_100kb_coverage.tsv"

# 1) 规范化：只取前三列，确保start<end，排序
awk 'BEGIN{OFS="\t"}
     NF>=3 && $2<$3 {print $1,$2,$3}' "$FEATURE_BED" \
  | sort -k1,1 -k2,2n > features.sorted.bed

awk 'BEGIN{OFS="\t"}
     NF>=3 && $2<$3 {print $1,$2,$3}' "$BIN_BED" \
  | sort -k1,1 -k2,2n > bins100k.sorted.bed

# 2) 计算 features 与 bin 的重叠长度（-wo 最后列是 overlap bp）
# 输出列：
# bin_chr bin_start bin_end feat_chr feat_start feat_end overlap
bedtools intersect -a bins100k.sorted.bed -b features.sorted.bed -wo \
  | awk 'BEGIN{OFS="\t"}
         {
           key=$1"\t"$2"\t"$3;
           ov[key]+=$7;
         }
         END{
           for(k in ov){
             print k, ov[k];
           }
         }' \
  | sort -k1,1 -k2,2n > bin_overlap_bp.tsv

# 3) 把没有任何 overlap 的 bin 补 0，并计算比例 = overlap / bin_len
# 输出：chr start end bin_len overlap_bp proportion
awk 'BEGIN{OFS="\t"}
     FNR==NR{
       ov[$1"\t"$2"\t"$3]=$4;
       next
     }
     {
       key=$1"\t"$2"\t"$3;
       len=$3-$2;
       bp=(key in ov)?ov[key]:0;
       prop=(len>0)?bp/len:0;
       print $1,$2,$3,len,bp,prop
     }' bin_overlap_bp.tsv bins100k.sorted.bed > "$OUT_TSV"

echo "Done: $OUT_TSV"
```

</details>
