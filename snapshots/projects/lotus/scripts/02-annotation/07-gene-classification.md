---
title: '基因分类'
sidebar:
  label: '基因分类'
  order: 7
---


## 中文

**作用**：`new_genes_finder.sh` 把更新后的 T2T 基因模型与上一版注释比对分类：在旧蛋白组中无合格 BLASTP 匹配的 T2T 基因，拆成旧组装未覆盖区域的基因（`new_genes_Gap`）、旧序列存在但旧模型需修正的基因（`genes_anno_correct`）和其余基因（`new_genes_anno`）。
**入口**：`new_genes_find.sh`（内部调用 `new_genes_finder.sh`）。
**输入**：7 个位置参数，1–3 为 T2T 集（脚本中的 ref），4–6 为上一版注释（query）：`pep1` T2T 蛋白（BLASTP 库）、`genome1` T2T 基因组（BLASTN 查询来源）、`gff1` T2T GFF3；`pep2` 旧版蛋白、`genome2` 旧版基因组、`gff2` 旧版 GFF3；`gap` 旧组装未覆盖区域 BED，即使该 BED 为空也必须给出。
**输出**：`all_genes.list`、`all_new_genes.list`/`.bed`、`new_genes_none_Gap.list`/`.bed`/`.fa`；分类 `new_genes_Gap/`、`genes_anno_correct/`、`new_genes_anno/` 下的 `.list`/`.bed`/`.tsv`；中间文件 `blastp/`、`blastn/`、`best_hit/` 与 `ref_genes.bed`、`query_genes.bed`；计数：Gifu v1.3 修正 2,104 / 新注释 3,033 / 补 gap 22 = 5,159 个更新基因，MG20 v3.0 1,256 / 2,491 / 5,548 = 9,295。
**运行**：`bash new_genes_find.sh`；核心调用 `sh new_genes_finder.sh Lotus_GifuT2T.pep.fa Lotus_GifuT2T_v1.0.fasta Lotus_GifuT2T.gff3 Lotus_Gifu_Old.pep.fa Lotus_Gifu_Old.fa Lotus_Gifu_Old.gff3 Old_region.uncovered.bed`（MG20 换成 `Lotus_MG20T2T.*`/`Lotus_MG20_Old.*`）；`blast_merge.py` 以仓库外绝对路径 ``{SOFTWARE_HX}`/Script/Homologous_gene_search/blast_merge.py` 调用，须改为本仓库 `../10_gene_correspondence/blast_merge.py`。
**工具**：BLAST+ `makeblastdb`/`blastp`/`blastn`（`-dbtype prot`/`nucl`、`-evalue 1e-5`、`-outfmt "6 ... qcovs"`、`-num_threads 88`）；Python 3 + pandas（`blast_merge.py`）；seqkit `fx2tab -l -n -i`；bedtools `intersect`、`getfasta -fi -bed -name -s`；gawk 取 BLASTN 最佳命中；过滤式 `awk '$3>50 && $7>0.5 && $8>0.5'`。

<details>
<summary>English</summary>

**Purpose**: `new_genes_finder.sh` classifies the updated T2T gene models against the previous annotation: T2T genes without an acceptable BLASTP match in the old proteome are split into genes in regions the old assembly did not cover (`new_genes_Gap`), genes whose old sequence is present but whose old model needs correcting (`genes_anno_correct`), and the remainder (`new_genes_anno`).
**Entry point**: `new_genes_find.sh` (which calls `new_genes_finder.sh`).
**Inputs**: seven positional arguments; 1–3 are the T2T set (the script's ref) and 4–6 the previous annotation (its query): `pep1` T2T proteins (BLASTP database), `genome1` T2T genome (source of the BLASTN queries), `gff1` T2T GFF3; `pep2` previous-version proteins, `genome2` previous-version genome, `gff2` previous-version GFF3; `gap` BED of regions the old assembly did not cover, which must be given even when the BED is empty.
**Outputs**: `all_genes.list`, `all_new_genes.list`/`.bed`, `new_genes_none_Gap.list`/`.bed`/`.fa`; the `.list`/`.bed`/`.tsv` files under `new_genes_Gap/`, `genes_anno_correct/` and `new_genes_anno/`; intermediates `blastp/`, `blastn/`, `best_hit/`, `ref_genes.bed`, `query_genes.bed`; totals: Gifu v1.3 2,104 corrected / 3,033 newly annotated / 22 gap-filling = 5,159 updated genes, MG20 v3.0 1,256 / 2,491 / 5,548 = 9,295.
**Run**: `bash new_genes_find.sh`; the core call is `sh new_genes_finder.sh Lotus_GifuT2T.pep.fa Lotus_GifuT2T_v1.0.fasta Lotus_GifuT2T.gff3 Lotus_Gifu_Old.pep.fa Lotus_Gifu_Old.fa Lotus_Gifu_Old.gff3 Old_region.uncovered.bed` (MG20 uses `Lotus_MG20T2T.*`/`Lotus_MG20_Old.*`); `blast_merge.py` is called through the absolute path ``{SOFTWARE_HX}`/Script/Homologous_gene_search/blast_merge.py`, which must be changed to this repository's `../10_gene_correspondence/blast_merge.py`.
**Tools**: BLAST+ `makeblastdb`/`blastp`/`blastn` (`-dbtype prot`/`nucl`, `-evalue 1e-5`, `-outfmt "6 ... qcovs"`, `-num_threads 88`); Python 3 + pandas (`blast_merge.py`); seqkit `fx2tab -l -n -i`; bedtools `intersect`, `getfasta -fi -bed -name -s`; gawk for the BLASTN best hit; the filter `awk '$3>50 && $7>0.5 && $8>0.5'`.

</details>

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>new_genes_find.sh</code> · Shell · 21 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/02_annotation/07_gene_classification/new_genes_find.sh)

```bash
#!/bin/bash
#CSUB -J FindNewGenes
#CSUB -q c01
#CSUB -o FindNewGenes.out
#CSUB -e FindNewGenes.error
#CSUB -n 64
#CSUB -R span[hosts=1]

sh new_genes_finder.sh Lotus_GifuT2T.pep.fa\
    Lotus_GifuT2T_v1.0.fasta Lotus_GifuT2T.gff3\
    Lotus_Gifu_Old.pep.fa\
    Lotus_Gifu_Old.fa \
    Lotus_Gifu_Old.gff3\
    Old_region.uncovered.bed

sh new_genes_finder.sh Lotus_MG20T2T.pep.fa\
    Lotus_MG20T2T_v1.0.fasta Lotus_MG20T2T.gff3\
    Lotus_MG20_Old.pep.fa\
    Lotus_MG20_Old.fa \
    Lotus_MG20_Old.gff3\
    Old_region.uncovered.bed
```

</details>

<details open>
<summary><code>new_genes_finder.sh</code> · Shell · 254 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/02_annotation/07_gene_classification/new_genes_finder.sh)

```bash
#!/bin/bash
set -euo pipefail

#=================== 参数设置 ===================
cul=$(pwd)
THREADS=88

#=================== 输入文件 ===================
pep1=$1
genome1=$2
gff1=$3
pep2=$4
genome2=$5
gff2=$6
gap=$7   # gap 可为空，但必须占位

#=================== 获取蛋白质序列一致性或覆盖率较低的基因列表 ===================
mkdir -p blastp

# 构建 BLASTP 数据库
makeblastdb -in "$pep1" -dbtype prot -out "blastp/db"

# 运行 BLASTP
blastp -query "$pep2" -db "blastp/db" \
    -out "blastp/result.blast" \
    -evalue 1e-5 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
    -num_threads $THREADS

cd blastp

cut -f1,2,4,3 result.blast > result.blast.res   # qid sid length pident

seqkit fx2tab -l -n -i "$cul/$pep1" > ref_protein_lengths.txt
seqkit fx2tab -l -n -i "$cul/$pep2" > query_protein_lengths.txt

python3 {SOFTWARE_HX}/Script/Homologous_gene_search/blast_merge.py \
    result.blast.res query_protein_lengths.txt ref_protein_lengths.txt blast.merge

# 筛选一致性差 或 覆盖率差的基因
awk -F'\t' 'BEGIN{OFS="\t"}
    $3=="gene" {
        match($9, /ID=([^;]+)/, a)
        if(a[1]!="") print a[1]
    }' "$cul/$gff1" > ../all_genes.list

awk -F "\t" '$3>50 && $7>0.5 && $8>0.5 {print $2}' blast.merge > ../blast_gene.tmp.list

cd "$cul"

grep -Fv -f blast_gene.tmp.list all_genes.list > all_new_genes.list

rm blast_gene.tmp.list

#=================== 提取基因所在位置 ===================
awk -F'\t' 'BEGIN{OFS="\t"}
    $3=="gene" {
        match($9, /ID=([^;]+)/, a)
        if(a[1]!="") print $1, $4-1, $5, a[1], $7
    }' "$gff1" > ref_genes.bed

awk -F'\t' 'BEGIN{OFS="\t"}
    $3=="gene" {
        match($9, /ID=([^;]+)/, a)
        if(a[1]!="") print $1, $4-1, $5, a[1], $7
    }' "$gff2" > query_genes.bed

grep -F -f all_new_genes.list ref_genes.bed > all_new_genes.bed

#=================== 分离 gap 区域基因 ===================
if [ -n "$gap" ] && [ -f "$gap" ]; then
    mkdir -p new_genes_Gap

    bedtools intersect \
        -a all_new_genes.bed \
        -b "$gap" \
        -wa \
        | awk '{print $4}' \
        | sort -u > new_genes_Gap/new_genes_Gap.list

    grep -F -f new_genes_Gap/new_genes_Gap.list ref_genes.bed \
        > new_genes_Gap/new_genes_Gap.bed

    grep -Fv -f new_genes_Gap/new_genes_Gap.list all_new_genes.list \
        > new_genes_none_Gap.list
else
    cp all_new_genes.list new_genes_none_Gap.list
fi

#=================== 获取非 gap 区的新注释基因序列 ===================
grep -F -f new_genes_none_Gap.list ref_genes.bed > new_genes_none_Gap.bed

bedtools getfasta \
    -fi "$genome1" \
    -bed new_genes_none_Gap.bed \
    -name -s \
    > new_genes_none_Gap.fa

#=================== BLASTN 到旧基因组 ===================
mkdir -p blastn

makeblastdb -in "$genome2" -dbtype nucl -out blastn/db

blastn -query new_genes_none_Gap.fa -db blastn/db \
    -out blastn/result.blast \
    -evalue 1e-5 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
    -num_threads $THREADS

#=================== 区分其余新基因的类型 ===================
echo "[INFO] 开始筛选 blastn 最优匹配..."

mkdir -p best_hit
mkdir -p genes_anno_correct
mkdir -p new_genes_anno

# ------------------------------
# 1. 从 qseqid 中提取染色体编号（Chr1 / Chr2 ...）
# ------------------------------
extract_chr() {
    echo "$1" | sed -E 's/.*Chr([0-9]+).*/\1/'
}

# ------------------------------
# 2. 按优先级挑选每个 qseqid 的最佳匹配
# evalue 最小 > pident 最大 > length 最大
# 并要求染色体一致
# ------------------------------
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$9,$10,$11}' blastn/result.blast > blastn/tmp.rearranged


gawk -v OFS="\t" '
{
    q=$1; 
    s=$2; 
    pident=$3; 
    len=$4; 
    start_pos=$5; 
    end_pos=$6; 
    evalue=$7;

    # 提取 qseqid 染色体 Chr1
    chr_q = q; 
    sub(/.*Chr/, "", chr_q); 
    sub(/:.*/, "", chr_q);

    # 提取 sseqid 染色体 chr1
    chr_s = s; 
    sub(/.*chr/, "", chr_s);

    if (chr_q != chr_s) next;

    key = q;

    if (!(key in best)) {
        best[key]=$0;
        evalue_min[key]=evalue;
        pident_max[key]=pident;
        len_max[key]=len;
    } else {
        if (evalue < evalue_min[key] ||
            (evalue == evalue_min[key] && pident > pident_max[key]) ||
            (evalue == evalue_min[key] && pident == pident_max[key] && len > len_max[key])) {

            best[key]=$0;
            evalue_min[key]=evalue;
            pident_max[key]=pident;
            len_max[key]=len;
        }
    }
}
END {
    for (k in best) print best[k];
}
' blastn/tmp.rearranged > best_hit/best_hit.raw

echo "[INFO] 提取最佳匹配完毕 → best_hit/best_hit.raw"

# ------------------------------
# 3. 只保留 qseqid, sseqid, sstart, send
# ------------------------------
awk -v OFS="\t" '{print $1, $2, $5, $6}' best_hit/best_hit.raw > best_hit/best_hit.simple

echo "[INFO] 最佳匹配简化文件 → best_hit/best_hit.simple"

# ------------------------------
# 4. 与 query_genes.bed 比对是否有重叠
# ------------------------------
# 将 best_hit.simple 转为 BED
awk -v OFS="\t" '
{
    chrom=$2;

    # sstart=$3; send=$4;
    sstart=$3; 
    send=$4;

    # BED 要求 start < end
    if (sstart <= send) {
        start = sstart - 1;
        end = send;
    } else {
        start = send - 1;
        end = sstart;
    }

    print chrom, start, end, $1;
}
' best_hit/best_hit.simple > best_hit/best_hit.bed

# 进行 overlap 检查
bedtools intersect -wa -wb \
    -a best_hit/best_hit.bed \
    -b query_genes.bed \
    > best_hit/overlap.txt

# ------------------------------
# 5. 输出分类结果
# ------------------------------
# 获取修正基因
# 有重叠 → 输出 qseqid 与 query gene
awk -v OFS="\t" '
{
    qseq=$4;
    query_gene=$8;
    print qseq, query_gene;
}
' best_hit/overlap.txt \
> genes_anno_correct/genes_anno_correct.tsv

# 清洗基因名称
awk -v OFS="\t" '{
    gsub(/::.*/, "", $1);
    print $1, $2;
}' genes_anno_correct/genes_anno_correct.tsv > genes_anno_correct/genes_anno_correct.cleaned.tsv \
    && mv genes_anno_correct/genes_anno_correct.cleaned.tsv genes_anno_correct/genes_anno_correct.tsv

cut -f1 genes_anno_correct/genes_anno_correct.tsv | sort -u > genes_anno_correct/genes_anno_correct.list

grep -F -f genes_anno_correct/genes_anno_correct.list ref_genes.bed \
    > genes_anno_correct/genes_anno_correct.bed

echo "[INFO] 修正的基因 → genes_anno_correct/genes_anno_correct.list"

# 获取新注释的基因
# 无重叠 → 输出 qseqid 列表
grep -Fv -f genes_anno_correct/genes_anno_correct.list new_genes_none_Gap.list > new_genes_anno/new_genes_anno.list

grep -F -f new_genes_anno/new_genes_anno.list ref_genes.bed \
    > new_genes_anno/new_genes_anno.bed

echo "[INFO] 新注释的基因 → new_genes_anno/new_genes_anno.list"

echo "[INFO] 任务完成！"
```

</details>
