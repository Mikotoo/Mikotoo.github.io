---
title: 'RNA-seq 比对'
sidebar:
  label: 'RNA-seq 比对'
  order: 2
---


## 中文

**作用**：用 StringTie 把 Gifu RNA-seq 文库分别对 T2T 组装注释与 BRAKER 过滤注释定量，得到下游表达证据。
**入口**：`01_align/align.sh`。
**输入**：``{PROJ_LOTUS_ZC}`/data/Gifu/RNA/clean/align.conf`（逐行 `sample fq1 fq2`）；`$dir/01.align/<sample>.bam`；`$ref1=`{PROJ_LOTUS_ZC}`/Gifu/00.assemble_result/Lotus_GifuT2T.gtf`、`$ref2=`{PROJ_LOTUS_ZC}`/Gifu/09.annotation/06.braker/Gifu_filtered.gtf`。
**输出**：`02.stringtie_geta/<sample>/<sample>.gtf`+`.tsv`、`03.stringtie_braker/<sample>/` 下同名文件；`02.tpm_geta/<sample>.tpm`、`03.tpm_braker/<sample>.tpm` 及各自 `all_samples.tpm.tsv`（gene_id + TPM，取 `.tsv` 第 9 列去掉表头后用 `join` 逐样本合并）。
**运行**：`bash 01_align/align.sh`；每样本 `stringtie -p 64 -G $ref1 -e -B -o $sample.gtf -A $sample.tsv $dir/01.align/$sample.bam`（`$ref2` 同样式）；索引构建与比对块被注释，须自备 `$dir/01.align/<sample>.bam`，且 `03.*braker` 分支依赖 `04_braker` 的输出。
**工具**：StringTie `-p 64 -G <gtf> -e -B -o <sample>.gtf -A <sample>.tsv`；HISAT2 `hisat2-build -p 64`、`hisat2 -p 64 -x <index> -1 <fq1> -2 <fq2> -S <sample>.sam --new-summary`；samtools `view -H`、`sort -@ 64`、`index`（比对块已注释）。

<details>
<summary>English</summary>

**Purpose**: Quantify the Gifu RNA-seq libraries with StringTie against the T2T assembly annotation and the BRAKER-filtered annotation, giving the expression evidence used downstream.
**Entry point**: `01_align/align.sh`.
**Inputs**: ``{PROJ_LOTUS_ZC}`/data/Gifu/RNA/clean/align.conf` (read line by line as `sample fq1 fq2`); `$dir/01.align/<sample>.bam`; `$ref1=`{PROJ_LOTUS_ZC}`/Gifu/00.assemble_result/Lotus_GifuT2T.gtf` and `$ref2=`{PROJ_LOTUS_ZC}`/Gifu/09.annotation/06.braker/Gifu_filtered.gtf`.
**Outputs**: `02.stringtie_geta/<sample>/<sample>.gtf`+`.tsv` and the same files under `03.stringtie_braker/<sample>/`; `02.tpm_geta/<sample>.tpm`, `03.tpm_braker/<sample>.tpm` and each `all_samples.tpm.tsv` (gene_id plus TPM, column 9 of the `.tsv` with the header stripped, merged across samples with `join`).
**Run**: `bash 01_align/align.sh`; per sample `stringtie -p 64 -G $ref1 -e -B -o $sample.gtf -A $sample.tsv $dir/01.align/$sample.bam` (same form with `$ref2`); the index build and alignment blocks are commented out, so `$dir/01.align/<sample>.bam` must already exist, and the `03.*braker` branch needs the `04_braker` output.
**Tools**: StringTie `-p 64 -G <gtf> -e -B -o <sample>.gtf -A <sample>.tsv`; HISAT2 `hisat2-build -p 64` and `hisat2 -p 64 -x <index> -1 <fq1> -2 <fq2> -S <sample>.sam --new-summary`; samtools `view -H`, `sort -@ 64`, `index` (the alignment blocks are commented out).

</details>

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>01_align/align.sh</code> · Shell · 76 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/02_annotation/02_rnaseq_align/01_align/align.sh)

```bash
#!/bin/bash
#CSUB -J hisat
#CSUB -q c01
#CSUB -o hisat.out
#CSUB -e hisat.error
#CSUB -n 64
#CSUB -R span[hosts=1]

#cd {PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/00.index
#hisat2-build -p 64 Gifu_v1.0.fasta Gifu
#index={PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/00.index/Gifu
dir={PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq
ref1={PROJ_LOTUS_ZC}/Gifu/00.assemble_result/Lotus_GifuT2T.gtf
ref2={PROJ_LOTUS_ZC}/Gifu/09.annotation/06.braker/Gifu_filtered.gtf
cat {PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/align.conf | while read id;
do
arr=($id)
sample=${arr[0]}
fq1=${arr[1]}
fq2=${arr[2]}

cd $dir
#hisat2 -p 64 -x $index -1 $fq1 -2 $fq2 -S $sample.sam --summary-file $sample_hisat.log --new-summary
#grep "NH:i:1" $sample.sam | grep "YT:Z:CP" > ${sample}_unique.sam
#samtools view -H $sample.sam > ${sample}_head.sam
#cat ${sample}_unique.sam >> ${sample}_head.sam
#samtools sort -@ 64 -o $sample.bam ${sample}_head.sam
#rm *sam
#samtools index $sample.bam
mkdir 02.stringtie_geta
mkdir 02.stringtie_geta/$sample && cd 02.stringtie_geta/$sample
stringtie -p 64 -G $ref1 -e -B -o $sample.gtf -A $sample.tsv $dir/01.align/$sample.bam

cd $dir
mkdir 03.stringtie_braker
mkdir 03.stringtie_braker/$sample && cd 03.stringtie_braker/$sample
stringtie -p 64 -G $ref2 -e -B -o $sample.gtf -A $sample.tsv $dir/01.align/$sample.bam

cd $dir
mkdir 02.tpm_geta
sed '1d' 02.stringtie_geta/$sample/$sample.tsv | awk '{print $1"\t"$9}'|sort > 02.tpm_geta/$sample.tpm

cd $dir
mkdir 03.tpm_braker
sed '1d' 03.stringtie_braker/$sample/$sample.tsv | awk '{print $1"\t"$9}'|sort > 03.tpm_braker/$sample.tpm

done

cd $dir/02.tpm_geta
ls *.tpm | sort > tpm.files
mapfile -t files < tpm.files
cp "${files[0]}" __merged.tmp
for f in "${files[@]:1}"; do
  join -t $'\t' -a1 -a2 -e 0 -o auto __merged.tmp "$f" > __merged.next
  mv __merged.next __merged.tmp
done
header="gene_id"
for f in "${files[@]}"; do header+="\t${f%.tpm}"; done
{ echo -e "$header"; cat __merged.tmp; } > all_samples.tpm.tsv

rm -f __merged.tmp tpm.files


cd $dir/03.tpm_braker
ls *.tpm | sort > tpm.files
mapfile -t files < tpm.files
cp "${files[0]}" __merged.tmp
for f in "${files[@]:1}"; do
  join -t $'\t' -a1 -a2 -e 0 -o auto __merged.tmp "$f" > __merged.next
  mv __merged.next __merged.tmp
done
header="gene_id"
for f in "${files[@]}"; do header+="\t${f%.tpm}"; done
{ echo -e "$header"; cat __merged.tmp; } > all_samples.tpm.tsv

rm -f __merged.tmp tpm.files
```

</details>
