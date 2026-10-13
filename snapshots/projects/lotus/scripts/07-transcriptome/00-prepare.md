---
title: '数据准备'
sidebar:
  label: '数据准备'
  order: 0
---


## 中文
**作用**：bulk RNA-seq 上游处理——fastp 质控、HISAT2 比对到 Gifu T2T 参考、StringTie 逐样本定量，并把各样本 TPM 列合并为一张矩阵。
**入口**：`fastp.sh`、`align.sh`
**输入**：`align.conf`（样本名、R1 FASTQ、R2 FASTQ 三列）及其指向的 `raw/<sample>-R1.fq.gz`、`-R2.fq.gz`；建索引用 `Gifu_v1.0.fasta`，StringTie 用 `Lotus_GifuT2T.gtf`。
**输出**：`clean/<sample>-R{1,2}.fastq.gz`、`report/<sample>.json|.html`、排序并建索引的 `<sample>.bam(.bai)`、`02.stringtie/<sample>/<sample>.gtf|.tsv`、`02.tpm/<sample>.tpm` 与合并矩阵 `02.tpm/all_samples.tpm.tsv`。
**运行**：`bash fastp.sh`、`bash align.sh`（`fastp.sh` 须在本目录内运行；`align.sh` 的 `$index` 从未赋值、`-x` 拿不到索引路径，且须以 `01.align` 为当前目录，运行前需修正脚本）。
**工具**：fastp（`-w 88`，JSON/HTML 报告）、HISAT2（`hisat2-build -p 64`；`-p 64 --new-summary`）、samtools（`view -H`、`sort -@ 64`、`index`）、StringTie（`-p 64 -G -e -B -o -A`）、Awk 与 `sort`/`join` 合并 TPM。

<details>
<summary>English</summary>

**Purpose**: Upstream bulk RNA-seq processing — fastp read QC, HISAT2 alignment to the Gifu T2T reference and StringTie per-sample quantification, joined into one TPM matrix.
**Entry point**: `fastp.sh`, `align.sh`
**Inputs**: `align.conf` (sample name, R1 FASTQ, R2 FASTQ) and the `raw/<sample>-R1.fq.gz` / `-R2.fq.gz` files it points to; `Gifu_v1.0.fasta` for the index and `Lotus_GifuT2T.gtf` for StringTie.
**Outputs**: `clean/<sample>-R{1,2}.fastq.gz`, `report/<sample>.json|.html`, the sorted and indexed `<sample>.bam(.bai)`, `02.stringtie/<sample>/<sample>.gtf|.tsv`, `02.tpm/<sample>.tpm` and the joined matrix `02.tpm/all_samples.tpm.tsv`.
**Run**: `bash fastp.sh`, `bash align.sh` (`fastp.sh` must run in this directory; in `align.sh` `$index` is never set, so `-x` receives no index path, and `01.align` must be the current directory — fix the script before running).
**Tools**: fastp (`-w 88`, JSON/HTML reports), HISAT2 (`hisat2-build -p 64`; `-p 64 --new-summary`), samtools (`view -H`, `sort -@ 64`, `index`), StringTie (`-p 64 -G -e -B -o -A`), Awk with `sort`/`join` for the TPM merge.

</details>

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>align.conf</code> · ini · 39 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/07_transcriptome/00_prepare/align.conf)

```ini
Gifu_flower_rep1	raw/Gifu_flower_rep1-R1.fq.gz	raw/Gifu_flower_rep1-R2.fq.gz
Gifu_flower_rep2	raw/Gifu_flower_rep2-R1.fq.gz	raw/Gifu_flower_rep2-R2.fq.gz
Gifu_flower_rep3	raw/Gifu_flower_rep3-R1.fq.gz	raw/Gifu_flower_rep3-R2.fq.gz
Gifu_leaf_rep1	raw/Gifu_leaf_rep1-R1.fq.gz	raw/Gifu_leaf_rep1-R2.fq.gz
Gifu_leaf_rep2	raw/Gifu_leaf_rep2-R1.fq.gz	raw/Gifu_leaf_rep2-R2.fq.gz
Gifu_leaf_rep3	raw/Gifu_leaf_rep3-R1.fq.gz	raw/Gifu_leaf_rep3-R2.fq.gz
Gifu_nodule_dpi10_rep1	raw/Gifu_nodule_dpi10_rep1-R1.fq.gz	raw/Gifu_nodule_dpi10_rep1-R2.fq.gz
Gifu_nodule_dpi10_rep2	raw/Gifu_nodule_dpi10_rep2-R1.fq.gz	raw/Gifu_nodule_dpi10_rep2-R2.fq.gz
Gifu_nodule_dpi10_rep3	raw/Gifu_nodule_dpi10_rep3-R1.fq.gz	raw/Gifu_nodule_dpi10_rep3-R2.fq.gz
Gifu_nodule_dpi21_rep1	raw/Gifu_nodule_dpi21_rep1-R1.fq.gz	raw/Gifu_nodule_dpi21_rep1-R2.fq.gz
Gifu_nodule_dpi21_rep2	raw/Gifu_nodule_dpi21_rep2-R1.fq.gz	raw/Gifu_nodule_dpi21_rep2-R2.fq.gz
Gifu_nodule_dpi21_rep3	raw/Gifu_nodule_dpi21_rep3-R1.fq.gz	raw/Gifu_nodule_dpi21_rep3-R2.fq.gz
Gifu_nodule_dpi45_rep1	raw/Gifu_nodule_dpi45_rep1-R1.fq.gz	raw/Gifu_nodule_dpi45_rep1-R2.fq.gz
Gifu_nodule_dpi45_rep2	raw/Gifu_nodule_dpi45_rep2-R1.fq.gz	raw/Gifu_nodule_dpi45_rep2-R2.fq.gz
Gifu_nodule_dpi45_rep3	raw/Gifu_nodule_dpi45_rep3-R1.fq.gz	raw/Gifu_nodule_dpi45_rep3-R2.fq.gz
Gifu_pod_rep1	raw/Gifu_pod_rep1-R1.fq.gz	raw/Gifu_pod_rep1-R2.fq.gz
Gifu_pod_rep2	raw/Gifu_pod_rep2-R1.fq.gz	raw/Gifu_pod_rep2-R2.fq.gz
Gifu_pod_rep3	raw/Gifu_pod_rep3-R1.fq.gz	raw/Gifu_pod_rep3-R2.fq.gz
Gifu_root_dpi02_rep1	raw/Gifu_root_dpi02_rep1-R1.fq.gz	raw/Gifu_root_dpi02_rep1-R2.fq.gz
Gifu_root_dpi02_rep2	raw/Gifu_root_dpi02_rep2-R1.fq.gz	raw/Gifu_root_dpi02_rep2-R2.fq.gz
Gifu_root_dpi02_rep3	raw/Gifu_root_dpi02_rep3-R1.fq.gz	raw/Gifu_root_dpi02_rep3-R2.fq.gz
Gifu_root_dpi10_rep1	raw/Gifu_root_dpi10_rep1-R1.fq.gz	raw/Gifu_root_dpi10_rep1-R2.fq.gz
Gifu_root_dpi10_rep2	raw/Gifu_root_dpi10_rep2-R1.fq.gz	raw/Gifu_root_dpi10_rep2-R2.fq.gz
Gifu_root_dpi10_rep3	raw/Gifu_root_dpi10_rep3-R1.fq.gz	raw/Gifu_root_dpi10_rep3-R2.fq.gz
Gifu_root_dpi21_rep1	raw/Gifu_root_dpi21_rep1-R1.fq.gz	raw/Gifu_root_dpi21_rep1-R2.fq.gz
Gifu_root_dpi21_rep2	raw/Gifu_root_dpi21_rep2-R1.fq.gz	raw/Gifu_root_dpi21_rep2-R2.fq.gz
Gifu_root_dpi21_rep3	raw/Gifu_root_dpi21_rep3-R1.fq.gz	raw/Gifu_root_dpi21_rep3-R2.fq.gz
Gifu_root_dpi45_rep1	raw/Gifu_root_dpi45_rep1-R1.fq.gz	raw/Gifu_root_dpi45_rep1-R2.fq.gz
Gifu_root_dpi45_rep2	raw/Gifu_root_dpi45_rep2-R1.fq.gz	raw/Gifu_root_dpi45_rep2-R2.fq.gz
Gifu_root_dpi45_rep3	raw/Gifu_root_dpi45_rep3-R1.fq.gz	raw/Gifu_root_dpi45_rep3-R2.fq.gz
Gifu_root_hpi48_rep1	raw/Gifu_root_hpi48_rep1-R1.fq.gz	raw/Gifu_root_hpi48_rep1-R2.fq.gz
Gifu_root_hpi48_rep2	raw/Gifu_root_hpi48_rep2-R1.fq.gz	raw/Gifu_root_hpi48_rep2-R2.fq.gz
Gifu_root_hpi48_rep3	raw/Gifu_root_hpi48_rep3-R1.fq.gz	raw/Gifu_root_hpi48_rep3-R2.fq.gz
Gifu_root_uni_rep1	raw/Gifu_root_uni_rep1-R1.fq.gz	raw/Gifu_root_uni_rep1-R2.fq.gz
Gifu_root_uni_rep2	raw/Gifu_root_uni_rep2-R1.fq.gz	raw/Gifu_root_uni_rep2-R2.fq.gz
Gifu_root_uni_rep3	raw/Gifu_root_uni_rep3-R1.fq.gz	raw/Gifu_root_uni_rep3-R2.fq.gz
Gifu_stem_rep1	raw/Gifu_stem_rep1-R1.fq.gz	raw/Gifu_stem_rep1-R2.fq.gz
Gifu_stem_rep2	raw/Gifu_stem_rep2-R1.fq.gz	raw/Gifu_stem_rep2-R2.fq.gz
Gifu_stem_rep3	raw/Gifu_stem_rep3-R1.fq.gz	raw/Gifu_stem_rep3-R2.fq.gz
```

</details>

<details open>
<summary><code>align.sh</code> · Shell · 51 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/07_transcriptome/00_prepare/align.sh)

```bash
#!/bin/bash
#CSUB -J stringtie
#CSUB -q c01
#CSUB -o stringtie.out
#CSUB -e stringtie.error
#CSUB -n 64
#CSUB -R span[hosts=1]

cd {PROJ_LOTUS_HX}/03_annotation/Gifu/04.rnasq/00.index
hisat2-build -p 64 Gifu_v1.0.fasta Gifu
#ndex={PROJ_LOTUS_HX}/03_annotation/Gifu/04.rnasq/00.index/Gifu
dir={PROJ_LOTUS_HX}/03_annotation/Gifu/04.rnasq
ref1={PROJ_LOTUS_HX}/00_assembly_result/GifuT2T/Lotus_GifuT2T.gtf
cat {PROJ_LOTUS_HX}/01_data/Gifu/RNA/clean/align.conf | while read id;
do
arr=($id)
sample=${arr[0]}
fq1=${arr[1]}
fq2=${arr[2]}

cd $dir
hisat2 -p 64 -x $index -1 $fq1 -2 $fq2 -S $sample.sam --summary-file $sample_hisat.log --new-summary
grep "NH:i:1" $sample.sam | grep "YT:Z:CP" > ${sample}_unique.sam
samtools view -H $sample.sam > ${sample}_head.sam
cat ${sample}_unique.sam >> ${sample}_head.sam
samtools sort -@ 64 -o $sample.bam ${sample}_head.sam
rm *sam
samtools index $sample.bam
mkdir 02.stringtie
mkdir 02.stringtie/$sample && cd 02.stringtie/$sample
stringtie -p 64 -G $ref1 -e -B -o $sample.gtf -A $sample.tsv $dir/01.align/$sample.bam

cd $dir
mkdir 02.tpm
sed '1d' 02.stringtie/$sample/$sample.tsv | awk '{print $1"\t"$9}'|sort > 02.tpm/$sample.tpm

done

cd $dir/02.tpm
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

<details open>
<summary><code>fastp.sh</code> · Shell · 18 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/07_transcriptome/00_prepare/fastp.sh)

```bash
#!/bin/bash
#CSUB -J fastp
#CSUB -q c01
#CSUB -o fastp.out
#CSUB -e fastp.error
#CSUB -n 88
#CSUB -R span[hosts=1]

cat align.conf|while read id;
do
arr=($id)
sample=${arr[0]}
fq1=${arr[1]}
fq2=${arr[2]}

fastp -i $fq1 -I $fq2 -o clean/${sample}-R1.fastq.gz -O clean/${sample}-R2.fastq.gz -w 88 -j report/${sample}.json -h report/${sample}.html

done
```

</details>
