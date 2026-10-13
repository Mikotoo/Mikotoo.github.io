---
title: '08 单核转录组'
sidebar:
  label: '08 单核转录组'
  order: -1
---


## 中文
**作用**：把 Gifu T2T 基因组上的单核 RNA-seq（根的两个阶段 Dpi00 与 Dpi02，各两个重复）从 Cell Ranger 定量做到聚类与细胞类型注释、簇/细胞类型 DEG 与 scTenifoldKnk 虚拟敲除。
**入口**：下表 7 个 shell 脚本；`Script/*.py` 中需手动运行的为 `03_1-cluster_statitcs.py`、`03_2-fliter_to_all_gene_h5ad.py`、`05_3-marker_cluster_DEGs_own_order.py` 与 `makesure_02dpi_samples_all_have_Early_signal_cell.py`。
**输入**：`sample.list`（`ProjectID CNName RunID SampleName`）、四个 `*_filtered.h5ad`、Gifu T2T FASTA/GTF 与 `.pep.fa`、`{Gmax,Gm2.1,GmZH,At}_marker_genes.pep.fa`、eggNOG 注释；后续步骤按当前目录读 `02_processData/` 与 `04_Anno/`，scTenifoldKnk 另需 `zz_00_pea.rds` 与 `All_species_merged_by_GifuT2T.SNF_Symbol.tsv`。
**输出**：`00_ref/`、`01_cellranger_out/<SampleName>/`、`02_processData/`（h5ad、`zz_*`、`zz_00_pea.rds`、`00_02_pseudobulk_fpkm.csv`）、`03_figures/`（降维图、`clustree_plot_*.pdf`、火山图）、`04_Anno/`（Seurat rds、marker CSV、`blast_*/`、富集结果）、`06_batch_multiple_scTenifoldKnk/<group>/<gene>/`。
**运行**：在本阶段根目录（含 `02_processData/`、`03_figures/`、`04_Anno/`）按表中顺序运行，编号脚本携带 `#CSUB` 作业指令、交由调度器提交；``{TOKEN}`` 路径占位符见 `scripts/README.md`；注意 `03_1-dimenRedu_adjusted.py` 写出的文件名与 `03_1-cluster_statitcs.py`、`03_2-fliter_to_all_gene_h5ad.py` 读取的文件名不一致，需先统一。
**工具**：Cell Ranger 9.0.1；Python 3.9.19 + scanpy/anndata/scvi-tools/leidenalg/rpy2；R 4.4.1 + Seurat/scTenifoldKnk/DESeq2/clustree；BLAST+ 与 seqkit；参数：`n_top_genes` 4000（扫描 1000–5000）、Leiden 0.4（扫描 0.2–1.2）、每基因 100 次 scTenifoldKnk 重复（K ∈ {10,30,50,70,100}）。
| 脚本 | 用途 |
|---|---|
| `00_mkref.sh` | `cellranger mkref --genome=Gifu`（Gifu T2T FASTA + GTF）→ `./00_ref` |
| `00_cellranger.sh` | 对 `sample.list` 中每个样本执行一次 `cellranger count` → `01_cellranger_out/<SampleName>` |
| `run_singlecell.sh` | 一个驱动依次跑 `Script/01_annaDataQC.py`、`02_integrate.py`、`03_1-dimenRedu_adjusted.py`、`05_3-markerAnno_Ps_v5_own.py`、`06_4-DEGs_cell_type.R` |
| `05_cluster_DEGs_finder.sh` | 先 `05_1-…py` 再 `05_2-…R`：Seurat 对象与簇 marker，随后按簇做富集 |
| `06_1.marker_gene_finder.sh` | DEG → 肽 FASTA → 对 Gmax/Gm2.1/GmZH/At marker 做 `blastp` 并合并结果 |
| `06_1.Individual_clusters_marker_gene_finder.sh` | 取各簇 top-20 DEG 序列并抽 eggNOG 注释（启用簇 `8`） |
| `07_batch_multiple_scTenifoldKnk.sh` | 每个候选基因生成一个 scTenifoldKnk 作业脚本与 `scTenifoldKnk_run.list` |

<details>
<summary>English</summary>

**Purpose**: Take the single-nucleus RNA-seq data on the Gifu T2T genome (the two root stages Dpi00 and Dpi02, two replicates each) from Cell Ranger quantification through clustering and cell-type annotation, cluster/cell-type DEGs and scTenifoldKnk virtual knockout.
**Entry point**: the seven shell scripts in the table below; the `Script/*.py` files with no entry script, run by hand, are `03_1-cluster_statitcs.py`, `03_2-fliter_to_all_gene_h5ad.py`, `05_3-marker_cluster_DEGs_own_order.py` and `makesure_02dpi_samples_all_have_Early_signal_cell.py`.
**Inputs**: `sample.list` (`ProjectID CNName RunID SampleName`), the four `*_filtered.h5ad` objects, Gifu T2T FASTA/GTF and `.pep.fa`, `{Gmax,Gm2.1,GmZH,At}_marker_genes.pep.fa`, eggNOG annotation; later steps read `02_processData/` and `04_Anno/` from the working directory, and scTenifoldKnk additionally needs `zz_00_pea.rds` and `All_species_merged_by_GifuT2T.SNF_Symbol.tsv`.
**Outputs**: `00_ref/`, `01_cellranger_out/<SampleName>/`, `02_processData/` (h5ad files, `zz_*`, `zz_00_pea.rds`, `00_02_pseudobulk_fpkm.csv`), `03_figures/` (dimension-reduction plots, `clustree_plot_*.pdf`, volcano plots), `04_Anno/` (Seurat rds, marker CSVs, `blast_*/`, enrichment results) and `06_batch_multiple_scTenifoldKnk/<group>/<gene>/`.
**Run**: from the stage root (holding `02_processData/`, `03_figures/`, `04_Anno/`) run the scripts in table order; the numbered scripts carry `#CSUB` directives and are submitted to the scheduler; for the ``{TOKEN}`` path placeholders see `scripts/README.md`; note that `03_1-dimenRedu_adjusted.py` writes a file name that `03_1-cluster_statitcs.py` and `03_2-fliter_to_all_gene_h5ad.py` do not read, so it must be unified first.
**Tools**: Cell Ranger 9.0.1; Python 3.9.19 with scanpy/anndata/scvi-tools/leidenalg/rpy2; R 4.4.1 with Seurat/scTenifoldKnk/DESeq2/clustree; BLAST+ and seqkit; settings as coded: `n_top_genes` 4000 (scanned 1000–5000), Leiden 0.4 (scanned 0.2–1.2), and 100 scTenifoldKnk replicates per gene (K ∈ {10,30,50,70,100}).
| Script | Purpose |
|---|---|
| `00_mkref.sh` | `cellranger mkref --genome=Gifu` from the Gifu T2T FASTA + GTF → `./00_ref` |
| `00_cellranger.sh` | one `cellranger count` per sample in `sample.list` → `01_cellranger_out/<SampleName>` |
| `run_singlecell.sh` | one driver running `Script/01_annaDataQC.py`, `02_integrate.py`, `03_1-dimenRedu_adjusted.py`, `05_3-markerAnno_Ps_v5_own.py`, `06_4-DEGs_cell_type.R` in order |
| `05_cluster_DEGs_finder.sh` | `05_1-…py` then `05_2-…R`: Seurat object and cluster markers, then per-cluster enrichment |
| `06_1.marker_gene_finder.sh` | DEGs → peptide FASTA → `blastp` against the Gmax/Gm2.1/GmZH/At markers, results merged |
| `06_1.Individual_clusters_marker_gene_finder.sh` | top-20 DEG sequences per cluster plus eggNOG annotation (active cluster `8`) |
| `07_batch_multiple_scTenifoldKnk.sh` | one scTenifoldKnk job script per candidate gene plus `scTenifoldKnk_run.list` |

</details>

### 本阶段步骤

| 步骤 | 内容 |
| --- | --- |
| [单核分析脚本](/_projects/lotus/scripts/08-single-nucleus/script/) | — |

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>00_cellranger.sh</code> · Shell · 41 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/08_single_nucleus/00_cellranger.sh)

```bash
#!/bin/bash
#CSUB -J cellranger
#CSUB -q c01
#CSUB -o cellranger.out
#CSUB -e cellranger.error
#CSUB -n 64
#CSUB -R span[hosts=1]

REF={PROJ_LOTUS_HX}/07_singleCell/00_ref/fasta/genome.fa
DATADIR={HOME_HX}/Single_cell_data
SAMPLELIST=sample.list

# module load cellranger/7.2.0   # 如果集群环境要求，请加载软件

while read ProjectID CNName RunID SampleName
do
    echo ">>> Processing $SampleName ($RunID)"

    R1=${DATADIR}/${ProjectID}/${RunID}/${RunID}_S1_L001_R1_001.fastq.gz
    R2=${DATADIR}/${ProjectID}/${RunID}/${RunID}_S1_L001_R2_001.fastq.gz

    if [[ ! -f $R1 ]]; then
        echo "❌ Missing R1 file: $R1"
        continue
    fi

    outdir=./01_cellranger_out/${SampleName}
    mkdir -p ${outdir}

    cellranger count \
        --id=${SampleName} \
        --transcriptome={PROJ_LOTUS_HX}/07_singleCell/00_ref \
        --fastqs=$(dirname $R1) \
        --sample=${RunID} \
        --create-bam=true \
        --localcores=32 \
        --localmem=120 \
        --chemistry=SC3Pv3 \
        --output-dir=${outdir}

done < ${SAMPLELIST}
```

</details>

<details open>
<summary><code>00_mkref.sh</code> · Shell · 17 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/08_single_nucleus/00_mkref.sh)

```bash
#!/bin/bash
#CSUB -J mkref
#CSUB -q c01
#CSUB -o mkref.out
#CSUB -e mkref.error
#CSUB -n 64
#CSUB -R span[hosts=1]

reffasta={PROJ_LOTUS_HX}/00_assembly_result/GifuT2T/Lotus_GifuT2T_v1.0.fasta
refgtf={PROJ_LOTUS_HX}/00_assembly_result/GifuT2T/Lotus_GifuT2T.gtf

cellranger mkref \
  --genome=Gifu \
  --fasta=$reffasta \
  --genes=$refgtf \
  --output-dir=./00_ref \
  --nthreads=64
```

</details>

<details open>
<summary><code>05_cluster_DEGs_finder.sh</code> · Shell · 76 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/08_single_nucleus/05_cluster_DEGs_finder.sh)

```bash
#!/bin/bash
#CSUB -J clusterDEG
#CSUB -q c01
#CSUB -o clusterDEG.out
#CSUB -e clusterDEG.error
#CSUB -n 64
#CSUB -R span[hosts=1]

# 加载 conda
source ~/anaconda3/etc/profile.d/conda.sh

# 激活环境
conda activate singlecell

# 运行脚本
python ./Script/05_1-Seurat_Exp_h5ad_to_rds.py
Rscript ./Script/05_2-Exp_h5ad_to_rds_cluster.R

# 对每簇的marker genes的功能进行富集
#========================
# 0. 基本路径设置
#========================
output=./04_Anno/cluster_makers_genes_enrichment
mkdir -p "$output"

script={SOFTWARE_HX}/Script/RNA-seq
anno={PROJ_LOTUS_HX}/00_assembly_result/GifuT2T/Lotus_Gifu.emapper.annotations
file=./04_Anno/cluster_06_cluster_markers_gene1.csv


#========================
# 1. 获取所有 cluster 编号
#========================
clusters=$(tail -n +2 "$file" \
  | awk -F',' '{gsub(/"/,""); print $6}' \
  | sort -u)

#========================
# 2. 逐 cluster 筛选基因 + 富集
#========================
for cl in $clusters; do
  echo "处理 cluster ${cl}..."

  # 固定文件名，不使用 mktemp，避免报错
  tmpfile="${output}/cluster_${cl}.genelist"

  # 筛选基因并替换 "-" 为 "_"
  awk -F',' -v c="$cl" '
    BEGIN { OFS="\t" }
    NR>1 {
      gsub(/"/, "", $0)
      if ($6 == c && $2 > 1 && $5 < 0.05) {
        gsub("-", "_", $7)
        print $7
      }
    }
  ' "$file" > "$tmpfile"


  # 若无基因，则跳过
  if [ ! -s "$tmpfile" ]; then
    echo "  cluster ${cl} 没有满足条件的基因，跳过富集分析"
    rm "$tmpfile"
    continue
  fi

  # cluster 单独输出目录
  outdir="${output}/cluster_${cl}"
  mkdir -p "$outdir"

  echo "  正在进行富集分析..."
  Rscript "${script}/Enrichment.R" "$anno" "$tmpfile" "$outdir"

done

echo "全部 cluster 处理完毕！"
```

</details>

<details open>
<summary><code>06_1.Individual_clusters_marker_gene_finder.sh</code> · Shell · 55 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/08_single_nucleus/06_1.Individual_clusters_marker_gene_finder.sh)

```bash
#!/bin/bash
#CSUB -J findmarker
#CSUB -q c01
#CSUB -o findmarker.out
#CSUB -e findmarker.error
#CSUB -n 64
#CSUB -R span[hosts=1]

# 输入文件
csv=./04_Anno/cluster_06_cluster_markers_gene.csv
fasta={PROJ_LOTUS_HX}/00_assembly_result/GifuT2T/Lotus_GifuT2T.pep.fa
annot=./Lotus_Gifu.emapper.annotations
outdir=./04_Anno

mkdir -p $outdir

# 需要处理的 cluster
#clusters=("3" "4" "6")
clusters=("8")

for n in "${clusters[@]}"; do
    echo ">>> Processing cluster $n ..."

    # Step1: 获取 top10 基因列表
    awk -F',' -v c="$n" '
        NR>1 {
            gsub(/"/, "", $0);                  
            if ($6 == c) {
                gsub(/-/, "_", $7);             
                print $2 "\t" $7                
            }
        }' "$csv" \
        | sort -k1,1nr \
        | head -n 20 \
        | cut -f2 \
        > ${outdir}/cluster_${n}_DEGs.list

    # Step2: 从 pep.fa 提取序列
    seqkit grep -f ${outdir}/cluster_${n}_DEGs.list "$fasta" \
        > ${outdir}/cluster_${n}_DEGs.pep.fa

    echo "Sequence extraction for cluster $n done."

    # Step3: 从 EggNOG 注释中提取对应基因的注释
    # 注释文件的第一列是 query ID
    grep -v "^#" "$annot" | awk -v OFS="\t" '
        NR==FNR {a[$1]=1; next}
        ($1 in a) {print}
    ' ${outdir}/cluster_${n}_DEGs.list - \
        > ${outdir}/cluster_${n}_DEGs.annotation

    echo "Annotation extraction for cluster $n done."
done

echo "All clusters finished."
```

</details>

<details open>
<summary><code>06_1.marker_gene_finder.sh</code> · Shell · 175 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/08_single_nucleus/06_1.marker_gene_finder.sh)

```bash
#!/bin/bash
#CSUB -J findmarker
#CSUB -q c01
#CSUB -o findmarker.out
#CSUB -e findmarker.error
#CSUB -n 64
#CSUB -R span[hosts=1]

cul=`pwd`

# 输入文件
csv=./04_Anno/cluster_06_cluster_markers_gene.csv
fasta={PROJ_LOTUS_HX}/00_assembly_result/GifuT2T/Lotus_GifuT2T.pep.fa
out=./04_Anno/cluster_DEGs.pep.fa

# 对DEGs进行过滤，提取 gene 列，去掉表头并去重，并将-替换为_
awk -F',' 'NR>1 {
    gsub(/"/, "", $0);
    if ( $5 < 0.05) {
        gsub(/-/, "_", $7);
        print $7
    }
}' "$csv" | sort -u > ./04_Anno/cluster_DEGs.list

# 用 seqkit 从 pep.fa 中提取对应序列
seqkit grep -f ./04_Anno/cluster_DEGs.list "$fasta" > "$out"

echo "✅ Done! Extracted sequences saved to: $out"

outdir=./04_Anno/blast_Gmax
mkdir -p "$outdir"

# 构建blast数据库
makeblastdb -in "./04_Anno/Gmax_marker_genes.pep.fa" -dbtype prot -out "$outdir/db"

# 运行blastp
blastp -query "$out" -db "$outdir/db" \
    -out "$outdir/result.blast" \
    -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
    -num_threads 32

cut -f1,2,4,3 $outdir/result.blast > $outdir/result.blast.res

seqkit fx2tab -l -n -i ./04_Anno/Gmax_marker_genes.pep.fa > $outdir/ref_protein_lengths.txt

seqkit fx2tab -l -n -i $out > $outdir/query_protein_lengths.txt

cd $outdir

python3 $cul/blast_merge.py result.blast.res query_protein_lengths.txt ref_protein_lengths.txt blast1.merge

awk '$3>50 && $7>0.5 && $8>0.5' blast1.merge | \
awk '{
    diff = (1-$7>=0?1-$7:$7-1) + (1-$8>=0?1-$8:$8-1)
    key = $2
    if (!(key in best) || $3 > best_id[key] || ($3 == best_id[key] && diff < best_diff[key])) {
        best[key] = $0
        best_id[key] = $3
        best_diff[key] = diff
    }
}
END {
    for (k in best) print best[k]
}' > blast1.filtered.best

cd $cul
outdir=./04_Anno/blast_Gm2.1
mkdir -p "$outdir"

# 构建blast数据库
makeblastdb -in "./04_Anno/Gm2.1_marker_genes.pep.fa" -dbtype prot -out "$outdir/db"

# 运行blastp
blastp -query "$out" -db "$outdir/db" \
    -out "$outdir/result.blast" \
    -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
    -num_threads 32

cut -f1,2,4,3 $outdir/result.blast > $outdir/result.blast.res

seqkit fx2tab -l -n -i ./04_Anno/Gm2.1_marker_genes.pep.fa > $outdir/ref_protein_lengths.txt

seqkit fx2tab -l -n -i $out > $outdir/query_protein_lengths.txt

cd $outdir

python3 $cul/blast_merge.py result.blast.res query_protein_lengths.txt ref_protein_lengths.txt blast1.merge

awk '$3>50 && $7>0.5 && $8>0.5' blast1.merge | \
awk '{
    diff = (1-$7>=0?1-$7:$7-1) + (1-$8>=0?1-$8:$8-1)
    key = $2
    if (!(key in best) || $3 > best_id[key] || ($3 == best_id[key] && diff < best_diff[key])) {
        best[key] = $0
        best_id[key] = $3
        best_diff[key] = diff
    }
}
END {
    for (k in best) print best[k]
}' > blast1.filtered.best

cd $cul
outdir=./04_Anno/blast_GmZH
mkdir -p "$outdir"

# 构建blast数据库
makeblastdb -in "./04_Anno/GmZH_marker_genes.pep.fa" -dbtype prot -out "$outdir/db"

# 运行blastp
blastp -query "$out" -db "$outdir/db" \
    -out "$outdir/result.blast" \
    -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
    -num_threads 32

cut -f1,2,4,3 $outdir/result.blast > $outdir/result.blast.res

seqkit fx2tab -l -n -i ./04_Anno/GmZH_marker_genes.pep.fa > $outdir/ref_protein_lengths.txt

seqkit fx2tab -l -n -i $out > $outdir/query_protein_lengths.txt

cd $outdir

python3 $cul/blast_merge.py result.blast.res query_protein_lengths.txt ref_protein_lengths.txt blast1.merge

awk '$3>50 && $7>0.5 && $8>0.5' blast1.merge | \
awk '{
    diff = (1-$7>=0?1-$7:$7-1) + (1-$8>=0?1-$8:$8-1)
    key = $2
    if (!(key in best) || $3 > best_id[key] || ($3 == best_id[key] && diff < best_diff[key])) {
        best[key] = $0
        best_id[key] = $3
        best_diff[key] = diff
    }
}
END {
    for (k in best) print best[k]
}' > blast1.filtered.best

cd $cul
outdir=./04_Anno/blast_At
mkdir -p "$outdir"

# 构建blast数据库
makeblastdb -in "./04_Anno/At_marker_genes.pep.fa" -dbtype prot -out "$outdir/db"

# 运行blastp
blastp -query "$out" -db "$outdir/db" \
    -out "$outdir/result.blast" \
    -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
    -num_threads 32

cut -f1,2,4,3 $outdir/result.blast > $outdir/result.blast.res

seqkit fx2tab -l -n -i ./04_Anno/At_marker_genes.pep.fa > $outdir/ref_protein_lengths.txt

seqkit fx2tab -l -n -i $out > $outdir/query_protein_lengths.txt

cd $outdir

python3 $cul/blast_merge.py result.blast.res query_protein_lengths.txt ref_protein_lengths.txt blast1.merge

awk '$3>50 && $7>0.5 && $8>0.5' blast1.merge | \
awk '{
    diff = (1-$7>=0?1-$7:$7-1) + (1-$8>=0?1-$8:$8-1)
    key = $2
    if (!(key in best) || $3 > best_id[key] || ($3 == best_id[key] && diff < best_diff[key])) {
        best[key] = $0
        best_id[key] = $3
        best_diff[key] = diff
    }
}
END {
    for (k in best) print best[k]
}' > blast1.filtered.best
```

</details>

<details open>
<summary><code>07_batch_multiple_scTenifoldKnk.sh</code> · Shell · 59 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/08_single_nucleus/07_batch_multiple_scTenifoldKnk.sh)

```bash
#!/bin/bash
set -euo pipefail

cul="$(pwd)"
runlist="$cul/scTenifoldKnk_run.list"
genelist=$1
celltype=$2

# 清空 run list
: > "$runlist"

while IFS=$'\t ' read -r gene group || [ -n "$gene" ]; do
  # 去掉 Windows CR
  gene="${gene//$'\r'/}"
  group="${group//$'\r'/}"

  # 跳过非法行
  [ -z "$gene" ] && continue
  [ -z "$group" ] && continue

  # 构建目录：06_batch_multiple_scTenifoldKnk / group / gene
  base_dir="$cul/06_batch_multiple_scTenifoldKnk/$group/$gene"
  mkdir -p "$base_dir"

  job="$base_dir/${gene}_scTenifoldKnk.sh"

  cat > "$job" <<EOF
#!/bin/bash
#CSUB -J ${gene}_Knk
#CSUB -q c02
#CSUB -o $base_dir/${gene}_Knk.out
#CSUB -e $base_dir/${gene}_Knk.error
#CSUB -n 16
#CSUB -R span[hosts=1]

set -eo pipefail   # conda 前不要 -u

source ~/anaconda3/etc/profile.d/conda.sh
conda activate singlecell

set -euo pipefail  # 激活后再启用 -u

# 进入 gene 专属目录
cd "$base_dir"

# 软链接数据（避免重复拷贝）
ln -sf $cul/02_processData/zz_00_pea.rds
ln -sf $cul/All_species_merged_by_GifuT2T.SNF_Symbol.tsv

# 执行 scTenifoldKnk
Rscript $cul/Script/07_batch_multiple_scTenifoldKnk.R "$gene" "$celltype"
EOF

  chmod +x "$job"
  echo "$job" >> "$runlist"

done < "$genelist"

echo "Generated run list: $runlist"
```

</details>

<details open>
<summary><code>run_singlecell.sh</code> · Shell · 43 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/08_single_nucleus/run_singlecell.sh)

```bash
#!/bin/bash
#CSUB -J singlecell
#CSUB -q c01
#CSUB -o singlecell.out
#CSUB -e singlecell.error
#CSUB -n 64
#CSUB -R span[hosts=1]

# Steps 01-03 and 06_2-06_3 of the single-nucleus pipeline, which were
# previously one six-line wrapper each. The longer steps keep their own
# scripts: 00_cellranger.sh, 00_mkref.sh, 05_cluster_DEGs_finder.sh,
# 06_1.marker_gene_finder.sh, 06_1.Individual_clusters_marker_gene_finder.sh
# and 07_batch_multiple_scTenifoldKnk.sh.
#
# Run from 08_single_nucleus/:
#   bash run_singlecell.sh

# 加载 conda
source ~/anaconda3/etc/profile.d/conda.sh

# 激活环境
conda activate singlecell

# ── 01 QC ────────────────────────────────────────────────────────────────
python ./Script/01_annaDataQC.py

# ── 02 integration ───────────────────────────────────────────────────────
python ./Script/02_integrate.py

# ── 03 dimensionality reduction and clustering ───────────────────────────
python ./Script/03_1-dimenRedu_adjusted.py

# ── 06_2 marker annotation ───────────────────────────────────────────────
# 已有的文献的marker gene效果不好，故改用先发掘各簇DEG，在比对到大豆marker gene的思路。
# python ./Script/05_3-marker_knowngene_lotus_own_order.py
python ./Script/05_3-markerAnno_Ps_v5_own.py

# ── 06_3 cell-type expression and DEGs ───────────────────────────────────
# The three calls below feed this one and are commented out in the original run.
# python ./Script/06_1-Seurat_Exp_h5ad_to_rds.py
# Rscript ./Script/06_2-Exp_h5ad_to_rds_cell_type.R
# Rscript ./Script/06_3-agg_to_fpkm.R
Rscript ./Script/06_4-DEGs_cell_type.R
```

</details>
