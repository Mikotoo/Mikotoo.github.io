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


