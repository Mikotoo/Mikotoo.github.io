---
title: '基因对应关系'
sidebar:
  label: '基因对应关系'
  order: 10
---


## 中文

**作用**：`blast.sh` 把 Gifu T2T 蛋白组分别对上一版 Gifu、上一版 MG20 与 MG20 T2T 蛋白组做 BLASTP，用 `blast_merge.py` 合并每次结果，并归约为每个 Gifu T2T 基因一条对应记录；三张配对表 `paste` 成 `All_species_merged_by_GifuT2T.tsv`，即四基因组之间的基因对应关系。
**入口**：`blast.sh`；`blast_merge.py` 可独立运行。
**输入**：工作目录中的 `Lotus_GifuT2T.pep.fa`（三次运行的查询，也是 MG20T2T 运行的主体）、`Lotus_GifuOld.pep.fa`、`Lotus_MG20Old.pep.fa`、`Lotus_MG20T2T.pep.fa`；`blast_merge.py`，以 `$cul/blast_merge.py` 调用（`$cul` 即工作目录）；写死路径的 `ALL_GENE_LIST`（只读第 1 列，为没有对应基因的条目输出 `NA` 行）。
**输出**：`GifuT2T_GifuOld/`、`GifuT2T_MG20Old/`、`GifuT2T_MG20T2T/`，各含 BLAST 库 `db.*`、`result.blast`、`result.blast.res`、`ref_len.txt`、`query_len.txt`、`blast1.merge`（脚本中使用的字面文件名）、`final_pairs.tsv` 等中间表；`All_species_merged_by_GifuT2T.tsv`（表头 `GifuT2T GifuOld MG20Old MG20T2T`，每个基因一行）。
**运行**：`bash blast.sh`；`ALL_GENE_LIST` 在脚本里是仓库外的绝对路径，须先改成实际文件；独立用法为 `python3 blast_merge.py <blast_result> <query_length_file> <ref_length_file> <output_file> [chunksize]`，chunksize 默认 2,000,000。
**工具**：BLAST+ `makeblastdb -dbtype prot` 与 `blastp`（`-evalue 1e-5`、`-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs"`、`-num_threads 88`）；Python 3 + pandas（`blast_merge.py`，分块 `chunksize=2_000_000`）；seqkit `fx2tab -l -n -i`；GNU Awk 过滤 `$3 > 50`（确信对为 `$3==100 && $7==1 && $8==1`）；`cut`、`sort`、`paste`、`find`。

<details>
<summary>English</summary>
**Purpose**: `blast.sh` runs BLASTP of the Gifu T2T proteome against the previous Gifu, previous MG20 and MG20 T2T proteomes, merges each result with `blast_merge.py`, and reduces it to one counterpart per Gifu T2T gene; the three per-genome pair tables are pasted into `All_species_merged_by_GifuT2T.tsv`, the gene correspondence among the four genomes.
**Entry point**: `blast.sh`; `blast_merge.py` also runs standalone.
**Inputs**: in the working directory, `Lotus_GifuT2T.pep.fa` (query of all three runs and subject of the MG20T2T run), `Lotus_GifuOld.pep.fa`, `Lotus_MG20Old.pep.fa`, `Lotus_MG20T2T.pep.fa`; `blast_merge.py`, invoked as `$cul/blast_merge.py` where `$cul` is the working directory; `ALL_GENE_LIST`, a hardcoded path read for column 1 only, used to emit an `NA` row for genes without a counterpart.
**Outputs**: `GifuT2T_GifuOld/`, `GifuT2T_MG20Old/`, `GifuT2T_MG20T2T/`, each holding the BLAST database `db.*`, `result.blast`, `result.blast.res`, `ref_len.txt`, `query_len.txt`, `blast1.merge` (the literal file name used by the script) and intermediate tables up to `final_pairs.tsv`; `All_species_merged_by_GifuT2T.tsv` (header `GifuT2T GifuOld MG20Old MG20T2T`, one row per gene).
**Run**: `bash blast.sh`; `ALL_GENE_LIST` is an absolute path outside this repository in the script and must be pointed at the real file first; standalone use is `python3 blast_merge.py <blast_result> <query_length_file> <ref_length_file> <output_file> [chunksize]` with chunksize defaulting to 2,000,000.
**Tools**: BLAST+ `makeblastdb -dbtype prot` and `blastp` (`-evalue 1e-5`, `-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs"`, `-num_threads 88`); Python 3 + pandas (`blast_merge.py`, chunked read `chunksize=2_000_000`); seqkit `fx2tab -l -n -i`; GNU Awk keeping `$3 > 50` (sure pairs `$3==100 && $7==1 && $8==1`); `cut`, `sort`, `paste`, `find`.

</details>


