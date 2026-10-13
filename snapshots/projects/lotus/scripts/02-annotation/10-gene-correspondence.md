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

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>blast.sh</code> · Shell · 163 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/02_annotation/10_gene_correspondence/blast.sh)

```bash
#!/bin/bash
#CSUB -J blast
#CSUB -q c01
#CSUB -o blast.out
#CSUB -e blast.error
#CSUB -n 88
#CSUB -R span[hosts=1]

set -euo pipefail

#=================== 参数设置 ===================
THREADS=88
cul=$(pwd)

#=================== 输入文件 ===================
GifuOld="Lotus_GifuOld.pep.fa"
GifuT2T="Lotus_GifuT2T.pep.fa"
MG20Old="Lotus_MG20Old.pep.fa"
MG20T2T="Lotus_MG20T2T.pep.fa"

ALL_GENE_LIST="{PROJ_LOTUS_HX}/03_annotation/Ortholog/04.blast_Ortholog3/GifuT2T_all_gene_list.txt"

#=================== 函数定义 ===================
run_blastp() {
    query=$1
    subject=$2
    outdir=$3

    mkdir -p "$outdir"

    makeblastdb -in "$subject" -dbtype prot -out "$outdir/db"

    blastp -query "$query" -db "$outdir/db" \
        -out "$outdir/result.blast" \
        -evalue 1e-5 \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
        -num_threads $THREADS
}

#=================== BLAST 比对 ===================
run_blastp "$GifuT2T" "$GifuOld"  "GifuT2T_GifuOld"
run_blastp "$GifuT2T" "$MG20Old"  "GifuT2T_MG20Old"
run_blastp "$GifuT2T" "$MG20T2T"  "GifuT2T_MG20T2T"

echo "✅ BLASTP 完成"

#=================== 结果处理 ===================
find ./GifuT2T_* -type f -name "*.blast" | while read blast; do

    blast1=$(basename "$blast")
    blast2=$(basename "$(dirname "$blast")")

    ref=${blast2##*_}
    query=${blast2%%_*}

    workdir="$cul/$blast2"
    cd "$workdir"

    # ------------------------------------------------
    # step0: 基础整理
    # ------------------------------------------------
    cut -f1,2,4,3 "$blast1" > "${blast1}.res"

    seqkit fx2tab -l -n -i "$cul/Lotus_${ref}.pep.fa"   > ref_len.txt
    seqkit fx2tab -l -n -i "$cul/Lotus_${query}.pep.fa" > query_len.txt

    python3 $cul/blast_merge.py \
        "${blast1}.res" \
        query_len.txt \
        ref_len.txt \
        blast1.merge

    # ------------------------------------------------
    # 1️⃣ 去版本号 + 基础过滤 → tmp1
    # ------------------------------------------------
    awk -F'\t' '
    BEGIN{OFS="\t"}
    {
        gsub(/\.[0-9]+$/, "", $1)
        gsub(/\.[0-9]+$/, "", $2)
        if ($3 > 50)
            print
    }' blast1.merge > tmp1.tsv

    # ------------------------------------------------
    # 2️⃣ 提取 sure pairs
    #   100% identity + 完全覆盖
    # ------------------------------------------------
    awk -F'\t' '
    BEGIN{OFS="\t"}
    $3==100 && $7==1 && $8==1 {print $1,$2}
    ' tmp1.tsv | sort -u > sure_pairs.tsv

    # ------------------------------------------------
    # 3️⃣ 去除已是 sure 的 query → tmp2
    # ------------------------------------------------
    awk -F'\t' '
    NR==FNR {a[$1]; next}
    !($1 in a)
    ' sure_pairs.tsv tmp1.tsv > tmp2.tsv

    # ------------------------------------------------
    # 4️⃣ 对剩余 query，按 $9 最大选 best hit
    # ------------------------------------------------
    awk -F'\t' '
    BEGIN{OFS="\t"}
    {
        if (!($1 in max) || $9 > max[$1]) {
            max[$1] = $9
            line[$1] = $1 OFS $2
        }
    }
    END{
        for (i in line) print line[i]
    }' tmp2.tsv > best_pairs.tsv

    # ------------------------------------------------
    # 5️⃣ 合并 sure + best
    # ------------------------------------------------
    cat sure_pairs.tsv best_pairs.tsv | sort -u > all_pairs.tsv

    # ------------------------------------------------
    # 6️⃣ 补全所有 GifuT2T 基因
    # ------------------------------------------------
    awk -F'\t' '
    BEGIN{OFS="\t"}
    NR==FNR {map[$1]=$2; next}
    {
        if ($1 in map)
            print $1, map[$1]
        else
            print $1, "NA"
    }' all_pairs.tsv "$ALL_GENE_LIST" > final_pairs.tsv

    echo "✅ [$blast2] finished → final_pairs.tsv"

    cd "$cul"
done

echo "🎉 所有比对与基因配对流程完成"

#================================================================
# 7️⃣ 合并所有 final_pairs.tsv → 四列表
#================================================================
echo "🔗 合并所有 final_pairs.tsv → All_species_merged_by_GifuT2T.tsv"

paste \
    "$cul/GifuT2T_GifuOld/final_pairs.tsv" \
    "$cul/GifuT2T_MG20Old/final_pairs.tsv" \
    "$cul/GifuT2T_MG20T2T/final_pairs.tsv" \
| awk -F'\t' '
BEGIN{
    OFS="\t";
    print "GifuT2T","GifuOld","MG20Old","MG20T2T"
}
{
    # 每个 final_pairs.tsv 都是两列：GifuT2T  Other
    # paste 后为：1 2 | 3 4 | 5 6
    print $1, $2, $4, $6
}
' > All_species_merged_by_GifuT2T.tsv

echo "🎉 所有物种合并完成 → All_species_merged_by_GifuT2T.tsv"
```

</details>

<details open>
<summary><code>blast_merge.py</code> · Python · 102 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/02_annotation/10_gene_correspondence/blast_merge.py)

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pandas as pd
from collections import defaultdict


def load_length_map(length_file: str, key_col: str, val_col: str) -> dict:
    df = pd.read_csv(length_file, sep="\t", header=None, usecols=[0, 1])
    df.columns = [key_col, val_col]
    df[key_col] = df[key_col].astype(str)
    df[val_col] = df[val_col].astype(float)
    return dict(zip(df[key_col], df[val_col]))


def process_blast_results(input_file, query_lengths_file, reference_lengths_file, output_file,
                          chunksize=2_000_000):

    qlen_map = load_length_map(query_lengths_file, "query_id", "query_length")
    rlen_map = load_length_map(reference_lengths_file, "reference_id", "reference_length")

    # (q, r) -> 累加值
    acc_sim_len = defaultdict(float)   # Σ(similarity * aln_len)
    acc_len = defaultdict(float)       # Σ(aln_len)

    dtypes = {0: "string", 1: "string", 2: "float32", 3: "float32"}

    reader = pd.read_csv(
        input_file,
        sep="\t",
        header=None,
        usecols=[0, 1, 2, 3],
        dtype=dtypes,
        chunksize=chunksize
    )

    for chunk in reader:
        chunk.columns = ["query_id", "reference_id", "similarity", "alignment_length"]
        chunk["sim_x_len"] = chunk["similarity"] * chunk["alignment_length"]

        g = chunk.groupby(["query_id", "reference_id"], sort=False, observed=True).agg(
            sum_sim_x_len=("sim_x_len", "sum"),
            sum_len=("alignment_length", "sum")
        )

        for (q, r), row in g.iterrows():
            key = (str(q), str(r))
            acc_sim_len[key] += float(row["sum_sim_x_len"])
            acc_len[key] += float(row["sum_len"])

    # ===== 输出 =====
    with open(output_file, "w") as out:
        for (q, r), total_len in acc_len.items():
            if total_len <= 0:
                continue

            sim_sum = acc_sim_len[(q, r)]
            sim = sim_sum / total_len

            qlen = qlen_map.get(q, float("nan"))
            rlen = rlen_map.get(r, float("nan"))

            qcov = total_len / qlen if qlen == qlen else float("nan")
            scov = total_len / rlen if rlen == rlen else float("nan")

            # ===== 新增两列 =====
            q_weighted_sim_cov = sim_sum / qlen if qlen == qlen else float("nan")
            r_weighted_sim_cov = sim_sum / rlen if rlen == rlen else float("nan")

            out.write(
                f"{q}\t{r}\t{sim:.6f}\t{total_len:.2f}\t"
                f"{qlen if qlen == qlen else 'NA'}\t"
                f"{rlen if rlen == rlen else 'NA'}\t"
                f"{qcov if qcov == qcov else 'NA'}\t"
                f"{scov if scov == scov else 'NA'}\t"
                f"{q_weighted_sim_cov if q_weighted_sim_cov == q_weighted_sim_cov else 'NA'}\t"
                f"{r_weighted_sim_cov if r_weighted_sim_cov == r_weighted_sim_cov else 'NA'}\n"
            )


if __name__ == "__main__":
    if len(sys.argv) not in (5, 6):
        print(
            "用法: python3 blast_merge.py "
            "<blast_result> <query_length_file> <ref_length_file> <output_file> [chunksize]"
        )
        sys.exit(1)

    input_file = sys.argv[1]
    query_lengths_file = sys.argv[2]
    reference_lengths_file = sys.argv[3]
    output_file = sys.argv[4]
    chunksize = int(sys.argv[5]) if len(sys.argv) == 6 else 2_000_000

    process_blast_results(
        input_file,
        query_lengths_file,
        reference_lengths_file,
        output_file,
        chunksize=chunksize
    )
```

</details>
