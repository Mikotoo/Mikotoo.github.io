---
title: '重复单元'
sidebar:
  label: '重复单元'
  order: 2
---


## 中文

**作用**：把 TRF 记录变成每条染色体的候选着丝粒单体序列并跨染色体连接：对每个重复单元做 canonical 化（旋转与反向互补）、选出最能解释该家族 k×L 长度结构的单体长度 L，再自比对单体得到染色体对连接数与连通分量。
**入口**：`run.sh`
**输入**：`gifu_cent.bed`、`mg20_cent.bed`（TRF 记录，单元序列在指定列）；`gifu_chr_top_monomer.tsv`、`mg20_chr_top_monomer.tsv`（chr、单元长度、单元序列、span）；`cluster.py` 用 `self.sim90.tsv`。
**输出**：`run.sh` 写 `lotus_chr_top_monomer.fa`、`lotus_monomer.*` BLAST 库、`self.allvsall.tsv`、`self.noSelf.tsv`、`self.sim80.tsv`、`self.sim90.tsv`、`top_hit_per_query.tsv`、`chr_pair_links.tsv`；`canon_unit_and_length.py` 写扩展 TSV 与 `length_summary.tsv`；`pick_monomer_from_clusters.py` 写 `cluster_monomer.tsv`；`cluster.py` 把簇打印到 stdout。
**运行**：`bash run.sh`；`python3 canon_unit_and_length.py -i mg20_cent.bed -o mg20.units.canonical.tsv --summary mg20.length_summary.tsv --seq-col 9`；`python3 pick_monomer_from_clusters.py -i clusters.tsv -o cluster_monomer.tsv --tol 6 --max-k 6 --top 5`；`python3 cluster.py`（`run.sh` 第 62 行以上为注释记录需单独运行；`cluster.py` 硬编码 `self.sim90.tsv`，须在本目录运行）。
**工具**：BLAST+ `makeblastdb -dbtype nucl`、`blastn -task blastn-short -dust no -soft_masking false -word_size 7 -evalue 1e-10 -max_target_seqs 10000`，`outfmt 6`；awk/sort 过滤 `pident>=80 && len>=100`（`self.sim80.tsv`）与 `pident>=90 && len>=140`（`self.sim90.tsv`）；三个 Python 脚本仅用标准库；`#CSUB` 队列 c01、88 槽。

<details>
<summary>English</summary>

**Purpose**: Turn TRF records into a per-chromosome set of candidate centromeric monomer sequences and link them across chromosomes: canonicalise each repeat unit (rotation and reverse complement), pick the monomer length L that best explains the family's k×L length structure, and self-align the monomers to build chromosome-pair link counts and connected components.
**Entry point**: `run.sh`
**Inputs**: `gifu_cent.bed`, `mg20_cent.bed` (TRF records, unit sequence in the configured column); `gifu_chr_top_monomer.tsv`, `mg20_chr_top_monomer.tsv` (chr, unit length, unit sequence, span); `self.sim90.tsv` for `cluster.py`.
**Outputs**: `run.sh` writes `lotus_chr_top_monomer.fa`, the `lotus_monomer.*` BLAST database, `self.allvsall.tsv`, `self.noSelf.tsv`, `self.sim80.tsv`, `self.sim90.tsv`, `top_hit_per_query.tsv` and `chr_pair_links.tsv`; `canon_unit_and_length.py` writes the extended TSV and `length_summary.tsv`; `pick_monomer_from_clusters.py` writes `cluster_monomer.tsv`; `cluster.py` prints clusters to stdout.
**Run**: `bash run.sh`; `python3 canon_unit_and_length.py -i mg20_cent.bed -o mg20.units.canonical.tsv --summary mg20.length_summary.tsv --seq-col 9`; `python3 pick_monomer_from_clusters.py -i clusters.tsv -o cluster_monomer.tsv --tol 6 --max-k 6 --top 5`; `python3 cluster.py` (everything above line 62 of `run.sh` is a commented record to be run separately; `cluster.py` hardcodes `self.sim90.tsv` and must run in this directory).
**Tools**: BLAST+ `makeblastdb -dbtype nucl`, `blastn -task blastn-short -dust no -soft_masking false -word_size 7 -evalue 1e-10 -max_target_seqs 10000` with `outfmt 6`; awk/sort filters `pident>=80 && len>=100` (`self.sim80.tsv`) and `pident>=90 && len>=140` (`self.sim90.tsv`); all three Python scripts are standard library only; `#CSUB` queue c01, 88 slots.

</details>

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>canon_unit_and_length.py</code> · Python · 120 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/05_centromere/02_repeat_units/canon_unit_and_length.py)

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from collections import defaultdict

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

def revcomp(seq: str) -> str:
    return seq.translate(COMP)[::-1]

def canonical_rotation(seq: str) -> str:
    """Return lexicographically smallest rotation of seq."""
    if not seq:
        return seq
    # O(n^2) but n is small (<~200) for satellite monomers, ok.
    rots = (seq[i:] + seq[:i] for i in range(len(seq)))
    return min(rots)

def canonical_unit_strand_invariant(seq: str) -> str:
    """Canonicalize by rotation, and also consider reverse-complement."""
    s = seq.upper()
    c1 = canonical_rotation(s)
    rc = revcomp(s)
    c2 = canonical_rotation(rc)
    return min(c1, c2)

def is_header(line: str) -> bool:
    # ignore TRF header-like lines if present
    low = line.lower()
    return (low.startswith("sequence:") or low.startswith("parameters:"))

def main():
    ap = argparse.ArgumentParser(
        description="Compute canonical unit (rotation + reverse-complement) and unit length from TRF-derived table."
    )
    ap.add_argument("-i", "--input", required=True, help="Input TSV (or whitespace-separated) file.")
    ap.add_argument("-o", "--output", required=True, help="Output TSV with canonical unit columns added.")
    ap.add_argument("--seq-col", type=int, default=4,
                    help="1-based column index for unit sequence (default: 4).")
    ap.add_argument("--start-col", type=int, default=2,
                    help="1-based column index for start (default: 2).")
    ap.add_argument("--end-col", type=int, default=3,
                    help="1-based column index for end (default: 3).")
    ap.add_argument("--summary", default="length_summary.tsv",
                    help="Output summary TSV filename (default: length_summary.tsv).")
    args = ap.parse_args()

    seq_idx = args.seq_col - 1
    s_idx = args.start_col - 1
    e_idx = args.end_col - 1

    # length -> count, total_span
    len_count = defaultdict(int)
    len_span = defaultdict(int)

    with open(args.input, "r", encoding="utf-8") as fin, \
         open(args.output, "w", encoding="utf-8") as fout:

        # 写表头
        fout.write("\t".join([
            "chrom", "start", "end", "unit",
            "unit_len", "canon_unit", "canon_len",
            "span_bp"
        ]) + "\t" + "raw_line" + "\n")

        for line in fin:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if is_header(line):
                continue

            # 允许 TSV 或空格分隔
            parts = line.split("\t")
            if len(parts) == 1:
                parts = line.split()

            # 跳过明显不是数据行的情况
            if len(parts) <= max(seq_idx, s_idx, e_idx):
                continue

            chrom = parts[0]
            try:
                start = int(parts[s_idx])
                end = int(parts[e_idx])
            except ValueError:
                # 可能是表头行
                continue

            unit = parts[seq_idx].strip()
            # 若 unit 中混入了非ACGTN字符（一般不会），这里简单清理
            unit_clean = "".join([c for c in unit.upper() if c in "ACGTN"])
            if not unit_clean:
                continue

            unit_len = len(unit_clean)
            canon = canonical_unit_strand_invariant(unit_clean)
            span = end - start + 1
            if span < 0:
                span = 0

            len_count[unit_len] += 1
            len_span[unit_len] += span

            # 输出：固定字段 + 原始整行（方便你回溯）
            fout.write("\t".join(map(str, [
                chrom, start, end, unit_clean,
                unit_len, canon, len(canon),
                span
            ])) + "\t" + line + "\n")

    # 输出长度统计
    with open(args.summary, "w", encoding="utf-8") as fs:
        fs.write("unit_len\tcount\ttotal_span_bp\n")
        for L in sorted(len_count.keys()):
            fs.write(f"{L}\t{len_count[L]}\t{len_span[L]}\n")

if __name__ == "__main__":
    main()
```

</details>

<details open>
<summary><code>cluster.py</code> · Python · 26 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/05_centromere/02_repeat_units/cluster.py)

```python
import sys
from collections import defaultdict, deque

edges = defaultdict(set)
for line in open("self.sim90.tsv"):
    q,s,*_ = line.rstrip("\n").split("\t")
    edges[q].add(s)
    edges[s].add(q)

seen=set()
cid=0
for node in edges:
    if node in seen: 
        continue
    cid += 1
    dq=deque([node])
    seen.add(node)
    comp=[]
    while dq:
        x=dq.popleft()
        comp.append(x)
        for y in edges[x]:
            if y not in seen:
                seen.add(y)
                dq.append(y)
    print(f"Cluster{cid}\t" + "\t".join(sorted(comp)))
```

</details>

<details open>
<summary><code>pick_monomer_from_clusters.py</code> · Python · 118 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/05_centromere/02_repeat_units/pick_monomer_from_clusters.py)

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
from collections import defaultdict, Counter

LEN_RE = re.compile(r"\|L=(\d+)\|")

def parse_len(seq_id: str) -> int:
    m = LEN_RE.search(seq_id)
    if not m:
        raise ValueError(f"Cannot parse length from id: {seq_id}")
    return int(m.group(1))

def score_candidate(L, lengths, tol=6, max_k=6):
    """
    Score how well a candidate monomer length L explains the observed lengths.
    A length x is considered explained if exists k s.t. |x - k*L| <= tol.
    Return: (explained_count, explained_weighted, details)
    """
    explained = 0
    explained_w = 0.0
    detail = []
    for x in lengths:
        best = None
        for k in range(1, max_k + 1):
            target = k * L
            d = abs(x - target)
            if d <= tol:
                if best is None or d < best[0]:
                    best = (d, k, target)
        if best:
            explained += 1
            # closer match contributes more
            explained_w += 1.0 / (1.0 + best[0])
            detail.append((x, best[1], best[2], best[0]))
        else:
            detail.append((x, None, None, None))
    return explained, explained_w, detail

def main():
    ap = argparse.ArgumentParser(
        description="Pick shortest monomer per cluster by explaining k×L length structure."
    )
    ap.add_argument("-i", "--input", required=True, help="clusters.tsv (cluster_id\\tmembers...)")
    ap.add_argument("-o", "--output", default="cluster_monomer.tsv", help="Output TSV")
    ap.add_argument("--tol", type=int, default=6, help="Tolerance in bp for k×L matching (default: 6)")
    ap.add_argument("--max-k", type=int, default=6, help="Max multiplier k to try (default: 6)")
    ap.add_argument("--top", type=int, default=5, help="Report top N monomer-length member IDs closest to L* (default: 5)")
    args = ap.parse_args()

    # cluster -> list of member IDs
    clusters = {}
    with open(args.input, "r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        # accept either: cluster_id size members  OR your pasted format cluster_id members...
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            cid = parts[0]
            if len(parts) == 2:
                members = parts[1].split(",") if "," in parts[1] else parts[1].split()
            elif len(parts) >= 3 and parts[1].isdigit():
                members = parts[2].split(",")
            else:
                # fallback: everything after first col is member IDs (space separated)
                members = []
                for p in parts[1:]:
                    members.extend(p.split(","))
            members = [m for m in members if m]
            clusters[cid] = members

    with open(args.output, "w", encoding="utf-8") as out:
        out.write("cluster\tbest_monomer_len\tn_members\tn_distinct_lens\texplained_count\texplained_weighted\tlength_histogram\ttop_monomer_like_ids\n")

        for cid, members in sorted(clusters.items()):
            lens = [parse_len(m) for m in members]
            lens_sorted = sorted(lens)
            hist = Counter(lens_sorted)

            # candidate monomer lengths: use observed lengths in a reasonable range
            # for centromere satellites, monomer often 80–250, but keep broad
            candidates = sorted({L for L in hist if 60 <= L <= 250})

            if not candidates:
                # fallback: take min observed as candidate
                candidates = [min(lens_sorted)]

            best = None
            best_detail = None
            for L in candidates:
                expl, expl_w, detail = score_candidate(L, lens_sorted, tol=args.tol, max_k=args.max_k)
                # primary: maximize explained_count, secondary: explained_weighted, tertiary: smaller L
                key = (expl, expl_w, -L)
                if best is None or key > best[0]:
                    best = (key, L, expl, expl_w)
                    best_detail = detail

            _, bestL, expl, expl_w = best

            # pick member IDs closest to bestL (k=1)
            # choose those with length within tol of bestL, or just closest
            monomer_like = [(abs(parse_len(m) - bestL), parse_len(m), m) for m in members]
            monomer_like.sort(key=lambda x: (x[0], x[1]))
            top_ids = [x[2] for x in monomer_like[:args.top]]

            # histogram string
            hist_str = ",".join([f"{L}:{c}" for L, c in sorted(hist.items())])

            out.write(f"{cid}\t{bestL}\t{len(members)}\t{len(hist)}\t{expl}\t{expl_w:.3f}\t{hist_str}\t" + ",".join(top_ids) + "\n")

    print(f"[DONE] Wrote: {args.output}")

if __name__ == "__main__":
    main()
```

</details>

<details open>
<summary><code>run.sh</code> · Shell · 99 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/05_centromere/02_repeat_units/run.sh)

```bash
#!/bin/bash
#CSUB -J process
#CSUB -q c01
#CSUB -o process.out
#CSUB -e process.error
#CSUB -n 88
#CSUB -R span[hosts=1]

#python3 canon_unit_and_length.py \
#  -i mg20_cent.bed \
#  -o mg20.units.canonical.tsv \
#  --summary mg20.length_summary.tsv \
#  --seq-col 9

#####
#in="gifu.units.canonical.tsv"   # 改成你的文件名
#out="mg20_chr_length_summary.tsv"

#awk -F'\t' 'NR>1{
#  chr=$1; L=$5; span=$8;
#  key=chr"\t"L;
#  cov[key]+=span; cnt[key]+=1
#}
#END{
#  print "chr\tunit_len\tcount\ttotal_span_bp";
#  for(k in cov) print k"\t"cnt[k]"\t"cov[k]
#}' "$in" | sort -k1,1 -k4,4nr > "$out"


#每条染色体覆盖最大的长度

#awk -F'\t' 'NR>1{key=$1"\t"$5; cov[key]+=$8}
#END{for(k in cov) print k"\t"cov[k]}' "$in" \
#| sort -k1,1 -k3,3nr \
#| awk 'BEGIN{FS=OFS="\t"} {if($1!=p){print; p=$1}}'

#in="mg20.units.canonical.tsv"
#out="mg20_chr_top_monomer.tsv"

#awk -F'\t' 'NR>1{
#  chr=$1; L=$5; unit=$6; span=$8;
#  if(L>=100 && L<=400){
#    key=chr"\t"unit;
#    cov[key]+=span;
#    len[key]=L;
#  }
#}
#END{
#  for(k in cov){
#    split(k,a,"\t");
#    chr=a[1]; unit=a[2];
#    print chr"\t"len[k]"\t"unit"\t"cov[k];
#  }
#}' "$in" | sort -k1,1 -k4,4nr \
#| awk 'BEGIN{FS=OFS="\t"}
#{
#  chr=$1;
#  top[chr]++
#  if(top[chr]<=10) print
#}' > "$out"

cat mg20_chr_top_monomer.tsv gifu_chr_top_monomer.tsv > lotus_chr_top_monomer.tsv

in="lotus_chr_top_monomer.tsv"         # 你的这个输出文件名
fa="lotus_chr_top_monomer.fa"

awk 'BEGIN{FS=OFS="\t"}
{
  chr=$1; L=$2; seq=$3; cov=$4;
  id=chr"|L="L"|cov="cov"|n="NR;
  print ">"id"\n"seq
}' "$in" > "$fa"


makeblastdb -in "$fa" -dbtype nucl -out lotus_monomer
blastn -task blastn-short \
  -query "$fa" -db lotus_monomer \
  -dust no -soft_masking false \
  -word_size 7 \
  -evalue 1e-10 \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
  -max_target_seqs 10000 \
  > self.allvsall.tsv

awk -F'\t' '$1!=$2' self.allvsall.tsv > self.noSelf.tsv
awk -F'\t' '($3>=80 && $4>=100)' self.noSelf.tsv > self.sim80.tsv
awk -F'\t' '($3>=90 && $4>=140)' self.noSelf.tsv > self.sim90.tsv
sort -k1,1 -k12,12nr self.sim90.tsv \
| awk -F'\t' '!seen[$1]++{print $1"\t"$2"\tpid="$3"\taln="$4"\tbits="$12}' \
> top_hit_per_query.tsv

awk -F'\t' '{
  split($1,a,"|"); split($2,b,"|");
  chr1=a[1]; chr2=b[1];
  key=chr1"\t"chr2;
  n[key]++
} END{
  for(k in n) print k"\t"n[k]
}' self.sim90.tsv | sort -k3,3nr > chr_pair_links.tsv
```

</details>
