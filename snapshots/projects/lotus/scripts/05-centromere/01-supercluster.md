---
title: '超级簇合并'
sidebar:
  label: '超级簇合并'
  order: 1
---


## 中文

**作用**：TRF 全基因组找串联重复，每个重复家族取一个代表单元用 CD-HIT-EST 聚类，代表序列用 BLASTN 回贴基因组，命中区间 Jaccard ≥ 0.7 的家族用并查集并入超级簇，只在 100 kb bin 中足够致密时才保留。
**入口**：`trf.sh` → `blast/blast.sh` → `blast/cluster_density.sh`（驱动）
**输入**：`Lotus_GifuT2T_v1.0.fasta`；`Gifu_cluster_reps.fa`；`Gifu_cluster_reps.blast`（第 1–4 列 = cluster、chr、start、end）；`Gifu_cluster_representatives.txt`（第 4 列 = cluster、第 5 列 = 单元序列）。
**输出**：驱动脚本第三个参数指定目录（默认 `supercluster`）下的 `windows.100000.bed`、`clusters.list`、`merged/*.bed`、`pairs.jaccard_ge0.7.tsv`、`cluster2super.tsv`、`super_summary.tsv`、`superbeds/*.bed`、`density/density.long.tsv`、`kept_supers.tsv`、`superunits/*.units.bed`。
**运行**：`bash trf.sh`；`bash blast/blast.sh`；`bash blast/cluster_density.sh Gifu_cluster_reps.blast {PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Lotus_GifuT2T_v1.0.fasta supercluster`（`trf.sh` 与 `blast.sh` 多为注释记录，仅末行有效；`cluster_density.sh` 须在 `blast/` 下运行）。
**工具**：TRF `2 7 7 80 10 50 500 -d -h`；CD-HIT-EST `-c 0.85 -n 3 -d 0 -T 8 -M 16000`；BLAST+ `makeblastdb -dbtype nucl`、`blastn -task blastn -evalue 1e-10 -perc_identity 80 -num_threads 64`；bedtools `makewindows -w 100000`、`merge -d 50`、`jaccard`；阈值 `WIN=100000`、`JAC=0.7`、`DENS_CUT=0.5`、`BIN_KEEP=5`。

<details>
<summary>English</summary>

**Purpose**: TRF calls tandem repeats genome-wide, one representative unit per repeat family is clustered with CD-HIT-EST, the representatives are mapped back with BLASTN, and families whose merged genomic hits reach a Jaccard ≥ 0.7 are merged with Union-Find into superclusters, kept only when dense enough in 100-kb bins.
**Entry point**: `trf.sh` → `blast/blast.sh` → `blast/cluster_density.sh` (driver)
**Inputs**: `Lotus_GifuT2T_v1.0.fasta`; `Gifu_cluster_reps.fa`; `Gifu_cluster_reps.blast` (columns 1–4 = cluster, chr, start, end); `Gifu_cluster_representatives.txt` (column 4 = cluster, column 5 = unit sequence).
**Outputs**: in the directory given as the third argument of the driver (default `supercluster`): `windows.100000.bed`, `clusters.list`, `merged/*.bed`, `pairs.jaccard_ge0.7.tsv`, `cluster2super.tsv`, `super_summary.tsv`, `superbeds/*.bed`, `density/density.long.tsv`, `kept_supers.tsv`, `superunits/*.units.bed`.
**Run**: `bash trf.sh`; `bash blast/blast.sh`; `bash blast/cluster_density.sh Gifu_cluster_reps.blast {PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Lotus_GifuT2T_v1.0.fasta supercluster` (`trf.sh` and `blast.sh` are mostly commented records with only their last line live; `cluster_density.sh` must run from `blast/`).
**Tools**: TRF `2 7 7 80 10 50 500 -d -h`; CD-HIT-EST `-c 0.85 -n 3 -d 0 -T 8 -M 16000`; BLAST+ `makeblastdb -dbtype nucl`, `blastn -task blastn -evalue 1e-10 -perc_identity 80 -num_threads 64`; bedtools `makewindows -w 100000`, `merge -d 50`, `jaccard`; thresholds `WIN=100000`, `JAC=0.7`, `DENS_CUT=0.5`, `BIN_KEEP=5`.

</details>

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>blast/blast.sh</code> · Shell · 29 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/05_centromere/01_supercluster/blast/blast.sh)

```bash
#!/bin/bash
#CSUB -J blast
#CSUB -q c02
#CSUB -o blast.out
#CSUB -e blast.error
#CSUB -n 64
#CSUB -R span[hosts=1]



GENOME={PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Lotus_GifuT2T_v1.0.fasta
DB={PROJ_04_LOTUS_GENOME}/03_annotation/Gifu/02.cent/blast/Gifu_genome_db

#makeblastdb -in "$GENOME" -dbtype nucl -out "$DB"


#blastn \
#  -task blastn \
#  -db "$DB" \
#  -query Gifu_cluster_reps.fa \
#  -evalue 1e-10 \
#  -perc_identity 80 \
#  -dust no \
# -outfmt "6 qseqid sseqid sstart send pident length bitscore" \
#  -num_threads 64 \
# > Gifu_cluster_reps.blast


bash cluster_density.sh Gifu_cluster_reps.blast {PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Lotus_GifuT2T_v1.0.fasta supercluster
```

</details>

<details open>
<summary><code>blast/build_supercluster_unionfind.py</code> · Python · 114 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/05_centromere/01_supercluster/blast/build_supercluster_unionfind.py)

```python
#!/usr/bin/env python3
import argparse
from pathlib import Path


def parse_args():
    ap = argparse.ArgumentParser(
        description="Build superclusters with Union-Find from a pairs table (A,B,Jaccard)."
    )
    ap.add_argument("--outdir", required=True, help="Working directory")
    ap.add_argument(
        "--pairs",
        default=None,
        help="Pairs TSV path. Default: <outdir>/pairs.jaccard_ge<JAC>.tsv if --jac is given, else <outdir>/pairs.jaccard.tsv",
    )
    ap.add_argument(
        "--clusters",
        default=None,
        help="Clusters list file. Default: <outdir>/clusters.list",
    )
    ap.add_argument(
        "--jac",
        default=None,
        help="Jaccard cutoff label used to infer default pairs filename (e.g. 0.7 -> pairs.jaccard_ge0.7.tsv)",
    )
    ap.add_argument(
        "--map-out",
        default="cluster2super.tsv",
        help="Output mapping filename under outdir",
    )
    ap.add_argument(
        "--summary-out",
        default="super_summary.tsv",
        help="Output summary filename under outdir",
    )
    ap.add_argument(
        "--prefix",
        default="super_",
        help="Supercluster ID prefix (default: super_)",
    )
    ap.add_argument(
        "--pad",
        type=int,
        default=4,
        help="Zero-padding width for supercluster IDs (default: 4 -> super_0001)",
    )
    return ap.parse_args()


args = parse_args()
outdir = Path(args.outdir)

if args.pairs:
    pairs = Path(args.pairs)
else:
    if args.jac is not None:
        pairs = outdir / f"pairs.jaccard_ge{args.jac}.tsv"
    else:
        pairs = outdir / "pairs.jaccard.tsv"

clusters_file = Path(args.clusters) if args.clusters else (outdir / "clusters.list")
clusters = [x.strip() for x in clusters_file.read_text().splitlines() if x.strip()]

parent = {c: c for c in clusters}
rank = {c: 0 for c in clusters}

def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x

def union(a, b):
    ra, rb = find(a), find(b)
    if ra == rb:
        return
    if rank[ra] < rank[rb]:
        parent[ra] = rb
    elif rank[ra] > rank[rb]:
        parent[rb] = ra
    else:
        parent[rb] = ra
        rank[ra] += 1

if pairs.exists():
    for line in pairs.read_text().splitlines():
        if not line.strip():
            continue
        a, b, _ = line.split("\t")
        union(a, b)

groups = {}
for c in clusters:
    r = find(c)
    groups.setdefault(r, []).append(c)

roots = sorted(groups.keys())

map_out = outdir / args.map_out
summary_out = outdir / args.summary_out

with open(map_out, "w") as f:
    for i, r in enumerate(roots, 1):
        sid = f"{args.prefix}{i:0{args.pad}d}"
        for c in sorted(groups[r]):
            f.write(f"{c}\t{sid}\n")

with open(summary_out, "w") as f:
    for i, r in enumerate(roots, 1):
        f.write(
            f"{args.prefix}{i:0{args.pad}d}\t{len(groups[r])}\t{','.join(sorted(groups[r]))}\n"
        )

print(f"OK: wrote {map_out.name} and {summary_out.name}")
```

</details>

<details open>
<summary><code>blast/cluster_density.sh</code> · Shell · 156 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/05_centromere/01_supercluster/blast/cluster_density.sh)

```bash
#!/usr/bin/env bash
set -euo pipefail

BLAST="$1"      # trf-blast 输出（至少包含：cluster_id chr start end ...）
FASTA="$2"      # 参考基因组 fasta
OUTDIR="${3:-supercluster_0.7_100kb}"

# ===================== 参数（都从这里控制，python 不写死） =====================
WIN=100000          # bin 大小（100kb）
JAC=0.7             # jaccard 阈值
MERGE_D=50          # merge 间隔（用于计算 super union 区间 & jaccard）
DENS_CUT=0.5        # density 阈值（>0.5 认为高密度）
BIN_KEEP=5          # 一个 super 中 density>DENS_CUT 的 bin 数量必须 > BIN_KEEP 才保留

mkdir -p "$OUTDIR"/{merged,superbeds,superunits,density,tmp}

# 取 blast 中某个 cluster 的所有命中区间（不 merge），输出 BED
blast_units_for_cluster() {
  local cl="$1"
  awk -v C="$cl" 'BEGIN{OFS="\t"}
    $1==C{
      chr=$2; s=$3; e=$4;
      if(s>e){t=s;s=e;e=t}
      s0=s-1; if(s0<0)s0=0;
      print chr,s0,e
    }' "$BLAST"
}

echo "[1/7] Prepare genome sizes & windows (${WIN}bp)"
samtools faidx "$FASTA"
cut -f1,2 "${FASTA}.fai" > "$OUTDIR/genome.sizes"
bedtools makewindows -g "$OUTDIR/genome.sizes" -w "$WIN" \
  | sort -k1,1 -k2,2n \
  > "$OUTDIR/windows.${WIN}.bed"

echo "[2/7] Split BLAST by cluster and merge hits (only for Jaccard/super union)"
cut -f1 "$BLAST" | sort -u > "$OUTDIR/clusters.list"

while read -r cl; do
  blast_units_for_cluster "$cl" \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - -d "$MERGE_D" \
  > "$OUTDIR/merged/${cl}.bed"
done < "$OUTDIR/clusters.list"

echo "[3/7] Compute pairwise Jaccard and keep >= ${JAC}"
: > "$OUTDIR/pairs.jaccard_ge${JAC}.tsv"
mapfile -t CLS < "$OUTDIR/clusters.list"
N=${#CLS[@]}

for ((i=0;i<N;i++)); do
  for ((j=i+1;j<N;j++)); do
    A="${CLS[i]}"
    B="${CLS[j]}"
    J=$(bedtools jaccard -a "$OUTDIR/merged/${A}.bed" -b "$OUTDIR/merged/${B}.bed" \
        | awk 'NR==2{print $3}')
    awk -v a="$A" -v b="$B" -v j="$J" -v cut="$JAC" 'BEGIN{
      if((j+0) >= cut) print a"\t"b"\t"j
    }' >> "$OUTDIR/pairs.jaccard_ge${JAC}.tsv"
  done
done

echo "[4/7] Union-Find to build superclusters"
python build_supercluster_unionfind.py \
  --outdir "$OUTDIR" \
  --jac "$JAC" \
  --pairs "$OUTDIR/pairs.jaccard_ge${JAC}.tsv"

echo "[5/7] Build supercluster union BEDs (merge across member clusters)"
cut -f2 "$OUTDIR/cluster2super.tsv" | sort -u > "$OUTDIR/supers.list"

while read -r sup; do
  awk -v S="$sup" '$2==S{print $1}' "$OUTDIR/cluster2super.tsv" \
  | while read -r cl; do
      cat "$OUTDIR/merged/${cl}.bed"
    done \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - -d "$MERGE_D" \
  > "$OUTDIR/superbeds/${sup}.bed"
done < "$OUTDIR/supers.list"

echo "[6/7] Compute density per bin (col5 = density)"
# 输出列：chr start end group densitybin_len
calc_one () {
  local NAME="$1"
  local BED="$2"
  bedtools coverage -a "$OUTDIR/windows.${WIN}.bed" -b "$BED" \
    | awk -v G="$NAME" 'BEGIN{OFS="\t"}
        {
          # bedtools coverage 默认输出：
          # $1 chr, $2 start, $3 end,
          # $4 numFeaturesInA, $5 lenA,
          # $6 basesCovered, $7 fractionCovered
          density=$7
          covered_bp=$6
          bin_len=$5
          print $1,$2,$3,G,density,covered_bp,bin_len
        }'
}


# ALL：所有 super union 再合并
cat "$OUTDIR"/superbeds/super_*.bed \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - -d "$MERGE_D" \
  > "$OUTDIR/superbeds/ALL.bed"

calc_one "ALL" "$OUTDIR/superbeds/ALL.bed" > "$OUTDIR/density/density.ALL.tsv"

: > "$OUTDIR/density/density.supers.tsv"
while read -r sup; do
  calc_one "$sup" "$OUTDIR/superbeds/${sup}.bed" >> "$OUTDIR/density/density.supers.tsv"
done < "$OUTDIR/supers.list"

cat "$OUTDIR/density/density.ALL.tsv" "$OUTDIR/density/density.supers.tsv" \
  > "$OUTDIR/density/density.long.tsv"

echo "[7/7] Filter supers; export unit beds from BLAST"
# 统计每个 super：density>DENS_CUT 的 bin 数
awk -v cut="$DENS_CUT" '
  $4!="ALL" && $5>cut {cnt[$4]++}
  END{
    for(g in cnt){
      if(cnt[g] > '"$BIN_KEEP"') print g"\t"cnt[g]
    }
  }' "$OUTDIR/density/density.supers.tsv" \
  | sort -k1,1 \
  > "$OUTDIR/kept_supers.tsv"

cut -f1 "$OUTDIR/kept_supers.tsv" > "$OUTDIR/kept_supers.list" || true

KEPT_N=$(wc -l < "$OUTDIR/kept_supers.list" 2>/dev/null || echo 0)
echo "[INFO] Kept supers: ${KEPT_N} (criterion: bins(density>${DENS_CUT}) > ${BIN_KEEP})"

if [[ "${KEPT_N}" -eq 0 ]]; then
  echo "[WARN] No super passed filter."
  echo "Done."
  exit 0
fi

# 导出 unit bed：从 BLAST 中提取（不使用 merge 后结果）
while read -r sup; do
  out="$OUTDIR/superunits/${sup}.units.bed"
  : > "$out"
  awk -v S="$sup" '$2==S{print $1}' "$OUTDIR/cluster2super.tsv" \
    | while read -r cl; do
        blast_units_for_cluster "$cl" >> "$out"
      done
  sort -k1,1 -k2,2n "$out" -o "$out"
done < "$OUTDIR/kept_supers.list"

echo "Done."
echo "  Mapping:       $OUTDIR/cluster2super.tsv"
echo "  Density long:  $OUTDIR/density/density.long.tsv"
echo "  Kept list:     $OUTDIR/kept_supers.list"
echo "  Units bed:     $OUTDIR/superunits/*.units.bed"
```

</details>

<details open>
<summary><code>trf.sh</code> · Shell · 145 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/05_centromere/01_supercluster/trf.sh)

```bash
#!/bin/bash
#CSUB -J trf
#CSUB -q c01
#CSUB -o trf.out
#CSUB -e trf.error
#CSUB -n 64
#CSUB -R span[hosts=1]

#trf Gifu_v1.0.fasta 2 7 7 80 10 50 500 -d -h

#awk -f trf_dat2bed.awk Gifu_v1.0.fasta.2.7.7.80.10.50.500.dat > Gifu_TRF_all.bed

#awk '($7 >= 10) && ($3-$2 >= 1000)' Gifu_TRF_all.bed > Gifu_TRF_satellite.bed

#sort -k1,1 -k2,2n -k7,7n Gifu_TRF_satellite.bed \
#| awk '!seen[$1":"$2":"$3]++' \
#> Gifu_TRF_minUnit.bed


# 提取 unique unit
#cut -f9 Gifu_TRF_minUnit.bed | sort -u > units.unique.txt

# 转 fasta
#awk '{print ">"$1"\n"$1}' units.unique.txt > units.fa

# cd-hit-est 聚类
#cd-hit-est -i units.fa -o units_clustered -c 0.85 -n 3 -d 0 -T 8 -M 16000

# 生成 unit -> cluster 映射
#awk '
#BEGIN{ OFS="\t" }
#/^>Cluster/{
#    cluster = $2;
#    next;
#}
#{
#    match($0, />[^ ]+/);
#    if (RSTART > 0) {
#        seq = substr($0, RSTART+1, RLENGTH-1);
#        sub(/\.\.\.$/, "", seq);
#        print seq, cluster;
#    }
#}' units_clustered.clstr > unit2cluster.txt

# 加回 TRF 记录
#awk '
#BEGIN{ OFS="\t" }
#NR==FNR { unit2cl[$1]=$2; next }
#{
#    unit = $9;
#    cl = unit2cl[unit];
#    if (cl=="") cl="NA";
#    print $0, cl;
#}
#' unit2cluster.txt Gifu_TRF_minUnit.bed \
#  > Gifu_TRF_satellite_unitCluster.bed



#awk '
#BEGIN{ OFS="\t" }
#{
#    chr  = $1
#    beg  = $2
#    end  = $3
#    cl   = $NF      # 最后一列视为 cluster
#    if (cl == "NA") next    # 不要 NA 的 cluster（如果有）

#    len = end - beg
#    key = chr "\t" cl
#    a[key] += len
#}
#END{
#    for (k in a) {
#        print k, a[k]   # chr cluster total_len
#    }
#}' Gifu_TRF_satellite_unitCluster.bed \
#  > Gifu_cluster_len_per_chr.txt

#sort -k1,1 -k 3,3nr Gifu_cluster_len_per_chr.txt > Gifu_cluster_len_per_chr.sorted.txt


#TOPN=10
#awk -v N=$TOPN '
#BEGIN{ OFS="\t" }
#{
#    chr = $1
#    if (chr != last_chr) {
#        last_chr = chr
#        count = 0
#    }
#    if (count < N) {
#        print
#        count++
#    }
#}' Gifu_cluster_len_per_chr.sorted.txt \
#  > Gifu_top${TOPN}_clusters_per_chr.txt


#cut -f2 Gifu_top${TOPN}_clusters_per_chr.txt \
#  | sort -u \
#  > Gifu_candidate_clusters.txt


#awk '
#BEGIN{ OFS="\t" }
#{
#    chr  = $1
#    beg  = $2
#    end  = $3
#    cl   = $10
#    unit = $9
#    len  = end - beg
#    print chr, beg, end, cl, unit, len
#}' Gifu_TRF_satellite_unitCluster.bed \
#  | sort -k4,4 -k6,6nr \
#  > Gifu_unit_by_cluster_len.txt


#awk '
#BEGIN{ OFS="\t" }
#NR==FNR{
#    cand[$1] = 1;   # 候选 cluster 集合
#    next;
#}
#{
#    cl = $4
#    if (!(cl in cand)) next
#    if (!(cl in seen)) {
#        seen[cl] = 1
#        print
#    }
#}' Gifu_candidate_clusters.txt \
#   Gifu_unit_by_cluster_len.txt \
# > Gifu_cluster_representatives.txt

awk '
BEGIN{ OFS="" }
{
    cl   = $4
    unit = $5
    print ">cluster_", cl
    print unit
}' Gifu_cluster_representatives.txt \
  > blast/Gifu_cluster_reps.fa
```

</details>
