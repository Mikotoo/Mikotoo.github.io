---
title: '一致性序列'
sidebar:
  label: '一致性序列'
  order: 3
---


## 中文

**作用**：为四个着丝粒重复单元 `unit_161`、`unit_172`、`unit_186`、`unit_330` 构建共识序列：从每个基因组提取该单元的命中、用 MAFFT 比对、折叠成共识，再把 Gifu 与 MG20 的共识相互比对得到跨材料共识。
**入口**：`run_unit_mafft_consensus.sh`
**输入**：`Lotus_GifuT2T_v1.0.fasta`、`Lotus_MG20T2T_v1.0.fasta`；每个单元的 BED `gifu_unit_<unit>.bed` 与 `mg20_unit_<unit>.bed`（BED6）。
**输出**：`separate_then_merge_consensus/` 下 `01_seqs/<tag>_<unit>.fa`、`02_msa/<tag>_<unit>.aln.fa`、`03_consensus/<tag>_<unit>.cons.fa`、`03_dist/<tag>_<unit>.dist.tsv`、`04_merge2cons_msa/<unit>.Gifu_MG20.cons.{fa,aln.fa}`、`05_final/<unit>.final.cons.fa` 与 `<unit>.final.dist.tsv`。
**运行**：`bash run_unit_mafft_consensus.sh`（须在含两个基因组 FASTA 与 `gifu_`/`mg20_` 前缀 BED 的目录中运行；`blast.sh` 除末行外均为注释）；单次调用：`python get_consensus.py -i aln.fa -o cons.fa -d dist.tsv -t 0.5 --include-gap --name <name>`。
**工具**：MAFFT `--thread 32 --auto`；bedtools `getfasta -name -s`、`jaccard`、`intersect`、`sort`、`merge`；BLAST+ `makeblastdb -dbtype nucl`、`blastn -task blastn-short -word_size 7 -evalue 1e-10`（`outfmt 6`）；`get_consensus.py` 用 Biopython `AlignIO`，`Levenshtein` 可选（内置纯 Python 回退）；`THRESH=0.5`。

<details>
<summary>English</summary>

**Purpose**: Build consensus sequences for the four centromeric repeat units `unit_161`, `unit_172`, `unit_186` and `unit_330`: hits of each unit are extracted from each genome, aligned with MAFFT and collapsed into a consensus, and the Gifu and MG20 consensuses are then aligned to each other to give a cross-accession consensus.
**Entry point**: `run_unit_mafft_consensus.sh`
**Inputs**: `Lotus_GifuT2T_v1.0.fasta`, `Lotus_MG20T2T_v1.0.fasta`; the per-unit BEDs `gifu_unit_<unit>.bed` and `mg20_unit_<unit>.bed` (BED6).
**Outputs**: under `separate_then_merge_consensus/`: `01_seqs/<tag>_<unit>.fa`, `02_msa/<tag>_<unit>.aln.fa`, `03_consensus/<tag>_<unit>.cons.fa`, `03_dist/<tag>_<unit>.dist.tsv`, `04_merge2cons_msa/<unit>.Gifu_MG20.cons.{fa,aln.fa}`, `05_final/<unit>.final.cons.fa` and `<unit>.final.dist.tsv`.
**Run**: `bash run_unit_mafft_consensus.sh` (must run from a directory holding the two genome FASTAs and the `gifu_`/`mg20_`-prefixed BEDs; `blast.sh` is commented out except its last line); a single call: `python get_consensus.py -i aln.fa -o cons.fa -d dist.tsv -t 0.5 --include-gap --name <name>`.
**Tools**: MAFFT `--thread 32 --auto`; bedtools `getfasta -name -s`, `jaccard`, `intersect`, `sort`, `merge`; BLAST+ `makeblastdb -dbtype nucl`, `blastn -task blastn-short -word_size 7 -evalue 1e-10` (`outfmt 6`); `get_consensus.py` uses Biopython `AlignIO` with `Levenshtein` optional (a pure-Python fallback is built in); `THRESH=0.5`.

</details>

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>blast.sh</code> · Shell · 86 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/05_centromere/03_consensus/blast.sh)

```bash
#!/bin/bash
#CSUB -J process
#CSUB -q c02
#CSUB -o process.out
#CSUB -e process.error
#CSUB -n 64
#CSUB -R span[hosts=1]

#makeblastdb -in Lotus_GifuT2T_v1.0.fasta -dbtype nucl -out gifu_genome
#makeblastdb -in Lotus_MG20T2T_v1.0.fasta -dbtype nucl -out mg20_genome

#blastn -task blastn-short \
#  -query cent_unit.fa -db gifu_genome \
#  -dust no -soft_masking false \
#  -word_size 7 \
#  -evalue 1e-10 \
#  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
#  > gifu_cent.blast

#blastn -task blastn-short \
#  -query cent_unit.fa -db mg20_genome \
#  -dust no -soft_masking false \
#  -word_size 7 \
#  -evalue 1e-10 \
#  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
#  > mg20_cent.blast


#blast=gifu_cent.blast

# 阈值可按需要调：例如 identity>=90 且 aln_len>=120
#PID=80
#ALN=100

#awk -v PID="$PID" -v ALN="$ALN" 'BEGIN{FS=OFS="\t"}
#{
#  q=$1; chr=$2; pid=$3; alen=$4;
#  s1=$9; s2=$10;

#  if(pid < PID) next;
#  if(alen < ALN) next;

  # strand: sstart < send -> + ; else -
#  strand = (s1 < s2 ? "+" : "-");

  # BED coordinates (0-based, half-open)
#  st = (s1 < s2 ? s1 : s2);
#  en = (s1 > s2 ? s1 : s2);

  # output BED6: chr start end name score strand
#  print chr, st-1, en, q, pid, strand
#}' "$blast" \
#| sort -k1,1 -k2,2n > gifu_units.all.bed

# 按 unit 拆成 4 个 bed
#for u in unit_161 unit_172 unit_186 unit_330; do
#  awk -v u=$u '$4==u' gifu_units.all.bed > gifu_${u}.bed
#done


#for a in mg20_unit_161 mg20_unit_172 mg20_unit_186 mg20_unit_330; do
#  for b in mg20_unit_161 mg20_unit_172 mg20_unit_186 mg20_unit_330; do
#    if [[ "$a" < "$b" ]]; then
#      echo "== $a vs $b (ROI) =="
#      bedtools jaccard -a ${a}.bed -b ${b}.bed
#    fi
#  done
#done


# 定义窗口
#echo -e "Chr5\t34000000\t34100000" > chr5_34.0_34.1.bed

# 提取窗口内命中
#for u in gifu_unit_161 gifu_unit_172 gifu_unit_186 gifu_unit_330; do
#  bedtools intersect -a ${u}.bed -b chr5_34.0_34.1.bed > ${u}.Chr5_34_34.1.bed
#done

# 可选：合并碎片，减少视觉“交错假象”
#for u in gifu_unit_161 gifu_unit_172 gifu_unit_186 gifu_unit_330; do
#  bedtools sort -i ${u}.Chr5_34_34.1.bed \
#  | bedtools merge -i - > ${u}.Chr5_34_34.1.merged.bed
#done


bash run_unit_mafft_consensus.sh
```

</details>

<details open>
<summary><code>get_consensus.py</code> · Python · 125 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/05_centromere/03_consensus/get_consensus.py)

```python
from Bio import AlignIO
from collections import Counter
import argparse
import os

try:
    from Levenshtein import distance as levenshtein_distance
except ImportError:
    def levenshtein_distance(s1, s2):
        """纯Python实现的Levenshtein距离计算（适用于短序列）"""
        if len(s1) < len(s2):
            return levenshtein_distance(s2, s1)
        if len(s2) == 0:
            return len(s1)
        
        previous_row = range(len(s2) + 1)
        for i, c1 in enumerate(s1):
            current_row = [i + 1]
            for j, c2 in enumerate(s2):
                insertions = previous_row[j + 1] + 1
                deletions = current_row[j] + 1
                substitutions = previous_row[j] + (c1 != c2)
                current_row.append(min(insertions, deletions, substitutions))
            previous_row = current_row
        return previous_row[-1]

def calculate_consensus(alignment, threshold=0.5, include_gap=False, ambiguity=False):
    """增强型共有序列计算"""
    consensus = []
    valid_bases = set('ACGTUacgtu-') if include_gap else set('ACGTUacgtu')
    iupac_code = {
        frozenset(['A']): 'A', frozenset(['C']): 'C', frozenset(['G']): 'G',
        frozenset(['T']): 'T', frozenset(['U']): 'U', frozenset(['A','G']): 'R',
        frozenset(['C','T']): 'Y', frozenset(['G','C']): 'S', frozenset(['A','T']): 'W',
        frozenset(['G','T']): 'K', frozenset(['A','C']): 'M', frozenset(['C','G','T']): 'B',
        frozenset(['A','G','T']): 'D', frozenset(['A','C','T']): 'H', frozenset(['A','C','G']): 'V',
        frozenset(['A','C','G','T']): 'N'
    }

    for col_idx in range(alignment.get_alignment_length()):
        column = alignment[:, col_idx].upper()
        filtered = [b for b in column if b in valid_bases] or ['N']  # 处理全无效列
        
        counts = Counter(filtered)
        total = len(filtered)
        
        candidates = []
        current_sum = 0.0
        for base, count in counts.most_common():
            freq = count / total
            if current_sum + freq > threshold:
                if current_sum == 0:  # 单个碱基超过阈值
                    candidates = [base]
                else:  # 多个碱基组合超过阈值
                    candidates.append(base)
                break
            candidates.append(base)
            current_sum += freq

        if ambiguity:
            key = frozenset(c.upper() for c in candidates)
            consensus.append(iupac_code.get(key, 'N'))
        else:
            consensus.append(candidates[0] if candidates else 'N')

    return ''.join(consensus)

def main():
    parser = argparse.ArgumentParser(
        description='Generate consensus sequence and calculate distances',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-i', '--input', required=True, help='Input aligned FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output consensus FASTA file')
    parser.add_argument('-d', '--distances', required=True, help='Output TSV file for distances')
    parser.add_argument('-t', '--threshold', type=float, default=0.5,
                        help='Consensus frequency threshold')
    parser.add_argument('--include-gap', action='store_true',
                        help='Include gaps in consensus calculation')
    parser.add_argument('--ambiguity', action='store_true',
                        help='Use IUPAC ambiguity codes')
    parser.add_argument('--name', help='Custom consensus name')

    args = parser.parse_args()

    # 读取比对文件并生成共识序列
    try:
        alignment = AlignIO.read(args.input, "fasta")
    except Exception as e:
        raise SystemExit(f"Error reading alignment: {str(e)}")

    consensus_name = args.name if args.name else f"Consensus_{os.path.basename(args.input)}"
    consensus_seq = calculate_consensus(
        alignment,
        threshold=args.threshold,
        include_gap=args.include_gap,
        ambiguity=args.ambiguity
    )

    # 写入共识序列

    consensus_clean = consensus_seq.replace('-', '').upper()

    with open(args.output, 'w') as f:
        f.write(f">{consensus_name}\n")
        for i in range(0, len(consensus_clean), 80):
            f.write(consensus_clean[i:i+80] + "\n")

    # 计算并保存距离数据
    with open(args.distances, 'w') as f_dist:
        f_dist.write("SequenceID\tEffectiveLength\tDistance\tNormalizedDistance\n")
        
        for record in alignment:
            # 去除gap并转为大写
            seq_clean = str(record.seq).upper().replace('-', '')
            
            # 计算实际参数
            effective_length = len(seq_clean)
            distance = levenshtein_distance(consensus_clean, seq_clean)
            normalized = distance / effective_length if effective_length > 0 else 0
            
            f_dist.write(f"{record.id}\t{effective_length}\t{distance}\t{normalized:.4f}\n")

if __name__ == "__main__":
    main()
```

</details>

<details open>
<summary><code>run_unit_mafft_consensus.sh</code> · Shell · 111 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/05_centromere/03_consensus/run_unit_mafft_consensus.sh)

```bash
#!/usr/bin/env bash
set -euo pipefail

# ========= 需要你改的 =========
GIFU_FA="Lotus_GifuT2T_v1.0.fasta"
MG20_FA="Lotus_MG20T2T_v1.0.fasta"

UNITS=(unit_161 unit_172 unit_186 unit_330)
GIFU_BED_PREFIX="gifu_"
MG20_BED_PREFIX="mg20_"

THREADS=32
OUTDIR="separate_then_merge_consensus"

CONS_PY="get_consensus.py"   # 你给的共识脚本文件名
THRESH=0.5
USE_IUPAC=""   # 不用IUPAC就置空：USE_IUPAC=""
INCLUDE_GAP="--include-gap"            # 要算gap就：INCLUDE_GAP="--include-gap"
# =================================

mkdir -p "$OUTDIR"/{01_seqs,02_msa,03_consensus,03_dist,04_merge2cons_msa,05_final}

command -v bedtools >/dev/null 2>&1 || { echo "[ERROR] bedtools not found"; exit 1; }
command -v mafft    >/dev/null 2>&1 || { echo "[ERROR] mafft not found"; exit 1; }
[[ -f "$CONS_PY" ]] || { echo "[ERROR] missing $CONS_PY"; exit 1; }

extract_and_consensus () {
  local genome="$1"
  local bed="$2"
  local tag="$3"     # Gifu or MG20
  local unit="$4"

  local seqfa="$OUTDIR/01_seqs/${tag}_${unit}.fa"
  local alnfa="$OUTDIR/02_msa/${tag}_${unit}.aln.fa"
  local consfa="$OUTDIR/03_consensus/${tag}_${unit}.cons.fa"
  local disttsv="$OUTDIR/03_dist/${tag}_${unit}.dist.tsv"

  if [[ ! -s "$bed" ]]; then
    echo "[WARN] $tag $unit bed missing/empty: $bed"
    return 1
  fi

  # 提取序列（按strand统一方向）
  bedtools getfasta -fi "$genome" -bed "$bed" -fo "$seqfa" -name -s

  # header 加前缀避免冲突（可选但建议）
  awk -v p="${tag}|" '/^>/{sub(/^>/,">"p);} {print}' "$seqfa" > "$seqfa.tmp" && mv "$seqfa.tmp" "$seqfa"

  local nseq
  nseq=$(grep -c "^>" "$seqfa" || true)
  if [[ "$nseq" -lt 2 ]]; then
    echo "[WARN] $tag $unit sequences <2 ($nseq), skip"
    return 1
  fi

  # MAFFT
  mafft --thread "$THREADS" --auto "$seqfa" > "$alnfa"

  # 共识 + 距离
  python "$CONS_PY" \
    -i "$alnfa" \
    -o "$consfa" \
    -d "$disttsv" \
    -t "$THRESH" \
    $USE_IUPAC \
    $INCLUDE_GAP \
    --name "${tag}_${unit}_consensus"

  echo "[OK] $tag $unit done: nseq=$nseq"
  return 0
}

for u in "${UNITS[@]}"; do
  echo "==== $u ===="

  gbed="${GIFU_BED_PREFIX}${u}.bed"
  mbed="${MG20_BED_PREFIX}${u}.bed"

  g_ok=0
  m_ok=0
  extract_and_consensus "$GIFU_FA" "$gbed" "Gifu" "$u" && g_ok=1 || true
  extract_and_consensus "$MG20_FA" "$mbed" "MG20" "$u" && m_ok=1 || true

  # 两个基因组都成功时：把两个共识再对齐，得到跨基因组比较/最终合并共识
  if [[ "$g_ok" -eq 1 && "$m_ok" -eq 1 ]]; then
    two_cons="$OUTDIR/04_merge2cons_msa/${u}.Gifu_MG20.cons.fa"
    cat "$OUTDIR/03_consensus/Gifu_${u}.cons.fa" \
        "$OUTDIR/03_consensus/MG20_${u}.cons.fa" > "$two_cons"

    two_aln="$OUTDIR/04_merge2cons_msa/${u}.Gifu_MG20.cons.aln.fa"
    mafft --thread "$THREADS" --auto "$two_cons" > "$two_aln"

    # 可选：再生成“最终跨基因组共识”（其实就是两条共识的共识）
    final_cons="$OUTDIR/05_final/${u}.final.cons.fa"
    final_dist="$OUTDIR/05_final/${u}.final.dist.tsv"
    python "$CONS_PY" \
      -i "$two_aln" \
      -o "$final_cons" \
      -d "$final_dist" \
      -t 0.5 \
      $USE_IUPAC \
      $INCLUDE_GAP \
      --name "${u}_Final_GifuMG20"

    echo "[OK] $u merged consensus done"
  else
    echo "[WARN] $u cannot merge 2-consensus (need both genomes succeed)"
  fi
done

echo "[DONE] outputs in: $OUTDIR"
```

</details>
