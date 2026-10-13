---
title: '易位核查'
sidebar:
  label: '易位核查'
  order: 3
---


## 中文

**作用**：用 read-pair 证据验证区分 Gifu 与 MG20 的 Chr1↔Chr2 易位是真实结构差异而非组装错误，并对多个样本批量分型：同一套 reads 比对到两个基因组，统计跨染色体 pair。
**入口**：`run.sh`
**输入**：`trans.conf`（3 列制表符分隔、无表头：`sample`、`bam_gifu`、`bam_mg`）及其中列出的 BAM；逐对计数另需坐标排序并建索引的 `gifu2gifu.bam`、`gifu2mg.bam`、`mg2gifu.bam`、`mg2mg.bam`；区域 `Chr1:22829936-22849936`（Gifu）与 `Chr2:99260160-99280160`（MG20）。
**输出**：`bam.list`；每样本 `{sample}.trans.gifu_anchored.merged_pairs.tsv`、`{sample}.trans.mg_anchored.merged_pairs.tsv` 与 4 个子集 BAM；汇总 `lotus_translocation_typing.tsv`（`sample`、`n_cross_in_gifu_bam`、`n_cross_in_mg_bam`、`call`、`ratio_MG20`、`note`）；逐对计数在 stdout 出计数表并写 `<bam>_chr1_chr2_pairs.txt` 等文件。
**运行**：`bash run.sh`（须在存放 `trans.conf`、BAM 与 TSV 的目录中运行；文件末尾的逐对计数命令以注释保留，需单独运行）。
**工具**：pysam（BAM 访问）、pandas、numpy；`samtools index`；Python 3（finder 经 `subprocess` 调用，`--min-support` 默认 5，两侧都达标为 `Mixed`，都低于为 `Uncertain`）；CSUB（`-q c01`、`-n 88`）。

<details>
<summary>English</summary>

**Purpose**: Read-pair evidence that the Chr1↔Chr2 translocation distinguishing Gifu from MG20 is a real structural difference rather than a mis-assembly, plus batch typing across samples: the same reads are mapped to both genomes and cross-chromosome pairs are counted.
**Entry point**: `run.sh`
**Inputs**: `trans.conf` (3 tab-separated columns, no header: `sample`, `bam_gifu`, `bam_mg`) and the BAMs it lists; per-pair counting also needs the coordinate-sorted, indexed `gifu2gifu.bam`, `gifu2mg.bam`, `mg2gifu.bam`, `mg2mg.bam`; regions `Chr1:22829936-22849936` (Gifu) and `Chr2:99260160-99280160` (MG20).
**Outputs**: `bam.list`; per sample `{sample}.trans.gifu_anchored.merged_pairs.tsv`, `{sample}.trans.mg_anchored.merged_pairs.tsv` and four subset BAMs; the summary `lotus_translocation_typing.tsv` (`sample`, `n_cross_in_gifu_bam`, `n_cross_in_mg_bam`, `call`, `ratio_MG20`, `note`); per-pair counting prints count tables to stdout and writes files such as `<bam>_chr1_chr2_pairs.txt`.
**Run**: `bash run.sh` (must run in the directory holding `trans.conf`, the BAMs and the TSVs; the per-pair counting commands are kept commented at the end of the file and must be run separately).
**Tools**: `pysam` (BAM access), pandas, numpy; `samtools index`; Python 3 (the finder is called through `subprocess`; `--min-support` defaults to 5, both sides at or above it give `Mixed`, below it on both sides gives `Uncertain`); CSUB (`-q c01`, `-n 88`).

</details>

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>batch_translocation_typing.py</code> · Python · 168 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/06_genome_comparison/03_translocation_check/batch_translocation_typing.py)

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
批量调用 find_cross_mapped_reads.py，对多个材料做“Gifu 型 / MG20 型”结构分型。

输入：
  一个 3 列的配置文件，比如 trans.conf.txt：
      sample   bam_gifu   bam_mg
  或你当前用的 lotus_translocation_typing.conf 那种格式。

含义：
  - bam_gifu : 该样本 reads 比对到 Gifu 基因组的 BAM
  - bam_mg   : 该样本 reads 比对到 MG20 基因组的 BAM

对于每个样本：
  调用 find_cross_mapped_reads.py 生成：
    {sample}.trans.gifu_anchored.merged_pairs.tsv
    {sample}.trans.mg_anchored.merged_pairs.tsv

  统计：
    n_gifu = gifu_anchored.merged_pairs.tsv 行数（在 Gifu 参考上 Chr1↔Chr2 跨染色体 pairs）
    n_mg   = mg_anchored.merged_pairs.tsv 行数（在 MG20 参考上 Chr1↔Chr2 跨染色体 pairs）

  判定逻辑（按你的要求）：
    - 如果 n_mg ≥ min_support 且 n_gifu < min_support  → Gifu-like
    - 如果 n_gifu ≥ min_support 且 n_mg < min_support  → MG20-like
    - 如果 n_gifu ≥ min_support 且 n_mg ≥ min_support → Mixed
    - 否则                                         → Uncertain

  默认 min_support = 5 （≥5 条就可以判断）
"""

import argparse
import subprocess
import os
import pandas as pd


def run_find_cross(bam_gifu, bam_mg, region_gifu, region_mg,
                   out_prefix, script_path):
    """调用 find_cross_mapped_reads.py"""
    cmd = [
        "python", script_path,
        "--bam-gifu", bam_gifu,
        "--bam-mg", bam_mg,
        "--region-gifu", region_gifu,
        "--region-mg", region_mg,
        "--out-prefix", out_prefix,
    ]
    print("[CMD]", " ".join(cmd))
    subprocess.run(cmd, check=True)


def count_rows_if_exists(path):
    """存在就读 TSV 数行，不存在返回 0"""
    if os.path.exists(path):
        df = pd.read_csv(path, sep="\t")
        return df.shape[0]
    else:
        return 0


def classify_sample(n_gifu, n_mg, min_support=5):
    """
    根据 GifuBAM / MG20BAM 中跨染色体 pair 数量给出分型（方向已经按你的要求反过来）：

      - n_mg >= min_support 且 n_gifu < min_support  → Gifu-like
      - n_gifu >= min_support 且 n_mg < min_support → MG20-like
      - 两边都 >= min_support                       → Mixed
      - 否则                                         → Uncertain
    """
    total = n_gifu + n_mg
    ratio_mg = (n_mg / total) if total > 0 else 0.0

    # 只在 MG20 参考中有跨染色体 pair：说明样本本身结构更像 Gifu
    if n_mg >= min_support and n_gifu < min_support:
        return "Gifu-like", ratio_mg, "cross only in MG20-BAM"

    # 只在 Gifu 参考中有跨染色体 pair：说明样本本身结构更像 MG20
    if n_gifu >= min_support and n_mg < min_support:
        return "MG20-like", ratio_mg, "cross only in Gifu-BAM"

    # 两边都有足够的跨染色体 pair：可能是混合群体 / 杂合结构
    if n_gifu >= min_support and n_mg >= min_support:
        return "Mixed", ratio_mg, "cross in both BAMs"

    # 支持数太少：不确定
    return "Uncertain", ratio_mg, f"support< {min_support} on both sides"


def main():
    ap = argparse.ArgumentParser(
        description="Batch typing of Chr1↔Chr2 translocation (Gifu vs MG20)"
    )
    ap.add_argument("--conf", required=True,
                    help="3 列: sample  bam_gifu  bam_mg")
    ap.add_argument("--region-gifu", required=True,
                    help="Gifu 断点区域，如 Chr1:22829936-22839936")
    ap.add_argument("--region-mg", required=True,
                    help="MG20 断点区域，如 Chr2:99270160-99280160")
    ap.add_argument("--find-script", default="find_cross_mapped_reads.py",
                    help="find_cross_mapped_reads.py 的路径")
    ap.add_argument("--out-summary", default="translocation_typing.summary.tsv",
                    help="输出 summary 表")
    ap.add_argument("--min-support", type=int, default=5,
                    help="最少跨染色体 pairs 支持数（默认 5）")
    args = ap.parse_args()

    # 读取配置文件
    conf = pd.read_csv(args.conf, sep="\t", header=None,
                       names=["sample", "bam_gifu", "bam_mg"])

    records = []

    for _, row in conf.iterrows():
        sample = row["sample"]
        bam_gifu = row["bam_gifu"]
        bam_mg = row["bam_mg"]

        out_prefix = f"{sample}.trans"
        print(f"\n===== Processing sample: {sample} =====")
        print(f"  bam_gifu: {bam_gifu}")
        print(f"  bam_mg  : {bam_mg}")

        # 1. 调用 find_cross_mapped_reads.py
        run_find_cross(
            bam_gifu=bam_gifu,
            bam_mg=bam_mg,
            region_gifu=args.region_gifu,
            region_mg=args.region_mg,
            out_prefix=out_prefix,
            script_path=args.find_script,
        )

        # 2. 统计 gifu-anchored / mg-anchored cross pairs
        gifu_tsv = f"{out_prefix}.gifu_anchored.merged_pairs.tsv"
        mg_tsv = f"{out_prefix}.mg_anchored.merged_pairs.tsv"

        n_gifu = count_rows_if_exists(gifu_tsv)
        n_mg = count_rows_if_exists(mg_tsv)

        call, ratio_mg, note = classify_sample(
            n_gifu=n_gifu,
            n_mg=n_mg,
            min_support=args.min_support,
        )

        print(f"  n_cross_in_gifu_bam = {n_gifu}")
        print(f"  n_cross_in_mg_bam   = {n_mg}")
        print(f"  => call = {call}  (ratio_MG20 = {ratio_mg:.3f})  note: {note}")

        records.append({
            "sample": sample,
            "n_cross_in_gifu_bam": n_gifu,
            "n_cross_in_mg_bam": n_mg,
            "call": call,
            "ratio_MG20": f"{ratio_mg:.4f}",
            "note": note,
        })

    summary_df = pd.DataFrame(records)
    summary_df.to_csv(args.out_summary, sep="\t", index=False)
    print(f"\n[Done] Summary written to: {args.out_summary}")


if __name__ == "__main__":
    main()
```

</details>

<details open>
<summary><code>bp_support_summary.py</code> · Python · 175 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/06_genome_comparison/03_translocation_check/bp_support_summary.py)

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pysam

def parse_bp(s):
    """
    解析形如: Chr1:22839936:-  或 Chr2:90335659:+
    返回: (chrom, pos(int), side)
    side 为 '+' 或 '-'
    """
    parts = s.split(":")
    if len(parts) != 3:
        raise ValueError("断点格式必须为 Chr:Pos:+/- ，例如 Chr1:22839936:-")
    chrom = parts[0]
    pos = int(parts[1])
    side = parts[2]
    if side not in ["+", "-"]:
        raise ValueError("方向必须是 + 或 -")
    return chrom, pos, side

def bp_to_region(chrom, pos, side, win):
    """
    根据断点和方向生成 1-based 闭区间 [start, end]
    side == '-'  -> 左侧窗口: [pos - (win-1), pos]
    side == '+'  -> 右侧窗口: [pos, pos + (win-1)]
    """
    if side == "-":
        start = max(1, pos - (win - 1))
        end = pos
    else:  # '+'
        start = pos
        end = pos + (win - 1)
    return chrom, start, end

def in_region(chrom, pos1based, target_chr, start, end):
    return (chrom == target_chr) and (start <= pos1based <= end)

def count_pairs_between_regions(bam, chr1, s1, e1, chr2, s2, e2):
    """
    统计：一端在 region1，一端在 region2 的成对 reads 数（按 read name 去重）
    region1 = chr1:s1-e1, region2 = chr2:s2-e2
    """
    pair_names = set()

    def scan_region(this_chr, this_start, this_end,
                    other_chr, other_start, other_end):
        start0 = this_start - 1  # 1-based -> 0-based
        end0 = this_end          # pysam.fetch: [start0, end0)
        for read in bam.fetch(this_chr, start0, end0):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if not read.is_paired or read.mate_is_unmapped:
                continue

            pos1based = read.reference_start + 1
            if not (this_start <= pos1based <= this_end):
                continue

            mate_id = read.next_reference_id
            if mate_id < 0:
                continue
            mate_chr = bam.get_reference_name(mate_id)
            mate_pos = read.next_reference_start + 1

            if in_region(mate_chr, mate_pos, other_chr, other_start, other_end):
                pair_names.add(read.query_name)

    scan_region(chr1, s1, e1, chr2, s2, e2)
    scan_region(chr2, s2, e2, chr1, s1, e1)

    return len(pair_names)

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "同时统计：一套 reads 在一个 BAM 中对真结构/假结构的支持 "
            "(true & false breakpoint pairs)."
        )
    )
    p.add_argument("--bam", required=True,
                   help="输入 BAM（排序并建立索引）")

    # 真结构的两个断点
    p.add_argument("--true-bp1", required=True,
                   help="真结构断点1: Chr:Pos:+/-  例如 Chr1:22839936:-")
    p.add_argument("--true-bp2", required=True,
                   help="真结构断点2: Chr:Pos:+/-")

    # 假结构的两个断点
    p.add_argument("--false-bp1", required=True,
                   help="假结构断点1: Chr:Pos:+/-")
    p.add_argument("--false-bp2", required=True,
                   help="假结构断点2: Chr:Pos:+/-")

    p.add_argument("--win", type=int, default=10000,
                   help="窗口大小（bp），默认 10000 = 10 kb")
    p.add_argument("--no-header", action="store_true",
                   help="不输出表头（用于 >> 追加到同一个 summary）")
    p.add_argument("--out", default="-",
                   help="输出文件名，默认 - 为 stdout")
    return p.parse_args()

def main():
    args = parse_args()
    bam = pysam.AlignmentFile(args.bam, "rb")

    # 解析真结构断点
    t_chr1_bp, t_pos1_bp, t_side1 = parse_bp(args.true_bp1)
    t_chr2_bp, t_pos2_bp, t_side2 = parse_bp(args.true_bp2)
    t_chr1, t_s1, t_e1 = bp_to_region(t_chr1_bp, t_pos1_bp, t_side1, args.win)
    t_chr2, t_s2, t_e2 = bp_to_region(t_chr2_bp, t_pos2_bp, t_side2, args.win)

    # 解析假结构断点
    f_chr1_bp, f_pos1_bp, f_side1 = parse_bp(args.false_bp1)
    f_chr2_bp, f_pos2_bp, f_side2 = parse_bp(args.false_bp2)
    f_chr1, f_s1, f_e1 = bp_to_region(f_chr1_bp, f_pos1_bp, f_side1, args.win)
    f_chr2, f_s2, f_e2 = bp_to_region(f_chr2_bp, f_pos2_bp, f_side2, args.win)

    # 统计真结构 / 假结构 read-pairs
    true_pairs = count_pairs_between_regions(
        bam, t_chr1, t_s1, t_e1, t_chr2, t_s2, t_e2
    )
    false_pairs = count_pairs_between_regions(
        bam, f_chr1, f_s1, f_e1, f_chr2, f_s2, f_e2
    )

    denom = true_pairs + false_pairs
    if denom > 0:
        true_ratio = true_pairs / denom
        false_ratio = false_pairs / denom
    else:
        true_ratio = float("nan")
        false_ratio = float("nan")

    header = "\t".join([
        "bam",
        "win_bp",
        "true_bp1", "true_bp2",
        "true_region1", "true_region2", "true_pairs",
        "false_bp1", "false_bp2",
        "false_region1", "false_region2", "false_pairs",
        "true_ratio", "false_ratio"
    ]) + "\n"

    true_region1 = f"{t_chr1}:{t_s1}-{t_e1}"
    true_region2 = f"{t_chr2}:{t_s2}-{t_e2}"
    false_region1 = f"{f_chr1}:{f_s1}-{f_e1}"
    false_region2 = f"{f_chr2}:{f_s2}-{f_e2}"

    line = "\t".join([
        args.bam,
        str(args.win),
        args.true_bp1, args.true_bp2,
        true_region1, true_region2, str(true_pairs),
        args.false_bp1, args.false_bp2,
        false_region1, false_region2, str(false_pairs),
        f"{true_ratio:.4f}", f"{false_ratio:.4f}",
    ]) + "\n"

    if args.out == "-" or args.out == "stdout":
        if not args.no_header:
            print(header, end="")
        print(line, end="")
    else:
        with open(args.out, "a") as f:
            if not args.no_header:
                f.write(header)
            f.write(line)

    bam.close()

if __name__ == "__main__":
    main()
```

</details>

<details open>
<summary><code>count_pairs_chr_chr2.py</code> · Python · 67 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/06_genome_comparison/03_translocation_check/count_pairs_chr_chr2.py)

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pysam

def parse_args():
    p = argparse.ArgumentParser(
        description="统计 BAM 中一对 reads 分别位于 chr1 和 chr2 的数量"
    )
    p.add_argument("--bam", required=True,
                   help="输入 BAM 文件（已排序并建立索引）")
    p.add_argument("--chr1", required=True,
                   help="染色体1名称，例如 Chr1")
    p.add_argument("--chr2", required=True,
                   help="染色体2名称，例如 Chr2")
    p.add_argument("--proper_only", action="store_true",
                   help="只统计 proper pair（read.is_proper_pair 为 True）")
    return p.parse_args()

def main():
    args = parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")

    chr1_name = args.chr1
    chr2_name = args.chr2

    total_pairs = 0
    pairs_chr1_chr2 = 0
    pairs_chr1_chr2_dir1 = 0  # read1 在 chr1，read2 在 chr2
    pairs_chr1_chr2_dir2 = 0  # read1 在 chr2，read2 在 chr1

    for read in bam.fetch(until_eof=True):
        # 只统计成对 reads 的 read1，避免重复计数
        if not read.is_paired or not read.is_read1:
            continue
        if read.is_unmapped or read.mate_is_unmapped:
            continue
        if args.proper_only and (not read.is_proper_pair):
            continue

        rchr = bam.get_reference_name(read.reference_id)
        mchr = bam.get_reference_name(read.next_reference_id)

        total_pairs += 1

        # 一端在 chr1，一端在 chr2
        if ((rchr == chr1_name and mchr == chr2_name) or
            (rchr == chr2_name and mchr == chr1_name)):
            pairs_chr1_chr2 += 1
            if rchr == chr1_name and mchr == chr2_name:
                pairs_chr1_chr2_dir1 += 1
            else:
                pairs_chr1_chr2_dir2 += 1

    bam.close()

    print(f"# BAM: {args.bam}")
    print(f"# chr1: {chr1_name}, chr2: {chr2_name}")
    print(f"total_read_pairs\t{total_pairs}")
    print(f"pairs_chr1_chr2_total\t{pairs_chr1_chr2}")
    print(f"pairs_chr1(chr1)-chr2(chr2)\t{pairs_chr1_chr2_dir1}")
    print(f"pairs_chr2(chr2)-chr1(chr1)\t{pairs_chr1_chr2_dir2}")

if __name__ == "__main__":
    main()
```

</details>

<details open>
<summary><code>count_pairs_region1_chr1tail_chr2.py</code> · Python · 127 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/06_genome_comparison/03_translocation_check/count_pairs_region1_chr1tail_chr2.py)

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pysam


def in_region(chrom, pos1based, target_chr, start, end):
    """判断 (chrom, pos) 是否落在 target_chr:start-end (1-based, 闭区间)"""
    return (chrom == target_chr) and (start <= pos1based <= end)


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "在参考基因组坐标下统计：一端在 Chr1[1-bp]，"
            "另一端在 Chr1[bp+1-end] 或 Chr2 任意位置的成对 reads 数量"
        )
    )
    p.add_argument("--bam", required=True,
                   help="比对到 Gifu 基因组的 BAM（已 sort + index）")
    p.add_argument("--chr1", default="Chr1",
                   help="Chr1 名称（默认 Chr1）")
    p.add_argument("--chr2", default="Chr2",
                   help="Chr2 名称（默认 Chr2）")
    p.add_argument("--bp", type=int, required=True,
                   help="Chr1 上的断点坐标（例如 22839936）")
    p.add_argument("--min-mapq", type=int, default=0,
                   help="最小比对质量（MAPQ），默认 0 不过滤")
    return p.parse_args()


def main():
    args = parse_args()
    bam = pysam.AlignmentFile(args.bam, "rb")

    chr1 = args.chr1
    chr2 = args.chr2
    bp   = args.bp

    # 1) 检查染色体名是否存在
    try:
        chr1_len = bam.get_reference_length(chr1)
    except ValueError:
        raise SystemExit(f"[ERROR] 在 BAM header 中找不到 {chr1}，"
                         f"实际染色体名可能是：{bam.references[:10]} ...")
    try:
        chr2_len = bam.get_reference_length(chr2)
    except ValueError:
        raise SystemExit(f"[ERROR] 在 BAM header 中找不到 {chr2}，"
                         f"实际染色体名可能是：{bam.references[:10]} ...")

    # 2) 定义区域（1-based 闭区间）
    region1_chr   = chr1
    region1_start = 1
    region1_end   = bp

    region2a_chr   = chr1
    region2a_start = bp + 1
    region2a_end   = chr1_len

    region2b_chr   = chr2
    region2b_start = 1
    region2b_end   = chr2_len

    total_pairs = 0
    n_1_to_1tail = 0   # 区域1 <-> Chr1 tail
    n_1_to_chr2  = 0   # 区域1 <-> Chr2
    n_read1_in_region1 = 0  # 有多少 read1 落在 region1（自检用）

    # 只看 read1，避免重复计数
    for read in bam.fetch(until_eof=True):
        if not read.is_paired or not read.is_read1:
            continue
        if read.is_unmapped or read.mate_is_unmapped:
            continue
        if read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < args.min_mapq:
            continue

        rchr = bam.get_reference_name(read.reference_id)
        rpos = read.reference_start + 1  # 1-based

        mate_id = read.next_reference_id
        if mate_id < 0:
            continue
        mchr = bam.get_reference_name(mate_id)
        mpos = read.next_reference_start + 1  # 1-based

        total_pairs += 1

        r_in_1     = in_region(rchr, rpos, region1_chr,   region1_start,  region1_end)
        r_in_1tail = in_region(rchr, rpos, region2a_chr,  region2a_start, region2a_end)
        r_in_chr2  = in_region(rchr, rpos, region2b_chr,  region2b_start, region2b_end)

        m_in_1     = in_region(mchr, mpos, region1_chr,   region1_start,  region1_end)
        m_in_1tail = in_region(mchr, mpos, region2a_chr,  region2a_start, region2a_end)
        m_in_chr2  = in_region(mchr, mpos, region2b_chr,  region2b_start, region2b_end)

        if r_in_1:
            n_read1_in_region1 += 1

        # 一端在区域1，另一端在 Chr1 tail
        if (r_in_1 and m_in_1tail) or (r_in_1tail and m_in_1):
            n_1_to_1tail += 1

        # 一端在区域1，另一端在 Chr2
        if (r_in_1 and m_in_chr2) or (r_in_chr2 and m_in_1):
            n_1_to_chr2 += 1

    bam.close()

    print(f"# BAM: {args.bam}")
    print(f"# Chr1: {chr1} (len={chr1_len}), Chr2: {chr2} (len={chr2_len})")
    print(f"# Chr1 break at: {bp}")
    print(f"# region1  = {chr1}:{region1_start}-{region1_end}")
    print(f"# region2a = {chr1}:{region2a_start}-{region2a_end} (Chr1 tail)")
    print(f"# region2b = {chr2}:{region2b_start}-{region2b_end} (Chr2 full)")
    print(f"total_read_pairs\t{total_pairs}")
    print(f"read1_in_region1\t{n_read1_in_region1}")
    print(f"pairs_region1_to_chr1tail\t{n_1_to_1tail}")
    print(f"pairs_region1_to_chr2\t{n_1_to_chr2}")


if __name__ == "__main__":
    main()
```

</details>

<details>
<summary><code>find_cross_mapped_reads.py</code> · Python · 318 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/06_genome_comparison/03_translocation_check/find_cross_mapped_reads.py)

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pysam
from collections import Counter


# ----------------- 基础工具函数 -----------------

def parse_region(s):
    """
    解析形如 Chr1:1-22839936 的字符串
    返回 (chrom, start, end) （1-based, 闭区间）
    """
    chrom, coords = s.split(":")
    start, end = coords.split("-")
    return chrom, int(start), int(end)


def in_region(chrom, pos1based, target_chr, start, end):
    """判断 (chrom, pos) 是否落在指定区域内"""
    return (chrom == target_chr) and (start <= pos1based <= end)


def collect_read1_names_in_region(bam_path, region_str):
    """
    在给定 BAM + 区域里收集所有落在区域内的 read1 名称集合
    只看 primary、成对、两端都 mapped
    """
    chrom, start, end = parse_region(region_str)
    names = set()

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        start0 = max(start - 1, 0)
        end0 = end
        for read in bam.fetch(chrom, start0, end0):
            if not read.is_paired or not read.is_read1:
                continue
            if read.is_unmapped or read.mate_is_unmapped:
                continue
            if read.is_secondary or read.is_supplementary:
                continue

            pos1based = read.reference_start + 1
            if in_region(chrom, pos1based, chrom, start, end):
                names.add(read.query_name)

    return names


def summarize_mate_positions(bam_path, region_str, read_names):
    """
    对于给定 BAM 和 区域 + read_name 集合：
    - 找到所有 read1 落在该区域且 read_name 属于 read_names
    - 统计 mate 比对到的染色体数量
    - 返回:
        mate_chr_counter: Counter(chr -> count)
        mate_records: list[(qname, read_chr, read_pos, mate_chr, mate_pos)]
    """
    chrom, start, end = parse_region(region_str)
    mate_chr_counter = Counter()
    mate_records = []

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        start0 = max(start - 1, 0)
        end0 = end
        for read in bam.fetch(chrom, start0, end0):
            if not read.is_paired or not read.is_read1:
                continue
            if read.is_unmapped or read.mate_is_unmapped:
                continue
            if read.is_secondary or read.is_supplementary:
                continue

            if read.query_name not in read_names:
                continue

            rchr = bam.get_reference_name(read.reference_id)
            rpos = read.reference_start + 1

            mate_id = read.next_reference_id
            if mate_id < 0:
                continue
            mchr = bam.get_reference_name(mate_id)
            mpos = read.next_reference_start + 1

            mate_chr_counter[mchr] += 1
            mate_records.append((read.query_name, rchr, rpos, mchr, mpos))

    return mate_chr_counter, mate_records


def write_subset_bam(in_bam_path, out_bam_path, read_name_set):
    """
    从原始 BAM 中抽取 query_name 属于 read_name_set 的所有 reads（两端都写出）
    """
    with pysam.AlignmentFile(in_bam_path, "rb") as bam_in:
        with pysam.AlignmentFile(out_bam_path, "wb", template=bam_in) as bam_out:
            for read in bam_in.fetch(until_eof=True):
                if read.query_name in read_name_set:
                    bam_out.write(read)


def fetch_positions_from_bam(bam_path, read_names):
    """
    在 BAM 中查找这些 read_names 的 read1 比对位置 + mate 位置
    返回字典: name -> (read_chr, read_pos, mate_chr, mate_pos)
    """
    name_set = set(read_names)
    result = {}

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.query_name not in name_set:
                continue
            if not read.is_paired or not read.is_read1:
                continue
            if read.is_unmapped or read.mate_is_unmapped:
                continue
            if read.is_secondary or read.is_supplementary:
                continue

            rchr = bam.get_reference_name(read.reference_id)
            rpos = read.reference_start + 1
            mate_id = read.next_reference_id
            if mate_id < 0:
                continue
            mchr = bam.get_reference_name(mate_id)
            mpos = read.next_reference_start + 1

            result[read.query_name] = (rchr, rpos, mchr, mpos)

            if len(result) == len(name_set):
                break

    return result


def pair_distance(chr_a, pos_a, chr_b, pos_b):
    """
    计算一对 reads 之间的距离：
    - 如果在同一条染色体: abs(pos_b - pos_a)
    - 如果跨染色体: 返回 -1 （可当 NA）
    """
    if chr_a == chr_b:
        return abs(pos_b - pos_a)
    else:
        return -1


# ----------------- 主程序 -----------------

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "同一套 reads 比对到两个基因组："
            "在 bam-gifu 中落在 region-gifu，且在 bam-mg 中落在 region-mg 的 read1 交集；\n"
            "在两边分别找出 Chr1<->Chr2 的跨染色体 pairs，"
            "输出 merged 表（两基因组上的位置 + pair 距离），"
            "并导出对应的小 BAM。"
        )
    )
    p.add_argument("--bam-gifu", required=True,
                   help="同一套 reads 比对到 Gifu 基因组的 BAM")
    p.add_argument("--bam-mg", required=True,
                   help="同一套 reads 比对到 MG20 基因组的 BAM")
    p.add_argument("--region-gifu", required=True,
                   help="Gifu 区域，如 Chr1:22829936-22839936")
    p.add_argument("--region-mg", required=True,
                   help="MG20 区域，如 Chr2:99270160-99280160")
    p.add_argument("--out-prefix", default="cross_reads",
                   help="输出文件前缀")
    return p.parse_args()


def main():
    args = parse_args()

    # 把 region 中的染色体名当作“Chr1”和“Chr2”的标签
    gifu_chr, _, _ = parse_region(args.region_gifu)   # 例如 Chr1
    mg_chr, _, _   = parse_region(args.region_mg)     # 例如 Chr2

    chr1_name = gifu_chr
    chr2_name = mg_chr

    # ---------- 1. 区域内 read1 名称 & 交集 ----------
    print("[*] 收集 Gifu BAM 区域内的 read1 名称...")
    names_gifu = collect_read1_names_in_region(args.bam_gifu, args.region_gifu)
    print(f"    Gifu 区域内 read1 数量: {len(names_gifu)}")

    print("[*] 收集 MG20 BAM 区域内的 read1 名称...")
    names_mg = collect_read1_names_in_region(args.bam_mg, args.region_mg)
    print(f"    MG20 区域内 read1 数量: {len(names_mg)}")

    cross_names = names_gifu & names_mg
    print(f"[*] 同时满足两个区域的 read1 数量（交集）: {len(cross_names)}")
    if not cross_names:
        print("[!] 交集为空，结束。")
        return

    # ---------- 2. 在两边 BAM 中统计 mate 位置 ----------
    print("[*] 分析这些 reads 在 Gifu BAM 中的 mate 位置...")
    gifu_mate_chr_counter, gifu_mate_records = summarize_mate_positions(
        args.bam_gifu, args.region_gifu, cross_names
    )

    print("[*] 分析这些 reads 在 MG20 BAM 中的 mate 位置...")
    mg_mate_chr_counter, mg_mate_records = summarize_mate_positions(
        args.bam_mg, args.region_mg, cross_names
    )

    print("\n=== Summary: mate 染色体分布（Gifu BAM） ===")
    for chr_name, count in gifu_mate_chr_counter.most_common():
        print(f"Gifu_BAM_mate_chr\t{chr_name}\t{count}")

    print("\n=== Summary: mate 染色体分布（MG20 BAM） ===")
    for chr_name, count in mg_mate_chr_counter.most_common():
        print(f"MG20_BAM_mate_chr\t{chr_name}\t{count}")

    # ---------- 3. 各自 BAM 中的 Chr1<->Chr2 跨染色体 pairs ----------
    gifu_cross_records = []
    gifu_cross_names = set()
    for qname, rchr, rpos, mchr, mpos in gifu_mate_records:
        if ((rchr == chr1_name and mchr == chr2_name) or
            (rchr == chr2_name and mchr == chr1_name)):
            gifu_cross_records.append((qname, rchr, rpos, mchr, mpos))
            gifu_cross_names.add(qname)

    mg_cross_records = []
    mg_cross_names = set()
    for qname, rchr, rpos, mchr, mpos in mg_mate_records:
        if ((rchr == chr1_name and mchr == chr2_name) or
            (rchr == chr2_name and mchr == chr1_name)):
            mg_cross_records.append((qname, rchr, rpos, mchr, mpos))
            mg_cross_names.add(qname)

    print("\n=== Chr1<->Chr2 跨染色体 pairs 数量 ===")
    print(f"Gifu BAM 中: {len(gifu_cross_records)} 条")
    print(f"MG20 BAM 中: {len(mg_cross_records)} 条")

    # ---------- 4A. 以 Gifu 为锚点的 merged 表 ----------
    if gifu_cross_records:
        print("[*] 生成 Gifu-anchored merged 表...")
        mg_pos_map_for_gifu = fetch_positions_from_bam(args.bam_mg, gifu_cross_names)

        merged_tsv_gifu = f"{args.out_prefix}.gifu_anchored.merged_pairs.tsv"
        with open(merged_tsv_gifu, "w") as out:
            out.write(
                "read_name\t"
                "gifu_read_chr\tgifu_read_pos\tgifu_mate_chr\tgifu_mate_pos\tgifu_pair_dist\t"
                "mg_read_chr\tmg_read_pos\tmg_mate_chr\tmg_mate_pos\tmg_pair_dist\n"
            )
            for qname, g_rchr, g_rpos, g_mchr, g_mpos in gifu_cross_records:
                g_dist = pair_distance(g_rchr, g_rpos, g_mchr, g_mpos)
                if qname in mg_pos_map_for_gifu:
                    m_rchr, m_rpos, m_mchr, m_mpos = mg_pos_map_for_gifu[qname]
                    m_dist = pair_distance(m_rchr, m_rpos, m_mchr, m_mpos)
                else:
                    m_rchr = m_mchr = "NA"
                    m_rpos = m_mpos = -1
                    m_dist = -1
                out.write(
                    f"{qname}\t"
                    f"{g_rchr}\t{g_rpos}\t{g_mchr}\t{g_mpos}\t{g_dist}\t"
                    f"{m_rchr}\t{m_rpos}\t{m_mchr}\t{m_mpos}\t{m_dist}\n"
                )
        print(f"    写入: {merged_tsv_gifu}")

        # 导出对应的小 BAM
        gifu_pairs_bam = f"{args.out_prefix}.gifu_chr1_chr2_pairs.bam"
        print(f"[*] 从 {args.bam_gifu} 中抽取 {len(gifu_cross_names)} 条 reads 写入 {gifu_pairs_bam}")
        write_subset_bam(args.bam_gifu, gifu_pairs_bam, gifu_cross_names)

        gifu_reads_in_mg_bam = f"{args.out_prefix}.gifu_reads_in_mg.bam"
        print(f"[*] 从 {args.bam_mg} 中抽取同名 reads 写入 {gifu_reads_in_mg_bam}")
        write_subset_bam(args.bam_mg, gifu_reads_in_mg_bam, gifu_cross_names)

    # ---------- 4B. 以 MG 为锚点的 merged 表（一般你更关注这一部分） ----------
    if mg_cross_records:
        print("[*] 生成 MG-anchored merged 表...")
        gifu_pos_map_for_mg = fetch_positions_from_bam(args.bam_gifu, mg_cross_names)

        merged_tsv_mg = f"{args.out_prefix}.mg_anchored.merged_pairs.tsv"
        with open(merged_tsv_mg, "w") as out:
            out.write(
                "read_name\t"
                "mg_read_chr\tmg_read_pos\tmg_mate_chr\tmg_mate_pos\tmg_pair_dist\t"
                "gifu_read_chr\tgifu_read_pos\tgifu_mate_chr\tgifu_mate_pos\tgifu_pair_dist\n"
            )
            for qname, m_rchr, m_rpos, m_mchr, m_mpos in mg_cross_records:
                m_dist = pair_distance(m_rchr, m_rpos, m_mchr, m_mpos)
                if qname in gifu_pos_map_for_mg:
                    g_rchr, g_rpos, g_mchr, g_mpos = gifu_pos_map_for_mg[qname]
                    g_dist = pair_distance(g_rchr, g_rpos, g_mchr, g_mpos)
                else:
                    g_rchr = g_mchr = "NA"
                    g_rpos = g_mpos = -1
                    g_dist = -1
                out.write(
                    f"{qname}\t"
                    f"{m_rchr}\t{m_rpos}\t{m_mchr}\t{m_mpos}\t{m_dist}\t"
                    f"{g_rchr}\t{g_rpos}\t{g_mchr}\t{g_mpos}\t{g_dist}\n"
                )
        print(f"    写入: {merged_tsv_mg}")

        # 导出对应的小 BAM
        mg_pairs_bam = f"{args.out_prefix}.mg_chr1_chr2_pairs.bam"
        print(f"[*] 从 {args.bam_mg} 中抽取 {len(mg_cross_names)} 条 reads 写入 {mg_pairs_bam}")
        write_subset_bam(args.bam_mg, mg_pairs_bam, mg_cross_names)

        mg_reads_in_gifu_bam = f"{args.out_prefix}.mg_reads_in_gifu.bam"
        print(f"[*] 从 {args.bam_gifu} 中抽取同名 reads 写入 {mg_reads_in_gifu_bam}")
        write_subset_bam(args.bam_gifu, mg_reads_in_gifu_bam, mg_cross_names)


if __name__ == "__main__":
    main()
```

</details>

<details open>
<summary><code>run.sh</code> · Shell · 102 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/06_genome_comparison/03_translocation_check/run.sh)

```bash
#!/bin/bash
#CSUB -J run
#CSUB -q c01
#CSUB -o run.out
#CSUB -e run.error
#CSUB -n 88
#CSUB -R span[hosts=1]

# Two entry points share this directory:
#   (a) batch structural typing of the Chr1<->Chr2 translocation across samples
#       (the `python batch_translocation_typing.py` block below);
#   (b) the per-pair counting commands used for the per-pair plots, recorded
#       at the end of the file. They were run interactively rather than from
#       this script, which is why they are commented out.

echo "[*] STEP 1: 从 trans.conf 提取 BAM 列表 ..."
cut -f2,3 trans.conf | tr '\t' '\n' | sort -u > bam.list

echo "[*] STEP 2: 自动检查/补齐 BAM 索引 ..."
while read bam; do
    if [ ! -f "$bam" ]; then
        echo "[WARN] 找不到 BAM: $bam"
        continue
    fi

    bai1="${bam}.bai"
    bai2="${bam%.bam}.bai"

    if [ -f "$bai1" ] || [ -f "$bai2" ]; then
        echo "[OK] 已存在索引: $bam"
    else
        echo "[INDEX] 开始建立索引: $bam"
        samtools index "$bam"
    fi
done < bam.list


python batch_translocation_typing.py \
  --conf trans.conf \
  --region-gifu Chr1:22829936-22849936 \
  --region-mg   Chr2:99260160-99280160 \
  --find-script find_cross_mapped_reads.py \
  --out-summary lotus_translocation_typing.tsv


# ── Per-pair counting for the per-pair plots (run interactively) ─────────
#samtools index gifu2gifu.bam
#samtools index gifu2mg.bam
#samtools index mg2gifu.bam
#samtools index mg2mg.bam

#python count_pairs_chr_chr2.py \
#  --bam gifu2gifu.bam \
#  --chr1 Chr1 \
#  --chr2 Chr2 \
#  > gifu2gifu_chr1_chr2_pairs.txt

#python count_pairs_chr_chr2.py \
#  --bam mg2gifu.bam \
#  --chr1 Chr1 \
#  --chr2 Chr2 \
#  > mg2gifu_chr1_chr2_pairs.txt

#python count_pairs_chr_chr2.py \
#  --bam gifu2mg.bam \
#  --chr1 Chr1 \
#  --chr2 Chr2 \
#  > gifu2mg_chr1_chr2_pairs.txt

#python count_pairs_chr_chr2.py \
#  --bam mg2mg.bam \
#  --chr1 Chr1 \
#  --chr2 Chr2 \
#  > mg2mg_chr1_chr2_pairs.txt

#python count_pairs_region1_chr1tail_chr2.py \
#  --bam gifu2gifu.bam \
#  --chr1 Chr1 \
#  --chr2 Chr2 \
#  --bp 22839936 \
#  > gifu2gifu_region_pairs.txt

#python count_pairs_region1_chr1tail_chr2.py \
#  --bam mg2gifu.bam \
#  --chr1 Chr1 \
#  --chr2 Chr2 \
#  --bp 22839936 \
#  > mg2gifu_region_pairs.txt

#python find_cross_mapped_reads.py \
#  --bam-gifu gifu2gifu.bam \
#  --bam-mg   gifu2mg.bam \
#  --region-gifu Chr1:22829936-22849936 \
#  --region-mg   Chr2:99260160-99280160 \
#  --out-prefix gifu_cross

#python find_cross_mapped_reads.py \
#  --bam-gifu mg2gifu.bam \
#  --bam-mg   mg2mg.bam \
#  --region-gifu Chr1:22829936-22849936 \
#  --region-mg   Chr2:99260160-99280160 \
#  --out-prefix mg_cross
```

</details>
