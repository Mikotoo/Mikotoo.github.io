---
title: 'rDNA 位点与阵列'
sidebar:
  label: 'rDNA 位点与阵列'
  order: 1
---


## 中文

**作用**：用 Barrnap 注释 Gifu T2T 的 rRNA 基因以确定 rDNA 位点，切出各阵列做 StainedGlass 自比对，并统计每个阵列的 18S/5.8S/28S/5S 特征数以得到完整 45S 重复单元数。
**入口**：`barrnap.sh`
**输入**：`Lotus_GifuT2T_v1.0.fasta`（含随附 `.fai`，Chr1–Chr6）；区段 BED `Chr2_0_3Mbp.bed`、`Chr2_24_26Mbp.bed`、`Chr5_14_16Mbp.bed`、`Chr6_8_10Mbp.bed`、`Chr6_25_28Mbp.bed`；Barrnap 的 `rRNA.gff`/`rRNA.fa`；生产 GFF `Gifu_rRNA.gff`、`MG20_rRNA.gff`。
**输出**：`rRNA.gff`/`rRNA.fa`；`result/*_output.tbl.gz`、`*_output.bed.gz` 及 `.fa`/`.fai`/`.bed`/`.fasta`/`.mmi`/`.bam` 中间文件；`results/` 下各阵列的 PDF/PNG 图；阵列组成表打印到标准输出并写入 `03_annotation/rDNA_array_summary.tsv`（45S 完整单元数 = 阵列内 `min(18S, 5.8S, 28S)`，阵列坐标硬编码在 `count_rDNA.py` 中）。
**运行**：`bash barrnap.sh` → `python3 count_rDNA.py` → `bash result/process.sh`（`process.sh` 从自身目录读 `../Lotus_GifuT2T_v1.0.fasta` 与区段 BED，结果写到当前目录，须在 `result/` 下运行）。
**工具**：Barrnap `--kingdom euk --threads 88`；minimap2 建索引 `-f 1000 -s 400 -ax ava-ont`、比对 `-t 88 -f 10000 -s 400 -ax ava-ont --dual=yes --eqx`；StainedGlass 0.6 `samIdentity.py --threads 88 --matches 400 --header`、`refmt.py --window 2000`、`aln_plot.R --threads 88`；bedtools `getfasta`、`makewindows -w 2000`；samtools `faidx`、`sort -m 4G`；`bgzip`、`Rscript`、`python3.13`（经 ``{SOFTWARE_LOCAL}``）；`count_rDNA.py` 仅用标准库。

<details>
<summary>English</summary>

**Purpose**: Annotate the Gifu T2T rRNA genes with Barrnap to define the rDNA loci, extract each array for StainedGlass self-identity analysis, and count the 18S/5.8S/28S/5S features per array to derive the number of complete 45S repeat units.
**Entry point**: `barrnap.sh`
**Inputs**: `Lotus_GifuT2T_v1.0.fasta` (with the shipped `.fai`, Chr1–Chr6); the region BEDs `Chr2_0_3Mbp.bed`, `Chr2_24_26Mbp.bed`, `Chr5_14_16Mbp.bed`, `Chr6_8_10Mbp.bed`, `Chr6_25_28Mbp.bed`; Barrnap's `rRNA.gff`/`rRNA.fa`; the production GFFs `Gifu_rRNA.gff` and `MG20_rRNA.gff`.
**Outputs**: `rRNA.gff`/`rRNA.fa`; `result/*_output.tbl.gz` and `*_output.bed.gz` plus the `.fa`/`.fai`/`.bed`/`.fasta`/`.mmi`/`.bam` intermediates; PDF/PNG plots per array under `results/`; the array composition table on stdout and in `03_annotation/rDNA_array_summary.tsv` (complete 45S units = `min(18S, 5.8S, 28S)` inside each array, with the array coordinates hardcoded in `count_rDNA.py`).
**Run**: `bash barrnap.sh` → `python3 count_rDNA.py` → `bash result/process.sh` (`process.sh` reads `../Lotus_GifuT2T_v1.0.fasta` and the region BEDs from its own directory and writes into the current directory, so it must run inside `result/`).
**Tools**: Barrnap `--kingdom euk --threads 88`; minimap2 index `-f 1000 -s 400 -ax ava-ont`, alignment `-t 88 -f 10000 -s 400 -ax ava-ont --dual=yes --eqx`; StainedGlass 0.6 `samIdentity.py --threads 88 --matches 400 --header`, `refmt.py --window 2000`, `aln_plot.R --threads 88`; bedtools `getfasta`, `makewindows -w 2000`; samtools `faidx`, `sort -m 4G`; `bgzip`, `Rscript`, `python3.13` (through ``{SOFTWARE_LOCAL}``); `count_rDNA.py` is standard library only.

</details>

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>barrnap.sh</code> · Shell · 3 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/04_rDNA/01_locus_and_arrays/barrnap.sh)

```bash
#!/bin/bash

barrnap --kingdom euk --threads 88 -o rRNA.fa Lotus_GifuT2T_v1.0.fasta > rRNA.gff
```

</details>

<details open>
<summary><code>count_rDNA.py</code> · Python · 186 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/04_rDNA/01_locus_and_arrays/count_rDNA.py)

```python
#!/usr/bin/env python3
"""
Parse barrnap rRNA GFF files to count 18S/5.8S/28S/5S features per rDNA array,
and compute the number of complete 45S repeat units.

Output: rDNA array composition table for Gifu and MG20.
"""

import sys
from collections import defaultdict
from pathlib import Path

# ── rDNA array definitions (1-based coordinates, matching GFF) ──────────────
# Format: label -> (chrom, start, end, expected_type)
# expected_type: "45S" or "5S"

GIFU_ARRAYS = [
    ("Array-1  (Chr2 45S)",   "Chr2", 537436,   949316,   "45S"),
    ("Array-2  (Chr2 5S)",    "Chr2", 24476203,  25027540, "5S"),
    ("Array-3  (Chr5 45S)",   "Chr5", 14943605,  15344676, "45S"),
    ("Array-4  (Chr6 45S)",   "Chr6", 8844951,   8943185,  "45S"),
    ("Array-5  (Chr6 45S)",   "Chr6", 26125599,  27223832, "45S"),
]

MG20_ARRAYS = [
    ("Array-1  (Chr2 45S)",   "Chr2", 1528120,   2857257,  "45S"),
    ("Array-2  (Chr2 5S)",    "Chr2", 26673364,  27897859, "5S"),
    ("Array-3  (Chr5 45S)",   "Chr5", 15013665,  15659889, "45S"),
    ("Array-4  (Chr6 45S)",   "Chr6", 15057397,  15451961, "45S"),
    ("Array-5  (Chr6 45S)",   "Chr6", 26255620,  26622326, "45S"),
]

# ── GFF paths ───────────────────────────────────────────────────────────────
GIFU_GFF = Path(
    "{PROJ_04_LOTUS_GENOME}/"
    "03_annotation/Gifu/rDNA/Gifu_rRNA.gff"
)
MG20_GFF = Path(
    "{PROJ_04_LOTUS_GENOME}/"
    "03_annotation/MG20/rDNA/MG20_rRNA.gff"
)


def parse_rrna_type(name_attr: str) -> str:
    """Extract rRNA subunit type from barrnap GFF Name attribute."""
    if "5S" in name_attr and "5_8S" not in name_attr and "5.8S" not in name_attr:
        return "5S"
    if "5_8S" in name_attr or "5.8S" in name_attr:
        return "5.8S"
    if "18S" in name_attr:
        return "18S"
    if "28S" in name_attr:
        return "28S"
    return "other"


def load_gff(gff_path: Path):
    """Load barrnap GFF, return list of (chrom, start, end, rrna_type, partial)."""
    features = []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            chrom, _, ftype, start, end, _, _, _, attr = parts
            if ftype != "rRNA":
                continue
            rrna_type = parse_rrna_type(attr)
            is_partial = "partial" in attr
            features.append((chrom, int(start), int(end), rrna_type, is_partial))
    return features


def count_in_region(features, chrom, region_start, region_end):
    """Count 18S/5.8S/28S/5S features with any overlap to the region."""
    counts = {"18S": 0, "5.8S": 0, "28S": 0, "5S": 0}
    partials = {"18S": 0, "5.8S": 0, "28S": 0, "5S": 0}
    for chrom_f, start, end, rrna_type, is_partial in features:
        if chrom_f != chrom:
            continue
        if start > region_end or end < region_start:
            continue
        if rrna_type in counts:
            counts[rrna_type] += 1
            if is_partial:
                partials[rrna_type] += 1
    return counts, partials


def fmt_complete(total, partial):
    """Format as 'total (N partial)' if partials exist, else just total."""
    if partial > 0:
        return f"{total} ({partial} partial)"
    return str(total)


def process_species(species_name, arrays, gff_path):
    """Print a table for one species."""
    features = load_gff(gff_path)
    print(f"## {species_name}")
    print()
    header = (
        f"{'Array':<28s} {'Chrom':>5s} {'Coord (bp)':>24s} "
        f"{'Length (bp)':>12s} {'18S':>8s} {'5.8S':>8s} "
        f"{'28S':>8s} {'5S':>8s} {'45S units':>10s}"
    )
    print(header)
    print("-" * len(header))

    for label, chrom, rs, re, expected in arrays:
        length = re - rs
        coord_str = f"{rs:,}-{re:,}"
        counts, partials = count_in_region(features, chrom, rs, re)

        if expected == "45S":
            # Complete 45S repeat units = min(18S, 5.8S, 28S)
            units_45s = min(counts["18S"], counts["5.8S"], counts["28S"])
            units_str = str(units_45s)
        else:
            units_str = "—"

        row = (
            f"{label:<28s} {chrom:>5s} {coord_str:>24s} "
            f"{length:>12,d} "
            f"{fmt_complete(counts['18S'], partials['18S']):>8s} "
            f"{fmt_complete(counts['5.8S'], partials['5.8S']):>8s} "
            f"{fmt_complete(counts['28S'], partials['28S']):>8s} "
            f"{fmt_complete(counts['5S'], partials['5S']):>8s} "
            f"{units_str:>10s}"
        )
        print(row)
    print()


def main():
    print("# rDNA array composition table for Lotus japonicus")
    print()
    process_species("Gifu", GIFU_ARRAYS, GIFU_GFF)
    process_species("MG20", MG20_ARRAYS, MG20_GFF)

    # ── TSV output for spreadsheet ──
    tsv_path = Path(
        "{PROJ_04_LOTUS_GENOME}/"
        "03_annotation/rDNA_array_summary.tsv"
    )
    with open(tsv_path, "w") as f:
        hdr = [
            "Species", "Array", "Chromosome", "Start", "End",
            "Length_bp", "18S_count", "18S_partial",
            "5.8S_count", "5.8S_partial",
            "28S_count", "28S_partial",
            "5S_count", "5S_partial",
            "Complete_45S_units",
        ]
        f.write("\t".join(hdr) + "\n")

        for species, arrays, gff in [
            ("Gifu", GIFU_ARRAYS, GIFU_GFF),
            ("MG20", MG20_ARRAYS, MG20_GFF),
        ]:
            feats = load_gff(gff)
            for label, chrom, rs, re, expected in arrays:
                counts, partials = count_in_region(feats, chrom, rs, re)
                if expected == "45S":
                    units = min(counts["18S"], counts["5.8S"], counts["28S"])
                else:
                    units = 0
                f.write("\t".join(map(str, [
                    species,
                    label.strip().split("(")[0].strip(),
                    chrom, rs, re,
                    re - rs,
                    counts["18S"], partials["18S"],
                    counts["5.8S"], partials["5.8S"],
                    counts["28S"], partials["28S"],
                    counts["5S"], partials["5S"],
                    units,
                ])) + "\n")

    print(f"TSV saved to: {tsv_path}")


if __name__ == "__main__":
    main()
```

</details>

<details open>
<summary><code>result/process.sh</code> · Shell · 78 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/04_rDNA/01_locus_and_arrays/result/process.sh)

```bash
#!/bin/bash

bedtools getfasta -fi ../Lotus_GifuT2T_v1.0.fasta -bed Chr2_0_3Mbp.bed > Chr2_0_3Mbp.fa
samtools faidx Chr2_0_3Mbp.fa
bedtools makewindows -g Chr2_0_3Mbp.fa.fai -w 2000 > Chr2_0_3Mbp_output.bed
bedtools getfasta -fi Chr2_0_3Mbp.fa -bed Chr2_0_3Mbp_output.bed > Chr2_0_3Mbp_output.2000.fasta

minimap2 -f 1000 -s 400 -ax ava-ont -d Chr2_0_3Mbp_output.fasta.mmi Chr2_0_3Mbp_output.2000.fasta
minimap2 -t 88 -f 10000 -s 400 -ax ava-ont --dual=yes --eqx Chr2_0_3Mbp_output.fasta.mmi Chr2_0_3Mbp_output.2000.fasta | samtools sort -m 4G -o Chr2_0_3Mbp_output.bam
python3.13 {SOFTWARE_LOCAL}/StainedGlass-0.6/workflow/scripts/samIdentity.py --threads 88 --matches 400 --header Chr2_0_3Mbp_output.bam > Chr2_0_3Mbp_output.tbl
bgzip -c Chr2_0_3Mbp_output.tbl > Chr2_0_3Mbp_output.tbl.gz
python3.13 {SOFTWARE_LOCAL}/StainedGlass-0.6/workflow/scripts/refmt.py --window 2000 --fai Chr2_0_3Mbp.fa.fai --full Chr2_0_3Mbp_output.full.tbl.gz Chr2_0_3Mbp_output.tbl.gz Chr2_0_3Mbp_output.bed.gz
mkdir -p results/Chr2_0_3Mbp_figures/pdfs
mkdir -p results/Chr2_0_3Mbp_figures/pngs
Rscript {SOFTWARE_LOCAL}/StainedGlass-0.6/workflow/scripts/aln_plot.R -b Chr2_0_3Mbp_output.bed.gz --threads 88 --prefix Chr2_0_3Mbp



bedtools getfasta -fi ../Lotus_GifuT2T_v1.0.fasta -bed Chr2_24_26Mbp.bed > Chr2_24_26Mbp.fa
samtools faidx Chr2_24_26Mbp.fa
bedtools makewindows -g Chr2_24_26Mbp.fa.fai -w 2000 > Chr2_24_26Mbp_output.bed
bedtools getfasta -fi Chr2_24_26Mbp.fa -bed Chr2_24_26Mbp_output.bed > Chr2_24_26Mbp_output.2000.fasta

minimap2 -f 1000 -s 400 -ax ava-ont -d Chr2_24_26Mbp_output.fasta.mmi Chr2_24_26Mbp_output.2000.fasta
minimap2 -t 88 -f 10000 -s 400 -ax ava-ont --dual=yes --eqx Chr2_24_26Mbp_output.fasta.mmi Chr2_24_26Mbp_output.2000.fasta | samtools sort -m 4G -o Chr2_24_26Mbp_output.bam
python3.13 {SOFTWARE_LOCAL}/StainedGlass-0.6/workflow/scripts/samIdentity.py --threads 88 --matches 400 --header Chr2_24_26Mbp_output.bam > Chr2_24_26Mbp_output.tbl
bgzip -c Chr2_24_26Mbp_output.tbl > Chr2_24_26Mbp_output.tbl.gz
python3.13 {SOFTWARE_LOCAL}/StainedGlass-0.6/workflow/scripts/refmt.py --window 2000 --fai Chr2_24_26Mbp.fa.fai --full Chr2_24_26Mbp_output.full.tbl.gz Chr2_24_26Mbp_output.tbl.gz Chr2_24_26Mbp_output.bed.gz
mkdir -p results/Chr2_24_26Mbp_figures/pdfs
mkdir -p results/Chr2_24_26Mbp_figures/pngs
Rscript {SOFTWARE_LOCAL}/StainedGlass-0.6/workflow/scripts/aln_plot.R -b Chr2_24_26Mbp_output.bed.gz --threads 88 --prefix Chr2_24_26Mbp



bedtools getfasta -fi ../Lotus_GifuT2T_v1.0.fasta -bed Chr5_14_16Mbp.bed > Chr5_14_16Mbp.fa
samtools faidx Chr5_14_16Mbp.fa
bedtools makewindows -g Chr5_14_16Mbp.fa.fai -w 2000 > Chr5_14_16Mbp_output.bed
bedtools getfasta -fi Chr5_14_16Mbp.fa -bed Chr5_14_16Mbp_output.bed > Chr5_14_16Mbp_output.2000.fasta

minimap2 -f 1000 -s 400 -ax ava-ont -d Chr5_14_16Mbp_output.fasta.mmi Chr5_14_16Mbp_output.2000.fasta
minimap2 -t 88 -f 10000 -s 400 -ax ava-ont --dual=yes --eqx Chr5_14_16Mbp_output.fasta.mmi Chr5_14_16Mbp_output.2000.fasta | samtools sort -m 4G -o Chr5_14_16Mbp_output.bam
python3.13 {SOFTWARE_LOCAL}/StainedGlass-0.6/workflow/scripts/samIdentity.py --threads 88 --matches 400 --header Chr5_14_16Mbp_output.bam > Chr5_14_16Mbp_output.tbl
bgzip -c Chr5_14_16Mbp_output.tbl > Chr5_14_16Mbp_output.tbl.gz
python3.13 {SOFTWARE_LOCAL}/StainedGlass-0.6/workflow/scripts/refmt.py --window 2000 --fai Chr5_14_16Mbp.fa.fai --full Chr5_14_16Mbp_output.full.tbl.gz Chr5_14_16Mbp_output.tbl.gz Chr5_14_16Mbp_output.bed.gz
mkdir -p results/Chr5_14_16Mbp_figures/pdfs
mkdir -p results/Chr5_14_16Mbp_figures/pngs
Rscript {SOFTWARE_LOCAL}/StainedGlass-0.6/workflow/scripts/aln_plot.R -b Chr5_14_16Mbp_output.bed.gz --threads 88 --prefix Chr5_14_16Mbp



bedtools getfasta -fi ../Lotus_GifuT2T_v1.0.fasta -bed Chr6_8_10Mbp.bed > Chr6_8_10Mbp.fa
samtools faidx Chr6_8_10Mbp.fa
bedtools makewindows -g Chr6_8_10Mbp.fa.fai -w 2000 > Chr6_8_10Mbp_output.bed
bedtools getfasta -fi Chr6_8_10Mbp.fa -bed Chr6_8_10Mbp_output.bed > Chr6_8_10Mbp_output.2000.fasta

minimap2 -f 1000 -s 400 -ax ava-ont -d Chr6_8_10Mbp_output.fasta.mmi Chr6_8_10Mbp_output.2000.fasta
minimap2 -t 88 -f 10000 -s 400 -ax ava-ont --dual=yes --eqx Chr6_8_10Mbp_output.fasta.mmi Chr6_8_10Mbp_output.2000.fasta | samtools sort -m 4G -o Chr6_8_10Mbp_output.bam
python3.13 {SOFTWARE_LOCAL}/StainedGlass-0.6/workflow/scripts/samIdentity.py --threads 88 --matches 400 --header Chr6_8_10Mbp_output.bam > Chr6_8_10Mbp_output.tbl
bgzip -c Chr6_8_10Mbp_output.tbl > Chr6_8_10Mbp_output.tbl.gz
python3.13 {SOFTWARE_LOCAL}/StainedGlass-0.6/workflow/scripts/refmt.py --window 2000 --fai Chr6_8_10Mbp.fa.fai --full Chr6_8_10Mbp_output.full.tbl.gz Chr6_8_10Mbp_output.tbl.gz Chr6_8_10Mbp_output.bed.gz
mkdir -p results/Chr6_8_10Mbp_figures/pdfs
mkdir -p results/Chr6_8_10Mbp_figures/pngs
Rscript {SOFTWARE_LOCAL}/StainedGlass-0.6/workflow/scripts/aln_plot.R -b Chr6_8_10Mbp_output.bed.gz --threads 88 --prefix Chr6_8_10Mbp


bedtools getfasta -fi ../Lotus_GifuT2T_v1.0.fasta -bed Chr6_25_28Mbp.bed > Chr6_25_28Mbp.fa
samtools faidx Chr6_25_28Mbp.fa
bedtools makewindows -g Chr6_25_28Mbp.fa.fai -w 2000 > Chr6_25_28Mbp_output.bed
bedtools getfasta -fi Chr6_25_28Mbp.fa -bed Chr6_25_28Mbp_output.bed > Chr6_25_28Mbp_output.2000.fasta

minimap2 -f 1000 -s 400 -ax ava-ont -d Chr6_25_28Mbp_output.fasta.mmi Chr6_25_28Mbp_output.2000.fasta
minimap2 -t 88 -f 10000 -s 400 -ax ava-ont --dual=yes --eqx Chr6_25_28Mbp_output.fasta.mmi Chr6_25_28Mbp_output.2000.fasta | samtools sort -m 4G -o Chr6_25_28Mbp_output.bam
python3.13 {SOFTWARE_LOCAL}/StainedGlass-0.6/workflow/scripts/samIdentity.py --threads 88 --matches 400 --header Chr6_25_28Mbp_output.bam > Chr6_25_28Mbp_output.tbl
bgzip -c Chr6_25_28Mbp_output.tbl > Chr6_25_28Mbp_output.tbl.gz
python3.13 {SOFTWARE_LOCAL}/StainedGlass-0.6/workflow/scripts/refmt.py --window 2000 --fai Chr6_25_28Mbp.fa.fai --full Chr6_25_28Mbp_output.full.tbl.gz Chr6_25_28Mbp_output.tbl.gz Chr6_25_28Mbp_output.bed.gz
mkdir -p results/Chr6_25_28Mbp_figures/pdfs
mkdir -p results/Chr6_25_28Mbp_figures/pngs
Rscript {SOFTWARE_LOCAL}/StainedGlass-0.6/workflow/scripts/aln_plot.R -b Chr6_25_28Mbp_output.bed.gz --threads 88 --prefix Chr6_25_28Mbp
```

</details>
