---
title: '基因预测（BRAKER）'
sidebar:
  label: '基因预测（BRAKER）'
  order: 4
---


## 中文

**作用**：在 EDTA 屏蔽后的 Gifu 基因组上用蛋白、Iso-Seq 与 RNA-seq 证据跑 BRAKER3，再按最小长度、N 端甲硫氨酸、无内部终止密码子过滤蛋白及其对应的 CDS 与 GFF3。
**入口**：`Gifu_braker.sh`（BRAKER3）；`makeblastdb.sh`、`filter.sh`、`bedtools.sh`、`filter_gff_overlaps.py` 建库、过滤与去冗余。
**输入**：屏蔽基因组 ``{PROJ_LOTUS_ZC}`/Gifu/09.annotation/03.edta/Gifu_v1.0.fasta.mod.MAKER.masked`（来自 `01_repeat_edta`）；蛋白证据 `.../06.braker/douke_pep.fasta`；Iso-Seq BAM `.../05.isoseq/Gifu_isoforms.sorted.bam`；RNA-seq BAM 在 `.../04.rnasq/01.align/`（逗号分隔的 `--bam` 列表在 `Gifu_braker.sh` 第 11 行）；`filter.sh` 另读 `braker.aa`、`braker.codingseq`、`braker.gff3`，`bedtools.sh` 读 `gene_filter.bed`。
**输出**：`.../06.braker` 内的 `braker.aa`、`braker.codingseq`、`braker.gff3`；`Gifu_filtered.pep.fa`、`Gifu_filtered.cds.fa`、`Gifu_filtered.gff3`（`02_rnaseq_align` 用的是 `Gifu_filtered.gtf`，需自行转换）；`overlaps.tsv`；`filter_gff_overlaps.py` 写其第二个位置参数。
**运行**：`bash makeblastdb.sh` → `bash Gifu_braker.sh` → `bash filter.sh` → `python3 filter_gff_overlaps.py <in.gff3> <out.gff3>`；`filter.sh` 调 `python AnnotationFilter.py -l 20 -i braker.aa -o Gifu_filtered.pep.fa -c braker.codingseq -C Gifu_filtered.cds.fa -g braker.gff3 -G Gifu_filtered.gff3`；`bedtools.sh` 跑 `bedtools intersect -a gene_filter.bed -b gene_filter.bed -wa -wb | awk -F'\t' '$4 < $9' > overlaps.tsv`。
**工具**：BRAKER3 `braker.pl --species=Gifu_1 --threads 88 --gff3`（`singularity exec braker3_sandbox`，镜像是本地镜像名而非路径）；AUGUSTUS `--AUGUSTUS_CONFIG_PATH=`{SOFTWARE_ZC}`/Augustus-master/config/`；BLAST+ `makeblastdb -in douke_pep.fasta -dbtype prot -out douke_pep.fasta`；bedtools `intersect -wa -wb`；`AnnotationFilter.py -l 20`；`filter_gff_overlaps.py --strand-aware`。

<details>
<summary>English</summary>

**Purpose**: Run BRAKER3 on the EDTA-masked Gifu genome with protein, Iso-Seq and RNA-seq evidence, then filter the predicted proteins (minimum length, N-terminal methionine, no internal stop codon) together with their matching CDS and GFF3 features.
**Entry point**: `Gifu_braker.sh` (BRAKER3); `makeblastdb.sh`, `filter.sh`, `bedtools.sh` and `filter_gff_overlaps.py` build the BLAST database, filter and reduce redundancy.
**Inputs**: masked genome ``{PROJ_LOTUS_ZC}`/Gifu/09.annotation/03.edta/Gifu_v1.0.fasta.mod.MAKER.masked` (from `01_repeat_edta`); protein evidence `.../06.braker/douke_pep.fasta`; Iso-Seq BAM `.../05.isoseq/Gifu_isoforms.sorted.bam`; RNA-seq BAMs in `.../04.rnasq/01.align/` (the comma-separated `--bam` list is on line 11 of `Gifu_braker.sh`); `filter.sh` also reads `braker.aa`, `braker.codingseq`, `braker.gff3`, and `bedtools.sh` reads `gene_filter.bed`.
**Outputs**: `braker.aa`, `braker.codingseq`, `braker.gff3` in `.../06.braker`; `Gifu_filtered.pep.fa`, `Gifu_filtered.cds.fa`, `Gifu_filtered.gff3` (note that `02_rnaseq_align` uses `Gifu_filtered.gtf`); `overlaps.tsv`; `filter_gff_overlaps.py` writes its second positional argument.
**Run**: `bash makeblastdb.sh` → `bash Gifu_braker.sh` → `bash filter.sh` → `python3 filter_gff_overlaps.py <in.gff3> <out.gff3>`; `filter.sh` calls `python AnnotationFilter.py -l 20 -i braker.aa -o Gifu_filtered.pep.fa -c braker.codingseq -C Gifu_filtered.cds.fa -g braker.gff3 -G Gifu_filtered.gff3`; `bedtools.sh` runs `bedtools intersect -a gene_filter.bed -b gene_filter.bed -wa -wb | awk -F'\t' '$4 < $9' > overlaps.tsv`.
**Tools**: BRAKER3 `braker.pl --species=Gifu_1 --threads 88 --gff3` (via `singularity exec braker3_sandbox`, a local image name rather than a path); AUGUSTUS `--AUGUSTUS_CONFIG_PATH=`{SOFTWARE_ZC}`/Augustus-master/config/`; BLAST+ `makeblastdb -in douke_pep.fasta -dbtype prot -out douke_pep.fasta`; bedtools `intersect -wa -wb`; `AnnotationFilter.py -l 20`; `filter_gff_overlaps.py --strand-aware`.

</details>

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>AnnotationFilter.py</code> · Python · 147 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/02_annotation/04_braker/AnnotationFilter.py)

```python
#!/usr/bin/env python3
"""
Filter predicted proteins and corresponding CDS and GFF3 features:
    - Keep proteins >= min length (default 50 aa), starting with M, and with no internal stop codons.
    - Output filtered proteins and (optionally) matching CDS.
    - Output GFF3 with removed features whose ID or Parent matches filtered-out proteins.
"""

import argparse
from Bio import SeqIO
import re
from collections import defaultdict

def parse_args():
    ap = argparse.ArgumentParser(
        description="Filter proteins by length, N-terminal methionine, and premature stop codon; filter corresponding CDS and GFF3."
    )
    ap.add_argument("-i", "--input", required=True, help="protein fasta")
    ap.add_argument("-o", "--out", required=True, help="output protein fasta")
    ap.add_argument("-l", "--length", type=int, default=20, help="minimum protein length")
    ap.add_argument("-c", "--cds", required=False, help="input CDS fasta (optional)")
    ap.add_argument("-C", "--out_cds", required=False, help="output filtered CDS fasta (optional)")
    ap.add_argument("-g", "--gff", required=False, help="input GFF3 (optional)")
    ap.add_argument("-G", "--out_gff", required=False, help="output filtered GFF3 (optional)")
    return ap.parse_args()

def protein_pass(protein_seq, min_length):
    """Return True if sequence passes all criteria."""
    s = str(protein_seq)
    if len(s) < min_length:
        return False
    if s[0] != "M":
        return False
    # Allow only a single terminal stop, no internal
    if "*" in s[:-1]:
        return False
    return True

def filter_proteins_and_write(input_fasta, min_length, output_fasta):
    """Return set of IDs of proteins that pass filters and write them."""
    passing_ids = set()
    with open(output_fasta, "w") as out:
        for rec in SeqIO.parse(input_fasta, "fasta"):
            if protein_pass(rec.seq, min_length):
                SeqIO.write(rec, out, "fasta")
                passing_ids.add(rec.id.split()[0])
    return passing_ids

def filter_cds_and_write(cds_fasta, passing_ids, output_fasta):
    n = 0
    with open(output_fasta, "w") as out:
        for rec in SeqIO.parse(cds_fasta, "fasta"):
            if rec.id.split()[0] in passing_ids:
                SeqIO.write(rec, out, "fasta")
                n += 1
    return n

def parse_attr(attr_field):
    """Parse GFF3 attribute field into a dict."""
    attrs = {}
    for part in attr_field.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k] = v
    return attrs

def collect_transcript_gene_relations(gff_file):
    """Return mappings: transcript → gene and gene → [transcripts]."""
    t2g = {}
    g2t = defaultdict(list)
    with open(gff_file) as fh:
        for ln in fh:
            if ln.startswith("#") or ln.strip() == "":
                continue
            flds = ln.rstrip().split("\t")
            if len(flds) < 9:
                continue
            ftype, attrs = flds[2], flds[8]
            a = parse_attr(attrs)
            if ftype.lower() in {"mrna", "transcript"}:
                tid = a.get("ID")
                gid = a.get("Parent")
                if tid and gid:
                    t2g[tid] = gid
                    g2t[gid].append(tid)
    return t2g, g2t

def decide_gene_removal(g2t, filtered_out_ids):
    """Genes to drop = every transcript filtered out."""
    drop_genes = set()
    for gid, tids in g2t.items():
        if all(tid in filtered_out_ids for tid in tids):
            drop_genes.add(gid)
    return drop_genes

def filter_gff(gff_in, gff_out, passing_ids, drop_genes):
    kept, removed = 0, 0
    with open(gff_in) as inp, open(gff_out, "w") as out:
        for ln in inp:
            if ln.startswith("#") or ln.strip() == "":
                out.write(ln)
                continue
            flds = ln.rstrip().split("\t")
            if len(flds) < 9:
                out.write(ln)
                continue
            ftype, attrs = flds[2], flds[8]
            a = parse_attr(attrs)
            gid = a.get("ID") if ftype.lower() == "gene" else None
            tid = None
            if ftype.lower() in {"mrna", "transcript"}:
                tid = a.get("ID")
            elif "Parent" in a:
                tid = a["Parent"].split(",")[0]  # first parent

            # drop rules
            # Remove gene if all transcripts filtered out
            if gid and gid in drop_genes:
                removed += 1
                continue
            # Remove transcript and all features whose Parent matches
            if tid and tid not in passing_ids:
                removed += 1
                continue

            out.write(ln)
            kept += 1
    return kept, removed

def main():
    args = parse_args()
    passing_ids = filter_proteins_and_write(args.input, args.length, args.out)
    print(f"{len(passing_ids)} proteins passed all filters and written to {args.out}")
    if args.cds and args.out_cds:
        n = filter_cds_and_write(args.cds, passing_ids, args.out_cds)
        print(f"{n} CDS records matched filtered proteins and written to {args.out_cds}")
    if args.gff and args.out_gff:
        t2g, g2t = collect_transcript_gene_relations(args.gff)
        filtered_out_ids = set(t2g.keys()) - passing_ids
        drop_genes = decide_gene_removal(g2t, filtered_out_ids)
        kept, removed = filter_gff(args.gff, args.out_gff, passing_ids, drop_genes)
        print(f"GFF filtered: kept {kept} lines; removed {removed} lines → {args.out_gff}")

if __name__ == "__main__":
    main()
```

</details>

<details open>
<summary><code>Gifu_braker.sh</code> · Shell · 11 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/02_annotation/04_braker/Gifu_braker.sh)

```bash
#!/bin/bash
#CSUB -J singularity
#CSUB -q c01
#CSUB -o singularity.out
#CSUB -e singularity.error
#CSUB -n 88
#CSUB -R span[hosts=1]

singularity exec braker3_sandbox braker.pl  --genome={PROJ_LOTUS_ZC}/Gifu/09.annotation/03.edta/Gifu_v1.0.fasta.mod.MAKER.masked --species=Gifu_1 \
 --threads 88 --gff3  --workingdir={PROJ_LOTUS_ZC}/Gifu/09.annotation/06.braker --AUGUSTUS_CONFIG_PATH={SOFTWARE_ZC}/Augustus-master/config/  \
 --prot_seq={PROJ_LOTUS_ZC}/Gifu/09.annotation/06.braker/douke_pep.fasta --bam={PROJ_LOTUS_ZC}/Gifu/09.annotation/05.isoseq/Gifu_isoforms.sorted.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_flower_rep1.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_flower_rep2.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_flower_rep3.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_leaf_rep1.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_leaf_rep2.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_leaf_rep3.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_nodule_dpi10_rep1.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_nodule_dpi10_rep2.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_nodule_dpi10_rep3.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_nodule_dpi21_rep1.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_nodule_dpi21_rep2.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_nodule_dpi21_rep3.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_nodule_dpi45_rep1.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_nodule_dpi45_rep2.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_nodule_dpi45_rep3.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_pod_rep1.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_pod_rep2.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_pod_rep3.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_root_dpi02_rep1.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_root_dpi02_rep2.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_root_dpi02_rep3.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_root_dpi10_rep1.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_root_dpi10_rep2.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_root_dpi10_rep3.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_root_dpi21_rep1.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_root_dpi21_rep2.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_root_dpi21_rep3.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_root_dpi45_rep1.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_root_dpi45_rep2.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_root_dpi45_rep3.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_root_hpi48_rep1.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_root_hpi48_rep2.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_root_hpi48_rep3.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_root_uni_rep1.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_root_uni_rep2.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_root_uni_rep3.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_stem_rep1.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_stem_rep2.bam,{PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/01.align/Gifu_stem_rep3.bam
```

</details>

<details open>
<summary><code>bedtools.sh</code> · Shell · 3 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/02_annotation/04_braker/bedtools.sh)

```bash
bedtools intersect -a gene_filter.bed -b gene_filter.bed -wa -wb \
| awk -F'\t' '$4 < $9 {print $0}' \
> overlaps.tsv
```

</details>

<details open>
<summary><code>filter.sh</code> · Shell · 12 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/02_annotation/04_braker/filter.sh)

```bash
#!/bin/bash
#CSUB -J filter
#CSUB -q c01
#CSUB -o filter.out
#CSUB -e filter.error
#CSUB -n 88
#CSUB -R span[hosts=1]

python AnnotationFilter.py -l 20 \
-i braker.aa -o Gifu_filtered.pep.fa \
-c braker.codingseq -C Gifu_filtered.cds.fa \
-g braker.gff3 -G Gifu_filtered.gff3
```

</details>

<details>
<summary><code>filter_gff_overlaps.py</code> · Python · 334 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/02_annotation/04_braker/filter_gff_overlaps.py)

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filter overlapping/duplicate genes from a GFF3, while preserving *all* child features
(mRNA/exon/CDS/UTR/...) of the kept genes and keeping output format/order close to input.

Usage:
  python filter_gff_overlaps_full.py in.gff3 out.gff3 \
    --strand-aware \
    --reciprocal-threshold 0.0 \
    --keep-criterion longest

Options:
  --strand-aware            Only treat same-strand genes as overlapping (default off).
  --reciprocal-threshold    Reciprocal overlap threshold in [0,1]; 0 = any overlap (default 0).
  --keep-criterion          Which gene to keep within a cluster: longest | id-min (default longest).
"""

import sys
import argparse
from collections import defaultdict, OrderedDict, deque

# -------------------- GFF utils --------------------
def parse_attrs(attr_field: str) -> OrderedDict:
    attrs = OrderedDict()
    for entry in attr_field.split(";"):
        entry = entry.strip()
        if not entry:
            continue
        if "=" in entry:
            k, v = entry.split("=", 1)
            attrs[k] = v
        else:
            attrs[entry] = ""
    return attrs

def attrs_to_str(attrs: OrderedDict) -> str:
    return ";".join(f"{k}={v}" if v != "" else k for k, v in attrs.items())

def length1(s: int, e: int) -> int:
    return e - s + 1  # GFF is 1-based inclusive

def overlap_len(a, b) -> int:
    # a,b are dicts with start,end
    return max(0, min(a["end"], b["end"]) - max(a["start"], b["start"]) + 1)

def reciprocal_overlap(a, b) -> float:
    inter = overlap_len(a, b)
    if inter == 0:
        return 0.0
    la = length1(a["start"], a["end"])
    lb = length1(b["start"], b["end"])
    return min(inter / la, inter / lb)

# -------------------- Parsing whole GFF --------------------
def read_gff(path):
    """
    Return:
      headers: list[str]                 # lines starting with '#' or malformed kept as-is
      records: dict[id]->rec             # every feature with ID (gene/mRNA/exon/CDS/... if they have ID)
      anon_lines: list[(lineno,line)]    # features without ID (会作为“匿名记录”保留，但仅在其 Parent 最终被保留时输出)
      parent_of: dict[id]->list[parent_ids]
      children_of: dict[id]->list[child_ids]
      genes: dict[gid]->rec              # only feature type == 'gene'
      order_index: dict[id]->lineno      # first appearance order for features with ID
    """
    headers = []
    records = {}          # id -> rec
    anon_lines = []       # lines without ID but may have Parent (exon 等常有 ID，也可能没有)
    parent_of = defaultdict(list)    # child_id -> [parent_id,...]
    children_of = defaultdict(list)  # parent_id -> [child_id,...]
    genes = {}
    order_index = {}

    with open(path, "r") as fh:
        lineno = 0
        for raw in fh:
            lineno += 1
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                headers.append(line)
                continue

            parts = line.split("\t")
            if len(parts) < 9:
                # 非标准行，当作 header 原样保留
                headers.append(line)
                continue

            seqid, source, ftype, start, end, score, strand, phase, attr = parts
            try:
                start_i = int(start); end_i = int(end)
            except Exception:
                headers.append(line)
                continue

            attrs = parse_attrs(attr)
            fid = attrs.get("ID", None)
            parents = []
            if "Parent" in attrs:
                parents = [p for p in attrs["Parent"].split(",") if p]

            rec = {
                "chr": seqid,
                "source": source,
                "type": ftype,
                "start": start_i,
                "end": end_i,
                "score": score,
                "strand": strand if strand in "+-" else ".",
                "phase": phase,
                "attrs": attrs,
                "line": line,
                "lineno": lineno,
            }

            if fid:
                records[fid] = rec
                if fid not in order_index:
                    order_index[fid] = lineno
                for p in parents:
                    parent_of[fid].append(p)
                    children_of[p].append(fid)
                if ftype.lower() == "gene":
                    genes[fid] = rec
            else:
                # 没有 ID 的行：先收集，等到输出阶段根据其 Parent 是否被保留决定是否写出
                anon_lines.append((lineno, rec, parents))

    return headers, records, anon_lines, parent_of, children_of, genes, order_index

# -------------------- Dedup & clustering on genes --------------------
def dedup_exact(genes_items):
    """同 chr,start,end,strand 完全一致视为重复，只保留一个（按最长、再按 ID 排序）。"""
    buckets = defaultdict(list)
    for gid, g in genes_items:
        key = (g["chr"], g["start"], g["end"], g["strand"])
        buckets[key].append((gid, g))

    kept, dropped = [], set()
    for key, lst in buckets.items():
        if len(lst) == 1:
            kept.append(lst[0])
        else:
            lst.sort(key=lambda x: (-(x[1]["end"] - x[1]["start"]), x[0]))
            kept.append(lst[0])
            for gid, _ in lst[1:]:
                dropped.add(gid)
    return kept, dropped

def build_gene_clusters(genes_items, strand_aware=False, reciprocal_threshold=0.0):
    """把基因连通为簇。"""
    groups = defaultdict(list)  # per (chr[,strand])
    for gid, g in genes_items:
        key = (g["chr"], g["strand"]) if strand_aware else (g["chr"], None)
        groups[key].append((gid, g))

    clusters = []
    for key, arr in groups.items():
        arr.sort(key=lambda x: x[1]["start"])
        n = len(arr)
        visited = [False]*n
        for i in range(n):
            if visited[i]:
                continue
            comp = []
            dq = deque([i]); visited[i] = True
            while dq:
                u = dq.popleft()
                comp.append(u)
                urec = arr[u][1]
                # 只需检查 start <= u.end 的右侧邻域
                v = u + 1
                while v < n and arr[v][1]["start"] <= urec["end"]:
                    if not visited[v]:
                        vrec = arr[v][1]
                        ok = False
                        if reciprocal_threshold <= 0:
                            ok = not (vrec["start"] > urec["end"] or vrec["end"] < urec["start"])
                        else:
                            ok = (reciprocal_overlap(urec, vrec) >= reciprocal_threshold)
                        if ok:
                            visited[v] = True
                            dq.append(v)
                    v += 1
            clusters.append([arr[idx] for idx in comp])
    return clusters

def choose_rep(cluster, keep_criterion="longest"):
    if keep_criterion == "id-min":
        cluster.sort(key=lambda x: x[0])
        return cluster[0]
    # longest
    cluster.sort(key=lambda x: (-(x[1]["end"] - x[1]["start"] + 1), x[0]))
    return cluster[0]

# -------------------- Collect descendants & rewrite --------------------
def collect_descendants(root_gene_id, children_of, records):
    """返回以 gene 为根的所有后代（包含 gene 自身）的 ID 集合。"""
    keep_ids = set()
    stack = [root_gene_id]
    while stack:
        pid = stack.pop()
        if pid in keep_ids:
            continue
        keep_ids.add(pid)
        for cid in children_of.get(pid, []):
            if cid in records:       # 只收有 ID 的
                stack.append(cid)
    return keep_ids

def prune_parent_attr(rec, kept_ids_set):
    """修剪 Parent= 使之只保留仍存在的父 ID；若清空则返回 False 表示应丢弃该行。"""
    attrs = rec["attrs"]
    if "Parent" not in attrs:
        return True
    parents = [p for p in attrs["Parent"].split(",") if p]
    new_parents = [p for p in parents if p in kept_ids_set]
    if not new_parents:
        return False
    if ",".join(new_parents) != attrs["Parent"]:
        attrs = attrs.copy()
        attrs["Parent"] = ",".join(new_parents)
        rec["attrs"] = attrs
    return True

def rec_to_line(rec):
    # 用当前 attrs 重建第 9 列，并保持其他列不变
    parts = [
        rec["chr"],
        rec["source"],
        rec["type"],
        str(rec["start"]),
        str(rec["end"]),
        rec["score"],
        rec["strand"],
        rec["phase"],
        attrs_to_str(rec["attrs"])
    ]
    return "\t".join(parts)

# -------------------- Main --------------------
def main():
    ap = argparse.ArgumentParser(description="Filter overlapping genes from GFF3 and keep full structures.")
    ap.add_argument("input")
    ap.add_argument("output")
    ap.add_argument("--strand-aware", action="store_true", help="Overlap only within same strand.")
    ap.add_argument("--reciprocal-threshold", type=float, default=0.0,
                    help="Reciprocal overlap threshold in [0,1]. 0: any overlap.")
    ap.add_argument("--keep-criterion", choices=["longest","id-min"], default="longest",
                    help="Representative in a cluster.")
    args = ap.parse_args()

    headers, records, anon_lines, parent_of, children_of, genes, order_index = read_gff(args.input)

    # 1) 去完全重复
    gene_items = [(gid, genes[gid]) for gid in genes]
    dedup_kept, dropped_exact = dedup_exact(gene_items)

    # 2) 基于重叠建立簇并选择代表
    clusters = build_gene_clusters(dedup_kept,
                                   strand_aware=args.strand_aware,
                                   reciprocal_threshold=args.reciprocal_threshold)
    kept_gene_ids = set()
    for cl in clusters:
        gid, _ = choose_rep(cl, keep_criterion=args.keep_criterion)
        kept_gene_ids.add(gid)

    # 3) 收集每个保留基因的所有后代（递归）
    kept_ids = set()
    for gid in kept_gene_ids:
        kept_ids |= collect_descendants(gid, children_of, records)
    # 顺带把基因自身也包含在 kept_ids 里（collect 已含）

    # 4) 对匿名行（无 ID，但有 Parent 的），仅当 Parent 之一被保留时写出，并修剪 Parent
    #    对有 ID 的行：如果在 kept_ids，则写出；否则丢弃
    output_items = []  # (lineno, line_str)
    # 4.1 有 ID 的
    for fid in kept_ids:
        if fid not in records:
            continue
        rec = records[fid]
        # 修剪 Parent（例如 mRNA 的 Parent=gene；若 gene 被保留则 OK；若还有其他被删的 Parent 则去掉）
        if not prune_parent_attr(rec, kept_ids):
            continue
        # 只对 Parent 变动的行重构；否则尽量用原行，减少改动
        line_out = rec_to_line(rec) if "Parent" in rec["attrs"] and line_parent_changed(rec["line"], rec["attrs"]) else rec["line"]
        output_items.append((rec["lineno"], line_out))

    # 4.2 无 ID 的（多数情况下也会有 Parent）
    for lineno, rec, parents in anon_lines:
        if not parents:
            # 无 Parent 的匿名行：通常与结构无关，保留原样
            output_items.append((lineno, rec["line"]))
        else:
            keep = [p for p in parents if p in kept_ids]
            if keep:
                # 修剪 Parent
                if "attrs" in rec and "Parent" in rec["attrs"]:
                    attrs = rec["attrs"].copy()
                    attrs["Parent"] = ",".join(keep)
                    rec2 = rec.copy(); rec2["attrs"] = attrs
                    output_items.append((lineno, rec_to_line(rec2)))
                else:
                    output_items.append((lineno, rec["line"]))

    # 5) 拼装输出（头注释 + 按原行号排序后的主体；若存在 ##FASTA 段在 headers 里也会被保留）
    output_items.sort(key=lambda x: x[0])
    with open(args.output, "w") as out:
        for h in headers:
            out.write(h + "\n")
        for _, ln in output_items:
            out.write(ln + "\n")

    # 统计
    sys.stderr.write(f"[INFO] Input genes: {len(genes)}\n")
    sys.stderr.write(f"[INFO] Removed exact duplicates: {len(dropped_exact)}\n")
    sys.stderr.write(f"[INFO] Kept genes (after overlap filtering): {len(kept_gene_ids)}\n")
    sys.stderr.write(f"[INFO] Total features written: {len(output_items)}\n")
    sys.stderr.write(f"[INFO] Output: {args.output}\n")

def line_parent_changed(orig_line: str, new_attrs: OrderedDict) -> bool:
    """轻量判断 Parent 是否被修改，决定是否重建第9列（尽量减少改动）。"""
    try:
        old_attr = orig_line.split("\t", 8)[8]
    except Exception:
        return True
    old = parse_attrs(old_attr)
    return old.get("Parent","") != new_attrs.get("Parent","")

if __name__ == "__main__":
    main()
```

</details>

<details open>
<summary><code>makeblastdb.sh</code> · Shell · 9 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/02_annotation/04_braker/makeblastdb.sh)

```bash
#!/bin/bash
#CSUB -J blast
#CSUB -q c01
#CSUB -o blast.out
#CSUB -e blast.error
#CSUB -n 88
#CSUB -R span[hosts=1]

makeblastdb -in douke_pep.fasta -dbtype prot -out douke_pep.fasta
```

</details>
