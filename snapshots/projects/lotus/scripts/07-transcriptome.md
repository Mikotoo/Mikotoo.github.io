---
title: '07 转录组'
sidebar:
  label: '07 转录组'
  order: -1
---


## 中文
**作用**：*Lotus japonicus* bulk RNA-seq——Gifu 与 MG20 组织图谱文库的预处理与定量、组织图谱与跨材料相关性、共表达分析（WGCNA、Mfuzz、Pearson 网络），以及突变体与公共数据定量。
**入口**：各步骤脚本见下表；阶段级 `run.sh` 调用 `Jaccard.py`。
**输入**：`00_prepare/align.conf` 与 Gifu、MG20 组织文库原始 FASTQ；`01_assembly`/`02_annotation` 的 GifuT2T FASTA 与 GTF；`07_coexpression_network/` 读取的基因、SNF 与同源列表；`01_tissue_atlas/` 与 `08_module_consistency/` 另需条件均值矩阵 `Gifu_all_samples_mean.tsv`、`MG20_all_samples_mean.tsv`。
**输出**：`00_prepare/all_samples.tpm.tsv` 供后续各步使用；各子目录的图表与表格见其 README。
**运行**：按 `00_prepare/` → `01_tissue_atlas/` → `05_wgcna/` → `06_mfuzz/` → `07_coexpression_network/` → `08_module_consistency/` 的顺序执行各目录脚本（各脚本按自身工作目录解析输入，不能写成一行）；路径占位符 ``{TOKEN}`` 见 `scripts/README.md`。
**工具**：`00_prepare/` 用 fastp、HISAT2、samtools、StringTie、GNU Awk；`01`、`02` 用 Python 3（pandas/numpy/scipy）与 R（tidyverse）；`05`–`08` 用 R 4.4.3（WGCNA、Mfuzz、Biobase、psych、data.table、igraph、ggraph、pheatmap、ggplot2）。
| 步骤 | 作用 |
|---|---|
| `00_prepare/` | fastp 质控、HISAT2 比对、StringTie 定量、TPM 矩阵 |
| `01_tissue_atlas/` | 组织图谱、SNF × 新基因相关、`spearman.R` 跨材料相关性 |
| `02_mutant_and_public_data/` | *Lotus* 结瘤突变体、大豆与 *Medicago truncatula* 定量 |
| `05_wgcna/` | 模块检测与合并、模块–性状相关、hub 基因 |
| `06_mfuzz/` | brown 与 lightcyan1 模块的时间序列软聚类 |
| `07_coexpression_network/` | Pearson 共表达网络 |
| `08_module_consistency/` | 模块基因的跨材料一致性 |
| `run.sh` + `Jaccard.py` | τ 组织特异性指数、Gifu/MG20 top-5000 基因集的 Jaccard 重叠 |

<details>
<summary>English</summary>

**Purpose**: Bulk RNA-seq of *Lotus japonicus* — preparation and quantification of the Gifu and MG20 tissue-atlas libraries, the atlas and its cross-accession correlation, co-expression analysis (WGCNA, Mfuzz, Pearson networks), and mutant and public-dataset quantification.
**Entry point**: the scripts listed in the table; the stage-level `run.sh` invokes `Jaccard.py`.
**Inputs**: `00_prepare/align.conf` and the raw FASTQ of the Gifu and MG20 tissue libraries; the GifuT2T FASTA/GTF from `01_assembly`/`02_annotation`; the gene, SNF and orthology lists read by `07_coexpression_network/`; `01_tissue_atlas/` and `08_module_consistency/` additionally need the condition-mean matrices `Gifu_all_samples_mean.tsv` and `MG20_all_samples_mean.tsv`.
**Outputs**: `00_prepare/all_samples.tpm.tsv`, consumed by every later step; the plots and tables of each sub-directory are listed in its README.
**Run**: run the scripts of `00_prepare/` → `01_tissue_atlas/` → `05_wgcna/` → `06_mfuzz/` → `07_coexpression_network/` → `08_module_consistency/` in that order (each script resolves its inputs relative to its own working directory, so the steps are not a one-liner); for the ``{TOKEN}`` path placeholders see `scripts/README.md`.
**Tools**: `00_prepare/` uses fastp, HISAT2, samtools, StringTie and GNU Awk; `01` and `02` use Python 3 (pandas/numpy/scipy) and R (tidyverse); `05`–`08` use R 4.4.3 (WGCNA, Mfuzz, Biobase, psych, data.table, igraph, ggraph, pheatmap, ggplot2).
| Step | Purpose |
|---|---|
| `00_prepare/` | fastp QC, HISAT2 alignment, StringTie quantification, TPM matrix |
| `01_tissue_atlas/` | tissue atlas, SNF × new-gene correlation, cross-accession correlation (`spearman.R`) |
| `02_mutant_and_public_data/` | *Lotus* nodulation mutants, soybean and *Medicago truncatula* quantification |
| `05_wgcna/` | module detection and merging, module–trait correlation, hub genes |
| `06_mfuzz/` | time-series soft clustering of the brown and lightcyan1 modules |
| `07_coexpression_network/` | Pearson co-expression network |
| `08_module_consistency/` | module gene consistency across accessions |
| `run.sh` + `Jaccard.py` | τ tissue-specificity index, Jaccard overlap of the Gifu/MG20 top-5000 gene sets |

</details>

### 本阶段步骤

| 步骤 | 内容 |
| --- | --- |
| [数据准备](/_projects/lotus/scripts/07-transcriptome/00-prepare/) | — |
| [组织表达图谱](/_projects/lotus/scripts/07-transcriptome/01-tissue-atlas/) | — |
| [突变体与公共数据](/_projects/lotus/scripts/07-transcriptome/02-mutant-and-public-data/) | — |
| [WGCNA 共表达模块](/_projects/lotus/scripts/07-transcriptome/05-wgcna/) | — |
| [Mfuzz 时序聚类](/_projects/lotus/scripts/07-transcriptome/06-mfuzz/) | — |
| [共表达网络](/_projects/lotus/scripts/07-transcriptome/07-coexpression-network/) | — |
| [模块一致性](/_projects/lotus/scripts/07-transcriptome/08-module-consistency/) | — |

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>Jaccard.py</code> · Python · 161 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/07_transcriptome/Jaccard.py)

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from pathlib import Path
import argparse
import sys

EPS = 1e-6

def calc_tau(row: pd.Series) -> float:
    """
    τ (tau) tissue-specificity index.
    row: expression across tissues (non-negative)
    """
    xmax = row.max()
    if xmax <= 0:
        return 0.0
    n = row.shape[0]
    return float(((1.0 - row / xmax).sum()) / (n - 1))

def build_tissue_matrix(df: pd.DataFrame, tissue_patterns: dict) -> pd.DataFrame:
    """
    df: gene_id as index, columns are sample names (TPM)
    tissue_patterns: {tissue: [substr1, substr2, ...]}  (OR match)
    returns: DataFrame (gene_id x tissue) with mean TPM across matched columns
    """
    tissue_expr = {}
    for tissue, pats in tissue_patterns.items():
        cols = [c for c in df.columns if any(p in c for p in pats)]
        if len(cols) == 0:
            raise ValueError(f"[ERROR] No columns matched for tissue '{tissue}' using patterns: {pats}")
        tissue_expr[tissue] = df[cols].mean(axis=1)
    return pd.DataFrame(tissue_expr)

def process_one(expr_tsv: str, label: str, outdir: str, topn: int,
                tissue_patterns: dict, min_tpm_max: float = 0.0) -> dict:
    """
    Read TPM matrix, compute tissue mean, tau, max_tissue, and export topN lists per tissue.
    Returns: dict {tissue: set(gene_id)}
    """
    df = pd.read_csv(expr_tsv, sep="\t", dtype={0: str})
    if "gene_id" not in df.columns:
        raise ValueError(f"[ERROR] '{expr_tsv}' has no column named 'gene_id'")
    df = df.set_index("gene_id")

    # Ensure numeric
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    # 1) tissue-level mean TPM
    tissue_df = build_tissue_matrix(df, tissue_patterns)

    # Optional filter: only keep genes expressed above some level in at least one tissue
    if min_tpm_max > 0:
        keep = tissue_df.max(axis=1) >= min_tpm_max
        tissue_df = tissue_df.loc[keep]

    # 2) tau + max tissue
    tau = tissue_df.apply(calc_tau, axis=1)
    max_tissue = tissue_df.idxmax(axis=1)

    res = tissue_df.copy()
    res["tau"] = tau
    res["max_tissue"] = max_tissue

    # 3) export topN per tissue
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    top_sets = {}
    for tissue in tissue_patterns.keys():
        sub = res[res["max_tissue"] == tissue].sort_values("tau", ascending=False)
        sub_top = sub.head(topn)
        genes = list(sub_top.index)
        top_sets[tissue] = set(genes)

        # write list
        (out / f"{label}_{tissue}_top{topn}.txt").write_text("\n".join(genes) + "\n", encoding="utf-8")

        # optional: write a table for debugging/traceability
        sub_top.to_csv(out / f"{label}_{tissue}_top{topn}.detail.tsv", sep="\t", index=True)

    # write a full summary table (tissue means + tau + max_tissue)
    res.to_csv(out / f"{label}.tissue_mean_tau.tsv", sep="\t", index=True)

    return top_sets

def jaccard(a: set, b: set) -> float:
    if len(a) == 0 and len(b) == 0:
        return 1.0
    if len(a) == 0 or len(b) == 0:
        return 0.0
    return len(a & b) / len(a | b)

def main():
    ap = argparse.ArgumentParser(
        description="Compute tau tissue-specificity, pick topN tissue-enriched genes, and compute Jaccard between Gifu and MG20."
    )
    ap.add_argument("--gifu", required=True, help="Gifu TPM matrix TSV (gene_id + samples)")
    ap.add_argument("--mg20", required=True, help="MG20 TPM matrix TSV (gene_id + samples)")
    ap.add_argument("--topn", type=int, default=500, help="Top N genes per tissue (default: 500)")
    ap.add_argument("--min_tpm_max", type=float, default=0.0,
                    help="Optional filter: keep genes with max tissue-mean TPM >= this value (default: 0, no filter)")
    ap.add_argument("--outprefix", default="tau_topN_jaccard", help="Output prefix (default: tau_topN_jaccard)")
    args = ap.parse_args()

    # ========= Tissue definitions (root_uni is separated) =========
    # NOTE: patterns are substrings; OR logic within each tissue
    tissue_patterns = {
        "flower":   ["flower_"],
        "leaf":     ["leaf_"],
        "nodule":   ["nodule_"],
        "pod":      ["pod_"],
        "stem":     ["stem_"],

        # root excludes root_uni explicitly by using dpi/hpi patterns
        "root":     ["root_dpi", "root_hpi"],

        # root_uni is its own group
        "root_uni": ["root_uni_"],
    }

    gifu_out = f"{args.outprefix}.Gifu"
    mg20_out = f"{args.outprefix}.MG20"

    print(f"[INFO] Processing Gifu: {args.gifu}")
    gifu_sets = process_one(
        args.gifu, "Gifu", gifu_out, args.topn, tissue_patterns, min_tpm_max=args.min_tpm_max
    )

    print(f"[INFO] Processing MG20: {args.mg20}")
    mg20_sets = process_one(
        args.mg20, "MG20", mg20_out, args.topn, tissue_patterns, min_tpm_max=args.min_tpm_max
    )

    # ========= Jaccard per tissue =========
    rows = []
    for tissue in tissue_patterns.keys():
        A = gifu_sets.get(tissue, set())
        B = mg20_sets.get(tissue, set())
        inter = len(A & B)
        union = len(A | B)
        jac = jaccard(A, B)
        rows.append([tissue, len(A), len(B), inter, union, jac])

    out = pd.DataFrame(rows, columns=["tissue", "Gifu_topN", "MG20_topN", "intersection", "union", "jaccard"])
    out_tsv = f"{args.outprefix}.jaccard_by_tissue.tsv"
    out.to_csv(out_tsv, sep="\t", index=False)

    print("[OK] Jaccard results:")
    print(out.to_string(index=False))
    print(f"[OUT] {out_tsv}")
    print(f"[OUT] TopN lists and details: {gifu_out}/ , {mg20_out}/")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(str(e), file=sys.stderr)
        sys.exit(1)
```

</details>

<details open>
<summary><code>run.sh</code> · Shell · 6 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/07_transcriptome/run.sh)

```bash
python3 Jaccard.py \
  --gifu 01.Gifu_tissue/Gifu_all_samples.tpm.tsv \
  --mg20 02.MG20.tissue/MG20_all_samples.tpm.tsv \
  --topn 5000 \
  --min_tpm_max 0 \
  --outprefix Lotus_tau_top500
```

</details>
