---
title: '组织表达图谱'
sidebar:
  label: '组织表达图谱'
  order: 1
---


## 中文
**作用**：Gifu 组织图谱分析——把 SNF 家族基因与全图谱中相关性最高的 top-5“新”基因配对，并计算 Gifu 与 MG20 的跨材料表达相关性。
**入口**：`snf_new_corr.sh`、`spearman.R`；`getCount.sh` 只有计数任务的投递头，没有命令体。
**输入**：`snf_new_corr.sh` 需要 `expr.tsv`（TPM 矩阵，首列为基因 ID）与 `genes.tsv`（GeneID、Category，Category 含 `SNF` 或 `New`）；`spearman.R` 需要 `Gifu_all_samples_mean.tsv` 与 `MG20_all_samples_mean.tsv`。
**输出**：`<outprefix>.snf_new_corr/` 下的 `*.SNF_vs_New.all_pairs.tsv`、`*.SNF_top5_newgenes.tsv`、`*.SNF_top5_matrix.tsv`、`*.summary.txt`、`*.SNF_top5_heatmap.pdf|.png`、`*.SNF_top5_dotplot.pdf|.png`；`spearman.R` 写 `Gifu_MG20.sample_expression_spearman.tsv` 与 `Gifu_MG20.organ_expression_spearman.tsv` 及各自的 PDF。
**运行**：`bash snf_new_corr.sh expr.tsv genes.tsv outprefix`（须以本目录为工作目录）；`Rscript spearman.R`
**工具**：Python 3 + pandas/numpy/scipy（Pearson r 与 t 检验 p）；R + ggplot2/readr/dplyr/tidyr 与 tidyverse；集群投递头 `#CSUB -q c01 -n 64`。

<details>
<summary>English</summary>

**Purpose**: Tissue-atlas analysis for the Gifu accession — pairing SNF-family genes with the top-5 most correlated "new" genes across the atlas, and computing the cross-accession expression correlation between Gifu and MG20.
**Entry point**: `snf_new_corr.sh`, `spearman.R`; `getCount.sh` holds only the job header of the counting step and has no command body.
**Inputs**: `snf_new_corr.sh` needs `expr.tsv` (TPM matrix, first column = gene ID) and `genes.tsv` (GeneID and Category; Category containing `SNF` or `New`); `spearman.R` needs `Gifu_all_samples_mean.tsv` and `MG20_all_samples_mean.tsv`.
**Outputs**: in `<outprefix>.snf_new_corr/` — `*.SNF_vs_New.all_pairs.tsv`, `*.SNF_top5_newgenes.tsv`, `*.SNF_top5_matrix.tsv`, `*.summary.txt`, `*.SNF_top5_heatmap.pdf|.png` and `*.SNF_top5_dotplot.pdf|.png`; `spearman.R` writes `Gifu_MG20.sample_expression_spearman.tsv` and `Gifu_MG20.organ_expression_spearman.tsv`, each with a PDF.
**Run**: `bash snf_new_corr.sh expr.tsv genes.tsv outprefix` (run it with this directory as the working directory); `Rscript spearman.R`
**Tools**: Python 3 with pandas/numpy/scipy (Pearson r and t-test p); R with ggplot2/readr/dplyr/tidyr and tidyverse; cluster job header `#CSUB -q c01 -n 64`.

</details>

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>getCount.sh</code> · Shell · 7 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/07_transcriptome/01_tissue_atlas/getCount.sh)

```bash
#!/bin/bash
#CSUB -J count
#CSUB -q c01
#CSUB -o count.out
#CSUB -e count.error
#CSUB -n 64
#CSUB -R span[hosts=1]
```

</details>

<details open>
<summary><code>snf_new_corr.sh</code> · Shell · 173 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/07_transcriptome/01_tissue_atlas/snf_new_corr.sh)

```bash
#!/bin/bash
set -euo pipefail

if [ $# -lt 3 ]; then
  echo "Usage: bash $0 expr.tsv genes.tsv outprefix"
  exit 1
fi

EXPR="$1"
GENES="$2"
OUT="$3"

ODIR="${OUT}.snf_new_corr"
mkdir -p "$ODIR"

# =========================
# 1) Python: corr + top5
# =========================
python3 - "$EXPR" "$GENES" "$OUT" <<'PY'
import sys, os
import pandas as pd
import numpy as np
from scipy.stats import t as tdist

expr, genes, outp = sys.argv[1:4]
odir = f"{outp}.snf_new_corr"
os.makedirs(odir, exist_ok=True)

df = pd.read_csv(expr, sep="\t", header=0)
df = df.rename(columns={df.columns[0]:"GeneID"})
df = df.drop_duplicates("GeneID").set_index("GeneID")

meta = pd.read_csv(genes, sep="\t", header=None)
if meta.shape[1] < 2:
    raise SystemExit("gene.list 至少2列：GeneID 与 Category")
meta.columns = ["GeneID","Category"] + [f"c{i}" for i in range(3, meta.shape[1]+1)]
meta = meta.drop_duplicates("GeneID").set_index("GeneID")

snf_ids = meta.index[meta["Category"].astype(str).str.contains("SNF", case=False, na=False)].tolist()
new_ids = meta.index[meta["Category"].astype(str).str.contains("New", case=False, na=False)].tolist()

snf_ids = [g for g in snf_ids if g in df.index]
new_ids = [g for g in new_ids if g in df.index]

if len(snf_ids) == 0 or len(new_ids) == 0:
    raise SystemExit(f"在表达矩阵中匹配到 SNF={len(snf_ids)}, New={len(new_ids)}；检查 GeneID 是否一致")

X = df.loc[snf_ids].to_numpy(float)
Y = df.loc[new_ids].to_numpy(float)
n = X.shape[1]
dfree = n - 2

def zscore(A):
    mu = A.mean(axis=1, keepdims=True)
    sd = A.std(axis=1, ddof=0, keepdims=True)
    sd[sd == 0] = np.nan
    return (A - mu) / sd

Xz = zscore(X)
Yz = zscore(Y)

R = (Xz @ Yz.T) / n
Rp = np.clip(R, -0.999999999, 0.999999999)
tval = Rp * np.sqrt(dfree / (1.0 - Rp**2))
P = 2.0 * tdist.sf(np.abs(tval), df=dfree)

out_long = pd.DataFrame({
    "SNF_gene": np.repeat(snf_ids, len(new_ids)),
    "New_gene": np.tile(new_ids, len(snf_ids)),
    "r": R.reshape(-1),
    "p": P.reshape(-1),
}).dropna(subset=["r"])

odir = f"{outp}.snf_new_corr"
out_long.to_csv(f"{odir}/{outp}.SNF_vs_New.all_pairs.tsv", sep="\t", index=False)

top_list = []
for i, snf in enumerate(snf_ids):
    rr = R[i, :]
    pp = P[i, :]
    idx = np.argsort(-np.abs(rr))[:5]
    for j in idx:
        if np.isnan(rr[j]): 
            continue
        top_list.append([snf, new_ids[j], rr[j], pp[j], abs(rr[j])])

top = pd.DataFrame(top_list, columns=["SNF_gene","New_gene","r","p","abs_r"])
top = top.sort_values(["SNF_gene","abs_r"], ascending=[True, False])
top.to_csv(f"{odir}/{outp}.SNF_top5_newgenes.tsv", sep="\t", index=False)

uniq_new = top["New_gene"].unique().tolist()
mat = pd.DataFrame(index=snf_ids, columns=uniq_new, dtype=float)
for _, row in top.iterrows():
    mat.loc[row["SNF_gene"], row["New_gene"]] = row["r"]
mat.to_csv(f"{odir}/{outp}.SNF_top5_matrix.tsv", sep="\t", na_rep="NA")

with open(f"{odir}/{outp}.summary.txt","w") as f:
    f.write(f"Samples: {n}\nSNF genes used: {len(snf_ids)}\nNew genes used: {len(new_ids)}\nAll pairs: {out_long.shape[0]}\nTop5 rows: {top.shape[0]}\n")

print("[OK] written:")
print(f"  {odir}/{outp}.SNF_vs_New.all_pairs.tsv")
print(f"  {odir}/{outp}.SNF_top5_newgenes.tsv")
print(f"  {odir}/{outp}.SNF_top5_matrix.tsv")
print(f"  {odir}/{outp}.summary.txt")
PY

# =========================
# 2) R: plot (write temp R script)
# =========================
RSCRIPT="${ODIR}/plot_top5.R"
cat > "$RSCRIPT" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
outp <- args[1]
odir <- paste0(outp, ".snf_new_corr")

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
})

top_fn <- file.path(odir, paste0(outp, ".SNF_top5_newgenes.tsv"))
mat_fn <- file.path(odir, paste0(outp, ".SNF_top5_matrix.tsv"))

if (!file.exists(top_fn)) stop(paste("Missing:", top_fn))
if (!file.exists(mat_fn)) stop(paste("Missing:", mat_fn))

top <- read_tsv(top_fn, show_col_types = FALSE)
mat <- read_tsv(mat_fn, show_col_types = FALSE)

# Heatmap data
mat_long <- mat %>%
  pivot_longer(-1, names_to="New_gene", values_to="r") %>%
  rename(SNF_gene = 1) %>%
  filter(!is.na(r))

p1 <- ggplot(mat_long, aes(x=New_gene, y=SNF_gene, fill=r)) +
  geom_tile(color="white", linewidth=0.2) +
  scale_fill_gradient2(low="#2b6cb0", mid="white", high="#c53030", midpoint=0) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=45, hjust=1, vjust=1),
    panel.grid = element_blank()
  ) +
  labs(x="Top New genes (union)", y="SNF genes", fill="Pearson r")

h <- max(3, 0.35*length(unique(mat_long$SNF_gene)))
ggsave(file.path(odir, paste0(outp, ".SNF_top5_heatmap.pdf")), p1, width=10, height=h)
ggsave(file.path(odir, paste0(outp, ".SNF_top5_heatmap.png")), p1, width=10, height=h, dpi=300)

# Dotplot
top2 <- top %>%
  group_by(SNF_gene) %>%
  mutate(rank = row_number()) %>%
  ungroup()

p2 <- ggplot(top2, aes(x=rank, y=SNF_gene, size=abs_r, color=r)) +
  geom_point() +
  scale_color_gradient2(low="#2b6cb0", mid="grey70", high="#c53030", midpoint=0) +
  scale_x_continuous(breaks=1:5, labels=paste0("Top",1:5)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x=NULL, y="SNF genes", size="|r|", color="r")

h2 <- max(3, 0.35*length(unique(top2$SNF_gene)))
ggsave(file.path(odir, paste0(outp, ".SNF_top5_dotplot.pdf")), p2, width=7, height=h2)
ggsave(file.path(odir, paste0(outp, ".SNF_top5_dotplot.png")), p2, width=7, height=h2, dpi=300)
RS

Rscript "$RSCRIPT" "$OUT"

echo "[DONE] outputs in: $ODIR/"
```

</details>

<details open>
<summary><code>spearman.R</code> · R · 187 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/07_transcriptome/01_tissue_atlas/spearman.R)

```r
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

options(stringsAsFactors = FALSE)

# ===============================
# 1. 输入文件
# ===============================
gifu_file <- "Gifu_all_samples_mean.tsv"
mg20_file <- "MG20_all_samples_mean.tsv"
out_prefix <- "Gifu_MG20"

# ===============================
# 2. 读入数据
# ===============================
gifu <- read_tsv(gifu_file, show_col_types = FALSE)
mg20 <- read_tsv(mg20_file, show_col_types = FALSE)

# ===============================
# 3. 转为 long format
#    并去掉 Gifu_root_dpi02
# ===============================
gifu_long <- gifu %>%
  pivot_longer(
    -gene_id,
    names_to = "sample",
    values_to = "expr_gifu"
  ) %>%
  filter(sample != "Gifu_root_dpi02") %>%
  mutate(
    sample = str_remove(sample, "^Gifu_")
  )

mg20_long <- mg20 %>%
  pivot_longer(
    -gene_id,
    names_to = "sample",
    values_to = "expr_mg20"
  ) %>%
  mutate(
    sample = str_remove(sample, "^MG20_")
  )

# ===============================
# 4. 合并 Gifu & MG20
# ===============================
expr_merge_raw <- gifu_long %>%
  inner_join(mg20_long, by = c("gene_id", "sample"))

# ===============================
# 5. 不添加器官分类：
#    按每个 sample 计算表达相似性
# ===============================
sample_corr <- expr_merge_raw %>%
  group_by(sample) %>%
  filter(
    sum(!is.na(expr_gifu) & !is.na(expr_mg20)) > 1,
    sd(expr_gifu, na.rm = TRUE) > 0,
    sd(expr_mg20, na.rm = TRUE) > 0
  ) %>%
  summarise(
    spearman_r = cor(
      expr_gifu,
      expr_mg20,
      method = "spearman",
      use = "pairwise.complete.obs"
    ),
    n_gene = n_distinct(gene_id),
    .groups = "drop"
  )

write_tsv(
  sample_corr,
  paste0(out_prefix, ".sample_expression_spearman.tsv")
)

# sample 排序
sample_corr <- sample_corr %>%
  mutate(
    sample = factor(sample, levels = sample)
  )

p_sample_corr <- ggplot(sample_corr, aes(x = sample, y = spearman_r)) +
  geom_col(fill = "#59A14F", width = 0.7) +
  geom_text(
    aes(label = sprintf("%.2f", spearman_r)),
    vjust = -0.5,
    size = 3.2
  ) +
  ylim(0, 1) +
  labs(
    x = NULL,
    y = "Spearman correlation (Gifu vs MG20)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  )

ggsave(
  paste0(out_prefix, ".sample_expression_correlation.pdf"),
  p_sample_corr,
  width = 9,
  height = 4.5
)

# ===============================
# 6. 器官分类函数
# ===============================
make_sample_info <- function(samples) {
  tibble(
    sample = samples,
    organ = case_when(
      str_detect(sample, "nodule_dpi10|nodule_dpi21|nodule_dpi45") ~ "nodule",
      str_detect(sample, "root_dpi10|root_dpi21|root_dpi45|root_hpi48") ~ "root_infect",
      str_detect(sample, "root_uni") ~ "root",
      str_detect(sample, "leaf") ~ "leaf",
      str_detect(sample, "stem") ~ "stem",
      str_detect(sample, "flower") ~ "flower",
      str_detect(sample, "pod") ~ "pod",
      TRUE ~ NA_character_
    )
  )
}

sample_info <- make_sample_info(
  intersect(gifu_long$sample, mg20_long$sample)
)

# ===============================
# 7. 添加器官分类后：
#    按 organ 计算表达相似性
# ===============================
expr_merge <- expr_merge_raw %>%
  left_join(sample_info, by = "sample") %>%
  filter(!is.na(organ))

organ_corr <- expr_merge %>%
  group_by(organ) %>%
  filter(
    sum(!is.na(expr_gifu) & !is.na(expr_mg20)) > 1,
    sd(expr_gifu, na.rm = TRUE) > 0,
    sd(expr_mg20, na.rm = TRUE) > 0
  ) %>%
  summarise(
    spearman_r = cor(
      expr_gifu,
      expr_mg20,
      method = "spearman",
      use = "pairwise.complete.obs"
    ),
    n_gene = n_distinct(gene_id),
    .groups = "drop"
  )

write_tsv(
  organ_corr,
  paste0(out_prefix, ".organ_expression_spearman.tsv")
)

p_corr <- ggplot(organ_corr, aes(x = organ, y = spearman_r)) +
  geom_col(fill = "#4E79A7", width = 0.7) +
  geom_text(
    aes(label = sprintf("%.2f", spearman_r)),
    vjust = -0.5,
    size = 3.5
  ) +
  ylim(0, 1) +
  labs(
    x = NULL,
    y = "Spearman correlation (Gifu vs MG20)"
  ) +
  theme_classic(base_size = 12)

ggsave(
  paste0(out_prefix, ".organ_expression_correlation.pdf"),
  p_corr,
  width = 7,
  height = 4
)
```

</details>
