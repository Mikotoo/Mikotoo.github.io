---
title: 'Contig ordering'
sidebar:
  label: 'Contig ordering'
  order: 6
---


## 中文
**作用**：用 nucmer 把组装 contig 比对到参考组装，推导每个 query contig 的主导参考 contig、比对长度、比例与链向，用于排序与定向。
**入口**：`05.1.nucmer.sh` → `05.2.contigOrder.py`
**输入**：`Gifu_ref.fa`（参考）与 `Gifu_contig.fa`（组装 contig），文件名硬编码并在工作目录解析。
**输出**：`Gifu_scaffold.delta`、`Gifu_scaffold.filter.delta`、`Gifu_scaffold.filter.delta.coords`、`Gifu_scaffold.delta.coords`、一张 `mummerplot` PNG；`mummer_analysis_results/dominant_mapping.csv` 与 `.../query_contig_on_dominant_ref_colored.pdf`（`+` 蓝、`-` 红）。
**运行**：`bash 05.1.nucmer.sh`，然后 `python 05.2.contigOrder.py`
**工具**：MUMmer `nucmer --prefix Gifu_scaffold`、`delta-filter -i 89 -l 1000 -1`、`show-coords -THrd`、`mummerplot --layout --png --large`（需 gnuplot）；Python 3 + pandas/matplotlib（`read_csv(sep='\t', skiprows=5)`）。

<details>
<summary>English</summary>

**Purpose**: Align the assembled contigs to a reference with nucmer and derive each query contig's dominant reference contig, aligned length, ratio and strand, for ordering and orientation.
**Entry point**: `05.1.nucmer.sh` → `05.2.contigOrder.py`
**Inputs**: `Gifu_ref.fa` (reference) and `Gifu_contig.fa` (contigs) — hardcoded filenames resolved in the working directory.
**Outputs**: `Gifu_scaffold.delta`, `Gifu_scaffold.filter.delta`, `Gifu_scaffold.filter.delta.coords`, `Gifu_scaffold.delta.coords`, a `mummerplot` PNG; `mummer_analysis_results/dominant_mapping.csv` and `.../query_contig_on_dominant_ref_colored.pdf` (`+` blue, `-` red).
**Run**: `bash 05.1.nucmer.sh`, then `python 05.2.contigOrder.py`
**Tools**: MUMmer `nucmer --prefix Gifu_scaffold`, `delta-filter -i 89 -l 1000 -1`, `show-coords -THrd`, `mummerplot --layout --png --large` (needs gnuplot); Python 3 with pandas and matplotlib (`read_csv(sep='\t', skiprows=5)`).

</details>

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>06.1.nucmer.sh</code> · Shell · 9 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/01_assembly/06_contig_ordering/06.1.nucmer.sh)

```bash


nucmer --prefix Gifu_scaffold Gifu_ref.fa Gifu_contig.fa
delta-filter -i 89 -l 1000 -1 Gifu_scaffold.delta > Gifu_scaffold.filter.delta
show-coords -THrd Gifu_scaffold.filter.delta > Gifu_scaffold.filter.delta.coords

show-coords -THrd Gifu_scaffold.delta > Gifu_scaffold.delta.coords

mummerplot Gifu_scaffold.filter.delta -R Gifu_ref.fa -Q Gifu_contig.fa --layout --png --large
```

</details>

<details open>
<summary><code>06.2.contigOrder.py</code> · Python · 105 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/01_assembly/06_contig_ordering/06.2.contigOrder.py)

```python
import pandas as pd
import matplotlib.pyplot as plt
import os
import seaborn as sns

# 输入输出设置
coords_path = "Gifu_scaffold.filter.delta.coords"
output_dir = "mummer_analysis_results"
os.makedirs(output_dir, exist_ok=True)

# 读取 coords 文件
df = pd.read_csv(coords_path, sep='\t', skiprows=5, header=None)
df.columns = [
    "ref_start", "ref_end", "qry_start", "qry_end",
    "ref_len", "qry_len", "percent_id", "ref_cov", "qry_cov",
    "ref_contig", "qry_contig"
]

# 方向、长度等辅助列
df["strand"] = df.apply(lambda x: "+" if x["qry_start"] <= x["qry_end"] else "-", axis=1)
df["alignment_len"] = abs(df["qry_end"] - df["qry_start"]) + 1
df["ref_start_clean"] = df[["ref_start", "ref_end"]].min(axis=1)
df["ref_end_clean"] = df[["ref_start", "ref_end"]].max(axis=1)

# 汇总每条query对ref的比对长度与方向
summary_records = []
for qry, subdf in df.groupby("qry_contig"):
    total_qry_len = subdf["qry_len"].iloc[0]
    ref_stats = subdf.groupby("ref_contig").agg({
        "alignment_len": "sum",
        "strand": lambda x: list(x)
    }).reset_index()

    for _, row in ref_stats.iterrows():
        ref = row["ref_contig"]
        aligned_len = row["alignment_len"]
        strands = row["strand"]
        ratio = aligned_len / total_qry_len
        plus = strands.count("+")
        minus = strands.count("-")
        dominant_strand = "+" if plus >= minus else "-"
        summary_records.append([qry, ref, aligned_len, ratio, dominant_strand])

summary_df = pd.DataFrame(summary_records, columns=[
    "Query_Contig", "Ref_Contig", "Aligned_Length", "Ratio", "Dominant_Strand"
])

# 找出每个 query 的 dominant ref
dominant_mapping = (
    summary_df.groupby("Query_Contig")
    .apply(lambda x: x.loc[x["Aligned_Length"].idxmax(), ["Query_Contig", "Ref_Contig", "Aligned_Length", "Dominant_Strand"]])
    .reset_index(drop=True)
)
dominant_mapping["strand"] = dominant_mapping["Dominant_Strand"]
dominant_mapping.to_csv(os.path.join(output_dir, "dominant_mapping.csv"), index=False)

# 按 dominant ref 分组，绘制每个 query 的合并线段
plot_data = []
for _, row in dominant_mapping.iterrows():
    qry = row["Query_Contig"]
    ref = row["Ref_Contig"]
    strand = row["strand"]

    sub = df[(df["qry_contig"] == qry) & (df["ref_contig"] == ref)]
    if sub.empty:
        continue
    start = sub["ref_start_clean"].min()
    end = sub["ref_end_clean"].max()
    plot_data.append([ref, qry, start, end, strand])

plot_df = pd.DataFrame(plot_data, columns=["ref_contig", "qry_contig", "start", "end", "strand"])

# 为每个 ref-query 分配独立的 y 坐标
plot_df["ref_qry"] = plot_df["ref_contig"] + "|" + plot_df["qry_contig"]
y_map = {}
y_val = 0
for ref in df["ref_contig"].unique():
    ref_df = dominant_mapping[dominant_mapping["Ref_Contig"] == ref]
    for qry in ref_df["Query_Contig"]:
        key = f"{ref}|{qry}"
        y_map[key] = y_val
        y_val += 1
plot_df["y"] = plot_df["ref_qry"].map(y_map)

# ytick 按照 ref 显示
yticks = plot_df[["y", "ref_contig"]].drop_duplicates().sort_values("y")

# 画图
fig, ax = plt.subplots(figsize=(14, 0.5 * len(y_map) + 2))
for _, row in plot_df.iterrows():
    color = "tab:blue" if row["strand"] == "+" else "tab:red"
    ax.hlines(y=row["y"], xmin=row["start"], xmax=row["end"], color=color, linewidth=2)
    ax.text((row["start"] + row["end"]) / 2, row["y"] + 0.3, row["qry_contig"], fontsize=6, ha="center", rotation=45)

ax.set_yticks(yticks["y"])
ax.set_yticklabels(yticks["ref_contig"])
ax.set_xlabel("Reference Genome Position")
ax.set_ylabel("Reference Contigs")
ax.set_title("Query Contig Mapping on Dominant Reference Contigs")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "query_contig_on_dominant_ref_colored.pdf"))
plt.close()

print("✔ dominant_mapping.csv 已保存")
print("✔ 主比对关系图 query_contig_on_dominant_ref_colored.pdf 已生成")
```

</details>
