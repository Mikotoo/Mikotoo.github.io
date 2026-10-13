---
title: '单核分析脚本'
sidebar:
  label: '单核分析脚本'
  order: 999
---


## 中文
**作用**：`08_single_nucleus` 流程的 Python 与 R 实现脚本，由上一级目录的 7 个 `*.sh` 按相对路径调用。
**入口**：无独立入口，由上一级入口脚本驱动（`python ./Script/<name>.py`、`Rscript ./Script/<name>.R`）。
**输入**：`../02_processData/` 下的 `*_filtered.h5ad`、`data_clusterd_all_gene_{n}.h5ad`、`zz_00_pea.rds`、`All_species_merged_by_GifuT2T.SNF_Symbol.tsv`；`../04_Anno/` 下的 marker CSV 与 Gifu T2T GTF。
**输出**：`../02_processData/` 的 h5ad/rds/CSV、`../03_figures/` 的图、`../04_Anno/` 的 Seurat 对象与 marker 表，以及每个基因目录下的 `01_runs/`、`02_summaries/`、`03_plots/`。
**运行**：`bash ../run_singlecell.sh`、`bash ../05_cluster_DEGs_finder.sh` 等；无入口脚本的成员手动运行，目录名 `Script` 不可更改。
**工具**：scanpy、anndata、scvi-tools、leidenalg、rpy2（Python 3.9.19）；Seurat、scTenifoldKnk、DESeq2、clustree（R 4.4.1）。

<details>
<summary>English</summary>

**Purpose**: The Python and R implementation scripts of the `08_single_nucleus` workflow, called by relative path from the seven `*.sh` entry scripts in the parent directory.
**Entry point**: none of their own; they are driven by the parent entry scripts (`python ./Script/<name>.py`, `Rscript ./Script/<name>.R`).
**Inputs**: `*_filtered.h5ad`, `data_clusterd_all_gene_{n}.h5ad`, `zz_00_pea.rds` and `All_species_merged_by_GifuT2T.SNF_Symbol.tsv` under `../02_processData/`; the marker CSVs and the Gifu T2T GTF under `../04_Anno/`.
**Outputs**: h5ad/rds/CSV files in `../02_processData/`, plot files in `../03_figures/`, Seurat objects and marker tables in `../04_Anno/`, and `01_runs/`, `02_summaries/`, `03_plots/` inside each gene directory.
**Run**: `bash ../run_singlecell.sh`, `bash ../05_cluster_DEGs_finder.sh`, and so on; the members without an entry script are run by hand, and the `Script` directory name must not change.
**Tools**: scanpy, anndata, scvi-tools, leidenalg, rpy2 (Python 3.9.19); Seurat, scTenifoldKnk, DESeq2, clustree (R 4.4.1).

</details>

## 代码

以下为该步骤目录中的实际文件（同步快照，与该目录内容一致）。

<details open>
<summary><code>01_annaDataQC.py</code> · Python · 40 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/08_single_nucleus/Script/01_annaDataQC.py)

```python
import sys
import anndata
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

Dpi00_1 =anndata.read_h5ad("{PROJ_LOTUS_HX}/07_singleCell/02_processData/Dpi00_1_filtered.h5ad")
Dpi00_2 =anndata.read_h5ad("{PROJ_LOTUS_HX}/07_singleCell/02_processData/Dpi00_2_filtered.h5ad")
Dpi02_1 =anndata.read_h5ad("{PROJ_LOTUS_HX}/07_singleCell/02_processData/Dpi02_1_filtered.h5ad")
Dpi02_2 =anndata.read_h5ad("{PROJ_LOTUS_HX}/07_singleCell/02_processData/Dpi02_2_filtered.h5ad")

Dpi00 = Dpi00_1.concatenate(Dpi00_2)
Dpi02 = Dpi02_1.concatenate(Dpi02_2)

samples=["Dpi00", "Dpi02"]
fig, aex = plt.subplots(ncols=2, nrows=2, figsize=(50, 30))
for idx, data in enumerate([Dpi00, Dpi02]):
    data_concatenated = data
    data.write_h5ad(f"02_processData/{samples[idx]}_concatenated.h5ad")
    data.layers['counts'] = data.X.copy()
    sc.pp.normalize_total(data,target_sum=1e4,inplace=True)
    sc.pp.log1p(data)
    data_scaled = data.copy()
    data.write_h5ad(f"02_processData/{samples[idx]}_scaled.h5ad")
    sns.histplot(data_concatenated.obs["total_counts"], bins=100,kde=True,ax=aex[idx,0])
    aex[idx,0].set_title(f"{samples[idx]}_Total counts")
    sns.histplot(data_scaled.X.sum(1), bins=100,kde=True,ax=aex[idx,1])
    aex[idx,1].set_title(f"{samples[idx]}_Shifted logarithm")
plt.savefig("03_figures/hist_counts.png")

adata_list =[]
for sample in samples:
    ad = anndata.read_h5ad(f"02_processData/{sample}_scaled.h5ad")
    ad.obs["Sample"] = sample
    adata_list.append(ad)

data_all = anndata.concat(adata_list)
data_all.write_h5ad("02_processData/data_all.h5ad")
```

</details>

<details open>
<summary><code>02_integrate.py</code> · Python · 56 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/08_single_nucleus/Script/02_integrate.py)

```python
import sys
import anndata
import scanpy as sc
import scvi

# 读取数据
data_scaled = anndata.read_h5ad("./02_processData/data_all.h5ad")

# 循环运行不同的 n_top_genes
for n in [1000, 2000, 3000, 4000, 5000]:

    print(f"==== Running n_top_genes = {n} ====")

    # 每次都复制数据，避免被上一次处理修改
    data = data_scaled.copy()

    # 选择高变基因
    sc.experimental.pp.highly_variable_genes(
        data, 
        flavor="pearson_residuals",
        layer="counts",
        batch_key="Sample",
        n_top_genes=n,
        subset=True,
        inplace=True,
    )

    # 设置 scVI anndata
    scvi.model.SCVI.setup_anndata(
        data,
        layer="counts",
        batch_key="Sample",
        categorical_covariate_keys=["Sample"],
        continuous_covariate_keys=["total_counts"]
    )

    # 固定随机种子
    scvi.settings.seed = 7

    # 初始化模型
    model = scvi.model.SCVI(data)

    # 训练模型
    model.train()

    # 保存 latent space
    data.obsm['X_scVI'] = model.get_latent_representation()

    # 保存 normalized expression
    data.layers['scvi_normalized'] = model.get_normalized_expression(library_size=1e4)

    # 输出文件名
    out_file = f"02_processData/data_integrated_{n}.h5ad"
    data.write_h5ad(out_file)

    print(f"==== Saved: {out_file} ====\n")
```

</details>

<details open>
<summary><code>03_1-cluster_statitcs.py</code> · Python · 123 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/08_single_nucleus/Script/03_1-cluster_statitcs.py)

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# =============================
# 1. 读取命令行参数
# =============================
if len(sys.argv) < 2:
    print("请提供 n 值，例如：python 03_1-cluster_statitcs.py 4000")
    sys.exit(1)

n = sys.argv[1]

# =============================
# 2. 读取 .h5ad 数据
# =============================
h5_file = f"02_processData/data_clusterd_{n}.h5ad"
adata = anndata.read_h5ad(h5_file)

df = adata.obs.copy()

# =============================
# 3. 提取 cluster & Sample
# =============================
df["cluster"] = df["leiden_0.4"].astype(str)
df["sample"] = df["Sample"].astype(str)

# =============================
# 4. cluster 指定顺序
# =============================
cluster_order = ['9', '3', '4', '7', '2', '6', '8', '0', '1', '5']
df["cluster"] = pd.Categorical(df["cluster"], categories=cluster_order, ordered=True)

# =============================
# 5. 每个 cluster 的总细胞数
# =============================
cluster_counts = df.groupby("cluster").size().reset_index(name="total_cells")

# =============================
# 6. 各 Sample 在每 cluster 的真实细胞数（raw counts）
# =============================
cluster_sample_count = df.groupby(["cluster", "sample"]).size().reset_index(name="count")

cluster_sample_count = cluster_sample_count.merge(cluster_counts, on="cluster", how="left")

# =============================
# 7. pivot 创建 raw count 堆叠矩阵
# =============================
plot_df = cluster_sample_count.pivot_table(
    index="cluster", columns="sample", values="count", fill_value=0
)

plot_df = plot_df.reindex(cluster_order)

# =============================
# 8. 绘图（横向 raw count stacked bar）
# =============================
fig, ax = plt.subplots(figsize=(10, 8), dpi=300)

# 自定义颜色
palette = ["#1F77B4", "#FF7F0E"]
samples = plot_df.columns.tolist()
colors = (palette * (len(samples) // len(palette) + 1))[:len(samples)]

# 横向堆叠条形图（使用真实数量）
plot_df.plot(
    kind="barh",
    stacked=True,
    ax=ax,
    color=colors,
    width=0.7
)

# =============================
# 9. 添加 cluster 总细胞数注释（右侧）
# =============================
for idx, row in cluster_counts.iterrows():
    cluster = row["cluster"]
    total_cells = row["total_cells"]
    y_pos = cluster_order.index(cluster)

    ax.text(
        plot_df.sum(axis=1).max() * 1.02,   # 在右侧一点点
        y_pos,
        str(total_cells),
        va="center",
        ha="left",
        fontsize=9,
        fontweight="bold"
    )

# =============================
# 10. 删除外边框，仅保留坐标轴
# =============================
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

ax.tick_params(axis='both', which='both', length=4)

# =============================
# 图形标签
# =============================
ax.set_title(f"Cluster Composition (raw cell counts, n_top_genes = {n})")
ax.set_xlabel("Number of Cells")
ax.set_ylabel("Cluster")
ax.legend(samples, title="Sample", bbox_to_anchor=(1.05, 1), loc="upper left")

plt.tight_layout()

# =============================
# 11. 保存 PDF
# =============================
os.makedirs("03_figures", exist_ok=True)
out_file = f"03_figures/cluster_stacked_raw_{n}.pdf"
plt.savefig(out_file, bbox_inches="tight")

print("完成绘图:", out_file)
```

</details>

<details open>
<summary><code>03_1-dimenRedu_adjusted.py</code> · Python · 122 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/08_single_nucleus/Script/03_1-dimenRedu_adjusted.py)

```python
import sys
sys.path.append('{CONDA_SINGLECELL}/lib/python3.9/site-packages/Scripts')
import anndata
import pandas as pd
import numpy as np
import scanpy as sc
import scipy.sparse as sp
import scvi
import Scripts
import matplotlib.pyplot as plt
import seaborn as sns
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

pandas2ri.activate()
clustree = importr('clustree')
ggplot2 = importr('ggplot2')

# 循环运行不同的 n_top_genes
for n in [1000, 2000, 3000, 4000, 5000]:

    print(f"==== Running n_top_genes = {n} ====")

    in_file = f"02_processData/data_integrated_{n}.h5ad"

    data_integrated = anndata.read_h5ad(in_file)
    data = data_integrated

    palette = ["#1F77B4","#FF7F0E"]


    categories = data.obs['Sample'].cat.categories
    color_map = {category: color for category, color in zip(categories, palette)}

    ###降维
    ## UMAP
    # 通过降低n_neighbors和min_dist让各类之间的边界更清晰
    sc.pp.neighbors(data, n_neighbors=10, use_rep='X_scVI')
    sc.tl.umap(data, min_dist=0.3)
    data.obsm["scvi_umap"]=data.obsm["X_umap"]

    fig,ax = plt.subplots(figsize=(10,6),dpi=300)
    sc.pl.embedding(data,basis='scvi_umap',color="Sample",palette=color_map,show=False, ax=ax)
    plt.subplots_adjust(left=0.1, right=0.8, top=0.9, bottom=0.1)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, title='Sample', bbox_to_anchor=(1.05, 1), loc='upper left')
    out_file = f"03_figures/dimenRedu_{n}_adjusted.png"
    plt.savefig(out_file, bbox_inches='tight')
    plt.close()

    ##聚类
    resolutions = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]

    for res in resolutions:
        sc.tl.leiden(data, resolution=res, key_added=f'leiden_{res}',random_state=0)


    fig,axes=plt.subplots(3,2,figsize=(30,15))
    plt.subplots_adjust(left=0.1,bottom=0.1,top=0.9,right=0.8,hspace=0.2,wspace=0.4)
    sc.pl.embedding(data,basis='scvi_umap',ax=axes[0,0],color="leiden_0.2")
    sc.pl.embedding(data,basis='scvi_umap',ax=axes[0,1],color="leiden_0.4")
    sc.pl.embedding(data,basis='scvi_umap',ax=axes[1,0],color="leiden_0.6")
    sc.pl.embedding(data,basis='scvi_umap',ax=axes[1,1],color="leiden_0.8")
    sc.pl.embedding(data,basis='scvi_umap',ax=axes[2,0],color="leiden_1.0")
    sc.pl.embedding(data,basis='scvi_umap',ax=axes[2,1],color="leiden_1.2")
    out_file = f"03_figures/cluster_{n}_adjusted.png"
    plt.savefig(out_file)
    plt.close()

    cluster_df = pd.DataFrame({f'leiden_{res}': data.obs[f'leiden_{res}'] for res in resolutions})
    r_cluster_df = pandas2ri.py2rpy(cluster_df)
    ro.globalenv['cluster_df'] = r_cluster_df
    ro.globalenv['n'] = n
    r_code = """
    library(clustree)
    library(ggplot2)
    p <- clustree(cluster_df, prefix = "leiden_")
    out_file = paste0("03_figures/clustree_plot_", n, ".pdf")
    ggsave(out_file, plot = p, device = "pdf", width=16, height=8)
    """
    ro.r(r_code)

    # final
    fig,axes=plt.subplots(2,1,figsize=(8,10),dpi=300)
    plt.subplots_adjust(left=0.1,bottom=0.1,top=0.9,right=0.8,hspace=0.2,wspace=0.25)
    sc.pl.embedding(data,basis='scvi_umap',ax=axes[0],color="Sample",palette=color_map,show=False)
    sc.pl.embedding(data,basis='scvi_umap',ax=axes[1],color="leiden_0.4",legend_loc="on data",show=False)
    axes[0].set_title("Samples")
    axes[1].set_title("Clusters")
    out_file = f"03_figures/final_cluster_{n}_adjusted.png"
    plt.savefig(out_file)
    plt.close()

    data.write_h5ad(f"02_processData/data_clusterd_{n}_adjusted.h5ad")

    # cluster_statitcs
    clusters = data.obs["leiden_0.4"]
    samples = data.obs['Sample']
    cluster_counts = clusters.value_counts().sort_index()
    cluster_sample_counts=pd.crosstab(clusters,samples)
    cluster_sample_props = cluster_sample_counts.div(cluster_sample_counts.sum(axis=1), axis=0)


    fig, axs = plt.subplots(nrows=2, figsize=(12, 16), dpi=300)
    cluster_counts.plot(kind='bar', ax=axs[0], color='grey')
    axs[0].set_title('Total Number of Cells in Each Cluster')
    axs[0].set_xlabel('Cluster')
    axs[0].set_ylabel('Number of Cells')

    cluster_sample_props.plot(kind='bar', stacked=True, color=[color_map[color] for color in categories], ax=axs[1])
    axs[1].set_title('Proportion of Samples in Each Cluster')
    axs[1].set_xlabel('Cluster')
    axs[1].set_ylabel('Proportion')

    handles, labels = axs[1].get_legend_handles_labels()
    axs[1].legend(handles, labels, title='Sample', bbox_to_anchor=(1.05, 1), loc='upper left')

    out_file = f"03_figures/cluster_statitcs_{n}_adjusted.png"
    plt.tight_layout()
    plt.savefig(out_file, bbox_inches='tight')
    plt.close()
```

</details>

<details open>
<summary><code>03_2-fliter_to_all_gene_h5ad.py</code> · Python · 20 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/08_single_nucleus/Script/03_2-fliter_to_all_gene_h5ad.py)

```python
import anndata

# 循环运行不同的 n_top_genes
for n in [1000, 2000, 3000, 4000, 5000]:

    data_concatenated = anndata.read_h5ad("./02_processData/data_all.h5ad")

    print(f"==== Running n_top_genes = {n} ====")

    in_file = f"02_processData/data_clusterd_{n}.h5ad"
    data_cluster = anndata.read_h5ad(in_file)

    # Copy 'obs' and 'var' if necessary
    data_concatenated.obs = data_cluster.obs
    data_concatenated.uns = data_cluster.uns
    data_concatenated.obsm = data_cluster.obsm
    #data_concatenated.layers = data_cluster.layers
    data_concatenated.obsp = data_cluster.obsp

    data_concatenated.write_h5ad(f"02_processData/data_clusterd_all_gene_{n}.h5ad")
```

</details>

<details open>
<summary><code>05_1-Seurat_Exp_h5ad_to_rds.py</code> · Python · 37 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/08_single_nucleus/Script/05_1-Seurat_Exp_h5ad_to_rds.py)

```python
import anndata
import os
import matplotlib.pyplot as plt
import scanpy as sc
import scipy.io as sio

os.makedirs("04_Anno", exist_ok=True)

data_cluster = anndata.read_h5ad("./02_processData/data_clusterd_all_gene_4000.h5ad")
selected_data = data_cluster
selected_data.obs_names_make_unique()
#############cell type
cluster2annotation = {
    '0': '0',
    '1': '1',
    '2': '2',
    '3': '3',
    '4': '4',
    '5': '5',
    '6': '6',
    '7': '7',
    '8': '8',
    '9': '9'
}

selected_data.obs['celltype'] = selected_data.obs['leiden_0.4'].map(cluster2annotation).astype('category')
#############################################

cellinfo=selected_data.obs
geneinfo=selected_data.var
mtx2=selected_data.layers['counts'].T
mtx1=selected_data.X.T

cellinfo.to_csv("04_Anno/cluster_01-cellinfo.csv")
geneinfo.to_csv("04_Anno/cluster_02-geneinfo.csv")
sio.mmwrite("04_Anno/cluster_03-sparse_matrix_counts.mtx",mtx2)
sio.mmwrite("04_Anno/cluster_04-sparse_matrix.mtx",mtx1)
```

</details>

<details open>
<summary><code>05_2-Exp_h5ad_to_rds_cluster.R</code> · R · 74 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/08_single_nucleus/Script/05_2-Exp_h5ad_to_rds_cluster.R)

```r
# 加载必要的包
# dplyr: 数据操作
# Seurat: 单细胞数据分析
# patchwork: 用于将多个图形组合成一个布局
# Biobase: 生物信息学分析工具
library(dplyr)
library(Seurat)
library(patchwork)
library(Biobase)

# 读取注释好的 h5ad 文件转换后的数据
# 读取细胞信息（cellinfo）、基因信息（geneinfo）和表达矩阵（counts）
cellinfo <- read.csv("04_Anno/cluster_01-cellinfo.csv", row.names = 1)  # 细胞元数据
geneinfo <- read.csv("04_Anno/cluster_02-geneinfo.csv", row.names = 1)  # 基因元数据
counts <- Matrix::readMM(file = "04_Anno/cluster_03-sparse_matrix_counts.mtx")  # 稀疏表达矩阵

# 为稀疏矩阵添加行名（基因）和列名（细胞）
rownames(counts) <- rownames(geneinfo)
colnames(counts) <- rownames(cellinfo)

# 创建 Seurat 对象
# 使用表达矩阵 counts 和细胞元数据 cellinfo
All.merge <- CreateSeuratObject(counts = counts, meta.data = cellinfo)

# 将 Seurat 对象保存为 pbmc（更常用的命名）
pbmc <- All.merge

# 数据归一化
# NormalizeData 对基因表达数据进行归一化处理
pbmc <- NormalizeData(pbmc)

# 筛选高变基因
# FindVariableFeatures 用于识别在不同细胞间具有高度可变性的基因
pbmc <- FindVariableFeatures(pbmc)

# 计算不同聚类的平均表达量
# group.by = "celltype" 按照列 celltype 分组
dat <- AverageExpression(pbmc, group.by = "celltype")$RNA

# 将计算的平均表达量保存为 CSV 文件
write.csv(dat, file = "04_Anno/cluster_05-AverageExpression_celltype.csv", row.names = TRUE)

### 计算两个cluster的差异表达情况
#Idents(pbmc) <- pbmc@meta.data$celltype
#cluster5.markers <- FindMarkers(pbmc, ident.1 = c(5,6), ident.2 = c(0,3))

saveRDS(pbmc, file = "04_Anno/cluster_00_pea.rds")

################################################################################
library(dplyr)
library(Seurat)
library(patchwork)
library(Biobase)
#读取注释好的 h5ad 文件转换后的数据
seurat_obj <- readRDS("04_Anno/cluster_00_pea.rds")
### 计算cluster的marker gene
Idents(pbmc) <- pbmc@meta.data$leiden_0.4
##所有cluster
all_markers <- FindAllMarkers(
  object = pbmc, 
  min.pct = 0.25, 
  logfc.threshold = 0.25, 
  only.pos = TRUE
)

#写出marker gene
write.csv(all_markers, file = "04_Anno/cluster_06_cluster_markers_gene.csv", row.names = FALSE)

all_markers <- FindAllMarkers(
  object = pbmc
)
all_markers <- all_markers[all_markers$p_val_adj<0.05&abs(all_markers$avg_log2FC)>1,]
#写出marker gene
write.csv(all_markers, file = "04_Anno/cluster_06_cluster_markers_gene1.csv", row.names = FALSE)
```

</details>

<details open>
<summary><code>05_3-markerAnno_Ps_v5_own.py</code> · Python · 103 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/08_single_nucleus/Script/05_3-markerAnno_Ps_v5_own.py)

```python
import sys
import anndata
import pandas as pd
import numpy as np
import scanpy as sc
import scipy.sparse as sp
import scvi
import scripts
import matplotlib.pyplot as plt

data_clustered = anndata.read_h5ad("./02_processData/data_clusterd_all_gene_4000.h5ad")
data=data_clustered

marker_genes = {
    "Cortex": [
        "Lotus_GifuT2T_01G0461100",
        "Lotus_GifuT2T_02G1487100",
    ],
    "Epidermis": [
        "Lotus_GifuT2T_05G2864900"
    ],
    "Meristem": [
        "Lotus_GifuT2T_03G1828600",
        "Lotus_GifuT2T_04G2253600",
        "Lotus_GifuT2T_04G2253800"
    ],
    "Pericycle": [
        "Lotus_GifuT2T_02G1146000",
        "Lotus_GifuT2T_04G2574000"
    ],
    "Phloem": [
        "Lotus_GifuT2T_03G2130600",
        "Lotus_GifuT2T_05G2958500",
    ],
    "Vascular": [
        "Lotus_GifuT2T_01G0541600", #
        "Lotus_GifuT2T_01G0093700",
        "Lotus_GifuT2T_01G0703300",
        "Lotus_GifuT2T_01G0626500"
    ],
    "Xylem": [
        "Lotus_GifuT2T_01G0273300",
        "Lotus_GifuT2T_04G2644400",
    ],
    "Infected cell": [
        "Lotus_GifuT2T_04G2692200",
        "Lotus_GifuT2T_01G0784100",
        "Lotus_GifuT2T_03G1515300",
        "Lotus_GifuT2T_05G3052700",
    ],
    "Nodule primordium": [
        "Lotus_GifuT2T_01G0025100",
        "Lotus_GifuT2T_04G2544700"
    ],
    "LjLb": [
        "Lotus_GifuT2T_05G2767100",
        "Lotus_GifuT2T_05G2767300",
        "Lotus_GifuT2T_05G2792500"
    ]
}

marker_genes_in_data = dict()
for ct, markers in marker_genes.items():
    markers_found = list()
    for marker in markers:
        if marker in data.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found


sc.tl.dendrogram(data,'leiden_0.4',use_rep="X_scVI")

print(marker_genes_in_data)


# 指定聚类顺序，按实际需求调整列表顺序
# cluster_order = ['3', '2', '6','12','10','7','8','0','1','4','5','9','11','13','14']  # 替换为你的实际顺序
# cluster_order = ['8', '6', '1', '0', '5', '4', '3', '2', '7', '9']
cluster_order = ['0', '1', '5', '8', '6', '2', '7', '3', '4', '9']

# 将 leiden_0.4 列转换为有序分类
data.obs['leiden_0.4'] = pd.Categorical(
    data.obs['leiden_0.4'],
    categories=cluster_order,  # 强制指定顺序
    ordered=True
)

# 创建画布
fig, axes = plt.subplots(2, 1, figsize=(16, 12), dpi=300)

# 绘制点图（无需 groups_order 参数）
sc.pl.dotplot(
    data,
    ax=axes[0],
    groupby="leiden_0.4",
    var_names=marker_genes_in_data,
    dendrogram=False,
    standard_scale="var"
)

axes[1].axis("off")

plt.savefig("03_figures/markerAnno_filter_dot_v5_own_0.4.pdf")
```

</details>

<details open>
<summary><code>05_3-marker_cluster_DEGs_own_order.py</code> · Python · 187 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/08_single_nucleus/Script/05_3-marker_cluster_DEGs_own_order.py)

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import font_manager
import os

# --------- 目录路径（你的 Arial 字体存放目录） ------------
font_dir = "{HOME_HX}/.fonts/Arial"
for f in os.listdir(font_dir):
    if f.endswith(".ttf") or f.endswith(".TTF"):
        font_manager.fontManager.addfont(os.path.join(font_dir, f))

# --------- 强制 matplotlib 使用 Arial 字体 ------------
plt.rcParams['font.family'] = 'Arial'
matplotlib.rcParams['pdf.fonttype'] = 42   # 可编辑文本
matplotlib.rcParams['ps.fonttype'] = 42

# ------------ 保证AI中可编辑的Arial字体 ------------
matplotlib.rcParams['pdf.fonttype'] = 42   # 不要转轮廓
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'
# ------------------------------------------------------

import os
import pandas as pd
import anndata
import scanpy as sc
import numpy as np

# ==========================
# 0. Read h5ad
# ==========================
data = anndata.read_h5ad("./02_processData/data_clusterd_all_gene_4000.h5ad")

# ==========================
# 1. Cluster order
# ==========================
# cluster_order = ['9', '3', '4', '7', '2', '6', '8', '0', '1', '5']
cluster_order = ['5', '1', '0', '8', '6', '2', '7', '4', '3', '9']

data.obs["leiden_0.4"] = pd.Categorical(
    data.obs["leiden_0.4"],
    categories=cluster_order,
    ordered=True
)

# ==========================
# 2. Read marker gene CSV
# ==========================
csv_path = "./04_Anno/cluster_06_cluster_markers_gene.csv"
df = pd.read_csv(csv_path)

df["cluster"] = df["cluster"].astype(str)

# ---- 统一 gene 名格式：替换 - 为 _ ----
df["gene"] = df["gene"].str.replace("-", "_")

# ==========================
# 3. special cluster 8 with corrected gene names
# ==========================
special_cluster8_genes_raw = [
    "Lotus_GifuT2T_04G2692200",
    "Lotus-GifuT2T-01G0002200",
    "Lotus-GifuT2T-02G1393700"
]

# 替换为与矩阵一致的格式
special_cluster8_genes = [g.replace("-", "_") for g in special_cluster8_genes_raw]

cluster_top_genes = {}

# ==========================
# 4. Process each cluster
# ==========================
for cl in cluster_order:
    sub = df[df["cluster"] == cl]

    # Cluster8 special handling
    if cl == "8":
        cluster_top_genes[cl] = special_cluster8_genes
        continue

    # Apply filtering rules
    sub_filtered = sub[
        (sub["avg_log2FC"] > 0.25) &
        (sub["pct.1"] > 0.25) &
        (sub["p_val_adj"] < 0.05)
    ]

    if sub_filtered.shape[0] == 0:
        print(f"⚠️ Cluster {cl} has no genes passing filter")
        continue

    # Sort by p_val_adj → avg_log2FC
    sub_filtered = sub_filtered.sort_values(
        by=["p_val_adj", "avg_log2FC"],
        ascending=[True, False]
    )

    # Select top3
    cluster_top_genes[cl] = list(sub_filtered["gene"].head(3))

# print final used genes
print("\n最终用于绘图的基因：")
for cl in cluster_order:
    if cl in cluster_top_genes:
        print(f"Cluster {cl}: {cluster_top_genes[cl]}")
    else:
        print(f"Cluster {cl}: 无基因符合条件")

# ==============================================================
# 5. DOTPLOT ALL CLUSTERS
# ==============================================================

all_genes = []
for cl in cluster_order:
    if cl in cluster_top_genes:
        all_genes.extend(cluster_top_genes[cl])

# remove genes not in var
all_genes = [g for g in all_genes if g in data.var.index]

# 如果仍然没有基因，直接退出防止报错
if len(all_genes) == 0:
    raise ValueError("❌ dotplot 无基因可绘制，检查基因名格式是否都正确匹配 data.var.index")

# 自动根据基因数调整画布宽度
n_genes = len(all_genes)
fig_width = max(12, n_genes * 0.8)

print(f"绘制 dotplot，基因数 {n_genes}，figsize=({fig_width}, 10)")

plt.figure(figsize=(fig_width, 10), dpi=300)

sc.pl.dotplot(
    data,
    var_names=all_genes,
    groupby="leiden_0.4",
    standard_scale="var",
    dendrogram=False,
    show=False
)

plt.savefig("./04_Anno/dotplot_all_clusters_top3.pdf", dpi=300)
plt.close()

print("\nDOTPLOT Done!")

# ==============================================================
# 6. UMAP for each cluster
# ==============================================================

out_root = "./04_Anno/UMAP_top/"
os.makedirs(out_root, exist_ok=True)

red_gray_cmap = plt.cm.colors.LinearSegmentedColormap.from_list(
    "redgray",
    ["lightgray", "red"]
)

for cl in cluster_order:
    if cl not in cluster_top_genes:
        continue

    cl_dir = os.path.join(out_root, f"cluster_{cl}")
    os.makedirs(cl_dir, exist_ok=True)

    for gene in cluster_top_genes[cl]:
        if gene not in data.var.index:
            print(f"⚠️ {gene} not found in data")
            continue

        sc.pl.umap(
            data,
            color=gene,
            show=False,
            color_map=red_gray_cmap,
            s=10
        )

        plt.title(f"Cluster {cl} - {gene}")
        plt.savefig(os.path.join(cl_dir, f"{gene}_UMAP.pdf"), dpi=300)
        plt.close()

print("\nUMAP All Done!")
```

</details>

<details open>
<summary><code>05_3-marker_knowngene_lotus_own_order.py</code> · Python · 123 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/08_single_nucleus/Script/05_3-marker_knowngene_lotus_own_order.py)

```python
import sys
import anndata
import pandas as pd
import numpy as np
import scanpy as sc
import scipy.sparse as sp
import scripts
import matplotlib.pyplot as plt

data_clustered = anndata.read_h5ad("./02_processData/data_clusterd_all_gene_4000.h5ad")
data=data_clustered


### marker基因注释
marker_genes = {"":[
"Lotus_GifuT2T_06G3179800",
"Lotus_GifuT2T_01G0061000",
"Lotus_GifuT2T_01G0087900",
"Lotus_GifuT2T_01G0128800"]}

# check if the markers are in the data
marker_genes_in_data = dict()
for ct, markers in marker_genes.items():
    markers_found = list()
    for marker in markers:
        if marker in data.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found


sc.tl.dendrogram(data,'leiden_0.4',use_rep="X_scVI")

print(marker_genes_in_data)

# 指定聚类顺序，按实际需求调整列表顺序
cluster_order = ['0', '1', '5', '8', '6', '2', '7', '3', '4', '9'] # 替换为你的实际顺序

# 将 leiden_0.4 列转换为有序分类
data.obs['leiden_0.4'] = pd.Categorical(
    data.obs['leiden_0.4'],
    categories=cluster_order,  # 强制指定顺序
    ordered=True
)

categories = data.obs['Sample'].cat.categories
for sample in categories:
    selected_cells = data.obs['Sample'] == sample
    sample_data = data[selected_cells].copy()

    marker_genes_in_data = dict()
    for ct, markers in marker_genes.items():
        markers_found = list()
        for marker in markers:
            if marker in sample_data.var.index:
                markers_found.append(marker)
        marker_genes_in_data[ct] = markers_found

    sc.tl.dendrogram(sample_data,'leiden_0.4',use_rep="X_scVI")

    fig, axes=plt.subplots(2,1,figsize=(16,12),dpi=300)
    sc.pl.dotplot(
        sample_data,
        ax=axes[0],
        groupby="leiden_0.4",
        var_names=marker_genes_in_data,
        dendrogram=False,
        standard_scale="var",
    )
    axes[1].axis("off")

    plt.savefig(f"03_figures/knowngene_dot_{sample}_1_own_0.4.pdf")

# ================================
#  绘制 marker 基因的 UMAP 分布图
# ================================
# 四个 marker 基因列表
genes_to_plot = [
    "Lotus_GifuT2T_06G3179800",
    "Lotus_GifuT2T_01G0061000",
    "Lotus_GifuT2T_01G0087900",
    "Lotus_GifuT2T_01G0128800"
]

# 逐个 sample 绘制 UMAP
for sample in categories:
    selected_cells = data.obs['Sample'] == sample
    sample_data = data[selected_cells].copy()

    # 筛选当前 sample 中确实存在的基因
    valid_genes = [g for g in genes_to_plot if g in sample_data.var.index]

    if len(valid_genes) == 0:
        print(f"⚠️ Sample {sample} 中没有任何 marker 基因")
        continue

    # ---- 创建一个 2×2 的 panel 图 ----
    n = len(valid_genes)
    ncols = 2
    nrows = int(np.ceil(n / ncols))

    fig, axs = plt.subplots(nrows, ncols, figsize=(10, 10), dpi=300)
    axs = axs.flatten()

    for i, gene in enumerate(valid_genes):
        sc.pl.umap(
            sample_data,
            color=gene,
            ax=axs[i],
            show=False,
            color_map="viridis",  # 可换 plasma, magma
            s=10
        )
        axs[i].set_title(gene, fontsize=10)

    # 隐藏多余子图
    for j in range(i+1, len(axs)):
        axs[j].axis('off')

    plt.tight_layout()
    plt.savefig(f"03_figures/knowngene_umap_{sample}_4genes.pdf")
    plt.close()

print("UMAP 分布图已全部绘制完成！")
```

</details>

<details open>
<summary><code>06_1-Seurat_Exp_h5ad_to_rds.py</code> · Python · 38 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/08_single_nucleus/Script/06_1-Seurat_Exp_h5ad_to_rds.py)

```python
import anndata
import os
import matplotlib.pyplot as plt
import scanpy as sc
import scipy.io as sio

data_cluster = anndata.read_h5ad("./02_processData/data_clusterd_all_gene_4000.h5ad")
data_cluster.obs_names_make_unique()
selected_data = data_cluster

#############cell type
#Xylem:1   Phloem:2   Pericycle:3   Vascular:4   Cortex:5   Meristem:6    Early symbiotic signaling cells:7 
cluster2annotation = {
    '0': '5',
    '1': '5',
    '2': '3',
    '3': '4',
    '4': '4',
    '5': '5',
    '6': '6',
    '7': '2',
    '8': '7',
    '9': '1'
}

selected_data.obs['celltype'] = selected_data.obs['leiden_0.4'].map(cluster2annotation).astype('category')

#############################################

cellinfo=selected_data.obs
geneinfo=selected_data.var
mtx2=selected_data.layers['counts'].T
mtx1=selected_data.X.T

cellinfo.to_csv("./02_processData/zz_01-cellinfo.csv")
geneinfo.to_csv("./02_processData/zz_02-geneinfo.csv")
sio.mmwrite("./02_processData/zz_03-sparse_matrix_counts.mtx",mtx2)
sio.mmwrite("./02_processData/zz_04-sparse_matrix.mtx",mtx1)
```

</details>

<details open>
<summary><code>06_2-Exp_h5ad_to_rds_cell_type.R</code> · R · 146 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/08_single_nucleus/Script/06_2-Exp_h5ad_to_rds_cell_type.R)

```r
# 加载必要的包
# dplyr: 数据操作
# Seurat: 单细胞数据分析
# patchwork: 用于将多个图形组合成一个布局
# Biobase: 生物信息学分析工具
library(dplyr)
library(Seurat)
library(patchwork)
library(Biobase)
library(tidyverse)

# 读取注释好的 h5ad 文件转换后的数据
# 读取细胞信息（cellinfo）、基因信息（geneinfo）和表达矩阵（counts）
cellinfo <- read.csv("./02_processData/zz_01-cellinfo.csv")  # 细胞元数据
geneinfo <- read.csv("./02_processData/zz_02-geneinfo.csv", row.names = 1)  # 基因元数据
counts <- Matrix::readMM(file = "./02_processData/zz_03-sparse_matrix_counts.mtx")  # 稀疏表达矩阵

#------------- Step 1. 读取 barcode -------------
cell00dpi1 <- read_tsv("00dpi_1_barcodes.tsv", col_names = FALSE)$X1
cell00dpi2 <- read_tsv("00dpi_2_barcodes.tsv", col_names = FALSE)$X1
cell02dpi1 <- read_tsv("02dpi_1_barcodes.tsv", col_names = FALSE)$X1
cell02dpi2 <- read_tsv("02dpi_2_barcodes.tsv", col_names = FALSE)$X1

# 合并为两个集合
set00dpi <- unique(c(cell00dpi1, cell00dpi2))
set02dpi <- unique(c(cell02dpi1, cell02dpi2))

#------------- Step 2. 清洗 cellinfo$X（去掉第二个 "-" 及后内容）-------------
# 原格式如：AAACCCACAACAGCCC-1-0
cellinfo$X_clean <- sub("^([^-]+-[^-]+).*", "\\1", cellinfo$X)

#------------- Step 3. 判断属于哪一类并加后缀 -------------
cellinfo$X_final <- case_when(
  cellinfo$X_clean %in% set00dpi ~ paste0(cellinfo$X_clean, "-0"),
  cellinfo$X_clean %in% set02dpi ~ paste0(cellinfo$X_clean, "-1"),
  TRUE ~ cellinfo$X_clean    # 不在任何列表就保持原样（一般不会出现）
)

#------------- Step 4. 处理重复并设置行名 -------------
cellinfo$X_final_unique <- make.unique(cellinfo$X_final, sep = "_")
rownames(cellinfo) <- cellinfo$X_final_unique
cellinfo <- cellinfo %>% select(-X, -X_clean, -X_final, -X_final_unique)

# 为稀疏矩阵添加行名（基因）和列名（细胞）
rownames(counts) <- rownames(geneinfo)
colnames(counts) <- rownames(cellinfo)

# 创建 Seurat 对象
# 使用表达矩阵 counts 和细胞元数据 cellinfo
All.merge <- CreateSeuratObject(counts = counts, meta.data = cellinfo)

# 将 Seurat 对象保存为 pbmc（更常用的命名）
pbmc <- All.merge

# 数据归一化
# NormalizeData 对基因表达数据进行归一化处理
pbmc <- NormalizeData(pbmc)

# 筛选高变基因
# FindVariableFeatures 用于识别在不同细胞间具有高度可变性的基因
pbmc <- FindVariableFeatures(pbmc)

# 计算不同聚类的平均表达量
# group.by = "celltype" 按照列 celltype 分组
dat <- AverageExpression(pbmc, group.by = "celltype")$RNA

# 将计算的平均表达量保存为 CSV 文件
write.csv(dat, file = "./02_processData/zz_05-AverageExpression_celltype.csv", row.names = TRUE)

### 计算两个cluster的差异表达情况
#Idents(pbmc) <- pbmc@meta.data$celltype
#cluster5.markers <- FindMarkers(pbmc, ident.1 = c(5,6), ident.2 = c(0,3))

### 将对象拆分为不同时期处理，比较不同时期的表达差异
# 假设 pbmc 是一个 Seurat 对象，提取元数据
metadata <- pbmc@meta.data

# 提取细胞名称，并获取最后一位数字
cell_names <- rownames(metadata)
last_digits <- as.numeric(substr(cell_names, nchar(cell_names), nchar(cell_names)))

# 动态调整 cell_type 的值 0-9 10-19 20-29 30-39 40-49
metadata$adjusted_celltype <- metadata$celltype + 10 * last_digits

# 将修改后的元数据更新到 Seurat 对象
pbmc@meta.data <- metadata

# 查看调整后的结果
#head(metadata[, c("leiden_0.6", "adjusted_leiden_0.6")])

# 计算基因表达
dat <- AverageExpression(pbmc, group.by = "adjusted_celltype")$RNA
write.csv(dat, file = "./02_processData/zz_06-AverageExpression_celltype_difftime.csv", row.names = TRUE)

saveRDS(pbmc, file = "./02_processData/zz_00_pea.rds")


################################################################################
library(dplyr)
library(Seurat)
library(patchwork)
library(Biobase)
#读取注释好的 h5ad 文件转换后的数据
seurat_obj <- readRDS("./02_processData/zz_00_pea.rds")
### 计算cluster的marker gene
Idents(pbmc) <- pbmc@meta.data$adjusted_celltype
##所有cluster
all_markers <- FindAllMarkers(
  object = pbmc, 
  min.pct = 0.25, 
  logfc.threshold = 0.25, 
  only.pos = TRUE
)
#写出marker gene
write.csv(all_markers, file = "./02_processData/zz_07_celltype_difftime_markers_gene.csv", row.names = FALSE)

### 计算cluster的marker gene
Idents(pbmc) <- pbmc@meta.data$celltype
##所有cluster
all_markers <- FindAllMarkers(
  object = pbmc, 
  min.pct = 0.25, 
  logfc.threshold = 0.25, 
  only.pos = TRUE
)
#写出marker gene
write.csv(all_markers, file = "./02_processData/zz_07_celltype_markers_gene.csv", row.names = FALSE)

################################################################################
library(dplyr)
library(Seurat)
library(patchwork)
library(Biobase)
#读取注释好的 h5ad 文件转换后的数据
seurat_obj <- readRDS("./02_processData/zz_00_pea.rds")
### 计算cluster的marker gene
Idents(pbmc) <- pbmc@meta.data$leiden_0.4
##所有cluster
all_markers <- FindAllMarkers(
  object = pbmc, 
  min.pct = 0.25, 
  logfc.threshold = 0.25, 
  only.pos = TRUE
)
#写出marker gene
write.csv(all_markers, file = "./02_processData/zz_08_cluster_markers_gene.csv", row.names = FALSE)
```

</details>

<details open>
<summary><code>06_3-agg_to_fpkm.R</code> · R · 119 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/08_single_nucleus/Script/06_3-agg_to_fpkm.R)

```r
# 加载必要的包
# dplyr: 数据操作
# Seurat: 单细胞数据分析
# patchwork: 用于将多个图形组合成一个布局
# Biobase: 生物信息学分析工具
library(dplyr)
library(Seurat)
library(patchwork)
library(Biobase)
library(tidyverse)

# 读取注释好的 h5ad 文件转换后的数据
# 读取细胞信息（cellinfo）、基因信息（geneinfo）和表达矩阵（counts）
cellinfo <- read.csv("./02_processData/zz_01-cellinfo.csv")  # 细胞元数据
geneinfo <- read.csv("./02_processData/zz_02-geneinfo.csv", row.names = 1)  # 基因元数据
counts <- Matrix::readMM(file = "./02_processData/zz_03-sparse_matrix_counts.mtx")  # 稀疏表达矩阵

#------------- Step 1. 读取 barcode -------------
cell00dpi1 <- read_tsv("00dpi_1_barcodes.tsv", col_names = FALSE)$X1
cell00dpi2 <- read_tsv("00dpi_2_barcodes.tsv", col_names = FALSE)$X1
cell02dpi1 <- read_tsv("02dpi_1_barcodes.tsv", col_names = FALSE)$X1
cell02dpi2 <- read_tsv("02dpi_2_barcodes.tsv", col_names = FALSE)$X1

# 合并为两个集合
set00dpi <- unique(c(cell00dpi1, cell00dpi2))
set02dpi <- unique(c(cell02dpi1, cell02dpi2))

#------------- Step 2. 清洗 cellinfo$X（去掉第二个 "-" 及后内容）-------------
# 原格式如：AAACCCACAACAGCCC-1-0
cellinfo$X_clean <- sub("^([^-]+-[^-]+).*", "\\1", cellinfo$X)

#------------- Step 3. 判断属于哪一类并加后缀 -------------
cellinfo$X_final <- case_when(
  cellinfo$X_clean %in% set00dpi ~ paste0(cellinfo$X_clean, "-0"),
  cellinfo$X_clean %in% set02dpi ~ paste0(cellinfo$X_clean, "-1"),
  TRUE ~ cellinfo$X_clean    # 不在任何列表就保持原样（一般不会出现）
)

#------------- Step 4. 处理重复并设置行名 -------------
cellinfo$X_final_unique <- make.unique(cellinfo$X_final, sep = "_")
rownames(cellinfo) <- cellinfo$X_final_unique
cellinfo <- cellinfo %>% select(-X, -X_clean, -X_final, -X_final_unique)

# 为稀疏矩阵添加行名（基因）和列名（细胞）
rownames(counts) <- rownames(geneinfo)
colnames(counts) <- rownames(cellinfo)

# 创建 Seurat 对象
# 使用表达矩阵 counts 和细胞元数据 cellinfo
All.merge <- CreateSeuratObject(counts = counts, meta.data = cellinfo)

# 将 Seurat 对象保存为 pbmc（更常用的命名）
pbmc <- All.merge

### 将对象拆分为不同时期处理，比较不同时期的表达差异
# 假设 pbmc 是一个 Seurat 对象，提取元数据
metadata <- pbmc@meta.data

# 提取细胞名称，并获取最后一位数字
cell_names <- rownames(metadata)
last_digits <- as.numeric(substr(cell_names, nchar(cell_names), nchar(cell_names)))

# 动态调整 cell_type 的值 0-8 9-17 18-26 27-35 36-44
# 动态调整 cell_type 的值 0-9 10-19 20-29 30-39 40-49
metadata$adjusted_celltype <- metadata$celltype + 10 * last_digits

# 将修改后的元数据更新到 Seurat 对象
pbmc@meta.data <- metadata


# 计算不同聚类的平均表达量
# group.by = "celltype" 按照列 celltype 分组
# 通过AggregateExpression生成伪样本计数矩阵
pseudobulk <- Seurat::AggregateExpression(
  pbmc,
  group.by = "adjusted_celltype",  # 分组字段
  slot = "counts",                     # 使用原始counts
  fun = "sum"                          # 求和聚合
)
count_matrix <- pseudobulk$RNA          # 提取计数矩阵

# 从伪样本列名解析元数据（假设列名为"sample_celltype"）
sample_info <- data.frame(
  celltype = colnames(count_matrix),  # 直接使用列名
  row.names = colnames(count_matrix)
)

# 从GTF文件加载基因长度（需替换为实际数据）
library(GenomicFeatures)
txdb <- txdbmaker::makeTxDbFromGFF("{PROJ_LOTUS_HX}/00_assembly_result/GifuT2T/Lotus_GifuT2T.gtf")
gene_lengths <- transcriptLengths(txdb, with.cds_len=TRUE) %>%
  group_by(gene_id) %>%
  summarise(length = max(cds_len))  # 取最长转录本长度

# 确保基因名称与计数矩阵一致
rownames(gene_lengths) <- gene_lengths$gene_id
gene_lengths <- gene_lengths[rownames(count_matrix), ]

library(DESeq2)

# 创建DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ celltype  # 根据分析目标设定design
)

# 添加基因长度信息
mcols(dds)$basepairs <- gene_lengths$length

# 方法1：直接计算（未标准化）
fpkm_matrix <- fpkm(dds, robust = TRUE)  # robust=TRUE对长度做中位数校正

# 方法2：先运行标准化（推荐）
dds <- estimateSizeFactors(dds)
fpkm_matrix <- fpkm(dds, robust = TRUE)

# 保存为CSV
write.csv(fpkm_matrix, "./02_processData/00_02_pseudobulk_fpkm.csv")
```

</details>

<details open>
<summary><code>06_4-DEGs_cell_type.R</code> · R · 183 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/08_single_nucleus/Script/06_4-DEGs_cell_type.R)

```r
library(dplyr)
library(Seurat)
library(patchwork)
library(Biobase)
library(ggplot2)

#读取注释好的 h5ad 文件转换后的数据
pbmc <- readRDS("./02_processData/zz_00_pea.rds")
### 计算cluster的marker gene
#Idents(pbmc) <- pbmc@meta.data$adjusted_celltype
##所有cluster
#all_markers <- FindAllMarkers(
#  object = pbmc
#)
#写出marker gene
#write.csv(all_markers, file = "zz_07_celltype_difftime_markers_gene.csv", row.names = FALSE)

### 计算cluster的marker gene
# 定义 cluster 组别
cluster_groups <- list(
  Xylem = c(1, 11),
  Phloem = c(2, 12),
  Pericycle = c(3, 13),
  Vascular = c(4, 14),
  Cortex = c(5, 15),
  Meristem = c(6, 16),
  Early = c(7, 17)
)

# 设置细胞标识
Idents(pbmc) <- pbmc@meta.data$adjusted_celltype

# 初始化结果列表
marker_results <- list()

# 计算组内差异基因，以dpi01为背景
group_names <- names(cluster_groups)
for (i in 1:(length(group_names))) {
    group_name <- names(cluster_groups)[i]
    clusters <- cluster_groups[[group_name]]
    markers <- FindMarkers(
      pbmc,
      ident.1 = clusters[1],
      ident.2 = clusters[2]
      )
    markers <- markers %>%
      mutate(Significance = case_when(
        p_val_adj < 0.05 & avg_log2FC > 1 ~ "Upregulated",
        p_val_adj < 0.05 & avg_log2FC < -1 ~ "Downregulated",
        TRUE ~ "Not Significant"
      ))

      # 保存结果
    comparison_name <- paste(group_name, "02dpivs00dpi", sep = "_")
    marker_results[[comparison_name]] <- markers

    ggplot(markers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
      geom_point(alpha = 0.8) +
      scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
      theme_minimal() +
      labs(
        title = "Volcano Plot",
        x = "Log2 Fold Change",
        y = "-Log10 Adjusted P-Value"
      )

    ggsave(
      filename = paste0("./03_figures/", comparison_name, "_volcano_plots.pdf"),
      width = 6, height = 5)

    # 输出结果到文件
    write.csv(markers, file = paste0("./02_processData/", comparison_name, "_markers.csv"))

    # 打印提示信息
    cat("已完成比较：", comparison_name, "\n")
}

# 初始化结果列表
marker_results <- list()

# 遍历所有组合
group_names <- names(cluster_groups)
for (i in 1:(length(group_names) - 1)) {
  for (j in (i + 1):length(group_names)) {
    # 获取当前组合的组名和 cluster 编号
    group1 <- group_names[i]
    group2 <- group_names[j]
    clusters1 <- cluster_groups[[group1]]
    clusters2 <- cluster_groups[[group2]]

    # 计算差异基因
    markers <- FindMarkers(
      pbmc,
      ident.1 = clusters1,
      ident.2 = clusters2 
    )
    markers <- markers %>%
      mutate(Significance = case_when(
        p_val_adj < 0.05 & avg_log2FC > 1 ~ "Upregulated",
        p_val_adj < 0.05 & avg_log2FC < -1 ~ "Downregulated",
        TRUE ~ "Not Significant"
      ))

    # 保存结果
    comparison_name <- paste(group1, "vs", group2, sep = "_")
    marker_results[[comparison_name]] <- markers

    # 绘制火山图
    ggplot(markers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
      geom_point(alpha = 0.8) +
      scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
      theme_minimal() +
      labs(
        title = "Volcano Plot",
        x = "Log2 Fold Change",
        y = "-Log10 Adjusted P-Value"
      )

    ggsave(
    filename = paste0("./03_figures/", comparison_name, "_volcano_plots.pdf"),
    width = 6, height = 5)

    # 输出结果到文件
    write.csv(markers, file = paste0("./02_processData/", comparison_name, "_markers.csv"))

    # 打印提示信息
    cat("已完成比较：", comparison_name, "\n")
  }
}

###############################################
### 额外分析：Early symbiotic signaling cells vs 非 Early symbiotic signaling cells
###############################################

cat("开始 Early symbiotic signaling cells vs Others 分析...\n")

# 定义 Early symbiotic signaling cells 的 clusters
Early_clusters <- cluster_groups$Early  # 即 c(7, 17)

# 所有 clusters
all_clusters <- unique(unlist(cluster_groups))

# 除 Early symbiotic signaling cells 外的其他 clusters
other_clusters <- setdiff(all_clusters, Early_clusters)

# 执行差异分析（Early symbiotic signaling cells 为 ident.1，上调表示 Early symbiotic signaling cells 相对 others 上调）
markers_inf <- FindMarkers(
  pbmc,
  ident.1 = Early_clusters,
  ident.2 = other_clusters
)

markers_inf <- markers_inf %>%
  mutate(Significance = case_when(
    p_val_adj < 0.05 & avg_log2FC > 1  ~ "Upregulated",
    p_val_adj < 0.05 & avg_log2FC < -1 ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

# 保存结果
comparison_name <- "Early_symbiotic_signaling_cells_vs_Others"
marker_results[[comparison_name]] <- markers_inf

# 绘制火山图
p_inf <- ggplot(markers_inf, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Early symbiotic signaling cells vs Others",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value"
  )

ggsave(
  filename = paste0("./03_figures/", comparison_name, "_volcano_plots.pdf"),
  width = 6, height = 5
)

# 输出到 CSV
write.csv(markers_inf, file = paste0("./02_processData/", comparison_name, "_markers.csv"))

cat("已完成比较：", comparison_name, "\n")
```

</details>

<details>
<summary><code>07_batch_multiple_scTenifoldKnk.R</code> · R · 579 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/08_single_nucleus/Script/07_batch_multiple_scTenifoldKnk.R)

```r
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(ggplot2)
  library(ggrepel)
  library(remotes)
  library(scTenifoldKnk)
  library(glue)
  library(gghalves)
})

options(stringsAsFactors = FALSE)

# ============================================================
# 0) 参数
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("用法：Rscript 10.scTenifoldKnk_replicates.R gKO", call. = FALSE)

gKO <- args[1]
celltype <- args[2]

# 运行次数
N_RUN <- 100

# 子集数量
K_SET <- c(10, 30, 50, 70, 100)

# hit阈值（用于次数统计）
LOG2FC_CUT <- 1
PADJ_CUT <- 0.05

# 细胞选择（与你当前脚本一致）
N_CELLS_USE <- 2000

# HVG数
N_HVG <- 1000

# ============================================================
# 1) 输出目录
# ============================================================
dir.create("01_runs", showWarnings = FALSE)
dir.create("02_summaries", showWarnings = FALSE)
dir.create("03_plots", showWarnings = FALSE)

# ============================================================
# 2) 数据准备：提取count矩阵 + SNF重命名
# ============================================================
obj5 <- readRDS("zz_00_pea.rds")
SNF  <- read.table("All_species_merged_by_GifuT2T.SNF_Symbol.tsv", header = TRUE)

SNF_fix <- SNF %>%
  mutate(GifuT2T_fix = gsub("_", "-", GifuT2T))

Idents(obj5) <- "celltype"
message("celltype table:")
print(table(obj5$celltype))

# 提取该细胞类型的特异表达基因
deg_ct <- FindMarkers(
  object = obj5,
  ident.1 = celltype,
  ident.2 = NULL
)

deg_ct_sig <- deg_ct %>%
  dplyr::filter(
    p_val_adj < 0.05,
    abs(avg_log2FC) > 1
  )

genes_ct_all <- rownames(deg_ct_sig)

# 早期信号响应细胞（与你当前一致：idents="7"）
sc.s <- subset(obj5, idents = celltype)

count <- GetAssayData(sc.s, slot = "counts")

# HVGs
hvg <- VariableFeatures(sc.s)
if (length(hvg) < N_HVG) {
  warning(glue("HVG数量不足({length(hvg)}), 将使用全部HVG。"))
  hvg_use <- hvg
} else {
  hvg_use <- hvg[1:N_HVG]
}

genes_use <- unique(c(hvg_use, gKO, genes_ct_all,SNF_fix$GifuT2T_fix))
cells_use <- colnames(count)[1:min(N_CELLS_USE, ncol(count))]

count_use <- count[genes_use, cells_use]

# --- 将SNF的基因号替换为GeneSymbol（处理同名重复加后缀） ---
genes_count <- rownames(count_use)

snf_map <- SNF_fix %>%
  filter(GifuT2T_fix %in% genes_count) %>%
  select(GeneSymbol, GifuT2T_fix) %>%
  group_by(GeneSymbol) %>%
  mutate(
    index = row_number(),
    NewGene = if (n() > 1) paste0(GeneSymbol, "-", index) else GeneSymbol
  ) %>%
  ungroup()

rename_table <- data.frame(OldGene = genes_count, stringsAsFactors = FALSE) %>%
  left_join(snf_map %>% select(GifuT2T_fix, NewGene),
            by = c("OldGene" = "GifuT2T_fix")) %>%
  mutate(FinalGene = ifelse(is.na(NewGene), OldGene, NewGene))

count_use_renamed <- count_use
rownames(count_use_renamed) <- rename_table$FinalGene

write.table(
  rename_table,
  file = "02_summaries/SNF_gene_rename_map.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

# 虚拟敲除的基因名称要和矩阵保持一致
if (gKO %in% rename_table$OldGene){
  gKO_raw <- gKO
  gKO <- rename_table[ rename_table$OldGene %in% gKO_raw, "FinalGene"]
}

# ============================================================
# 3) scTenifoldKnk重复运行：100次并保留结果
# ============================================================
run_one <- function(seed_i, count_mat, gKO) {
  set.seed(seed_i)

  res <- scTenifoldKnk(
    countMatrix = count_mat,
    qc = TRUE,
    gKO = gKO,
    qc_mtThreshold = 0.1,
    qc_minLSize = 1000,
    nc_lambda = 0,
    nc_nNet = 10,
    nc_nCells = min(500, ncol(count_mat)),
    nc_nComp = 3,
    nc_scaleScores = TRUE,
    nc_symmetric = FALSE,
    nc_q = 0.9,
    td_K = 3,
    td_maxIter = 1000,
    td_maxError = 1e-5,
    td_nDecimal = 3,
    ma_nDim = 2,
    nCores = parallel::detectCores()
  )

  df <- res[["diffRegulation"]] %>% as.data.frame()
  # 统一列名：确保 gene / FC / p.adj 存在
  if (!("gene" %in% colnames(df))) stop("diffRegulation中缺少 gene 列")
  if (!("FC" %in% colnames(df))) stop("diffRegulation中缺少 FC 列")
  if (!("p.adj" %in% colnames(df))) stop("diffRegulation中缺少 p.adj 列")

  df <- df %>%
    mutate(
      log2FC = log2(FC),
      logp.adj = -log10(pmax(p.adj, 1e-300))
    ) %>%
    select(gene, distance, Z, FC, p.value, p.adj, log2FC, logp.adj, everything())

  df
}

# 主循环（带失败保护）
all_runs <- vector("list", N_RUN)
run_status <- tibble(
  run = 1:N_RUN,
  seed = 1000 + 1:N_RUN,
  ok = FALSE,
  msg = NA_character_
)

message(glue("Start {N_RUN} runs ..."))
for (i in 1:N_RUN) {
  seed_i <- 1000 + i
  message(glue("[Run {i}/{N_RUN}] seed={seed_i}"))

  out_rds <- file.path("01_runs", glue("run_{sprintf('%03d', i)}.rds"))
  out_csv <- file.path("01_runs", glue("run_{sprintf('%03d', i)}.csv"))

  df_i <- tryCatch(
    {
      df <- run_one(seed_i, count_use_renamed, gKO)

      saveRDS(df, out_rds)
      write.csv(df, out_csv, row.names = FALSE)

      run_status$ok[i] <- TRUE
      run_status$msg[i] <- "OK"
      df
    },
    error = function(e) {
      run_status$ok[i] <- FALSE
      run_status$msg[i] <- as.character(e$message)
      warning(glue("[Run {i}] FAILED: {e$message}"))
      NULL
    }
  )

  all_runs[[i]] <- df_i
}

all_runs <- lapply(
  all_runs,
  function(df) {
    if (is.null(df)) return(NULL)
    df %>% dplyr::filter(gene != gKO)
  }
)

write.csv(run_status, "02_summaries/run_status.csv", row.names = FALSE)
saveRDS(all_runs, "02_summaries/all_runs_list.rds")

ok_idx <- which(run_status$ok)
if (length(ok_idx) < max(K_SET)) {
  warning(glue(
    "成功运行次数只有 {length(ok_idx)}，不足以覆盖最大K={max(K_SET)}。将只对可用范围内的K做分析。"
  ))
}
# 按顺序取成功的前若干次
all_runs_ok <- all_runs[ok_idx]
N_OK <- length(all_runs_ok)

K_SET2 <- K_SET[K_SET <= N_OK]
if (length(K_SET2) == 0) stop("没有足够的成功运行结果用于后续分析。")

# ============================================================
# 4) (2) 相关性：各K下，K次结果两两FC Spearman相关 + 箱线图 + 均值拟合
# ============================================================
# 将每次结果转为 named vector: gene -> FC
to_fc_vec <- function(df) {
  v <- df$FC
  names(v) <- df$gene
  v
}

# 计算两次run的spearman相关（按gene对齐）
spearman_fc <- function(v1, v2) {
  g <- intersect(names(v1), names(v2))
  if (length(g) < 10) return(NA_real_)
  suppressWarnings(cor(v1[g], v2[g], method = "spearman", use = "pairwise.complete.obs"))
}

cor_records <- list()
mean_records <- list()

for (K in K_SET2) {
  runsK <- all_runs_ok[1:K]
  vecs <- lapply(runsK, to_fc_vec)

  cors <- c()
  if (K >= 2) {
    for (i in 1:(K - 1)) {
      for (j in (i + 1):K) {
        cors <- c(cors, spearman_fc(vecs[[i]], vecs[[j]]))
      }
    }
  }

  cor_df <- tibble(
    K = K,
    pair_id = seq_along(cors),
    rho = cors
  ) %>% filter(!is.na(rho))

  cor_records[[as.character(K)]] <- cor_df

  mean_records[[as.character(K)]] <- tibble(
    K = K,
    mean_rho = mean(cor_df$rho, na.rm = TRUE)
  )
}

cor_all <- bind_rows(cor_records)
mean_all <- bind_rows(mean_records)

write.csv(cor_all, "02_summaries/spearman_pairwise_all.csv", row.names = FALSE)
write.csv(mean_all, "02_summaries/spearman_mean_byK.csv", row.names = FALSE)

cor_all <- cor_all %>%
  mutate(K = factor(K, levels = c(10, 30, 50, 70, 100)))

mean_all <- mean_all %>%
  mutate(K = factor(K, levels = c(10, 30, 50, 70, 100)))

mean_all <- mean_all %>%
  mutate(K_index = as.numeric(K))


## -----------------------------
## 1. 数据准备（再次确认）
## -----------------------------
cor_all <- cor_all %>%
  mutate(K = factor(K, levels = c("10", "30", "50", "70", "100")))

mean_all <- mean_all %>%
  mutate(
    K = factor(K, levels = c("10", "30", "50", "70", "100")),
    K_index = as.numeric(K)   # 仅用于趋势线定位
  )

## -----------------------------
## 2. 作图（CNS / PNAS 风格）
## -----------------------------
p_cor <- ggplot() +

  ## 箱线图（主体，黑色）
  geom_boxplot(
    data = cor_all,
    aes(x = K, y = rho),
    width = 0.55,
    outlier.size = 0.6,
    outlier.alpha = 0.25,
    linewidth = 0.4,
    color = "black",
    fill = NA
  ) +

  ## 均值点（小而克制）
  stat_summary(
    data = cor_all,
    aes(x = K, y = rho),
    fun = mean,
    geom = "point",
    size = 1.8,
    color = "black"
  ) +

  ## 均值趋势线（强调色，非回归）
  geom_line(
    data = mean_all,
    aes(x = K_index, y = mean_rho, group = 1),
    linewidth = 0.9,
    color = "#B2182B"   # CNS 常用深红
  ) +

  ## 均值锚点
  geom_point(
    data = mean_all,
    aes(x = K_index, y = mean_rho),
    size = 2.2,
    color = "#B2182B"
  ) +

  ## X 轴：离散标签
  scale_x_discrete(
    limits = c("10", "30", "50", "70", "100"),
    name = "Number of scTenifoldKnk runs (K)"
  ) +

  ## Y 轴
  scale_y_continuous(
    name = "Spearman correlation of FC (pairwise)",
    expand = expansion(mult = c(0.05, 0.1))
  ) +

  ## 主题（PNAS / CNS 通用）
  theme_classic(base_size = 11) +
  theme(
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(
      hjust = 0.5,
      face = "plain",
      size = 11
    )
  ) +

  labs(
    title = glue("Reproducibility of inferred FC across repeated scTenifoldKnk runs (gKO = {gKO})")
  )

ggsave("03_plots/01_FC_spearman_boxplot_with_mean_fit.pdf", p_cor, width = 8, height = 5)

# ============================================================
# 5) (3) 各K下：每个基因 log2FC>1 且 p.adj<0.05 的次数
# ============================================================
# 以及 p.adj<0.05 的次数（给第4步上色用）
count_hits_for_K <- function(runs_list, K, log2fc_cut, padj_cut) {
  runsK <- runs_list[1:K]

  # gene全集（通常一致，但稳妥起见取union）
  gene_all <- sort(unique(unlist(lapply(runsK, function(df) df$gene))))

  hit_counts <- setNames(rep(0L, length(gene_all)), gene_all)
  sigp_counts <- setNames(rep(0L, length(gene_all)), gene_all)

  for (df in runsK) {
    df2 <- df %>% select(gene, log2FC, p.adj)

    hit_gene <- df2 %>% filter(log2FC < log2fc_cut, p.adj < padj_cut) %>% pull(gene)
    sig_gene <- df2 %>% filter(p.adj < padj_cut) %>% pull(gene)

    hit_counts[hit_gene] <- hit_counts[hit_gene] + 1L
    sigp_counts[sig_gene] <- sigp_counts[sig_gene] + 1L
  }

  tibble(
    gene = gene_all,
    K = K,
    hit_n = as.integer(hit_counts[gene_all]),
    sigp_n = as.integer(sigp_counts[gene_all])
  )
}

hit_tables <- lapply(K_SET2, function(K) {
  count_hits_for_K(all_runs_ok, K, LOG2FC_CUT, PADJ_CUT)
})
hit_all <- bind_rows(hit_tables)

write.csv(hit_all, "02_summaries/gene_hit_counts_byK.csv", row.names = FALSE)

# 同时输出宽表（便于查看）
hit_wide <- hit_all %>%
  select(gene, K, hit_n) %>%
  pivot_wider(names_from = K, values_from = hit_n, names_prefix = "K_")

sigp_wide <- hit_all %>%
  select(gene, K, sigp_n) %>%
  pivot_wider(names_from = K, values_from = sigp_n, names_prefix = "K_")

write.csv(hit_wide, "02_summaries/gene_hit_counts_byK_wide.csv", row.names = FALSE)
write.csv(sigp_wide, "02_summaries/gene_sigp_counts_byK_wide.csv", row.names = FALSE)

# ============================================================
# 6) (4) 各K下：Top20基因云雨图（雨=散点，云=小提琴）
#     - 横轴：FC
#     - 纵轴：gene
#     - 基因排序：按FC均值排序
#     - 颜色：由“p.adj<0.05的次数(sigp_n)”控制（渐变）
# ============================================================
get_fc_long_for_K <- function(runs_list, K) {
  runsK <- runs_list[1:K]
  fc_long <- bind_rows(lapply(seq_along(runsK), function(i) {
    df <- runsK[[i]]
    df %>%
      transmute(
        run = i,
        gene = gene,
        FC = FC,
        log2FC = log2FC,
        p.adj = p.adj
      )
  }))
  fc_long
}

make_raincloud_plot <- function(K, fc_long, hit_tbl_K, out_pdf) {
  # Top20：按 hit_n 降序
  top_genes <- hit_tbl_K %>%
    arrange(desc(hit_n), desc(sigp_n)) %>%
    slice_head(n = 20) %>%
    pull(gene)

  dat <- fc_long %>% filter(gene %in% top_genes)

  # gene排序：按FC均值（对这K次）
  gene_order <- dat %>%
    group_by(gene) %>%
    summarise(mean_FC = mean(FC, na.rm = TRUE), .groups = "drop") %>%
    arrange(mean_FC) %>%
    pull(gene)

  # 合并颜色信息：sigp_n（p.adj<0.05次数）
  sigp_map <- hit_tbl_K %>% select(gene, sigp_n, hit_n)

  dat2 <- dat %>%
    left_join(sigp_map, by = "gene") %>%
    mutate(gene = factor(gene, levels = gene_order))

  # 云雨图：violin + boxplot + jitter
  p <- ggplot(dat2, aes(x = gene, y = log2FC)) +

    ## 上半部小提琴（云）
    geom_half_violin(
      aes(fill = sigp_n),
      side = "r",
      trim = TRUE,
      width = 0.9,
      alpha = 0.85,
      color = NA
    ) +

    ## 箱线（雨的“核心”）
    geom_boxplot(
      width = 0.12,
      outlier.shape = NA,
      alpha = 0.9,
      linewidth = 0.35
    ) +

    ## 散点（雨）
    geom_jitter(
      aes(color = sigp_n),
      width = 0.08,
      height = 0,
      size = 1.2,
      alpha = 0.6
    ) +

    coord_flip() +

    ## CNS 风格连续渐变（小提琴 & 点共用）
    scale_fill_gradient(
      low = "grey85",
      high = "#B2182B"   # Nature / CNS 常用暗红
    ) +
    scale_color_gradient(
      low = "grey60",
      high = "#B2182B"
    ) +

    labs(
      title = glue(
        "Top 20 genes ranked by hit frequency (log2FC > {LOG2FC_CUT}, adj. P < {PADJ_CUT})\n",
        "K = {K}, gKO = {gKO}"
      ),
      x = "Gene (ordered by mean log2FC)",
      y = "log2(FC)",
      fill  = "Count (adj. P < 0.05)",
      color = "Count (adj. P < 0.05)"
    ) +

    theme_classic(base_size = 11) +
    theme(
      axis.line = element_line(linewidth = 0.4),
      axis.ticks = element_line(linewidth = 0.4),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 9),
      legend.position = "right",
      plot.title = element_text(
        hjust = 0,
        size = 11,
        face = "plain"
      )
    )

  ggsave(out_pdf, p, width = 8.5, height = 6)

  p
}

for (K in K_SET2) {
  message(glue("Make raincloud plot for K={K}"))

  fc_long_K <- get_fc_long_for_K(all_runs_ok, K)

  hit_tbl_K <- hit_all %>% filter(K == !!K)

  out_pdf <- file.path("03_plots", glue("02_raincloud_Top20_byHit_K{K}.pdf"))
  make_raincloud_plot(K, fc_long_K, hit_tbl_K, out_pdf)
}

# ============================================================
# 7) 额外输出：每个K下 Top20基因列表（方便检查）
# ============================================================
top20_list <- bind_rows(lapply(K_SET2, function(K) {
  (hit_all %>% filter(K == !!K) %>%
     arrange(desc(hit_n), desc(sigp_n)) %>%
     slice_head(n = 20) %>%
     mutate(rank = row_number()) %>%
     select(K, rank, gene, hit_n, sigp_n))
}))

write.csv(top20_list, "02_summaries/top20_genes_byHit_eachK.csv", row.names = FALSE)

message("All done.")
```

</details>

<details open>
<summary><code>makesure_02dpi_samples_all_have_Early_signal_cell.py</code> · Python · 141 行</summary>

[在 GitHub 上查看此文件](https://github.com/Mikotoo/Lotus_genome/blob/main/scripts/08_single_nucleus/Script/makesure_02dpi_samples_all_have_Early_signal_cell.py)

```python
#!/usr/bin/env python3

import anndata
import pandas as pd
import matplotlib.pyplot as plt

# ===============================
# 1. 读取 h5ad
# ===============================
adata = anndata.read_h5ad("./02_processData/data_clusterd_all_gene_4000.h5ad")
adata.obs_names_make_unique()

# ===============================
# 2. 注释 celltype
# ===============================
cluster2annotation = {
    "0": "5",
    "1": "5",
    "2": "3",
    "3": "4",
    "4": "4",
    "5": "5",
    "6": "6",
    "7": "2",
    "8": "7",
    "9": "1"
}

adata.obs["celltype"] = (
    adata.obs["leiden_0.4"]
    .astype(str)
    .map(cluster2annotation)
    .astype("category")
)

# ===============================
# 3. 读取 barcode
# ===============================
cell02dpi1 = pd.read_csv(
    "02dpi_1_barcodes.tsv",
    header=None,
    sep="\t"
)[0].astype(str).tolist()

cell02dpi2 = pd.read_csv(
    "02dpi_2_barcodes.tsv",
    header=None,
    sep="\t"
)[0].astype(str).tolist()

# ===============================
# 4. 处理 obs_names
# adata.obs_names: AAAC...-1-0
# barcode 文件: AAAC...-1
# 去掉最后的 -0 / -1 / -2 ...
# ===============================
obs = adata.obs.copy()
obs["cell_id"] = obs.index
obs["barcode"] = obs["cell_id"].str.replace(r"-[0-9]+$", "", regex=True)

obs["group"] = "other"
obs.loc[obs["barcode"].isin(cell02dpi1), "group"] = "cell02dpi1"
obs.loc[obs["barcode"].isin(cell02dpi2), "group"] = "cell02dpi2"

# ===============================
# 5. 加入 UMAP 坐标
# ===============================
obs["UMAP_1"] = adata.obsm["X_umap"][:, 0]
obs["UMAP_2"] = adata.obsm["X_umap"][:, 1]

# ===============================
# 6. 统计 celltype == 7 数量
# ===============================
n_02dpi1_total = (obs["group"] == "cell02dpi1").sum()
n_02dpi2_total = (obs["group"] == "cell02dpi2").sum()

n_02dpi1_ct7 = ((obs["group"] == "cell02dpi1") & (obs["celltype"].astype(str) == "7")).sum()
n_02dpi2_ct7 = ((obs["group"] == "cell02dpi2") & (obs["celltype"].astype(str) == "7")).sum()

print(f"cell02dpi1 total cells: {n_02dpi1_total}")
print(f"cell02dpi1 celltype==7 cells: {n_02dpi1_ct7}")

print(f"cell02dpi2 total cells: {n_02dpi2_total}")
print(f"cell02dpi2 celltype==7 cells: {n_02dpi2_ct7}")

# ===============================
# 7. 绘图函数
# ===============================
def plot_highlight_umap(
    obs,
    group_name,
    outfile_prefix
):
    target = (obs["group"] == group_name) & (obs["celltype"].astype(str) == "7")
    other = ~target

    plt.figure(figsize=(5, 4.5))

    plt.scatter(
        obs.loc[other, "UMAP_1"],
        obs.loc[other, "UMAP_2"],
        s=2,
        c="lightgrey",
        alpha=0.5,
        linewidths=0
    )

    plt.scatter(
        obs.loc[target, "UMAP_1"],
        obs.loc[target, "UMAP_2"],
        s=6,
        c="#D7B76A",
        alpha=0.9,
        linewidths=0
    )

    plt.title(f"{group_name}: celltype 7", fontsize=14, fontweight="bold")
    plt.xlabel("UMAP 1")
    plt.ylabel("UMAP 2")
    plt.tight_layout()

    plt.savefig(f"{outfile_prefix}.pdf")
    plt.savefig(f"{outfile_prefix}.png", dpi=300)
    plt.close()

# ===============================
# 8. 分别绘制两个 UMAP
# ===============================
plot_highlight_umap(
    obs,
    group_name="cell02dpi1",
    outfile_prefix="cell02dpi1_celltype7_umap"
)

plot_highlight_umap(
    obs,
    group_name="cell02dpi2",
    outfile_prefix="cell02dpi2_celltype7_umap"
)

print("Done.")
```

</details>
