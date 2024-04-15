"""
3个样本分别导入
step1. 去除双胞
step2. 过滤 每个基因在至少10细胞中，每个细胞基因数[400:4000]，每个细胞DMI[600:6000]
step3. 归一化
"""

import sys
import anndata
import pandas as pd
import numpy as np
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns
import scvi
import scripts
import scripts.scDblFinder
from scripts.scDblFinder import run_ScDblFinder
import matplotlib.pyplot as plt


data1 = sc.read_10x_mtx("/share/home/yzwl_zhangchao/Project/soybean_sn/01_cellRanger/resultDir/step1_cellRanger/nodule_large/nodule_large/outs/filtered_feature_bc_matrix", cache=True)
data2 = sc.read_10x_mtx("/share/home/yzwl_zhangchao/Project/soybean_sn/01_cellRanger/resultDir/step1_cellRanger/nodule_small/nodule_small/outs/filtered_feature_bc_matrix", cache=True)
data3 = sc.read_10x_mtx("/share/home/yzwl_zhangchao/Project/soybean_sn/01_cellRanger/resultDir/step1_cellRanger/root/root/outs/filtered_feature_bc_matrix", cache=True)


def run_plot_scatter(data, x, y, sample, ax):
    """
    绘制散点图并设置标题、x轴和y轴的范围。
    
    参数：
    - data: 数据集
    - x: x轴变量
    - y: y轴变量
    - sample: 样本名
    - ax: 子图对象
    """
    sc.pl.scatter(data, x=x, y=y, ax=ax)
    ax.set_title(sample)
    ax.set_xlim(0, 50000)
    ax.set_ylim(0, 16000)


datasets = {}
samples = ["nodule_large", "nodule_small", "root"]
fig, aex = plt.subplots(ncols=3, nrows=3, figsize=(12, 15))
for idx, data in enumerate([data1, data2, data3]):
    filename = f"_processData/data{idx + 1}_filtered.h5ad"
    datasets[str(data) + "_raw"] = data
    data.obs["Sample"] = samples[idx]
    sc.pp.calculate_qc_metrics(data, inplace=True, percent_top=None, log1p=False, )
    
    # 绘制原始数据散点图
    run_plot_scatter(data, 'total_counts', 'n_genes_by_counts', samples[idx], aex[0, idx])
    
    # 去除双胞
    run_ScDblFinder(data, copy=False, doubletRatio=0.1)
    
    # 绘制去除双胞后的散点图
    run_plot_scatter(data, 'total_counts', 'n_genes_by_counts', samples[idx], aex[1, idx])
    
    # 基因/细胞过滤
    sc.pp.filter_genes(data, min_cells=10)
    sc.pp.filter_cells(data, min_genes=400)
    sc.pp.filter_cells(data, max_genes=4000)
    sc.pp.filter_cells(data, min_counts=600)
    sc.pp.filter_cells(data, max_counts=6000)
    
    # 绘制过滤后的散点图
    run_plot_scatter(data, 'total_counts', 'n_genes_by_counts', samples[idx], aex[2, idx])
    
    # 将数据写入文件
    data.write_h5ad(filename)

# 保存图形
plt.savefig("figures/scatter.png")

#合并数据并保存
data_concatenated = data1.concatenate(data2,data3)
data_concatenated.write_h5ad("_processData/data_concatenated.h5ad")

## data_concatenated = anndata.read_h5ad("/share/home/yzwl_zhangchao/Project/soybean_sn/02_QC/_processData/data_concatenated.h5ad")

data = data_concatenated
plt.figure(figsize=(12, 12))
sc.pl.scatter(data_concatenated, x='total_counts', y='n_genes_by_counts')
plt.savefig("figures/concatenated_scatter.png")

#归一化
data.layers['counts'] = data.X.copy()
sc.pp.normalize_total(data,target_sum=1e4,inplace=True)
data_scaled = data
sc.pp.log1p(data_scaled)
data_scaled.write_h5ad("_processData/data_scaled.h5ad")

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
sns.histplot(data_concatenated.obs["total_counts"], bins=100,kde=True,ax=axes[0])
axes[0].set_title("Total counts")
sns.histplot(data_scaled.X.sum(1), bins=100,kde=True,ax=axes[1])
axes[1].set_title("Shifted logarithm")
plt.savefig("figures/hist_counts.png")




