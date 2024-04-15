import sys
import anndata
import pandas as pd
import numpy as np
import scanpy as sc
import scipy.sparse as sp
import scvi
import scripts
import scripts.scDblFinder
from scripts.scDblFinder import run_ScDblFinder
import matplotlib.pyplot as plt

data_integrated = anndata.read_h5ad("/share/home/yzwl_zhangchao/Project/soybean_sn/02_QC/_processData/data_integrated.h5ad")
data = data_integrated


sc.pp.pca(data,layer='scvi_normalized')
data.obsm["scvi_pca"]=data.obsm["X_pca"]

##### tsne
sc.tl.tsne(data,use_rep="X_scVI")
data.obsm["scvi_tsne"] = data.obsm["X_tsne"]

sc.tl.tsne(data,use_rep="scvi_pca")
data.obsm["scvi_pca_tsne"] = data.obsm["X_tsne"]


#### UMAP
sc.pp.neighbors(data,n_neighbors=15,use_rep='X_scVI')
sc.tl.umap(data)
data.obsm["scvi_umap"]=data.obsm["X_umap"]

sc.pp.neighbors(data,n_neighbors=15,use_rep='scvi_pca')
sc.tl.umap(data)
data.obsm["scvi_pca_umap"]=data.obsm["X_umap"]


fig,axes=plt.subplots(3,2,figsize=(30,20))

sc.pl.embedding(data,basis='scvi_pca',frameon='small',ax=axes[0,1],color="Sample")
sc.pl.embedding(data,basis='scvi_tsne',frameon='small',ax=axes[1,0],color="Sample")
sc.pl.embedding(data,basis='scvi_pca_tsne',frameon='small',ax=axes[1,1],color="Sample")
sc.pl.embedding(data,basis='scvi_umap',frameon='small',ax=axes[2,0],color="Sample")
sc.pl.embedding(data,basis='scvi_pca_umap',frameon='small',ax=axes[2,1],color="Sample")
plt.savefig("figures/dimenRedu.png")

data.write_h5ad("_processData/data_dimenRedu.h5ad")

