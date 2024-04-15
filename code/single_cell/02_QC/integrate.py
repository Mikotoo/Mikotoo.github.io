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


data_scaled = anndata.read_h5ad("/share/home/yzwl_zhangchao/Project/soybean_sn/02_QC/_processData/data_scaled.h5ad")
data = data_scaled
sc.experimental.pp.highly_variable_genes(
    data, 
    flavor="pearson_residuals",
    layer='counts',
    batch_key = "Sample",
    n_top_genes=5000,
    subset=True,
    inplace=True,
)

scvi.model.SCVI.setup_anndata(data, layer = "counts",
                             categorical_covariate_keys=["Sample"],
                             continuous_covariate_keys=['total_counts'])

model = scvi.model.SCVI(data)

model.train()

latent = model.get_latent_representation()
latent.shape
data.obsm['X_scVI'] = latent
data.layers['scvi_normalized'] = model.get_normalized_expression(library_size = 1e4)

data.write_h5ad("_processData/data_integrated.h5ad")

