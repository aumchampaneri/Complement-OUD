
# Prepare the environment
import gget
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from numpy.lib.npyio import savez

gget.setup("cellxgene")

# Uncomment the following line to see the documentation
# help(gget.cellxgene)

'''
Query the cellxgene database for the genes of interest
Filters:
disease: normal
tissue_general: brain

> Query Time: 3m 30s
'''

adata = gget.cellxgene(
    ensembl=True,
    verbose=True,
    species='mus_musculus',
    disease='normal',
    tissue_general='brain'
)

all_genes = {
"""
Dictionary mapping a subset of gene names to their corresponding Ensembl IDs.

Keys:
   str: Gene names (e.g., 'C3', 'C3AR1', etc.)

Values:
   str: Ensembl IDs (e.g., 'ENSG00000125730', 'ENSG00000171860', etc.)
"""
    "C3": "ENSMUSG00000024164",
    "C3AR1": "ENSMUSG00000040552",
    "C5": "ENSMUSG00000026874",
    "C5ar1": "ENSMUSG00000049130",
    "C5ar2": "ENSMUSG00000074361",
    "Cfh": "ENSMUSG00000026365",
    "C1qa": "ENSMUSG00000036887",
    "C1qb": "ENSMUSG00000036905",
    "C1qc": "ENSMUSG00000036896",
    "C4a": "ENSMUSG00000015451",
    "C4b": "ENSMUSG00000073418"
}


"""
Preprocess the kidney dataset by filtering genes, normalizing counts,
log-transforming the data, and scaling the data. Then, perform PCA
and plot an overview of the PCA results.

Steps:
1. Filter genes with at least 1 count.
2. Normalize the total counts to 1e6.
3. Log-transform the data.
4. Scale the data.
5. Perform PCA.
6. Plot an overview of the PCA results.
"""

adata.write_h5ad(filename="/Users/aumchampaneri/PycharmProjects/Complement-OUD/Mm Census gget Query/Data Files/Mouse_Census_Brain.h5ad")

sc.pp.filter_genes(adata, min_counts=3)
sc.pp.filter_cells(adata, min_genes=500)
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)
# sc.pp.scale(adata) # Uses too much RAM -> swap to a zero center pca

# Save the raw (normalized + log-transformed, but unscaled) data
adata.raw = adata.copy()

# Calculate PCA
sc.tl.pca(adata, zero_center=False)

# Algorithically integrate multiple experiments
sc.external.pp.harmony_integrate(adata, key='dataset_id')

# Calculate UMAP
sc.pp.neighbors(adata, use_rep='X_pca_harmony')
sc.tl.umap(adata)

adata.write_h5ad(
    filename="/Users/aumchampaneri/PycharmProjects/Complement-OUD/Mm Census gget Query/Data Files/Mouse_Census_Brain_PP.h5ad",
)