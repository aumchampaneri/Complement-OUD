# Prepare the environment
# !pip install gget --quiet
# !pip install harmonypy --quiet
# !pip install scanpy --quiet

import gget
import scanpy as sc
from scanpy.experimental.pp import recipe_pearson_residuals

# Setup CELLxGENE
# gget.setup("cellxgene")

# Download mouse brain dataset from CELLxGENE (~3.5 million cells)
adata = gget.cellxgene(
    census_version="2025-01-30",
    species="mus_musculus",
    ensembl=True,
    dataset_id=[
        "1229ecc2-b067-4664-91da-0251aec31574",
        "98e5ea9f-16d6-47ec-a529-686e76515e39",
        "58b01044-c5e5-4b0f-8a2d-6ebf951e01ff",
    ],
    disease="normal",
    tissue_general="brain"
)

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

# Save the processed file
adata_processed_path = "/Super Folder - Mus musculus/gget census UMAP/Data Files/Mouse_Census_Brain_PP.h5ad"
adata.write_h5ad(adata_processed_path)