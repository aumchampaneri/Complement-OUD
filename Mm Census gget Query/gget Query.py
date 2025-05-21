import gget
import scanpy as sc
import harmonypy  # Needed for Harmony integration
from scanpy.experimental.pp import recipe_pearson_residuals

# Download mouse brain dataset from CELLxGENE (~3.5 million cells)
gget.setup("cellxgene")
adata = gget.cellxgene(
    species="mus_musculus",
    ensembl=True,
    dataset_id=['e0ed3c55-aff6-4bb7-b6ff-98a2d90b890c',
            '35081d47-99bf-4507-9541-735428df9a9f',
            'dbb4e1ed-d820-4e83-981f-88ef7eb55a35',
            '79a2344d-eddd-45b1-b376-39eddfab1899',
            '1229ecc2-b067-4664-91da-0251aec31574',
            '72eb2332-b308-4014-8d25-95233a9aff1e',
            '3bbb6cf9-72b9-41be-b568-656de6eb18b5',
            '98e5ea9f-16d6-47ec-a529-686e76515e39',
            '58b01044-c5e5-4b0f-8a2d-6ebf951e01ff',
            ],
    disease="normal",
    tissue_general="brain"
)

# Basic filtering
sc.pp.filter_genes(adata, min_counts=3)
sc.pp.filter_cells(adata, min_genes=500)

# Check for required metadata
if 'dataset_id' not in adata.obs.columns:
    raise ValueError("Missing 'dataset_id' in adata.obs, which is needed for Harmony.")

# Run Pearson residuals (replaces normalize, log1p, scale, and hvg)
recipe_pearson_residuals(adata, batch_key='dataset_id')

# Integrate batches using Harmony
adata.obsm["X_pca"] = adata.obsm["X_pca"]  # Harmony expects 'X_pca'
sc.external.pp.harmony_integrate(adata, key='dataset_id')

# Neighbors, clustering, and UMAP
sc.pp.neighbors(adata, use_rep='X_pca_harmony')
sc.tl.umap(adata)

# Save processed file
adata.write_h5ad("/Users/aumchampaneri/PycharmProjects/Complement-OUD/Mm Census gget Query/Data Files/Mouse_Census_Brain_PP.h5ad", compression="gzip")
