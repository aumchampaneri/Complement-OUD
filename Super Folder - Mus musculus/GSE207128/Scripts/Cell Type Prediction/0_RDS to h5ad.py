import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
import os

# Activate automatic conversion between R and pandas
pandas2ri.activate()

# Import R packages
base = importr('base')
seurat = importr('Seurat')

# Set paths
input_path = "/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus/GSE207128/Outputs/Processed Data/processed_seurat.rds"
output_dir = "/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus/GSE207128/Outputs/Processed Data/"
output_path = os.path.join(output_dir, "processed_seurat.h5ad")

print("Loading Seurat object...")
# Load the RDS file
seurat_obj = base.readRDS(input_path)

print("Extracting data from Seurat object...")
# Extract count matrix (try different methods for v5 compatibility)
try:
    # Try v5 method first
    counts = seurat.LayerData(seurat_obj, assay="RNA", layer="counts")
except:
    try:
        # Fall back to v4 method
        counts = seurat.GetAssayData(seurat_obj, assay="RNA", slot="counts")
    except:
        # Last resort - direct access
        counts = seurat_obj.slots['assays'].slots['RNA'].slots['counts']

# Convert to pandas DataFrame
counts_df = pandas2ri.rpy2py(counts)

# Get cell names and gene names
cell_names = list(seurat_obj.colnames)
gene_names = list(seurat_obj.rownames)

# Create DataFrame with proper index and columns
counts_df = pd.DataFrame(counts_df.toarray() if hasattr(counts_df, 'toarray') else counts_df,
                        index=gene_names, columns=cell_names)

print(f"Expression matrix shape: {counts_df.shape}")

# Extract metadata
try:
    metadata = pandas2ri.rpy2py(seurat_obj.slots['meta.data'])
    metadata.index = cell_names
except:
    print("Could not extract metadata, creating minimal metadata")
    metadata = pd.DataFrame(index=cell_names)

print(f"Metadata shape: {metadata.shape}")

# Create AnnData object
print("Creating AnnData object...")
adata = ad.AnnData(X=counts_df.T)  # Transpose so cells are rows, genes are columns
adata.obs = metadata
adata.var_names = gene_names
adata.obs_names = cell_names

# Try to extract embeddings if they exist
try:
    # Extract UMAP
    umap_embeddings = pandas2ri.rpy2py(seurat_obj.slots['reductions'].slots['umap'].slots['cell.embeddings'])
    adata.obsm['X_umap'] = umap_embeddings
    print("Added UMAP embeddings")
except:
    print("No UMAP embeddings found")

try:
    # Extract PCA
    pca_embeddings = pandas2ri.rpy2py(seurat_obj.slots['reductions'].slots['pca'].slots['cell.embeddings'])
    adata.obsm['X_pca'] = pca_embeddings
    print("Added PCA embeddings")
except:
    print("No PCA embeddings found")

# Save as H5AD
print(f"Saving H5AD file to: {output_path}")
adata.write_h5ad(output_path)

print("Conversion completed successfully!")
print(f"Final AnnData object:")
print(f"  - Cells: {adata.n_obs}")
print(f"  - Genes: {adata.n_vars}")
print(f"  - Embeddings: {list(adata.obsm.keys())}")