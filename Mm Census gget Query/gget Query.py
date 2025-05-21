# Prepare the environment
import gget
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc

gget.setup("cellxgene")

# Query the cellxgene database for glial and immune cells from normal mouse brain
# This pulls only cells relevant to memory-constrained glia/immune-focused analyses
adata = gget.cellxgene(
    ensembl=True,
    verbose=True,
    species='mus_musculus',
    disease='normal',
    tissue_general='brain',
    cell_type=[
        # Glial cells
        "astrocyte",
        "oligodendrocyte",
        "oligodendrocyte precursor cell",
        "microglial cell",
        "Bergmann glial cell",
        "astrocyte of the cerebellum",
        "glioblast",
        "macroglial cell",
        "ependymal cell",
        "tanycyte",
        "immature astrocyte",
        "hypendymal cell",

        # Immune cells
        "macrophage",
        "CD4-positive, alpha-beta T cell",
        "CD8-positive, alpha-beta T cell",
        "B cell",
        "natural killer cell",
        "monocyte",
        "innate lymphoid cell",
        "meningeal macrophage",
        "perivascular macrophage",
        "central nervous system macrophage",
        "plasma cell",
        "mature NK T cell",
        "mast cell",
        "conventional dendritic cell",
        "plasmacytoid dendritic cell",
        "T cell"
    ]
)

# Optional: Add a minimal subset of neurons for reference or QC
# Downsample neuron cell types post-query if needed
if "cell_type" in adata.obs.columns:
    neuron_like = adata[adata.obs["cell_type"].str.contains("neuron", case=False, na=False)]
    neuron_sample = neuron_like[np.random.choice(neuron_like.shape[0], size=min(5000, neuron_like.shape[0]), replace=False)]
    adata = adata.concatenate(neuron_sample, batch_key="subset", batch_categories=["glia_immune", "neurons_subsampled"])

# Save unprocessed raw data
adata.write_h5ad("/Users/aumchampaneri/PycharmProjects/Complement-OUD/Mm Census gget Query/Data Files/Mouse_Census_Brain.h5ad")

# Basic filtering and preprocessing
sc.pp.filter_genes(adata, min_counts=3)
sc.pp.filter_cells(adata, min_genes=500)
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)

# Save normalized but unscaled version
adata.raw = adata.copy()

# PCA without zero-centering (to save RAM)
sc.tl.pca(adata, zero_center=False)

# Harmony integration across datasets
sc.external.pp.harmony_integrate(adata, key='dataset_id')

# Neighbors + UMAP
sc.pp.neighbors(adata, use_rep='X_pca_harmony')
sc.tl.umap(adata)

# Save final processed object
adata.write_h5ad("/Users/aumchampaneri/PycharmProjects/Complement-OUD/Mm Census gget Query/Data Files/Mouse_Census_Brain_PP.h5ad")
