'''
🔬 1. Preprocessing & Quality Control (QC)
A. Load data
    Tool: scanpy (.h5ad, .mtx, 10x output, etc.)
    Considerations: If available, use spliced/unspliced data for future RNA velocity.
        import scanpy as sc
        adata = sc.read_10x_h5("filtered_feature_bc_matrix.h5")
B. Initial QC
    Metrics:
        n_genes_by_counts
        total_counts
        percent_mito (and optionally ribosomal, nuclear, etc.)
        Use adata.var['feature_types'] to restrict to nuclear RNA
            import numpy as np
            adata.var['mt'] = adata.var_names.str.startswith('MT-')
            sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
    Thresholding (with robust stats):
        Median absolute deviation (MAD) or quantile-based filtering
        Consider sex- and tissue-aware filtering using miQC (R package)

🔁 2. Normalization & Feature Selection

    A. Ambient RNA correction (important for snRNA-seq)
        Tool: CellBender (deep learning)
        Alternative: SoupX (if using droplet-based data)
    B. Normalization
        Pearson residuals normalization from SCTransform (Seurat) – best variance stabilization
            In Python: pySCTransform or fallback to:
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
        Optional: sc.experimental.pp.highly_variable_genes() (newer method)
    C. Batch Correction (if applicable)
        Best in class:
            scVI/scANVI (deep generative model)
            Harmony (faster, interpretable)
            scGen or trVAE (if conditioning on experimental covariates)
                # scVI example
                import scvi
                scvi.model.SCVI.setup_anndata(adata, batch_key='sample')
                model = scvi.model.SCVI(adata)
                model.train()
                adata.obsm["X_scVI"] = model.get_latents()
🔍 3. Dimensionality Reduction & Clustering

    A. PCA → UMAP/tSNE
            Consider diffusion maps or PHATE for developmental trajectories
                sc.pp.pca(adata)
                sc.pp.neighbors(adata)
                sc.tl.umap(adata)

    B. Clustering
        Leiden or Louvain, tune resolution
        For stability: clustree (R), consensus clustering, or resampling-based metrics
'''

# Import necessary libraries
import scanpy as sc
import pandas as pd
import numpy as np

