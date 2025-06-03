'''
🔬 Single-nucleus RNA-seq Processing Pipeline - scVI Version
GSE225158 - OUD Striatum snRNA-seq data reprocessing with scVI

Alternative to Pearson residuals using probabilistic modeling
- Better batch integration across 22 subjects
- Improved statistical framework for clinical comparisons
- State-of-the-art normalization and dimensionality reduction
'''

import scanpy as sc
import scvi
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings('ignore')

# Use same paths as main pipeline
from scripts.snrna.preprocessing_qc import (
    INPUT_PATH, PROCESSED_DATA_DIR, RESULTS_DIR, PLOTS_DIR,
    calculate_mad_outliers
)

# Output files for scVI version
OUTPUT_H5AD_SCVI = "GSE225158_reprocessed_scvi.h5ad"
OUTPUT_H5AD_PATH_SCVI = f"{PROCESSED_DATA_DIR}/{OUTPUT_H5AD_SCVI}"

def main_scvi():
    """scVI-based processing pipeline"""
    
    print("=" * 70)
    print("🔬 scVI-BASED PREPROCESSING PIPELINE")
    print("=" * 70)
    
    # ...existing data loading and QC code...
    # (Same as main pipeline through Step 2)
    
    # =======================================
    # 🔁 STEP 3: scVI SETUP & TRAINING
    # =======================================
    print("\n" + "=" * 70)
    print("🔁 STEP 3: scVI MODEL SETUP & TRAINING")
    print("=" * 70)
    
    # Set up scVI
    scvi.settings.seed = 42
    raw_adata.raw = raw_adata  # Save raw counts
    
    # Setup anndata for scVI
    batch_key = 'ID'  # Use subject ID as batch
    if batch_key in raw_adata.obs.columns:
        scvi.model.SCVI.setup_anndata(
            raw_adata,
            batch_key=batch_key,
            continuous_covariate_keys=['pct_counts_mt', 'total_counts']
        )
        print(f"✅ scVI setup complete with batch_key='{batch_key}'")
    else:
        scvi.model.SCVI.setup_anndata(raw_adata)
        print("✅ scVI setup complete (no batch correction)")
    
    # Train scVI model
    print("🏋️  Training scVI model...")
    model = scvi.model.SCVI(raw_adata, n_latent=50, n_hidden=256, n_layers=2)
    model.train(max_epochs=300, early_stopping=True, batch_size=1024)
    
    # Get latent representation
    print("🔍 Extracting latent representation...")
    raw_adata.obsm["X_scvi"] = model.get_latent_representation()
    
    # Compute neighborhood graph from scVI latent space
    sc.pp.neighbors(raw_adata, use_rep="X_scvi", n_neighbors=15)
    sc.tl.umap(raw_adata, random_state=42)
    
    # Clustering
    sc.tl.leiden(raw_adata, resolution=0.5, random_state=42)
    
    print("✅ scVI pipeline complete!")
    
    # Save results
    raw_adata.write(OUTPUT_H5AD_PATH_SCVI)
    print(f"💾 scVI results saved to: {OUTPUT_H5AD_PATH_SCVI}")
    
    return raw_adata

if __name__ == "__main__":
    scvi_adata = main_scvi()
