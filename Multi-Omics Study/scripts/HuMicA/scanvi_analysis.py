#!/usr/bin/env python3
"""
scANVI Analysis Pipeline for HuMicA Dataset

This script performs scANVI analysis on the HuMicA single-cell RNA-seq dataset,
following best practices for reproducibility, GPU utilization, and model training.

Author: Generated for Complement-OUD project
Date: 2025
"""

import os
import sys
import warnings
import logging
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('scanvi_analysis.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# Configuration
class Config:
    """Configuration parameters for the analysis"""
    # File paths
    INPUT_H5AD = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study/data/raw/HuMicA/HuMicA_for_scVI.h5ad"
    OUTPUT_DIR = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study/results/scanvi"

    # Model parameters
    BATCH_KEY = "Sample_ID"
    LABELS_KEY = "seurat_clusters"
    N_LATENT = 30
    N_LAYERS = 2
    N_HIDDEN = 128

    # Training parameters
    MAX_EPOCHS = 400
    EARLY_STOPPING = True
    EARLY_STOPPING_PATIENCE = 45
    LEARNING_RATE = 1e-3
    BATCH_SIZE = 128

    # Quality control
    MIN_GENES = 200
    MIN_CELLS = 3
    MAX_GENES_BY_COUNTS = 5000
    MITO_THRESHOLD = 20

    # Reproducibility
    RANDOM_SEED = 42

    # Visualization
    FIGSIZE = (10, 8)
    DPI = 300

def setup_environment():
    """Setup analysis environment and check requirements"""
    logger.info("Setting up analysis environment...")

    # Set random seeds for reproducibility
    np.random.seed(Config.RANDOM_SEED)
    torch.manual_seed(Config.RANDOM_SEED)
    scvi.settings.seed = Config.RANDOM_SEED

    # Configure scanpy
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=Config.DPI, facecolor='white')

    # Check GPU availability
    if torch.cuda.is_available():
        device = torch.cuda.get_device_name(0)
        logger.info(f"GPU available: {device}")
        scvi.settings.dl_pin_memory_gpu_training = True
    else:
        logger.info("GPU not available, using CPU")

    # Create output directory
    output_path = Path(Config.OUTPUT_DIR)
    output_path.mkdir(parents=True, exist_ok=True)
    logger.info(f"Output directory: {output_path}")

    # Log package versions
    logger.info(f"scanpy=={sc.__version__}")
    logger.info(f"scvi-tools=={scvi.__version__}")
    logger.info(f"torch=={torch.__version__}")

def load_and_inspect_data():
    """Load H5AD file and perform initial inspection"""
    logger.info("Loading data...")

    if not os.path.exists(Config.INPUT_H5AD):
        raise FileNotFoundError(f"Input file not found: {Config.INPUT_H5AD}")

    adata = sc.read_h5ad(Config.INPUT_H5AD)
    logger.info(f"Loaded data: {adata.n_obs} cells × {adata.n_vars} genes")

    # Log metadata info
    logger.info(f"Available metadata columns: {list(adata.obs.columns)}")
    logger.info(f"Batch key '{Config.BATCH_KEY}' unique values: {adata.obs[Config.BATCH_KEY].nunique()}")
    logger.info(f"Labels key '{Config.LABELS_KEY}' unique values: {adata.obs[Config.LABELS_KEY].nunique()}")

    # Check if raw counts are in X
    logger.info(f"Data matrix type: {type(adata.X)}")
    logger.info(f"Data range: {adata.X.min():.2f} to {adata.X.max():.2f}")

    return adata

def quality_control(adata):
    """Perform quality control filtering"""
    logger.info("Performing quality control...")

    # Store original dimensions
    n_cells_orig, n_genes_orig = adata.n_obs, adata.n_vars

    # Calculate QC metrics
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

    # Add mitochondrial gene percentage
    adata.obs['pct_counts_mt'] = (adata.obs['total_counts_mt'] / adata.obs['total_counts']) * 100

    # Plot QC metrics
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                 jitter=0.4, multi_panel=True, ax=axes)
    plt.tight_layout()
    plt.savefig(f"{Config.OUTPUT_DIR}/qc_metrics_before_filtering.png", dpi=Config.DPI, bbox_inches='tight')
    plt.close()

    # Filter cells and genes
    sc.pp.filter_cells(adata, min_genes=Config.MIN_GENES)
    sc.pp.filter_genes(adata, min_cells=Config.MIN_CELLS)

    # Filter cells with too many genes or high mitochondrial content
    adata = adata[adata.obs.n_genes_by_counts < Config.MAX_GENES_BY_COUNTS, :]
    adata = adata[adata.obs.pct_counts_mt < Config.MITO_THRESHOLD, :]

    logger.info(f"After QC filtering: {adata.n_obs} cells ({n_cells_orig - adata.n_obs} removed), "
                f"{adata.n_vars} genes ({n_genes_orig - adata.n_vars} removed)")

    return adata

def setup_scanvi_model(adata):
    """Setup scANVI model with proper configuration"""
    logger.info("Setting up scANVI model...")

    # Ensure categorical data types
    adata.obs[Config.BATCH_KEY] = adata.obs[Config.BATCH_KEY].astype('category')
    adata.obs[Config.LABELS_KEY] = adata.obs[Config.LABELS_KEY].astype('category')

    # Setup AnnData for scANVI
    scvi.model.SCANVI.setup_anndata(
        adata,
        layer=None,  # Use adata.X (raw counts)
        batch_key=Config.BATCH_KEY,
        labels_key=Config.LABELS_KEY,
        categorical_covariate_keys=None,
        continuous_covariate_keys=None
    )

    # Create model
    model = scvi.model.SCANVI(
        adata,
        n_latent=Config.N_LATENT,
        n_layers=Config.N_LAYERS,
        n_hidden=Config.N_HIDDEN,
        dropout_rate=0.1,
        dispersion='gene',
        gene_likelihood='nb'
    )

    logger.info(f"Model created with {Config.N_LATENT} latent dimensions")
    logger.info(f"Model architecture: {Config.N_LAYERS} layers, {Config.N_HIDDEN} hidden units")

    return model

def train_model(model, adata):
    """Train the scANVI model with monitoring"""
    logger.info("Training scANVI model...")

    # Train the model
    model.train(
        max_epochs=Config.MAX_EPOCHS,
        batch_size=Config.BATCH_SIZE,
        early_stopping=Config.EARLY_STOPPING,
        early_stopping_patience=Config.EARLY_STOPPING_PATIENCE,
        lr=Config.LEARNING_RATE,
        use_gpu=torch.cuda.is_available()
    )

    # Plot training history
    train_elbo = model.history["elbo_train"]
    val_elbo = model.history["elbo_validation"]

    plt.figure(figsize=Config.FIGSIZE)
    plt.plot(train_elbo, label='Training ELBO')
    plt.plot(val_elbo, label='Validation ELBO')
    plt.xlabel('Epoch')
    plt.ylabel('ELBO')
    plt.title('scANVI Training History')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(f"{Config.OUTPUT_DIR}/training_history.png", dpi=Config.DPI, bbox_inches='tight')
    plt.close()

    logger.info(f"Training completed after {len(train_elbo)} epochs")

    return model

def extract_results(model, adata):
    """Extract latent representations and predictions"""
    logger.info("Extracting results...")

    # Get latent representation
    adata.obsm["X_scANVI"] = model.get_latent_representation()

    # Get normalized expression
    adata.obsm["X_scANVI_normalized"] = model.get_normalized_expression(library_size=1e4)

    # Get predictions for unlabeled cells (if any)
    predictions = model.predict()
    adata.obs["scanvi_predictions"] = predictions
    adata.obs["scanvi_predictions"] = adata.obs["scanvi_predictions"].astype('category')

    # Get prediction probabilities
    pred_probs = model.predict(soft=True)
    adata.obsm["scanvi_pred_probs"] = pred_probs

    logger.info("Results extracted successfully")

    return adata

def run_downstream_analysis(adata):
    """Perform downstream analysis on latent space"""
    logger.info("Running downstream analysis...")

    # Compute neighborhood graph on latent space
    sc.pp.neighbors(adata, use_rep="X_scANVI", n_neighbors=15, n_pcs=None)

    # Compute UMAP
    sc.tl.umap(adata, min_dist=0.3, spread=1.0)

    # Leiden clustering on latent space
    sc.tl.leiden(adata, resolution=0.5, key_added="leiden_scanvi")

    logger.info("Downstream analysis completed")

def create_visualizations(adata):
    """Create comprehensive visualizations"""
    logger.info("Creating visualizations...")

    # UMAP plots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # Original clusters
    sc.pl.umap(adata, color=Config.LABELS_KEY, legend_loc='on data',
               title='Original Seurat Clusters', ax=axes[0,0], show=False)

    # scANVI predictions
    sc.pl.umap(adata, color="scanvi_predictions", legend_loc='on data',
               title='scANVI Predictions', ax=axes[0,1], show=False)

    # Batch effect
    sc.pl.umap(adata, color=Config.BATCH_KEY, title='Batch (Sample ID)',
               ax=axes[1,0], show=False)

    # New clustering
    sc.pl.umap(adata, color="leiden_scanvi", legend_loc='on data',
               title='Leiden Clustering (scANVI latent)', ax=axes[1,1], show=False)

    plt.tight_layout()
    plt.savefig(f"{Config.OUTPUT_DIR}/umap_overview.png", dpi=Config.DPI, bbox_inches='tight')
    plt.close()

    # Batch integration assessment
    plt.figure(figsize=Config.FIGSIZE)
    sc.pl.umap(adata, color=Config.BATCH_KEY, title='Batch Integration Assessment')
    plt.savefig(f"{Config.OUTPUT_DIR}/batch_integration.png", dpi=Config.DPI, bbox_inches='tight')
    plt.close()

    logger.info("Visualizations saved")

def save_results(model, adata):
    """Save model and processed data"""
    logger.info("Saving results...")

    # Save model
    model_path = f"{Config.OUTPUT_DIR}/scanvi_model"
    model.save(model_path, overwrite=True)
    logger.info(f"Model saved to: {model_path}")

    # Save processed AnnData
    adata_path = f"{Config.OUTPUT_DIR}/adata_scanvi_processed.h5ad"
    adata.write_h5ad(adata_path)
    logger.info(f"Processed data saved to: {adata_path}")

    # Save latent embeddings as CSV
    latent_df = pd.DataFrame(
        adata.obsm["X_scANVI"],
        index=adata.obs.index,
        columns=[f"scANVI_{i+1}" for i in range(adata.obsm["X_scANVI"].shape[1])]
    )
    latent_df.to_csv(f"{Config.OUTPUT_DIR}/scanvi_latent_embeddings.csv")

    # Save metadata with predictions
    metadata_cols = [Config.BATCH_KEY, Config.LABELS_KEY, "scanvi_predictions", "leiden_scanvi"]
    metadata_df = adata.obs[metadata_cols].copy()
    metadata_df.to_csv(f"{Config.OUTPUT_DIR}/scanvi_metadata.csv")

    logger.info("All results saved successfully")

def generate_summary_report(adata, model):
    """Generate a summary report of the analysis"""
    logger.info("Generating summary report...")

    report = f"""
scANVI Analysis Summary Report
Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}

Dataset Information:
- Input file: {Config.INPUT_H5AD}
- Final dimensions: {adata.n_obs} cells × {adata.n_vars} genes
- Batch key: {Config.BATCH_KEY} ({adata.obs[Config.BATCH_KEY].nunique()} batches)
- Labels key: {Config.LABELS_KEY} ({adata.obs[Config.LABELS_KEY].nunique()} clusters)

Model Configuration:
- Latent dimensions: {Config.N_LATENT}
- Architecture: {Config.N_LAYERS} layers, {Config.N_HIDDEN} hidden units
- Training epochs: {len(model.history["elbo_train"])}
- Batch size: {Config.BATCH_SIZE}

Results:
- scANVI latent representation: X_scANVI ({adata.obsm["X_scANVI"].shape})
- New Leiden clusters: {adata.obs["leiden_scanvi"].nunique()}
- GPU used: {torch.cuda.is_available()}

Output Files:
- Processed data: adata_scanvi_processed.h5ad
- Model: scanvi_model/
- Latent embeddings: scanvi_latent_embeddings.csv
- Metadata: scanvi_metadata.csv
- Visualizations: umap_overview.png, batch_integration.png
- Training history: training_history.png
"""

    with open(f"{Config.OUTPUT_DIR}/analysis_summary.txt", "w") as f:
        f.write(report)

    logger.info("Summary report generated")
    print(report)

def main():
    """Main analysis pipeline"""
    try:
        # Setup
        setup_environment()

        # Load and inspect data
        adata = load_and_inspect_data()

        # Quality control (optional - data might already be filtered)
        # Uncomment if you want additional QC filtering
        # adata = quality_control(adata)

        # Setup and train model
        model = setup_scanvi_model(adata)
        model = train_model(model, adata)

        # Extract results and run downstream analysis
        adata = extract_results(model, adata)
        run_downstream_analysis(adata)

        # Create visualizations
        create_visualizations(adata)

        # Save results
        save_results(model, adata)

        # Generate summary
        generate_summary_report(adata, model)

        logger.info("scANVI analysis completed successfully!")

    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        raise

if __name__ == "__main__":
    main()
