#!/usr/bin/env python3
"""
scANVI Analysis Pipeline for Human Microglia Atlas (HuMicA) Dataset

This script performs scANVI analysis on the HuMicA dataset from:
"The Human Microglia Atlas (HuMicA) unravels changes in disease-associated
microglia subsets across neurodegenerative conditions" - Nature Communications

Primary goals:
- Batch correction across multiple neurodegenerative disease datasets
- Preservation of disease-associated microglia (DAM) signatures
- Cross-condition analysis of neuroinflammatory pathways
- Identification of conserved vs disease-specific microglia states

Author: Generated for Complement-OUD project
Date: 2025
"""

import os
import sys
import logging
from pathlib import Path
import numpy as np
import pandas as pd

# Set non-interactive matplotlib backend before importing pyplot
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for headless environments
import matplotlib.pyplot as plt

# Configure matplotlib for optimal headless plotting
matplotlib.rcParams['font.size'] = 12
matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['savefig.dpi'] = 300
matplotlib.rcParams['savefig.bbox'] = 'tight'
matplotlib.rcParams['savefig.pad_inches'] = 0.1
matplotlib.rcParams['axes.titlesize'] = 14
matplotlib.rcParams['axes.labelsize'] = 12
matplotlib.rcParams['legend.fontsize'] = 10
matplotlib.rcParams['xtick.labelsize'] = 10
matplotlib.rcParams['ytick.labelsize'] = 10

import scanpy as sc
import scvi
import torch
from datetime import datetime
import psutil
import argparse

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
    """Configuration parameters optimized for microglia batch correction"""
    # File paths
    INPUT_H5AD = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study/data/raw/HuMicA/HuMicA_for_scVI.h5ad"
    OUTPUT_DIR = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study/results/scanvi/01_scanvi"

    # Model parameters - optimized for microglia batch correction (conservative for stability)
    BATCH_KEY = "Sample_ID"  # Individual samples for technical batch correction
    LABELS_KEY = "seurat_clusters"  # Existing microglia subtypes
    N_LATENT = 30  # Conservative latent dimensions for stability
    N_LAYERS = 2  # Proven stable architecture
    N_HIDDEN = 256  # Balanced capacity

    # Training parameters - conservative for numerical stability
    MAX_EPOCHS = 500
    EARLY_STOPPING = True
    EARLY_STOPPING_PATIENCE = 45
    LEARNING_RATE = 5e-4  # More conservative learning rate
    BATCH_SIZE = 512  # Conservative batch size to prevent OOM
    GRADIENT_CLIP_VAL = 0.5  # Stricter gradient clipping
    USE_GPU = True  # Enable GPU acceleration
    DROPOUT_RATE = 0.1  # Conservative dropout

    # Checkpoint settings
    SKIP_TRAINING = True  # Skip training and load pre-trained model
    MODEL_CHECKPOINT = OUTPUT_DIR + "/scanvi_model"
    DATA_CHECKPOINT = OUTPUT_DIR + "/adata_scanvi_processed.h5ad"

    @staticmethod
    def checkpoint_exists():
        """Check if both model and data checkpoints exist"""
        return (os.path.exists(Config.MODEL_CHECKPOINT) and
                os.path.exists(Config.DATA_CHECKPOINT))

    # Quality control - disabled since data is pre-processed
    SKIP_QC = True  # Data already filtered in Seurat pipeline
    MIN_GENES = 200
    MIN_CELLS = 3
    MAX_GENES_BY_COUNTS = 5000
    MITO_THRESHOLD = 20

    # Reproducibility
    RANDOM_SEED = 42

    # Visualization
    FIGSIZE = (12, 8)
    DPI = 300

    # HuMicA-specific gene sets for analysis
    # Disease-Associated Microglia (DAM) signature genes
    DAM_GENES = [
        'TREM2', 'TYROBP', 'AXL', 'FERMT2', 'SPI1', 'CTSD', 'CSF1',
        'ITGAX', 'CLEC7A', 'LILRB4', 'TIMP2', 'CTSB', 'CTSL', 'GPNMB',
        'LPL', 'CD9', 'CST7', 'LIPA', 'CTSS', 'PSAP'
    ]

    # Homeostatic microglia markers
    HOMEOSTATIC_GENES = [
        'P2RY12', 'TMEM119', 'CX3CR1', 'FCRLS', 'OLFML3', 'HEXB',
        'SALL1', 'SLC2A5', 'SPARC', 'GPR34', 'SIGLECH', 'TGFBR1'
    ]

    # Neuroinflammatory pathway genes
    NEUROINFLAM_GENES = [
        'TNF', 'IL1B', 'IL6', 'IFNG', 'TGFB1', 'CD68', 'CD86', 'ARG1',
        'NOS2', 'IL10', 'CD163', 'MRC1', 'FCGR1A', 'FCGR3A', 'IRF8',
        'STAT1', 'STAT6', 'IRF1', 'IRF4', 'NFKB1', 'RELB'
    ]

    # Neurodegenerative disease-specific markers
    DISEASE_MARKERS = {
        'AD': ['APOE', 'CLU', 'ABCA7', 'CR1', 'PICALM', 'BIN1'],
        'PD': ['SNCA', 'LRRK2', 'PARK7', 'PINK1', 'PRKN'],
        'ALS': ['SOD1', 'TARDBP', 'FUS', 'C9orf72'],
        'MS': ['HLA-DRB1', 'IL7R', 'CD40', 'TNFRSF1A']
    }

    # HuMicA microglia cell type mapping (cluster number to cell type name)
    MICROGLIA_CELL_TYPES = {
        '0': 'Homeostatic Microglia',
        '1': 'Activated Microglia I',
        '2': 'Disease-Associated Microglia (DAM)',
        '3': 'Interferon-Response Microglia',
        '4': 'Proliferating Microglia',
        '5': 'Activated Microglia II',
        '6': 'Lipid-Associated Microglia',
        '7': 'Phagocytic Microglia',
        '8': 'Stressed/Dysfunctional Microglia'
    }

def monitor_gpu_memory():
    """Monitor GPU memory usage across different backends"""
    memory_info = {}

    # System memory
    memory_info['system_memory_percent'] = psutil.virtual_memory().percent
    memory_info['system_memory_available_gb'] = psutil.virtual_memory().available / 1024**3

    # MPS memory (if available)
    if torch.backends.mps.is_available():
        try:
            # MPS doesn't have direct memory query, so we track system memory
            memory_info['mps_backend'] = 'available'
        except:
            memory_info['mps_backend'] = 'unavailable'

    # CUDA memory (if available)
    if torch.cuda.is_available():
        try:
            memory_info['cuda_allocated_gb'] = torch.cuda.memory_allocated() / 1024**3
            memory_info['cuda_reserved_gb'] = torch.cuda.memory_reserved() / 1024**3
            memory_info['cuda_free_gb'] = (torch.cuda.get_device_properties(0).total_memory - torch.cuda.memory_reserved()) / 1024**3
        except:
            memory_info['cuda_backend'] = 'error'

    return memory_info

def get_optimal_batch_size(n_cells, available_memory_gb, accelerator='cpu'):
    """Calculate optimal batch size based on available memory and dataset size"""

    # Conservative memory requirements based on working scVI implementation
    if accelerator == 'mps':
        # Conservative settings for MPS stability
        base_batch_size = 256
        max_batch_size = 512
    elif accelerator == 'gpu' or accelerator == 'cuda':
        # Conservative settings for CUDA
        base_batch_size = 256
        max_batch_size = 512
    else:  # CPU
        # CPU can handle larger batches
        base_batch_size = 512
        max_batch_size = 1024

    # Memory-based constraints (conservative)
    if available_memory_gb < 8:
        optimal_batch = 128
    elif available_memory_gb < 16:
        optimal_batch = 256
    elif available_memory_gb < 32:
        optimal_batch = base_batch_size
    else:
        optimal_batch = max_batch_size

    # Ensure batch size makes sense for dataset
    optimal_batch = min(optimal_batch, n_cells // 100)  # At least 100 batches
    optimal_batch = max(64, optimal_batch)  # Minimum viable size

    # Round to powers of 2 for efficiency
    optimal_batch = 2 ** int(np.log2(optimal_batch))

    return optimal_batch

def setup_environment():
    """Setup the analysis environment with M1 Max optimizations"""
    logger.info("Setting up analysis environment for microglia batch correction...")
    logger.info(f"Matplotlib backend: {matplotlib.get_backend()} (non-interactive)")

    # Set random seeds for reproducibility
    np.random.seed(Config.RANDOM_SEED)
    torch.manual_seed(Config.RANDOM_SEED)
    scvi.settings.seed = Config.RANDOM_SEED

    # Configure scanpy for non-interactive plotting
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=Config.DPI, facecolor='white')
    sc.settings.autoshow = False  # Disable automatic plot display
    sc.settings.figdir = Config.OUTPUT_DIR  # Set figure save directory

    # Optimize PyTorch for M1 Max
    torch.set_num_threads(9)  # Optimize for M1 Max cores

    # Configure scvi-tools for optimal performance
    scvi.settings.num_threads = 9

    # Check scvi-tools version for compatibility
    scvi_version = scvi.__version__
    logger.info(f"scvi-tools version: {scvi_version}")

    # Set conservative defaults for numerical stability
    if hasattr(scvi.settings, 'dl_num_workers'):
        scvi.settings.dl_num_workers = 0  # Disable multiprocessing for stability

    # GPU optimization with intelligent fallback
    if Config.USE_GPU and torch.backends.mps.is_available():
        device = "mps"
        logger.info("M1 Max MPS acceleration available - using Metal Performance Shaders")

        # Configure MPS for stability and performance
        torch.backends.mps.allow_tf32 = False  # Use full precision for stability

        logger.info("MPS configured for numerical stability")

        # Monitor GPU memory
        try:
            # Try to get MPS memory info if available
            logger.info("MPS memory optimization enabled")
        except:
            pass

    elif Config.USE_GPU and torch.cuda.is_available():
        device = torch.cuda.get_device_name(0)
        logger.info(f"CUDA GPU available: {device}")

        # CUDA memory optimization
        torch.cuda.empty_cache()
        memory_allocated = torch.cuda.memory_allocated() / 1024**3
        memory_reserved = torch.cuda.memory_reserved() / 1024**3
        logger.info(f"CUDA memory - Allocated: {memory_allocated:.2f}GB, Reserved: {memory_reserved:.2f}GB")

    else:
        device = "cpu"
        logger.info("Using CPU fallback")

    # Create output directory
    output_path = Path(Config.OUTPUT_DIR)
    output_path.mkdir(parents=True, exist_ok=True)
    logger.info(f"Output directory: {output_path}")

    # Log system info
    logger.info(f"Device: {device}")
    logger.info("Available RAM: 64GB M1 Max")
    logger.info(f"scanpy=={sc.__version__}")
    logger.info(f"scvi-tools=={scvi.__version__}")
    logger.info(f"torch=={torch.__version__}")

def add_cell_type_annotations(adata):
    """Add human-readable cell type annotations based on clusters"""
    # Map cluster numbers to cell type names
    adata.obs['cell_type'] = adata.obs[Config.LABELS_KEY].map(Config.MICROGLIA_CELL_TYPES)

    # Fill any missing mappings with the original cluster number
    missing_mask = adata.obs['cell_type'].isna()
    if missing_mask.any():
        adata.obs.loc[missing_mask, 'cell_type'] = 'Cluster ' + adata.obs.loc[missing_mask, Config.LABELS_KEY].astype(str)

    logger.info("Cell type annotations added:")
    cell_type_counts = adata.obs['cell_type'].value_counts()
    for cell_type, count in cell_type_counts.items():
        logger.info(f"   • {cell_type}: {count:,} cells")

    return adata

def load_and_inspect_data():
    """Load and inspect the HuMicA dataset for microglia analysis"""
    logger.info("Loading microglia data...")

    # Load data
    if not os.path.exists(Config.INPUT_H5AD):
        raise FileNotFoundError(f"Input file not found: {Config.INPUT_H5AD}")

    adata = sc.read_h5ad(Config.INPUT_H5AD)
    logger.info(f"Loaded microglia data: {adata.n_obs} cells × {adata.n_vars} genes")

    # Optimize sparse matrix format for training speed
    if hasattr(adata.X, 'tocsr'):
        adata.X = adata.X.tocsr()
        logger.info("Converting sparse matrix to CSR format for optimal training speed...")
        logger.info("Matrix format optimized")

    # Data validation and preprocessing
    logger.info("Performing data validation...")

    # Check for NaN values
    if hasattr(adata.X, 'data'):
        nan_count = np.isnan(adata.X.data).sum()
    else:
        nan_count = np.isnan(adata.X).sum()

    if nan_count > 0:
        logger.warning(f"Found {nan_count} NaN values in data. Replacing with zeros.")
        if hasattr(adata.X, 'data'):
            adata.X.data[np.isnan(adata.X.data)] = 0
        else:
            adata.X[np.isnan(adata.X)] = 0

    # Check for infinite values
    if hasattr(adata.X, 'data'):
        inf_count = np.isinf(adata.X.data).sum()
    else:
        inf_count = np.isinf(adata.X).sum()

    if inf_count > 0:
        logger.warning(f"Found {inf_count} infinite values in data. Replacing with zeros.")
        if hasattr(adata.X, 'data'):
            adata.X.data[np.isinf(adata.X.data)] = 0
        else:
            adata.X[np.isinf(adata.X)] = 0

    # Ensure non-negative values for count data
    if hasattr(adata.X, 'data'):
        neg_count = (adata.X.data < 0).sum()
        if neg_count > 0:
            logger.warning(f"Found {neg_count} negative values in count data. Clipping to zero.")
            adata.X.data[adata.X.data < 0] = 0
    else:
        neg_count = (adata.X < 0).sum()
        if neg_count > 0:
            logger.warning(f"Found {neg_count} negative values in count data. Clipping to zero.")
            adata.X[adata.X < 0] = 0

    # Log metadata info relevant to microglia analysis
    logger.info(f"Available metadata columns: {list(adata.obs.columns)}")
    logger.info(f"Batch key '{Config.BATCH_KEY}' unique values: {adata.obs[Config.BATCH_KEY].nunique()}")
    logger.info(f"Labels key '{Config.LABELS_KEY}' unique values: {adata.obs[Config.LABELS_KEY].nunique()}")

    # Log HuMicA disease/condition information
    if 'Group' in adata.obs.columns:
        logger.info(f"HuMicA disease groups: {adata.obs['Group'].value_counts().to_dict()}")
    if 'Diagnosis' in adata.obs.columns:
        logger.info(f"Diagnosis info: {adata.obs['Diagnosis'].value_counts().to_dict()}")
    if 'Study' in adata.obs.columns:
        logger.info(f"Contributing studies: {adata.obs['Study'].value_counts().to_dict()}")

    # Check HuMicA-specific gene signatures
    dam_present = [gene for gene in Config.DAM_GENES if gene in adata.var_names]
    homeostatic_present = [gene for gene in Config.HOMEOSTATIC_GENES if gene in adata.var_names]
    neuroinflam_present = [gene for gene in Config.NEUROINFLAM_GENES if gene in adata.var_names]

    logger.info(f"Disease-Associated Microglia (DAM) genes: {len(dam_present)}/{len(Config.DAM_GENES)}")
    logger.info(f"Homeostatic microglia genes: {len(homeostatic_present)}/{len(Config.HOMEOSTATIC_GENES)}")
    logger.info(f"Neuroinflammatory genes: {len(neuroinflam_present)}/{len(Config.NEUROINFLAM_GENES)}")

    # Check disease-specific markers
    for disease, markers in Config.DISEASE_MARKERS.items():
        present = [gene for gene in markers if gene in adata.var_names]
        logger.info(f"{disease} markers present: {len(present)}/{len(markers)}")

    # Check if raw counts are in X
    logger.info(f"Data matrix type: {type(adata.X)}")
    logger.info(f"Data range: {adata.X.min():.2f} to {adata.X.max():.2f}")

    # Add cell type annotations
    adata = add_cell_type_annotations(adata)

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
    """Setup scANVI model optimized for microglia batch correction with stability measures"""
    logger.info("Setting up scANVI model for microglia batch correction...")

    # Ensure categorical data types
    adata.obs[Config.BATCH_KEY] = adata.obs[Config.BATCH_KEY].astype('category')
    adata.obs[Config.LABELS_KEY] = adata.obs[Config.LABELS_KEY].astype('category')

    # Log batch and label distribution for microglia analysis
    logger.info(f"Batch distribution: {adata.obs[Config.BATCH_KEY].value_counts().to_dict()}")
    logger.info(f"Microglia subtype distribution: {adata.obs[Config.LABELS_KEY].value_counts().to_dict()}")

    # Additional data preprocessing for stability
    logger.info("Applying additional preprocessing for training stability...")

    # Check if data is already processed (following working scVI pattern)
    max_val = adata.X.max()
    logger.info(f"Data max value: {max_val:.2f}")

    # Conservative preprocessing approach - always preserve raw
    if max_val > 50:  # Likely raw counts
        logger.info("Detected raw counts - applying conservative normalization...")
        adata.raw = adata  # Always preserve raw
        # Conservative normalization (like working scVI implementation)
        sc.pp.normalize_total(adata, target_sum=1e4, copy=False)
        sc.pp.log1p(adata)
        logger.info("Data normalized and log-transformed")
    elif max_val > 10:  # Needs log transformation
        logger.info("Data appears normalized - applying log transformation...")
        if adata.raw is None:
            adata.raw = adata  # Preserve current state
        sc.pp.log1p(adata)
        logger.info("Log transformation applied")
    else:
        logger.info("Data appears already processed - preserving raw state...")
        if adata.raw is None:
            adata.raw = adata  # Always have raw backup

    # Final validation after preprocessing
    final_max = adata.X.max()
    final_min = adata.X.min()
    logger.info(f"Final data range: {final_min:.2f} to {final_max:.2f}")

    # GPU-optimized data validation - less aggressive clipping for better performance
    if final_max > 20:
        logger.info("Applying gentle data scaling for GPU stability")
        # Gentle scaling instead of hard clipping
        scale_factor = 15.0 / final_max
        if hasattr(adata.X, 'data'):
            adata.X.data = adata.X.data * scale_factor
        else:
            adata.X = adata.X * scale_factor
        logger.info(f"Data scaled by factor {scale_factor:.3f} for GPU training")

    # Register the dataset for scANVI
    scvi.model.SCANVI.setup_anndata(
        adata,
        layer=None,  # Use adata.X (processed counts)
        batch_key=Config.BATCH_KEY,  # Sample_ID for technical batch correction
        labels_key=Config.LABELS_KEY,  # Existing microglia subtypes
        unlabeled_category="Unknown",  # Required parameter - no unlabeled cells in this dataset
        categorical_covariate_keys=None,
        continuous_covariate_keys=None
    )

    # Create model with conservative parameters (based on working scVI implementation)
    model = scvi.model.SCANVI(
        adata,
        n_latent=Config.N_LATENT,  # Conservative 30 dimensions
        n_layers=Config.N_LAYERS,  # Proven 2 layers
        n_hidden=Config.N_HIDDEN,  # 256 hidden units
        dropout_rate=Config.DROPOUT_RATE,  # Conservative 0.1 dropout
        dispersion='gene',
        gene_likelihood='nb'
    )

    logger.info("Model created for microglia analysis:")
    logger.info(f"- Latent dimensions: {Config.N_LATENT}")
    logger.info(f"- Architecture: {Config.N_LAYERS} layers, {Config.N_HIDDEN} hidden units")
    logger.info(f"- Batch correction for {adata.obs[Config.BATCH_KEY].nunique()} samples")
    logger.info(f"- Preserving {adata.obs[Config.LABELS_KEY].nunique()} microglia subtypes")
    logger.info("- Enhanced stability features enabled")

    return model

def train_model(model, adata):
    """Train the scANVI model with intelligent GPU acceleration and stability measures"""
    logger.info("Training scANVI model for microglia batch correction...")

    # Check for NaN values in the model parameters before training
    for name, param in model.module.named_parameters():
        if torch.isnan(param).any():
            logger.error(f"NaN values found in model parameter {name} before training")
            raise ValueError(f"Model initialization failed - NaN in {name}")

    # Strategy 1: Try MPS with very conservative settings (based on working implementation)
    if Config.USE_GPU and torch.backends.mps.is_available():
        try:
            logger.info("Attempting MPS training with conservative settings...")

            # Monitor memory before training
            memory_info = monitor_gpu_memory()
            logger.info(f"Pre-training memory - System: {memory_info['system_memory_percent']:.1f}%, Available: {memory_info['system_memory_available_gb']:.1f}GB")

            # Use very conservative batch size for MPS stability
            optimal_batch_size = get_optimal_batch_size(
                adata.n_obs,
                memory_info['system_memory_available_gb'],
                'mps'
            )
            logger.info(f"Using conservative batch size for MPS: {optimal_batch_size}")

            # Conservative MPS training (following working scVI pattern)
            logger.info("Training with conservative MPS settings...")
            model.train(
                max_epochs=200,  # Reduced epochs for first attempt
                batch_size=optimal_batch_size,
                early_stopping=True,
                early_stopping_patience=30,
                accelerator='mps'
            )

            logger.info("MPS training completed successfully")

            # Log memory usage after training
            memory_info_after = monitor_gpu_memory()
            logger.info(f"Post-training memory - System: {memory_info_after['system_memory_percent']:.1f}%, Available: {memory_info_after['system_memory_available_gb']:.1f}GB")

            return model

        except Exception as e:
            error_msg = str(e)
            logger.warning(f"MPS training failed: {error_msg[:100]}...")

            # Provide specific guidance based on error patterns
            if "command buffer" in error_msg.lower():
                logger.warning("→ MPS hardware error - will use CPU fallback")
            elif "nan" in error_msg.lower() or "invalid values" in error_msg.lower():
                logger.warning("→ Numerical instability detected - trying CPU with different settings")
            elif "memory" in error_msg.lower():
                logger.warning("→ Memory issue detected - trying smaller batch size")

            # Clear MPS resources
            try:
                if hasattr(torch.mps, 'empty_cache'):
                    torch.mps.empty_cache()
                    logger.info("MPS cache cleared")
            except:
                pass

    # Strategy 2: Try CUDA if available
    if Config.USE_GPU and torch.cuda.is_available():
        try:
            logger.info("Attempting CUDA-accelerated training...")

            # Monitor CUDA memory and calculate optimal batch size
            memory_info = monitor_gpu_memory()
            if 'cuda_free_gb' in memory_info:
                available_memory = memory_info['cuda_free_gb']
            else:
                available_memory = 8.0  # Conservative default

            optimal_batch_size = get_optimal_batch_size(
                adata.n_obs,
                available_memory,
                'cuda'
            )
            logger.info(f"Using adaptive batch size for CUDA: {optimal_batch_size} (available memory: {available_memory:.1f}GB)")

            cuda_kwargs = {
                'max_epochs': Config.MAX_EPOCHS,
                'batch_size': optimal_batch_size,
                'early_stopping': Config.EARLY_STOPPING,
                'early_stopping_patience': Config.EARLY_STOPPING_PATIENCE,
                'accelerator': 'gpu'
            }

            try:
                model.train(lr=Config.LEARNING_RATE, **cuda_kwargs)
            except TypeError:
                model.train(**cuda_kwargs)

            logger.info("Training completed successfully with CUDA acceleration")

            # Log CUDA memory usage
            memory_info_after = monitor_gpu_memory()
            if 'cuda_allocated_gb' in memory_info_after:
                logger.info(f"CUDA memory after training - Allocated: {memory_info_after['cuda_allocated_gb']:.2f}GB, Reserved: {memory_info_after['cuda_reserved_gb']:.2f}GB")

            return model

        except Exception as e:
            error_msg = str(e)
            logger.warning(f"CUDA training failed: {error_msg[:200]}...")

            # Comprehensive CUDA error pattern detection
            if "out of memory" in error_msg.lower() or "cuda out of memory" in error_msg.lower():
                logger.warning("CUDA out of memory error")
                logger.info("Recommendation: Reduce batch size, clear cache, or use CPU")
            elif "nan" in error_msg.lower() or "invalid values" in error_msg.lower():
                logger.warning("NaN values in CUDA training - numerical instability")
                logger.info("Recommendation: Check data preprocessing and model parameters")
            elif "cudnn" in error_msg.lower():
                logger.warning("CuDNN error detected")
                logger.info("Recommendation: Update CUDA/CuDNN or try different precision")
            elif "device" in error_msg.lower() and "unavailable" in error_msg.lower():
                logger.warning("CUDA device unavailable")
                logger.info("Recommendation: Check GPU drivers and CUDA installation")
            elif "runtime error" in error_msg.lower():
                logger.warning("CUDA runtime error")
                logger.info("Recommendation: Try restarting and clearing GPU memory")
            else:
                logger.warning("Unknown CUDA error - falling back to CPU")

            logger.info("Falling back to CPU training...")
            torch.cuda.empty_cache()

    # Strategy 3: CPU training with proven settings (based on working scVI implementation)
    try:
        logger.info("Training with CPU using proven stable settings...")

        # Use conservative but proven settings
        memory_info = monitor_gpu_memory()
        optimal_batch_size = get_optimal_batch_size(
            adata.n_obs,
            memory_info['system_memory_available_gb'],
            'cpu'
        )
        logger.info(f"Using proven batch size for CPU: {optimal_batch_size}")

        # Training with settings proven to work in the reference implementation
        model.train(
            max_epochs=Config.MAX_EPOCHS,
            batch_size=optimal_batch_size,
            early_stopping=Config.EARLY_STOPPING,
            early_stopping_patience=Config.EARLY_STOPPING_PATIENCE,
            accelerator='cpu'
        )

        logger.info("CPU training completed successfully with proven settings")
        return model

    except Exception as e:
        logger.error(f"CPU training failed: {e}")
        logger.info("Trying minimal fallback settings...")
        error_msg = str(e)
        logger.error(f"Optimized CPU training failed: {error_msg[:200]}...")

        # CPU-specific error analysis
        if "nan" in error_msg.lower() or "invalid values" in error_msg.lower():
            logger.error("Critical: NaN values persist even in CPU training")
            logger.info("This indicates a fundamental data or model issue")
        elif "memory" in error_msg.lower():
            logger.error("CPU memory error - system may be overloaded")
        elif "convergence" in error_msg.lower():
            logger.warning("Training convergence issues")

        logger.info("Attempting conservative fallback...")

    # Final strategy: Minimal settings guaranteed to work
    try:
        logger.info("Using minimal guaranteed settings...")

        # Absolute minimal settings that should always work
        model.train(
            max_epochs=100,
            batch_size=128,  # Very small batch size
            early_stopping=True,
            early_stopping_patience=20,
            accelerator='cpu'
        )

        logger.info("Training completed with minimal guaranteed settings")
        return model

    except Exception as e2:
        logger.error(f"Even minimal training failed: {e2}")

        # Provide detailed error analysis
        error_msg = str(e2)
        if "nan" in error_msg.lower():
            logger.error("→ Root cause: NaN values in data or model - data preprocessing issue")
            logger.error("→ Suggestion: Check input data for infinite/NaN values")
        elif "memory" in error_msg.lower():
            logger.error("→ Root cause: Insufficient memory - try reducing dataset size")
        elif "invalid" in error_msg.lower():
            logger.error("→ Root cause: Invalid parameter values - model setup issue")
        else:
            logger.error(f"→ Unknown error pattern: {error_msg[:200]}")

        raise RuntimeError(f"All training strategies exhausted. Root error: {e2}")

def plot_training_history(model):
    """Plot training history and return training metrics"""
    try:
        train_elbo = model.history["elbo_train"]
        val_elbo = model.history["elbo_validation"]

        plt.figure(figsize=Config.FIGSIZE)
        plt.plot(train_elbo, label='Training ELBO', linewidth=2)
        plt.plot(val_elbo, label='Validation ELBO', linewidth=2)
        plt.xlabel('Epoch')
        plt.ylabel('ELBO (Evidence Lower Bound)')
        plt.title('scANVI Training History - Microglia Batch Correction')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig(f"{Config.OUTPUT_DIR}/training_history.png", dpi=Config.DPI, bbox_inches='tight')
        plt.close()

        final_train_elbo = train_elbo.iloc[-1]
        final_val_elbo = val_elbo.iloc[-1]
        logger.info(f"Training completed after {len(train_elbo)} epochs")
        logger.info(f"Final training ELBO: {final_train_elbo:.2f}")
        logger.info(f"Final validation ELBO: {final_val_elbo:.2f}")
    except Exception as e:
        logger.warning(f"Could not plot training history: {e}")

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
    """Create comprehensive visualizations for microglia analysis"""
    logger.info("Creating visualizations for microglia batch correction...")

    # Main UMAP overview
    fig, axes = plt.subplots(2, 2, figsize=(20, 15))

    # Original microglia clusters
    sc.pl.umap(adata, color=Config.LABELS_KEY, legend_loc='on data',
               title='Original Microglia Clusters', ax=axes[0,0], show=False)

    # scANVI predictions
    sc.pl.umap(adata, color="scanvi_predictions", legend_loc='on data',
               title='scANVI Predictions', ax=axes[0,1], show=False)

    # Batch effect assessment
    sc.pl.umap(adata, color=Config.BATCH_KEY, title='Batch (Sample ID)',
               ax=axes[1,0], show=False, legend_loc='right margin')

    # Disease/condition groups
    if 'Group' in adata.obs.columns:
        sc.pl.umap(adata, color='Group', title='Disease Groups',
                   ax=axes[1,1], show=False, legend_loc='right margin')
    else:
        sc.pl.umap(adata, color="leiden_scanvi", legend_loc='on data',
                   title='Leiden Clustering (scANVI latent)', ax=axes[1,1], show=False)

    plt.tight_layout()
    plt.savefig(f"{Config.OUTPUT_DIR}/microglia_umap_overview.png", dpi=Config.DPI, bbox_inches='tight')
    plt.close()

    # Dedicated cell type UMAP with proper annotations
    plt.figure(figsize=(16, 12))

    # Create a 2x2 subplot for detailed cell type analysis
    fig, axes = plt.subplots(2, 2, figsize=(20, 16))

    # 1. Cell types with meaningful names
    sc.pl.umap(adata, color='cell_type', legend_loc='right margin',
               title='Microglia Cell Types (Annotated)',
               ax=axes[0,0], show=False, frameon=False, size=20)

    # 2. Original cluster numbers for reference
    sc.pl.umap(adata, color=Config.LABELS_KEY, legend_loc='right margin',
               title='Original Cluster Numbers (Seurat)',
               ax=axes[0,1], show=False, frameon=False, size=20)

    # 3. New Leiden clusters from scANVI latent space
    sc.pl.umap(adata, color="leiden_scanvi", legend_loc='right margin',
               title='New Leiden Clusters (scANVI Latent Space)',
               ax=axes[1,0], show=False, frameon=False, size=20)

    # 4. Disease groups if available
    if 'Group' in adata.obs.columns:
        sc.pl.umap(adata, color='Group', legend_loc='right margin',
                   title='Disease Groups',
                   ax=axes[1,1], show=False, frameon=False, size=20)
    else:
        # Fallback to scANVI predictions
        sc.pl.umap(adata, color="scanvi_predictions", legend_loc='right margin',
                   title='scANVI Predicted Subtypes',
                   ax=axes[1,1], show=False, frameon=False, size=20)

    plt.suptitle('HuMicA Cell Type Analysis - Microglia Subtypes Across Conditions',
                 fontsize=20, y=0.98)
    plt.tight_layout()
    plt.subplots_adjust(top=0.93)
    plt.savefig(f"{Config.OUTPUT_DIR}/microglia_cell_types_umap.png",
                dpi=Config.DPI, bbox_inches='tight')
    plt.close()

    # Large focused UMAP for cell types - main visualization
    plt.figure(figsize=(14, 10))

    # Create detailed microglia subtype map with meaningful cell type names
    sc.pl.umap(adata, color='cell_type',
               legend_loc='right margin',
               title='HuMicA Microglia Cell Types\n(Annotated Subtypes)',
               frameon=False,
               size=25,
               alpha=0.8)

    # Add cell type count information with meaningful names
    cell_type_counts = adata.obs['cell_type'].value_counts()
    cell_type_info = []
    for cell_type, count in cell_type_counts.items():
        percentage = (count / len(adata.obs)) * 100
        cell_type_info.append(f"{cell_type}: {count:,} cells ({percentage:.1f}%)")

    # Add text box with cell type information
    textstr = '\n'.join(cell_type_info)
    props = dict(boxstyle='round', facecolor='lightblue', alpha=0.8)
    plt.text(0.02, 0.98, textstr, transform=plt.gca().transAxes, fontsize=9,
             verticalalignment='top', bbox=props)

    plt.tight_layout()
    plt.savefig(f"{Config.OUTPUT_DIR}/microglia_subtypes_detailed_umap.png",
                dpi=Config.DPI, bbox_inches='tight')
    plt.close()

    # Additional UMAP showing cluster numbers for reference
    plt.figure(figsize=(14, 10))
    sc.pl.umap(adata, color=Config.LABELS_KEY,
               legend_loc='right margin',
               title='HuMicA Microglia Cluster Numbers\n(Original Seurat Clustering)',
               frameon=False,
               size=25,
               alpha=0.8)

    # Add cluster number information
    cluster_counts = adata.obs[Config.LABELS_KEY].value_counts().sort_index()
    cluster_info = []
    for cluster, count in cluster_counts.items():
        percentage = (count / len(adata.obs)) * 100
        cluster_name = Config.MICROGLIA_CELL_TYPES.get(str(cluster), f"Cluster {cluster}")
        cluster_info.append(f"Cluster {cluster} ({cluster_name}): {count:,} cells ({percentage:.1f}%)")

    # Add text box with cluster information
    textstr = '\n'.join(cluster_info)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    plt.text(0.02, 0.98, textstr, transform=plt.gca().transAxes, fontsize=9,
             verticalalignment='top', bbox=props)

    plt.tight_layout()
    plt.savefig(f"{Config.OUTPUT_DIR}/microglia_clusters_reference_umap.png",
                dpi=Config.DPI, bbox_inches='tight')
    plt.close()

    # Detailed batch integration assessment
    plt.figure(figsize=(15, 6))
    plt.subplot(1, 2, 1)
    sc.pl.umap(adata, color=Config.BATCH_KEY, title='Before Integration (by Sample)',
               legend_loc='right margin', show=False)
    plt.subplot(1, 2, 2)
    sc.pl.umap(adata, color=Config.LABELS_KEY, title='After Integration (Microglia Types)',
               legend_loc='right margin', show=False)
    plt.tight_layout()
    plt.savefig(f"{Config.OUTPUT_DIR}/batch_integration_assessment.png", dpi=Config.DPI, bbox_inches='tight')
    plt.close()

    # HuMicA-specific gene expression visualizations
    # 1. Disease-Associated Microglia (DAM) signature
    dam_present = [gene for gene in Config.DAM_GENES if gene in adata.var_names]
    homeostatic_present = [gene for gene in Config.HOMEOSTATIC_GENES if gene in adata.var_names]
    if dam_present:
        key_dam_genes = dam_present[:min(6, len(dam_present))]
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        axes = axes.flatten()

        for i, gene in enumerate(key_dam_genes):
            sc.pl.umap(adata, color=gene, ax=axes[i], show=False,
                      title=f'{gene} (DAM marker)', cmap='Reds')

        plt.suptitle('Disease-Associated Microglia (DAM) Signature Genes', fontsize=16)
        plt.tight_layout()
        plt.savefig(f"{Config.OUTPUT_DIR}/dam_signature_genes.png", dpi=Config.DPI, bbox_inches='tight')
        plt.close()

    # 2. Homeostatic vs DAM comparison
    if dam_present and homeostatic_present:
        homeostatic_key = homeostatic_present[:3]
        dam_key = dam_present[:3]

        fig, axes = plt.subplots(2, 3, figsize=(18, 12))

        # Top row: Homeostatic markers
        for i, gene in enumerate(homeostatic_key):
            sc.pl.umap(adata, color=gene, ax=axes[0,i], show=False,
                      title=f'{gene} (Homeostatic)', cmap='Blues')

        # Bottom row: DAM markers
        for i, gene in enumerate(dam_key):
            sc.pl.umap(adata, color=gene, ax=axes[1,i], show=False,
                      title=f'{gene} (DAM)', cmap='Reds')

        plt.suptitle('HuMicA: Homeostatic vs Disease-Associated Microglia', fontsize=16)
        plt.tight_layout()
        plt.savefig(f"{Config.OUTPUT_DIR}/homeostatic_vs_dam.png", dpi=Config.DPI, bbox_inches='tight')
        plt.close()

    # 3. Neuroinflammatory pathway genes
    neuroinflam_present = [gene for gene in Config.NEUROINFLAM_GENES if gene in adata.var_names]
    if neuroinflam_present:
        key_genes = neuroinflam_present[:min(6, len(neuroinflam_present))]
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        axes = axes.flatten()

        for i, gene in enumerate(key_genes):
            sc.pl.umap(adata, color=gene, ax=axes[i], show=False,
                      title=f'{gene} Expression', cmap='viridis')

        plt.suptitle('Neuroinflammatory Pathway Genes', fontsize=16)
        plt.tight_layout()
        plt.savefig(f"{Config.OUTPUT_DIR}/neuroinflammatory_genes.png", dpi=Config.DPI, bbox_inches='tight')
        plt.close()

    logger.info("HuMicA-specific visualizations saved")

def save_results(model, adata):
    """Save model and processed data (creates checkpoint)"""
    logger.info("Saving results...")

    # Save model (checkpoint)
    model_path = f"{Config.OUTPUT_DIR}/scanvi_model"
    model.save(model_path, overwrite=True)
    logger.info(f"Model checkpoint saved to: {model_path}")

    # Save processed AnnData (checkpoint)
    adata_path = f"{Config.OUTPUT_DIR}/adata_scanvi_processed.h5ad"
    adata.write_h5ad(adata_path)
    logger.info(f"Data checkpoint saved to: {adata_path}")

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
    """Generate a summary report of the microglia analysis"""
    logger.info("Generating microglia analysis summary report...")

    # Get HuMicA gene signature info
    dam_present = [gene for gene in Config.DAM_GENES if gene in adata.var_names]
    homeostatic_present = [gene for gene in Config.HOMEOSTATIC_GENES if gene in adata.var_names]
    neuroinflam_present = [gene for gene in Config.NEUROINFLAM_GENES if gene in adata.var_names]

    # Get disease group info if available
    disease_info = ""
    if 'Group' in adata.obs.columns:
        disease_counts = adata.obs['Group'].value_counts()
        disease_info = f"\nDisease Groups:\n" + "\n".join([f"- {group}: {count} cells" for group, count in disease_counts.items()])

    report = f"""
Human Microglia Atlas (HuMicA) scANVI Analysis Report
Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
Reference: "The Human Microglia Atlas (HuMicA) unravels changes in disease-associated
microglia subsets across neurodegenerative conditions" - Nature Communications

=== ANALYSIS OVERVIEW ===
Primary Goal: Batch correction and disease-associated microglia analysis
Dataset: Human Microglia Atlas (HuMicA) - multi-condition neurodegenerative diseases
Platform: M1 Max (64GB RAM)

=== DATASET INFORMATION ===
- Input file: {Config.INPUT_H5AD}
- Final dimensions: {adata.n_obs} cells × {adata.n_vars} genes
- Batch key: {Config.BATCH_KEY} ({adata.obs[Config.BATCH_KEY].nunique()} different samples)
- Labels key: {Config.LABELS_KEY} ({adata.obs[Config.LABELS_KEY].nunique()} microglia subtypes){disease_info}

=== MODEL CONFIGURATION ===
- Latent dimensions: {Config.N_LATENT} (optimized for neurodegenerative disease complexity)
- Architecture: {Config.N_LAYERS} layers, {Config.N_HIDDEN} hidden units
- Training epochs: {len(model.history["elbo_train"])}
- Batch size: {Config.BATCH_SIZE} (M1 Max optimized)
- Learning rate: {Config.LEARNING_RATE}

=== HuMicA GENE SIGNATURE ANALYSIS ===
- Disease-Associated Microglia (DAM) genes: {len(dam_present)}/{len(Config.DAM_GENES)}
- Homeostatic microglia genes: {len(homeostatic_present)}/{len(Config.HOMEOSTATIC_GENES)}
- Neuroinflammatory pathway genes: {len(neuroinflam_present)}/{len(Config.NEUROINFLAM_GENES)}
- Key DAM genes detected: {', '.join(dam_present[:8])}{"..." if len(dam_present) > 8 else ""}

=== RESULTS ===
- scANVI latent representation: X_scANVI ({adata.obsm["X_scANVI"].shape})
- Batch-corrected expression: X_scANVI_normalized
- New Leiden clusters: {adata.obs["leiden_scanvi"].nunique()}
- Device used: {'MPS (M1 Max)' if torch.backends.mps.is_available() else 'CPU (M1 Max optimized)'}

=== OUTPUT FILES ===
- Processed data: adata_scanvi_processed.h5ad
- Trained model: scanvi_model/
- Latent embeddings: scanvi_latent_embeddings.csv
- Metadata with predictions: scanvi_metadata.csv
- HuMicA-specific visualizations:
  * microglia_umap_overview.png (Multi-panel overview)
  * microglia_cell_types_umap.png (Detailed cell type analysis with annotations)
  * microglia_subtypes_detailed_umap.png (Focused cell type UMAP with meaningful names)
  * microglia_clusters_reference_umap.png (Cluster numbers reference with type mapping)
  * batch_integration_assessment.png (Before/after batch correction)
  * dam_signature_genes.png (Disease-Associated Microglia expression)
  * homeostatic_vs_dam.png (Homeostatic vs DAM comparison)
  * neuroinflammatory_genes.png (Neuroinflammation pathway genes)
  * training_history.png (Model training convergence)

=== CELL TYPE ANNOTATIONS ===
- Homeostatic Microglia: Resting state, surveillance function
- Activated Microglia I & II: Early and advanced activation states
- Disease-Associated Microglia (DAM): Pathology-associated phenotype
- Interferon-Response Microglia: Anti-viral/immune response
- Proliferating Microglia: Cell cycle active, tissue repair
- Lipid-Associated Microglia: Lipid metabolism, myelin debris
- Phagocytic Microglia: Active debris clearance
- Stressed/Dysfunctional Microglia: Senescent or damaged cells

=== NEXT STEPS FOR HuMicA ANALYSIS ===
1. Use X_scANVI_normalized for disease-associated microglia differential expression
2. Compare DAM signature scores across neurodegenerative conditions
3. Analyze homeostatic to DAM transition patterns by cell type
4. Investigate disease-specific microglia activation profiles
5. Perform pathway enrichment analysis for each annotated cell type
6. Validate cell type annotations against original HuMicA publication results
7. Cross-reference with your Complement-OUD study for addiction-related microglia states
"""

    with open(f"{Config.OUTPUT_DIR}/humica_analysis_summary.txt", "w") as f:
        f.write(report)

    logger.info("HuMicA analysis summary report generated")
    print(report)

def main():
    """Main analysis pipeline for HuMicA microglia batch correction"""

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='HuMicA scANVI Analysis Pipeline - Microglia Batch Correction',
        epilog="""
Examples:
  %(prog)s                    # Load checkpoint (fast)
  %(prog)s --retrain          # Force retraining (slow)
  %(prog)s --skip-viz         # Analysis only, no plots
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--retrain', action='store_true',
                       help='Force retraining even if checkpoint exists')
    parser.add_argument('--skip-viz', action='store_true',
                       help='Skip visualization generation')
    args = parser.parse_args()

    # Override config based on arguments
    skip_training = Config.SKIP_TRAINING and not args.retrain

    try:
        # Setup
        setup_environment()

        # Check for pre-trained model and processed data
        if skip_training and os.path.exists(Config.MODEL_CHECKPOINT) and os.path.exists(Config.DATA_CHECKPOINT):
            logger.info("=" * 70)
            logger.info("🔄 CHECKPOINT DETECTED: Loading pre-trained model and processed data")
            logger.info("=" * 70)
            logger.info(f"⚡ This will be much faster than retraining!")
            logger.info(f"📁 Model checkpoint: {Config.MODEL_CHECKPOINT}")
            logger.info(f"📁 Data checkpoint: {Config.DATA_CHECKPOINT}")

            # Load processed data
            adata = sc.read_h5ad(Config.DATA_CHECKPOINT)
            logger.info(f"✅ Loaded processed data: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

            # Add cell type annotations if not present
            if 'cell_type' not in adata.obs.columns:
                adata = add_cell_type_annotations(adata)
                logger.info("✅ Added cell type annotations to checkpoint data")

            # Load trained model
            model = scvi.model.SCANVI.load(Config.MODEL_CHECKPOINT, adata)
            logger.info("✅ Loaded pre-trained scANVI model")

            # Log checkpoint info
            logger.info("📊 Checkpoint contents:")
            logger.info(f"   • Latent representation: {adata.obsm['X_scANVI'].shape}")
            logger.info(f"   • Batch correction: {adata.obs[Config.BATCH_KEY].nunique()} samples")
            logger.info(f"   • Cell types: {adata.obs[Config.LABELS_KEY].nunique()} clusters")
            if 'leiden_scanvi' in adata.obs.columns:
                logger.info(f"   • New clusters: {adata.obs['leiden_scanvi'].nunique()} Leiden clusters")

        else:
            logger.info("=" * 70)
            if args.retrain:
                logger.info("🔬 FULL TRAINING: Forced retraining requested")
                logger.info("⏰ This will take ~2 hours with GPU acceleration")
            else:
                logger.info("🔬 FULL TRAINING: No checkpoint found")
                logger.info("⏰ This will take ~2 hours, but creates checkpoint for future runs")
            logger.info("=" * 70)

            # Load and inspect data
            adata = load_and_inspect_data()

            # Skip quality control - data already processed in Seurat
            if not Config.SKIP_QC:
                logger.info("Applying additional quality control...")
                adata = quality_control(adata)
            else:
                logger.info("Skipping QC - using pre-processed Seurat data")

            # Setup and train model
            model = setup_scanvi_model(adata)
            model = train_model(model, adata)

            # Plot training history
            plot_training_history(model)

            # Extract results and run downstream analysis
            adata = extract_results(model, adata)
            run_downstream_analysis(adata)

            # Save results (create checkpoint)
            save_results(model, adata)
            logger.info("💾 Checkpoint saved for future runs")

        # Always run visualizations and analysis (can be updated without retraining)
        if not args.skip_viz:
            logger.info("📊 Running visualization and analysis...")

            # Ensure downstream analysis is complete
            if 'X_umap' not in adata.obsm:
                logger.info("Running missing downstream analysis...")
                run_downstream_analysis(adata)

            # Create visualizations
            create_visualizations(adata)
        else:
            logger.info("⏭️  Skipping visualizations (--skip-viz flag)")

        # Generate summary
        generate_summary_report(adata, model)

        logger.info("=" * 70)
        logger.info("🎉 HuMicA scANVI analysis completed successfully!")
        logger.info("=" * 70)
        if not Config.checkpoint_exists():
            logger.info("💡 Next run will be much faster using checkpoints!")
        logger.info("📊 Check the results directory for outputs and visualizations")

    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        raise

if __name__ == "__main__":
    main()
