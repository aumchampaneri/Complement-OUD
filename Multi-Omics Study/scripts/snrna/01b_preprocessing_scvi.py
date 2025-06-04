'''
🔬 Single-nucleus RNA-seq Processing Pipeline - scVI Version
GSE225158 - OUD Striatum snRNA-seq data reprocessing with scVI

Alternative to Pearson residuals using probabilistic modeling
- Better batch integration across 22 subjects
- Improved statistical framework for clinical comparisons
- State-of-the-art normalization and dimensionality reduction
'''

# =======================================
# 📁 INPUT/OUTPUT PATHS CONFIGURATION
# =======================================

# Input data path
INPUT_PATH = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study/data/raw/snrna/GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5ad"

# Output directory structure
BASE_OUTPUT_DIR = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study"
PROCESSED_DATA_DIR = f"{BASE_OUTPUT_DIR}/data/processed/snrna_scvi"
RESULTS_DIR = f"{BASE_OUTPUT_DIR}/results/snrna_scvi"
PLOTS_DIR = f"{RESULTS_DIR}/qc_plots"

# Output files for scVI version
OUTPUT_H5AD_SCVI = "GSE225158_reprocessed_scvi.h5ad"
OUTPUT_H5AD_PATH_SCVI = f"{PROCESSED_DATA_DIR}/{OUTPUT_H5AD_SCVI}"
SUMMARY_PATH_SCVI = f"{PROCESSED_DATA_DIR}/scvi_summary.txt"

# Checkpoint files
CHECKPOINT_DIR = f"{PROCESSED_DATA_DIR}/checkpoints"
CHECKPOINT_QC = f"{CHECKPOINT_DIR}/01_after_qc.h5ad"
CHECKPOINT_SCVI_MODEL = f"{CHECKPOINT_DIR}/02_scvi_model"
CHECKPOINT_ANALYSIS = f"{CHECKPOINT_DIR}/03_after_analysis.h5ad"
CHECKPOINT_CLUSTERING = f"{CHECKPOINT_DIR}/04_after_clustering.h5ad"

# =======================================
# 📚 IMPORT LIBRARIES
# =======================================

import scanpy as sc
import scvi
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend to prevent plot display
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import sparse
import warnings
import os
import torch
warnings.filterwarnings('ignore')

# Set scanpy to not show plots automatically
sc.settings.autoshow = False

# =======================================
# 🚀 M1 MAX OPTIMIZATIONS
# =======================================

# Configure for M1 Max with 64GB RAM
print("🚀 Configuring for M1 Max with 64GB RAM...")

# Set PyTorch to use Metal Performance Shaders (MPS) for M1
if torch.backends.mps.is_available():
    device = torch.device("mps")
    print(f"   • Using Metal Performance Shaders (MPS): {device}")
    # Don't set default device to avoid generator issues
    # torch.set_default_device(device)
else:
    device = torch.device("cpu")
    print(f"   • MPS not available, using CPU: {device}")

# Optimize memory usage for large dataset (with version compatibility)
try:
    torch.backends.mps.enable_fallback_warnings(False)  # Suppress MPS warnings
except AttributeError:
    print("   • Note: PyTorch version doesn't support MPS fallback warning control")

# Set memory optimization environment variable
os.environ['PYTORCH_MPS_HIGH_WATERMARK_RATIO'] = '0.0'  # Use available GPU memory

# Set scanpy settings optimized for large memory
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=300, facecolor='white')
sc.settings.n_jobs = -1  # Use all CPU cores
sc.settings.max_memory = 32  # Use up to 32GB for scanpy operations

# Configure scvi-tools for M1 Max
scvi.settings.seed = 42
scvi.settings.dl_pin_memory_gpu_training = False  # Better for MPS
scvi.settings.batch_size = 2048  # Larger batch size for your RAM

print("✅ M1 Max optimizations configured")

# =======================================
# 🔥 MEMORY USAGE OPTIMIZATIONS
# =======================================

# Reduce memory usage for large data processing
print("🔥 Applying memory usage optimizations...")

# Remove CUDA-specific settings that don't work with MPS
# torch.backends.cudnn.benchmark = True  # Not needed for MPS
# torch.backends.cudnn.enabled = True    # Not needed for MPS

# Optimize NumPy settings
np.set_printoptions(precision=4, suppress=True)  # Better float formatting

# Optimize Pandas settings
pd.options.mode.chained_assignment = None  # Disable SettingWithCopyWarning

print("✅ Memory usage optimizations applied")

# =======================================
# 🔬 CUSTOM FUNCTIONS
# =======================================

def calculate_mad_outliers(values, n_mads=3):
    """Calculate outliers using Median Absolute Deviation (MAD)"""
    median = np.median(values)
    mad = np.median(np.abs(values - median))
    lower = median - n_mads * mad
    upper = median + n_mads * mad
    return lower, upper

def save_checkpoint(adata, checkpoint_path, step_name):
    """Save checkpoint with error handling"""
    try:
        os.makedirs(os.path.dirname(checkpoint_path), exist_ok=True)
        adata.write(checkpoint_path)
        print(f"✅ Checkpoint saved: {step_name}")
        return True
    except Exception as e:
        print(f"⚠️  Failed to save checkpoint {step_name}: {e}")
        return False

def load_checkpoint(checkpoint_path, step_name):
    """Load checkpoint with error handling"""
    try:
        if os.path.exists(checkpoint_path):
            adata = sc.read_h5ad(checkpoint_path)
            print(f"✅ Checkpoint loaded: {step_name}")
            return adata
        else:
            print(f"📝 No checkpoint found for: {step_name}")
            return None
    except Exception as e:
        print(f"⚠️  Failed to load checkpoint {step_name}: {e}")
        return None

def main_scvi(resume_from_step=None):
    """scVI-based processing pipeline optimized for M1 Max with checkpointing
    
    Args:
        resume_from_step (str): Step to resume from ('qc', 'model', 'analysis', 'clustering', None)
    """
    
    print("=" * 70)
    print("🔬 scVI-BASED PREPROCESSING PIPELINE (M1 MAX OPTIMIZED)")
    print("=" * 70)
    
    if resume_from_step:
        print(f"🔄 Resuming from step: {resume_from_step}")
    
    # Create checkpoint directory
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    
    # =======================================
    # 🔬 STEP 1: LOAD DATA & EXTRACT RAW COUNTS
    # =======================================
    
    # Check if we can skip QC step
    if resume_from_step in ['model', 'analysis', 'clustering']:
        print("🔄 RESUMING: Loading QC checkpoint...")
        raw_adata = load_checkpoint(CHECKPOINT_QC, "Quality Control")
        if raw_adata is not None:
            print(f"✅ Resumed from QC: {raw_adata.n_obs:,} cells, {raw_adata.n_vars:,} genes")
            
            # 🧬 CRITICAL FIX: Ensure gene names are properly set even when resuming
            print("🔧 Verifying gene names after checkpoint loading...")
            if 'features' in raw_adata.var.columns:
                # Check if gene names are numeric (indicating they need fixing)
                if raw_adata.var_names[0].isdigit():
                    print("   ⚠️  Gene names are numeric indices, fixing...")
                    raw_adata.var_names = raw_adata.var['features'].values
                    raw_adata.var_names_make_unique()  # Ensure unique names
                    print(f"   ✅ Gene names fixed from features column")
                    print(f"   📋 Sample gene names: {raw_adata.var_names[:5].tolist()}")
                else:
                    print(f"   ✅ Gene names already correct: {raw_adata.var_names[:5].tolist()}")
            else:
                print("   ⚠️  Warning: 'features' column not found in var")
                print(f"   📋 Available var columns: {raw_adata.var.columns.tolist()}")
            
            goto_step = resume_from_step
        else:
            print("❌ QC checkpoint not found, starting from beginning...")
            goto_step = None
    else:
        goto_step = None
    
    if goto_step is None:
        print("🔬 STEP 1: DATA LOADING & RAW EXTRACTION")
        print("=" * 70)
        
        # Load and extract raw data (same as main pipeline)
        adata = sc.read_h5ad(INPUT_PATH)
        print(f"📊 Loaded dataset: {adata.shape}")
        
        # Extract raw counts
        if adata.raw is not None:
            raw_adata = adata.raw.to_adata()
            
            # 🧬 CRITICAL FIX: Set gene names from 'features' column as var index
            print("🔧 Setting gene names from 'features' column as var index...")
            if 'features' in raw_adata.var.columns:
                # Set the gene names as the index
                raw_adata.var_names = raw_adata.var['features'].values
                raw_adata.var_names_make_unique()  # Ensure unique names
                print(f"✅ Gene names set from features column")
                print(f"📋 Sample gene names: {raw_adata.var_names[:5].tolist()}")
            else:
                print("⚠️  Warning: 'features' column not found in var")
                print(f"📋 Available var columns: {raw_adata.var.columns.tolist()}")
            
            # Copy clinical metadata - FOCUSED LIST
            clinical_cols = [
                # Demographics & General (essential)
                'ID', 'Sex', 'Race', 'Age', 'BMI', 'PMI', 'pH', 'RIN',
                
                # Diagnosis & Mental Health (core for OUD study)
                'Dx_OUD', 'Dx_Substances', 'Dx_Comorbid',
                
                # Brain/Anatomical (essential for striatum study)
                'Region', 'Case',
                
                # Cell Type Annotations (analysis-derived but useful)
                'celltype1', 'celltype2', 'celltype3'
            ]
            available_clinical = [col for col in clinical_cols if col in adata.obs.columns]
            for col in available_clinical:
                raw_adata.obs[col] = adata.obs[col].copy()
            
            print(f"✅ Copied {len(available_clinical)} clinical metadata columns")
            
            # Add specific check for brain region
            if 'Region' in raw_adata.obs.columns:
                unique_regions = raw_adata.obs['Region'].unique()
                print(f"🧠 Brain regions identified: {list(unique_regions)}")
                print(f"   • Cells per region:")
                region_counts = raw_adata.obs['Region'].value_counts()
                for region, count in region_counts.items():
                    print(f"     - {region}: {count:,} cells")
            else:
                print("⚠️  Warning: No brain region information found")
            
            del adata
            print(f"✅ Extracted raw data: {raw_adata.shape}")
        else:
            print("❌ No raw data found")
            return None
        
        # =======================================
        # 🔬 STEP 2: QUALITY CONTROL (Same as main pipeline)
        # =======================================
        print("\n🔬 STEP 2: QUALITY CONTROL")
        print("=" * 70)
        
        # Calculate gene annotations
        raw_adata.var['mt'] = raw_adata.var_names.str.startswith('MT-')
        raw_adata.var['ribo'] = raw_adata.var_names.str.startswith(('RPS', 'RPL'))
        
        # Calculate QC metrics with proper qc_vars specification
        sc.pp.calculate_qc_metrics(
            raw_adata, 
            qc_vars=['mt', 'ribo'], 
            percent_top=None, 
            log1p=False, 
            inplace=True
        )
        
        # Add percentage metrics with optimized column handling (no debug prints)
        mt_counts_col = 'total_counts_mt' if 'total_counts_mt' in raw_adata.obs.columns else 'n_counts_mt'
        ribo_counts_col = 'total_counts_ribo' if 'total_counts_ribo' in raw_adata.obs.columns else 'n_counts_ribo'
        
        if mt_counts_col in raw_adata.obs.columns:
            raw_adata.obs['pct_counts_mt'] = (raw_adata.obs[mt_counts_col] / raw_adata.obs['total_counts']) * 100
        else:
            print("⚠️  Warning: Mitochondrial counts column not found")
            raw_adata.obs['pct_counts_mt'] = 0
        
        if ribo_counts_col in raw_adata.obs.columns:
            raw_adata.obs['pct_counts_ribo'] = (raw_adata.obs[ribo_counts_col] / raw_adata.obs['total_counts']) * 100
        else:
            print("⚠️  Warning: Ribosomal counts column not found")
            raw_adata.obs['pct_counts_ribo'] = 0
        
        print(f"📊 Starting QC: {raw_adata.n_obs:,} cells, {raw_adata.n_vars:,} genes")
        
        # Apply same QC filtering as main pipeline
        n_genes_lower, n_genes_upper = calculate_mad_outliers(raw_adata.obs['n_genes_by_counts'], n_mads=2.5)
        total_counts_lower, total_counts_upper = calculate_mad_outliers(raw_adata.obs['total_counts'], n_mads=2.5)
        mt_lower, mt_upper = calculate_mad_outliers(raw_adata.obs['pct_counts_mt'], n_mads=2.5)
        
        cell_filter = (
            (raw_adata.obs['n_genes_by_counts'] >= max(200, n_genes_lower)) &
            (raw_adata.obs['n_genes_by_counts'] <= n_genes_upper) &
            (raw_adata.obs['total_counts'] >= total_counts_lower) &
            (raw_adata.obs['total_counts'] <= total_counts_upper) &
            (raw_adata.obs['pct_counts_mt'] <= min(25, mt_upper))
        )
        
        # Apply cell filters with detailed logging
        print(f"\n🗑️  Detailed filtering results:")
        print(f"   • Cells before filtering: {len(cell_filter):,}")
        print(f"   • Cells passing QC: {cell_filter.sum():,}")
        print(f"   • Cells removed: {(~cell_filter).sum():,}")
        print(f"   • Retention rate: {cell_filter.sum()/len(cell_filter)*100:.1f}%")
        
        # Show what's being removed by QC criteria
        gene_count_fail = (raw_adata.obs['n_genes_by_counts'] < max(200, n_genes_lower)) | (raw_adata.obs['n_genes_by_counts'] > n_genes_upper)
        umi_count_fail = (raw_adata.obs['total_counts'] < total_counts_lower) | (raw_adata.obs['total_counts'] > total_counts_upper)
        mt_fail = raw_adata.obs['pct_counts_mt'] > min(25, mt_upper)
        
        print(f"\n📊 Cells failing each QC criterion:")
        print(f"   • Gene count outliers: {gene_count_fail.sum():,} cells")
        print(f"   • UMI count outliers: {umi_count_fail.sum():,} cells")
        print(f"   • High MT% (>{min(25, mt_upper):.1f}%): {mt_fail.sum():,} cells")
        
        raw_adata = raw_adata[cell_filter, :].copy()
        
        # Filter genes with logging
        genes_before = raw_adata.n_vars
        sc.pp.filter_genes(raw_adata, min_cells=3)
        genes_after = raw_adata.n_vars
        print(f"\n🧬 Gene filtering:")
        print(f"   • Genes before filtering: {genes_before:,}")
        print(f"   • Genes after filtering: {genes_after:,}")
        print(f"   • Genes removed: {genes_before - genes_after:,}")
        
        print(f"✅ After QC: {raw_adata.n_obs:,} cells, {raw_adata.n_vars:,} genes")
        
        # Save QC checkpoint
        save_checkpoint(raw_adata, CHECKPOINT_QC, "Quality Control")
        goto_step = 'model' if resume_from_step != 'qc' else None
    
    # =======================================
    # 🔁 STEP 3: scVI SETUP & TRAINING (OPTIMIZED)
    # =======================================
    
    # Check if we can skip model training
    if goto_step in ['analysis', 'clustering']:
        print("🔄 RESUMING: Loading scVI model checkpoint...")
        if os.path.exists(CHECKPOINT_SCVI_MODEL):
            try:
                model = scvi.model.SCVI.load(CHECKPOINT_SCVI_MODEL, raw_adata)
                print("✅ scVI model loaded from checkpoint")
                goto_step = resume_from_step
            except Exception as e:
                print(f"⚠️  Failed to load scVI model: {e}")
                print("🔄 Will retrain model...")
                goto_step = 'model'
        else:
            print("📝 No scVI model checkpoint found")
            goto_step = 'model'
    
    if goto_step in [None, 'model']:
        print("\n🔁 STEP 3: scVI MODEL SETUP & TRAINING (M1 MAX OPTIMIZED)")
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
            print(f"   • Batches: {raw_adata.obs[batch_key].nunique()}")
        else:
            scvi.model.SCVI.setup_anndata(raw_adata)
            print("✅ scVI setup complete (no batch correction)")
        
        # Train scVI model with M1 Max optimizations
        print("🏋️  Training scVI model (M1 Max optimized)...")
        
        # Optimized model architecture for your system
        model = scvi.model.SCVI(
            raw_adata, 
            n_latent=64,        # Increased from 50 (more capacity with your RAM)
            n_hidden=512,       # Increased from 256 (better representation)
            n_layers=3,         # Increased from 2 (deeper model)
            dropout_rate=0.1,   # Add dropout for regularization
            gene_likelihood="nb"  # Negative binomial for count data
        )
        
        # Force model to MPS device if available
        if torch.backends.mps.is_available():
            try:
                print(f"🔧 Moving model to MPS device...")
                model.module = model.module.to(device)
                print(f"   • Model device after move: {next(model.module.parameters()).device}")
            except Exception as e:
                print(f"   ⚠️  Could not move to MPS: {e}")
                print(f"   • Model device: {next(model.module.parameters()).device}")
        else:
            print(f"   • Model device: {next(model.module.parameters()).device}")
        
        model.train(
            max_epochs=400,              # More epochs with faster training
            early_stopping=True,
            early_stopping_patience=15,  # More patience
            batch_size=1024,            # Reduced batch size for memory stability
            train_size=0.9,             # Use more data for training
            validation_size=0.1,
            check_val_every_n_epoch=10, # Less frequent validation to save time
            plan_kwargs={
                'lr': 5e-4,             # Slightly higher learning rate
                'weight_decay': 1e-6,   # L2 regularization
                'eps': 1e-8,            # Numerical stability
                'reduce_lr_on_plateau': True,
                'lr_scheduler_metric': 'elbo_validation',
                'lr_patience': 8,
                'lr_factor': 0.6
            },
            # Force MPS usage
            accelerator='mps' if torch.backends.mps.is_available() else 'auto',
            devices=1 if torch.backends.mps.is_available() else 'auto'
        )
        
        print("✅ scVI training complete!")
        
        # Save model checkpoint
        try:
            model.save(CHECKPOINT_SCVI_MODEL, overwrite=True)
            print(f"✅ scVI model saved to: {CHECKPOINT_SCVI_MODEL}")
        except Exception as e:
            print(f"⚠️  Failed to save scVI model: {e}")
        
        goto_step = 'analysis' if resume_from_step not in [None, 'model'] else None
    
    # =======================================
    # 🔍 STEP 4: scVI-SPECIFIC ANALYSIS (OPTIMIZED)
    # =======================================
    
    # Check if we can skip analysis
    if goto_step == 'clustering':
        print("🔄 RESUMING: Loading analysis checkpoint...")
        temp_adata = load_checkpoint(CHECKPOINT_ANALYSIS, "scVI Analysis")
        if temp_adata is not None:
            raw_adata = temp_adata
            print("✅ Resumed from scVI analysis")
        else:
            print("📝 No analysis checkpoint found, running analysis...")
            goto_step = 'analysis'
    
    if goto_step in [None, 'analysis']:
        print("\n🔍 STEP 4: scVI-SPECIFIC ANALYSIS (M1 MAX OPTIMIZED)")
        print("=" * 70)
        
        # Get latent representation with optimized batch size
        print("🔍 Extracting latent representation...")
        raw_adata.obsm["X_scvi"] = model.get_latent_representation(batch_size=1024)  # Reduced from 4096
        
        # Get normalized expression with memory-optimized parameters
        print("🔍 Extracting normalized expression (memory optimized)...")
        try:
            # Use smaller batch size and enable garbage collection
            import gc
            raw_adata.layers["scvi_normalized"] = model.get_normalized_expression(
                library_size=10000,
                batch_size=512,         # Much smaller batch for memory efficiency
                return_mean=True
            )
            gc.collect()  # Force garbage collection
            print("   ✅ Normalized expression extracted successfully")
        except Exception as e:
            print(f"   ⚠️  Error with full normalized expression: {e}")
            print("   🔧 Skipping normalized expression layer to save memory...")
            # Continue without this layer if memory issues persist
        
        # Enable uncertainty calculation with 64GB RAM optimizations
        CALCULATE_UNCERTAINTY = True  # Enable for 64GB RAM system
        
        if CALCULATE_UNCERTAINTY:
            print("🔍 Computing gene expression uncertainty (chunked for 64GB RAM)...")
            try:
                # Optimized chunked parameters for 64GB system
                print("   • Using chunked processing for memory efficiency...")
                
                n_samples = 25          # Full samples for better uncertainty estimates
                batch_size = 512        # Conservative batch size
                chunk_size = 10000      # Process 10K cells at a time
                n_cells = raw_adata.n_obs
                n_genes = raw_adata.n_vars
                
                print(f"   • Processing {n_cells:,} cells in chunks of {chunk_size:,}")
                print(f"   • {n_samples} samples per cell, batch size {batch_size}")
                
                # Initialize arrays to store results
                uncertainty_matrix = np.zeros((n_cells, n_genes), dtype=np.float32)
                cv_matrix = np.zeros((n_cells, n_genes), dtype=np.float32)
                
                # Process data in chunks
                n_chunks = (n_cells + chunk_size - 1) // chunk_size
                
                for chunk_idx in range(n_chunks):
                    start_idx = chunk_idx * chunk_size
                    end_idx = min((chunk_idx + 1) * chunk_size, n_cells)
                    
                    print(f"   • Processing chunk {chunk_idx + 1}/{n_chunks}: cells {start_idx:,}-{end_idx:,}")
                    
                    # Create subset for this chunk
                    chunk_adata = raw_adata[start_idx:end_idx, :].copy()
                    
                    try:
                        # Get uncertainty samples for this chunk
                        chunk_samples = model.get_normalized_expression(
                            adata=chunk_adata,
                            n_samples=n_samples,
                            batch_size=batch_size,
                            return_mean=False
                        )
                        
                        # Calculate uncertainty metrics for this chunk
                        chunk_uncertainty = chunk_samples.var(axis=0)
                        chunk_mean = chunk_samples.mean(axis=0)
                        chunk_std = chunk_samples.std(axis=0)
                        chunk_cv = chunk_std / (chunk_mean + 1e-8)
                        
                        # Store results
                        uncertainty_matrix[start_idx:end_idx, :] = chunk_uncertainty
                        cv_matrix[start_idx:end_idx, :] = chunk_cv
                        
                        # Clean up chunk data
                        del chunk_samples, chunk_uncertainty, chunk_mean, chunk_std, chunk_cv
                        del chunk_adata
                        
                        # Force garbage collection every few chunks
                        if (chunk_idx + 1) % 3 == 0:
                            import gc
                            gc.collect()
                            print(f"     ✓ Completed chunk {chunk_idx + 1}, memory cleaned")
                        
                    except Exception as chunk_error:
                        print(f"     ⚠️  Error in chunk {chunk_idx + 1}: {chunk_error}")
                        print(f"     🔧 Filling with NaN values for this chunk...")
                        uncertainty_matrix[start_idx:end_idx, :] = np.nan
                        cv_matrix[start_idx:end_idx, :] = np.nan
                
                # Store results in adata
                raw_adata.layers["scvi_uncertainty"] = uncertainty_matrix
                raw_adata.layers["scvi_cv"] = cv_matrix
                
                # Calculate summary statistics (excluding NaN values)
                valid_uncertainty = uncertainty_matrix[~np.isnan(uncertainty_matrix)]
                valid_cv = cv_matrix[~np.isnan(cv_matrix)]
                
                if len(valid_uncertainty) > 0:
                    mean_uncertainty = valid_uncertainty.mean()
                    mean_cv = valid_cv.mean()
                    nan_cells = np.isnan(uncertainty_matrix).any(axis=1).sum()
                    
                    print(f"   ✅ Chunked uncertainty calculation completed successfully")
                    print(f"   • Mean gene uncertainty: {mean_uncertainty:.4f}")
                    print(f"   • Mean coefficient of variation: {mean_cv:.4f}")
                    print(f"   • Cells with valid uncertainty: {n_cells - nan_cells:,}/{n_cells:,}")
                    print(f"   • Added layers: 'scvi_uncertainty', 'scvi_cv'")
                else:
                    print(f"   ⚠️  No valid uncertainty values calculated")
                
                # Clean up final arrays
                del uncertainty_matrix, cv_matrix
                import gc
                gc.collect()
                
            except Exception as e:
                print(f"   ⚠️  Chunked uncertainty calculation failed: {e}")
                print("   🔧 Continuing without uncertainty layers...")
        else:
            print("🔍 Skipping uncertainty calculation (CALCULATE_UNCERTAINTY=False)")
            print("   💡 Set CALCULATE_UNCERTAINTY=True to enable if needed")
        
        # Optimized batch entropy calculation (keep this as it's useful)
        if 'ID' in raw_adata.obs.columns:
            print("📊 Computing batch integration metrics (optimized)...")
            
            # Optimized batch entropy calculation using vectorization
            print("   • Computing batch entropy for visualization (vectorized)...")
            
            # Use vectorized operations for better performance
            from sklearn.neighbors import NearestNeighbors
            import gc
            
            try:
                # Fit k-NN model on latent space with smaller neighborhood
                nbrs = NearestNeighbors(
                    n_neighbors=10,     # Reduced from 16 for memory
                    algorithm='auto',
                    n_jobs=4            # Limit parallel jobs
                ).fit(raw_adata.obsm["X_scvi"])
                
                # Get all nearest neighbors at once
                distances, indices = nbrs.kneighbors(raw_adata.obsm["X_scvi"])
                
                # Vectorized entropy calculation
                batch_entropy = []
                batch_ids = raw_adata.obs['ID'].values
                
                # Process in smaller chunks to save memory
                chunk_size = 1000
                for chunk_start in range(0, len(indices), chunk_size):
                    chunk_end = min(chunk_start + chunk_size, len(indices))
                    
                    for i in range(chunk_start, chunk_end):
                        neighbor_batches = batch_ids[indices[i][1:]]  # Exclude self
                        unique, counts = np.unique(neighbor_batches, return_counts=True)
                        probs = counts / counts.sum()
                        entropy = -np.sum(probs * np.log2(probs + 1e-10))
                        batch_entropy.append(entropy)
                    
                    # Force garbage collection every chunk
                    if chunk_start % 5000 == 0:
                        gc.collect()
                
                raw_adata.obs['batch_entropy'] = batch_entropy
                print(f"   • Mean batch entropy: {np.mean(batch_entropy):.2f} (higher = better mixing)")
                
                # Clean up
                del distances, indices, nbrs
                gc.collect()
                
            except Exception as e:
                print(f"   ⚠️  Error computing batch entropy: {e}")
                print("   🔧 Skipping batch entropy calculation...")
        
        # Save analysis checkpoint
        save_checkpoint(raw_adata, CHECKPOINT_ANALYSIS, "scVI Analysis")
        goto_step = 'clustering' if resume_from_step not in [None, 'analysis'] else None
    
    # =======================================
    # 🕸️ STEP 5: NEIGHBORHOOD GRAPH & CLUSTERING (OPTIMIZED)
    # =======================================
    print("\n🕸️ STEP 5: NEIGHBORHOOD GRAPH & CLUSTERING (M1 MAX OPTIMIZED)")
    print("=" * 70)
    
    # Optimized neighborhood graph computation
    print("🕸️  Computing neighborhood graph (optimized)...")
    sc.pp.neighbors(
        raw_adata, 
        use_rep="X_scvi", 
        n_neighbors=20,         # Slightly more neighbors
        n_pcs=50,              # Use more PCs from latent space
        method='umap',         # UMAP method for better performance
        metric='euclidean',
        random_state=42
    )
    
    # Optimized UMAP with more iterations for better embedding
    print("🗺️  Computing UMAP embedding (optimized)...")
    sc.tl.umap(
        raw_adata, 
        random_state=42, 
        min_dist=0.3, 
        spread=1.0,
        n_components=2,
        maxiter=500,           # More iterations for stability
        alpha=1.0,
        gamma=1.0,
        init_pos='spectral'    # Better initialization
    )
    
    # Multi-resolution clustering for robustness
    print("🎯 Multi-resolution clustering analysis...")
    resolutions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0]
    
    for res in resolutions:
        sc.tl.leiden(raw_adata, resolution=res, key_added=f'leiden_res_{res}', random_state=42)
        n_clusters = len(raw_adata.obs[f'leiden_res_{res}'].unique())
        print(f"   • Resolution {res}: {n_clusters} clusters")
    
    # Set primary clustering
    raw_adata.obs['leiden'] = raw_adata.obs['leiden_res_0.3']
    primary_clusters = len(raw_adata.obs['leiden'].unique())
    print(f"✅ Primary clustering: {primary_clusters} clusters (resolution 0.3)")
    
    # Compute cluster stability across resolutions
    print("📊 Analyzing cluster stability...")
    stability_scores = []
    for i in range(len(resolutions)-1):
        res1, res2 = resolutions[i], resolutions[i+1]
        # Adjusted Rand Index between consecutive resolutions
        from sklearn.metrics import adjusted_rand_score
        ari = adjusted_rand_score(
            raw_adata.obs[f'leiden_res_{res1}'], 
            raw_adata.obs[f'leiden_res_{res2}']
        )
        stability_scores.append(ari)
        print(f"   • ARI {res1}-{res2}: {ari:.3f}")
    
    # =======================================
    # 🧬 STEP 6: DIFFERENTIAL EXPRESSION & MARKER GENES
    # =======================================
    print("\n🧬 STEP 6: DIFFERENTIAL EXPRESSION & MARKER GENES")
    print("=" * 70)
    
    # Use scVI-corrected expression for DE analysis
    print("🔍 Finding marker genes using scVI-corrected expression...")
    
    # Ensure we have proper gene names before DE analysis
    if 'scvi_normalized' in raw_adata.layers:
        # Set the normalized layer as the main expression matrix for DE
        raw_adata.X = raw_adata.layers["scvi_normalized"].copy()
        print("   ✅ Using scVI normalized expression for marker gene analysis")
    else:
        print("   ⚠️  scVI normalized layer not found, using raw counts")
    
    # Verify gene names are properly set
    print(f"📋 Gene name verification:")
    print(f"   • var_names type: {type(raw_adata.var_names)}")
    print(f"   • First 5 gene names: {raw_adata.var_names[:5].tolist()}")
    print(f"   • Are gene names numeric? {raw_adata.var_names[0].isdigit() if len(raw_adata.var_names) > 0 else 'Unknown'}")
    
    # Find marker genes for each cluster
    sc.tl.rank_genes_groups(
        raw_adata, 
        groupby='leiden', 
        method='wilcoxon',
        use_raw=False,  # Use scVI-corrected data
        n_genes=50,
        pts=True
    )
    
    # Create marker gene dataframe with proper gene names
    marker_df = sc.get.rank_genes_groups_df(raw_adata, group=None)
    print(f"✅ Identified {len(marker_df)} marker genes across {primary_clusters} clusters")
    
    # Verify marker gene names are correct
    print(f"📋 Marker gene verification:")
    print(f"   • Sample marker genes from cluster 0: {marker_df[marker_df['group'] == '0']['names'].head(5).tolist()}")
    
    # Save top markers per cluster with gene name verification
    top_markers = {}
    for cluster in raw_adata.obs['leiden'].unique():
        cluster_markers = marker_df[marker_df['group'] == cluster].head(10)
        top_markers[cluster] = cluster_markers['names'].tolist()
        # Display actual gene names, not indices
        gene_names = cluster_markers['names'].head(5).tolist()
        print(f"   • Cluster {cluster}: {', '.join(gene_names)}")
    
    # =======================================
    # 📊 STEP 7: COMPREHENSIVE VISUALIZATION
    # =======================================
    print("\n📊 STEP 7: COMPREHENSIVE VISUALIZATION")
    print("=" * 70)
    
    os.makedirs(PLOTS_DIR, exist_ok=True)
    
    # Figure 1: scVI Training and Integration
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Training progress
    if hasattr(model, 'history'):
        train_elbo = model.history.get("elbo_train", [])
        val_elbo = model.history.get("elbo_validation", [])
        if len(train_elbo) > 0:
            axes[0,0].plot(train_elbo, label='Training ELBO', alpha=0.7)
            if len(val_elbo) > 0:
                axes[0,0].plot(val_elbo, label='Validation ELBO', alpha=0.7)
            axes[0,0].set_xlabel('Epoch')
            axes[0,0].set_ylabel('ELBO')
            axes[0,0].legend()
            axes[0,0].set_title('scVI Training Progress')
    
    # Batch integration
    if 'batch_entropy' in raw_adata.obs.columns:
        sc.pl.umap(raw_adata, color='batch_entropy', ax=axes[0,1], show=False, size=2)
        axes[0,1].set_title('Batch Mixing Entropy\n(Higher = Better Integration)')
    
    # Clusters
    sc.pl.umap(raw_adata, color='leiden', ax=axes[1,0], show=False, 
               legend_loc='on data', legend_fontsize=8, size=2)
    axes[1,0].set_title(f'scVI Leiden Clusters ({primary_clusters} clusters)')
    
    # Subject batch
    if 'ID' in raw_adata.obs.columns:
        sc.pl.umap(raw_adata, color='ID', ax=axes[1,1], show=False, size=1)
        axes[1,1].set_title('Subject ID (Batch Effect)')
    
    plt.suptitle('scVI Processing - Training & Integration', fontsize=16, y=0.98)
    plt.tight_layout()
    plt.savefig(f"{PLOTS_DIR}/scvi_training_integration.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Figure 2: Clinical Associations
    if 'Dx_OUD' in raw_adata.obs.columns:
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        sc.pl.umap(raw_adata, color='Dx_OUD', ax=axes[0,0], show=False, size=2)
        axes[0,0].set_title('OUD Diagnosis')
        
        if 'Sex' in raw_adata.obs.columns:
            sc.pl.umap(raw_adata, color='Sex', ax=axes[0,1], show=False, size=2)
            axes[0,1].set_title('Sex')
        
        if 'Age' in raw_adata.obs.columns:
            sc.pl.umap(raw_adata, color='Age', ax=axes[1,0], show=False, size=2)
            axes[1,0].set_title('Age')
        
        sc.pl.umap(raw_adata, color='total_counts', ax=axes[1,1], show=False, size=2)
        axes[1,1].set_title('Total UMI Counts')
        
        plt.suptitle('Clinical Metadata Overview', fontsize=16, y=0.98)
        plt.tight_layout()
        plt.savefig(f"{PLOTS_DIR}/scvi_clinical_overview.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    print(f"✅ Comprehensive plots saved to: {PLOTS_DIR}/")

    # =======================================
    # 💾 STEP 8: SAVE RESULTS
    # =======================================
    print("\n💾 STEP 8: SAVING RESULTS")
    print("=" * 70)
    
    os.makedirs(PROCESSED_DATA_DIR, exist_ok=True)
    
    # Save processed data
    raw_adata.write(OUTPUT_H5AD_PATH_SCVI)
    print(f"💾 scVI results saved to: {OUTPUT_H5AD_PATH_SCVI}")
    
    # Save marker genes
    marker_df_path = f"{PROCESSED_DATA_DIR}/scvi_marker_genes.csv"
    marker_df.to_csv(marker_df_path, index=False)
    print(f"💾 Marker genes saved to: {marker_df_path}")
    
    # Determine batch correction used
    batch_key = 'ID' if 'ID' in raw_adata.obs.columns else 'None'
    
    # Save summary
    with open(SUMMARY_PATH_SCVI, 'w') as f:
        f.write("GSE225158 OUD Striatum snRNA-seq scVI Processing Summary\n")
        f.write("=" * 60 + "\n")
        f.write(f"Processing date: {pd.Timestamp.now()}\n")
        f.write(f"Method: scVI probabilistic modeling\n\n")
        f.write("DATASET OVERVIEW:\n")
        f.write(f"• Final cells: {raw_adata.n_obs:,}\n")
        f.write(f"• Final genes: {raw_adata.n_vars:,}\n")
        f.write(f"• Clusters: {primary_clusters}\n")
        f.write(f"• Batch correction: {batch_key}\n\n")
        f.write("OUTPUTS:\n")
        f.write(f"• Processed data: {OUTPUT_H5AD_SCVI}\n")
        f.write(f"• Marker genes: scvi_marker_genes.csv\n")
        f.write(f"• Plots: {PLOTS_DIR}/\n\n")
        f.write("NEXT STEPS:\n")
        f.write("• Run 02_cell_type_annotation.py for cell type assignment\n")
    
    print(f"💾 Summary saved to: {SUMMARY_PATH_SCVI}")
    
    print("\n" + "=" * 70)
    print("✅ scVI PREPROCESSING COMPLETE!")
    print("=" * 70)
    print(f"📊 Dataset: {raw_adata.n_obs:,} cells × {raw_adata.n_vars:,} genes")
    print(f"🎯 Clusters: {primary_clusters}")
    print(f"🔬 Method: scVI with batch correction")
    print(f"📁 Ready for cell type annotation!")
    print("=" * 70)
    
    return raw_adata

if __name__ == "__main__":
    import sys
    
    # Check for resume argument
    resume_step = None
    if len(sys.argv) > 1:
        resume_step = sys.argv[1]
        valid_steps = ['qc', 'model', 'analysis', 'clustering']
        if resume_step not in valid_steps:
            print(f"❌ Invalid resume step. Valid options: {valid_steps}")
            sys.exit(1)
        print(f"🔄 Resuming from: {resume_step}")
    
    scvi_adata = main_scvi(resume_from_step=resume_step)
