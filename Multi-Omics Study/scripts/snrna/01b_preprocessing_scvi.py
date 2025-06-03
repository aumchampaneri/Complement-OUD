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

# =======================================
# 📚 IMPORT LIBRARIES
# =======================================

import scanpy as sc
import scvi
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import sparse
import warnings
import os
import torch
warnings.filterwarnings('ignore')

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

def main_scvi():
    """scVI-based processing pipeline optimized for M1 Max"""
    
    print("=" * 70)
    print("🔬 scVI-BASED PREPROCESSING PIPELINE (M1 MAX OPTIMIZED)")
    print("=" * 70)
    
    # =======================================
    # 🔬 STEP 1: LOAD DATA & EXTRACT RAW COUNTS
    # =======================================
    print("🔬 STEP 1: DATA LOADING & RAW EXTRACTION")
    print("=" * 70)
    
    # Load and extract raw data (same as main pipeline)
    adata = sc.read_h5ad(INPUT_PATH)
    print(f"📊 Loaded dataset: {adata.shape}")
    
    # Extract raw counts
    if adata.raw is not None:
        raw_adata = adata.raw.to_adata()
        
        # Copy clinical metadata
        clinical_cols = ['ID', 'Region', 'Case', 'Sex', 'Race', 'Age', 'BMI', 'PMI', 'pH', 'RIN', 
                        'Dx_OUD', 'Dx_Substances', 'Dx_Comorbid', 'celltype1', 'celltype2', 'celltype3']
        available_clinical = [col for col in clinical_cols if col in adata.obs.columns]
        for col in available_clinical:
            raw_adata.obs[col] = adata.obs[col].copy()
        
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
    
    # Add percentage metrics with robust column name handling
    # Handle potential column name variations
    mt_counts_col = 'n_counts_mt' if 'n_counts_mt' in raw_adata.obs.columns else 'mt_counts'
    ribo_counts_col = 'n_counts_ribo' if 'n_counts_ribo' in raw_adata.obs.columns else 'ribo_counts'
    
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
    
    raw_adata = raw_adata[cell_filter, :].copy()
    sc.pp.filter_genes(raw_adata, min_cells=3)
    
    print(f"✅ After QC: {raw_adata.n_obs:,} cells, {raw_adata.n_vars:,} genes")
    
    # =======================================
    # 🔁 STEP 3: scVI SETUP & TRAINING (OPTIMIZED)
    # =======================================
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
    
    # Optimized training parameters for M1 Max
    print(f"   • Model device: {next(model.module.parameters()).device}")
    
    model.train(
        max_epochs=400,              # More epochs with faster training
        early_stopping=True,
        early_stopping_patience=15,  # More patience
        batch_size=2048,            # Larger batch size for your RAM
        train_size=0.9,             # Use more data for training
        validation_size=0.1,
        check_val_every_n_epoch=5,  # More frequent validation
        plan_kwargs={
            'lr': 5e-4,             # Slightly higher learning rate
            'weight_decay': 1e-6,   # L2 regularization
            'eps': 1e-8,            # Numerical stability
            'reduce_lr_on_plateau': True,
            'lr_scheduler_metric': 'elbo_validation',
            'lr_patience': 8,
            'lr_factor': 0.6
        },
        # Remove use_gpu parameter - scvi-tools auto-detects MPS
        accelerator='auto',         # Let scvi-tools auto-detect MPS/GPU/CPU
        devices='auto'              # Auto-select available devices
    )
    
    # Create plots directory before saving
    os.makedirs(PLOTS_DIR, exist_ok=True)
    
    # Plot training history
    train_elbo = model.history["elbo_train"]
    val_elbo = model.history["elbo_validation"] 
    
    plt.figure(figsize=(10, 6))
    plt.plot(train_elbo, label='Training ELBO')
    plt.plot(val_elbo, label='Validation ELBO')
    plt.xlabel('Epoch')
    plt.ylabel('ELBO')
    plt.legend()
    plt.title('scVI Training Progress')
    plt.savefig(f"{PLOTS_DIR}/scvi_training_progress.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print("✅ scVI training complete!")
    
    # =======================================
    # 🔍 STEP 4: scVI-SPECIFIC ANALYSIS (OPTIMIZED)
    # =======================================
    print("\n🔍 STEP 4: scVI-SPECIFIC ANALYSIS (M1 MAX OPTIMIZED)")
    print("=" * 70)
    
    # Get latent representation with larger batch size
    print("🔍 Extracting latent representation...")
    raw_adata.obsm["X_scvi"] = model.get_latent_representation(batch_size=4096)
    
    # Get normalized expression with optimized parameters
    print("🔍 Extracting normalized expression...")
    raw_adata.layers["scvi_normalized"] = model.get_normalized_expression(
        library_size=10000,
        batch_size=4096,        # Larger batch for your RAM
        return_mean=True
    )
    
    # Optimized uncertainty calculation
    print("🔍 Computing gene expression uncertainty (optimized)...")
    # Use fewer samples but larger batches for efficiency
    uncertainty_samples = model.get_normalized_expression(
        n_samples=10,           # Reduced from 25 for speed
        batch_size=4096,        # Larger batch size
        return_mean=False
    )
    raw_adata.layers["scvi_uncertainty"] = uncertainty_samples.var(axis=0)
    
    # Optimized batch entropy calculation
    if 'ID' in raw_adata.obs.columns:
        print("📊 Computing batch integration metrics (optimized)...")
        
        # ...existing scib_metrics code...
        
        # Optimized batch entropy calculation using vectorization
        print("   • Computing batch entropy for visualization (vectorized)...")
        
        # Use vectorized operations for better performance
        from sklearn.neighbors import NearestNeighbors
        
        # Fit k-NN model on latent space
        nbrs = NearestNeighbors(
            n_neighbors=16, 
            algorithm='auto',
            n_jobs=-1  # Use all cores
        ).fit(raw_adata.obsm["X_scvi"])
        
        # Get all nearest neighbors at once
        distances, indices = nbrs.kneighbors(raw_adata.obsm["X_scvi"])
        
        # Vectorized entropy calculation
        batch_entropy = []
        batch_ids = raw_adata.obs['ID'].values
        
        for i in range(len(indices)):
            neighbor_batches = batch_ids[indices[i][1:]]  # Exclude self
            unique, counts = np.unique(neighbor_batches, return_counts=True)
            probs = counts / counts.sum()
            entropy = -np.sum(probs * np.log2(probs + 1e-10))
            batch_entropy.append(entropy)
        
        raw_adata.obs['batch_entropy'] = batch_entropy
        print(f"   • Mean batch entropy: {np.mean(batch_entropy):.2f} (higher = better mixing)")

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
    
    # Set the normalized layer as the main expression matrix for DE
    raw_adata.X = raw_adata.layers["scvi_normalized"].copy()
    
    # Find marker genes for each cluster
    sc.tl.rank_genes_groups(
        raw_adata, 
        groupby='leiden', 
        method='wilcoxon',
        use_raw=False,  # Use scVI-corrected data
        n_genes=50,
        pts=True
    )
    
    # Create marker gene dataframe
    marker_df = sc.get.rank_genes_groups_df(raw_adata, group=None)
    print(f"✅ Identified {len(marker_df)} marker genes across {primary_clusters} clusters")
    
    # Save top markers per cluster
    top_markers = {}
    for cluster in raw_adata.obs['leiden'].unique():
        cluster_markers = marker_df[marker_df['group'] == cluster].head(10)
        top_markers[cluster] = cluster_markers['names'].tolist()
        print(f"   • Cluster {cluster}: {', '.join(cluster_markers['names'].head(5).tolist())}")
    
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
    scvi_adata = main_scvi()
