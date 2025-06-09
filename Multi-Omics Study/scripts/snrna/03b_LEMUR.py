#!/usr/bin/env python3
"""
🌊 LEMUR Analysis - Latent Expression Modeling for Unified Representation
GSE225158 - OUD vs Control using pyLEMUR for continuous differential expression

This script performs continuous differential expression analysis using LEMUR,
which models cell states in latent space and captures gradual transitions
between conditions while preserving neighborhood relationships.

Author: Research Team
Date: 2024
"""

import os
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for server environments
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests

# Configure scanpy settings
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=300, facecolor='white')
warnings.filterwarnings('ignore')

# LEMUR dependencies
try:
    import pylemur
    import anndata as ad
    PYLEMUR_AVAILABLE = True
    print("✅ pyLEMUR available")
except ImportError:
    PYLEMUR_AVAILABLE = False
    print("❌ pyLEMUR not available. Install with: pip install pylemur")

# ============================================================================
# 📁 CONFIGURATION
# ============================================================================

# File paths
BASE_DIR = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study"
INPUT_H5AD = f"{BASE_DIR}/data/processed/snrna_scvi/GSE225158_annotated_scvi.h5ad"
OUTPUT_DIR = f"{BASE_DIR}/results/snrna_scvi/lemur_analysis"
PLOTS_DIR = f"{OUTPUT_DIR}/plots"
TABLES_DIR = f"{OUTPUT_DIR}/tables"

# Analysis parameters
ANALYSIS_CONFIG = {
    'n_embedding': 15,           # Number of LEMUR embedding dimensions
    'max_cells': 50000,          # Maximum cells for memory management
    'max_cells_fallback': 40000, # Maximum cells for fallback analysis
    'n_hvg': 3000,              # Number of highly variable genes
    'random_seed': 42,          # Random seed for reproducibility
    
    # Significance thresholds
    'fdr_lenient': 0.1,         # Lenient FDR threshold
    'lfc_lenient': 0.1,         # Lenient log fold change threshold
    'fdr_strict': 0.05,         # Strict FDR threshold
    'lfc_strict': 0.25,         # Strict log fold change threshold
    
    # Plotting parameters
    'figure_dpi': 300,          # Plot resolution
    'point_size': 5,            # Scatter plot point size
    'alpha': 0.6,               # Point transparency
}

# Create output directories
for directory in [OUTPUT_DIR, PLOTS_DIR, TABLES_DIR]:
    os.makedirs(directory, exist_ok=True)

print("🌊 LEMUR CONTINUOUS STATE ANALYSIS")
print("==================================")
print(f"Input: {INPUT_H5AD}")
print(f"Output: {OUTPUT_DIR}")

# ============================================================================
# 📊 DATA LOADING AND PREPARATION
# ============================================================================

def load_and_prepare_data():
    """Load and prepare data for LEMUR analysis"""
    print(f"\n📁 LOADING ANNOTATED DATA")
    print("=" * 30)
    
    if not os.path.exists(INPUT_H5AD):
        print(f"   ❌ File not found: {INPUT_H5AD}")
        print("   Make sure to run 02_scvi_annotation.py first")
        return None
    
    # Load the data
    adata = sc.read_h5ad(INPUT_H5AD)
    print(f"   Loaded data: {adata.n_obs} cells × {adata.n_vars} genes")
    
    # Display available columns for inspection
    available_cols = list(adata.obs.columns)
    if len(available_cols) > 10:
        print(f"   Available columns: {available_cols[:10]}... ({len(available_cols)} total)")
    else:
        print(f"   Available columns: {available_cols}")
    
    # Create condition column from Dx_OUD if needed
    if 'condition' not in adata.obs.columns:
        if 'Dx_OUD' in adata.obs.columns:
            print("   Creating condition column from Dx_OUD...")
            adata.obs['condition'] = adata.obs['Dx_OUD'].map({
                'OUD': 'OUD',
                'None': 'Control'
            }).fillna('Control')
        else:
            print("   ❌ No condition information found")
            return None
    
    print(f"   Conditions: {dict(adata.obs['condition'].value_counts())}")
    
    # Create sample_id if needed
    if 'sample_id' not in adata.obs.columns and 'ID' in adata.obs.columns:
        print("   Creating sample_id from ID column...")
        adata.obs['sample_id'] = adata.obs['ID'].astype(str)
    
    # Display metadata information
    metadata_info = []
    for col in ['sex', 'region', 'cell_type']:
        if col in adata.obs.columns:
            unique_vals = adata.obs[col].unique()
            metadata_info.append(f"{col.title()}: {len(unique_vals)} categories")
    
    if metadata_info:
        print("   " + " | ".join(metadata_info))
    
    print(f"   Cell types: {adata.obs['cell_type'].nunique()}")
    print(f"   Conditions: {adata.obs['condition'].unique().tolist()}")
    if 'sample_id' in adata.obs.columns:
        print(f"   Samples: {adata.obs['sample_id'].nunique()}")
    
    return adata

# ============================================================================
# 🌊 LEMUR ANALYSIS FUNCTIONS
# ============================================================================

def run_pylemur_analysis(adata, n_embedding=None, max_cells=None):
    """Run pyLEMUR analysis for continuous state modeling"""
    print("\n🌊 RUNNING pyLEMUR ANALYSIS")
    print("=" * 40)
    
    # Use config defaults if not provided
    if n_embedding is None:
        n_embedding = ANALYSIS_CONFIG['n_embedding']
    if max_cells is None:
        max_cells = ANALYSIS_CONFIG['max_cells']
    
    if not PYLEMUR_AVAILABLE:
        print("   pyLEMUR not available - running enhanced continuous analysis")
        return run_enhanced_continuous_analysis(adata, n_embedding)
    
    try:
        # Memory optimization for large datasets
        if adata.n_obs > max_cells:
            print(f"   Large dataset detected. Subsampling to {max_cells} cells for memory efficiency...")
            adata_sample = subsample_balanced(adata, max_cells)
        else:
            adata_sample = adata.copy()
        
        print("   Selecting highly variable genes...")
        sc.pp.highly_variable_genes(adata_sample, n_top_genes=ANALYSIS_CONFIG['n_hvg'], flavor='seurat_v3')
        adata_hvg = adata_sample[:, adata_sample.var.highly_variable].copy()
        
        # Prepare expression matrix
        if adata_hvg.raw is not None:
            X = adata_hvg.raw.X
            gene_names = adata_hvg.raw.var_names
        else:
            X = adata_hvg.X
            gene_names = adata_hvg.var_names
        
        # Convert to dense and normalize
        if sparse.issparse(X):
            X = X.toarray()
        
        print("   Normalizing expression data...")
        X_norm = normalize_expression(X)
        
        print(f"   Expression matrix: {X_norm.shape}")
        print(f"   Conditions: OUD={np.sum(adata_hvg.obs['condition'] == 'OUD')}, "
              f"Control={np.sum(adata_hvg.obs['condition'] == 'Control')}")
        
        # Create LEMUR model
        print("   Creating LEMUR model...")
        adata_lemur = ad.AnnData(X=X_norm)
        adata_lemur.obs = adata_hvg.obs.copy()
        adata_lemur.var_names = gene_names
        
        lemur_model = pylemur.tl.LEMUR(
            adata_lemur, 
            design="~ condition",
            n_embedding=n_embedding
        )
        
        # Fit the model
        print("   Fitting LEMUR model...")
        lemur_model.fit(verbose=False)
        
        # Perform differential expression testing
        print("   Running differential expression test...")
        de_results = perform_de_testing(lemur_model, X_norm, adata_hvg.obs, gene_names)
        
        print(f"   ✅ Analysis completed!")
        print(f"   Significant genes (FDR < 0.1, |LFC| > 0.1): {de_results['n_sig_lenient']}")
        print(f"   Significant genes (FDR < 0.05, |LFC| > 0.25): {de_results['n_sig_strict']}")
        
        return {
            'embedding': lemur_model.embedding,
            'differential_results': de_results['de_df'],
            'lemur_result': lemur_model,
            'adata_used': adata_sample,
            'method': 'pylemur'
        }
        
    except Exception as e:
        print(f"   ❌ pyLEMUR analysis failed: {str(e)}")
        print("   Falling back to enhanced continuous analysis...")
        return run_enhanced_continuous_analysis(adata, n_embedding)

def subsample_balanced(adata, n_sample):
    """Perform balanced subsampling maintaining condition proportions"""
    # Check if we have valid conditions
    if 'condition' not in adata.obs.columns:
        sampled_cells = np.random.choice(adata.obs.index, min(n_sample, len(adata.obs)), replace=False)
        return adata[sampled_cells].copy()
    
    # Remove any NaN conditions
    valid_mask = ~adata.obs['condition'].isna()
    if not valid_mask.all():
        adata = adata[valid_mask].copy()
    
    control_cells = adata.obs[adata.obs['condition'] == 'Control'].index
    oud_cells = adata.obs[adata.obs['condition'] == 'OUD'].index
    
    if len(control_cells) == 0 or len(oud_cells) == 0:
        sampled_cells = np.random.choice(adata.obs.index, min(n_sample, len(adata.obs)), replace=False)
        return adata[sampled_cells].copy()
    
    control_prop = len(control_cells) / len(adata.obs)
    n_control = int(n_sample * control_prop)
    n_oud = n_sample - n_control
    
    sampled_control = np.random.choice(control_cells, min(n_control, len(control_cells)), replace=False)
    sampled_oud = np.random.choice(oud_cells, min(n_oud, len(oud_cells)), replace=False)
    
    sampled_cells = np.concatenate([sampled_control, sampled_oud])
    return adata[sampled_cells].copy()

def normalize_expression(X):
    """Normalize expression matrix"""
    X_norm = X.copy()
    # Library size normalization
    lib_sizes = np.sum(X_norm, axis=1)
    X_norm = X_norm / lib_sizes[:, np.newaxis] * np.median(lib_sizes)
    # Log transformation
    return np.log1p(X_norm)

def perform_de_testing(lemur_model, X_norm, obs_data, gene_names):
    """Perform differential expression testing using LEMUR linear coefficients"""
    # Extract linear coefficients
    linear_coeffs = lemur_model.linear_coefficients
    
    # Get condition effect (second row for OUD vs Control)
    if linear_coeffs.shape[0] > 1:
        oud_effects = linear_coeffs[1, :]
    else:
        oud_effects = linear_coeffs[0, :]
    
    # Calculate condition masks
    oud_mask = obs_data['condition'] == 'OUD'
    ctrl_mask = obs_data['condition'] == 'Control'
    
    # Calculate mean expression and variance for statistical testing
    mean_oud = np.mean(X_norm[oud_mask, :], axis=0)
    mean_ctrl = np.mean(X_norm[ctrl_mask, :], axis=0)
    
    # Calculate pooled variance for standard error estimation
    var_oud = np.var(X_norm[oud_mask, :], axis=0)
    var_ctrl = np.var(X_norm[ctrl_mask, :], axis=0)
    pooled_var = (var_oud + var_ctrl) / 2
    
    # Standard error calculation
    n_oud = np.sum(oud_mask)
    n_ctrl = np.sum(ctrl_mask)
    stderr = np.sqrt(pooled_var / n_oud + pooled_var / n_ctrl)
    
    # T-statistics and p-values
    t_stats = oud_effects / (stderr + 1e-8)
    df = min(n_oud, n_ctrl) - 1
    p_vals = 2 * (1 - stats.t.cdf(np.abs(t_stats), df))
    
    # FDR correction
    valid_pvals = ~np.isnan(p_vals) & (p_vals >= 0) & (p_vals <= 1)
    adj_pvals = np.full(len(p_vals), 1.0)
    
    if valid_pvals.sum() > 0:
        _, adj_pvals[valid_pvals], _, _ = multipletests(
            p_vals[valid_pvals], method='fdr_bh'
        )
    
    # Calculate effect size (Cohen's d)
    pooled_std = np.sqrt(pooled_var)
    effect_size = (mean_oud - mean_ctrl) / (pooled_std + 1e-8)
    
    # Create results DataFrame
    de_df = pd.DataFrame({
        'gene': gene_names,
        'lfc': oud_effects,
        'pval': p_vals,
        'adj_pval': adj_pvals,
        'stat': t_stats,
        'mean_oud': mean_oud,
        'mean_ctrl': mean_ctrl,
        'effect_size': effect_size
    })
    
    # Add significance flags
    de_df['significant'] = (de_df['adj_pval'] < ANALYSIS_CONFIG['fdr_lenient']) & (abs(de_df['lfc']) > ANALYSIS_CONFIG['lfc_lenient'])
    de_df['significant_strict'] = (de_df['adj_pval'] < ANALYSIS_CONFIG['fdr_strict']) & (abs(de_df['lfc']) > ANALYSIS_CONFIG['lfc_strict'])
    
    return {
        'de_df': de_df,
        'n_sig_lenient': de_df['significant'].sum(),
        'n_sig_strict': de_df['significant_strict'].sum()
    }

def run_enhanced_continuous_analysis(adata, n_embedding=None):
    """Enhanced continuous analysis using neighborhood-preserving methods"""
    print("\n🔧 ENHANCED CONTINUOUS ANALYSIS")
    print("=" * 40)
    
    # Use config defaults if not provided
    if n_embedding is None:
        n_embedding = ANALYSIS_CONFIG['n_embedding']
    
    # Memory optimization for large datasets
    if adata.n_obs > ANALYSIS_CONFIG['max_cells_fallback']:
        print(f"   Large dataset detected. Subsampling for memory efficiency...")
        adata_work = subsample_balanced(adata, ANALYSIS_CONFIG['max_cells_fallback'])
    else:
        adata_work = adata.copy()
    
    # Use highly variable genes
    print("   Identifying highly variable genes...")
    sc.pp.highly_variable_genes(adata_work, n_top_genes=ANALYSIS_CONFIG['n_hvg'], flavor='seurat_v3')
    adata_hvg = adata_work[:, adata_work.var.highly_variable].copy()
    
    # Preprocessing
    print("   Preprocessing expression data...")
    sc.pp.normalize_total(adata_hvg, target_sum=1e4)
    sc.pp.log1p(adata_hvg)
    sc.pp.scale(adata_hvg, max_value=10)
    
    # Create neighborhood-preserving embedding
    print("   Creating neighborhood-preserving embedding...")
    sc.tl.pca(adata_hvg, svd_solver='arpack', n_comps=50)
    sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=40)
    sc.tl.umap(adata_hvg, min_dist=0.5, spread=1.0)
    
    # Create multi-dimensional embedding
    pca_coords = adata_hvg.obsm['X_pca'][:, :min(n_embedding-2, 48)]
    umap_coords = adata_hvg.obsm['X_umap']
    
    if pca_coords.shape[1] + 2 == n_embedding:
        embedding_coords = np.column_stack([umap_coords, pca_coords])
    else:
        needed_dims = n_embedding - 2
        if pca_coords.shape[1] >= needed_dims:
            embedding_coords = np.column_stack([umap_coords, pca_coords[:, :needed_dims]])
        else:
            extra_dims = needed_dims - pca_coords.shape[1]
            padding = np.zeros((adata_hvg.n_obs, extra_dims))
            embedding_coords = np.column_stack([umap_coords, pca_coords, padding])
    
    # Perform neighborhood-aware differential expression
    print("   Running neighborhood-aware differential analysis...")
    de_results = perform_neighborhood_de(adata, adata_hvg)
    
    print(f"   ✅ Analysis completed!")
    print(f"   Significant genes: {de_results['n_sig_lenient']}")
    
    return {
        'embedding': embedding_coords,
        'differential_results': de_results['de_df'],
        'adata_processed': adata_hvg,
        'adata_used': adata_work,
        'method': 'enhanced_continuous'
    }

def perform_neighborhood_de(adata, adata_hvg):
    """Perform neighborhood-aware differential expression analysis"""
    # Get original expression matrix
    if adata.raw is not None:
        expr_matrix = adata.raw.X
        gene_names = adata.raw.var_names
    else:
        expr_matrix = adata.X
        gene_names = adata.var_names
    
    if sparse.issparse(expr_matrix):
        expr_matrix = expr_matrix.toarray()
    
    condition_numeric = (adata.obs['condition'] == 'OUD').astype(int)
    neighbors_graph = adata_hvg.obsp['connectivities']
    
    de_results = []
    
    print("   Computing neighborhood statistics...")
    for i, gene in enumerate(gene_names):
        if i % 10000 == 0 and i > 0:
            print(f"     Processing gene {i+1}/{len(gene_names)}")
        
        gene_expr = expr_matrix[:, i]
        
        # Smooth expression using neighborhood graph
        if sparse.issparse(neighbors_graph):
            smoothed_expr = neighbors_graph.dot(gene_expr) / neighbors_graph.sum(axis=1).A1
        else:
            smoothed_expr = neighbors_graph.dot(gene_expr) / neighbors_graph.sum(axis=1)
        
        # Test differential expression
        oud_expr = smoothed_expr[condition_numeric == 1]
        ctrl_expr = smoothed_expr[condition_numeric == 0]
        
        if len(oud_expr) > 0 and len(ctrl_expr) > 0 and np.var(smoothed_expr) > 0:
            t_stat, p_val = stats.ttest_ind(oud_expr, ctrl_expr)
            mean_oud = np.mean(oud_expr)
            mean_ctrl = np.mean(ctrl_expr)
            lfc = np.log2((mean_oud + 1e-8) / (mean_ctrl + 1e-8))
            effect_size = (mean_oud - mean_ctrl) / np.sqrt(np.var(smoothed_expr))
        else:
            t_stat, p_val, lfc, effect_size = 0, 1, 0, 0
            mean_oud = mean_ctrl = 0
        
        de_results.append({
            'gene': gene,
            'lfc': lfc,
            'stat': t_stat,
            'pval': p_val,
            'effect_size': effect_size,
            'mean_oud': mean_oud,
            'mean_ctrl': mean_ctrl
        })
    
    # Convert to DataFrame and adjust p-values
    de_df = pd.DataFrame(de_results)
    
    valid_pvals = ~np.isnan(de_df['pval'])
    adj_pvals = np.full(len(de_df), 1.0)
    
    if valid_pvals.sum() > 0:
        _, adj_pvals[valid_pvals], _, _ = multipletests(
            de_df.loc[valid_pvals, 'pval'], method='fdr_bh'
        )
    
    de_df['adj_pval'] = adj_pvals
    de_df['significant'] = (de_df['adj_pval'] < ANALYSIS_CONFIG['fdr_lenient']) & (abs(de_df['lfc']) > ANALYSIS_CONFIG['lfc_lenient'])
    de_df['significant_strict'] = (de_df['adj_pval'] < ANALYSIS_CONFIG['fdr_strict']) & (abs(de_df['lfc']) > ANALYSIS_CONFIG['lfc_strict'])
    
    return {
        'de_df': de_df,
        'n_sig_lenient': de_df['significant'].sum(),
        'n_sig_strict': de_df['significant_strict'].sum()
    }

# ============================================================================
# 📊 VISUALIZATION FUNCTIONS
# ============================================================================

def plot_lemur_results(results, adata, save_plots=True):
    """Create comprehensive visualizations for LEMUR results"""
    print(f"\n📊 CREATING LEMUR VISUALIZATIONS")
    print("=" * 40)
    
    embedding = results['embedding']
    de_df = results['differential_results']
    
    # Create figure with multiple subplots
    fig = plt.figure(figsize=(20, 15))
    
    # 1. Embedding colored by condition
    ax1 = plt.subplot(3, 4, 1)
    condition_colors = {'OUD': 'red', 'Control': 'blue'}
    for condition in adata.obs['condition'].unique():
        mask = adata.obs['condition'] == condition
        plt.scatter(embedding[mask, 0], embedding[mask, 1], 
                   c=condition_colors.get(condition, 'gray'),
                   label=condition, alpha=0.6, s=5)
    plt.xlabel('LEMUR Embedding 1')
    plt.ylabel('LEMUR Embedding 2')
    plt.title('LEMUR Embedding - Condition')
    plt.legend()
    
    # 2. Embedding colored by cell type
    ax2 = plt.subplot(3, 4, 2)
    cell_types = adata.obs['cell_type'].unique()
    colors = plt.cm.tab20(np.linspace(0, 1, len(cell_types)))
    for i, cell_type in enumerate(cell_types):
        mask = adata.obs['cell_type'] == cell_type
        plt.scatter(embedding[mask, 0], embedding[mask, 1], 
                   c=[colors[i]], label=cell_type, alpha=0.6, s=5)
    plt.xlabel('LEMUR Embedding 1')
    plt.ylabel('LEMUR Embedding 2')
    plt.title('LEMUR Embedding - Cell Type')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # 3. Volcano plot
    ax3 = plt.subplot(3, 4, 3)
    x = de_df['lfc']
    y = -np.log10(de_df['adj_pval'].clip(lower=1e-300))
    
    sig_mask = (de_df['adj_pval'] < ANALYSIS_CONFIG['fdr_lenient']) & (abs(de_df['lfc']) > ANALYSIS_CONFIG['lfc_lenient'])
    colors = np.where(sig_mask, 
                     np.where(de_df['lfc'] > 0, 'red', 'blue'), 'lightgray')
    
    plt.scatter(x, y, c=colors, alpha=ANALYSIS_CONFIG['alpha'], s=10)
    plt.axhline(-np.log10(ANALYSIS_CONFIG['fdr_lenient']), color='red', linestyle='--', alpha=0.7, label=f"FDR={ANALYSIS_CONFIG['fdr_lenient']}")
    plt.axhline(-np.log10(ANALYSIS_CONFIG['fdr_strict']), color='black', linestyle='--', alpha=0.5, label=f"FDR={ANALYSIS_CONFIG['fdr_strict']}")
    plt.axvline(ANALYSIS_CONFIG['lfc_lenient'], color='black', linestyle='--', alpha=0.5)
    plt.axvline(-ANALYSIS_CONFIG['lfc_lenient'], color='black', linestyle='--', alpha=0.5)
    plt.xlabel('Log2 Fold Change (OUD vs Control)')
    plt.ylabel('-log10(Adjusted p-value)')
    plt.title('Differential Expression\n(Volcano Plot)')
    plt.legend(fontsize=8)
    
    # 4. MA plot
    ax4 = plt.subplot(3, 4, 4)
    if 'mean_oud' in de_df.columns and 'mean_ctrl' in de_df.columns:
        mean_expr = (de_df['mean_oud'] + de_df['mean_ctrl']) / 2
        plt.scatter(np.log10(mean_expr + 1e-8), de_df['lfc'], c=colors, alpha=ANALYSIS_CONFIG['alpha'], s=10)
        plt.xlabel('log10(Mean Expression)')
    else:
        plt.scatter(range(len(de_df)), de_df['lfc'], c=colors, alpha=ANALYSIS_CONFIG['alpha'], s=10)
        plt.xlabel('Gene Index')
    plt.axhline(0, color='black', linestyle='-', alpha=0.5)
    plt.axhline(ANALYSIS_CONFIG['lfc_lenient'], color='red', linestyle='--', alpha=0.5)
    plt.axhline(-ANALYSIS_CONFIG['lfc_lenient'], color='red', linestyle='--', alpha=0.5)
    plt.ylabel('Log2 Fold Change')
    plt.title('MA Plot')
    
    # 5. Top upregulated genes
    ax5 = plt.subplot(3, 4, 5)
    top_up = de_df[de_df['significant'] & (de_df['lfc'] > 0)].nlargest(15, 'lfc')
    if len(top_up) > 0:
        y_pos = np.arange(len(top_up))
        plt.barh(y_pos, top_up['lfc'], color='red', alpha=0.7)
        plt.yticks(y_pos, top_up['gene'])
        plt.xlabel('Log2 Fold Change')
        plt.title('Top Upregulated Genes')
        plt.gca().invert_yaxis()
    
    # 6. Top downregulated genes
    ax6 = plt.subplot(3, 4, 6)
    top_down = de_df[de_df['significant'] & (de_df['lfc'] < 0)].nsmallest(15, 'lfc')
    if len(top_down) > 0:
        y_pos = np.arange(len(top_down))
        plt.barh(y_pos, top_down['lfc'], color='blue', alpha=0.7)
        plt.yticks(y_pos, top_down['gene'])
        plt.xlabel('Log2 Fold Change')
        plt.title('Top Downregulated Genes')
        plt.gca().invert_yaxis()
    
    # 7. P-value distribution
    ax7 = plt.subplot(3, 4, 7)
    plt.hist(de_df['pval'].dropna(), bins=50, alpha=0.7, color='skyblue', edgecolor='black')
    plt.xlabel('P-value')
    plt.ylabel('Frequency')
    plt.title('P-value Distribution')
    
    # 8. Effect size vs significance
    ax8 = plt.subplot(3, 4, 8)
    effect_col = 'effect_size' if 'effect_size' in de_df.columns else 'lfc'
    effect_data = de_df[effect_col].abs() if effect_col == 'lfc' else de_df[effect_col]
    plt.scatter(effect_data, -np.log10(de_df['adj_pval'].clip(lower=1e-300)), 
               c=colors, alpha=ANALYSIS_CONFIG['alpha'], s=10)
    plt.xlabel('Effect Size' if effect_col == 'effect_size' else '|Log2 Fold Change|')
    plt.ylabel('-log10(Adjusted p-value)')
    plt.title('Effect Size vs Significance')
    
    # 9-12. Additional metadata plots
    plot_metadata_distributions(embedding, adata, fig)
    
    # Summary statistics
    ax12 = plt.subplot(3, 4, 12)
    ax12.axis('off')
    
    n_total = len(de_df)
    n_sig = de_df['significant'].sum()
    n_up = ((de_df['significant']) & (de_df['lfc'] > 0)).sum()
    n_down = ((de_df['significant']) & (de_df['lfc'] < 0)).sum()
    
    summary_text = f"""
    LEMUR Analysis Summary
    =====================
    
    Method: {results['method']}
    
    Total genes: {n_total:,}
    Significant: {n_sig:,} ({n_sig/n_total*100:.1f}%)
    
    Upregulated: {n_up:,}
    Downregulated: {n_down:,}
    
    Embedding dimensions: {embedding.shape[1]}
    Cells analyzed: {embedding.shape[0]:,}
    """
    
    ax12.text(0.05, 0.95, summary_text, transform=ax12.transAxes, 
             verticalalignment='top', fontsize=10, family='monospace')
    
    plt.tight_layout()
    
    if save_plots:
        plot_filename = f"{PLOTS_DIR}/lemur_comprehensive_analysis.png"
        plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"   Plots saved: {plot_filename}")
    
    return fig

def plot_metadata_distributions(embedding, adata, fig):
    """Plot additional metadata distributions"""
    # Sample distribution
    ax9 = plt.subplot(3, 4, 9)
    if 'sample_id' in adata.obs.columns:
        samples = adata.obs['sample_id'].unique()[:10]
        sample_colors = plt.cm.tab10(np.linspace(0, 1, len(samples)))
        for i, sample in enumerate(samples):
            mask = adata.obs['sample_id'] == sample
            plt.scatter(embedding[mask, 0], embedding[mask, 1], 
                       c=[sample_colors[i]], label=sample, alpha=0.6, s=5)
        plt.xlabel('LEMUR Embedding 1')
        plt.ylabel('LEMUR Embedding 2')
        plt.title('Sample Distribution')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Sex distribution
    ax10 = plt.subplot(3, 4, 10)
    if 'sex' in adata.obs.columns:
        sex_colors = {'M': 'lightblue', 'F': 'pink'}
        for sex in adata.obs['sex'].unique():
            mask = adata.obs['sex'] == sex
            plt.scatter(embedding[mask, 0], embedding[mask, 1], 
                       c=sex_colors.get(sex, 'gray'), label=sex, alpha=0.6, s=5)
        plt.xlabel('LEMUR Embedding 1')
        plt.ylabel('LEMUR Embedding 2')
        plt.title('Sex Distribution')
        plt.legend()
    
    # Region distribution
    ax11 = plt.subplot(3, 4, 11)
    if 'region' in adata.obs.columns:
        region_colors = {'Putamen': 'orange', 'Caudate': 'purple'}
        for region in adata.obs['region'].unique():
            mask = adata.obs['region'] == region
            plt.scatter(embedding[mask, 0], embedding[mask, 1], 
                       c=region_colors.get(region, 'gray'), label=region, alpha=0.6, s=5)
        plt.xlabel('LEMUR Embedding 1')
        plt.ylabel('LEMUR Embedding 2')
        plt.title('Brain Region Distribution')
        plt.legend()

# ============================================================================
# 💾 EXPORT FUNCTIONS
# ============================================================================

def export_lemur_results(results, adata):
    """Export LEMUR results to files"""
    print(f"\n💾 EXPORTING LEMUR RESULTS")
    print("=" * 40)
    
    de_df = results['differential_results']
    embedding = results['embedding']
    
    # Export differential expression results
    de_filename = f"{TABLES_DIR}/lemur_differential_expression.csv"
    de_df.to_csv(de_filename, index=False)
    print(f"   DE results saved: {de_filename}")
    
    # Export embedding coordinates
    embedding_df = pd.DataFrame(
        embedding,
        columns=[f'LEMUR_{i+1}' for i in range(embedding.shape[1])],
        index=adata.obs.index
    )
    
    # Add metadata
    metadata_cols = ['condition', 'cell_type']
    for col in ['sample_id', 'sex', 'region']:
        if col in adata.obs.columns:
            metadata_cols.append(col)
    
    for col in metadata_cols:
        if col in adata.obs.columns:
            embedding_df[col] = adata.obs[col].values
    
    embedding_filename = f"{TABLES_DIR}/lemur_embedding_coordinates.csv"
    embedding_df.to_csv(embedding_filename)
    print(f"   Embedding saved: {embedding_filename}")
    
    # Export significant genes
    sig_genes = de_df[de_df['significant']].copy().sort_values('adj_pval')
    sig_filename = f"{TABLES_DIR}/lemur_significant_genes.csv"
    sig_genes.to_csv(sig_filename, index=False)
    print(f"   Significant genes saved: {sig_filename}")
    
    # Export top genes by effect size
    abs_lfc = de_df['lfc'].abs()
    top_genes = de_df.loc[abs_lfc.nlargest(100).index]
    top_filename = f"{TABLES_DIR}/lemur_top_genes_by_effect.csv"
    top_genes.to_csv(top_filename, index=False)
    print(f"   Top genes by effect saved: {top_filename}")
    
    # Create summary report
    n_sig = de_df['significant'].sum()
    n_sig_strict = de_df.get('significant_strict', pd.Series([False]*len(de_df))).sum()
    
    summary_data = {
        'Analysis': 'LEMUR',
        'Method': results['method'],
        'Total_Genes': len(de_df),
        'Significant_Genes': n_sig,
        'Significant_Genes_Strict': n_sig_strict,
        'Percent_Significant': n_sig / len(de_df) * 100,
        'Upregulated': ((de_df['significant']) & (de_df['lfc'] > 0)).sum(),
        'Downregulated': ((de_df['significant']) & (de_df['lfc'] < 0)).sum(),
        'Embedding_Dimensions': embedding.shape[1],
        'Cells_Analyzed': embedding.shape[0],
        'Mean_Abs_LFC': de_df['lfc'].abs().mean(),
        'Max_Abs_LFC': de_df['lfc'].abs().max()
    }
    
    summary_df = pd.DataFrame([summary_data])
    summary_filename = f"{TABLES_DIR}/lemur_analysis_summary.csv"
    summary_df.to_csv(summary_filename, index=False)
    print(f"   Summary saved: {summary_filename}")
    
    return {
        'de_results': de_filename,
        'embedding': embedding_filename, 
        'significant_genes': sig_filename,
        'top_genes': top_filename,
        'summary': summary_filename
    }

# ============================================================================
# 🚀 MAIN EXECUTION
# ============================================================================

def main():
    """Main execution function"""
    print("🚀 STARTING LEMUR ANALYSIS PIPELINE")
    print("===================================")
    
    try:
        # Load and prepare data
        adata = load_and_prepare_data()
        if adata is None:
            return None
        
        print(f"\n📊 DATASET COMPOSITION:")
        print(f"   Total cells: {adata.n_obs:,}")
        print(f"   Total genes: {adata.n_vars:,}")
        print(f"   Conditions: {dict(adata.obs['condition'].value_counts())}")
        print(f"   Cell types: {adata.obs['cell_type'].nunique()}")
        if 'sample_id' in adata.obs.columns:
            print(f"   Samples: {adata.obs['sample_id'].nunique()}")
        
        # Run LEMUR analysis
        print(f"\n{'#'*60}")
        print(f"RUNNING LEMUR CONTINUOUS ANALYSIS")
        print(f"{'#'*60}")
        
        lemur_results = run_pylemur_analysis(adata)
        
        if lemur_results is not None:
            # Use the actual adata that was used for analysis
            adata_for_plotting = lemur_results.get('adata_used', adata)
            
            # Create visualizations
            plot_lemur_results(lemur_results, adata_for_plotting)
            
            # Export results
            export_files = export_lemur_results(lemur_results, adata_for_plotting)
            
            # Print final summary
            print("\n✅ LEMUR ANALYSIS COMPLETED!")
            print("=" * 40)
            print(f"Method used: {lemur_results['method']}")
            print(f"Results saved to: {OUTPUT_DIR}")
            
            de_df = lemur_results['differential_results']
            n_sig = de_df['significant'].sum()
            n_sig_strict = de_df.get('significant_strict', pd.Series([False]*len(de_df))).sum()
            
            print(f"Total significant genes (lenient): {n_sig}")
            print(f"Total significant genes (strict): {n_sig_strict}")
            
            # Show top genes
            if n_sig > 0:
                print(f"\nTop significant genes:")
                top_genes = de_df[de_df['significant']].nsmallest(5, 'adj_pval')
                for _, gene in top_genes.iterrows():
                    direction = "↑" if gene['lfc'] > 0 else "↓"
                    print(f"  {direction} {gene['gene']}: LFC={gene['lfc']:.3f}, padj={gene['adj_pval']:.2e}")
            
            print(f"\nTop genes by absolute effect size:")
            abs_lfc = de_df['lfc'].abs()
            top_by_effect = de_df.loc[abs_lfc.nlargest(5).index]
            for _, gene in top_by_effect.iterrows():
                direction = "↑" if gene['lfc'] > 0 else "↓"
                sig_status = "*" if gene['significant'] else ""
                print(f"  {direction} {gene['gene']}: LFC={gene['lfc']:.3f}, padj={gene['adj_pval']:.2e} {sig_status}")
        
        return lemur_results
        
    except Exception as e:
        print(f"\n❌ ERROR: {str(e)}")
        import traceback
        traceback.print_exc()
        return None

if __name__ == "__main__":
    # Set random seed for reproducibility
    np.random.seed(ANALYSIS_CONFIG['random_seed'])
    
    # Run main analysis
    results = main()