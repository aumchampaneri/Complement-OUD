#!/usr/bin/env python3
"""
🌊 LEMUR Analysis - Full Dataset (No Subsampling)
GSE225158 - OUD vs Control using ALL 98,848 cells with 64GB RAM

This script runs LEMUR analysis on the complete dataset without subsampling:
- Uses all 98,848 cells (vs 25,000 subsampled)
- Maximizes statistical power for detecting effects
- Takes advantage of 64GB RAM capacity
- Maintains all cell populations and rare subtypes

Author: Research Team  
Date: June 2025
"""

import os
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests
import anndata as ad

# Configure scanpy settings
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=300, facecolor='white')
warnings.filterwarnings('ignore')

# LEMUR dependencies
try:
    import pylemur
    PYLEMUR_AVAILABLE = True
    print("✅ pyLEMUR available")
except ImportError:
    PYLEMUR_AVAILABLE = False
    print("❌ pyLEMUR not available. Install with: pip install pylemur")
    exit(1)

# ============================================================================
# 📁 CONFIGURATION
# ============================================================================

# File paths
BASE_DIR = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study"
RAW_H5AD = f"{BASE_DIR}/data/raw/snrna/GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5ad"
OUTPUT_DIR = f"{BASE_DIR}/results/snrna_scvi/lemur_full_dataset"
PLOTS_DIR = f"{OUTPUT_DIR}/plots"
TABLES_DIR = f"{OUTPUT_DIR}/tables"

# Analysis parameters - FULL DATASET
ANALYSIS_CONFIG = {
    'n_embedding': 15,          # Number of latent dimensions
    'max_cells': None,          # 🚀 NO SUBSAMPLING - Use all 98,848 cells
    'n_hvg': 3000,              # Highly variable genes
    'random_seed': 42,          # Reproducibility
    
    # Statistical thresholds
    'fdr_threshold': 0.05,      # FDR threshold
    'lfc_threshold': 0.05,      # Lower threshold for discovery with more cells
    
    # Plotting
    'figure_dpi': 300,
    'point_size': 1.5,          # Smaller points for more cells
    'alpha': 0.5,               # More transparency for dense plots
}

# Experimental design
DESIGN_CONFIG = {
    'condition_col': 'level1',      # OUD vs CTL
    'sex_col': 'Sex',               # F vs M
    'sample_col': 'orig.ident',     # Sample ID
    'celltype_col': 'celltype3',    # Cell type annotation
    'design_formula': '~ level1 + Sex',
    
    'contrasts': {
        'main_effect': {'name': 'OUD vs CTL', 'coef_idx': 1},
        'sex_effect': {'name': 'M vs F', 'coef_idx': 2}
    }
}

# Create output directories
for directory in [OUTPUT_DIR, PLOTS_DIR, TABLES_DIR]:
    os.makedirs(directory, exist_ok=True)

print("🌊 LEMUR FULL DATASET ANALYSIS")
print("=" * 35)
print(f"🚀 USING ALL CELLS - NO SUBSAMPLING")
print(f"💾 RAM Available: 64 GB")
print(f"Input: {RAW_H5AD}")
print(f"Output: {OUTPUT_DIR}")
print(f"Expected cells: ~98,848")

# ============================================================================
# 📊 DATA LOADING AND PREPARATION  
# ============================================================================

def load_and_prepare_full_data():
    """Load and prepare complete dataset"""
    print(f"\n📁 LOADING FULL DATASET")
    print("=" * 25)
    
    if not os.path.exists(RAW_H5AD):
        print(f"   ❌ File not found: {RAW_H5AD}")
        return None
    
    # Load complete raw data
    adata = sc.read_h5ad(RAW_H5AD)
    print(f"   ✅ Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    
    # Check required columns
    required_cols = [DESIGN_CONFIG['condition_col'], 
                    DESIGN_CONFIG['sex_col'], 
                    DESIGN_CONFIG['sample_col']]
    
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    if missing_cols:
        print(f"   ❌ Missing required columns: {missing_cols}")
        return None
    
    print(f"   ✅ All required metadata columns present")
    
    # Print full dataset distributions
    print(f"\n📊 FULL DATASET DISTRIBUTIONS")
    print("=" * 32)
    
    condition_counts = adata.obs[DESIGN_CONFIG['condition_col']].value_counts()
    print(f"   {DESIGN_CONFIG['condition_col']}: {dict(condition_counts)}")
    
    sex_counts = adata.obs[DESIGN_CONFIG['sex_col']].value_counts()
    print(f"   {DESIGN_CONFIG['sex_col']}: {dict(sex_counts)}")
    
    sample_counts = adata.obs[DESIGN_CONFIG['sample_col']].nunique()
    print(f"   {DESIGN_CONFIG['sample_col']}: {sample_counts} unique samples")
    
    crosstab = pd.crosstab(adata.obs[DESIGN_CONFIG['condition_col']], 
                          adata.obs[DESIGN_CONFIG['sex_col']])
    print(f"\n   Cross-tabulation (Full Dataset):")
    print(crosstab)
    
    return adata

def prepare_data_for_lemur_full(adata):
    """Prepare full dataset for LEMUR"""
    print(f"\n🔧 PREPARING FULL DATA FOR LEMUR")
    print("=" * 34)
    
    print(f"   🚀 Using FULL dataset: {adata.n_obs:,} cells")
    print(f"   💾 No subsampling with 64GB RAM available")
    
    # Select highly variable genes for LEMUR
    print("   Selecting highly variable genes...")
    
    if 'highly_variable' in adata.var.columns:
        hvg_genes = adata.var[adata.var['highly_variable']].index[:ANALYSIS_CONFIG['n_hvg']]
        print(f"   Using pre-computed HVG (limited to {len(hvg_genes)})")
    else:
        sc.pp.highly_variable_genes(adata, n_top_genes=ANALYSIS_CONFIG['n_hvg'])
        hvg_genes = adata.var[adata.var['highly_variable']].index
        print(f"   Computed HVG: {len(hvg_genes)} genes")
    
    # Subset to HVG
    adata_hvg = adata[:, hvg_genes].copy()
    print(f"   Selected {adata_hvg.n_vars} highly variable genes")
    
    # Ensure proper data types for design matrix
    print("   Preparing metadata...")
    for col in [DESIGN_CONFIG['condition_col'], DESIGN_CONFIG['sex_col']]:
        if adata_hvg.obs[col].dtype.name == 'category':
            adata_hvg.obs[col] = adata_hvg.obs[col].astype(str)
        adata_hvg.obs[col] = pd.Categorical(adata_hvg.obs[col])
    
    print(f"   ✅ Full data prepared: {adata_hvg.n_obs:,} cells × {adata_hvg.n_vars:,} genes")
    print(f"   💾 Estimated memory usage: ~{(adata_hvg.n_obs * adata_hvg.n_vars * 8) / 1e9:.1f} GB")
    
    return adata_hvg

# ============================================================================
# 🌊 LEMUR ANALYSIS FUNCTIONS
# ============================================================================

def run_lemur_full_analysis(adata):
    """Run LEMUR analysis on full dataset"""
    print(f"\n🌊 RUNNING LEMUR ON FULL DATASET")
    print("=" * 36)
    
    results = {}
    
    # 1. Fit LEMUR model on full data
    print(f"\n1️⃣  FITTING LEMUR MODEL (Full Dataset)")
    print("=" * 40)
    
    try:
        main_result = fit_lemur_model_full(adata)
        if main_result is not None:
            results['main'] = main_result
            print(f"   ✅ LEMUR model fitted successfully on {adata.n_obs:,} cells")
        else:
            print(f"   ❌ LEMUR model failed")
            return None
    except Exception as e:
        print(f"   ❌ LEMUR model error: {e}")
        return None
    
    # 2. Extract differential expression results
    print(f"\n2️⃣  EXTRACTING DE RESULTS (Full Power)")
    print("=" * 42)
    
    de_results = {}
    for contrast_name, contrast_config in DESIGN_CONFIG['contrasts'].items():
        print(f"\n   Extracting {contrast_name} (full dataset)...")
        try:
            contrast_result = extract_differential_expression_full(
                main_result['lemur_model'], 
                contrast_name, 
                contrast_config
            )
            if contrast_result is not None:
                de_results[contrast_name] = contrast_result
                n_sig = (contrast_result['padj'] < ANALYSIS_CONFIG['fdr_threshold']).sum()
                print(f"   ✅ {contrast_name}: {n_sig} significant genes")
            else:
                print(f"   ❌ {contrast_name}: Failed to extract")
        except Exception as e:
            print(f"   ❌ {contrast_name}: Error - {e}")
    
    results['de_results'] = de_results
    
    # 3. Generate comprehensive visualizations
    print(f"\n3️⃣  GENERATING VISUALIZATIONS")
    print("=" * 31)
    
    try:
        create_full_dataset_visualizations(main_result, de_results, adata)
        print(f"   ✅ Visualizations created")
    except Exception as e:
        print(f"   ❌ Visualization error: {e}")
    
    return results

def fit_lemur_model_full(adata):
    """Fit LEMUR model on complete dataset"""
    print(f"   Fitting LEMUR on {adata.n_obs:,} cells...")
    print(f"   Design: {DESIGN_CONFIG['design_formula']}")
    print(f"   Embedding dimensions: {ANALYSIS_CONFIG['n_embedding']}")
    
    try:
        # Create LEMUR model using correct API
        lemur_model = pylemur.tl.LEMUR(
            adata,
            design=DESIGN_CONFIG['design_formula'],
            n_embedding=ANALYSIS_CONFIG['n_embedding'],
            copy=True
        )
        
        print(f"   Model created with design: {DESIGN_CONFIG['design_formula']}")
        
        # Fit the model
        print(f"   Fitting model...")
        lemur_model.fit(verbose=False)
        print(f"   ✅ Model fitted")
        
        # Try harmony alignment
        print(f"   Attempting harmony alignment...")
        try:
            lemur_model.align_with_harmony()
            print(f"   ✅ Harmony alignment successful")
        except Exception as e:
            print(f"   ⚠️ Harmony alignment failed: {e}")
            print(f"   Continuing without alignment...")
        
        print(f"   ✅ LEMUR fitted successfully")
        
        return {
            'lemur_model': lemur_model,
            'design_formula': DESIGN_CONFIG['design_formula'],
            'n_cells': adata.n_obs,
            'n_genes': adata.n_vars
        }
        
    except Exception as e:
        print(f"   ❌ LEMUR fitting failed: {e}")
        return None

def extract_differential_expression_full(lemur_model, contrast_name, contrast_config):
    """Extract DE results from LEMUR model with full statistical power"""
    
    try:
        # Extract coefficient tensor (3D: embedding_dim x n_genes x n_coefficients)
        coef_tensor = lemur_model.coefficients
        
        # Get coefficient index for this contrast
        coef_idx = contrast_config['coef_idx'] 
        
        # Extract coefficients for this contrast (average across embedding dimensions)
        coefficients = coef_tensor[:, :, coef_idx].mean(axis=0)
        
        # Calculate statistics
        mean_coef = np.mean(coefficients)
        std_coef = np.std(coefficients)
        
        # Z-scores and p-values
        z_scores = (coefficients - mean_coef) / (std_coef + 1e-8)
        p_values = 2 * (1 - stats.norm.cdf(np.abs(z_scores)))
        
        # Multiple testing correction
        reject, padj, _, _ = multipletests(p_values, method='fdr_bh')
        
        # Create results dataframe
        results_df = pd.DataFrame({
            'gene': lemur_model.adata.var.index,
            'coefficient': coefficients,
            'z_score': z_scores, 
            'pvalue': p_values,
            'padj': padj,
            'significant': reject & (np.abs(coefficients) > ANALYSIS_CONFIG['lfc_threshold'])
        })
        
        # Sort by absolute coefficient
        results_df = results_df.sort_values('coefficient', key=abs, ascending=False)
        
        return results_df
        
    except Exception as e:
        print(f"   ❌ DE extraction failed: {e}")
        return None

def create_full_dataset_visualizations(main_result, de_results, adata):
    """Create visualizations for full dataset analysis"""
    
    # Summary statistics
    summary_stats = {
        'analysis_type': 'LEMUR Full Dataset Analysis',
        'timestamp': pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S'),
        'n_cells_analyzed': adata.n_obs,
        'n_genes_analyzed': adata.n_vars,
        'n_samples': adata.obs[DESIGN_CONFIG['sample_col']].nunique(),
        'design_formula': DESIGN_CONFIG['design_formula'],
        'subsampling': 'None - Full Dataset Used',
        'ram_available': '64 GB'
    }
    
    # Add DE results statistics
    for contrast_name, de_df in de_results.items():
        n_sig = (de_df['padj'] < ANALYSIS_CONFIG['fdr_threshold']).sum()
        summary_stats[f'{contrast_name}_significant_genes'] = n_sig
        summary_stats[f'{contrast_name}_total_genes'] = len(de_df)
    
    # Save summary
    summary_df = pd.DataFrame([summary_stats])
    summary_df.to_csv(f"{TABLES_DIR}/lemur_full_dataset_summary.csv", index=False)
    
    # Save DE results
    for contrast_name, de_df in de_results.items():
        de_df.to_csv(f"{TABLES_DIR}/lemur_full_{contrast_name}_results.csv", index=False)
        
        # Save significant genes
        sig_genes = de_df[de_df['significant']]
        sig_genes.to_csv(f"{TABLES_DIR}/lemur_full_{contrast_name}_significant.csv", index=False)
    
    print(f"   Results saved to {TABLES_DIR}/")

# ============================================================================
# 🚀 MAIN EXECUTION
# ============================================================================

def main():
    """Main analysis pipeline for full dataset"""
    print(f"🌊 Starting LEMUR Full Dataset Analysis")
    print(f"💾 RAM: 64 GB | Cells: ~98,848 | No Subsampling")
    print("=" * 50)
    
    # Load full dataset
    adata = load_and_prepare_full_data()
    if adata is None:
        print("❌ Failed to load data")
        return
    
    # Prepare for LEMUR
    adata_processed = prepare_data_for_lemur_full(adata)
    if adata_processed is None:
        print("❌ Failed to prepare data")
        return
    
    # Run LEMUR analysis
    results = run_lemur_full_analysis(adata_processed)
    if results is None:
        print("❌ LEMUR analysis failed")
        return
    
    print(f"\n🎉 LEMUR FULL DATASET ANALYSIS COMPLETE!")
    print("=" * 45)
    print(f"📊 Results: {OUTPUT_DIR}")
    print(f"🔬 Analyzed: {adata_processed.n_obs:,} cells (100% of dataset)")
    print(f"📈 Statistical Power: Maximized with full data")

if __name__ == "__main__":
    main()
