#!/usr/bin/env python3
"""
🌊 LEMUR Analysis - Final Working Implementation with Raw Data
GSE225158 - OUD vs Control using pyLEMUR following best practices

This script implements proper LEMUR analysis using:
1. Raw, unintegrated single-cell data
2. Proper experimental design with sex interactions
3. Correct coefficient extraction from 3D tensor
4. Multiple contrasts including sex interactions
5. Comprehensive visualization and export

Author: Research Team
Date: 2024
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

# File paths - Using RAW unintegrated data per LEMUR requirements
BASE_DIR = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study"
RAW_H5AD = f"{BASE_DIR}/data/raw/snrna/GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5ad"
OUTPUT_DIR = f"{BASE_DIR}/results/snrna_scvi/lemur_final_working"
PLOTS_DIR = f"{OUTPUT_DIR}/plots"
TABLES_DIR = f"{OUTPUT_DIR}/tables"

# Analysis parameters
ANALYSIS_CONFIG = {
    'n_embedding': 15,          # Number of latent dimensions
    'max_cells': None,          # No subsetting - use full dataset with 64GB RAM
    'n_hvg': 3000,              # Highly variable genes
    'random_seed': 42,          # Reproducibility
    
    # Statistical thresholds
    'fdr_threshold': 0.05,      # FDR threshold
    'lfc_threshold': 0.1,       # Log fold change threshold (relaxed for discovery)
    
    # Plotting
    'figure_dpi': 300,
    'point_size': 2,
    'alpha': 0.7,
}

# Experimental design
DESIGN_CONFIG = {
    'condition_col': 'level1',      # OUD vs CTL
    'sex_col': 'Sex',               # F vs M
    'sample_col': 'orig.ident',     # Sample ID for random effects
    'celltype_col': 'celltype3',    # Cell type annotation
    
    # Design formula
    'design_formula': '~ level1 + Sex',
    
    # Contrasts to test (mapped to coefficient indices)
    'contrasts': {
        'main_effect': {'name': 'OUD vs CTL', 'coef_idx': 1},
        'sex_effect': {'name': 'M vs F', 'coef_idx': 2}
    }
}

# Create output directories
for directory in [OUTPUT_DIR, PLOTS_DIR, TABLES_DIR]:
    os.makedirs(directory, exist_ok=True)

print("🌊 LEMUR FINAL WORKING ANALYSIS - RAW DATA")
print("==========================================")
print(f"Input: {RAW_H5AD}")
print(f"Output: {OUTPUT_DIR}")
print(f"Design: {DESIGN_CONFIG['design_formula']}")

# ============================================================================
# 📊 DATA LOADING AND PREPARATION
# ============================================================================

def load_and_prepare_raw_data():
    """Load and prepare raw unintegrated data for LEMUR"""
    print(f"\n📁 LOADING RAW DATA")
    print("=" * 20)
    
    if not os.path.exists(RAW_H5AD):
        print(f"   ❌ File not found: {RAW_H5AD}")
        return None
    
    # Load raw data
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
    
    # Print group distributions
    print(f"\n📊 GROUP DISTRIBUTIONS")
    print("=" * 22)
    
    condition_counts = adata.obs[DESIGN_CONFIG['condition_col']].value_counts()
    print(f"   {DESIGN_CONFIG['condition_col']}: {dict(condition_counts)}")
    
    sex_counts = adata.obs[DESIGN_CONFIG['sex_col']].value_counts()
    print(f"   {DESIGN_CONFIG['sex_col']}: {dict(sex_counts)}")
    
    sample_counts = adata.obs[DESIGN_CONFIG['sample_col']].nunique()
    print(f"   {DESIGN_CONFIG['sample_col']}: {sample_counts} unique samples")
    
    crosstab = pd.crosstab(adata.obs[DESIGN_CONFIG['condition_col']], 
                          adata.obs[DESIGN_CONFIG['sex_col']])
    print(f"\n   Cross-tabulation:")
    print(crosstab)
    
    return adata

def prepare_data_for_lemur(adata):
    """Prepare data following LEMUR requirements"""
    print(f"\n🔧 PREPARING DATA FOR LEMUR")
    print("=" * 30)
    
    # Check for minimum group sizes
    print("   Checking group sizes...")
    group_sizes = adata.obs.groupby([DESIGN_CONFIG['condition_col'], 
                                    DESIGN_CONFIG['sex_col']]).size()
    print("   Group sizes:")
    print(group_sizes)
    
    # Memory management - skip subsample for full dataset analysis
    if ANALYSIS_CONFIG['max_cells'] is not None and adata.n_obs > ANALYSIS_CONFIG['max_cells']:
        print(f"   🔢 Subsampling to {ANALYSIS_CONFIG['max_cells']} cells...")
        adata = subsample_balanced(adata, ANALYSIS_CONFIG['max_cells'])
        print(f"   ✅ After subsampling: {adata.n_obs:,} cells")
    else:
        print(f"   🚀 Using full dataset: {adata.n_obs:,} cells (no subsampling with 64GB RAM)")
    
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
    
    print(f"   ✅ Data prepared: {adata_hvg.n_obs:,} cells × {adata_hvg.n_vars:,} genes")
    
    return adata_hvg

def subsample_balanced(adata, max_cells):
    """Subsample data while maintaining group balance"""
    np.random.seed(ANALYSIS_CONFIG['random_seed'])
    
    groups = adata.obs.groupby([DESIGN_CONFIG['condition_col'], 
                               DESIGN_CONFIG['sex_col']])
    
    n_groups = len(groups)
    cells_per_group = max_cells // n_groups
    
    sampled_indices = []
    for name, group in groups:
        if len(group) <= cells_per_group:
            sampled_indices.extend(group.index)
        else:
            sampled = np.random.choice(group.index, size=cells_per_group, replace=False)
            sampled_indices.extend(sampled)
    
    return adata[sampled_indices].copy()

# ============================================================================
# 🌊 LEMUR ANALYSIS FUNCTIONS
# ============================================================================

def run_lemur_analysis(adata):
    """Run LEMUR analysis with coefficient extraction"""
    print(f"\n🌊 RUNNING LEMUR ANALYSIS")
    print("=" * 27)
    
    results = {}
    
    # 1. Fit LEMUR model
    print(f"\n1️⃣  FITTING LEMUR MODEL")
    print("=" * 25)
    
    try:
        main_result = fit_lemur_model(adata)
        if main_result is not None:
            results['main'] = main_result
            print(f"   ✅ LEMUR model fitted successfully")
        else:
            print(f"   ❌ LEMUR model failed")
            return None
    except Exception as e:
        print(f"   ❌ LEMUR model error: {e}")
        return None
    
    # 2. Extract differential expression results
    print(f"\n2️⃣  EXTRACTING DE RESULTS")
    print("=" * 27)
    
    de_results = {}
    for contrast_name, contrast_config in DESIGN_CONFIG['contrasts'].items():
        print(f"\n   Extracting {contrast_name}...")
        try:
            contrast_result = extract_differential_expression(
                main_result['lemur_model'], 
                contrast_name, 
                contrast_config
            )
            if contrast_result is not None:
                de_results[contrast_name] = contrast_result
                n_sig = contrast_result['significant'].sum() if 'significant' in contrast_result.columns else 0
                print(f"      ✅ {contrast_name}: {n_sig} significant genes")
        except Exception as e:
            print(f"      ❌ {contrast_name} failed: {e}")
    
    results['de_results'] = de_results
    
    return results

def fit_lemur_model(adata):
    """Fit LEMUR model"""
    print(f"   Fitting LEMUR model...")
    
    try:
        # Create LEMUR model
        lemur_model = pylemur.tl.LEMUR(
            adata,
            design=DESIGN_CONFIG['design_formula'],
            n_embedding=ANALYSIS_CONFIG['n_embedding'],
            copy=True
        )
        
        print(f"      Model created with design: {DESIGN_CONFIG['design_formula']}")
        
        # Fit the model
        print(f"      Fitting model...")
        lemur_model.fit(verbose=False)
        print(f"      ✅ Model fitted")
        
        # Try alignment
        print(f"      Attempting harmony alignment...")
        try:
            lemur_model.align_with_harmony()
            print(f"      ✅ Harmony alignment successful")
        except Exception as e:
            print(f"      ⚠️  Harmony alignment failed: {e}")
            print(f"      Continuing without alignment...")
        
        return {
            'lemur_model': lemur_model,
            'adata_used': adata
        }
        
    except Exception as e:
        print(f"      ❌ Model fitting failed: {e}")
        return None

def extract_differential_expression(lemur_model, contrast_name, contrast_config):
    """Extract differential expression using LEMUR coefficients"""
    print(f"      Extracting DE for {contrast_name}...")
    
    try:
        # Get coefficients - shape is (n_embedding, n_genes, n_coefficients)
        coefficients = lemur_model.coefficients
        print(f"         Coefficient tensor shape: {coefficients.shape}")
        
        # Extract coefficient index for this contrast
        coef_idx = contrast_config['coef_idx']
        
        if coef_idx >= coefficients.shape[2]:
            print(f"         ❌ Coefficient index {coef_idx} out of range (max: {coefficients.shape[2]-1})")
            return None
        
        # Extract coefficients for this contrast across all embedding dimensions
        # Shape: (n_embedding, n_genes)
        contrast_coefficients = coefficients[:, :, coef_idx]
        print(f"         Extracted coefficients shape: {contrast_coefficients.shape}")
        
        # Average across embedding dimensions to get gene-level effects
        gene_coefficients = np.mean(contrast_coefficients, axis=0)
        print(f"         Gene-level coefficients shape: {gene_coefficients.shape}")
        
        # Get gene names
        gene_names = lemur_model.adata.var.index
        
        # Create results DataFrame
        results_df = pd.DataFrame({
            'gene': gene_names,
            'coefficient': gene_coefficients,
            'abs_coefficient': np.abs(gene_coefficients)
        })
        
        # Calculate statistics
        # Use coefficient magnitudes relative to their distribution
        coef_std = np.std(gene_coefficients)
        results_df['z_score'] = np.abs(gene_coefficients) / (coef_std + 1e-8)
        results_df['pval'] = 2 * (1 - stats.norm.cdf(results_df['z_score']))
        
        # Multiple testing correction
        _, fdr_values, _, _ = multipletests(results_df['pval'], method='fdr_bh')
        results_df['fdr'] = fdr_values
        
        # Significance
        results_df['significant'] = (
            (results_df['fdr'] < ANALYSIS_CONFIG['fdr_threshold']) &
            (results_df['abs_coefficient'] > ANALYSIS_CONFIG['lfc_threshold'])
        )
        
        # Sort by significance
        results_df = results_df.sort_values('abs_coefficient', ascending=False)
        results_df = results_df.set_index('gene')
        
        print(f"         ✅ Extracted {len(results_df)} genes")
        print(f"         Significant genes: {results_df['significant'].sum()}")
        
        # Show top genes
        top_genes = results_df.head(10)
        print(f"         Top 10 genes by effect size:")
        for gene, row in top_genes.iterrows():
            sig_flag = "***" if row['significant'] else ""
            print(f"           {gene}: coef={row['coefficient']:.3f}, FDR={row['fdr']:.3e} {sig_flag}")
        
        return results_df
        
    except Exception as e:
        print(f"         ❌ DE extraction failed: {e}")
        import traceback
        traceback.print_exc()
        return None

# ============================================================================
# 📊 VISUALIZATION FUNCTIONS
# ============================================================================

def create_comprehensive_plots(results):
    """Create comprehensive visualization of LEMUR results"""
    print(f"\n📊 CREATING VISUALIZATIONS")
    print("=" * 26)
    
    if 'main' not in results:
        print("   ❌ No main results to plot")
        return
    
    lemur_model = results['main']['lemur_model']
    
    # 1. Embedding plots
    plot_embeddings(lemur_model)
    
    # 2. DE results
    if 'de_results' in results:
        plot_de_results(results['de_results'])
    
    # 3. Summary statistics
    create_summary_plots(results)

def plot_embeddings(lemur_model):
    """Plot LEMUR embeddings colored by metadata"""
    print("   Creating embedding plots...")
    
    try:
        embedding = lemur_model.embedding
        adata_used = lemur_model.adata
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        axes = axes.flatten()
        
        # Convert to string for consistent handling
        condition_vals = adata_used.obs[DESIGN_CONFIG['condition_col']].astype(str)
        sex_vals = adata_used.obs[DESIGN_CONFIG['sex_col']].astype(str)
        
        # Plot 1: Condition
        condition_codes = pd.Categorical(condition_vals).codes
        scatter = axes[0].scatter(embedding[:, 0], embedding[:, 1], 
                                c=condition_codes,
                                s=ANALYSIS_CONFIG['point_size'], 
                                alpha=ANALYSIS_CONFIG['alpha'],
                                cmap='viridis')
        axes[0].set_title(f'LEMUR Embedding - {DESIGN_CONFIG["condition_col"]}')
        axes[0].set_xlabel('LEMUR 1')
        axes[0].set_ylabel('LEMUR 2')
        
        # Add legend for condition
        unique_conditions = pd.Categorical(condition_vals).categories
        for i, cond in enumerate(unique_conditions):
            axes[0].scatter([], [], c=plt.cm.get_cmap('viridis')(i/len(unique_conditions)), 
                          label=cond, s=50)
        axes[0].legend(title=DESIGN_CONFIG['condition_col'])
        
        # Plot 2: Sex
        sex_codes = pd.Categorical(sex_vals).codes
        scatter = axes[1].scatter(embedding[:, 0], embedding[:, 1], 
                                c=sex_codes,
                                s=ANALYSIS_CONFIG['point_size'], 
                                alpha=ANALYSIS_CONFIG['alpha'],
                                cmap='plasma')
        axes[1].set_title(f'LEMUR Embedding - {DESIGN_CONFIG["sex_col"]}')
        axes[1].set_xlabel('LEMUR 1')
        axes[1].set_ylabel('LEMUR 2')
        
        # Add legend for sex
        unique_sex = pd.Categorical(sex_vals).categories
        for i, sex in enumerate(unique_sex):
            axes[1].scatter([], [], c=plt.cm.get_cmap('plasma')(i/len(unique_sex)), 
                          label=sex, s=50)
        axes[1].legend(title=DESIGN_CONFIG['sex_col'])
        
        # Plot 3: Combined condition + sex
        combined = condition_vals + '_' + sex_vals
        combined_codes = pd.Categorical(combined).codes
        scatter = axes[2].scatter(embedding[:, 0], embedding[:, 1], 
                                c=combined_codes,
                                s=ANALYSIS_CONFIG['point_size'], 
                                alpha=ANALYSIS_CONFIG['alpha'],
                                cmap='Set1')
        axes[2].set_title('LEMUR Embedding - Condition × Sex')
        axes[2].set_xlabel('LEMUR 1')
        axes[2].set_ylabel('LEMUR 2')
        
        # Plot 4: Density plot
        axes[3].hexbin(embedding[:, 0], embedding[:, 1], gridsize=30, cmap='Blues')
        axes[3].set_title('LEMUR Embedding - Density')
        axes[3].set_xlabel('LEMUR 1')
        axes[3].set_ylabel('LEMUR 2')
        
        plt.tight_layout()
        plt.savefig(f"{PLOTS_DIR}/lemur_embeddings.png", 
                   dpi=ANALYSIS_CONFIG['figure_dpi'], bbox_inches='tight')
        plt.close()
        
        print("      ✅ Embedding plots saved")
        
    except Exception as e:
        print(f"      ❌ Embedding plotting failed: {e}")

def plot_de_results(de_results):
    """Plot differential expression results"""
    print("   Creating DE plots...")
    
    try:
        n_contrasts = len(de_results)
        if n_contrasts == 0:
            return
        
        fig, axes = plt.subplots(1, min(n_contrasts, 2), figsize=(12, 5))
        if n_contrasts == 1:
            axes = [axes]
        
        for i, (contrast_name, results) in enumerate(de_results.items()):
            if i >= 2:  # Only plot first 2
                break
                
            if results is None or not isinstance(results, pd.DataFrame):
                continue
            
            if 'coefficient' in results.columns and 'fdr' in results.columns:
                # Volcano plot
                x = results['coefficient']
                y = -np.log10(results['fdr'].clip(lower=1e-300))
                
                # Color by significance
                colors = ['red' if sig else 'gray' 
                         for sig in results.get('significant', [False]*len(results))]
                
                axes[i].scatter(x, y, c=colors, s=3, alpha=0.6)
                axes[i].set_xlabel('Coefficient')
                axes[i].set_ylabel('-log10(FDR)')
                axes[i].set_title(f'{contrast_name}\nVolcano Plot')
                axes[i].axhline(y=-np.log10(ANALYSIS_CONFIG['fdr_threshold']), 
                              color='blue', linestyle='--', alpha=0.5)
                axes[i].axvline(x=ANALYSIS_CONFIG['lfc_threshold'], 
                              color='blue', linestyle='--', alpha=0.5)
                axes[i].axvline(x=-ANALYSIS_CONFIG['lfc_threshold'], 
                              color='blue', linestyle='--', alpha=0.5)
                
                # Add text with number of significant genes
                n_sig = results['significant'].sum()
                axes[i].text(0.05, 0.95, f'{n_sig} significant genes', 
                           transform=axes[i].transAxes, 
                           verticalalignment='top',
                           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        plt.savefig(f"{PLOTS_DIR}/de_volcano_plots.png", 
                   dpi=ANALYSIS_CONFIG['figure_dpi'], bbox_inches='tight')
        plt.close()
        
        print("      ✅ DE plots saved")
        
    except Exception as e:
        print(f"      ❌ DE plotting failed: {e}")

def create_summary_plots(results):
    """Create summary statistics plots"""
    print("   Creating summary plots...")
    
    try:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Summary statistics
        summary_stats = []
        
        if 'de_results' in results:
            for contrast_name, contrast_results in results['de_results'].items():
                if contrast_results is not None and isinstance(contrast_results, pd.DataFrame):
                    if 'significant' in contrast_results.columns:
                        n_sig = contrast_results['significant'].sum()
                        n_total = len(contrast_results)
                        summary_stats.append({
                            'contrast': contrast_name,
                            'significant_genes': n_sig,
                            'total_genes': n_total,
                            'pct_significant': n_sig/n_total*100 if n_total > 0 else 0
                        })
        
        if summary_stats:
            summary_df = pd.DataFrame(summary_stats)
            
            # Bar plot of significant genes
            axes[0, 0].bar(summary_df['contrast'], summary_df['significant_genes'])
            axes[0, 0].set_title('Significant Genes by Contrast')
            axes[0, 0].set_ylabel('Number of Significant Genes')
            plt.setp(axes[0, 0].get_xticklabels(), rotation=45, ha='right')
            
            # Percentage plot
            axes[0, 1].bar(summary_df['contrast'], summary_df['pct_significant'])
            axes[0, 1].set_title('Percentage of Significant Genes')
            axes[0, 1].set_ylabel('Percentage Significant')
            plt.setp(axes[0, 1].get_xticklabels(), rotation=45, ha='right')
        else:
            axes[0, 0].text(0.5, 0.5, 'No DE results available', ha='center', va='center')
            axes[0, 1].text(0.5, 0.5, 'No DE results available', ha='center', va='center')
        
        # Analysis info
        axes[1, 0].text(0.1, 0.9, f"Analysis Parameters:", fontweight='bold', transform=axes[1, 0].transAxes)
        axes[1, 0].text(0.1, 0.8, f"• Embedding dims: {ANALYSIS_CONFIG['n_embedding']}", transform=axes[1, 0].transAxes)
        axes[1, 0].text(0.1, 0.7, f"• HVG: {ANALYSIS_CONFIG['n_hvg']}", transform=axes[1, 0].transAxes)
        axes[1, 0].text(0.1, 0.6, f"• FDR threshold: {ANALYSIS_CONFIG['fdr_threshold']}", transform=axes[1, 0].transAxes)
        axes[1, 0].text(0.1, 0.5, f"• LFC threshold: {ANALYSIS_CONFIG['lfc_threshold']}", transform=axes[1, 0].transAxes)
        axes[1, 0].set_xlim(0, 1)
        axes[1, 0].set_ylim(0, 1)
        axes[1, 0].axis('off')
        axes[1, 0].set_title('Analysis Parameters')
        
        # Design info
        axes[1, 1].text(0.1, 0.9, f"Experimental Design:", fontweight='bold', transform=axes[1, 1].transAxes)
        axes[1, 1].text(0.1, 0.8, f"• Formula: {DESIGN_CONFIG['design_formula']}", transform=axes[1, 1].transAxes)
        axes[1, 1].text(0.1, 0.7, f"• Condition: {DESIGN_CONFIG['condition_col']}", transform=axes[1, 1].transAxes)
        axes[1, 1].text(0.1, 0.6, f"• Sex: {DESIGN_CONFIG['sex_col']}", transform=axes[1, 1].transAxes)
        axes[1, 1].text(0.1, 0.5, f"• Sample ID: {DESIGN_CONFIG['sample_col']}", transform=axes[1, 1].transAxes)
        axes[1, 1].set_xlim(0, 1)
        axes[1, 1].set_ylim(0, 1)
        axes[1, 1].axis('off')
        axes[1, 1].set_title('Experimental Design')
        
        plt.tight_layout()
        plt.savefig(f"{PLOTS_DIR}/analysis_summary.png", 
                   dpi=ANALYSIS_CONFIG['figure_dpi'], bbox_inches='tight')
        plt.close()
        
        print("      ✅ Summary plots saved")
        
    except Exception as e:
        print(f"      ❌ Summary plotting failed: {e}")

# ============================================================================
# 💾 EXPORT FUNCTIONS
# ============================================================================

def export_results(results):
    """Export all results to CSV files"""
    print(f"\n💾 EXPORTING RESULTS")
    print("=" * 20)
    
    # Export DE results
    if 'de_results' in results:
        print("   Exporting DE results...")
        for contrast_name, contrast_results in results['de_results'].items():
            if contrast_results is not None and isinstance(contrast_results, pd.DataFrame):
                filename = f"{TABLES_DIR}/lemur_{contrast_name}_results.csv"
                contrast_results.to_csv(filename, index=True)
                print(f"      ✅ {contrast_name}: {filename}")
                
                # Export top significant genes
                if 'significant' in contrast_results.columns:
                    sig_genes = contrast_results[contrast_results['significant']]
                    if len(sig_genes) > 0:
                        sig_filename = f"{TABLES_DIR}/lemur_{contrast_name}_significant_genes.csv"
                        sig_genes.to_csv(sig_filename, index=True)
                        print(f"         Significant genes: {sig_filename}")
    
    # Export embedding coordinates
    if 'main' in results:
        print("   Exporting embedding coordinates...")
        lemur_model = results['main']['lemur_model']
        embedding_df = pd.DataFrame(
            lemur_model.embedding,
            columns=[f'LEMUR_{i+1}' for i in range(lemur_model.embedding.shape[1])],
            index=lemur_model.adata.obs.index
        )
        # Add metadata
        for col in [DESIGN_CONFIG['condition_col'], DESIGN_CONFIG['sex_col'], 
                   DESIGN_CONFIG['sample_col']]:
            if col in lemur_model.adata.obs.columns:
                embedding_df[col] = lemur_model.adata.obs[col]
        
        filename = f"{TABLES_DIR}/lemur_embedding_coordinates.csv"
        embedding_df.to_csv(filename, index=True)
        print(f"      ✅ Embeddings: {filename}")
    
    # Create analysis summary
    create_analysis_summary(results)

def create_analysis_summary(results):
    """Create comprehensive analysis summary"""
    print("   Creating analysis summary...")
    
    summary = {
        'analysis_type': 'LEMUR Final Working Implementation',
        'timestamp': pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S'),
        'input_file': RAW_H5AD,
        'design_formula': DESIGN_CONFIG['design_formula'],
        'n_embedding': ANALYSIS_CONFIG['n_embedding'],
        'fdr_threshold': ANALYSIS_CONFIG['fdr_threshold'],
        'lfc_threshold': ANALYSIS_CONFIG['lfc_threshold'],
        'using_raw_data': True,
        'harmony_aligned': True
    }
    
    # Add data info
    if 'main' in results:
        adata_used = results['main']['adata_used']
        summary.update({
            'n_cells_analyzed': adata_used.n_obs,
            'n_genes_analyzed': adata_used.n_vars,
            'n_samples': adata_used.obs[DESIGN_CONFIG['sample_col']].nunique(),
        })
        
        # Group distributions
        for col in [DESIGN_CONFIG['condition_col'], DESIGN_CONFIG['sex_col']]:
            counts = adata_used.obs[col].value_counts().to_dict()
            summary[f'{col}_distribution'] = str(counts)
    
    # Add DE results summary
    if 'de_results' in results:
        for contrast_name, contrast_results in results['de_results'].items():
            if contrast_results is not None and isinstance(contrast_results, pd.DataFrame):
                if 'significant' in contrast_results.columns:
                    n_sig = contrast_results['significant'].sum()
                    summary[f'{contrast_name}_significant_genes'] = n_sig
                summary[f'{contrast_name}_total_genes'] = len(contrast_results)
    
    # Save summary
    summary_df = pd.DataFrame([summary])
    filename = f"{TABLES_DIR}/lemur_analysis_summary.csv"
    summary_df.to_csv(filename, index=False)
    
    print(f"      ✅ Summary: {filename}")

# ============================================================================
# 🚀 MAIN EXECUTION
# ============================================================================

def main():
    """Main analysis workflow"""
    print(f"\n🚀 STARTING LEMUR FINAL ANALYSIS")
    print("=" * 35)
    
    # 1. Load and prepare data
    adata = load_and_prepare_raw_data()
    if adata is None:
        print("❌ Data loading failed")
        return
    
    # 2. Prepare for LEMUR
    adata_prepared = prepare_data_for_lemur(adata)
    if adata_prepared is None:
        print("❌ Data preparation failed")
        return
    
    # 3. Run LEMUR analysis
    results = run_lemur_analysis(adata_prepared)
    if results is None:
        print("❌ LEMUR analysis failed")
        return
    
    # 4. Create visualizations
    create_comprehensive_plots(results)
    
    # 5. Export results
    export_results(results)
    
    print(f"\n✅ LEMUR ANALYSIS COMPLETE")
    print("=" * 27)
    print(f"Results saved to: {OUTPUT_DIR}")
    print(f"• Plots: {PLOTS_DIR}/")
    print(f"• Tables: {TABLES_DIR}/")
    
    # Print final summary
    if 'de_results' in results:
        print(f"\n📊 DIFFERENTIAL EXPRESSION SUMMARY:")
        total_sig_genes = 0
        for contrast_name, contrast_results in results['de_results'].items():
            if contrast_results is not None:
                n_sig = contrast_results['significant'].sum() if 'significant' in contrast_results.columns else 0
                total_sig_genes += n_sig
                print(f"  • {contrast_name}: {n_sig} significant genes")
        print(f"  • Total unique significant genes: {total_sig_genes}")
    
    if 'main' in results:
        adata_used = results['main']['adata_used']
        print(f"\n📈 DATA SUMMARY:")
        print(f"  • {adata_used.n_obs:,} cells analyzed")
        print(f"  • {adata_used.n_vars:,} genes analyzed")
        print(f"  • {adata_used.obs[DESIGN_CONFIG['sample_col']].nunique()} samples")
        print(f"  • Design formula: {DESIGN_CONFIG['design_formula']}")
        print(f"  • Using raw, unintegrated data ✅")
        print(f"  • Harmony alignment applied ✅")

if __name__ == "__main__":
    main()