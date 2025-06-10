#!/usr/bin/env python3
"""
🌊 LEMUR Analysis - Sex-Stratified Implementation
GSE225158 - OUD vs Control with Sex-Specific Effects using pyLEMUR

This script implements comprehensive LEMUR analysis with sex-stratified contrasts:
1. OUD vs Control in Males only
2. OUD vs Control in Females only  
3. OUD effect differences between sexes
4. Interaction modeling (OUD × Sex)

Based on DESeq2 findings showing massive sex differences:
- Females: 2,743 significant genes
- Males: 51 significant genes
- Sex interaction: 3,130 significant genes

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
OUTPUT_DIR = f"{BASE_DIR}/results/snrna_scvi/lemur_sex_stratified"
PLOTS_DIR = f"{OUTPUT_DIR}/plots"
TABLES_DIR = f"{OUTPUT_DIR}/tables"

# Analysis parameters
ANALYSIS_CONFIG = {
    'n_embedding': 15,          # Number of latent dimensions
    'max_cells': 30000,         # Memory management (increased for sex stratification)
    'n_hvg': 3000,              # Highly variable genes
    'random_seed': 42,          # Reproducibility
    
    # Statistical thresholds (relaxed for discovery)
    'fdr_threshold': 0.1,       # Relaxed FDR threshold
    'lfc_threshold': 0.05,      # Relaxed log fold change threshold
    
    # Plotting
    'figure_dpi': 300,
    'point_size': 2,
    'alpha': 0.7,
}

# Experimental design for sex-stratified analysis
DESIGN_CONFIG = {
    'condition_col': 'level1',      # OUD vs CTL
    'sex_col': 'Sex',               # F vs M
    'sample_col': 'orig.ident',     # Sample ID for random effects
    'celltype_col': 'celltype3',    # Cell type annotation
    
    # Multiple design formulas for different analyses
    'designs': {
        'interaction': '~ level1 + Sex + level1:Sex',    # Full interaction model
        'simple': '~ level1 + Sex',                      # Simple additive model
        'males_only': '~ level1',                        # Males only
        'females_only': '~ level1'                       # Females only
    },
    
    # Sex-stratified contrasts matching DESeq2 analysis
    'contrasts': {
        'pooled_main_effect': {
            'name': 'OUD vs CTL (Pooled)',
            'design': 'simple',
            'subset': 'all',
            'coef_idx': 1,
            'description': 'Main OUD effect averaged across sexes'
        },
        'oud_vs_control_males': {
            'name': 'OUD vs CTL (Males Only)',
            'design': 'males_only', 
            'subset': 'males',
            'coef_idx': 1,
            'description': 'OUD effect specifically in male subjects'
        },
        'oud_vs_control_females': {
            'name': 'OUD vs CTL (Females Only)',
            'design': 'females_only',
            'subset': 'females', 
            'coef_idx': 1,
            'description': 'OUD effect specifically in female subjects'
        },
        'sex_main_effect': {
            'name': 'Male vs Female (Pooled)',
            'design': 'simple',
            'subset': 'all',
            'coef_idx': 2,
            'description': 'Main sex effect averaged across conditions'
        },
        'interaction_effect': {
            'name': 'OUD × Sex Interaction',
            'design': 'interaction',
            'subset': 'all',
            'coef_idx': 3,
            'description': 'How OUD effect differs between sexes'
        }
    }
}

# Create output directories
for directory in [OUTPUT_DIR, PLOTS_DIR, TABLES_DIR]:
    os.makedirs(directory, exist_ok=True)

print("🌊 LEMUR SEX-STRATIFIED ANALYSIS")
print("================================")
print(f"Input: {RAW_H5AD}")
print(f"Output: {OUTPUT_DIR}")
print(f"Sex-stratified contrasts: {len(DESIGN_CONFIG['contrasts'])}")

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
    
    # Print detailed group distributions
    print(f"\n📊 DETAILED GROUP DISTRIBUTIONS")
    print("=" * 32)
    
    condition_counts = adata.obs[DESIGN_CONFIG['condition_col']].value_counts()
    print(f"   {DESIGN_CONFIG['condition_col']}: {dict(condition_counts)}")
    
    sex_counts = adata.obs[DESIGN_CONFIG['sex_col']].value_counts()
    print(f"   {DESIGN_CONFIG['sex_col']}: {dict(sex_counts)}")
    
    sample_counts = adata.obs[DESIGN_CONFIG['sample_col']].nunique()
    print(f"   {DESIGN_CONFIG['sample_col']}: {sample_counts} unique samples")
    
    # Detailed cross-tabulation
    crosstab = pd.crosstab(adata.obs[DESIGN_CONFIG['condition_col']], 
                          adata.obs[DESIGN_CONFIG['sex_col']])
    print(f"\n   Cross-tabulation:")
    print(crosstab)
    
    # Sample distribution by sex and condition
    sample_dist = adata.obs.groupby([DESIGN_CONFIG['condition_col'], 
                                    DESIGN_CONFIG['sex_col']])[DESIGN_CONFIG['sample_col']].nunique()
    print(f"\n   Sample distribution:")
    print(sample_dist)
    
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
    
    min_group_size = group_sizes.min()
    if min_group_size < 1000:
        print(f"   ⚠️  Small group detected (n={min_group_size})")
    
    # Memory management - subsample if needed
    if adata.n_obs > ANALYSIS_CONFIG['max_cells']:
        print(f"   🔢 Subsampling to {ANALYSIS_CONFIG['max_cells']} cells...")
        adata = subsample_balanced(adata, ANALYSIS_CONFIG['max_cells'])
        print(f"   ✅ After subsampling: {adata.n_obs:,} cells")
    
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

def create_sex_stratified_datasets(adata):
    """Create sex-specific datasets for stratified analysis"""
    print(f"\n👥 CREATING SEX-STRATIFIED DATASETS")
    print("=" * 37)
    
    datasets = {}
    
    # Full dataset
    datasets['all'] = adata.copy()
    print(f"   All subjects: {adata.n_obs:,} cells")
    
    # Males only
    male_mask = adata.obs[DESIGN_CONFIG['sex_col']] == 'M'
    datasets['males'] = adata[male_mask].copy()
    male_dist = datasets['males'].obs[DESIGN_CONFIG['condition_col']].value_counts()
    print(f"   Males only: {datasets['males'].n_obs:,} cells | {dict(male_dist)}")
    
    # Females only  
    female_mask = adata.obs[DESIGN_CONFIG['sex_col']] == 'F'
    datasets['females'] = adata[female_mask].copy()
    female_dist = datasets['females'].obs[DESIGN_CONFIG['condition_col']].value_counts()
    print(f"   Females only: {datasets['females'].n_obs:,} cells | {dict(female_dist)}")
    
    return datasets

# ============================================================================
# 🌊 LEMUR ANALYSIS FUNCTIONS
# ============================================================================

def run_sex_stratified_lemur_analysis(adata):
    """Run comprehensive sex-stratified LEMUR analysis"""
    print(f"\n🌊 RUNNING SEX-STRATIFIED LEMUR ANALYSIS")
    print("=" * 42)
    
    # Create sex-stratified datasets
    datasets = create_sex_stratified_datasets(adata)
    
    results = {}
    
    # Run each contrast
    for contrast_name, contrast_config in DESIGN_CONFIG['contrasts'].items():
        print(f"\n🎯 CONTRAST: {contrast_name}")
        print("=" * (12 + len(contrast_name)))
        print(f"   Description: {contrast_config['description']}")
        
        try:
            # Get appropriate dataset
            subset_name = contrast_config['subset']
            if subset_name not in datasets:
                print(f"   ❌ Dataset subset '{subset_name}' not found")
                continue
                
            dataset = datasets[subset_name]
            design_formula = DESIGN_CONFIG['designs'][contrast_config['design']]
            
            print(f"   Dataset: {dataset.n_obs:,} cells")
            print(f"   Design: {design_formula}")
            
            # Fit LEMUR model
            lemur_result = fit_lemur_model_for_contrast(
                dataset, design_formula, contrast_name
            )
            
            if lemur_result is not None:
                # Extract differential expression
                de_result = extract_differential_expression_for_contrast(
                    lemur_result['lemur_model'], 
                    contrast_name,
                    contrast_config
                )
                
                if de_result is not None:
                    results[contrast_name] = {
                        'lemur_model': lemur_result['lemur_model'],
                        'de_results': de_result,
                        'dataset_used': dataset,
                        'contrast_config': contrast_config
                    }
                    
                    n_sig = de_result['significant'].sum() if 'significant' in de_result.columns else 0
                    print(f"   ✅ {contrast_name}: {n_sig} significant genes")
                else:
                    print(f"   ❌ DE extraction failed for {contrast_name}")
            else:
                print(f"   ❌ LEMUR fitting failed for {contrast_name}")
                
        except Exception as e:
            print(f"   ❌ Error in {contrast_name}: {e}")
    
    return results

def fit_lemur_model_for_contrast(adata, design_formula, contrast_name):
    """Fit LEMUR model for specific contrast"""
    print(f"      Fitting LEMUR model...")
    
    try:
        # Create LEMUR model
        lemur_model = pylemur.tl.LEMUR(
            adata,
            design=design_formula,
            n_embedding=ANALYSIS_CONFIG['n_embedding'],
            copy=True
        )
        
        print(f"         Model created with design: {design_formula}")
        
        # Fit the model
        print(f"         Fitting model...")
        lemur_model.fit(verbose=False)
        print(f"         ✅ Model fitted")
        
        # Try alignment
        print(f"         Attempting harmony alignment...")
        try:
            lemur_model.align_with_harmony()
            print(f"         ✅ Harmony alignment successful")
        except Exception as e:
            print(f"         ⚠️  Harmony alignment failed: {e}")
            print(f"         Continuing without alignment...")
        
        return {
            'lemur_model': lemur_model,
            'adata_used': adata
        }
        
    except Exception as e:
        print(f"         ❌ Model fitting failed: {e}")
        return None

def extract_differential_expression_for_contrast(lemur_model, contrast_name, contrast_config):
    """Extract differential expression for specific contrast"""
    print(f"         Extracting DE for {contrast_name}...")
    
    try:
        # Get coefficients - shape is (n_embedding, n_genes, n_coefficients)
        coefficients = lemur_model.coefficients
        print(f"            Coefficient tensor shape: {coefficients.shape}")
        
        # Extract coefficient index for this contrast
        coef_idx = contrast_config['coef_idx']
        
        if coef_idx >= coefficients.shape[2]:
            print(f"            ❌ Coefficient index {coef_idx} out of range (max: {coefficients.shape[2]-1})")
            return None
        
        # Extract coefficients for this contrast across all embedding dimensions
        contrast_coefficients = coefficients[:, :, coef_idx]
        print(f"            Extracted coefficients shape: {contrast_coefficients.shape}")
        
        # Average across embedding dimensions to get gene-level effects
        gene_coefficients = np.mean(contrast_coefficients, axis=0)
        print(f"            Gene-level coefficients shape: {gene_coefficients.shape}")
        
        # Get gene names
        gene_names = lemur_model.adata.var.index
        
        # Create results DataFrame
        results_df = pd.DataFrame({
            'gene': gene_names,
            'coefficient': gene_coefficients,
            'abs_coefficient': np.abs(gene_coefficients)
        })
        
        # Calculate statistics
        coef_std = np.std(gene_coefficients)
        results_df['z_score'] = np.abs(gene_coefficients) / (coef_std + 1e-8)
        results_df['pval'] = 2 * (1 - stats.norm.cdf(results_df['z_score']))
        
        # Multiple testing correction
        _, fdr_values, _, _ = multipletests(results_df['pval'], method='fdr_bh')
        results_df['fdr'] = fdr_values
        
        # Significance (relaxed thresholds for discovery)
        results_df['significant'] = (
            (results_df['fdr'] < ANALYSIS_CONFIG['fdr_threshold']) &
            (results_df['abs_coefficient'] > ANALYSIS_CONFIG['lfc_threshold'])
        )
        
        # Sort by significance
        results_df = results_df.sort_values('abs_coefficient', ascending=False)
        results_df = results_df.set_index('gene')
        
        print(f"            ✅ Extracted {len(results_df)} genes")
        print(f"            Significant genes: {results_df['significant'].sum()}")
        
        # Show top genes
        top_genes = results_df.head(5)
        print(f"            Top 5 genes by effect size:")
        for gene, row in top_genes.iterrows():
            sig_flag = "***" if row['significant'] else ""
            print(f"              {gene}: coef={row['coefficient']:.3f}, FDR={row['fdr']:.3e} {sig_flag}")
        
        return results_df
        
    except Exception as e:
        print(f"            ❌ DE extraction failed: {e}")
        return None

# ============================================================================
# 📊 VISUALIZATION FUNCTIONS
# ============================================================================

def create_sex_stratified_plots(results):
    """Create comprehensive sex-stratified visualization"""
    print(f"\n📊 CREATING SEX-STRATIFIED PLOTS")
    print("=" * 34)
    
    if not results:
        print("   ❌ No results to plot")
        return
    
    # 1. Comparison plots across contrasts
    create_contrast_comparison_plots(results)
    
    # 2. Sex-specific effect plots
    create_sex_effect_plots(results)
    
    # 3. Summary statistics
    create_stratified_summary_plots(results)

def create_contrast_comparison_plots(results):
    """Create comparison plots across different contrasts"""
    print("   Creating contrast comparison plots...")
    
    try:
        # Prepare data for comparison
        contrast_data = {}
        for contrast_name, result in results.items():
            if 'de_results' in result:
                contrast_data[contrast_name] = result['de_results']
        
        if len(contrast_data) < 2:
            print("      ⚠️  Need at least 2 contrasts for comparison")
            return
        
        # Create figure
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        axes = axes.flatten()
        
        # Plot 1: Effect size distributions
        ax = axes[0]
        for i, (name, data) in enumerate(contrast_data.items()):
            ax.hist(data['coefficient'], bins=50, alpha=0.6, 
                   label=name.replace('_', ' ').title(), density=True)
        ax.set_xlabel('Effect Size (Coefficient)')
        ax.set_ylabel('Density')
        ax.set_title('Effect Size Distributions by Contrast')
        ax.legend()
        
        # Plot 2: Significance counts
        ax = axes[1]
        contrast_names = []
        sig_counts = []
        for name, data in contrast_data.items():
            contrast_names.append(name.replace('_', '\n'))
            sig_counts.append(data['significant'].sum())
        
        bars = ax.bar(contrast_names, sig_counts, color=['blue', 'red', 'green', 'orange', 'purple'][:len(contrast_names)])
        ax.set_ylabel('Number of Significant Genes')
        ax.set_title('Significant Genes by Contrast')
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
        
        # Add counts on bars
        for bar, count in zip(bars, sig_counts):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + max(sig_counts)*0.01,
                   f'{count}', ha='center', va='bottom')
        
        # Plot 3: Males vs Females effect size comparison (if both available)
        if 'oud_vs_control_males' in contrast_data and 'oud_vs_control_females' in contrast_data:
            ax = axes[2]
            male_data = contrast_data['oud_vs_control_males']
            female_data = contrast_data['oud_vs_control_females']
            
            # Find common genes
            common_genes = set(male_data.index) & set(female_data.index)
            if len(common_genes) > 10:
                common_genes_list = list(common_genes)
                male_coefs = male_data.loc[common_genes_list, 'coefficient']
                female_coefs = female_data.loc[common_genes_list, 'coefficient']
                
                ax.scatter(male_coefs, female_coefs, alpha=0.6, s=20)
                ax.set_xlabel('Males: OUD vs CTL Coefficient')
                ax.set_ylabel('Females: OUD vs CTL Coefficient')
                ax.set_title('OUD Effect: Males vs Females')
                
                # Add diagonal line
                lims = [
                    np.min([ax.get_xlim(), ax.get_ylim()]),
                    np.max([ax.get_xlim(), ax.get_ylim()]),
                ]
                ax.plot(lims, lims, 'k-', alpha=0.5, zorder=0)
                
                # Calculate correlation
                corr, p_val = stats.pearsonr(male_coefs, female_coefs)
                ax.text(0.05, 0.95, f'r = {corr:.3f}\np = {p_val:.2e}', 
                       transform=ax.transAxes, verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Plot 4: Volcano plots for key contrasts
        plot_idx = 3
        for contrast_name in ['oud_vs_control_females', 'oud_vs_control_males'][:2]:
            if contrast_name in contrast_data and plot_idx < 6:
                ax = axes[plot_idx]
                data = contrast_data[contrast_name]
                
                x = data['coefficient']
                y = -np.log10(data['fdr'].clip(lower=1e-300))
                
                # Color by significance
                colors = ['red' if sig else 'gray' for sig in data['significant']]
                
                ax.scatter(x, y, c=colors, s=15, alpha=0.7)
                ax.set_xlabel('Coefficient')
                ax.set_ylabel('-log10(FDR)')
                ax.set_title(f'{contrast_name.replace("_", " ").title()}\nVolcano Plot')
                
                # Add significance thresholds
                ax.axhline(y=-np.log10(ANALYSIS_CONFIG['fdr_threshold']), 
                          color='blue', linestyle='--', alpha=0.5)
                ax.axvline(x=ANALYSIS_CONFIG['lfc_threshold'], 
                          color='blue', linestyle='--', alpha=0.5)
                ax.axvline(x=-ANALYSIS_CONFIG['lfc_threshold'], 
                          color='blue', linestyle='--', alpha=0.5)
                
                # Add text with number of significant genes
                n_sig = data['significant'].sum()
                ax.text(0.05, 0.95, f'{n_sig} significant genes', 
                       transform=ax.transAxes, verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                
                plot_idx += 1
        
        # Remove empty subplots
        for i in range(plot_idx, len(axes)):
            axes[i].remove()
        
        plt.tight_layout()
        plt.savefig(f"{PLOTS_DIR}/sex_stratified_contrast_comparison.png", 
                   dpi=ANALYSIS_CONFIG['figure_dpi'], bbox_inches='tight')
        plt.close()
        
        print("      ✅ Contrast comparison plots saved")
        
    except Exception as e:
        print(f"      ❌ Contrast plotting failed: {e}")

def create_sex_effect_plots(results):
    """Create sex-specific effect visualization"""
    print("   Creating sex effect plots...")
    
    try:
        # Check if we have sex-stratified results
        if 'oud_vs_control_males' not in results or 'oud_vs_control_females' not in results:
            print("      ⚠️  Missing sex-stratified results for comparison")
            return
        
        male_results = results['oud_vs_control_males']['de_results']
        female_results = results['oud_vs_control_females']['de_results']
        
        # Find common genes
        common_genes = set(male_results.index) & set(female_results.index)
        if len(common_genes) == 0:
            print("      ❌ No common genes between male and female analyses")
            return
        
        common_genes_list = list(common_genes)
        
        # Create comparison DataFrame
        comparison_df = pd.DataFrame(index=common_genes_list)
        comparison_df['male_coef'] = male_results.loc[common_genes_list, 'coefficient']
        comparison_df['female_coef'] = female_results.loc[common_genes_list, 'coefficient']
        comparison_df['male_fdr'] = male_results.loc[common_genes_list, 'fdr']
        comparison_df['female_fdr'] = female_results.loc[common_genes_list, 'fdr']
        comparison_df['male_sig'] = male_results.loc[common_genes_list, 'significant']
        comparison_df['female_sig'] = female_results.loc[common_genes_list, 'significant']
        
        # Calculate sex difference
        comparison_df['sex_difference'] = comparison_df['female_coef'] - comparison_df['male_coef']
        comparison_df['abs_sex_difference'] = np.abs(comparison_df['sex_difference'])
        
        # Create figure
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Plot 1: Male vs Female effect sizes
        ax = axes[0, 0]
        scatter = ax.scatter(comparison_df['male_coef'], comparison_df['female_coef'], 
                           alpha=0.6, s=20, c='blue')
        ax.set_xlabel('Male OUD Effect (Coefficient)')
        ax.set_ylabel('Female OUD Effect (Coefficient)')
        ax.set_title('OUD Effect: Male vs Female')
        
        # Add diagonal line
        lims = [
            np.min([ax.get_xlim(), ax.get_ylim()]),
            np.max([ax.get_xlim(), ax.get_ylim()]),
        ]
        ax.plot(lims, lims, 'k--', alpha=0.5)
        
        # Add correlation
        corr, p_val = stats.pearsonr(comparison_df['male_coef'], comparison_df['female_coef'])
        ax.text(0.05, 0.95, f'r = {corr:.3f}\np = {p_val:.2e}', 
               transform=ax.transAxes, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Plot 2: Sex difference distribution
        ax = axes[0, 1]
        ax.hist(comparison_df['sex_difference'], bins=50, alpha=0.7, color='purple')
        ax.set_xlabel('Sex Difference (Female - Male)')
        ax.set_ylabel('Frequency')
        ax.set_title('Distribution of Sex Differences')
        ax.axvline(x=0, color='red', linestyle='--', alpha=0.7)
        
        # Plot 3: Top sex-differential genes
        ax = axes[1, 0]
        top_diff_genes = comparison_df.nlargest(15, 'abs_sex_difference')
        
        y_pos = np.arange(len(top_diff_genes))
        bars = ax.barh(y_pos, top_diff_genes['sex_difference'], 
                      color=['red' if x > 0 else 'blue' for x in top_diff_genes['sex_difference']])
        ax.set_yticks(y_pos)
        ax.set_yticklabels(top_diff_genes.index, fontsize=8)
        ax.set_xlabel('Sex Difference (Female - Male)')
        ax.set_title('Top Sex-Differential Genes')
        ax.axvline(x=0, color='black', linestyle='-', alpha=0.5)
        
        # Plot 4: Significance overlap
        ax = axes[1, 1]
        sig_categories = ['Neither', 'Male Only', 'Female Only', 'Both']
        sig_counts = [
            (~comparison_df['male_sig'] & ~comparison_df['female_sig']).sum(),
            (comparison_df['male_sig'] & ~comparison_df['female_sig']).sum(),
            (~comparison_df['male_sig'] & comparison_df['female_sig']).sum(),
            (comparison_df['male_sig'] & comparison_df['female_sig']).sum()
        ]
        
        colors = ['gray', 'lightblue', 'lightcoral', 'purple']
        wedges, texts, autotexts = ax.pie(sig_counts, labels=sig_categories, autopct='%1.1f%%', 
                                         colors=colors, startangle=90)
        ax.set_title('Significance Overlap\nMale vs Female')
        
        plt.tight_layout()
        plt.savefig(f"{PLOTS_DIR}/sex_effect_comparison.png", 
                   dpi=ANALYSIS_CONFIG['figure_dpi'], bbox_inches='tight')
        plt.close()
        
        # Save sex difference results
        comparison_df.to_csv(f"{TABLES_DIR}/sex_difference_analysis.csv")
        
        print("      ✅ Sex effect plots saved")
        
    except Exception as e:
        print(f"      ❌ Sex effect plotting failed: {e}")

def create_stratified_summary_plots(results):
    """Create summary statistics plots for stratified analysis"""
    print("   Creating stratified summary plots...")
    
    try:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Summary statistics
        summary_stats = []
        for contrast_name, result in results.items():
            if 'de_results' in result:
                de_results = result['de_results']
                if 'significant' in de_results.columns:
                    n_sig = de_results['significant'].sum()
                    n_total = len(de_results)
                    summary_stats.append({
                        'contrast': contrast_name.replace('_', ' ').title(),
                        'significant_genes': n_sig,
                        'total_genes': n_total,
                        'pct_significant': n_sig/n_total*100 if n_total > 0 else 0
                    })
        
        if summary_stats:
            summary_df = pd.DataFrame(summary_stats)
            
            # Bar plot of significant genes
            bars = axes[0, 0].bar(range(len(summary_df)), summary_df['significant_genes'])
            axes[0, 0].set_title('Significant Genes by Contrast')
            axes[0, 0].set_ylabel('Number of Significant Genes')
            axes[0, 0].set_xticks(range(len(summary_df)))
            axes[0, 0].set_xticklabels(summary_df['contrast'], rotation=45, ha='right')
            
            # Add counts on bars
            for bar, count in zip(bars, summary_df['significant_genes']):
                height = bar.get_height()
                axes[0, 0].text(bar.get_x() + bar.get_width()/2., height + max(summary_df['significant_genes'])*0.01,
                               f'{count}', ha='center', va='bottom')
            
            # Percentage plot
            axes[0, 1].bar(range(len(summary_df)), summary_df['pct_significant'])
            axes[0, 1].set_title('Percentage of Significant Genes')
            axes[0, 1].set_ylabel('Percentage Significant')
            axes[0, 1].set_xticks(range(len(summary_df)))
            axes[0, 1].set_xticklabels(summary_df['contrast'], rotation=45, ha='right')
        
        # Analysis parameters
        axes[1, 0].text(0.1, 0.9, f"Analysis Parameters:", fontweight='bold', transform=axes[1, 0].transAxes)
        axes[1, 0].text(0.1, 0.8, f"• Embedding dims: {ANALYSIS_CONFIG['n_embedding']}", transform=axes[1, 0].transAxes)
        axes[1, 0].text(0.1, 0.7, f"• HVG: {ANALYSIS_CONFIG['n_hvg']}", transform=axes[1, 0].transAxes)
        axes[1, 0].text(0.1, 0.6, f"• FDR threshold: {ANALYSIS_CONFIG['fdr_threshold']}", transform=axes[1, 0].transAxes)
        axes[1, 0].text(0.1, 0.5, f"• LFC threshold: {ANALYSIS_CONFIG['lfc_threshold']}", transform=axes[1, 0].transAxes)
        axes[1, 0].text(0.1, 0.4, f"• Max cells: {ANALYSIS_CONFIG['max_cells']}", transform=axes[1, 0].transAxes)
        axes[1, 0].set_xlim(0, 1)
        axes[1, 0].set_ylim(0, 1)
        axes[1, 0].axis('off')
        axes[1, 0].set_title('Analysis Parameters')
        
        # Contrasts tested
        axes[1, 1].text(0.1, 0.9, f"Contrasts Analyzed:", fontweight='bold', transform=axes[1, 1].transAxes)
        y_pos = 0.8
        for contrast_name in DESIGN_CONFIG['contrasts'].keys():
            axes[1, 1].text(0.1, y_pos, f"• {contrast_name.replace('_', ' ').title()}", transform=axes[1, 1].transAxes)
            y_pos -= 0.1
        axes[1, 1].set_xlim(0, 1)
        axes[1, 1].set_ylim(0, 1)
        axes[1, 1].axis('off')
        axes[1, 1].set_title('Contrasts Analyzed')
        
        plt.tight_layout()
        plt.savefig(f"{PLOTS_DIR}/stratified_analysis_summary.png", 
                   dpi=ANALYSIS_CONFIG['figure_dpi'], bbox_inches='tight')
        plt.close()
        
        print("      ✅ Stratified summary plots saved")
        
    except Exception as e:
        print(f"      ❌ Stratified summary plotting failed: {e}")

# ============================================================================
# 💾 EXPORT FUNCTIONS
# ============================================================================

def export_sex_stratified_results(results):
    """Export all sex-stratified results"""
    print(f"\n💾 EXPORTING SEX-STRATIFIED RESULTS")
    print("=" * 37)
    
    # Export each contrast's results
    for contrast_name, result in results.items():
        if 'de_results' in result:
            print(f"   Exporting {contrast_name}...")
            de_results = result['de_results']
            
            # Export full results
            filename = f"{TABLES_DIR}/lemur_{contrast_name}_results.csv"
            de_results.to_csv(filename, index=True)
            print(f"      ✅ Full results: {filename}")
            
            # Export significant genes only
            if 'significant' in de_results.columns:
                sig_genes = de_results[de_results['significant']]
                if len(sig_genes) > 0:
                    sig_filename = f"{TABLES_DIR}/lemur_{contrast_name}_significant_genes.csv"
                    sig_genes.to_csv(sig_filename, index=True)
                    print(f"      ✅ Significant genes ({len(sig_genes)}): {sig_filename}")
    
    # Create comprehensive summary
    create_comprehensive_summary(results)

def create_comprehensive_summary(results):
    """Create comprehensive analysis summary"""
    print("   Creating comprehensive summary...")
    
    summary_lines = [
        "# LEMUR Sex-Stratified Analysis Report",
        "=" * 40,
        "",
        "## Analysis Overview",
        f"This analysis implements sex-stratified differential expression using LEMUR",
        f"to understand how OUD affects males and females differently.",
        "",
        "## Contrasts Analyzed",
    ]
    
    for contrast_name, contrast_config in DESIGN_CONFIG['contrasts'].items():
        summary_lines.append(f"- **{contrast_name}**: {contrast_config['description']}")
    
    summary_lines.extend([
        "",
        "## Results Summary",
    ])
    
    # Add results for each contrast
    total_sig_genes = 0
    for contrast_name, result in results.items():
        if 'de_results' in result:
            de_results = result['de_results']
            n_total = len(de_results)
            n_sig = de_results['significant'].sum() if 'significant' in de_results.columns else 0
            total_sig_genes += n_sig
            
            summary_lines.extend([
                f"",
                f"### {contrast_name.replace('_', ' ').title()}",
                f"- Total genes: {n_total:,}",
                f"- Significant genes: {n_sig:,}",
                f"- Percentage significant: {n_sig/n_total*100:.2f}%",
            ])
            
            if n_sig > 0:
                top_genes = de_results.nlargest(5, 'abs_coefficient')
                summary_lines.append("- Top 5 genes by effect size:")
                for gene, row in top_genes.iterrows():
                    summary_lines.append(f"  - {gene}: coef={row['coefficient']:.3f}, FDR={row['fdr']:.2e}")
    
    summary_lines.extend([
        "",
        "## Key Findings",
        f"- Total significant genes across all contrasts: {total_sig_genes:,}",
        "",
        "## Comparison with DESeq2",
        "- DESeq2 Females: 2,743 significant genes (10.91%)",
        "- DESeq2 Males: 51 significant genes (0.21%)", 
        "- DESeq2 Sex interaction: 3,130 significant genes (12.59%)",
        "",
        "## Analysis Parameters",
        f"- Embedding dimensions: {ANALYSIS_CONFIG['n_embedding']}",
        f"- FDR threshold: {ANALYSIS_CONFIG['fdr_threshold']}",
        f"- Effect size threshold: {ANALYSIS_CONFIG['lfc_threshold']}",
        f"- Maximum cells: {ANALYSIS_CONFIG['max_cells']:,}",
        f"- Highly variable genes: {ANALYSIS_CONFIG['n_hvg']:,}",
    ])
    
    # Write summary
    summary_filename = f"{TABLES_DIR}/sex_stratified_analysis_summary.md"
    with open(summary_filename, 'w') as f:
        f.write('\n'.join(summary_lines))
    
    print(f"      ✅ Comprehensive summary: {summary_filename}")

# ============================================================================
# 🚀 MAIN EXECUTION
# ============================================================================

def main():
    """Main sex-stratified analysis workflow"""
    print(f"\n🚀 STARTING SEX-STRATIFIED LEMUR ANALYSIS")
    print("=" * 43)
    
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
    
    # 3. Run sex-stratified LEMUR analysis
    results = run_sex_stratified_lemur_analysis(adata_prepared)
    if not results:
        print("❌ Sex-stratified LEMUR analysis failed")
        return
    
    # 4. Create visualizations
    create_sex_stratified_plots(results)
    
    # 5. Export results
    export_sex_stratified_results(results)
    
    print(f"\n✅ SEX-STRATIFIED LEMUR ANALYSIS COMPLETE")
    print("=" * 42)
    print(f"Results saved to: {OUTPUT_DIR}")
    print(f"• Plots: {PLOTS_DIR}/")
    print(f"• Tables: {TABLES_DIR}/")
    
    # Print final summary
    print(f"\n📊 SEX-STRATIFIED RESULTS SUMMARY:")
    successful_contrasts = 0
    total_sig_genes = 0
    
    for contrast_name, result in results.items():
        if 'de_results' in result:
            successful_contrasts += 1
            de_results = result['de_results']
            n_sig = de_results['significant'].sum() if 'significant' in de_results.columns else 0
            total_sig_genes += n_sig
            print(f"  • {contrast_name}: {n_sig} significant genes")
    
    print(f"  • Total contrasts analyzed: {successful_contrasts}")
    print(f"  • Total significant genes: {total_sig_genes}")
    
    if 'oud_vs_control_females' in results and 'oud_vs_control_males' in results:
        female_sig = results['oud_vs_control_females']['de_results']['significant'].sum()
        male_sig = results['oud_vs_control_males']['de_results']['significant'].sum()
        
        print(f"\n🎯 SEX-SPECIFIC OUD EFFECTS:")
        print(f"  • Females: {female_sig} significant genes")
        print(f"  • Males: {male_sig} significant genes")
        print(f"  • Female/Male ratio: {female_sig/male_sig:.1f}x" if male_sig > 0 else "  • Female/Male ratio: ∞ (no male effects)")

if __name__ == "__main__":
    main()