#!/usr/bin/env python3
"""
🌊 LEMUR Analysis - Regional Contrasts Implementation
GSE225158 - OUD vs Control with Regional Stratification (Caudate vs Putamen)

This script implements LEMUR analysis with regional contrasts to complete
the methodological comparison with DESeq2:

Core Regional Contrasts:
1. OUD vs Control in Caudate only
2. OUD vs Control in Putamen only
3. Regional differences in OUD effects
4. Regional × Sex interactions

Regional Biology:
- Caudate: Goal-directed behavior, cognitive control
- Putamen: Habit formation, motor control, addiction maintenance

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
OUTPUT_DIR = f"{BASE_DIR}/results/snrna_scvi/lemur_regional_contrasts"
PLOTS_DIR = f"{OUTPUT_DIR}/plots"
TABLES_DIR = f"{OUTPUT_DIR}/tables"

# Analysis parameters
ANALYSIS_CONFIG = {
    'n_embedding': 15,          # Number of latent dimensions
    'max_cells': 35000,         # Increased for regional analysis
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

# Regional experimental design
DESIGN_CONFIG = {
    'condition_col': 'level1',      # OUD vs CTL
    'sex_col': 'Sex',               # F vs M
    'region_col': 'Region',         # Caudate vs Putamen
    'sample_col': 'orig.ident',     # Sample ID for random effects
    'celltype_col': 'celltype3',    # Cell type annotation
    
    # Regional design formulas
    'designs': {
        'simple': '~ level1',                                # Simple OUD effect
        'with_sex': '~ level1 + Sex',                       # OUD + sex effects
        'interaction': '~ level1 + Sex + level1:Sex',       # Full interaction
        'regional': '~ level1 + Region',                    # OUD + region effects
        'full': '~ level1 + Sex + Region + level1:Sex'      # Complete model
    },
    
    # Regional contrasts matching DESeq2 analysis
    'contrasts': {
        'oud_vs_control_caudate': {
            'name': 'OUD vs CTL (Caudate Only)',
            'design': 'simple',
            'subset': 'caudate',
            'coef_idx': 1,
            'description': 'OUD effect specifically in Caudate region'
        },
        'oud_vs_control_putamen': {
            'name': 'OUD vs CTL (Putamen Only)',
            'design': 'simple',
            'subset': 'putamen',
            'coef_idx': 1,
            'description': 'OUD effect specifically in Putamen region'
        },
        'oud_vs_control_caudate_males': {
            'name': 'OUD vs CTL (Caudate Males)',
            'design': 'simple',
            'subset': 'caudate_males',
            'coef_idx': 1,
            'description': 'OUD effect in Caudate males only'
        },
        'oud_vs_control_caudate_females': {
            'name': 'OUD vs CTL (Caudate Females)',
            'design': 'simple',
            'subset': 'caudate_females',
            'coef_idx': 1,
            'description': 'OUD effect in Caudate females only'
        },
        'oud_vs_control_putamen_males': {
            'name': 'OUD vs CTL (Putamen Males)',
            'design': 'simple',
            'subset': 'putamen_males',
            'coef_idx': 1,
            'description': 'OUD effect in Putamen males only'
        },
        'oud_vs_control_putamen_females': {
            'name': 'OUD vs CTL (Putamen Females)',
            'design': 'simple',
            'subset': 'putamen_females',
            'coef_idx': 1,
            'description': 'OUD effect in Putamen females only'
        },
        'sex_effect_caudate': {
            'name': 'Male vs Female (Caudate)',
            'design': 'with_sex',
            'subset': 'caudate',
            'coef_idx': 2,
            'description': 'Sex effect in Caudate region'
        },
        'sex_effect_putamen': {
            'name': 'Male vs Female (Putamen)',
            'design': 'with_sex',
            'subset': 'putamen',
            'coef_idx': 2,
            'description': 'Sex effect in Putamen region'
        },
        'regional_main_effect': {
            'name': 'Caudate vs Putamen (Pooled)',
            'design': 'regional',
            'subset': 'all',
            'coef_idx': 2,
            'description': 'Main regional differences'
        }
    }
}

# Create output directories
for directory in [OUTPUT_DIR, PLOTS_DIR, TABLES_DIR]:
    os.makedirs(directory, exist_ok=True)

print("🌊 LEMUR REGIONAL CONTRASTS ANALYSIS")
print("====================================")
print(f"Input: {RAW_H5AD}")
print(f"Output: {OUTPUT_DIR}")
print(f"Regional contrasts: {len(DESIGN_CONFIG['contrasts'])}")

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
                    DESIGN_CONFIG['region_col'],
                    DESIGN_CONFIG['sample_col']]
    
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    if missing_cols:
        print(f"   ❌ Missing required columns: {missing_cols}")
        return None
    
    print(f"   ✅ All required metadata columns present")
    
    # Print detailed regional distributions
    print(f"\n📊 REGIONAL DISTRIBUTIONS")
    print("=" * 26)
    
    condition_counts = adata.obs[DESIGN_CONFIG['condition_col']].value_counts()
    print(f"   {DESIGN_CONFIG['condition_col']}: {dict(condition_counts)}")
    
    sex_counts = adata.obs[DESIGN_CONFIG['sex_col']].value_counts()
    print(f"   {DESIGN_CONFIG['sex_col']}: {dict(sex_counts)}")
    
    region_counts = adata.obs[DESIGN_CONFIG['region_col']].value_counts()
    print(f"   {DESIGN_CONFIG['region_col']}: {dict(region_counts)}")
    
    sample_counts = adata.obs[DESIGN_CONFIG['sample_col']].nunique()
    print(f"   {DESIGN_CONFIG['sample_col']}: {sample_counts} unique samples")
    
    # Regional cross-tabulations
    print(f"\n   Region × Condition:")
    region_condition = pd.crosstab(adata.obs[DESIGN_CONFIG['region_col']], 
                                  adata.obs[DESIGN_CONFIG['condition_col']])
    print(region_condition)
    
    print(f"\n   Region × Sex:")
    region_sex = pd.crosstab(adata.obs[DESIGN_CONFIG['region_col']], 
                            adata.obs[DESIGN_CONFIG['sex_col']])
    print(region_sex)
    
    # Three-way cross-tabulation
    print(f"\n   Region × Condition × Sex:")
    for region in adata.obs[DESIGN_CONFIG['region_col']].unique():
        region_data = adata.obs[adata.obs[DESIGN_CONFIG['region_col']] == region]
        region_crosstab = pd.crosstab(region_data[DESIGN_CONFIG['condition_col']], 
                                     region_data[DESIGN_CONFIG['sex_col']])
        print(f"   {region}:")
        print(region_crosstab)
        print()
    
    return adata

def prepare_data_for_lemur(adata):
    """Prepare data following LEMUR requirements"""
    print(f"\n🔧 PREPARING DATA FOR LEMUR")
    print("=" * 30)
    
    # Check for minimum group sizes
    print("   Checking regional group sizes...")
    group_sizes = adata.obs.groupby([DESIGN_CONFIG['region_col'],
                                    DESIGN_CONFIG['condition_col'], 
                                    DESIGN_CONFIG['sex_col']]).size()
    print("   Group sizes:")
    print(group_sizes)
    
    min_group_size = group_sizes.min()
    if min_group_size < 1000:
        print(f"   ⚠️  Small group detected (n={min_group_size})")
    
    # Memory management - subsample if needed
    if adata.n_obs > ANALYSIS_CONFIG['max_cells']:
        print(f"   🔢 Subsampling to {ANALYSIS_CONFIG['max_cells']} cells...")
        adata = subsample_balanced_regional(adata, ANALYSIS_CONFIG['max_cells'])
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
    for col in [DESIGN_CONFIG['condition_col'], DESIGN_CONFIG['sex_col'], DESIGN_CONFIG['region_col']]:
        if adata_hvg.obs[col].dtype.name == 'category':
            adata_hvg.obs[col] = adata_hvg.obs[col].astype(str)
        adata_hvg.obs[col] = pd.Categorical(adata_hvg.obs[col])
    
    print(f"   ✅ Data prepared: {adata_hvg.n_obs:,} cells × {adata_hvg.n_vars:,} genes")
    
    return adata_hvg

def subsample_balanced_regional(adata, max_cells):
    """Subsample data while maintaining regional and group balance"""
    np.random.seed(ANALYSIS_CONFIG['random_seed'])
    
    groups = adata.obs.groupby([DESIGN_CONFIG['region_col'],
                               DESIGN_CONFIG['condition_col'], 
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

def create_regional_datasets(adata):
    """Create region-specific datasets for stratified analysis"""
    print(f"\n🏗️  CREATING REGIONAL DATASETS")
    print("=" * 32)
    
    datasets = {}
    
    # Full dataset
    datasets['all'] = adata.copy()
    print(f"   All regions: {adata.n_obs:,} cells")
    
    # Regional datasets
    for region in adata.obs[DESIGN_CONFIG['region_col']].unique():
        region_mask = adata.obs[DESIGN_CONFIG['region_col']] == region
        datasets[region.lower()] = adata[region_mask].copy()
        
        region_data = datasets[region.lower()]
        condition_dist = region_data.obs[DESIGN_CONFIG['condition_col']].value_counts()
        sex_dist = region_data.obs[DESIGN_CONFIG['sex_col']].value_counts()
        
        print(f"   {region}: {region_data.n_obs:,} cells")
        print(f"      Condition: {dict(condition_dist)}")
        print(f"      Sex: {dict(sex_dist)}")
    
    # Region × Sex combinations
    for region in adata.obs[DESIGN_CONFIG['region_col']].unique():
        for sex in adata.obs[DESIGN_CONFIG['sex_col']].unique():
            subset_name = f"{region.lower()}_{sex.lower()}s"
            
            mask = (adata.obs[DESIGN_CONFIG['region_col']] == region) & \
                   (adata.obs[DESIGN_CONFIG['sex_col']] == sex)
            datasets[subset_name] = adata[mask].copy()
            
            subset_data = datasets[subset_name]
            if subset_data.n_obs > 0:
                condition_dist = subset_data.obs[DESIGN_CONFIG['condition_col']].value_counts()
                print(f"   {region} {sex}s: {subset_data.n_obs:,} cells | {dict(condition_dist)}")
    
    return datasets

# ============================================================================
# 🌊 LEMUR ANALYSIS FUNCTIONS
# ============================================================================

def run_regional_lemur_analysis(adata):
    """Run comprehensive regional LEMUR analysis"""
    print(f"\n🌊 RUNNING REGIONAL LEMUR ANALYSIS")
    print("=" * 36)
    
    # Create regional datasets
    datasets = create_regional_datasets(adata)
    
    results = {}
    
    # Run each regional contrast
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
            
            # Check minimum cell count
            if dataset.n_obs < 1000:
                print(f"   ⚠️  Warning: Low cell count ({dataset.n_obs})")
            
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
            import traceback
            traceback.print_exc()
    
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

def create_regional_plots(results):
    """Create comprehensive regional visualization"""
    print(f"\n📊 CREATING REGIONAL PLOTS")
    print("=" * 28)
    
    if not results:
        print("   ❌ No results to plot")
        return
    
    # 1. Regional comparison plots
    create_regional_comparison_plots(results)
    
    # 2. Caudate vs Putamen plots
    create_caudate_vs_putamen_plots(results)
    
    # 3. Regional summary statistics
    create_regional_summary_plots(results)

def create_regional_comparison_plots(results):
    """Create comparison plots across regional contrasts"""
    print("   Creating regional comparison plots...")
    
    try:
        # Filter results for main regional contrasts
        regional_results = {}
        for key in ['oud_vs_control_caudate', 'oud_vs_control_putamen']:
            if key in results and 'de_results' in results[key]:
                regional_results[key] = results[key]['de_results']
        
        if len(regional_results) < 2:
            print("      ⚠️  Need both Caudate and Putamen results for comparison")
            return
        
        # Create figure
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        axes = axes.flatten()
        
        # Plot 1: Effect size distributions
        ax = axes[0]
        for i, (name, data) in enumerate(regional_results.items()):
            region_name = name.replace('oud_vs_control_', '').title()
            ax.hist(data['coefficient'], bins=50, alpha=0.6, 
                   label=region_name, density=True)
        ax.set_xlabel('Effect Size (Coefficient)')
        ax.set_ylabel('Density')
        ax.set_title('OUD Effect Size by Region')
        ax.legend()
        
        # Plot 2: Significance counts
        ax = axes[1]
        region_names = []
        sig_counts = []
        for name, data in regional_results.items():
            region_names.append(name.replace('oud_vs_control_', '').title())
            sig_counts.append(data['significant'].sum())
        
        bars = ax.bar(region_names, sig_counts, color=['lightcoral', 'lightblue'])
        ax.set_ylabel('Number of Significant Genes')
        ax.set_title('Significant Genes by Region')
        
        # Add counts on bars
        for bar, count in zip(bars, sig_counts):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + max(sig_counts)*0.01,
                   f'{count}', ha='center', va='bottom')
        
        # Plot 3: Caudate vs Putamen scatter plot
        if 'oud_vs_control_caudate' in regional_results and 'oud_vs_control_putamen' in regional_results:
            ax = axes[2]
            caudate_data = regional_results['oud_vs_control_caudate']
            putamen_data = regional_results['oud_vs_control_putamen']
            
            # Find common genes
            common_genes = set(caudate_data.index) & set(putamen_data.index)
            if len(common_genes) > 10:
                common_genes_list = list(common_genes)
                caudate_coefs = caudate_data.loc[common_genes_list, 'coefficient']
                putamen_coefs = putamen_data.loc[common_genes_list, 'coefficient']
                
                ax.scatter(caudate_coefs, putamen_coefs, alpha=0.6, s=20)
                ax.set_xlabel('Caudate: OUD Effect Coefficient')
                ax.set_ylabel('Putamen: OUD Effect Coefficient')
                ax.set_title('Regional OUD Effects: Caudate vs Putamen')
                
                # Add diagonal line
                lims = [
                    np.min([ax.get_xlim(), ax.get_ylim()]),
                    np.max([ax.get_xlim(), ax.get_ylim()]),
                ]
                ax.plot(lims, lims, 'k-', alpha=0.5, zorder=0)
                
                # Calculate correlation
                corr, p_val = stats.pearsonr(caudate_coefs, putamen_coefs)
                ax.text(0.05, 0.95, f'r = {corr:.3f}\np = {p_val:.2e}', 
                       transform=ax.transAxes, verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Plot 4-5: Volcano plots
        plot_idx = 3
        for contrast_name in ['oud_vs_control_caudate', 'oud_vs_control_putamen']:
            if contrast_name in regional_results and plot_idx < 6:
                ax = axes[plot_idx]
                data = regional_results[contrast_name]
                
                x = data['coefficient']
                y = -np.log10(data['fdr'].clip(lower=1e-300))
                
                # Color by significance
                colors = ['red' if sig else 'gray' for sig in data['significant']]
                
                ax.scatter(x, y, c=colors, s=15, alpha=0.7)
                ax.set_xlabel('Coefficient')
                ax.set_ylabel('-log10(FDR)')
                region_name = contrast_name.replace('oud_vs_control_', '').title()
                ax.set_title(f'{region_name} OUD Effect\nVolcano Plot')
                
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
        
        # Remove empty subplot
        if plot_idx < len(axes):
            axes[plot_idx].remove()
        
        plt.tight_layout()
        plt.savefig(f"{PLOTS_DIR}/regional_comparison_plots.png", 
                   dpi=ANALYSIS_CONFIG['figure_dpi'], bbox_inches='tight')
        plt.close()
        
        print("      ✅ Regional comparison plots saved")
        
    except Exception as e:
        print(f"      ❌ Regional comparison plotting failed: {e}")

def create_caudate_vs_putamen_plots(results):
    """Create detailed Caudate vs Putamen comparison"""
    print("   Creating Caudate vs Putamen detailed plots...")
    
    try:
        # Collect all regional results
        regional_data = {}
        for contrast_name, result in results.items():
            if 'de_results' in result:
                regional_data[contrast_name] = result['de_results']
        
        if len(regional_data) == 0:
            return
        
        # Create comprehensive comparison figure
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Plot 1: Regional effect summary
        ax = axes[0, 0]
        contrast_names = []
        sig_counts = []
        
        for name, data in regional_data.items():
            if 'significant' in data.columns:
                display_name = name.replace('oud_vs_control_', '').replace('_', ' ').title()
                contrast_names.append(display_name)
                sig_counts.append(data['significant'].sum())
        
        if contrast_names:
            bars = ax.bar(range(len(contrast_names)), sig_counts, 
                         color=['lightcoral' if 'caudate' in name.lower() else 'lightblue' 
                               for name in contrast_names])
            ax.set_ylabel('Significant Genes')
            ax.set_title('Significant Genes by Regional Contrast')
            ax.set_xticks(range(len(contrast_names)))
            ax.set_xticklabels(contrast_names, rotation=45, ha='right')
            
            # Add counts
            for bar, count in zip(bars, sig_counts):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height + max(sig_counts)*0.01,
                       f'{count}', ha='center', va='bottom')
        
        # Plot 2: Sex-specific regional effects
        ax = axes[0, 1]
        sex_regional_results = {}
        for key in ['oud_vs_control_caudate_males', 'oud_vs_control_caudate_females',
                   'oud_vs_control_putamen_males', 'oud_vs_control_putamen_females']:
            if key in regional_data:
                sex_regional_results[key] = regional_data[key]['significant'].sum()
        
        if sex_regional_results:
            labels = [k.replace('oud_vs_control_', '').replace('_', ' ').title() 
                     for k in sex_regional_results.keys()]
            values = list(sex_regional_results.values())
            colors = ['lightcoral' if 'caudate' in label.lower() else 'lightblue' 
                     for label in labels]
            
            bars = ax.bar(range(len(labels)), values, color=colors)
            ax.set_ylabel('Significant Genes')
            ax.set_title('Sex-Specific Regional Effects')
            ax.set_xticks(range(len(labels)))
            ax.set_xticklabels(labels, rotation=45, ha='right')
        
        # Plot 3: Regional biology summary
        ax = axes[1, 0]
        ax.text(0.1, 0.9, "Regional Biology:", fontweight='bold', transform=ax.transAxes)
        ax.text(0.1, 0.8, "• Caudate: Goal-directed behavior", transform=ax.transAxes)
        ax.text(0.1, 0.7, "• Caudate: Cognitive control", transform=ax.transAxes)
        ax.text(0.1, 0.6, "• Putamen: Habit formation", transform=ax.transAxes)
        ax.text(0.1, 0.5, "• Putamen: Motor control", transform=ax.transAxes)
        ax.text(0.1, 0.4, "• Putamen: Addiction maintenance", transform=ax.transAxes)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        ax.set_title('Functional Roles')
        
        # Plot 4: Analysis summary
        ax = axes[1, 1]
        ax.text(0.1, 0.9, "Analysis Summary:", fontweight='bold', transform=ax.transAxes)
        ax.text(0.1, 0.8, f"• Total contrasts: {len(regional_data)}", transform=ax.transAxes)
        total_sig = sum(data['significant'].sum() for data in regional_data.values() 
                       if 'significant' in data.columns)
        ax.text(0.1, 0.7, f"• Total significant genes: {total_sig}", transform=ax.transAxes)
        ax.text(0.1, 0.6, f"• FDR threshold: {ANALYSIS_CONFIG['fdr_threshold']}", transform=ax.transAxes)
        ax.text(0.1, 0.5, f"• Effect threshold: {ANALYSIS_CONFIG['lfc_threshold']}", transform=ax.transAxes)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        ax.set_title('Analysis Parameters')
        
        plt.tight_layout()
        plt.savefig(f"{PLOTS_DIR}/caudate_vs_putamen_detailed.png", 
                   dpi=ANALYSIS_CONFIG['figure_dpi'], bbox_inches='tight')
        plt.close()
        
        print("      ✅ Caudate vs Putamen plots saved")
        
    except Exception as e:
        print(f"      ❌ Caudate vs Putamen plotting failed: {e}")

def create_regional_summary_plots(results):
    """Create summary statistics plots for regional analysis"""
    print("   Creating regional summary plots...")
    
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
            axes[0, 0].set_title('Significant Genes by Regional Contrast')
            axes[0, 0].set_ylabel('Number of Significant Genes')
            axes[0, 0].set_xticks(range(len(summary_df)))
            axes[0, 0].set_xticklabels(summary_df['contrast'], rotation=45, ha='right')
            
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
        
        # Regional contrasts
        axes[1, 1].text(0.1, 0.9, f"Regional Contrasts:", fontweight='bold', transform=axes[1, 1].transAxes)
        y_pos = 0.8
        for contrast_name in list(DESIGN_CONFIG['contrasts'].keys())[:8]:  # Show first 8
            display_name = contrast_name.replace('_', ' ').title()
            axes[1, 1].text(0.1, y_pos, f"• {display_name}", transform=axes[1, 1].transAxes, fontsize=8)
            y_pos -= 0.08
        axes[1, 1].set_xlim(0, 1)
        axes[1, 1].set_ylim(0, 1)
        axes[1, 1].axis('off')
        axes[1, 1].set_title('Contrasts Analyzed')
        
        plt.tight_layout()
        plt.savefig(f"{PLOTS_DIR}/regional_analysis_summary.png", 
                   dpi=ANALYSIS_CONFIG['figure_dpi'], bbox_inches='tight')
        plt.close()
        
        print("      ✅ Regional summary plots saved")
        
    except Exception as e:
        print(f"      ❌ Regional summary plotting failed: {e}")

# ============================================================================
# 💾 EXPORT FUNCTIONS
# ============================================================================

def export_regional_results(results):
    """Export all regional results"""
    print(f"\n💾 EXPORTING REGIONAL RESULTS")
    print("=" * 31)
    
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
    create_regional_summary_report(results)

def create_regional_summary_report(results):
    """Create comprehensive regional analysis summary"""
    print("   Creating regional summary report...")
    
    summary_lines = [
        "# LEMUR Regional Contrasts Analysis Report",
        "=" * 42,
        "",
        "## Analysis Overview",
        f"This analysis implements regional stratification in LEMUR to understand",
        f"how OUD affects different striatal regions (Caudate vs Putamen).",
        "",
        "## Regional Biology",
        "- **Caudate**: Goal-directed behavior, cognitive control, executive function",
        "- **Putamen**: Habit formation, motor control, addiction maintenance",
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
    regional_summary = {}
    
    for contrast_name, result in results.items():
        if 'de_results' in result:
            de_results = result['de_results']
            n_total = len(de_results)
            n_sig = de_results['significant'].sum() if 'significant' in de_results.columns else 0
            total_sig_genes += n_sig
            
            # Group by region
            if 'caudate' in contrast_name:
                region = 'Caudate'
            elif 'putamen' in contrast_name:
                region = 'Putamen'
            else:
                region = 'Overall'
            
            if region not in regional_summary:
                regional_summary[region] = []
            
            regional_summary[region].append({
                'contrast': contrast_name,
                'n_sig': n_sig,
                'n_total': n_total
            })
    
    # Report by region
    for region, contrasts in regional_summary.items():
        summary_lines.extend([
            f"",
            f"### {region} Region",
        ])
        
        for contrast_info in contrasts:
            contrast_name = contrast_info['contrast']
            n_sig = contrast_info['n_sig']
            n_total = contrast_info['n_total']
            
            summary_lines.extend([
                f"- **{contrast_name.replace('_', ' ').title()}**:",
                f"  - Total genes: {n_total:,}",
                f"  - Significant genes: {n_sig:,}",
                f"  - Percentage: {n_sig/n_total*100:.2f}%",
            ])
            
            if n_sig > 0 and 'de_results' in results[contrast_name]:
                top_genes = results[contrast_name]['de_results'].nlargest(3, 'abs_coefficient')
                summary_lines.append("  - Top genes:")
                for gene, row in top_genes.iterrows():
                    summary_lines.append(f"    - {gene}: coef={row['coefficient']:.3f}")
    
    summary_lines.extend([
        "",
        "## Key Findings",
        f"- Total significant genes across all regional contrasts: {total_sig_genes:,}",
        "",
        "## Regional Comparison",
        "Compare with DESeq2 regional findings for methodological validation.",
        "",
        "## Analysis Parameters",
        f"- Embedding dimensions: {ANALYSIS_CONFIG['n_embedding']}",
        f"- FDR threshold: {ANALYSIS_CONFIG['fdr_threshold']}",
        f"- Effect size threshold: {ANALYSIS_CONFIG['lfc_threshold']}",
        f"- Maximum cells: {ANALYSIS_CONFIG['max_cells']:,}",
        f"- Highly variable genes: {ANALYSIS_CONFIG['n_hvg']:,}",
    ])
    
    # Write summary
    summary_filename = f"{TABLES_DIR}/regional_analysis_summary.md"
    with open(summary_filename, 'w') as f:
        f.write('\n'.join(summary_lines))
    
    print(f"      ✅ Regional summary report: {summary_filename}")

# ============================================================================
# 🚀 MAIN EXECUTION
# ============================================================================

def main():
    """Main regional analysis workflow"""
    print(f"\n🚀 STARTING REGIONAL LEMUR ANALYSIS")
    print("=" * 37)
    
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
    
    # 3. Run regional LEMUR analysis
    results = run_regional_lemur_analysis(adata_prepared)
    if not results:
        print("❌ Regional LEMUR analysis failed")
        return
    
    # 4. Create visualizations
    create_regional_plots(results)
    
    # 5. Export results
    export_regional_results(results)
    
    print(f"\n✅ REGIONAL LEMUR ANALYSIS COMPLETE")
    print("=" * 36)
    print(f"Results saved to: {OUTPUT_DIR}")
    print(f"• Plots: {PLOTS_DIR}/")
    print(f"• Tables: {TABLES_DIR}/")
    
    # Print final summary
    print(f"\n📊 REGIONAL RESULTS SUMMARY:")
    successful_contrasts = 0
    total_sig_genes = 0
    
    caudate_effects = 0
    putamen_effects = 0
    
    for contrast_name, result in results.items():
        if 'de_results' in result:
            successful_contrasts += 1
            de_results = result['de_results']
            n_sig = de_results['significant'].sum() if 'significant' in de_results.columns else 0
            total_sig_genes += n_sig
            
            if 'caudate' in contrast_name:
                caudate_effects += n_sig
            elif 'putamen' in contrast_name:
                putamen_effects += n_sig
            
            print(f"  • {contrast_name}: {n_sig} significant genes")
    
    print(f"  • Total contrasts analyzed: {successful_contrasts}")
    print(f"  • Total significant genes: {total_sig_genes}")
    
    print(f"\n🏗️  REGIONAL EFFECTS SUMMARY:")
    print(f"  • Caudate effects: {caudate_effects} significant genes")
    print(f"  • Putamen effects: {putamen_effects} significant genes")
    
    if caudate_effects > 0 and putamen_effects > 0:
        ratio = caudate_effects / putamen_effects
        print(f"  • Caudate/Putamen ratio: {ratio:.1f}x")
    elif caudate_effects > 0:
        print(f"  • Caudate dominance: exclusive effects")
    elif putamen_effects > 0:
        print(f"  • Putamen dominance: exclusive effects")

if __name__ == "__main__":
    main()