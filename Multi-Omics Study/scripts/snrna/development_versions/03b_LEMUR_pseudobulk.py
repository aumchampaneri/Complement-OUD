#!/usr/bin/env python3
"""
🌊 LEMUR Pseudobulk Analysis - Latent Expression Modeling for Unified Representation
GSE225158 - OUD vs Control using pyLEMUR with pseudobulk aggregation and sex interactions

This script performs LEMUR analysis on pseudobulked single-cell data:
1. Aggregates single-cell data by sample/donor and condition
2. Includes sex interaction modeling
3. Performs multiple contrasts: condition, sex, and interaction effects
4. Uses proper biological replicates for statistical modeling

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
OUTPUT_DIR = f"{BASE_DIR}/results/snrna_scvi/lemur_pseudobulk_analysis"
PLOTS_DIR = f"{OUTPUT_DIR}/plots"
TABLES_DIR = f"{OUTPUT_DIR}/tables"

# Analysis parameters
ANALYSIS_CONFIG = {
    'n_embedding': 15,           # Number of LEMUR embedding dimensions
    'n_hvg': 3000,              # Number of highly variable genes
    'min_cells_per_sample': 50,  # Minimum cells per sample for pseudobulk
    'random_seed': 42,          # Random seed for reproducibility
    
    # Significance thresholds
    'fdr_lenient': 0.1,         # Lenient FDR threshold
    'lfc_lenient': 0.1,         # Lenient log fold change threshold
    'fdr_strict': 0.05,         # Strict FDR threshold
    'lfc_strict': 0.25,         # Strict log fold change threshold
    
    # Plotting parameters
    'figure_dpi': 300,          # Plot resolution
    'point_size': 10,           # Scatter plot point size
    'alpha': 0.7,               # Point transparency
}

# Create output directories
for directory in [OUTPUT_DIR, PLOTS_DIR, TABLES_DIR]:
    os.makedirs(directory, exist_ok=True)

print("🌊 LEMUR PSEUDOBULK ANALYSIS")
print("===========================")
print(f"Input: {INPUT_H5AD}")
print(f"Output: {OUTPUT_DIR}")

# ============================================================================
# 📊 DATA LOADING AND PSEUDOBULK AGGREGATION
# ============================================================================

def load_and_prepare_data():
    """Load scVI-annotated data and prepare for pseudobulk analysis"""
    print(f"\n📁 LOADING ANNOTATED DATA")
    print("=" * 30)
    
    if not os.path.exists(INPUT_H5AD):
        print(f"   ❌ File not found: {INPUT_H5AD}")
        print("   Make sure to run 02_scvi_annotation.py first")
        return None
    
    # Load the data
    adata = sc.read_h5ad(INPUT_H5AD)
    print(f"   Loaded data: {adata.n_obs} cells × {adata.n_vars} genes")
    
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
    
    # Standardize sex column
    if 'Sex' in adata.obs.columns:
        adata.obs['sex'] = adata.obs['Sex']
    
    # Create sample_id if needed
    if 'sample_id' not in adata.obs.columns and 'ID' in adata.obs.columns:
        print("   Creating sample_id from ID column...")
        adata.obs['sample_id'] = adata.obs['ID'].astype(str)
    
    # Display data composition
    print(f"   Conditions: {dict(adata.obs['condition'].value_counts())}")
    if 'sex' in adata.obs.columns:
        print(f"   Sex distribution: {dict(adata.obs['sex'].value_counts())}")
        
        # Cross-tabulation
        print("   Condition x Sex crosstab:")
        crosstab = pd.crosstab(adata.obs['condition'], adata.obs['sex'])
        print(crosstab)
    
    print(f"   Cell types: {adata.obs['cell_type'].nunique()}")
    if 'sample_id' in adata.obs.columns:
        print(f"   Samples: {adata.obs['sample_id'].nunique()}")
    
    return adata

def create_pseudobulk_data(adata):
    """Create pseudobulk aggregation by sample, condition, and sex"""
    print(f"\n🔄 CREATING PSEUDOBULK DATA")
    print("=" * 30)
    
    # Select highly variable genes first to reduce computational burden
    print("   Identifying highly variable genes...")
    sc.pp.highly_variable_genes(adata, n_top_genes=ANALYSIS_CONFIG['n_hvg'], flavor='seurat_v3')
    adata_hvg = adata[:, adata.var.highly_variable].copy()
    print(f"   Using {adata_hvg.n_vars} highly variable genes")
    
    # Get raw counts for aggregation
    if adata_hvg.raw is not None:
        counts_matrix = adata_hvg.raw.X
        gene_names = adata_hvg.raw.var_names
    else:
        counts_matrix = adata_hvg.X
        gene_names = adata_hvg.var_names
    
    # Convert to dense if sparse
    if sparse.issparse(counts_matrix):
        counts_matrix = counts_matrix.toarray()
    
    # Create aggregation groups
    print("   Aggregating counts by sample, condition, and sex...")
    
    # Required columns for pseudobulk
    required_cols = ['sample_id', 'condition', 'sex']
    for col in required_cols:
        if col not in adata_hvg.obs.columns:
            print(f"   ❌ Missing required column: {col}")
            return None, None
    
    # Create grouping key
    adata_hvg.obs['pseudobulk_group'] = (
        adata_hvg.obs['sample_id'].astype(str) + '_' +
        adata_hvg.obs['condition'].astype(str) + '_' +
        adata_hvg.obs['sex'].astype(str)
    )
    
    # Aggregate counts
    pseudobulk_counts = []
    pseudobulk_metadata = []
    
    for group in adata_hvg.obs['pseudobulk_group'].unique():
        mask = adata_hvg.obs['pseudobulk_group'] == group
        
        # Check minimum cell count
        n_cells = mask.sum()
        if n_cells < ANALYSIS_CONFIG['min_cells_per_sample']:
            print(f"   ⚠️  Skipping {group}: only {n_cells} cells")
            continue
        
        # Aggregate counts
        group_counts = np.sum(counts_matrix[mask, :], axis=0)
        pseudobulk_counts.append(group_counts)
        
        # Extract metadata
        group_obs = adata_hvg.obs[mask].iloc[0]
        pseudobulk_metadata.append({
            'pseudobulk_group': group,
            'sample_id': group_obs['sample_id'],
            'condition': group_obs['condition'],
            'sex': group_obs['sex'],
            'n_cells': n_cells
        })
    
    # Create pseudobulk AnnData object
    pseudobulk_counts = np.array(pseudobulk_counts)
    pseudobulk_metadata = pd.DataFrame(pseudobulk_metadata)
    
    print(f"   Created pseudobulk data: {pseudobulk_counts.shape[0]} samples × {pseudobulk_counts.shape[1]} genes")
    print(f"   Sample composition:")
    print(pseudobulk_metadata.groupby(['condition', 'sex']).size())
    
    # Create AnnData object
    adata_pseudobulk = ad.AnnData(
        X=pseudobulk_counts,
        obs=pseudobulk_metadata,
        var=pd.DataFrame(index=gene_names)
    )
    
    return adata_pseudobulk, adata_hvg

# ============================================================================
# 🌊 LEMUR ANALYSIS FUNCTIONS
# ============================================================================

def run_lemur_contrasts(adata_pseudobulk):
    """Run multiple LEMUR contrasts for comprehensive analysis"""
    print("\n🌊 RUNNING LEMUR CONTRASTS")
    print("=" * 40)
    
    if not PYLEMUR_AVAILABLE:
        print("   ❌ pyLEMUR not available")
        return None
    
    # Normalize pseudobulk data
    print("   Normalizing pseudobulk data...")
    X_norm = normalize_pseudobulk_expression(adata_pseudobulk.X)
    
    # Create normalized AnnData
    adata_norm = ad.AnnData(X=X_norm, obs=adata_pseudobulk.obs.copy(), var=adata_pseudobulk.var.copy())
    
    # Define contrasts to run
    contrasts = [
        {
            'name': 'condition_effect',
            'description': 'OUD vs Control (main effect)',
            'design': '~ condition',
            'coefficient_index': 1  # condition effect
        },
        {
            'name': 'sex_effect', 
            'description': 'Male vs Female (main effect)',
            'design': '~ sex',
            'coefficient_index': 1  # sex effect
        },
        {
            'name': 'condition_sex_interaction',
            'description': 'Condition + Sex + Interaction',
            'design': '~ condition + sex + condition:sex',
            'coefficient_index': 3  # interaction term
        },
        {
            'name': 'condition_in_males',
            'description': 'OUD vs Control in Males only',
            'design': '~ condition',
            'coefficient_index': 1,
            'filter': {'sex': 'M'}
        },
        {
            'name': 'condition_in_females',
            'description': 'OUD vs Control in Females only', 
            'design': '~ condition',
            'coefficient_index': 1,
            'filter': {'sex': 'F'}
        }
    ]
    
    results = {}
    
    for contrast in contrasts:
        print(f"\n   🔬 Running contrast: {contrast['description']}")
        
        # Apply filter if specified
        if 'filter' in contrast:
            filter_col, filter_val = list(contrast['filter'].items())[0]
            mask = adata_norm.obs[filter_col] == filter_val
            adata_filtered = adata_norm[mask].copy()
            print(f"      Filtered to {adata_filtered.n_obs} samples")
        else:
            adata_filtered = adata_norm
        
        # Skip if too few samples
        if adata_filtered.n_obs < 4:
            print(f"      ⚠️  Skipping: too few samples ({adata_filtered.n_obs})")
            continue
        
        try:
            # Create LEMUR model
            lemur_model = pylemur.tl.LEMUR(
                adata_filtered,
                design=contrast['design'],
                n_embedding=ANALYSIS_CONFIG['n_embedding']
            )
            
            # Fit model
            lemur_model.fit(verbose=False)
            
            # Extract results
            result = extract_lemur_results(
                lemur_model, 
                adata_filtered, 
                contrast['coefficient_index'],
                contrast['name']
            )
            
            if result is not None:
                result['contrast_info'] = contrast
                results[contrast['name']] = result
                print(f"      ✅ Completed: {result['n_sig']} significant genes")
            
        except Exception as e:
            print(f"      ❌ Failed: {str(e)}")
            continue
    
    return results

def normalize_pseudobulk_expression(X):
    """Normalize pseudobulk expression matrix"""
    # Library size normalization
    lib_sizes = np.sum(X, axis=1)
    X_norm = X / lib_sizes[:, np.newaxis] * np.median(lib_sizes)
    
    # Log transformation
    X_norm = np.log1p(X_norm)
    
    return X_norm

def extract_lemur_results(lemur_model, adata, coeff_index, contrast_name):
    """Extract differential expression results from LEMUR model"""
    try:
        # Get linear coefficients
        linear_coeffs = lemur_model.linear_coefficients
        
        if coeff_index >= linear_coeffs.shape[0]:
            print(f"      ⚠️  Coefficient index {coeff_index} out of range")
            return None
        
        # Extract effects for the specified coefficient
        effects = linear_coeffs[coeff_index, :]
        
        # Calculate statistics
        gene_names = adata.var.index
        
        # Simple t-test approach for p-values
        stderr = np.std(effects) + 1e-8
        t_stats = effects / stderr
        p_vals = 2 * (1 - stats.norm.cdf(np.abs(t_stats)))
        
        # FDR correction
        valid_pvals = ~np.isnan(p_vals) & (p_vals >= 0) & (p_vals <= 1)
        adj_pvals = np.full(len(p_vals), 1.0)
        
        if valid_pvals.sum() > 0:
            _, adj_pvals[valid_pvals], _, _ = multipletests(
                p_vals[valid_pvals], method='fdr_bh'
            )
        
        # Create results DataFrame
        de_df = pd.DataFrame({
            'gene': gene_names,
            'lfc': effects,
            'stat': t_stats,
            'pval': p_vals,
            'adj_pval': adj_pvals
        })
        
        # Add significance flags
        de_df['significant'] = (
            (de_df['adj_pval'] < ANALYSIS_CONFIG['fdr_lenient']) & 
            (abs(de_df['lfc']) > ANALYSIS_CONFIG['lfc_lenient'])
        )
        de_df['significant_strict'] = (
            (de_df['adj_pval'] < ANALYSIS_CONFIG['fdr_strict']) & 
            (abs(de_df['lfc']) > ANALYSIS_CONFIG['lfc_strict'])
        )
        
        return {
            'de_df': de_df,
            'embedding': lemur_model.embedding,
            'lemur_model': lemur_model,
            'adata_used': adata,
            'n_sig': de_df['significant'].sum(),
            'n_sig_strict': de_df['significant_strict'].sum(),
            'contrast_name': contrast_name
        }
        
    except Exception as e:
        print(f"      ❌ Result extraction failed: {str(e)}")
        return None

# ============================================================================
# 📊 VISUALIZATION FUNCTIONS
# ============================================================================

def create_comprehensive_plots(results):
    """Create comprehensive visualization of all contrasts"""
    print(f"\n📊 CREATING COMPREHENSIVE VISUALIZATIONS")
    print("=" * 40)
    
    n_contrasts = len(results)
    if n_contrasts == 0:
        print("   No results to plot")
        return
    
    # Create figure with subplots for each contrast
    fig, axes = plt.subplots(n_contrasts, 4, figsize=(20, 5*n_contrasts))
    if n_contrasts == 1:
        axes = axes.reshape(1, -1)
    
    for i, (contrast_name, result) in enumerate(results.items()):
        de_df = result['de_df']
        embedding = result['embedding']
        adata_used = result['adata_used']
        contrast_info = result['contrast_info']
        
        # Embedding plot
        ax = axes[i, 0]
        scatter_embedding_by_groups(ax, embedding, adata_used, contrast_info['description'])
        
        # Volcano plot
        ax = axes[i, 1]
        create_volcano_plot(ax, de_df, contrast_info['description'])
        
        # Top genes plot
        ax = axes[i, 2]
        plot_top_genes(ax, de_df, contrast_info['description'])
        
        # Summary statistics
        ax = axes[i, 3]
        plot_summary_stats(ax, de_df, result, contrast_info['description'])
    
    plt.tight_layout()
    plot_filename = f"{PLOTS_DIR}/lemur_pseudobulk_comprehensive.png"
    plt.savefig(plot_filename, dpi=ANALYSIS_CONFIG['figure_dpi'], bbox_inches='tight')
    plt.close()
    print(f"   Comprehensive plots saved: {plot_filename}")

def scatter_embedding_by_groups(ax, embedding, adata, title):
    """Create embedding scatter plot colored by groups"""
    # Determine coloring strategy based on available columns
    if 'condition' in adata.obs.columns and 'sex' in adata.obs.columns:
        # Create combined condition-sex groups
        groups = adata.obs['condition'] + '_' + adata.obs['sex']
        colors = {'OUD_M': 'darkred', 'OUD_F': 'lightcoral', 
                 'Control_M': 'darkblue', 'Control_F': 'lightblue'}
    elif 'condition' in adata.obs.columns:
        groups = adata.obs['condition']
        colors = {'OUD': 'red', 'Control': 'blue'}
    elif 'sex' in adata.obs.columns:
        groups = adata.obs['sex']
        colors = {'M': 'blue', 'F': 'pink'}
    else:
        groups = ['All'] * len(adata.obs)
        colors = {'All': 'gray'}
    
    for group in groups.unique():
        mask = groups == group
        color = colors.get(group, 'gray')
        ax.scatter(embedding[mask, 0], embedding[mask, 1], 
                  c=color, label=group, alpha=ANALYSIS_CONFIG['alpha'], 
                  s=ANALYSIS_CONFIG['point_size'])
    
    ax.set_xlabel('LEMUR Embedding 1')
    ax.set_ylabel('LEMUR Embedding 2')
    ax.set_title(f'Embedding - {title}')
    ax.legend()

def create_volcano_plot(ax, de_df, title):
    """Create volcano plot"""
    x = de_df['lfc']
    y = -np.log10(de_df['adj_pval'].clip(lower=1e-300))
    
    # Color by significance
    colors = np.where(de_df['significant'],
                     np.where(de_df['lfc'] > 0, 'red', 'blue'), 'lightgray')
    
    ax.scatter(x, y, c=colors, alpha=ANALYSIS_CONFIG['alpha'], s=5)
    ax.axhline(-np.log10(ANALYSIS_CONFIG['fdr_lenient']), color='red', 
               linestyle='--', alpha=0.7, label=f"FDR={ANALYSIS_CONFIG['fdr_lenient']}")
    ax.axvline(ANALYSIS_CONFIG['lfc_lenient'], color='black', linestyle='--', alpha=0.5)
    ax.axvline(-ANALYSIS_CONFIG['lfc_lenient'], color='black', linestyle='--', alpha=0.5)
    
    ax.set_xlabel('Log2 Fold Change')
    ax.set_ylabel('-log10(Adjusted p-value)')
    ax.set_title(f'Volcano - {title}')

def plot_top_genes(ax, de_df, title):
    """Plot top genes by effect size"""
    # Get top genes by absolute LFC
    abs_lfc = de_df['lfc'].abs()
    top_indices = abs_lfc.nlargest(10).index
    top_genes = de_df.loc[top_indices].copy()
    
    # Create barplot
    colors = ['red' if lfc > 0 else 'blue' for lfc in top_genes['lfc']]
    bars = ax.barh(range(len(top_genes)), top_genes['lfc'], color=colors, alpha=0.7)
    
    ax.set_yticks(range(len(top_genes)))
    ax.set_yticklabels(top_genes['gene'])
    ax.set_xlabel('Log2 Fold Change')
    ax.set_title(f'Top Genes - {title}')
    ax.axvline(0, color='black', linestyle='-', alpha=0.3)

def plot_summary_stats(ax, de_df, result, title):
    """Plot summary statistics"""
    ax.axis('off')
    
    n_total = len(de_df)
    n_sig = result['n_sig']
    n_sig_strict = result['n_sig_strict']
    n_up = ((de_df['significant']) & (de_df['lfc'] > 0)).sum()
    n_down = ((de_df['significant']) & (de_df['lfc'] < 0)).sum()
    
    summary_text = f"""
    {title}
    {'='*40}
    
    Total genes: {n_total:,}
    Significant (lenient): {n_sig:,} ({n_sig/n_total*100:.1f}%)
    Significant (strict): {n_sig_strict:,} ({n_sig_strict/n_total*100:.1f}%)
    
    Upregulated: {n_up:,}
    Downregulated: {n_down:,}
    
    LFC range: {de_df['lfc'].min():.3f} to {de_df['lfc'].max():.3f}
    """
    
    ax.text(0.05, 0.95, summary_text, transform=ax.transAxes, 
           verticalalignment='top', fontsize=10, family='monospace')

# ============================================================================
# 💾 EXPORT FUNCTIONS  
# ============================================================================

def export_all_results(results):
    """Export all contrast results"""
    print(f"\n💾 EXPORTING RESULTS")
    print("=" * 40)
    
    exported_files = {}
    
    for contrast_name, result in results.items():
        print(f"   Exporting {contrast_name}...")
        
        de_df = result['de_df']
        
        # Export full results
        full_filename = f"{TABLES_DIR}/{contrast_name}_differential_expression.csv"
        de_df.to_csv(full_filename, index=False)
        
        # Export significant genes
        sig_genes = de_df[de_df['significant']].copy().sort_values('adj_pval')
        sig_filename = f"{TABLES_DIR}/{contrast_name}_significant_genes.csv"
        sig_genes.to_csv(sig_filename, index=False)
        
        # Export embedding if available
        if 'embedding' in result:
            embedding_df = pd.DataFrame(
                result['embedding'],
                columns=[f'LEMUR_{i+1}' for i in range(result['embedding'].shape[1])]
            )
            # Add metadata
            for col in result['adata_used'].obs.columns:
                embedding_df[col] = result['adata_used'].obs[col].values
            
            emb_filename = f"{TABLES_DIR}/{contrast_name}_embedding.csv"
            embedding_df.to_csv(emb_filename, index=False)
        
        exported_files[contrast_name] = {
            'full_results': full_filename,
            'significant_genes': sig_filename,
            'embedding': emb_filename if 'embedding' in result else None
        }
    
    # Create master summary
    create_master_summary(results)
    
    return exported_files

def create_master_summary(results):
    """Create master summary of all contrasts"""
    summary_data = []
    
    for contrast_name, result in results.items():
        de_df = result['de_df']
        contrast_info = result['contrast_info']
        
        summary_data.append({
            'Contrast': contrast_name,
            'Description': contrast_info['description'],
            'Design_Formula': contrast_info['design'],
            'Total_Genes': len(de_df),
            'Significant_Lenient': result['n_sig'],
            'Significant_Strict': result['n_sig_strict'],
            'Percent_Significant': result['n_sig'] / len(de_df) * 100,
            'Max_Abs_LFC': de_df['lfc'].abs().max(),
            'Samples_Used': result['adata_used'].n_obs if 'adata_used' in result else 0
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_filename = f"{TABLES_DIR}/master_summary_all_contrasts.csv"
    summary_df.to_csv(summary_filename, index=False)
    print(f"   Master summary saved: {summary_filename}")

# ============================================================================
# 🚀 MAIN EXECUTION
# ============================================================================

def main():
    """Main execution function"""
    print("🚀 STARTING LEMUR PSEUDOBULK ANALYSIS PIPELINE")
    print("=" * 50)
    
    try:
        # Load and prepare data
        adata = load_and_prepare_data()
        if adata is None:
            return None
        
        # Create pseudobulk data
        adata_pseudobulk, adata_hvg = create_pseudobulk_data(adata)
        if adata_pseudobulk is None:
            return None
        
        print(f"\n📊 PSEUDOBULK DATASET COMPOSITION:")
        print(f"   Total samples: {adata_pseudobulk.n_obs}")
        print(f"   Total genes: {adata_pseudobulk.n_vars}")
        print(f"   Condition distribution: {dict(adata_pseudobulk.obs['condition'].value_counts())}")
        print(f"   Sex distribution: {dict(adata_pseudobulk.obs['sex'].value_counts())}")
        
        # Run LEMUR contrasts
        results = run_lemur_contrasts(adata_pseudobulk)
        
        if results and len(results) > 0:
            # Create visualizations
            create_comprehensive_plots(results)
            
            # Export results
            exported_files = export_all_results(results)
            
            # Print summary
            print("\n✅ LEMUR PSEUDOBULK ANALYSIS COMPLETED!")
            print("=" * 50)
            print(f"Results saved to: {OUTPUT_DIR}")
            
            for contrast_name, result in results.items():
                print(f"\n{contrast_name}:")
                print(f"  Description: {result['contrast_info']['description']}")
                print(f"  Significant genes: {result['n_sig']} (lenient), {result['n_sig_strict']} (strict)")
                
                # Show top genes
                de_df = result['de_df']
                if result['n_sig'] > 0:
                    top_genes = de_df[de_df['significant']].nsmallest(3, 'adj_pval')
                    print("  Top genes:")
                    for _, gene in top_genes.iterrows():
                        direction = "↑" if gene['lfc'] > 0 else "↓"
                        print(f"    {direction} {gene['gene']}: LFC={gene['lfc']:.3f}, padj={gene['adj_pval']:.2e}")
        
        return results
        
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