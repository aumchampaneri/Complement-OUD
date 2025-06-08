#!/usr/bin/env python3
"""
🔬 pyDESeq2 Differential Expression Analysis
GSE225158 - OUD vs Control using pyDESeq2 pseudobulk approach

Advantages:
- Same statistical framework as bulk RNA-seq
- Proper sample-level modeling with biological replicates
- Direct comparability with bulk DESeq2 results
- Cell-type specific insights with robust statistics

Pseudobulk Contrasts Generated (Separate Results by Contrast):
1. OUD vs. Control (Pooled)
2. OUD vs. Control <- Putamen
3. OUD vs. Control <- Caudate
4. OUD vs. Control <- Male
5. OUD vs. Control <- Female
6. OUD Effect <- Male vs. Female
7. OUD Effect <- Male vs. Female (Putamen)
8. OUD Effect <- Male vs. Female (Caudate)
9. Control Effect <- Male vs. Female
10. Control Effect <- Male vs. Female (Putamen)
11. Control Effect <- Male vs. Female (Caudate)
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from scipy import sparse
import warnings
warnings.filterwarnings('ignore')

# pyDESeq2 dependencies
try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    PYDESEQ2_AVAILABLE = True
    print("✅ pyDESeq2 available")
except ImportError:
    PYDESEQ2_AVAILABLE = False
    DeseqDataSet = None
    DeseqStats = None
    print("❌ pyDESeq2 not installed. Install with: pip install pydeseq2")

# ============================================================================
# 📁 CONFIGURATION
# ============================================================================

BASE_DIR = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study"
INPUT_H5AD = f"{BASE_DIR}/data/processed/snrna_scvi/GSE225158_annotated_scvi.h5ad"
OUTPUT_DIR = f"{BASE_DIR}/results/snrna_scvi/pydeseq2_analysis"
PLOTS_DIR = f"{OUTPUT_DIR}/plots"
TABLES_DIR = f"{OUTPUT_DIR}/tables"

# Create directories
for directory in [OUTPUT_DIR, PLOTS_DIR, TABLES_DIR]:
    os.makedirs(directory, exist_ok=True)

print("🔬 pyDESeq2 DIFFERENTIAL EXPRESSION ANALYSIS")
print("===========================================")
print(f"Input: {INPUT_H5AD}")
print(f"Output: {OUTPUT_DIR}")

# ============================================================================
# 📂 DATA LOADING AND PREPARATION
# ============================================================================

def load_annotated_data():
    """Load scVI-annotated data and prepare for analysis"""
    print(f"\n📁 LOADING ANNOTATED DATA")
    print("=" * 30)
    
    if not os.path.exists(INPUT_H5AD):
        print(f"   ❌ File not found: {INPUT_H5AD}")
        print("   Make sure to run 02_scvi_annotation.py first")
        return None
    
    # Load the data
    adata = sc.read_h5ad(INPUT_H5AD)
    print(f"   Loaded data: {adata.n_obs} cells × {adata.n_vars} genes")
    
    # Check required columns
    required_cols = ['cell_type', 'condition', 'sample_id']
    for col in required_cols:
        if col not in adata.obs.columns:
            print(f"   ❌ Missing column: {col}")
            return None
    
    # Show data composition
    print(f"   Cell types: {adata.obs['cell_type'].nunique()}")
    print(f"   Conditions: {list(adata.obs['condition'].unique())}")
    print(f"   Samples: {adata.obs['sample_id'].nunique()}")
    
    # Check for sex and region info
    if 'Sex' in adata.obs.columns:
        print(f"   Sex info: {list(adata.obs['Sex'].unique())}")
    if 'Region' in adata.obs.columns:
        print(f"   Region info: {list(adata.obs['Region'].unique())}")
    
    return adata

def create_pseudobulk_data(adata, subset_conditions=None):
    """Create pseudobulk count matrices by aggregating all cells by sample"""
    print(f"\n🧮 CREATING PSEUDOBULK DATA")
    print("   All cell types combined")
    if subset_conditions:
        print(f"   Subset conditions: {subset_conditions}")
    print("=" * 40)
    
    # Use all cell types - this is pseudobulk across the full dataset
    adata_subset = adata.copy()
    
    # Filter to specific conditions if specified
    if subset_conditions:
        condition_mask = adata_subset.obs['condition'].isin(subset_conditions)
        adata_subset = adata_subset[condition_mask].copy()
        if adata_subset.n_obs == 0:
            print(f"   No cells found for conditions: {subset_conditions}")
            return None, None
    
    print(f"   Cells in analysis: {adata_subset.n_obs}")
    
    # Get raw counts (should be in X if from scVI processing)
    if hasattr(adata_subset, 'raw') and adata_subset.raw is not None:
        counts_matrix = adata_subset.raw.X
        gene_names = adata_subset.raw.var_names
    else:
        counts_matrix = adata_subset.X
        gene_names = adata_subset.var_names
    
    # Convert to dense if sparse
    if sparse.issparse(counts_matrix):
        counts_matrix = counts_matrix.toarray()
    
    # Create sample-level aggregation
    sample_ids = adata_subset.obs['sample_id'].unique()
    print(f"   Unique samples: {len(sample_ids)}")
    
    # Aggregate counts by sample
    pseudobulk_counts = []
    sample_metadata = []
    
    for sample_id in sample_ids:
        sample_mask = adata_subset.obs['sample_id'] == sample_id
        sample_cells = adata_subset[sample_mask]
        
        if sample_cells.n_obs > 0:
            # Sum counts across cells in this sample
            if hasattr(sample_cells, 'raw') and sample_cells.raw is not None:
                sample_counts = np.array(sample_cells.raw.X.sum(axis=0)).flatten()
            else:
                sample_counts = np.array(sample_cells.X.sum(axis=0)).flatten()
            
            pseudobulk_counts.append(sample_counts)
            
            # Get sample metadata (take first cell's metadata as representative)
            sample_meta = sample_cells.obs.iloc[0]
            sample_metadata.append({
                'sample_id': sample_id,
                'condition': sample_meta['condition'],
                'donor_id': sample_meta.get('donor_id', 'Unknown'),
                'cell_count': sample_cells.n_obs,
                'sex': sample_meta.get('Sex', 'Unknown'),
                'region': sample_meta.get('Region', 'Unknown')
            })
    
    # Create count matrix and metadata DataFrame
    if len(pseudobulk_counts) > 0:
        count_matrix = np.column_stack(pseudobulk_counts)
        metadata_df = pd.DataFrame(sample_metadata)
        
        # Create count DataFrame - FIXED: use .values to get array
        count_df = pd.DataFrame(
            count_matrix,
            index=gene_names,
            columns=metadata_df['sample_id'].values
        )
        
        print(f"   Pseudobulk matrix shape: {count_df.shape}")
        print(f"   Samples per condition: {dict(metadata_df['condition'].value_counts())}")
        
        return count_df, metadata_df
    else:
        print("   No valid samples found")
        return None, None

# ============================================================================
# 🔬 DESEQ2 ANALYSIS FUNCTIONS
# ============================================================================

def run_single_contrast_analysis(count_df, metadata_df, contrast_name, contrast_spec, design_formula=None):
    """Run differential expression analysis for a single contrast"""
    print(f"\n🔬 RUNNING CONTRAST: {contrast_name}")
    print("=" * 60)
    
    if not PYDESEQ2_AVAILABLE:
        print("   pyDESeq2 not available, skipping analysis")
        return None
    
    # Check sample sizes
    condition_counts = metadata_df['condition'].value_counts()
    print(f"   Sample sizes: {dict(condition_counts)}")
    
    if len(condition_counts) < 2 or any(condition_counts < 2):
        print("   Insufficient samples for DESeq2 analysis (need ≥2 per condition)")
        return None
    
    # Filter low-expressed genes
    print("   Filtering low-expressed genes...")
    min_count = 10
    min_samples = 2
    
    expressed_genes = (count_df >= min_count).sum(axis=1) >= min_samples
    count_df_filtered = count_df[expressed_genes]
    
    print(f"   Genes before filtering: {len(count_df)}")
    print(f"   Genes after filtering: {len(count_df_filtered)}")
    
    if len(count_df_filtered) < 100:
        print("   Too few genes remaining after filtering")
        return None
    
    try:
        # Prepare data for pyDESeq2
        # Note: pyDESeq2 expects samples as rows, genes as columns
        count_matrix_t = count_df_filtered.T
        
        # Ensure metadata order matches count matrix
        sample_names = count_matrix_t.index.tolist()
        metadata_ordered = metadata_df.set_index('sample_id').loc[sample_names]
        
        # Use provided design formula or create default
        if design_formula is None:
            design_factors = ['condition']
            if 'sex' in metadata_ordered.columns and metadata_ordered['sex'].nunique() > 1:
                if not metadata_ordered['sex'].isin(['Unknown']).all():
                    design_factors.append('sex')
            
            if 'region' in metadata_ordered.columns and metadata_ordered['region'].nunique() > 1:
                if not metadata_ordered['region'].isin(['Unknown']).all():
                    design_factors.append('region')
        else:
            # Parse design formula to get factors
            design_factors = [f.strip() for f in design_formula.replace('~', '').split('+')]
        
        design_formula_str = "~ " + " + ".join(design_factors)
        print(f"   Design formula: {design_formula_str}")
        
        # Create DESeq2 dataset
        print("   Creating DESeq2 dataset...")
        if DeseqDataSet is None:
            raise ImportError("pyDESeq2 not available")
        dds = DeseqDataSet(
            counts=count_matrix_t,
            metadata=metadata_ordered,
            design_factors=design_factors,
            refit_cooks=True
        )
        
        # Run DESeq2 analysis
        print("   Running DESeq2 analysis...")
        dds.deseq2()
        
        # Get results based on contrast specification - FIXED: use list instead of tuple
        print("   Computing statistics...")
        if DeseqStats is None:
            raise ImportError("pyDESeq2 not available")
        
        if isinstance(contrast_spec, (tuple, list)) and len(contrast_spec) == 3:
            # Standard contrast: [factor, numerator, denominator]
            contrast_list = list(contrast_spec)  # Convert tuple to list
            stat_res = DeseqStats(dds, contrast=contrast_list)
        else:
            # Complex contrast or interaction - default to condition contrast
            if 'condition' in design_factors:
                stat_res = DeseqStats(dds, contrast=['condition', 'OUD', 'Control'])
            else:
                stat_res = DeseqStats(dds, contrast=['condition', 'OUD', 'Control'])
        
        stat_res.summary()
        
        # Get results DataFrame
        results_df = stat_res.results_df.copy()
        results_df['gene'] = results_df.index
        results_df = results_df.reset_index(drop=True)
        
        # Add additional statistics
        results_df['significant'] = (results_df['padj'] < 0.05) & (abs(results_df['log2FoldChange']) > 0.25)
        
        # Determine direction based on contrast
        if 'Male_vs_Female' in contrast_name:
            results_df['direction'] = np.where(
                results_df['significant'],
                np.where(results_df['log2FoldChange'] > 0, 'Up_in_Male', 'Down_in_Male'),
                'Not_significant'
            )
        else:
            results_df['direction'] = np.where(
                results_df['significant'],
                np.where(results_df['log2FoldChange'] > 0, 'Up_in_OUD', 'Down_in_OUD'),
                'Not_significant'
            )
        
        # Summary statistics
        n_total = len(results_df)
        n_significant = results_df['significant'].sum()
        direction_counts = results_df['direction'].value_counts()
        
        print(f"   Results summary:")
        print(f"     Total genes: {n_total}")
        print(f"     Significant: {n_significant} ({n_significant/n_total*100:.1f}%)")
        print(f"     Direction counts: {dict(direction_counts)}")
        
        return {
            'results': results_df,
            'dds': dds,
            'stats': stat_res,
            'contrast_name': contrast_name,
            'n_samples': len(metadata_ordered),
            'n_genes': n_total,
            'n_significant': n_significant,
            'metadata': metadata_ordered
        }
        
    except Exception as e:
        print(f"   Error in DESeq2 analysis: {str(e)}")
        return None

def filter_metadata_for_contrast(metadata_df, sex_filter=None, region_filter=None, condition_filter=None):
    """Filter metadata based on specified criteria"""
    filtered_df = metadata_df.copy()
    
    if sex_filter:
        filtered_df = filtered_df[filtered_df['sex'] == sex_filter]
    
    if region_filter:
        filtered_df = filtered_df[filtered_df['region'] == region_filter]
        
    if condition_filter:
        filtered_df = filtered_df[filtered_df['condition'].isin(condition_filter)]
    
    return filtered_df

def run_all_contrasts_for_dataset(adata):
    """Run all 11 contrasts for the full dataset (all cell types combined)"""
    print(f"\n🎯 ANALYZING ALL CONTRASTS FOR FULL DATASET")
    print("=" * 80)
    
    all_contrast_results = {}
    
    # Create base pseudobulk data from all cell types
    count_df, metadata_df = create_pseudobulk_data(adata)
    
    if count_df is None:
        print(f"   Skipping analysis - insufficient data")
        return None
    
    # Define all contrasts - FIXED: use lists instead of tuples for contrasts
    contrasts = [
        {
            'name': '01_OUD_vs_Control_Pooled',
            'description': 'OUD vs. Control (Pooled)',
            'metadata_filter': {},
            'contrast': ['condition', 'OUD', 'Control']
        },
        {
            'name': '02_OUD_vs_Control_Putamen',
            'description': 'OUD vs. Control <- Putamen',
            'metadata_filter': {'region_filter': 'Putamen'},
            'contrast': ['condition', 'OUD', 'Control']
        },
        {
            'name': '03_OUD_vs_Control_Caudate',
            'description': 'OUD vs. Control <- Caudate',
            'metadata_filter': {'region_filter': 'Caudate'},
            'contrast': ['condition', 'OUD', 'Control']
        },
        {
            'name': '04_OUD_vs_Control_Male',
            'description': 'OUD vs. Control <- Male',
            'metadata_filter': {'sex_filter': 'Male'},
            'contrast': ['condition', 'OUD', 'Control']
        },
        {
            'name': '05_OUD_vs_Control_Female',
            'description': 'OUD vs. Control <- Female',
            'metadata_filter': {'sex_filter': 'Female'},
            'contrast': ['condition', 'OUD', 'Control']
        },
        {
            'name': '06_OUD_Effect_Male_vs_Female',
            'description': 'OUD Effect <- Male vs. Female',
            'metadata_filter': {'condition_filter': ['OUD']},
            'contrast': ['sex', 'Male', 'Female']
        },
        {
            'name': '07_OUD_Effect_Male_vs_Female_Putamen',
            'description': 'OUD Effect <- Male vs. Female (Putamen)',
            'metadata_filter': {'condition_filter': ['OUD'], 'region_filter': 'Putamen'},
            'contrast': ['sex', 'Male', 'Female']
        },
        {
            'name': '08_OUD_Effect_Male_vs_Female_Caudate',
            'description': 'OUD Effect <- Male vs. Female (Caudate)',
            'metadata_filter': {'condition_filter': ['OUD'], 'region_filter': 'Caudate'},
            'contrast': ['sex', 'Male', 'Female']
        },
        {
            'name': '09_Control_Effect_Male_vs_Female',
            'description': 'Control Effect <- Male vs. Female',
            'metadata_filter': {'condition_filter': ['Control']},
            'contrast': ['sex', 'Male', 'Female']
        },
        {
            'name': '10_Control_Effect_Male_vs_Female_Putamen',
            'description': 'Control Effect <- Male vs. Female (Putamen)',
            'metadata_filter': {'condition_filter': ['Control'], 'region_filter': 'Putamen'},
            'contrast': ['sex', 'Male', 'Female']
        },
        {
            'name': '11_Control_Effect_Male_vs_Female_Caudate',
            'description': 'Control Effect <- Male vs. Female (Caudate)',
            'metadata_filter': {'condition_filter': ['Control'], 'region_filter': 'Caudate'},
            'contrast': ['sex', 'Male', 'Female']
        }
    ]
    
    # Run each contrast
    for contrast_info in contrasts:
        print(f"\n{'='*60}")
        print(f"CONTRAST {contrast_info['name']}: {contrast_info['description']}")
        print(f"{'='*60}")
        
        # Filter metadata and counts for this contrast
        filtered_metadata = filter_metadata_for_contrast(metadata_df, **contrast_info['metadata_filter'])
        
        # Check if we have enough samples
        if len(filtered_metadata) < 4:  # Need at least 2 samples per group
            print(f"   Insufficient samples after filtering: {len(filtered_metadata)}")
            continue
            
        # Filter count matrix to match filtered metadata
        sample_ids = filtered_metadata['sample_id'].tolist()
        filtered_counts = count_df[sample_ids]
        
        # Run the analysis
        result = run_single_contrast_analysis(
            filtered_counts,
            filtered_metadata,
            contrast_info['name'],
            contrast_info['contrast']
        )
        
        if result is not None:
            all_contrast_results[contrast_info['name']] = result
            
            # Save individual contrast results
            save_contrast_results(result, "Full_Dataset", contrast_info['name'])
    
    return all_contrast_results

# ============================================================================
# 📊 PLOTTING AND VISUALIZATION
# ============================================================================

def plot_contrast_results(result, analysis_type, contrast_name):
    """Create plots for contrast results"""
    if result is None:
        return
    
    results_df = result['results']
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle(f'{analysis_type} - {contrast_name}', fontsize=16, fontweight='bold')
    
    # 1. Volcano plot
    ax1 = axes[0, 0]
    x = results_df['log2FoldChange']
    y = -np.log10(results_df['padj'].fillna(1))
    
    # Color points
    colors = []
    for _, row in results_df.iterrows():
        if row['significant']:
            if row['log2FoldChange'] > 0:
                colors.append('red')
            else:
                colors.append('blue')
        else:
            colors.append('gray')
    
    ax1.scatter(x, y, c=colors, alpha=0.6, s=20)
    ax1.axhline(-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
    ax1.axvline(0.25, color='black', linestyle='--', alpha=0.5)
    ax1.axvline(-0.25, color='black', linestyle='--', alpha=0.5)
    ax1.set_xlabel('log2(Fold Change)')
    ax1.set_ylabel('-log10(padj)')
    ax1.set_title('Volcano Plot')
    
    # 2. MA plot
    ax2 = axes[0, 1]
    x_ma = results_df['baseMean']
    y_ma = results_df['log2FoldChange']
    ax2.scatter(np.log10(x_ma + 1), y_ma, c=colors, alpha=0.6, s=20)
    ax2.axhline(0, color='black', linestyle='-', alpha=0.5)
    ax2.set_xlabel('log10(baseMean)')
    ax2.set_ylabel('log2(Fold Change)')
    ax2.set_title('MA Plot')
    
    # 3. Significance barplot
    ax3 = axes[1, 0]
    direction_counts = results_df['direction'].value_counts()
    bars = ax3.bar(range(len(direction_counts)), direction_counts.values)
    ax3.set_xticks(range(len(direction_counts)))
    ax3.set_xticklabels(direction_counts.index, rotation=45, ha='right')
    ax3.set_ylabel('Number of Genes')
    ax3.set_title('Gene Direction Summary')
    
    # Add value labels on bars
    for bar, value in zip(bars, direction_counts.values):
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                str(value), ha='center', va='bottom')
    
    # 4. P-value distribution
    ax4 = axes[1, 1]
    pvals = results_df['pvalue'].dropna()
    ax4.hist(pvals, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
    ax4.set_xlabel('P-value')
    ax4.set_ylabel('Frequency')
    ax4.set_title('P-value Distribution')
    
    plt.tight_layout()
    
    # Save plot
    plot_filename = f"{PLOTS_DIR}/{analysis_type}_{contrast_name}_analysis.png"
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"   Plot saved: {plot_filename}")

# ============================================================================
# 💾 EXPORT AND SAVE FUNCTIONS
# ============================================================================

def save_contrast_results(result, analysis_type, contrast_name):
    """Save individual contrast results"""
    if result is None:
        return
    
    # Create analysis type directory
    analysis_dir = f"{TABLES_DIR}/{analysis_type}"
    os.makedirs(analysis_dir, exist_ok=True)
    
    # Save results table
    results_filename = f"{analysis_dir}/{contrast_name}_results.csv"
    result['results'].to_csv(results_filename, index=False)
    
    # Save summary stats
    summary_filename = f"{analysis_dir}/{contrast_name}_summary.txt"
    with open(summary_filename, 'w') as f:
        f.write(f"Contrast: {contrast_name}\n")
        f.write(f"Analysis Type: {analysis_type}\n")
        f.write(f"Total Samples: {result['n_samples']}\n")
        f.write(f"Total Genes: {result['n_genes']}\n")
        f.write(f"Significant Genes: {result['n_significant']}\n")
        f.write(f"Percent Significant: {result['n_significant']/result['n_genes']*100:.2f}%\n\n")
        
        # Direction summary
        direction_counts = result['results']['direction'].value_counts()
        f.write("Direction Summary:\n")
        for direction, count in direction_counts.items():
            f.write(f"  {direction}: {count}\n")
    
    print(f"   Results saved: {results_filename}")
    print(f"   Summary saved: {summary_filename}")

def export_all_results(all_contrast_results):
    """Export comprehensive summary of all results"""
    print(f"\n💾 EXPORTING COMPREHENSIVE RESULTS")
    print("=" * 40)
    
    # Create master summary
    summary_data = []
    
    for contrast_name, result in all_contrast_results.items():
        if result is not None:
            direction_counts = result['results']['direction'].value_counts()
            
            summary_data.append({
                'contrast': contrast_name,
                'n_samples': result['n_samples'],
                'n_genes': result['n_genes'],
                'n_significant': result['n_significant'],
                'percent_significant': result['n_significant']/result['n_genes']*100,
                'n_upregulated': direction_counts.get('Up_in_OUD', 0) + direction_counts.get('Up_in_Male', 0),
                'n_downregulated': direction_counts.get('Down_in_OUD', 0) + direction_counts.get('Down_in_Male', 0),
                'n_not_significant': direction_counts.get('Not_significant', 0)
            })
    
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_filename = f"{TABLES_DIR}/master_summary_all_contrasts.csv"
        summary_df.to_csv(summary_filename, index=False)
        
        print(f"   Master summary saved: {summary_filename}")
        
        # Print summary table
        print("\n📊 MASTER SUMMARY:")
        print(summary_df.to_string(index=False))
        
        return summary_df
    
    return None

# ============================================================================
# 🚀 MAIN EXECUTION FUNCTIONS
# ============================================================================

def run_full_dataset_analysis(adata):
    """Run complete analysis for the full dataset with all contrasts"""
    return run_all_contrasts_for_dataset(adata)

def main():
    """Main execution function"""
    print("🚀 STARTING COMPREHENSIVE pyDESeq2 ANALYSIS PIPELINE")
    print("====================================================")
    print("Will run all 11 contrasts for each viable cell type")
    
    if not PYDESEQ2_AVAILABLE:
        print("\n❌ pyDESeq2 not available!")
        print("Install with: pip install pydeseq2")
        return None
    
    try:
        # Load and prepare data
        adata = load_annotated_data()
        if adata is None:
            return None
        
        # Check data composition
        cell_type_counts = adata.obs['cell_type'].value_counts()
        sample_counts = adata.obs['sample_id'].nunique()
        condition_counts = adata.obs['condition'].value_counts()
        
        print(f"\n📊 DATASET COMPOSITION:")
        print(f"   Total cells: {adata.n_obs}")
        print(f"   Total samples: {sample_counts}")
        print(f"   Conditions: {dict(condition_counts)}")
        print(f"   Cell types: {len(cell_type_counts)}")
        
        # Run analysis for full dataset (all cell types combined)
        print(f"\n{'#'*80}")
        print(f"PROCESSING FULL DATASET (ALL CELL TYPES COMBINED)")
        print(f"{'#'*80}")
        
        all_contrast_results = run_full_dataset_analysis(adata)
        if all_contrast_results is not None and len(all_contrast_results) > 0:
            # Create plots for each contrast
            for contrast_name, result in all_contrast_results.items():
                plot_contrast_results(result, "Full_Dataset", contrast_name)
        
        # Export comprehensive results
        summary_df = export_all_results(all_contrast_results)
        
        # Print final summary
        print("\n✅ COMPREHENSIVE pyDESeq2 ANALYSIS COMPLETED!")
        print("=" * 50)
        print(f"Analyzed full dataset with all contrasts")
        print(f"Results organized by contrast in: {TABLES_DIR}")
        print(f"Plots saved to: {PLOTS_DIR}")
        
        if summary_df is not None:
            print(f"\nTotal contrasts run: {len(summary_df)}")
            print("Results saved with complete contrast breakdown!")
        
        return all_contrast_results
        
    except Exception as e:
        print(f"\n❌ ERROR: {str(e)}")
        raise e

# ============================================================================
# 🎯 SCRIPT EXECUTION
# ============================================================================

if __name__ == "__main__":
    main()