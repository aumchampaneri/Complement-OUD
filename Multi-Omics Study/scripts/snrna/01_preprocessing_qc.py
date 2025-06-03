'''
🔬 Single-nucleus RNA-seq Processing Pipeline
GSE225158 - OUD Striatum snRNA-seq data reprocessing

Dataset structure:
- 98,848 cells × 31,393 genes
- Rich clinical metadata (OUD diagnosis, demographics, etc.)
- Existing cell type annotations
- Raw counts stored in .raw

Workflow:
1. Extract raw counts and reset processing
2. Quality Control with robust filtering
3. Pearson residuals normalization pipeline
4. Dimensionality reduction and clustering
5. Integration with existing metadata
'''

# =======================================
# 📁 INPUT/OUTPUT PATHS CONFIGURATION
# =======================================

# Input data path
INPUT_PATH = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study/data/raw/snrna/GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5ad"

# Output directory structure
BASE_OUTPUT_DIR = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study"
PROCESSED_DATA_DIR = f"{BASE_OUTPUT_DIR}/data/processed/snrna"
RESULTS_DIR = f"{BASE_OUTPUT_DIR}/results/snrna"
PLOTS_DIR = f"{RESULTS_DIR}/qc_plots"

# Output file names
OUTPUT_H5AD = "GSE225158_reprocessed_pearson.h5ad"
SUMMARY_FILE = "reprocessing_summary.txt"

# Full output paths
OUTPUT_H5AD_PATH = f"{PROCESSED_DATA_DIR}/{OUTPUT_H5AD}"
SUMMARY_PATH = f"{PROCESSED_DATA_DIR}/{SUMMARY_FILE}"

# =======================================
# 📚 IMPORT LIBRARIES
# =======================================

# Import necessary libraries
import scanpy as sc
import scanpy.experimental as sce
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import sparse
import warnings
import os
warnings.filterwarnings('ignore')

# Set scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=300, facecolor='white')

def calculate_mad_outliers(values, n_mads=3):
    """Calculate outliers using Median Absolute Deviation (MAD)"""
    median = np.median(values)
    mad = np.median(np.abs(values - median))
    lower = median - n_mads * mad
    upper = median + n_mads * mad
    return lower, upper

def main():
    """Main processing pipeline"""
    
    # =======================================
    # 🔬 STEP 1: LOAD DATA & EXTRACT RAW COUNTS
    # =======================================
    print("=" * 70)
    print("🔬 STEP 1: DATA LOADING & RAW EXTRACTION")
    print("=" * 70)
    
    # Display paths
    print(f"📂 Input file: {INPUT_PATH}")
    print(f"📂 Output directory: {PROCESSED_DATA_DIR}")
    print(f"📂 Plots directory: {PLOTS_DIR}")
    
    # Load the data
    adata = sc.read_h5ad(INPUT_PATH)
    print(f"📊 Loaded dataset: {adata.shape}")
    
    # Show available clinical metadata
    clinical_cols = ['ID', 'Region', 'Case', 'Sex', 'Race', 'Age', 'BMI', 'PMI', 'pH', 'RIN', 
                     'Dx_OUD', 'Dx_Substances', 'Dx_Comorbid', 'celltype1', 'celltype2', 'celltype3']
    available_clinical = [col for col in clinical_cols if col in adata.obs.columns]
    print(f"📋 Available clinical metadata: {available_clinical}")
    
    # Extract raw counts
    print("\n🔍 Extracting raw count data...")
    if adata.raw is not None:
        raw_adata = adata.raw.to_adata()
        print(f"✅ Extracted raw data from .raw - Shape: {raw_adata.shape}")
        
        # Copy over the clinical metadata
        print("📋 Copying clinical metadata...")
        for col in available_clinical:
            if col in adata.obs.columns:
                raw_adata.obs[col] = adata.obs[col].copy()
        
        print(f"✅ Copied {len(available_clinical)} metadata columns")
    else:
        print("❌ No raw data found in .raw slot")
        return None
    
    # Clean up memory
    del adata
    
    # Validate raw counts
    print(f"\n🔍 Raw data validation:")
    print(f"   • Data type: {raw_adata.X.dtype}")
    print(f"   • Min value: {raw_adata.X.min()}")
    print(f"   • Max value: {raw_adata.X.max()}")
    print(f"   • Sparsity: {(raw_adata.X == 0).sum() / (raw_adata.X.shape[0] * raw_adata.X.shape[1]):.1%}")
    
    # =======================================
    # 🔬 STEP 2: QUALITY CONTROL
    # =======================================
    print("\n" + "=" * 70)
    print("🔬 STEP 2: QUALITY CONTROL")
    print("=" * 70)
    
    print(f"📊 Starting QC with {raw_adata.n_obs:,} cells and {raw_adata.n_vars:,} genes")
    
    # Calculate gene annotations
    print("🧬 Calculating gene annotations...")
    raw_adata.var['mt'] = raw_adata.var_names.str.startswith('MT-')
    raw_adata.var['ribo'] = raw_adata.var_names.str.startswith(('RPS', 'RPL'))
    raw_adata.var['hb'] = raw_adata.var_names.str.contains('^HB[^(P)]')
    
    # Calculate QC metrics with proper qc_vars specification
    print("📈 Calculating QC metrics...")
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
    
    # Print initial QC summary
    print("\n📋 Initial QC Summary:")
    print(f"   • Cells: {raw_adata.n_obs:,}")
    print(f"   • Genes: {raw_adata.n_vars:,}")
    print(f"   • MT genes: {raw_adata.var['mt'].sum()}")
    print(f"   • Ribosomal genes: {raw_adata.var['ribo'].sum()}")
    print(f"   • Median genes/cell: {raw_adata.obs['n_genes_by_counts'].median():.0f}")
    print(f"   • Median UMI/cell: {raw_adata.obs['total_counts'].median():.0f}")
    print(f"   • Median MT%: {raw_adata.obs['pct_counts_mt'].median():.2f}%")
    
    # Show QC by OUD diagnosis if available
    if 'Dx_OUD' in raw_adata.obs.columns:
        print(f"\n📊 QC by OUD diagnosis:")
        oud_summary = raw_adata.obs.groupby('Dx_OUD').agg({
            'total_counts': 'median',
            'n_genes_by_counts': 'median',
            'pct_counts_mt': 'median'
        }).round(1)
        print(oud_summary)
    
    # Robust QC filtering using MAD
    print("\n🎯 Robust QC Filtering (MAD-based)...")
    
    # Calculate thresholds using MAD
    n_genes_lower, n_genes_upper = calculate_mad_outliers(raw_adata.obs['n_genes_by_counts'], n_mads=2.5)
    total_counts_lower, total_counts_upper = calculate_mad_outliers(raw_adata.obs['total_counts'], n_mads=2.5)
    mt_lower, mt_upper = calculate_mad_outliers(raw_adata.obs['pct_counts_mt'], n_mads=2.5)
    
    print(f"   • Gene count thresholds: {n_genes_lower:.0f} - {n_genes_upper:.0f}")
    print(f"   • UMI count thresholds: {total_counts_lower:.0f} - {total_counts_upper:.0f}")
    print(f"   • MT% thresholds: {mt_lower:.2f}% - {mt_upper:.2f}%")
    
    # Apply conservative filters
    cell_filter = (
        (raw_adata.obs['n_genes_by_counts'] >= max(200, n_genes_lower)) &
        (raw_adata.obs['n_genes_by_counts'] <= n_genes_upper) &
        (raw_adata.obs['total_counts'] >= total_counts_lower) &
        (raw_adata.obs['total_counts'] <= total_counts_upper) &
        (raw_adata.obs['pct_counts_mt'] <= min(25, mt_upper))  # Cap MT at 25% for snRNA-seq
    )
    
    print(f"   • Cells passing QC: {cell_filter.sum():,} / {len(cell_filter):,} ({cell_filter.sum()/len(cell_filter)*100:.1f}%)")
    
    # Show filtering impact by group
    if 'Dx_OUD' in raw_adata.obs.columns:
        print(f"\n📊 QC impact by OUD diagnosis:")
        filter_impact = pd.crosstab(raw_adata.obs['Dx_OUD'], cell_filter, margins=True)
        print(filter_impact)
    
    # Apply cell filters
    raw_adata = raw_adata[cell_filter, :].copy()
    
    # Filter genes (expressed in at least 3 cells)
    sc.pp.filter_genes(raw_adata, min_cells=3)
    
    print(f"✅ After QC filtering: {raw_adata.n_obs:,} cells, {raw_adata.n_vars:,} genes")
    
    # =======================================
    # 🔁 STEP 3: PEARSON RESIDUALS PIPELINE
    # =======================================
    print("\n" + "=" * 70)
    print("🔁 STEP 3: PEARSON RESIDUALS NORMALIZATION")
    print("=" * 70)
    
    # Store raw counts
    raw_adata.raw = raw_adata
    print("💾 Saved raw counts to .raw")
    
    # Check for batch variables
    batch_key = None
    potential_batch_keys = ['ID', 'Region', 'Case']
    for key in potential_batch_keys:
        if key in raw_adata.obs.columns:
            n_batches = raw_adata.obs[key].nunique()
            if n_batches > 1 and n_batches < raw_adata.n_obs * 0.8:  # Reasonable batch size
                batch_key = key
                print(f"🔄 Using '{key}' as batch key ({n_batches} batches)")
                break
    
    if batch_key is None:
        print("🔄 No batch correction (no suitable batch variable found)")
    
    # Apply Pearson residuals recipe
    print("\n🔬 Applying Pearson residuals recipe...")
    print("   • Method: Full Pearson residuals pipeline (Lause et al., 2021)")
    print(f"   • Parameters: theta=100, n_top_genes=3000, n_comps=50")
    print(f"   • Batch key: {batch_key}")
    
    sce.pp.recipe_pearson_residuals(
        raw_adata,
        theta=100,
        clip=None,
        n_top_genes=3000,  # More genes for complex tissue
        batch_key=batch_key,
        n_comps=50,
        random_state=42,
        inplace=True
    )
    
    print("✅ Pearson residuals pipeline complete!")
    print(f"   • Final dataset: {raw_adata.n_obs:,} cells × {raw_adata.n_vars:,} genes")
    
    # =======================================
    # 🔍 STEP 4: NEIGHBORHOOD & EMBEDDING
    # =======================================
    print("\n" + "=" * 70)
    print("🔍 STEP 4: NEIGHBORHOOD GRAPH & EMBEDDING")
    print("=" * 70)
    
    # Compute neighborhood graph
    print("🕸️  Computing neighborhood graph...")
    sc.pp.neighbors(raw_adata, n_neighbors=15, n_pcs=40, random_state=42)
    
    # UMAP embedding
    print("🗺️  Computing UMAP embedding...")
    sc.tl.umap(raw_adata, random_state=42)
    
    print("✅ Neighborhood graph and UMAP complete")
    
    # =======================================
    # 🎯 STEP 5: CLUSTERING
    # =======================================
    print("\n" + "=" * 70)
    print("🎯 STEP 5: CLUSTERING")
    print("=" * 70)
    
    # Leiden clustering with multiple resolutions
    resolutions = [0.1, 0.3, 0.5, 0.7, 1.0, 1.5]
    print(f"🔍 Testing Leiden clustering at resolutions: {resolutions}")
    
    for res in resolutions:
        sc.tl.leiden(raw_adata, resolution=res, key_added=f'leiden_res_{res}', random_state=42)
        n_clusters = len(raw_adata.obs[f'leiden_res_{res}'].unique())
        print(f"   • Resolution {res}: {n_clusters} clusters")
    
    # Set default clustering
    raw_adata.obs['leiden'] = raw_adata.obs['leiden_res_0.5']
    print(f"✅ Using resolution 0.5 as default ({len(raw_adata.obs['leiden'].unique())} clusters)")
    
    # =======================================
    # 📊 STEP 6: VISUALIZATION & QC PLOTS
    # =======================================
    print("\n" + "=" * 70)
    print("📊 STEP 6: GENERATING VISUALIZATIONS")
    print("=" * 70)
    
    # Create output directory
    os.makedirs(PLOTS_DIR, exist_ok=True)
    
    # Main UMAP plots
    print("🎨 Creating UMAP visualizations...")
    
    # Figure 1: Basic clustering and QC
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    sc.pl.umap(raw_adata, color='leiden', ax=axes[0,0], show=False, legend_loc='on data', 
               legend_fontsize=8, size=2)
    axes[0,0].set_title('New Leiden Clusters (res=0.5)', fontsize=14)
    
    sc.pl.umap(raw_adata, color='pct_counts_mt', ax=axes[0,1], show=False, size=2)
    axes[0,1].set_title('Mitochondrial Gene %', fontsize=14)
    
    sc.pl.umap(raw_adata, color='total_counts', ax=axes[1,0], show=False, size=2)
    axes[1,0].set_title('Total UMI Counts', fontsize=14)
    
    sc.pl.umap(raw_adata, color='n_genes_by_counts', ax=axes[1,1], show=False, size=2)
    axes[1,1].set_title('Number of Genes', fontsize=14)
    
    plt.suptitle('Fresh Processing - Quality Control Overview', fontsize=16, y=0.98)
    plt.tight_layout()
    plt.savefig(f"{PLOTS_DIR}/umap_qc_overview.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Figure 2: Clinical metadata (if available)
    if 'Dx_OUD' in raw_adata.obs.columns:
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        sc.pl.umap(raw_adata, color='Dx_OUD', ax=axes[0,0], show=False, size=2)
        axes[0,0].set_title('OUD Diagnosis', fontsize=14)
        
        if 'Sex' in raw_adata.obs.columns:
            sc.pl.umap(raw_adata, color='Sex', ax=axes[0,1], show=False, size=2)
            axes[0,1].set_title('Sex', fontsize=14)
        
        if 'celltype1' in raw_adata.obs.columns:
            sc.pl.umap(raw_adata, color='celltype1', ax=axes[1,0], show=False, 
                      legend_loc='right margin', legend_fontsize=6, size=2)
            axes[1,0].set_title('Original Cell Types (Level 1)', fontsize=14)
        
        if 'Region' in raw_adata.obs.columns:
            sc.pl.umap(raw_adata, color='Region', ax=axes[1,1], show=False, size=2)
            axes[1,1].set_title('Brain Region', fontsize=14)
        
        plt.suptitle('Clinical Metadata Overview', fontsize=16, y=0.98)
        plt.tight_layout()
        plt.savefig(f"{PLOTS_DIR}/umap_clinical_overview.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    print(f"   • Saved UMAP plots to: {PLOTS_DIR}/")
    
    # =======================================
    # 💾 STEP 7: SAVE RESULTS
    # =======================================
    print("\n" + "=" * 70)
    print("💾 STEP 7: SAVING RESULTS")
    print("=" * 70)
    
    # Create output directory
    os.makedirs(PROCESSED_DATA_DIR, exist_ok=True)
    
    # Save processed data
    print(f"💾 Saving processed data to: {OUTPUT_H5AD_PATH}")
    raw_adata.write(OUTPUT_H5AD_PATH)
    
    # Save processing summary
    with open(SUMMARY_PATH, 'w') as f:
        f.write("GSE225158 OUD Striatum snRNA-seq Reprocessing Summary\n")
        f.write("=" * 60 + "\n")
        f.write(f"Processing date: {pd.Timestamp.now()}\n")
        f.write(f"Method: Fresh processing with Pearson residuals\n\n")
        f.write("PATHS:\n")
        f.write(f"• Input: {INPUT_PATH}\n")
        f.write(f"• Output: {OUTPUT_H5AD_PATH}\n")
        f.write(f"• Plots: {PLOTS_DIR}\n\n")
        f.write("DATASET OVERVIEW:\n")
        f.write(f"• Final cells: {raw_adata.n_obs:,}\n")
        f.write(f"• Final genes: {raw_adata.n_vars:,}\n")
        f.write(f"• New clusters: {len(raw_adata.obs['leiden'].unique())}\n")
        f.write(f"• Normalization: Pearson residuals (theta=100)\n")
        f.write(f"• Batch correction: {batch_key if batch_key else 'None'}\n\n")
        f.write("CLINICAL METADATA PRESERVED:\n")
        for col in available_clinical:
            if col in raw_adata.obs.columns:
                f.write(f"• {col}: {raw_adata.obs[col].nunique()} unique values\n")
    
    print(f"💾 Saved summary to: {SUMMARY_PATH}")
    
    # Final summary
    print("\n" + "=" * 70)
    print("✅ REPROCESSING COMPLETE!")
    print("=" * 70)
    print(f"📊 Dataset: {raw_adata.n_obs:,} cells × {raw_adata.n_vars:,} genes")
    print(f"🎯 New clusters: {len(raw_adata.obs['leiden'].unique())}")
    print(f"🔬 Method: Pearson residuals normalization")
    print(f"📋 Clinical metadata: {len(available_clinical)} variables preserved")
    print(f"💾 Data saved to: {OUTPUT_H5AD_PATH}")
    print(f"📊 Plots saved to: {PLOTS_DIR}/")
    print("=" * 70)
    
    return raw_adata

if __name__ == "__main__":
    # Run the reprocessing pipeline
    reprocessed_adata = main()