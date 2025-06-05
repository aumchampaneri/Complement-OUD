'''
🧠 Cell Type Annotation Pipeline - Clean & Simple
GSE225158 - OUD Striatum snRNA-seq using validated original annotations

Strategy: Use original author annotations (celltype3) as ground truth
Rationale: Original annotations are well-balanced and biologically validated

OUTPUT RECOMMENDATION:
Use scVI outputs for downstream analysis because:
✅ Batch correction - Removes technical batch effects
✅ Denoising - Reduces technical noise while preserving biological signal  
✅ Uncertainty quantification - Provides confidence estimates
✅ Better DE detection - More sensitive differential expression analysis
✅ Latent space - Clean embedding for visualization and clustering
'''

# =======================================
# 📁 INPUT/OUTPUT PATHS CONFIGURATION
# =======================================

# Base directory for data and results
BASE_DIR = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study"

# Input from scVI preprocessing (01b script)
INPUT_H5AD = f"{BASE_DIR}/data/processed/snrna_scvi/GSE225158_reprocessed_scvi.h5ad"

# Output files
OUTPUT_H5AD = f"{BASE_DIR}/data/processed/snrna_scvi/GSE225158_annotated_scvi.h5ad"
PLOTS_DIR = f"{BASE_DIR}/results/snrna_scvi/annotation_plots"
SUMMARY_FILE = f"{BASE_DIR}/data/processed/snrna_scvi/cell_type_annotation_summary.txt"

# =======================================
# 📚 IMPORT LIBRARIES
# =======================================

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

# Configure scanpy
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=300, facecolor='white')
sc.settings.autoshow = False

def create_clean_annotations(adata):
    """Map original celltype3 to standardized cell type names"""
    print("🎯 Creating clean cell type annotations...")
    
    if 'celltype3' not in adata.obs.columns:
        raise ValueError("celltype3 column not found in data!")
    
    # Standardized mapping
    mapping = {
        'Oligos': 'Oligodendrocytes',
        'Oligos_Pre': 'OPCs',
        'Microglia': 'Microglia', 
        'Astrocytes': 'Astrocytes',
        'Endothelial': 'Endothelial',
        'Mural': 'Pericytes',
        'D1-Matrix': 'D1_MSN',
        'D1-Striosome': 'D1_MSN',
        'D2-Matrix': 'D2_MSN', 
        'D2-Striosome': 'D2_MSN',
        'D1/D2-Hybrid': 'Mixed_MSN',
        'Int-TH': 'Cholinergic_Interneurons',
        'Int-CCK': 'CCK_Interneurons',
        'Int-PTHLH': 'PTHLH_Interneurons',
        'Int-SST': 'SST_Interneurons'
    }
    
    # Apply mapping
    adata.obs['cell_type'] = adata.obs['celltype3'].map(mapping).fillna(adata.obs['celltype3'])
    adata.obs['annotation_method'] = 'Original_Authors'
    
    # Show results
    counts = adata.obs['cell_type'].value_counts()
    print("   Cell type distribution:")
    for cell_type, count in counts.items():
        pct = count / len(adata) * 100
        print(f"   • {cell_type}: {count:,} cells ({pct:.1f}%)")
    
    return adata

def validate_markers(adata):
    """Quick validation using key marker genes"""
    print("\n🧬 Validating with marker genes...")
    
    markers = {
        'Oligodendrocytes': ['MBP', 'MOG', 'PLP1'],
        'Astrocytes': ['GFAP', 'AQP4', 'S100B'],
        'Microglia': ['AIF1', 'CX3CR1', 'P2RY12']
    }
    
    expr_data = adata.layers.get('scvi_normalized', adata.X)
    
    for cell_type, marker_list in markers.items():
        available = [m for m in marker_list if m in adata.var_names]
        if available and cell_type in adata.obs['cell_type'].values:
            cell_mask = adata.obs['cell_type'] == cell_type
            
            # Calculate enrichment
            indices = [adata.var_names.get_loc(m) for m in available]
            if hasattr(expr_data, 'toarray'):
                marker_expr = expr_data[:, indices].toarray()
            else:
                marker_expr = expr_data[:, indices]
            
            in_type = np.mean(marker_expr[cell_mask])
            out_type = np.mean(marker_expr[~cell_mask])
            enrichment = in_type / (out_type + 1e-8)
            
            status = "✅" if enrichment > 1.5 else "⚠️"
            print(f"   {status} {cell_type}: {enrichment:.1f}x enriched")

def create_plots(adata):
    """Create final annotation plots"""
    print("\n📊 Creating plots...")
    
    os.makedirs(PLOTS_DIR, exist_ok=True)
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Final annotations
    sc.pl.umap(adata, color='cell_type', ax=axes[0,0], show=False, 
               legend_loc='right margin', legend_fontsize=8)
    axes[0,0].set_title('Final Cell Types', fontweight='bold')
    
    # Original annotations
    sc.pl.umap(adata, color='celltype3', ax=axes[0,1], show=False,
               legend_loc='right margin', legend_fontsize=6)
    axes[0,1].set_title('Original celltype3')
    
    # Clinical variables
    if 'Region' in adata.obs.columns:
        sc.pl.umap(adata, color='Region', ax=axes[1,0], show=False)
        axes[1,0].set_title('Brain Region')
    
    if 'Dx_OUD' in adata.obs.columns:
        sc.pl.umap(adata, color='Dx_OUD', ax=axes[1,1], show=False)
        axes[1,1].set_title('OUD Diagnosis')
    
    plt.suptitle('GSE225158 Cell Type Annotation', fontsize=16, y=0.98)
    plt.tight_layout()
    plt.savefig(f"{PLOTS_DIR}/cell_type_annotation.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"   Plots saved to: {PLOTS_DIR}/")

def save_results(adata):
    """Save annotated data and summary"""
    print("\n💾 Saving results...")
    
    # Save data
    adata.write(OUTPUT_H5AD)
    print(f"   Data: {OUTPUT_H5AD}")
    
    # Create summary
    counts = adata.obs['cell_type'].value_counts()
    
    with open(SUMMARY_FILE, 'w') as f:
        f.write("GSE225158 Cell Type Annotation Summary\n")
        f.write("=" * 50 + "\n")
        f.write(f"Date: {pd.Timestamp.now()}\n")
        f.write(f"Method: Original author annotations\n")
        f.write(f"Total cells: {len(adata):,}\n")
        f.write(f"Cell types: {len(counts)}\n\n")
        
        f.write("RECOMMENDED FOR DOWNSTREAM ANALYSIS:\n")
        f.write("Use scVI outputs (this file) because:\n")
        f.write("• Batch corrected and denoised\n")
        f.write("• Better for differential expression\n") 
        f.write("• Includes uncertainty quantification\n")
        f.write("• Clean latent space representation\n\n")
        
        f.write("Cell Type Counts:\n")
        f.write("-" * 30 + "\n")
        for cell_type, count in counts.items():
            pct = count / len(adata) * 100
            f.write(f"{cell_type:25s}: {count:6,} ({pct:5.1f}%)\n")
        
        # Add data layer information
        f.write(f"\nDATA LAYERS AVAILABLE:\n")
        for layer in adata.layers.keys():
            f.write(f"• {layer}: {adata.layers[layer].shape}\n")
        f.write(f"• X (main): {adata.X.shape}\n")
    
    print(f"   Summary: {SUMMARY_FILE}")
    
    # Save cross-tabs if clinical data available
    if 'Region' in adata.obs.columns:
        region_file = f"{BASE_DIR}/data/processed/snrna_scvi/region_by_celltype.csv"
        pd.crosstab(adata.obs['Region'], adata.obs['cell_type']).to_csv(region_file)
        print(f"   Region analysis: {region_file}")
    
    if 'Dx_OUD' in adata.obs.columns:
        oud_file = f"{BASE_DIR}/data/processed/snrna_scvi/oud_by_celltype.csv"
        pd.crosstab(adata.obs['Dx_OUD'], adata.obs['cell_type']).to_csv(oud_file)
        print(f"   OUD analysis: {oud_file}")

def main():
    """Main annotation workflow"""
    print("🧠 CELL TYPE ANNOTATION - GSE225158")
    print("=" * 50)
    
    # Load data
    print("📂 Loading data...")
    if not os.path.exists(INPUT_H5AD):
        raise FileNotFoundError(f"Input file not found: {INPUT_H5AD}")
    
    adata = sc.read_h5ad(INPUT_H5AD)
    print(f"   Loaded: {adata.shape}")
    
    # Process
    adata = create_clean_annotations(adata)
    validate_markers(adata)
    create_plots(adata)
    save_results(adata)
    
    print(f"\n✅ ANNOTATION COMPLETE!")
    print(f"📊 {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    print(f"🧠 {len(adata.obs['cell_type'].unique())} cell types")
    print("🚀 Ready for downstream analysis!")
    
    return adata

if __name__ == "__main__":
    annotated_adata = main()
