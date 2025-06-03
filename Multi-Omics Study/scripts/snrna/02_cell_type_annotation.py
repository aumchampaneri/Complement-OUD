'''
🧠 Cell Type Annotation Pipeline
GSE225158 - OUD Striatum snRNA-seq cell type identification

Workflow:
1. Load preprocessed scVI data
2. Marker-based cell type prediction
3. Manual annotation refinement
4. Cell type validation and QC
5. Save annotated results
'''

# =======================================
# 📁 INPUT/OUTPUT PATHS CONFIGURATION
# =======================================

# Input from preprocessing
BASE_OUTPUT_DIR = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study"
PROCESSED_DATA_DIR = f"{BASE_OUTPUT_DIR}/data/processed/snrna"
RESULTS_DIR = f"{BASE_OUTPUT_DIR}/results/snrna"
PLOTS_DIR = f"{RESULTS_DIR}/annotation_plots"

# Input files
INPUT_H5AD = f"{PROCESSED_DATA_DIR}/GSE225158_reprocessed_scvi.h5ad"
MARKER_GENES_CSV = f"{PROCESSED_DATA_DIR}/scvi_marker_genes.csv"

# Output files
OUTPUT_H5AD = f"{PROCESSED_DATA_DIR}/GSE225158_annotated.h5ad"
ANNOTATION_SUMMARY = f"{PROCESSED_DATA_DIR}/cell_type_annotation_summary.txt"

# =======================================
# 📚 IMPORT LIBRARIES
# =======================================

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import os
warnings.filterwarnings('ignore')

# Set scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=300, facecolor='white')

def main_annotation():
    """Cell type annotation pipeline"""
    
    print("=" * 70)
    print("🧠 CELL TYPE ANNOTATION PIPELINE")
    print("=" * 70)
    
    # =======================================
    # 📂 STEP 1: LOAD PREPROCESSED DATA
    # =======================================
    print("📂 STEP 1: LOADING PREPROCESSED DATA")
    print("=" * 70)
    
    # Load scVI-processed data
    adata = sc.read_h5ad(INPUT_H5AD)
    print(f"📊 Loaded preprocessed data: {adata.shape}")
    
    # Load marker genes
    if os.path.exists(MARKER_GENES_CSV):
        marker_df = pd.read_csv(MARKER_GENES_CSV)
        print(f"📋 Loaded {len(marker_df)} marker genes")
    else:
        print("⚠️  Marker genes file not found, will compute fresh")
        marker_df = None
    
    # =======================================
    # 🧬 STEP 2: CELL TYPE SIGNATURE SCORING
    # =======================================
    print("\n🧬 STEP 2: CELL TYPE SIGNATURE SCORING")
    print("=" * 70)
    
    # Known brain cell type markers (striatum-focused)
    brain_markers = {
        'Excitatory_Neurons': ['SLC17A7', 'CAMK2A', 'GRIN1', 'GRIA1', 'NEUROD6'],
        'Inhibitory_Neurons': ['GAD1', 'GAD2', 'SLC32A1', 'PVALB', 'SST', 'VIP', 'LAMP5'],
        'Medium_Spiny_Neurons_D1': ['DRD1', 'FOXP1', 'PPP1R1B', 'DARPP32', 'PDYN', 'BCAM'],
        'Medium_Spiny_Neurons_D2': ['DRD2', 'FOXP2', 'PPP1R1B', 'DARPP32', 'PENK', 'GPR6'],
        'Cholinergic_Interneurons': ['CHAT', 'SLC5A7', 'SLC18A3', 'ISL1'],
        'Parvalbumin_Interneurons': ['PVALB', 'GPR149', 'ERBB4'],
        'Somatostatin_Interneurons': ['SST', 'CHODL', 'MEIS2'],
        'Astrocytes': ['GFAP', 'AQP4', 'S100B', 'ALDH1L1', 'SLC1A2', 'APOE'],
        'Oligodendrocytes': ['MBP', 'MOG', 'PLP1', 'CNP', 'MOBP', 'OPALIN'],
        'OPCs': ['PDGFRA', 'CSPG4', 'SOX10', 'OLIG2', 'GPR17'],
        'Microglia': ['AIF1', 'CX3CR1', 'P2RY12', 'TMEM119', 'HEXB', 'CSF1R'],
        'Endothelial': ['PECAM1', 'VWF', 'FLT1', 'CDH5', 'CLDN5'],
        'Pericytes': ['PDGFRB', 'RGS5', 'ACTA2', 'MCAM', 'ABCC9'],
        'VLMC': ['COLEC12', 'DCN', 'GSN', 'PDGFRA'],  # Vascular Leptomeningeal Cells
    }
    
    # Calculate signature scores
    print("🔍 Computing cell type signature scores...")
    signature_scores = {}
    
    for cell_type, markers in brain_markers.items():
        available_markers = [m for m in markers if m in adata.var_names]
        
        if len(available_markers) >= 2:  # Need at least 2 markers
            sc.tl.score_genes(
                adata, 
                gene_list=available_markers, 
                score_name=f'{cell_type}_score',
                use_raw=False
            )
            signature_scores[cell_type] = available_markers
            print(f"   • {cell_type}: {len(available_markers)}/{len(markers)} markers")
        else:
            print(f"   • {cell_type}: Insufficient markers ({len(available_markers)}/{len(markers)})")
    
    # =======================================
    # 🎯 STEP 3: CELL TYPE PREDICTION
    # =======================================
    print("\n🎯 STEP 3: CELL TYPE PREDICTION")
    print("=" * 70)
    
    # Get all score columns
    score_columns = [col for col in adata.obs.columns if col.endswith('_score')]
    print(f"🔢 Using {len(score_columns)} signature scores for prediction")
    
    if len(score_columns) > 0:
        # Predict cell types based on highest score
        score_matrix = adata.obs[score_columns].values
        predicted_types = []
        prediction_confidence = []
        prediction_scores = []
        
        for i in range(len(score_matrix)):
            scores = score_matrix[i]
            max_idx = np.argmax(scores)
            max_score = scores[max_idx]
            
            # Calculate confidence as difference from mean
            mean_score = np.mean(scores)
            confidence = max_score - mean_score
            
            predicted_type = score_columns[max_idx].replace('_score', '')
            predicted_types.append(predicted_type)
            prediction_confidence.append(confidence)
            prediction_scores.append(max_score)
        
        adata.obs['predicted_cell_type'] = predicted_types
        adata.obs['prediction_confidence'] = prediction_confidence
        adata.obs['prediction_score'] = prediction_scores
        
        # Summary
        print("📊 Initial cell type predictions:")
        pred_counts = adata.obs['predicted_cell_type'].value_counts()
        for cell_type, count in pred_counts.items():
            pct = count / len(adata.obs) * 100
            print(f"   • {cell_type}: {count:,} cells ({pct:.1f}%)")
    
    # =======================================
    # 🔍 STEP 4: ANNOTATION REFINEMENT
    # =======================================
    print("\n🔍 STEP 4: ANNOTATION REFINEMENT")
    print("=" * 70)
    
    # Refine annotations based on cluster consistency
    print("🔧 Refining annotations using cluster information...")
    
    # For each cluster, assign the most common predicted type
    cluster_annotations = {}
    refined_annotations = []
    
    for cluster in adata.obs['leiden'].unique():
        cluster_mask = adata.obs['leiden'] == cluster
        cluster_predictions = adata.obs.loc[cluster_mask, 'predicted_cell_type']
        
        # Get most common prediction for this cluster
        most_common = cluster_predictions.value_counts().index[0]
        cluster_annotations[cluster] = most_common
        
        # Count cells that agree with cluster consensus
        agreement = (cluster_predictions == most_common).sum()
        total = len(cluster_predictions)
        print(f"   • Cluster {cluster}: {most_common} ({agreement}/{total} = {agreement/total*100:.1f}% agreement)")
    
    # Apply cluster-based refinement
    for i, cluster in enumerate(adata.obs['leiden']):
        refined_annotations.append(cluster_annotations[cluster])
    
    adata.obs['refined_cell_type'] = refined_annotations
    
    # Calculate refinement confidence
    refinement_confidence = []
    for i, (original, refined) in enumerate(zip(adata.obs['predicted_cell_type'], refined_annotations)):
        if original == refined:
            confidence = adata.obs.iloc[i]['prediction_confidence']
        else:
            confidence = adata.obs.iloc[i]['prediction_confidence'] * 0.5  # Reduce confidence for changed annotations
        refinement_confidence.append(confidence)
    
    adata.obs['final_confidence'] = refinement_confidence
    
    # Final summary
    print("\n📊 Final cell type annotations:")
    final_counts = adata.obs['refined_cell_type'].value_counts()
    for cell_type, count in final_counts.items():
        pct = count / len(adata.obs) * 100
        print(f"   • {cell_type}: {count:,} cells ({pct:.1f}%)")
    
    # =======================================
    # 📊 STEP 5: VISUALIZATION
    # =======================================
    print("\n📊 STEP 5: GENERATING VISUALIZATIONS")
    print("=" * 70)
    
    os.makedirs(PLOTS_DIR, exist_ok=True)
    
    # Main annotation plot
    fig, axes = plt.subplots(2, 2, figsize=(20, 16))
    
    # Refined cell types
    sc.pl.umap(adata, color='refined_cell_type', ax=axes[0,0], show=False, 
               legend_loc='right margin', legend_fontsize=8, size=3)
    axes[0,0].set_title('Final Cell Type Annotations', fontsize=16)
    
    # Clusters
    sc.pl.umap(adata, color='leiden', ax=axes[0,1], show=False, 
               legend_loc='on data', legend_fontsize=10, size=3)
    axes[0,1].set_title('Leiden Clusters', fontsize=16)
    
    # Prediction confidence
    sc.pl.umap(adata, color='final_confidence', ax=axes[1,0], show=False, size=3)
    axes[1,0].set_title('Annotation Confidence', fontsize=16)
    
    # Compare with original if available
    if 'celltype1' in adata.obs.columns:
        sc.pl.umap(adata, color='celltype1', ax=axes[1,1], show=False, 
                   legend_loc='right margin', legend_fontsize=6, size=3)
        axes[1,1].set_title('Original Annotations', fontsize=16)
    else:
        # OUD diagnosis
        if 'Dx_OUD' in adata.obs.columns:
            sc.pl.umap(adata, color='Dx_OUD', ax=axes[1,1], show=False, size=3)
            axes[1,1].set_title('OUD Diagnosis', fontsize=16)
    
    plt.suptitle('Cell Type Annotation Results', fontsize=20, y=0.98)
    plt.tight_layout()
    plt.savefig(f"{PLOTS_DIR}/cell_type_annotation_overview.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✅ Annotation plots saved to: {PLOTS_DIR}/")
    
    # =======================================
    # 💾 STEP 6: SAVE RESULTS
    # =======================================
    print("\n💾 STEP 6: SAVING RESULTS")
    print("=" * 70)
    
    # Set final annotation as default
    adata.obs['cell_type'] = adata.obs['refined_cell_type']
    
    # Save annotated data
    adata.write(OUTPUT_H5AD)
    print(f"💾 Annotated data saved to: {OUTPUT_H5AD}")
    
    # Save annotation summary
    with open(ANNOTATION_SUMMARY, 'w') as f:
        f.write("GSE225158 Cell Type Annotation Summary\n")
        f.write("=" * 50 + "\n")
        f.write(f"Date: {pd.Timestamp.now()}\n\n")
        
        f.write("CELL TYPE COUNTS:\n")
        for cell_type, count in final_counts.items():
            pct = count / len(adata.obs) * 100
            f.write(f"• {cell_type}: {count:,} cells ({pct:.1f}%)\n")
        
        f.write(f"\nTOTAL CELLS: {len(adata.obs):,}\n")
        f.write(f"CELL TYPES: {len(final_counts)}\n")
        
        f.write("\nMARKER GENES USED:\n")
        for cell_type, markers in signature_scores.items():
            f.write(f"• {cell_type}: {', '.join(markers)}\n")
    
    print(f"💾 Summary saved to: {ANNOTATION_SUMMARY}")
    
    print("\n" + "=" * 70)
    print("✅ CELL TYPE ANNOTATION COMPLETE!")
    print("=" * 70)
    print(f"📊 Dataset: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    print(f"🧠 Cell types: {len(final_counts)}")
    print(f"📁 Ready for downstream analysis!")
    print("=" * 70)
    
    return adata

if __name__ == "__main__":
    annotated_adata = main_annotation()

