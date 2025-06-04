'''
🧠 Cell Type Annotation Pipeline
GSE225158 - OUD Striatum snRNA-seq cell type identification

Workflow:
1. Load preprocessed Pearson residuals data (01 output)
2. Marker-based cell type prediction
3. Manual annotation refinement
4. Cell type validation and QC
5. Save annotated results

Note: For scVI data, use 02b_cell_type_annotation_scvi.py
'''

# =======================================
# 📁 INPUT/OUTPUT PATHS CONFIGURATION
# =======================================

# Input from preprocessing (01 script only)
BASE_OUTPUT_DIR = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study"
PROCESSED_DATA_DIR_QC = f"{BASE_OUTPUT_DIR}/data/processed/snrna_qc"
RESULTS_DIR = f"{BASE_OUTPUT_DIR}/results/snrna"
PLOTS_DIR = f"{RESULTS_DIR}/annotation_plots"

# Input files (Pearson residuals only)
INPUT_H5AD = f"{PROCESSED_DATA_DIR_QC}/GSE225158_reprocessed_pearson.h5ad"

# Output files
OUTPUT_H5AD = f"{PROCESSED_DATA_DIR_QC}/GSE225158_annotated_pearson.h5ad"
ANNOTATION_SUMMARY = f"{PROCESSED_DATA_DIR_QC}/cell_type_annotation_summary.txt"

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
    """Cell type annotation pipeline for Pearson residuals data"""
    
    print("=" * 70)
    print("🧠 CELL TYPE ANNOTATION PIPELINE (PEARSON RESIDUALS)")
    print("=" * 70)
    
    # =======================================
    # 📂 STEP 1: LOAD PREPROCESSED DATA
    # =======================================
    print("📂 STEP 1: LOADING PREPROCESSED DATA")
    print("=" * 70)
    
    # Load Pearson residuals data
    if os.path.exists(INPUT_H5AD):
        adata = sc.read_h5ad(INPUT_H5AD)
        print(f"📊 Loaded Pearson residuals data: {adata.shape}")
        
        # Handle gene names stored in 'features' column
        if 'features' in adata.var.columns:
            print("🔧 Setting gene names from 'features' column as var index...")
            # Set the features column as the gene names index
            adata.var_names = adata.var['features'].values
            adata.var_names_make_unique()  # Ensure unique gene names - FIXED METHOD NAME
            print(f"✅ Gene names set from features column")
        
        # Verify gene names are now accessible
        print(f"📋 Sample gene names: {adata.var_names[:5].tolist()}")
        
    else:
        raise FileNotFoundError(f"Preprocessed data not found! Run 01_preprocessing_qc.py first.\nExpected file: {INPUT_H5AD}")
    
    print(f"✅ Using Pearson residuals processed data")

    # Remove the marker genes loading section since MARKER_GENES_CSV is not defined
    # This was left over from the original script
    print("⚠️  Will compute fresh marker genes from defined brain markers")

    # =======================================
    # 🧬 STEP 2: CELL TYPE SIGNATURE SCORING
    # =======================================
    print("\n🧬 STEP 2: CELL TYPE SIGNATURE SCORING")
    print("=" * 70)
    
    # Known brain cell type markers (ENHANCED STRIATUM-SPECIFIC)
    brain_markers = {
        # STRIATAL PROJECTION NEURONS (MSNs)
        'D1_MSN': ['DRD1', 'TAC1', 'PDYN', 'PPP1R1B', 'FOXP1', 'RGS9', 'CALB1'],
        'D2_MSN': ['DRD2', 'PENK', 'ADORA2A', 'PPP1R1B', 'FOXP2', 'GPR6', 'EPHA4'],

        # STRIATAL INTERNEURONS
        'Cholinergic_Interneurons': ['CHAT', 'SLC18A3', 'SLC5A7', 'ISL1', 'LHX8', 'TNFAIP8L3'],
        'PV_Interneurons': ['PVALB', 'TAC1', 'ERBB4', 'LHX6', 'GAD1'],
        'SST_Interneurons': ['SST', 'NPY', 'RELN', 'NOS1', 'CHODL'],
        'Calretinin_Interneurons': ['CALB2', 'VIP', 'LAMP5', 'NDNF', 'NTNG1'],

        # OTHER NEURONS
        'GABA_Interneurons': ['GAD1', 'GAD2', 'SLC32A1', 'LHX6', 'LAMP5'],
        'Cortical_Inputs': ['SLC17A7', 'CAMK2A', 'TBR1', 'SATB2'],  # Renamed from 'Excitatory_Inputs'

        # GLIA
        'Astrocytes': ['GFAP', 'AQP4', 'S100B', 'ALDH1L1', 'SLC1A2', 'CLU', 'GJA1'],
        'Oligodendrocytes': ['PLP1', 'MOG', 'MBP', 'CNP', 'MAG', 'OPALIN'],
        'OPCs': ['PDGFRA', 'CSPG4', 'OLIG1', 'OLIG2', 'SOX10', 'GPR17', 'TNR'],
        'Microglia': ['CX3CR1', 'TMEM119', 'P2RY12', 'AIF1', 'TREM2', 'SALL1', 'CSF1R'],

        # VASCULAR
        'Endothelial': ['CLDN5', 'FLT1', 'VWF', 'PECAM1', 'CDH5', 'TIE1'],
        'Pericytes': ['PDGFRB', 'RGS5', 'ACTA2', 'MCAM', 'NOTCH3', 'ABCC9'],
        'VLMC': ['COL1A2', 'LUM', 'DCN', 'EMCN', 'MRC1'],

        # RARE CELLS
        'Ependymal': ['FOXJ1', 'RSPH1', 'PIFO', 'CCDC153'],
        'Tanycytes': ['RAX', 'FXYD1', 'DIO2', 'LHX1'],
        
        # OPTIONAL: Add MSN subtypes if you want more granularity
        'Patch_MSN': ['OPRM1', 'MOR1', 'PDYN'],  # Patch compartment
        'Matrix_MSN': ['CALB1', 'TAC1'],         # Matrix compartment
    }
    
    # Calculate signature scores with improved method
    print("🔍 Computing cell type signature scores...")
    signature_scores = {}
    
    for cell_type, markers in brain_markers.items():
        # Now check against the properly set var_names (from features column)
        available_markers = [m for m in markers if m in adata.var_names]
        
        if len(available_markers) >= 2:  # Need at least 2 markers
            sc.tl.score_genes(
                adata, 
                gene_list=available_markers, 
                score_name=f'{cell_type}_score',
                use_raw=False,
                ctrl_size=min(50, len(available_markers)*5),
                n_bins=25
            )
            signature_scores[cell_type] = available_markers
            print(f"   • {cell_type}: {len(available_markers)}/{len(markers)} markers")
            if len(available_markers) < len(markers):
                missing_markers = [m for m in markers if m not in adata.var_names]
                print(f"     Missing: {missing_markers[:3]}{'...' if len(missing_markers) > 3 else ''}")
        else:
            print(f"   ⚠️  {cell_type}: Insufficient markers ({len(available_markers)}/{len(markers)})")
            missing_markers = [m for m in markers if m not in adata.var_names]
            print(f"     Missing: {missing_markers}")

    # =======================================
    # 🎯 STEP 3: CELL TYPE PREDICTION (MARKER-BASED ONLY)
    # =======================================
    print("\n🎯 STEP 3: CELL TYPE PREDICTION")
    print("=" * 70)
    
    # Get all score columns
    score_columns = [col for col in adata.obs.columns if col.endswith('_score')]
    print(f"🔢 Using {len(score_columns)} signature scores for prediction")
    
    # IMPROVED PREDICTION WITH THRESHOLDING (remove reference mapping code)
    if len(score_columns) > 0:
        # Calculate prediction with confidence thresholds
        score_matrix = adata.obs[score_columns].values
        predicted_types = []
        prediction_confidence = []
        prediction_scores = []
        
        # Calculate thresholds for each cell type
        score_thresholds = {}
        for i, col in enumerate(score_columns):
            scores = score_matrix[:, i]
            threshold = np.percentile(scores, 75)  # 75th percentile as threshold
            score_thresholds[col] = threshold
        
        for i in range(len(score_matrix)):
            scores = score_matrix[i]
            
            # Find scores above threshold
            above_threshold = []
            for j, score in enumerate(scores):
                if score > score_thresholds[score_columns[j]]:
                    above_threshold.append((j, score))
            
            if above_threshold:
                # Choose highest score among those above threshold
                max_idx, max_score = max(above_threshold, key=lambda x: x[1])
                confidence = max_score - np.mean(scores)
            else:
                # No strong prediction - assign based on highest score but low confidence
                max_idx = np.argmax(scores)
                max_score = scores[max_idx]
                confidence = 0.1  # Low confidence
            
            predicted_type = score_columns[max_idx].replace('_score', '')
            predicted_types.append(predicted_type)
            prediction_confidence.append(confidence)
            prediction_scores.append(max_score)
        
        adata.obs['predicted_cell_type'] = predicted_types
        adata.obs['prediction_confidence'] = prediction_confidence
        adata.obs['prediction_score'] = prediction_scores
        
        # Summary with confidence filtering
        print("📊 Initial cell type predictions:")
        pred_counts = adata.obs['predicted_cell_type'].value_counts()
        high_conf_mask = adata.obs['prediction_confidence'] > np.percentile(adata.obs['prediction_confidence'], 50)
        
        for cell_type, count in pred_counts.items():
            high_conf_count = ((adata.obs['predicted_cell_type'] == cell_type) & high_conf_mask).sum()
            pct = count / len(adata.obs) * 100
            pct_high_conf = high_conf_count / count * 100 if count > 0 else 0
            print(f"   • {cell_type}: {count:,} cells ({pct:.1f}%) - {pct_high_conf:.1f}% high confidence")

    # =======================================
    # 🔍 STEP 4: ANNOTATION REFINEMENT
    # =======================================
    print("\n🔍 STEP 4: ANNOTATION REFINEMENT")
    print("=" * 70)
    
    # ENHANCED REFINEMENT WITH MARKER EXPRESSION VALIDATION
    print("🔧 Refining annotations using cluster information and marker validation...")
    
    cluster_annotations = {}
    refined_annotations = []
    
    for cluster in adata.obs['leiden'].unique():
        cluster_mask = adata.obs['leiden'] == cluster
        cluster_predictions = adata.obs.loc[cluster_mask, 'predicted_cell_type']
        
        # Get most common prediction for this cluster
        most_common = cluster_predictions.value_counts().index[0]
        
        # Validate with marker expression if available
        if most_common in signature_scores:
            cluster_cells = adata[cluster_mask]
            marker_genes = signature_scores[most_common]
            
            # Check if cluster expresses the markers
            marker_expression = cluster_cells[:, marker_genes].X.mean(axis=0)
            mean_marker_expr = np.mean(marker_expression)
            
            # If low marker expression, mark as uncertain
            if mean_marker_expr < np.percentile(adata.X.toarray(), 25):
                most_common = f"Uncertain_{most_common}"
        
        cluster_annotations[cluster] = most_common
        
        # Count cells that agree with cluster consensus
        agreement = (cluster_predictions == most_common.replace("Uncertain_", "")).sum()
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
    
    # Main annotation plot (simplified for marker-based only)
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
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
    elif 'Dx_OUD' in adata.obs.columns:
        sc.pl.umap(adata, color='Dx_OUD', ax=axes[1,1], show=False, size=3)
        axes[1,1].set_title('OUD Diagnosis', fontsize=16)
    
    plt.suptitle('Cell Type Annotation Results (Pearson Residuals)', fontsize=20, y=0.98)
    plt.tight_layout()
    plt.savefig(f"{PLOTS_DIR}/cell_type_annotation_pearson.png", dpi=300, bbox_inches='tight')
    plt.close()

    print(f"✅ Annotation plots saved to: {PLOTS_DIR}/")

    # Save annotation summary (simplified)
    with open(ANNOTATION_SUMMARY, 'w') as f:
        f.write("GSE225158 Cell Type Annotation Summary (Pearson Residuals)\n")
        f.write("=" * 60 + "\n")
        f.write(f"Date: {pd.Timestamp.now()}\n")
        f.write(f"Data source: Pearson residuals normalization\n\n")
        
        f.write("ANNOTATION METHODS USED:\n")
        f.write("• Marker-based signature scoring\n")
        f.write("• Cluster-based refinement\n\n")
        
        f.write("FINAL CELL TYPE COUNTS:\n")
        for cell_type, count in final_counts.items():
            pct = count / len(adata.obs) * 100
            f.write(f"• {cell_type}: {count:,} cells ({pct:.1f}%)\n")
        
        f.write(f"\nTOTAL CELLS: {len(adata.obs):,}\n")
        f.write(f"CELL TYPES: {len(final_counts)}\n\n")
        
        f.write("MARKER GENES USED:\n")
        for cell_type, markers in signature_scores.items():
            f.write(f"• {cell_type}: {', '.join(markers)}\n")

    print(f"💾 Summary saved to: {ANNOTATION_SUMMARY}")
    
    # Fix the region analysis path
    if 'Region' in adata.obs.columns:
        print("\n🧠 Brain region-specific cell type analysis:")
        region_cell_type_counts = pd.crosstab(adata.obs['Region'], adata.obs['refined_cell_type'])
        print(region_cell_type_counts)
        
        # Save region-specific analysis
        region_analysis_path = f"{PROCESSED_DATA_DIR_QC}/region_celltype_analysis.csv"
        region_cell_type_counts.to_csv(region_analysis_path)
        print(f"💾 Region analysis saved to: {region_analysis_path}")
    
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

