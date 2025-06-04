'''
🧠 Cell Type Annotation Pipeline - scVI + Multi-Method Validation
GSE225158 - OUD Striatum snRNA-seq cell type identification using multiple methods

Workflow:
1. Load scVI-processed data
2. scVI-optimized marker scoring
3. Reference-based annotation (if available)
4. Cluster-based consensus
5. Multi-method validation
'''

# =======================================
# 📁 INPUT/OUTPUT PATHS CONFIGURATION
# =======================================

# Input from scVI preprocessing (01b script)
BASE_OUTPUT_DIR = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study"
PROCESSED_DATA_DIR_SCVI = f"{BASE_OUTPUT_DIR}/data/processed/snrna_scvi"
RESULTS_DIR = f"{BASE_OUTPUT_DIR}/results/snrna_scvi"
PLOTS_DIR = f"{RESULTS_DIR}/annotation_plots"

# Input files
INPUT_H5AD = f"{PROCESSED_DATA_DIR_SCVI}/GSE225158_reprocessed_scvi.h5ad"
MARKER_GENES_CSV = f"{PROCESSED_DATA_DIR_SCVI}/scvi_marker_genes.csv"

# Output files
OUTPUT_H5AD = f"{PROCESSED_DATA_DIR_SCVI}/GSE225158_annotated_scvi.h5ad"
ANNOTATION_SUMMARY = f"{PROCESSED_DATA_DIR_SCVI}/cell_type_annotation_summary_scvi.txt"

# =======================================
# 📚 IMPORT LIBRARIES
# =======================================

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend to prevent plot display
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import os

# scVI dependencies
try:
    import torch
    import scvi
    SCVI_AVAILABLE = True
    print("✅ scVI dependencies available")
except ImportError:
    SCVI_AVAILABLE = False
    print("⚠️  scVI dependencies not available")

warnings.filterwarnings('ignore')

# Set scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=300, facecolor='white')
sc.settings.autoshow = False

def run_reference_based_annotation(adata):
    """Run reference-based annotation using existing annotations if available"""
    print("🔍 Checking for reference-based annotation...")
    
    # Check if original annotations exist
    original_annotations = ['celltype1', 'celltype2', 'celltype3']
    available_refs = [col for col in original_annotations if col in adata.obs.columns]
    
    if len(available_refs) > 0:
        print(f"   ✅ Found reference annotations: {available_refs}")
        
        # Use the most detailed annotation available
        ref_col = available_refs[-1]  # Usually celltype3 is most detailed
        adata.obs['reference_annotation'] = adata.obs[ref_col]
        
        # Show reference distribution
        ref_counts = adata.obs['reference_annotation'].value_counts()
        print(f"   📊 Reference annotation distribution ({ref_col}):")
        for cell_type, count in ref_counts.head(10).items():
            pct = count / len(adata.obs) * 100
            print(f"     • {cell_type}: {count:,} cells ({pct:.1f}%)")
        
        return True
    else:
        print("   ⚠️  No reference annotations found")
        return False

def run_scanpy_leiden_annotation(adata, brain_markers):
    """Use scanpy's standard workflow for annotation validation"""
    print("🎯 Running scanpy-based annotation for validation...")
    
    # Calculate marker scores using standard approach
    adata = calculate_scvi_marker_scores(adata, brain_markers)
    
    # Find top expressed genes per cluster for validation
    if 'leiden' in adata.obs.columns:
        try:
            # Use scVI normalized expression for ranking
            if 'scvi_normalized' in adata.layers:
                adata.X = adata.layers['scvi_normalized'].copy()
            
            # Rank genes for each cluster
            sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', n_genes=25)
            
            # Extract top markers per cluster
            cluster_markers = {}
            for cluster in adata.obs['leiden'].unique():
                markers = sc.get.rank_genes_groups_df(adata, group=cluster)['names'].head(10).tolist()
                cluster_markers[f'Cluster_{cluster}'] = markers
                print(f"   • Cluster {cluster} top markers: {', '.join(markers[:5])}")
            
            # Store cluster-based predictions
            adata.obs['scanpy_prediction'] = [f'Cluster_{c}' for c in adata.obs['leiden']]
            
            return True, cluster_markers
            
        except Exception as e:
            print(f"   ⚠️  Scanpy annotation failed: {e}")
            return False, {}
    else:
        print("   ⚠️  No leiden clustering found")
        return False, {}

def create_multi_method_consensus(adata, methods_used):
    """Create consensus from multiple annotation methods"""
    print("🔀 Creating multi-method consensus...")
    
    consensus_annotations = []
    consensus_confidence = []
    method_agreement = []
    
    for cluster in adata.obs['leiden'].unique():
        cluster_mask = adata.obs['leiden'] == cluster
        cluster_cells = adata.obs.loc[cluster_mask]
        
        # Collect predictions from all methods
        cluster_predictions = {}
        
        # scVI method
        if 'scvi_prediction' in adata.obs.columns:
            scvi_pred = cluster_cells['scvi_prediction'].value_counts()
            if len(scvi_pred) > 0:
                cluster_predictions['scvi'] = scvi_pred.index[0]
        
        # Reference method
        if 'reference_annotation' in adata.obs.columns:
            ref_pred = cluster_cells['reference_annotation'].value_counts()
            if len(ref_pred) > 0:
                cluster_predictions['reference'] = ref_pred.index[0]
        
        # Scanpy method
        if 'scanpy_prediction' in adata.obs.columns:
            scanpy_pred = cluster_cells['scanpy_prediction'].value_counts()
            if len(scanpy_pred) > 0:
                cluster_predictions['scanpy'] = scanpy_pred.index[0]
        
        # Determine consensus
        if len(cluster_predictions) == 0:
            final_type = f"Unknown_Cluster_{cluster}"
            confidence = 'Low'
            agreement = 0
        elif len(cluster_predictions) == 1:
            final_type = list(cluster_predictions.values())[0]
            confidence = 'Medium'
            agreement = 1
        else:
            # Check agreement between methods
            unique_predictions = list(set(cluster_predictions.values()))
            agreement_score = len([p for p in cluster_predictions.values() if p == unique_predictions[0]]) / len(cluster_predictions)
            
            if agreement_score >= 0.67:  # 2/3 or more agree
                final_type = unique_predictions[0]
                confidence = 'High' if agreement_score == 1.0 else 'Medium'
            else:
                # Prioritize scVI if available, then reference
                if 'scvi' in cluster_predictions:
                    final_type = cluster_predictions['scvi']
                elif 'reference' in cluster_predictions:
                    final_type = cluster_predictions['reference']
                else:
                    final_type = list(cluster_predictions.values())[0]
                
                final_type = f"Mixed_{final_type}"
                confidence = 'Low'
            
            agreement = agreement_score
        
        print(f"   • Cluster {cluster}: {final_type} (confidence: {confidence}, agreement: {agreement:.2f})")
        print(f"     Methods: {', '.join(cluster_predictions.keys())}")
        
        # Apply to all cells in cluster
        for idx in cluster_cells.index:
            consensus_annotations.append(final_type)
            consensus_confidence.append(confidence)
            method_agreement.append(agreement)
    
    adata.obs['final_cell_type'] = consensus_annotations
    adata.obs['consensus_confidence'] = consensus_confidence
    adata.obs['method_agreement'] = method_agreement
    
    # Show final results
    print("\n📊 Final multi-method consensus:")
    final_counts = pd.Series(consensus_annotations).value_counts()
    for cell_type, count in final_counts.items():
        pct = count / len(consensus_annotations) * 100
        print(f"   • {cell_type}: {count:,} cells ({pct:.1f}%)")
    
    return adata

def calculate_scvi_marker_scores(adata, brain_markers):
    """Calculate marker scores using scVI's denoised expression"""
    print("🧬 Calculating scVI-optimized marker scores...")
    
    # Use scVI normalized expression (denoised and batch-corrected)
    if 'scvi_normalized' not in adata.layers:
        print("   ⚠️  scVI normalized layer not found, using X")
        expression_data = adata.X
    else:
        expression_data = adata.layers['scvi_normalized']
        print("   ✅ Using scVI denoised expression")
    
    # Calculate scores for each cell type
    for cell_type, markers in brain_markers.items():
        available_markers = [m for m in markers if m in adata.var_names]
        
        if len(available_markers) > 0:
            # Get marker expression
            marker_indices = [adata.var_names.get_loc(m) for m in available_markers]
            
            if hasattr(expression_data, 'toarray'):
                marker_expr = expression_data[:, marker_indices].toarray()
            else:
                marker_expr = expression_data[:, marker_indices]
            
            # Calculate robust score (mean + consideration of uncertainty if available)
            marker_score = np.mean(marker_expr, axis=1)
            
            # If uncertainty is available, weight by confidence
            if 'scvi_uncertainty' in adata.layers:
                uncertainty_data = adata.layers['scvi_uncertainty']
                if hasattr(uncertainty_data, 'toarray'):
                    marker_uncertainty = uncertainty_data[:, marker_indices].toarray()
                else:
                    marker_uncertainty = uncertainty_data[:, marker_indices]
                
                # Lower uncertainty = higher confidence weight
                confidence_weight = 1 / (1 + np.mean(marker_uncertainty, axis=1))
                marker_score = marker_score * confidence_weight
                print(f"   • {cell_type}: {len(available_markers)} markers (uncertainty-weighted)")
            else:
                print(f"   • {cell_type}: {len(available_markers)} markers")
            
            adata.obs[f'{cell_type}_score'] = marker_score
        else:
            print(f"   ⚠️  {cell_type}: No markers found")
            adata.obs[f'{cell_type}_score'] = 0
    
    return adata

def run_scvi_latent_annotation(adata, brain_markers):
    """Use scVI latent space for annotation with marker guidance"""
    print("🎯 Running scVI latent space annotation...")
    
    # First calculate marker scores
    adata = calculate_scvi_marker_scores(adata, brain_markers)
    
    # Get all score columns
    score_cols = [f'{ct}_score' for ct in brain_markers.keys()]
    score_cols = [col for col in score_cols if col in adata.obs.columns]
    
    if len(score_cols) == 0:
        print("   ⚠️  No marker scores available")
        return False
    
    # Use both marker scores and latent space clustering for annotation
    score_matrix = adata.obs[score_cols].values
    
    # Normalize scores within each cell (softmax-like)
    score_matrix_norm = score_matrix / (score_matrix.sum(axis=1, keepdims=True) + 1e-8)
    
    # Get top predictions
    max_indices = np.argmax(score_matrix_norm, axis=1)
    max_scores = np.max(score_matrix_norm, axis=1)
    second_max = np.partition(score_matrix_norm, -2, axis=1)[:, -2]
    
    # Calculate confidence as difference between top two scores
    confidence_scores = max_scores - second_max
    
    cell_types = [col.replace('_score', '') for col in score_cols]
    predicted_types = [cell_types[i] for i in max_indices]
    
    # Apply adaptive thresholds based on data distribution
    high_conf_threshold = np.percentile(confidence_scores, 75)
    med_conf_threshold = np.percentile(confidence_scores, 50)
    
    # Assign confidence levels
    confidence_levels = []
    final_predictions = []
    
    for i, (pred, conf, max_score) in enumerate(zip(predicted_types, confidence_scores, max_scores)):
        if conf > high_conf_threshold and max_score > 0.4:
            confidence_levels.append('High')
            final_predictions.append(pred)
        elif conf > med_conf_threshold and max_score > 0.25:
            confidence_levels.append('Medium')
            final_predictions.append(pred)
        else:
            confidence_levels.append('Low')
            final_predictions.append(f"Uncertain_{pred}")
    
    adata.obs['scvi_prediction'] = final_predictions
    adata.obs['scvi_confidence'] = confidence_levels
    adata.obs['scvi_score'] = max_scores
    adata.obs['scvi_conf_score'] = confidence_scores
    
    print("   ✅ scVI latent annotation completed")
    
    # Show results
    pred_counts = pd.Series(final_predictions).value_counts()
    for cell_type, count in pred_counts.items():
        pct = count / len(final_predictions) * 100
        print(f"   • {cell_type}: {count:,} cells ({pct:.1f}%)")
    
    return True

def main_annotation_scvi():
    """Multi-method cell type annotation pipeline for scVI data"""
    
    print("=" * 70)
    print("🧠 MULTI-METHOD CELL TYPE ANNOTATION PIPELINE")
    print("=" * 70)
    
    # =======================================
    # 📂 STEP 1: LOAD PREPROCESSED DATA
    # =======================================
    print("📂 STEP 1: LOADING PREPROCESSED DATA")
    print("=" * 70)
    
    # Load scVI data
    if os.path.exists(INPUT_H5AD):
        adata = sc.read_h5ad(INPUT_H5AD)
        print(f"📊 Loaded scVI-processed data: {adata.shape}")
    else:
        raise FileNotFoundError(f"scVI data not found! Run 01b_preprocessing_scvi.py first.\nExpected file: {INPUT_H5AD}")
    
    # Load marker genes if available
    if os.path.exists(MARKER_GENES_CSV):
        marker_df = pd.read_csv(MARKER_GENES_CSV)
        print(f"📋 Loaded {len(marker_df)} marker genes from scVI")
    else:
        print("⚠️  scVI marker genes not found, will compute fresh")
        marker_df = None

    # =======================================
    # 🧬 STEP 2: ENHANCED MARKER DEFINITION
    # =======================================
    print("\n🧬 STEP 2: ENHANCED MARKER DEFINITION")
    print("=" * 70)
    
    # Enhanced striatum-specific markers for CellAssign
    brain_markers = {
        'D1_MSN': ['DRD1', 'FOXP1', 'PPP1R1B', 'PDYN', 'BCAM', 'CRYM', 'EBF1', 'RXRG'],
        'D2_MSN': ['DRD2', 'FOXP2', 'PPP1R1B', 'PENK', 'GPR6', 'SP9', 'CRYM', 'ANKFN1'],
        'Cholinergic_Interneurons': ['CHAT', 'SLC5A7', 'SLC18A3', 'ISL1', 'ERBB4'],
        'PV_Interneurons': ['PVALB', 'GPR149', 'ERBB4', 'SULF1'],
        'SST_Interneurons': ['SST', 'CHODL', 'MEIS2', 'NPY', 'GRIK1'],
        'Astrocytes': ['GFAP', 'AQP4', 'S100B', 'ALDH1L1', 'SLC1A2', 'GJA1'],
        'Oligodendrocytes': ['MBP', 'MOG', 'PLP1', 'CNP', 'MOBP', 'MAG'],
        'OPCs': ['PDGFRA', 'CSPG4', 'SOX10', 'OLIG2', 'GPR17'],
        'Microglia': ['AIF1', 'CX3CR1', 'P2RY12', 'TMEM119', 'CSF1R'],
        'Endothelial': ['PECAM1', 'VWF', 'FLT1', 'CDH5', 'CLDN5'],
        'Pericytes': ['PDGFRB', 'RGS5', 'ACTA2', 'NOTCH3']
    }
    
    print(f"🧬 Using {len(brain_markers)} striatum cell type definitions")
    
    # =======================================
    # 🎯 STEP 3: MULTI-METHOD ANNOTATION
    # =======================================
    print("\n🎯 STEP 3: MULTI-METHOD ANNOTATION")
    print("=" * 70)
    
    methods_used = []
    
    # Method 1: scVI-optimized annotation
    scvi_success = run_scvi_latent_annotation(adata, brain_markers)
    if scvi_success:
        methods_used.append('scVI-optimized')
    
    # Method 2: Reference-based annotation
    ref_success = run_reference_based_annotation(adata)
    if ref_success:
        methods_used.append('Reference-based')
    
    # Method 3: Scanpy standard workflow
    scanpy_success, cluster_markers = run_scanpy_leiden_annotation(adata, brain_markers)
    if scanpy_success:
        methods_used.append('Scanpy-standard')
    
    print(f"\n✅ Active annotation methods: {', '.join(methods_used)}")
    
    # =======================================
    # 🔀 STEP 4: MULTI-METHOD CONSENSUS
    # =======================================
    print("\n🔀 STEP 4: MULTI-METHOD CONSENSUS")
    print("=" * 70)
    
    adata = create_multi_method_consensus(adata, methods_used)
    
    # =======================================
    # 📊 STEP 5: ENHANCED VALIDATION PLOTS
    # =======================================
    print("\n📊 STEP 5: ENHANCED VALIDATION PLOTS")
    print("=" * 70)
    
    os.makedirs(PLOTS_DIR, exist_ok=True)
    
    # Create validation plots showing all methods
    n_methods = len(methods_used)
    fig, axes = plt.subplots(2, 3, figsize=(24, 16))
    
    # Final consensus
    sc.pl.umap(adata, color='final_cell_type', ax=axes[0,0], show=False, 
               legend_loc='right margin', legend_fontsize=8, size=3)
    axes[0,0].set_title('Final Consensus Annotation', fontsize=16)
    
    # Method agreement
    sc.pl.umap(adata, color='method_agreement', ax=axes[0,1], show=False, size=3)
    axes[0,1].set_title('Method Agreement Score', fontsize=16)
    
    # Consensus confidence
    sc.pl.umap(adata, color='consensus_confidence', ax=axes[0,2], show=False, size=3)
    axes[0,2].set_title('Consensus Confidence', fontsize=16)
    
    # Individual method results
    if 'scvi_prediction' in adata.obs.columns:
        sc.pl.umap(adata, color='scvi_prediction', ax=axes[1,0], show=False, 
                   legend_loc='right margin', legend_fontsize=6, size=3)
        axes[1,0].set_title('scVI Method', fontsize=16)
    
    if 'reference_annotation' in adata.obs.columns:
        sc.pl.umap(adata, color='reference_annotation', ax=axes[1,1], show=False, 
                   legend_loc='right margin', legend_fontsize=6, size=3)
        axes[1,1].set_title('Reference Annotation', fontsize=16)
    
    # OUD diagnosis for clinical context
    if 'Dx_OUD' in adata.obs.columns:
        sc.pl.umap(adata, color='Dx_OUD', ax=axes[1,2], show=False, size=3)
        axes[1,2].set_title('OUD Diagnosis', fontsize=16)
    
    plt.suptitle('Multi-Method Cell Type Annotation Validation', fontsize=20, y=0.98)
    plt.tight_layout()
    plt.savefig(f"{PLOTS_DIR}/multi_method_annotation_validation.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✅ Enhanced plots saved to: {PLOTS_DIR}/")

    # =======================================
    # 💾 STEP 6: SAVE RESULTS
    # =======================================
    print("\n💾 STEP 6: SAVING RESULTS")
    print("=" * 70)
    
    # Set final annotation
    adata.obs['cell_type'] = adata.obs['final_cell_type']
    
    # Save annotated data
    adata.write(OUTPUT_H5AD)
    print(f"💾 scVI annotated data saved to: {OUTPUT_H5AD}")
    
    # Save comprehensive summary
    final_counts = adata.obs['final_cell_type'].value_counts()
    
    with open(ANNOTATION_SUMMARY, 'w') as f:
        f.write("GSE225158 Cell Type Annotation Summary (Multi-Method Validation)\n")
        f.write("=" * 70 + "\n")
        f.write(f"Date: {pd.Timestamp.now()}\n")
        f.write(f"Data source: scVI processed data\n\n")
        
        f.write("ANNOTATION METHODS USED:\n")
        for method in methods_used:
            f.write(f"• {method}\n")
        f.write("• Multi-method consensus validation\n\n")
        
        f.write("FINAL CELL TYPE COUNTS:\n")
        for cell_type, count in final_counts.items():
            pct = count / len(adata.obs) * 100
            f.write(f"• {cell_type}: {count:,} cells ({pct:.1f}%)\n")
        
        f.write(f"\nTOTAL CELLS: {len(adata.obs):,}\n")
        f.write(f"CELL TYPES: {len(final_counts)}\n")

    print(f"💾 Summary saved to: {ANNOTATION_SUMMARY}")
    
    # Brain region analysis
    if 'Region' in adata.obs.columns:
        print("\n🧠 Brain region-specific analysis:")
        region_analysis = pd.crosstab(adata.obs['Region'], adata.obs['final_cell_type'])
        print(region_analysis)
        
        region_path = f"{PROCESSED_DATA_DIR_SCVI}/region_celltype_analysis_scvi.csv"
        region_analysis.to_csv(region_path)
        print(f"💾 Region analysis saved to: {region_path}")
    
    print("\n" + "=" * 70)
    print("✅ MULTI-METHOD ANNOTATION COMPLETE!")
    print("=" * 70)
    print(f"📊 Dataset: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    print(f"🧠 Cell types: {len(final_counts)}")
    print(f"🔬 Methods used: {', '.join(methods_used)}")
    print("=" * 70)
    
    return adata

if __name__ == "__main__":
    annotated_adata = main_annotation_scvi()
