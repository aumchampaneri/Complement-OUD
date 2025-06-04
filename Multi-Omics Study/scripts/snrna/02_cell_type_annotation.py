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
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend to prevent plot display
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import os

# Optional CellTypist for reference mapping
try:
    import celltypist
    CELLTYPIST_AVAILABLE = True
    print("✅ CellTypist available for reference mapping")
except ImportError:
    CELLTYPIST_AVAILABLE = False
    print("⚠️  CellTypist not available - using marker-based only")

warnings.filterwarnings('ignore')

# Set scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=300, facecolor='white')
# Ensure scanpy doesn't show plots
sc.settings.autoshow = False

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
    
    # STATISTICAL VALIDATION: Find cluster marker genes
    print("🔬 Computing cluster-specific marker genes for validation...")
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', use_raw=False)
    
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
            
            # MEMORY-OPTIMIZED: Check if cluster expresses the markers
            # Use sparse-safe operations instead of .toarray()
            from scipy import sparse
            if sparse.issparse(cluster_cells[:, marker_genes].X):
                # For sparse matrices, use .A (equivalent to .toarray() but more memory efficient)
                marker_expression = cluster_cells[:, marker_genes].X.mean(axis=0).A1
                mean_marker_expr = np.mean(marker_expression)
                
                # Calculate background expression more efficiently
                if sparse.issparse(adata.X):
                    background_threshold = np.percentile(adata.X.data, 25)  # Only use non-zero values
                else:
                    background_threshold = np.percentile(adata.X, 25)
            else:
                # For dense matrices
                marker_expression = cluster_cells[:, marker_genes].X.mean(axis=0)
                mean_marker_expr = np.mean(marker_expression)
                background_threshold = np.percentile(adata.X, 25)
            
            # STATISTICAL VALIDATION: Check if cluster markers support prediction
            cluster_markers = adata.uns['rank_genes_groups']['names'][str(cluster)][:10]  # Top 10 markers
            marker_overlap = len(set(marker_genes) & set(cluster_markers))
            
            print(f"     Cluster {cluster} -> {most_common}: marker overlap {marker_overlap}/{len(marker_genes)}")
            
            # If low marker expression OR low statistical support, mark as uncertain
            if mean_marker_expr < background_threshold or marker_overlap < max(1, len(marker_genes) * 0.3):
                most_common = f"Uncertain_{most_common}"
                print(f"     ⚠️  Low confidence: expr={mean_marker_expr:.3f}, overlap={marker_overlap}")
        
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
    # 🎯 STEP 2.5: OPTIONAL CELLTYPIST REFERENCE MAPPING
    # =======================================
    if CELLTYPIST_AVAILABLE:
        print("\n🎯 STEP 2.5: CELLTYPIST REFERENCE MAPPING")
        print("=" * 70)
        
        celltypist_success = False  # Track if CellTypist works
        
        try:
            print("🔍 Loading CellTypist models...")
            models = celltypist.models.get_all_models()
            print(f"   Available models: {len(models)}")
            
            # Prioritized models for human striatum (adult brain tissue)
            preferred_models = [
                'Adult_Human_PrefrontalCortex.pkl',    # Best - adult human brain cortex (113 types)
                'Adult_Human_MTG.pkl',                 # Excellent - adult human middle temporal gyrus (127 types)
                'Human_AdultAged_Hippocampus.pkl',     # Good - adult human hippocampus (15 types)
                'Human_Longitudinal_Hippocampus.pkl',  # Good - adult human hippocampus (24 types)
                'Developing_Human_Brain.pkl',          # Developmental but human brain (129 types)
                'Adult_CynomolgusMacaque_Hippocampus.pkl', # Non-human primate brain (34 types)
                'Mouse_Isocortex_Hippocampus.pkl',     # Mouse brain but relevant (42 types)
                'Mouse_Whole_Brain.pkl',               # Comprehensive mouse brain (334 types)
                'Pan_Fetal_Human.pkl',                 # General human fetal (138 types)
                'Immune_All_Low.pkl'                   # Fallback general immune (98 types)
            ]
            
            selected_model = None
            for model in preferred_models:
                if model in models:
                    selected_model = model
                    print(f"   🧠 Selected optimal brain model: {selected_model}")
                    
                    # Provide context about the selected model
                    model_info = {
                        'Adult_Human_PrefrontalCortex.pkl': "Adult human prefrontal cortex - ideal for cortical inputs to striatum",
                        'Adult_Human_MTG.pkl': "Adult human middle temporal gyrus - excellent cortical reference",
                        'Human_AdultAged_Hippocampus.pkl': "Adult human hippocampus - good for limbic connections",
                        'Human_Longitudinal_Hippocampus.pkl': "Adult human hippocampus - longitudinal study",
                        'Developing_Human_Brain.pkl': "Developing human brain - comprehensive but developmental",
                        'Adult_CynomolgusMacaque_Hippocampus.pkl': "Non-human primate brain - close evolutionary match",
                        'Mouse_Isocortex_Hippocampus.pkl': "Mouse brain - well-characterized reference",
                        'Mouse_Whole_Brain.pkl': "Comprehensive mouse brain atlas",
                        'Pan_Fetal_Human.pkl': "General human fetal tissue reference",
                        'Immune_All_Low.pkl': "General immune cell reference"
                    }
                    
                    print(f"      📋 {model_info.get(selected_model, 'Brain-relevant model')}")
                    break
            
            if not selected_model:
                print("   ⚠️  No preferred brain models found. Available models:")
                for i, model in enumerate(models[:15]):
                    print(f"      {i}: {model}")
                selected_model = models[0] if models else None
                print(f"   Using first available model: {selected_model}")
            
            if selected_model:
                # PREPARE DATA FOR CELLTYPIST
                print("🔧 Preparing data for CellTypist...")
                
                # CellTypist expects log1p normalized data to 10,000 counts per cell
                # Our Pearson residuals data needs to be converted
                
                # Create a copy for CellTypist with proper normalization
                adata_ct = adata.copy()
                
                # Check if we have raw counts available
                if adata_ct.raw is not None:
                    print("   Using raw counts for CellTypist normalization")
                    # Use raw counts and normalize
                    adata_ct.X = adata_ct.raw.X.copy()
                    sc.pp.normalize_total(adata_ct, target_sum=1e4)
                    sc.pp.log1p(adata_ct)
                else:
                    print("   ⚠️  No raw counts available - converting Pearson residuals")
                    # Convert Pearson residuals to pseudo-normalized data
                    # This is a workaround since we don't have raw counts
                    
                    # Shift and scale Pearson residuals to positive range
                    from scipy import sparse
                    if sparse.issparse(adata_ct.X):
                        X_min = adata_ct.X.data.min() if adata_ct.X.data.size > 0 else 0
                    else:
                        X_min = adata_ct.X.min()
                    
                    # Shift to make all values positive
                    if X_min < 0:
                        adata_ct.X = adata_ct.X - X_min
                    
                    # Scale to reasonable range (simulate normalized counts)
                    if sparse.issparse(adata_ct.X):
                        scale_factor = 10.0 / (adata_ct.X.data.max() if adata_ct.X.data.size > 0 else 1.0)
                    else:
                        scale_factor = 10.0 / adata_ct.X.max()
                    
                    adata_ct.X = adata_ct.X * scale_factor
                    
                    # Apply log1p
                    sc.pp.log1p(adata_ct)
                
                print("   ✅ Data prepared for CellTypist")
                
                # Run CellTypist prediction
                print("🔬 Running CellTypist annotation...")
                predictions = celltypist.annotate(
                    adata_ct, 
                    model=selected_model, 
                    majority_voting=True,  # Use majority voting for robustness
                    over_clustering='leiden'  # Use existing leiden clusters
                )
                
                # Validate CellTypist predictions before using them
                print("🔍 Validating CellTypist predictions...")
                predicted_labels = predictions.predicted_labels['predicted_labels']
                
                # Clean CellTypist predictions - remove appended marker genes
                print("🧹 Cleaning CellTypist predictions (removing appended markers)...")
                cleaned_labels = []
                for label in predicted_labels:
                    # Common patterns to clean
                    clean_label = str(label)
                    
                    # Remove common marker genes appended to cell type names
                    marker_patterns = [
                        r'\s+MOG\s*.*$',           # Remove " MOG ..." 
                        r'\s+OPALIN\s*.*$',        # Remove " OPALIN ..."
                        r'\s+PDGFRA\s*.*$',        # Remove " PDGFRA ..."
                        r'\s+P2RY12\s*.*$',        # Remove " P2RY12 ..."
                        r'\s+GFAP\s*.*$',          # Remove " GFAP ..."
                        r'\s+AQP4\s*.*$',          # Remove " AQP4 ..."
                        r'\s+[A-Z0-9]+\s+[A-Z0-9]+.*$',  # Remove " GENE1 GENE2 ..."
                    ]
                    
                    import re
                    for pattern in marker_patterns:
                        clean_label = re.sub(pattern, '', clean_label)
                    
                    # Clean up cell type names to standard format
                    clean_label = clean_label.strip()
                    
                    # Standardize common brain cell type names
                    cell_type_mapping = {
                        'Oligo': 'Oligodendrocyte',
                        'OPC': 'OPC',
                        'Micro': 'Microglia',
                        'Astro': 'Astrocyte',
                        'Neuron': 'Neuron',
                        'Endo': 'Endothelial',
                        'Pericyte': 'Pericyte'
                    }
                    
                    # Apply standardization
                    for short_name, full_name in cell_type_mapping.items():
                        if clean_label.startswith(short_name):
                            clean_label = full_name
                            break
                    
                    cleaned_labels.append(clean_label)
                
                # Replace original predictions with cleaned ones
                predicted_labels = pd.Series(cleaned_labels)
                
                # Show before/after cleaning example
                original_unique = predictions.predicted_labels['predicted_labels'].unique()
                cleaned_unique = predicted_labels.unique()
                print(f"   📋 Cleaned predictions:")
                print(f"   Before: {original_unique[:3]}...")
                print(f"   After:  {cleaned_unique[:3]}...")
                
                # Check for corrupted predictions (repeated strings, invalid characters)
                unique_predictions = predicted_labels.unique()
                corrupted_predictions = []
                for pred in unique_predictions:
                    if (len(str(pred)) > 50 or  # Abnormally long cell type names
                        str(pred).count('OG') > 3 or  # Repeated patterns
                        str(pred).count('OPALIN') > 1 or
                        any(char * 3 in str(pred) for char in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ')):  # Triple repeated characters
                        corrupted_predictions.append(pred)
                
                if corrupted_predictions:
                    print(f"⚠️  Detected corrupted CellTypist predictions: {len(corrupted_predictions)} types")
                    print(f"   Examples: {corrupted_predictions[:3]}")
                    print("   This may be due to memory issues or model corruption")
                    
                    # Try to clean predictions or skip CellTypist
                    valid_predictions = [pred for pred in unique_predictions if pred not in corrupted_predictions]
                    
                    if len(valid_predictions) < 3:  # Too few valid predictions
                        print("   ❌ Too few valid predictions - skipping CellTypist")
                        raise Exception("Corrupted CellTypist predictions detected")
                    else:
                        print(f"   🔧 Using only valid predictions: {len(valid_predictions)} types")
                        # Replace corrupted predictions with 'Unknown'
                        cleaned_predictions = predicted_labels.copy()
                        for corrupt_pred in corrupted_predictions:
                            cleaned_predictions = cleaned_predictions.replace(corrupt_pred, 'Unknown_CellType')
                        predicted_labels = cleaned_predictions
                
                # Extract predictions back to original adata with error handling
                adata.obs['celltypist_prediction'] = predicted_labels
                
                # Handle confidence scores - different CellTypist versions may use different column names
                if 'conf_score' in predictions.predicted_labels.columns:
                    confidence_scores = predictions.predicted_labels['conf_score']
                elif 'confidence' in predictions.predicted_labels.columns:
                    confidence_scores = predictions.predicted_labels['confidence']
                elif 'max_prob' in predictions.predicted_labels.columns:
                    confidence_scores = predictions.predicted_labels['max_prob']
                else:
                    # Calculate confidence from probability columns if available
                    prob_cols = [col for col in predictions.predicted_labels.columns if col not in ['predicted_labels']]
                    if prob_cols:
                        prob_matrix = predictions.predicted_labels[prob_cols].values
                        confidence_scores = pd.Series(np.max(prob_matrix, axis=1))
                        print(f"   📊 Calculated confidence from {len(prob_cols)} probability columns")
                    else:
                        # Fallback: assign medium confidence
                        confidence_scores = pd.Series([0.5] * len(predicted_labels))
                        print("   ⚠️  No confidence scores found - using default 0.5")
                
                # Robust confidence score validation and cleaning
                print("   🔧 Cleaning and validating confidence scores...")
                try:
                    # Convert to pandas Series if not already
                    if not isinstance(confidence_scores, pd.Series):
                        confidence_scores = pd.Series(confidence_scores)
                    
                    # Convert to numeric, forcing errors to NaN
                    confidence_scores = pd.to_numeric(confidence_scores, errors='coerce')
                    
                    # Fill NaN values with 0.5 (medium confidence)
                    confidence_scores = confidence_scores.fillna(0.5)
                    
                    # Clip values to valid range [0, 1]
                    confidence_scores = confidence_scores.clip(0, 1)
                    
                    # Final validation - ensure all values are numeric
                    if not confidence_scores.dtype.kind in 'biufc':  # numeric types
                        print("   ⚠️  Non-numeric confidence scores detected - using default")
                        confidence_scores = pd.Series([0.5] * len(predicted_labels))
                    
                    print(f"   ✅ Confidence scores cleaned: range [{confidence_scores.min():.3f}, {confidence_scores.max():.3f}]")
                    
                except Exception as conf_error:
                    print(f"   ⚠️  Error processing confidence scores: {conf_error}")
                    confidence_scores = pd.Series([0.5] * len(predicted_labels))
                    print("   Using default confidence scores of 0.5")
                
                adata.obs['celltypist_confidence'] = confidence_scores
                adata.obs['celltypist_model_used'] = selected_model  # Track which model was used
                
                print("✅ CellTypist predictions completed and validated")
                celltypist_success = True
                
                # Summary of CellTypist results with validation
                ct_counts = adata.obs['celltypist_prediction'].value_counts()
                print("📊 CellTypist predictions (validated):")
                
                # Check if we have reasonable predictions
                if len(ct_counts) < 2:
                    print("   ⚠️  Too few cell types predicted - may indicate issues")
                    celltypist_success = False
                elif 'Unknown_CellType' in ct_counts.index and ct_counts['Unknown_CellType'] > len(adata.obs) * 0.5:
                    print("   ⚠️  Too many unknown predictions - may indicate issues")
                    celltypist_success = False
                else:
                    for cell_type, count in ct_counts.head(10).items():
                        pct = count / len(adata.obs) * 100
                        try:
                            # Safe confidence calculation
                            mask = adata.obs['celltypist_prediction'] == cell_type
                            if mask.sum() > 0:
                                avg_conf = adata.obs.loc[mask, 'celltypist_confidence'].mean()
                                print(f"   • {cell_type}: {count:,} cells ({pct:.1f}%) - avg conf: {avg_conf:.3f}")
                            else:
                                print(f"   • {cell_type}: {count:,} cells ({pct:.1f}%) - conf: N/A")
                        except Exception as summary_error:
                            print(f"   • {cell_type}: {count:,} cells ({pct:.1f}%) - conf: Error")
                
                # Validate brain relevance
                brain_types = ct_counts.index.tolist()
                brain_keywords = ['neuron', 'astro', 'oligo', 'microglia', 'brain', 'neural', 'glia', 'cortical', 'pyramidal', 'interneuron']
                brain_relevant = any(keyword in ' '.join(brain_types).lower() for keyword in brain_keywords)
                
                if brain_relevant:
                    print("   ✅ Model appears highly brain-relevant based on predicted cell types")
                else:
                    print("   ⚠️  Model may not be optimal for brain tissue")
                    print(f"   📋 Top predictions: {', '.join(ct_counts.head(5).index.tolist())}")
            else:
                raise Exception("No CellTypist models available")
                
        except Exception as e:
            print(f"⚠️  CellTypist failed: {e}")
            print("   Continuing with marker-based annotation only...")
            celltypist_success = False

    # =======================================
    # 🔀 STEP 3.5: CONSENSUS ANNOTATION (IF CELLTYPIST AVAILABLE)
    # =======================================
    if CELLTYPIST_AVAILABLE and celltypist_success and 'celltypist_prediction' in adata.obs.columns:
        print("\n🔀 STEP 3.5: CREATING CONSENSUS ANNOTATIONS")
        print("=" * 70)
        
        # Enhanced mapping for brain tissue
        celltypist_to_marker_mapping = {
            # Neurons - more specific for brain
            'Neuron': ['D1_MSN', 'D2_MSN', 'Cortical_Inputs'],
            'GABAergic neuron': ['D1_MSN', 'D2_MSN', 'PV_Interneurons', 'SST_Interneurons', 'GABA_Interneurons'],
            'Glutamatergic neuron': ['Cortical_Inputs'],
            'Inhibitory neuron': ['PV_Interneurons', 'SST_Interneurons', 'GABA_Interneurons'],
            'Excitatory neuron': ['Cortical_Inputs'],
            'Projection neuron': ['D1_MSN', 'D2_MSN'],
            'Interneuron': ['PV_Interneurons', 'SST_Interneurons', 'Cholinergic_Interneurons'],
            
            # Glia - brain specific
            'Astrocyte': ['Astrocytes'],
            'Oligodendrocyte': ['Oligodendrocytes'], 
            'OPC': ['OPCs'],
            'Oligodendrocyte precursor cell': ['OPCs'],
            'Microglia': ['Microglia'],
            'Microglial cell': ['Microglia'],
            
            # Vascular - brain specific
            'Endothelial cell': ['Endothelial'],
            'Brain endothelial cell': ['Endothelial'],
            'Pericyte': ['Pericytes'],
            'Vascular smooth muscle cell': ['Pericytes'],
            'VLMC': ['VLMC'],
            
            # Rare brain cells
            'Ependymal cell': ['Ependymal'],
            'Choroidal epithelial cell': ['Ependymal'],
            'Tanycyte': ['Tanycytes'],
            
            # General fallbacks
            'Neural cell': ['D1_MSN', 'D2_MSN', 'Cortical_Inputs'],
            'Brain cell': ['D1_MSN', 'D2_MSN', 'Astrocytes'],
        }
        
        print("🔀 Computing consensus between marker-based and CellTypist predictions...")
        consensus_predictions = []
        consensus_confidence = []
        
        for i in range(len(adata.obs)):
            marker_pred = adata.obs.iloc[i]['predicted_cell_type']
            marker_conf = adata.obs.iloc[i]['prediction_confidence']
            ct_pred = adata.obs.iloc[i]['celltypist_prediction']
            ct_conf = adata.obs.iloc[i]['celltypist_confidence']
            
            # Check if predictions are compatible
            compatible_markers = celltypist_to_marker_mapping.get(ct_pred, [])
            is_compatible = marker_pred in compatible_markers
            
            # Consensus logic
            if is_compatible and ct_conf > 0.7 and marker_conf > 0.3:
                # High confidence agreement
                consensus_predictions.append(marker_pred)  # Use more specific marker prediction
                consensus_confidence.append(max(marker_conf, ct_conf))
            elif ct_conf > 0.8 and marker_conf < 0.3:
                # High confidence CellTypist, low confidence marker
                # Map to best compatible marker type
                if compatible_markers:
                    consensus_predictions.append(compatible_markers[0])
                else:
                    consensus_predictions.append(f"CT_{ct_pred}")
                consensus_confidence.append(ct_conf * 0.9)
            elif marker_conf > 0.5:
                # Higher confidence marker prediction
                consensus_predictions.append(marker_pred)
                consensus_confidence.append(marker_conf * 0.8)
            else:
                # Low confidence both - mark as uncertain
                consensus_predictions.append(f"Uncertain_{marker_pred}")
                consensus_confidence.append(min(marker_conf, ct_conf))
        
        # Update predictions with consensus
        adata.obs['consensus_prediction'] = consensus_predictions
        adata.obs['consensus_confidence'] = consensus_confidence
        adata.obs['predicted_cell_type'] = consensus_predictions  # Use consensus as main prediction
        adata.obs['prediction_confidence'] = consensus_confidence
        
        print("✅ Consensus predictions created")

    # =======================================
    # 📊 STEP 5: ENHANCED VISUALIZATION
    # =======================================
    print("\n📊 STEP 5: GENERATING VISUALIZATIONS")
    print("=" * 70)
    
    os.makedirs(PLOTS_DIR, exist_ok=True)
    
    # Enhanced visualization comparing all annotation methods
    if CELLTYPIST_AVAILABLE and celltypist_success and 'celltypist_prediction' in adata.obs.columns:
        # Create comprehensive 3x2 layout for all methods
        fig, axes = plt.subplots(3, 2, figsize=(20, 24))
        
        # ROW 1: Original vs Final annotations
        if 'celltype1' in adata.obs.columns:
            sc.pl.umap(adata, color='celltype1', ax=axes[0,0], show=False, 
                       legend_loc='right margin', legend_fontsize=8, size=3)
            axes[0,0].set_title('Original Paper Annotations', fontsize=16, weight='bold')
        else:
            sc.pl.umap(adata, color='leiden', ax=axes[0,0], show=False, 
                       legend_loc='on data', legend_fontsize=8, size=3)
            axes[0,0].set_title('Leiden Clusters (No Original Annotations)', fontsize=16, weight='bold')
        
        sc.pl.umap(adata, color='refined_cell_type', ax=axes[0,1], show=False, 
                   legend_loc='right margin', legend_fontsize=8, size=3)
        axes[0,1].set_title('Final Refined Annotations', fontsize=16, weight='bold')
        
        # ROW 2: Marker-based vs CellTypist
        sc.pl.umap(adata, color='predicted_cell_type', ax=axes[1,0], show=False, 
                   legend_loc='right margin', legend_fontsize=8, size=3)
        axes[1,0].set_title('Marker-Based Predictions', fontsize=16, weight='bold')
        
        sc.pl.umap(adata, color='celltypist_prediction', ax=axes[1,1], show=False, 
                   legend_loc='right margin', legend_fontsize=8, size=3)
        axes[1,1].set_title('CellTypist Predictions\n(Adult Human Prefrontal Cortex)', fontsize=16, weight='bold')
        
        # ROW 3: Confidence scores and additional info
        sc.pl.umap(adata, color='final_confidence', ax=axes[2,0], show=False, size=3,
                   color_map='viridis')
        axes[2,0].set_title('Final Annotation Confidence', fontsize=16, weight='bold')
        
        sc.pl.umap(adata, color='celltypist_confidence', ax=axes[2,1], show=False, size=3,
                   color_map='plasma')
        axes[2,1].set_title('CellTypist Confidence Scores', fontsize=16, weight='bold')
        
        plt.suptitle('Comprehensive Cell Type Annotation Comparison\n(Marker-Based + CellTypist)', 
                     fontsize=20, weight='bold', y=0.98)
        plt.tight_layout()
        plt.savefig(f"{PLOTS_DIR}/comprehensive_annotation_comparison.png", dpi=300, bbox_inches='tight')
        plt.close()  # Explicitly close to prevent display
        
        # Create focused comparison plot (2x2)
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        # Original vs marker-based vs CellTypist vs final
        if 'celltype1' in adata.obs.columns:
            sc.pl.umap(adata, color='celltype1', ax=axes[0,0], show=False, 
                       legend_loc='right margin', legend_fontsize=6, size=3)
            axes[0,0].set_title('Original Annotations', fontsize=14)
        else:
            sc.pl.umap(adata, color='leiden', ax=axes[0,0], show=False, 
                       legend_loc='on data', legend_fontsize=8, size=3)
            axes[0,0].set_title('Leiden Clusters', fontsize=14)
        
        sc.pl.umap(adata, color='predicted_cell_type', ax=axes[0,1], show=False, 
                   legend_loc='right margin', legend_fontsize=6, size=3)
        axes[0,1].set_title('Marker-Based Predictions', fontsize=14)
        
        sc.pl.umap(adata, color='celltypist_prediction', ax=axes[1,0], show=False, 
                   legend_loc='right margin', legend_fontsize=6, size=3)
        axes[1,0].set_title('CellTypist Predictions', fontsize=14)
        
        sc.pl.umap(adata, color='refined_cell_type', ax=axes[1,1], show=False, 
                   legend_loc='right margin', legend_fontsize=6, size=3)
        axes[1,1].set_title('Final Refined Annotations', fontsize=14)
        
        plt.suptitle('Cell Type Annotation Methods Comparison', fontsize=16, y=0.98)
        
    else:
        # Fallback visualization without CellTypist (2x2 layout)
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        # Original annotations
        if 'celltype1' in adata.obs.columns:
            sc.pl.umap(adata, color='celltype1', ax=axes[0,0], show=False, 
                       legend_loc='right margin', legend_fontsize=8, size=3)
            axes[0,0].set_title('Original Paper Annotations', fontsize=16, weight='bold')
        else:
            sc.pl.umap(adata, color='leiden', ax=axes[0,0], show=False, 
                       legend_loc='on data', legend_fontsize=10, size=3)
            axes[0,0].set_title('Leiden Clusters', fontsize=16, weight='bold')
        
        # Marker-based predictions
        sc.pl.umap(adata, color='predicted_cell_type', ax=axes[0,1], show=False, 
                   legend_loc='right margin', legend_fontsize=8, size=3)
        axes[0,1].set_title('Marker-Based Predictions', fontsize=16, weight='bold')
        
        # Final refined annotations
        sc.pl.umap(adata, color='refined_cell_type', ax=axes[1,0], show=False, 
                   legend_loc='right margin', legend_fontsize=8, size=3)
        axes[1,0].set_title('Final Refined Annotations', fontsize=16, weight='bold')
        
        # Confidence scores
        sc.pl.umap(adata, color='final_confidence', ax=axes[1,1], show=False, size=3,
                   color_map='viridis')
        axes[1,1].set_title('Annotation Confidence', fontsize=16, weight='bold')
        
        plt.suptitle('Cell Type Annotation Results (Marker-Based Only)', 
                     fontsize=18, weight='bold', y=0.98)
    
    plt.tight_layout()
    plt.savefig(f"{PLOTS_DIR}/cell_type_annotation_pearson.png", dpi=300, bbox_inches='tight')
    plt.close()  # Explicitly close to prevent display

    # Create individual method plots for detailed inspection
    if CELLTYPIST_AVAILABLE and celltypist_success and 'celltypist_prediction' in adata.obs.columns:
        print("📊 Creating individual annotation method plots...")
        
        # Individual plot for marker-based
        plt.figure(figsize=(10, 8))
        sc.pl.umap(adata, color='predicted_cell_type', show=False, 
                   legend_loc='right margin', legend_fontsize=10, size=5)
        plt.title('Marker-Based Cell Type Predictions\n(Striatum-Specific Gene Signatures)', 
                  fontsize=16, weight='bold', pad=20)
        plt.savefig(f"{PLOTS_DIR}/marker_based_annotations.png", dpi=300, bbox_inches='tight')
        plt.close()  # Explicitly close to prevent display
        
        # Individual plot for CellTypist
        plt.figure(figsize=(10, 8))
        sc.pl.umap(adata, color='celltypist_prediction', show=False, 
                   legend_loc='right margin', legend_fontsize=10, size=5)
        plt.title('CellTypist Reference-Based Predictions\n(Adult Human Prefrontal Cortex Model)', 
                  fontsize=16, weight='bold', pad=20)
        plt.savefig(f"{PLOTS_DIR}/celltypist_annotations.png", dpi=300, bbox_inches='tight')
        plt.close()  # Explicitly close to prevent display
        
        # Individual plot for original (if available)
        if 'celltype1' in adata.obs.columns:
            plt.figure(figsize=(10, 8))
            sc.pl.umap(adata, color='celltype1', show=False, 
                       legend_loc='right margin', legend_fontsize=10, size=5)
            plt.title('Original Paper Annotations\n(GSE225158 Published Cell Types)', 
                      fontsize=16, weight='bold', pad=20)
            plt.savefig(f"{PLOTS_DIR}/original_annotations.png", dpi=300, bbox_inches='tight')
            plt.close()  # Explicitly close to prevent display

    # =======================================
    # 💾 STEP 6: SAVE ANNOTATED DATA
    # =======================================
    print("\n💾 STEP 6: SAVING ANNOTATED DATA")
    print("=" * 70)
    
    # Save the final annotated dataset
    adata.write(OUTPUT_H5AD)
    print(f"✅ Annotated data saved to: {OUTPUT_H5AD}")
    
    # Summary of saved annotations
    print(f"📊 Saved annotations include:")
    annotation_columns = [col for col in adata.obs.columns if any(term in col for term in 
                         ['predicted', 'refined', 'confidence', 'celltypist', 'consensus'])]
    for col in annotation_columns:
        print(f"   • {col}")

    # Enhanced summary with CellTypist info
    with open(ANNOTATION_SUMMARY, 'w') as f:
        f.write("GSE225158 Cell Type Annotation Summary (Pearson Residuals)\n")
        f.write("=" * 60 + "\n")
        f.write(f"Date: {pd.Timestamp.now()}\n")
        f.write(f"Data source: Pearson residuals normalization\n\n")
        
        f.write("ANNOTATION METHODS USED:\n")
        f.write("• Marker-based signature scoring\n")
        if CELLTYPIST_AVAILABLE and celltypist_success and 'celltypist_prediction' in adata.obs.columns:
            f.write("• CellTypist reference mapping\n")
            f.write("• Consensus prediction combining both methods\n")
        f.write("• Cluster-based refinement\n")
        f.write("• Statistical validation with cluster markers\n\n")
        
        f.write("FINAL CELL TYPE COUNTS:\n")
        for cell_type, count in final_counts.items():
            pct = count / len(adata.obs) * 100
            f.write(f"• {cell_type}: {count:,} cells ({pct:.1f}%)\n")
        
        f.write(f"\nTOTAL CELLS: {len(adata.obs):,}\n")
        f.write(f"CELL TYPES: {len(final_counts)}\n\n")
        
        f.write("MARKER GENES USED:\n")
        for cell_type, markers in signature_scores.items():
            f.write(f"• {cell_type}: {', '.join(markers)}\n")
        
        # Add cluster validation results
        f.write("\nCLUSTER VALIDATION:\n")
        for cluster in sorted(adata.obs['leiden'].unique()):
            cluster_size = (adata.obs['leiden'] == cluster).sum()
            annotation = cluster_annotations[cluster]
            f.write(f"• Cluster {cluster}: {annotation} ({cluster_size} cells)\n")

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

