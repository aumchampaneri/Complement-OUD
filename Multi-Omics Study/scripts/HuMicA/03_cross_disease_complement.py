"""
HuMicA Cross-Disease Complement Comparison Suite
==============================================

Option B: Cross-Disease Complement Comparison
- Statistical analysis of complement differences across diseases
- Volcano plots and differential expression analysis
- Disease classification based on complement profiles
- Advanced statistical modeling and machine learning

Author: Generated for Complement-OUD project
Date: 2025
"""

import os
import sys
import logging
from pathlib import Path
import numpy as np
import pandas as pd
import argparse
import warnings
warnings.filterwarnings('ignore')

# Set non-interactive matplotlib backend
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# Configure matplotlib for optimal plotting
matplotlib.rcParams['font.size'] = 12
matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['savefig.dpi'] = 300
matplotlib.rcParams['savefig.bbox'] = 'tight'
matplotlib.rcParams['axes.titlesize'] = 14
matplotlib.rcParams['axes.labelsize'] = 12
matplotlib.rcParams['legend.fontsize'] = 10

import scanpy as sc
from datetime import datetime
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.metrics import classification_report, confusion_matrix
from statsmodels.stats.multitest import multipletests
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('cross_disease_complement.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

class CrossDiseaseConfig:
    """Configuration for cross-disease complement analysis"""

    # Paths
    DATA_CHECKPOINT = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study/results/scanvi/01_scanvi/adata_scanvi_processed.h5ad"
    COMPLEMENT_DATA = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study/results/scanvi/02_complement/adata_with_complement_scores.h5ad"
    OUTPUT_DIR = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study/results/scanvi/03_cross_disease"

    # Analysis parameters
    FIGSIZE = (12, 8)
    DPI = 300
    ALPHA = 0.05  # Significance threshold
    FDR_METHOD = 'fdr_bh'  # Benjamini-Hochberg FDR correction

    # Disease groups for analysis
    PRIMARY_DISEASES = ['No Neuropathology', 'AD', 'PD', 'MS', 'Epilepsy']
    SECONDARY_DISEASES = ['LBD', 'COVID-19', 'ASD']
    ALL_DISEASES = PRIMARY_DISEASES + SECONDARY_DISEASES

    # Complement pathways (should match 02_complement_analysis.py)
    COMPLEMENT_PATHWAYS = [
        'complement_classical_pathway',
        'complement_alternative_pathway',
        'complement_lectin_pathway',
        'complement_terminal_pathway',
        'complement_regulators',
        'complement_anaphylatoxins_receptors',
        'complement_overall'
    ]

    # Key complement genes for detailed analysis
    KEY_COMPLEMENT_GENES = [
        'C1QA', 'C1QB', 'C1QC',  # Classical pathway initiators
        'C3', 'C5',              # Central components
        'C3AR1', 'C5AR1',        # Anaphylatoxin receptors
        'CD55', 'CD59',          # Regulators
        'CFH', 'CFI'             # Alternative pathway regulators
    ]

    # Cell types of interest
    KEY_CELL_TYPES = [
        'Homeostatic Microglia',
        'Disease-Associated Microglia (DAM)',
        'Activated Microglia I',
        'Activated Microglia II',
        'Phagocytic Microglia'
    ]

    # Statistical thresholds
    MIN_CELLS_PER_GROUP = 50
    LOG2FC_THRESHOLD = 0.5
    PVALUE_THRESHOLD = 0.05

def setup_environment():
    """Setup analysis environment"""
    logger.info("Setting up cross-disease complement analysis environment...")
    logger.info(f"Matplotlib backend: {matplotlib.get_backend()} (non-interactive)")

    # Create output directory
    os.makedirs(CrossDiseaseConfig.OUTPUT_DIR, exist_ok=True)

    # Configure scanpy
    sc.settings.verbosity = 2
    sc.settings.set_figure_params(dpi=CrossDiseaseConfig.DPI, facecolor='white')
    sc.settings.autoshow = False
    sc.settings.figdir = CrossDiseaseConfig.OUTPUT_DIR

    logger.info(f"Output directory: {CrossDiseaseConfig.OUTPUT_DIR}")
    logger.info(f"Analysis focus: {len(CrossDiseaseConfig.PRIMARY_DISEASES)} primary diseases")

def load_complement_data():
    """Load processed data with complement scores"""
    logger.info("Loading HuMicA data with complement pathway scores...")

    # First try to load complement-scored data from Option A
    if os.path.exists(CrossDiseaseConfig.COMPLEMENT_DATA):
        logger.info("Loading complement-scored data from Option A...")
        adata = sc.read_h5ad(CrossDiseaseConfig.COMPLEMENT_DATA)
        logger.info(f"✅ Loaded complement-scored data: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    elif os.path.exists(CrossDiseaseConfig.DATA_CHECKPOINT):
        logger.info("Loading original checkpoint data...")
        adata = sc.read_h5ad(CrossDiseaseConfig.DATA_CHECKPOINT)
        logger.info(f"✅ Loaded data: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    else:
        raise FileNotFoundError(f"No data checkpoint found")

    # Check if complement scores are available
    complement_scores = [col for col in adata.obs.columns if col.startswith('complement_')]

    if not complement_scores:
        logger.warning("No complement scores found! Run 02_complement_analysis.py first")
        logger.info("Calculating basic complement scores...")
        # Basic fallback scoring
        from scipy.sparse import issparse

        key_genes = ['C1QA', 'C1QB', 'C1QC', 'C3', 'C5AR1']
        available_genes = [gene for gene in key_genes if gene in adata.var_names]

        if available_genes:
            sc.tl.score_genes(adata, available_genes, score_name='complement_overall')
            complement_scores = ['complement_overall']
            logger.info(f"Created basic complement score with {len(available_genes)} genes")

    logger.info(f"Available complement scores: {len(complement_scores)}")
    for score in complement_scores:
        logger.info(f"   • {score}")

    # Filter for diseases with sufficient cells
    logger.info("Filtering diseases by cell count...")
    disease_counts = adata.obs['Group'].value_counts()

    valid_diseases = []
    for disease in CrossDiseaseConfig.ALL_DISEASES:
        if disease in disease_counts and disease_counts[disease] >= CrossDiseaseConfig.MIN_CELLS_PER_GROUP:
            valid_diseases.append(disease)
            logger.info(f"   ✅ {disease}: {disease_counts[disease]:,} cells")
        elif disease in disease_counts:
            logger.info(f"   ⚠️  {disease}: {disease_counts[disease]:,} cells (too few, skipping)")

    # Filter data to valid diseases
    adata = adata[adata.obs['Group'].isin(valid_diseases)].copy()
    logger.info(f"Final dataset: {adata.n_obs:,} cells from {len(valid_diseases)} diseases")

    return adata, valid_diseases

def perform_differential_complement_analysis(adata, valid_diseases):
    """Perform differential expression analysis for complement genes across diseases"""
    logger.info("Performing differential complement gene analysis...")

    # Get available complement genes
    available_complement_genes = [gene for gene in CrossDiseaseConfig.KEY_COMPLEMENT_GENES
                                 if gene in adata.var_names]

    if not available_complement_genes:
        logger.warning("No key complement genes found in dataset")
        return None

    logger.info(f"Analyzing {len(available_complement_genes)} complement genes")

    # Prepare results storage
    all_results = []

    # Control group
    control_group = 'No Neuropathology'
    if control_group not in valid_diseases:
        logger.warning(f"Control group '{control_group}' not available")
        return None

    control_mask = adata.obs['Group'] == control_group

    # Compare each disease to control
    for disease in valid_diseases:
        if disease == control_group:
            continue

        logger.info(f"   Analyzing {disease} vs {control_group}...")
        disease_mask = adata.obs['Group'] == disease

        for gene in available_complement_genes:
            # Get expression data with proper indexing
            gene_idx = list(adata.var_names).index(gene)

            # Convert boolean mask to numpy array for proper indexing
            control_indices = np.where(control_mask)[0]
            disease_indices = np.where(disease_mask)[0]

            if hasattr(adata.X, 'toarray'):
                # For sparse matrices
                control_expr = adata.X[control_indices, :][:, gene_idx].toarray().flatten()
                disease_expr = adata.X[disease_indices, :][:, gene_idx].toarray().flatten()
            else:
                # For dense matrices
                control_expr = adata.X[control_indices, gene_idx]
                disease_expr = adata.X[disease_indices, gene_idx]

            # Statistical test
            try:
                stat, pval = stats.mannwhitneyu(disease_expr, control_expr, alternative='two-sided')

                # Calculate fold change (add pseudocount to avoid division by zero)
                control_mean = np.mean(control_expr) + 1e-8
                disease_mean = np.mean(disease_expr) + 1e-8
                fold_change = disease_mean / control_mean
                log2fc = np.log2(fold_change)

                all_results.append({
                    'Gene': gene,
                    'Disease': disease,
                    'Control_Mean': np.mean(control_expr),
                    'Disease_Mean': np.mean(disease_expr),
                    'Fold_Change': fold_change,
                    'Log2FC': log2fc,
                    'P_Value': pval,
                    'Control_N': len(control_expr),
                    'Disease_N': len(disease_expr)
                })

            except Exception as e:
                logger.warning(f"Error analyzing {gene} in {disease}: {e}")
                continue

    if not all_results:
        logger.warning("No differential expression results generated")
        return None

    # Convert to DataFrame and apply multiple testing correction
    results_df = pd.DataFrame(all_results)

    # FDR correction within each disease
    for disease in results_df['Disease'].unique():
        disease_mask = results_df['Disease'] == disease
        pvals = results_df.loc[disease_mask, 'P_Value'].values

        if len(pvals) > 0:
            rejected, pvals_corrected, _, _ = multipletests(pvals, method=CrossDiseaseConfig.FDR_METHOD)
            results_df.loc[disease_mask, 'FDR_P_Value'] = pvals_corrected
            results_df.loc[disease_mask, 'Significant'] = rejected

    # Add significance categories
    results_df['Significance_Category'] = 'Not Significant'
    sig_mask = (results_df['FDR_P_Value'] < CrossDiseaseConfig.PVALUE_THRESHOLD)
    results_df.loc[sig_mask & (results_df['Log2FC'] > CrossDiseaseConfig.LOG2FC_THRESHOLD), 'Significance_Category'] = 'Upregulated'
    results_df.loc[sig_mask & (results_df['Log2FC'] < -CrossDiseaseConfig.LOG2FC_THRESHOLD), 'Significance_Category'] = 'Downregulated'

    # Save results
    results_df.to_csv(f"{CrossDiseaseConfig.OUTPUT_DIR}/complement_differential_expression.csv", index=False)

    logger.info(f"Differential expression analysis completed:")
    logger.info(f"   • Total comparisons: {len(results_df)}")
    logger.info(f"   • Significant (FDR < {CrossDiseaseConfig.PVALUE_THRESHOLD}): {results_df['Significant'].sum()}")
    logger.info(f"   • Upregulated: {(results_df['Significance_Category'] == 'Upregulated').sum()}")
    logger.info(f"   • Downregulated: {(results_df['Significance_Category'] == 'Downregulated').sum()}")

    return results_df

def create_volcano_plots(results_df):
    """Create volcano plots for complement gene differential expression"""
    logger.info("Creating volcano plots...")

    if results_df is None or results_df.empty:
        logger.warning("No results available for volcano plots")
        return

    diseases = results_df['Disease'].unique()
    n_diseases = len(diseases)

    # Multi-panel volcano plot
    n_cols = 3
    n_rows = (n_diseases + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    if n_rows == 1:
        axes = axes.reshape(1, -1)

    for i, disease in enumerate(diseases):
        row, col = i // n_cols, i % n_cols
        ax = axes[row, col]

        disease_data = results_df[results_df['Disease'] == disease]

        # Create volcano plot
        x = disease_data['Log2FC']
        y = -np.log10(disease_data['FDR_P_Value'])

        # Color by significance
        colors = []
        for _, row in disease_data.iterrows():
            if row['Significance_Category'] == 'Upregulated':
                colors.append('red')
            elif row['Significance_Category'] == 'Downregulated':
                colors.append('blue')
            else:
                colors.append('gray')

        ax.scatter(x, y, c=colors, alpha=0.7, s=50)

        # Add significance thresholds
        ax.axhline(-np.log10(CrossDiseaseConfig.PVALUE_THRESHOLD), color='black', linestyle='--', alpha=0.5)
        ax.axvline(CrossDiseaseConfig.LOG2FC_THRESHOLD, color='black', linestyle='--', alpha=0.5)
        ax.axvline(-CrossDiseaseConfig.LOG2FC_THRESHOLD, color='black', linestyle='--', alpha=0.5)

        # Label significant genes
        sig_genes = disease_data[disease_data['Significant']]
        for _, gene_row in sig_genes.iterrows():
            ax.annotate(gene_row['Gene'],
                       (gene_row['Log2FC'], -np.log10(gene_row['FDR_P_Value'])),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=8, alpha=0.8)

        ax.set_xlabel('Log2 Fold Change')
        ax.set_ylabel('-Log10 FDR P-Value')
        ax.set_title(f'{disease} vs Control')
        ax.grid(True, alpha=0.3)

    # Hide unused subplots
    for i in range(n_diseases, n_rows * n_cols):
        row, col = i // n_cols, i % n_cols
        axes[row, col].set_visible(False)

    plt.suptitle('Complement Gene Differential Expression Across Diseases', fontsize=16, y=0.98)
    plt.tight_layout()
    plt.subplots_adjust(top=0.93)
    plt.savefig(f"{CrossDiseaseConfig.OUTPUT_DIR}/complement_volcano_plots.png",
                dpi=CrossDiseaseConfig.DPI, bbox_inches='tight')
    plt.close()

    # Individual high-quality volcano plots for key diseases
    key_diseases = ['AD', 'PD', 'MS']

    for disease in key_diseases:
        if disease in diseases:
            plt.figure(figsize=(10, 8))

            disease_data = results_df[results_df['Disease'] == disease]

            # Create scatter plot with better styling
            upregulated = disease_data[disease_data['Significance_Category'] == 'Upregulated']
            downregulated = disease_data[disease_data['Significance_Category'] == 'Downregulated']
            not_sig = disease_data[disease_data['Significance_Category'] == 'Not Significant']

            plt.scatter(not_sig['Log2FC'], -np.log10(not_sig['FDR_P_Value']),
                       c='lightgray', alpha=0.6, s=50, label='Not Significant')
            plt.scatter(upregulated['Log2FC'], -np.log10(upregulated['FDR_P_Value']),
                       c='red', alpha=0.8, s=60, label='Upregulated')
            plt.scatter(downregulated['Log2FC'], -np.log10(downregulated['FDR_P_Value']),
                       c='blue', alpha=0.8, s=60, label='Downregulated')

            # Add threshold lines
            plt.axhline(-np.log10(CrossDiseaseConfig.PVALUE_THRESHOLD), color='black', linestyle='--', alpha=0.5)
            plt.axvline(CrossDiseaseConfig.LOG2FC_THRESHOLD, color='black', linestyle='--', alpha=0.5)
            plt.axvline(-CrossDiseaseConfig.LOG2FC_THRESHOLD, color='black', linestyle='--', alpha=0.5)

            # Label genes
            for _, gene_row in pd.concat([upregulated, downregulated]).iterrows():
                plt.annotate(gene_row['Gene'],
                           (gene_row['Log2FC'], -np.log10(gene_row['FDR_P_Value'])),
                           xytext=(5, 5), textcoords='offset points',
                           fontsize=10, fontweight='bold')

            plt.xlabel('Log2 Fold Change vs Control')
            plt.ylabel('-Log10 FDR Adjusted P-Value')
            plt.title(f'Complement Gene Expression: {disease} vs No Neuropathology')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(f"{CrossDiseaseConfig.OUTPUT_DIR}/volcano_{disease.lower()}_detailed.png",
                        dpi=CrossDiseaseConfig.DPI, bbox_inches='tight')
            plt.close()

def analyze_pathway_clustering(adata, valid_diseases):
    """Perform hierarchical clustering of diseases based on complement pathway scores"""
    logger.info("Performing pathway-based disease clustering...")

    # Get complement pathway scores
    pathway_scores = [col for col in adata.obs.columns if col.startswith('complement_')]

    if not pathway_scores:
        logger.warning("No pathway scores found for clustering")
        return

    # Calculate mean pathway scores per disease
    disease_pathway_matrix = []
    disease_names = []

    for disease in valid_diseases:
        disease_mask = adata.obs['Group'] == disease
        if disease_mask.sum() < CrossDiseaseConfig.MIN_CELLS_PER_GROUP:
            continue

        disease_scores = adata.obs.loc[disease_mask, pathway_scores].mean().values
        disease_pathway_matrix.append(disease_scores)
        disease_names.append(disease)

    if len(disease_pathway_matrix) < 2:
        logger.warning("Insufficient diseases for clustering")
        return

    # Create clustering matrix
    pathway_matrix = np.array(disease_pathway_matrix)
    pathway_df = pd.DataFrame(pathway_matrix,
                             index=disease_names,
                             columns=[col.replace('complement_', '').replace('_', ' ').title()
                                    for col in pathway_scores])

    # Standardize data
    scaler = StandardScaler()
    pathway_matrix_scaled = scaler.fit_transform(pathway_matrix)

    # Hierarchical clustering
    linkage_matrix = linkage(pathway_matrix_scaled, method='ward')

    # Create dendrogram
    plt.figure(figsize=(12, 8))

    dendr = dendrogram(linkage_matrix, labels=disease_names, orientation='top')
    plt.title('Disease Clustering Based on Complement Pathway Activity')
    plt.xlabel('Diseases')
    plt.ylabel('Distance')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f"{CrossDiseaseConfig.OUTPUT_DIR}/disease_complement_clustering.png",
                dpi=CrossDiseaseConfig.DPI, bbox_inches='tight')
    plt.close()

    # Create heatmap with clustering
    plt.figure(figsize=(12, 8))

    # Get cluster order from dendrogram
    cluster_order = dendr['ivl']
    reordered_df = pathway_df.loc[cluster_order]

    sns.clustermap(reordered_df, method='ward', cmap='RdBu_r', center=0,
                   figsize=(12, 8), annot=True, fmt='.3f',
                   cbar_kws={'label': 'Mean Pathway Score'})
    plt.title('Complement Pathway Activity Heatmap (Clustered)')
    plt.tight_layout()
    plt.savefig(f"{CrossDiseaseConfig.OUTPUT_DIR}/complement_pathway_clustered_heatmap.png",
                dpi=CrossDiseaseConfig.DPI, bbox_inches='tight')
    plt.close()

    # PCA analysis
    pca = PCA(n_components=min(3, len(disease_names)-1))
    pca_result = pca.fit_transform(pathway_matrix_scaled)

    # PCA plot
    plt.figure(figsize=(10, 8))

    colors = plt.cm.Set1(np.linspace(0, 1, len(disease_names)))
    for i, disease in enumerate(disease_names):
        plt.scatter(pca_result[i, 0], pca_result[i, 1],
                   c=[colors[i]], s=100, alpha=0.8, label=disease)
        plt.annotate(disease, (pca_result[i, 0], pca_result[i, 1]),
                    xytext=(5, 5), textcoords='offset points', fontsize=10)

    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
    plt.title('Disease Clustering: Complement Pathway PCA')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{CrossDiseaseConfig.OUTPUT_DIR}/complement_pathway_pca.png",
                dpi=CrossDiseaseConfig.DPI, bbox_inches='tight')
    plt.close()

    # Save clustering results
    cluster_results = pd.DataFrame({
        'Disease': disease_names,
        'PC1': pca_result[:, 0],
        'PC2': pca_result[:, 1] if pca_result.shape[1] > 1 else [0]*len(disease_names)
    })
    cluster_results.to_csv(f"{CrossDiseaseConfig.OUTPUT_DIR}/disease_clustering_results.csv", index=False)

    logger.info(f"Clustering analysis completed for {len(disease_names)} diseases")
    logger.info(f"PCA explained variance: PC1={pca.explained_variance_ratio_[0]:.1%}, PC2={pca.explained_variance_ratio_[1]:.1%}")

def perform_disease_classification(adata, valid_diseases):
    """Perform machine learning classification of diseases based on complement profiles"""
    logger.info("Performing disease classification analysis...")

    # Get complement pathway scores
    pathway_scores = [col for col in adata.obs.columns if col.startswith('complement_')]

    if len(pathway_scores) < 2:
        logger.warning("Insufficient pathway scores for classification")
        return

    # Prepare data
    X = adata.obs[pathway_scores].values
    y = adata.obs['Group'].values

    # Filter to valid diseases
    valid_mask = np.isin(y, valid_diseases)
    X = X[valid_mask]
    y = y[valid_mask]

    # Remove diseases with too few samples
    disease_counts = pd.Series(y).value_counts()
    sufficient_diseases = disease_counts[disease_counts >= CrossDiseaseConfig.MIN_CELLS_PER_GROUP].index

    sufficient_mask = np.isin(y, sufficient_diseases)
    X = X[sufficient_mask]
    y = y[sufficient_mask]

    if len(sufficient_diseases) < 2:
        logger.warning("Insufficient diseases with enough samples for classification")
        return

    logger.info(f"Classification dataset: {X.shape[0]} cells, {len(sufficient_diseases)} diseases")

    # Random Forest Classification
    rf_classifier = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1)

    # Cross-validation
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    cv_scores = cross_val_score(rf_classifier, X, y, cv=cv, scoring='accuracy')

    logger.info(f"Cross-validation accuracy: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")

    # Train full model for feature importance
    rf_classifier.fit(X, y)

    # Feature importance analysis
    feature_importance = pd.DataFrame({
        'Pathway': [col.replace('complement_', '').replace('_', ' ').title() for col in pathway_scores],
        'Importance': rf_classifier.feature_importances_
    }).sort_values('Importance', ascending=False)

    # Plot feature importance
    plt.figure(figsize=(10, 6))
    sns.barplot(data=feature_importance, x='Importance', y='Pathway')
    plt.title('Complement Pathway Importance for Disease Classification')
    plt.xlabel('Feature Importance')
    plt.tight_layout()
    plt.savefig(f"{CrossDiseaseConfig.OUTPUT_DIR}/disease_classification_feature_importance.png",
                dpi=CrossDiseaseConfig.DPI, bbox_inches='tight')
    plt.close()

    # Save classification results
    classification_results = {
        'cv_accuracy_mean': cv_scores.mean(),
        'cv_accuracy_std': cv_scores.std(),
        'cv_scores': cv_scores.tolist(),
        'feature_importance': feature_importance.to_dict('records'),
        'diseases_analyzed': sufficient_diseases.tolist(),
        'n_samples': X.shape[0],
        'n_features': X.shape[1]
    }

    import json
    with open(f"{CrossDiseaseConfig.OUTPUT_DIR}/disease_classification_results.json", "w") as f:
        json.dump(classification_results, f, indent=2)

    feature_importance.to_csv(f"{CrossDiseaseConfig.OUTPUT_DIR}/pathway_feature_importance.csv", index=False)

    logger.info("Disease classification completed:")
    logger.info(f"   • Accuracy: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
    logger.info(f"   • Top pathway: {feature_importance.iloc[0]['Pathway']} ({feature_importance.iloc[0]['Importance']:.3f})")

def create_disease_comparison_summary(adata, valid_diseases, results_df):
    """Create comprehensive summary of disease comparisons"""
    logger.info("Creating disease comparison summary...")

    # Disease similarity matrix based on complement profiles
    pathway_scores = [col for col in adata.obs.columns if col.startswith('complement_')]

    if pathway_scores:
        # Calculate pairwise correlations between diseases
        disease_correlations = []
        disease_pairs = []

        for i, disease1 in enumerate(valid_diseases):
            for j, disease2 in enumerate(valid_diseases):
                if i <= j:
                    continue

                mask1 = adata.obs['Group'] == disease1
                mask2 = adata.obs['Group'] == disease2

                if mask1.sum() < 20 or mask2.sum() < 20:
                    continue

                # Calculate correlation of pathway scores
                scores1 = adata.obs.loc[mask1, pathway_scores].mean()
                scores2 = adata.obs.loc[mask2, pathway_scores].mean()

                correlation = np.corrcoef(scores1, scores2)[0, 1]

                disease_correlations.append(correlation)
                disease_pairs.append(f"{disease1} vs {disease2}")

        # Create correlation summary plot
        if disease_correlations:
            plt.figure(figsize=(12, 8))

            correlation_df = pd.DataFrame({
                'Disease_Pair': disease_pairs,
                'Correlation': disease_correlations
            }).sort_values('Correlation', ascending=False)

            sns.barplot(data=correlation_df, x='Correlation', y='Disease_Pair')
            plt.title('Disease Similarity Based on Complement Pathway Profiles')
            plt.xlabel('Pathway Score Correlation')
            plt.tight_layout()
            plt.savefig(f"{CrossDiseaseConfig.OUTPUT_DIR}/disease_similarity_correlations.png",
                        dpi=CrossDiseaseConfig.DPI, bbox_inches='tight')
            plt.close()

            correlation_df.to_csv(f"{CrossDiseaseConfig.OUTPUT_DIR}/disease_similarity_matrix.csv", index=False)

    # Summary statistics by disease
    if results_df is not None:
        summary_stats = []

        for disease in valid_diseases:
            if disease == 'No Neuropathology':
                continue

            disease_results = results_df[results_df['Disease'] == disease]

            n_significant = disease_results['Significant'].sum()
            n_upregulated = (disease_results['Significance_Category'] == 'Upregulated').sum()
            n_downregulated = (disease_results['Significance_Category'] == 'Downregulated').sum()

            # Get most significant genes
            top_genes = disease_results.nsmallest(3, 'FDR_P_Value')['Gene'].tolist()

            summary_stats.append({
                'Disease': disease,
                'Total_Genes_Tested': len(disease_results),
                'Significant_Genes': n_significant,
                'Upregulated_Genes': n_upregulated,
                'Downregulated_Genes': n_downregulated,
                'Percent_Significant': (n_significant / len(disease_results)) * 100,
                'Top_Significant_Genes': ', '.join(top_genes)
            })

        summary_df = pd.DataFrame(summary_stats)
        summary_df.to_csv(f"{CrossDiseaseConfig.OUTPUT_DIR}/disease_complement_summary.csv", index=False)

        # Plot summary statistics
        fig, axes = plt.subplots(2, 2, figsize=(15, 10))

        # Significant genes per disease
        sns.barplot(data=summary_df, x='Disease', y='Significant_Genes', ax=axes[0,0])
        axes[0,0].set_title('Significant Complement Genes per Disease')
        axes[0,0].tick_params(axis='x', rotation=45)

        # Percentage significant
        sns.barplot(data=summary_df, x='Disease', y='Percent_Significant', ax=axes[0,1])
        axes[0,1].set_title('Percentage of Complement Genes Significantly Altered')
        axes[0,1].tick_params(axis='x', rotation=45)

        # Up vs down regulation
        regulation_data = summary_df.melt(id_vars=['Disease'],
                                        value_vars=['Upregulated_Genes', 'Downregulated_Genes'],
                                        var_name='Regulation', value_name='Count')
        sns.barplot(data=regulation_data, x='Disease', y='Count', hue='Regulation', ax=axes[1,0])
        axes[1,0].set_title('Up vs Down Regulated Complement Genes')
        axes[1,0].tick_params(axis='x', rotation=45)

        # Cell count per disease
        disease_cell_counts = adata.obs['Group'].value_counts()
        disease_cell_counts = disease_cell_counts.reindex(valid_diseases)

        axes[1,1].bar(range(len(disease_cell_counts)), disease_cell_counts.values)
        axes[1,1].set_xticks(range(len(disease_cell_counts)))
        axes[1,1].set_xticklabels(disease_cell_counts.index, rotation=45)
        axes[1,1].set_title('Cell Count per Disease')
        axes[1,1].set_ylabel('Number of Cells')

        plt.tight_layout()
        plt.savefig(f"{CrossDiseaseConfig.OUTPUT_DIR}/disease_comparison_summary.png",
                    dpi=CrossDiseaseConfig.DPI, bbox_inches='tight')
        plt.close()

def generate_comprehensive_report(adata, valid_diseases, results_df):
    """Generate comprehensive analysis report"""
    logger.info("Generating comprehensive cross-disease analysis report...")

    # Calculate summary statistics
    total_cells = adata.n_obs
    total_genes = adata.n_vars

    pathway_scores = [col for col in adata.obs.columns if col.startswith('complement_')]
    n_pathways = len(pathway_scores)

    # Disease-specific statistics
    disease_stats = []
    for disease in valid_diseases:
        mask = adata.obs['Group'] == disease
        n_cells = mask.sum()
        disease_stats.append(f"• {disease}: {n_cells:,} cells")

    # Differential expression summary
    de_summary = ""
    if results_df is not None:
        total_tests = len(results_df)
        significant_tests = results_df['Significant'].sum()
        de_summary = f"""
=== DIFFERENTIAL EXPRESSION SUMMARY ===
• Total comparisons: {total_tests:,}
• Significant results (FDR < 0.05): {significant_tests:,} ({significant_tests/total_tests*100:.1f}%)
• Upregulated findings: {(results_df['Significance_Category'] == 'Upregulated').sum():,}
• Downregulated findings: {(results_df['Significance_Category'] == 'Downregulated').sum():,}

Top Dysregulated Genes by Disease:
"""
        for disease in valid_diseases:
            if disease == 'No Neuropathology':
                continue
            disease_results = results_df[results_df['Disease'] == disease]
            if not disease_results.empty:
                top_gene = disease_results.nsmallest(1, 'FDR_P_Value').iloc[0]
                de_summary += f"• {disease}: {top_gene['Gene']} (Log2FC: {top_gene['Log2FC']:.2f}, FDR: {top_gene['FDR_P_Value']:.2e})\n"

    report = f"""
HuMicA Cross-Disease Complement Comparison Report
===============================================
Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}

=== ANALYSIS OVERVIEW ===
Objective: Compare complement system activation across neurological diseases
Approach: Statistical analysis, machine learning, and pathway-based clustering
Dataset: Human Microglia Atlas (HuMicA) with complement pathway scores

=== DATASET CHARACTERISTICS ===
• Total cells analyzed: {total_cells:,}
• Total genes: {total_genes:,}
• Complement pathways scored: {n_pathways}
• Diseases analyzed: {len(valid_diseases)}

Disease Breakdown:
{chr(10).join(disease_stats)}

=== COMPLEMENT PATHWAYS ANALYZED ===
1. Classical Pathway - C1q-mediated activation (synaptic pruning)
2. Alternative Pathway - Spontaneous amplification loop
3. Lectin Pathway - Carbohydrate pattern recognition
4. Terminal Pathway - Membrane attack complex formation
5. Complement Regulators - CD55, CD59, CFH system
6. Anaphylatoxin Receptors - C3AR1, C5AR1 signaling
7. Overall Complement Activity - Integrated score

{de_summary}

=== KEY ANALYTICAL APPROACHES ===
1. Differential Expression Analysis
   - Mann-Whitney U tests with FDR correction
   - Volcano plots for visualization
   - Log2 fold change thresholds

2. Disease Clustering Analysis
   - Hierarchical clustering based on pathway scores
   - Principal component analysis (PCA)
   - Disease similarity correlations

3. Machine Learning Classification
   - Random Forest classification
   - Cross-validated accuracy assessment
   - Feature importance ranking

=== OUTPUT FILES ===
Statistical Results:
• complement_differential_expression.csv - Gene-level DE results
• disease_clustering_results.csv - Clustering analysis results
• disease_classification_results.json - ML classification metrics
• pathway_feature_importance.csv - Pathway importance for classification
• disease_similarity_matrix.csv - Cross-disease correlations
• disease_complement_summary.csv - Summary statistics

Visualizations:
• complement_volcano_plots.png - Multi-disease volcano plots
• volcano_[disease]_detailed.png - Individual disease volcano plots
• disease_complement_clustering.png - Hierarchical clustering
• complement_pathway_clustered_heatmap.png - Clustered heatmap
• complement_pathway_pca.png - PCA analysis
• disease_classification_feature_importance.png - ML feature importance
• disease_similarity_correlations.png - Disease similarity
• disease_comparison_summary.png - Multi-panel summary

=== RESEARCH IMPLICATIONS ===
1. Disease-Specific Complement Signatures
   - Unique complement activation patterns per disease
   - Shared vs disease-specific pathway dysregulation
   - Therapeutic target prioritization

2. Mechanistic Insights
   - Classical pathway (C1q) role in synaptic elimination
   - Alternative pathway amplification in chronic inflammation
   - Complement regulator dysfunction patterns

3. Clinical Translation Potential
   - Biomarker development opportunities
   - Therapeutic intervention strategies
   - Patient stratification approaches

=== COMPLEMENT-OUD RESEARCH CONNECTIONS ===
• Neuroinflammation overlap between addiction and neurodegenerative diseases
• C5AR1 as potential target for addiction-related neuroinflammation
• Complement-mediated synaptic pruning in reward circuitry
• Shared complement dysregulation patterns across CNS diseases

=== NEXT STEPS FOR RESEARCH ===
1. Validate key findings in independent cohorts
2. Integrate with OUD complement data
3. Functional validation of top dysregulated genes
4. Drug target prioritization and screening
5. Biomarker validation studies
6. Cross-species validation in animal models

=== METHODOLOGICAL NOTES ===
• Statistical significance: FDR < 0.05
• Fold change threshold: |Log2FC| > 0.5
• Minimum cells per group: {CrossDiseaseConfig.MIN_CELLS_PER_GROUP}
• Multiple testing correction: Benjamini-Hochberg FDR
• Machine learning: 5-fold cross-validation

This analysis provides a comprehensive foundation for understanding
complement system dysfunction across neurological diseases and its
potential relevance to addiction neurobiology.
"""

    # Save report
    with open(f"{CrossDiseaseConfig.OUTPUT_DIR}/cross_disease_complement_report.txt", "w") as f:
        f.write(report)

    print(report)

def main():
    """Main analysis pipeline for cross-disease complement comparison"""

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='HuMicA Cross-Disease Complement Comparison Suite',
        epilog="""
Examples:
  %(prog)s                    # Full cross-disease analysis
  %(prog)s --quick            # Skip ML classification
  %(prog)s --diseases AD PD MS # Analyze specific diseases only
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--quick', action='store_true',
                       help='Skip computationally intensive analyses')
    parser.add_argument('--diseases', nargs='+',
                       help='Specific diseases to analyze')
    args = parser.parse_args()

    try:
        logger.info("=" * 70)
        logger.info("🔬 HuMicA CROSS-DISEASE COMPLEMENT COMPARISON SUITE")
        logger.info("=" * 70)

        # Setup
        setup_environment()

        # Load data
        adata, valid_diseases = load_complement_data()

        # Filter diseases if specified
        if args.diseases:
            specified_diseases = [d for d in args.diseases if d in valid_diseases]
            if specified_diseases:
                valid_diseases = specified_diseases
                adata = adata[adata.obs['Group'].isin(valid_diseases)].copy()
                logger.info(f"Analysis restricted to: {valid_diseases}")

        # Differential expression analysis
        results_df = perform_differential_complement_analysis(adata, valid_diseases)

        # Create volcano plots
        if results_df is not None:
            create_volcano_plots(results_df)

        # Pathway clustering analysis
        analyze_pathway_clustering(adata, valid_diseases)

        # Disease classification (unless quick mode)
        if not args.quick:
            perform_disease_classification(adata, valid_diseases)
        else:
            logger.info("⏭️  Skipping ML classification (--quick mode)")

        # Create comprehensive comparison summary
        create_disease_comparison_summary(adata, valid_diseases, results_df)

        # Generate report
        generate_comprehensive_report(adata, valid_diseases, results_df)

        logger.info("=" * 70)
        logger.info("🎉 Cross-disease complement analysis completed successfully!")
        logger.info("=" * 70)
        logger.info("📊 Check the results directory for comprehensive outputs")
        logger.info("🧬 Ready for integration with OUD complement research")

    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        raise

if __name__ == "__main__":
    main()
