#!/usr/bin/env python3
"""
LEMUR Latent Space Analysis - Comprehensive Embedding Diagnostics
================================================================================

This script performs comprehensive latent space diagnostics for LEMUR embeddings,
including variance analysis, covariate correlations, and visualization overlays.

Key Features:
1. Variance Explained Analysis - Component-wise variance decomposition
2. Embedding-Covariate Correlations - Statistical associations with metadata
3. Latent Space Visualization - UMAP and 2D scatter plots with overlays
4. Robust error handling and professional visualizations

Input: LEMUR results from 03b_LEMUR.py analysis
Output: Diagnostic plots and statistical summaries

Dataset: GSE225158 - Human striatal single-cell RNA-seq (OUD vs Control)
Author: Research Team
Date: June 2025
Version: 1.0
"""

import warnings
import logging
from pathlib import Path
from datetime import datetime
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)

# Configure plotting
plt.style.use('default')
sns.set_palette("husl")
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=300, facecolor='white')

# ================================================================================
# 📁 CONFIGURATION
# ================================================================================

class LatentSpaceConfig:
    """Configuration for latent space analysis"""

    def __init__(self):
        # Paths
        self.BASE_DIR = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study"
        self.LEMUR_RESULTS_DIR = f"{self.BASE_DIR}/results/snrna_scvi/lemur_comprehensive"
        self.OUTPUT_DIR = f"{self.BASE_DIR}/results/snrna_scvi/lemur_latent_space"
        self.PLOTS_DIR = f"{self.OUTPUT_DIR}/plots/latent_structure"

        # Required metadata columns (with flexible cell type naming)
        self.REQUIRED_METADATA = ['Sex', 'level1']
        self.CELLTYPE_CANDIDATES = ['celltype', 'celltype3', 'cell_type', 'celltype1', 'celltype2']

        # Plotting parameters
        self.FIGURE_SIZE = (10, 8)
        self.DPI = 300
        self.PALETTE = 'viridis'

        # Statistical parameters
        self.CORRELATION_METHOD = 'spearman'
        self.SIGNIFICANCE_THRESHOLD = 0.05

    def setup_directories(self):
        """Create output directories"""
        dirs_to_create = [
            self.OUTPUT_DIR,
            self.PLOTS_DIR,
            f"{self.OUTPUT_DIR}/tables"
        ]

        for directory in dirs_to_create:
            Path(directory).mkdir(parents=True, exist_ok=True)

        logger.info(f"📁 Output directories created: {self.OUTPUT_DIR}")


# ================================================================================
# 🔧 UTILITY FUNCTIONS
# ================================================================================

def load_lemur_results(results_dir):
    """Load LEMUR results and associated data"""

    logger.info("📊 Loading LEMUR results...")

    try:
        # Load the processed AnnData object - prioritize LEMUR-processed files
        h5ad_file = None

        # Try LEMUR-processed files first (these have embeddings)
        lemur_paths = [
            f"{results_dir}/lemur_processed_data.h5ad",
            "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study/data/processed/snrna_scvi/GSE225158_lemur_processed.h5ad",
            f"{results_dir}/../lemur_comprehensive/lemur_processed_data.h5ad"
        ]

        for path in lemur_paths:
            if Path(path).exists():
                h5ad_file = Path(path)
                logger.info(f"   🌊 Found LEMUR-processed file: {path}")
                break

        # If no LEMUR files found, try other locations
        if h5ad_file is None:
            fallback_paths = [
                Path(results_dir).parent / "lemur_analysis" / "processed_data.h5ad",
                f"{results_dir}/processed_data.h5ad",
                f"{results_dir}/../processed_data.h5ad",
                "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study/data/processed/snrna_scvi/GSE225158_annotated_scvi.h5ad"
            ]

            for path in fallback_paths:
                if Path(path).exists():
                    h5ad_file = Path(path)
                    logger.info(f"   📊 Found fallback file: {path}")
                    break

        if not h5ad_file.exists():
            raise FileNotFoundError("Could not find AnnData file with LEMUR embeddings")

        adata = sc.read_h5ad(str(h5ad_file))
        logger.info(f"   ✅ Loaded AnnData: {adata.n_obs:,} cells, {adata.n_vars:,} genes")

        # Check for LEMUR embeddings
        lemur_embeddings = None
        embedding_key = None

        # Check for LEMUR embeddings in priority order
        embedding_candidates = ['X_lemur_harmony', 'X_lemur', 'lemur_embedding']
        for key in embedding_candidates:
            if key in adata.obsm:
                lemur_embeddings = adata.obsm[key]
                embedding_key = key
                logger.info(f"   ✅ Found LEMUR embeddings: {lemur_embeddings.shape} ({embedding_key})")
                break

        if lemur_embeddings is None:
            available_embeddings = list(adata.obsm.keys())
            raise ValueError(f"No LEMUR embeddings found in AnnData.obsm. Available: {available_embeddings}")

        # Create lemur_result dictionary
        lemur_result = {
            'embedding': lemur_embeddings,
            'embedding_key': embedding_key,
            'n_components': lemur_embeddings.shape[1],
            'n_cells': lemur_embeddings.shape[0]
        }

        return lemur_result, adata

    except Exception as e:
        logger.error(f"❌ Failed to load LEMUR results: {e}")
        return None, None


def prepare_metadata_for_analysis(adata, required_columns, celltype_candidates):
    """Prepare and validate metadata for analysis with proper type consistency"""

    logger.info("🔍 Preparing metadata...")

    # Check for required columns
    missing_columns = [col for col in required_columns if col not in adata.obs.columns]
    if missing_columns:
        logger.warning(f"⚠️ Missing required metadata columns: {missing_columns}")
        # Use available columns only
        available_columns = [col for col in required_columns if col in adata.obs.columns]
        if not available_columns:
            raise ValueError("No required metadata columns found")
    else:
        available_columns = required_columns.copy()

    # Find cell type column
    celltype_col = None
    for candidate in celltype_candidates:
        if candidate in adata.obs.columns:
            celltype_col = candidate
            available_columns.append(candidate)
            logger.info(f"   ✅ Using '{candidate}' as cell type column")
            break

    if celltype_col is None:
        logger.warning(f"⚠️ No cell type column found. Tried: {celltype_candidates}")
        logger.warning("   Available columns with 'cell' or 'type': " +
                      str([col for col in adata.obs.columns if 'cell' in col.lower() or 'type' in col.lower()]))

    metadata_df = adata.obs[available_columns].copy()

    # Handle missing values and ensure proper categorical typing
    for col in available_columns:
        if metadata_df[col].isna().any():
            logger.warning(f"⚠️ Column '{col}' has missing values - filling with 'Unknown'")
            metadata_df[col] = metadata_df[col].fillna('Unknown')

        # Ensure categorical variables are properly typed for downstream analysis
        if metadata_df[col].dtype == 'object' or not pd.api.types.is_numeric_dtype(metadata_df[col]):
            # Convert to categorical type first for consistent handling
            metadata_df[col] = metadata_df[col].astype('category')
            logger.info(f"   🔧 Converted '{col}' to categorical type")

    logger.info(f"   ✅ Prepared metadata for {len(available_columns)} columns")
    for col in available_columns:
        unique_vals = metadata_df[col].nunique()
        dtype_info = metadata_df[col].dtype
        logger.info(f"      {col}: {unique_vals} unique values, dtype: {dtype_info}")

    return metadata_df, available_columns


def one_hot_encode_metadata(metadata_df, columns):
    """One-hot encode categorical metadata with proper type handling"""

    logger.info("🔢 One-hot encoding metadata...")

    encoded_data = {}
    feature_names = []

    for col in columns:
        # Ensure column is properly typed as categorical before operations
        if not pd.api.types.is_categorical_dtype(metadata_df[col]):
            metadata_df[col] = metadata_df[col].astype('category')
            logger.debug(f"   🔧 Converted {col} to categorical for encoding")

        # Get unique values from categorical data
        if pd.api.types.is_categorical_dtype(metadata_df[col]):
            unique_vals = sorted(metadata_df[col].cat.categories.tolist())
        else:
            unique_vals = sorted(metadata_df[col].unique())

        # Create one-hot encoding with proper type handling
        for val in unique_vals:
            feature_name = f"{col}_{val}"
            # Use categorical comparison which is more robust
            try:
                encoded_data[feature_name] = (metadata_df[col] == val).astype(int)
                feature_names.append(feature_name)
            except Exception as e:
                logger.warning(f"   ⚠️ Error encoding {feature_name}: {e}")
                # Fallback: create zero column
                encoded_data[feature_name] = np.zeros(len(metadata_df), dtype=int)
                feature_names.append(feature_name)

    encoded_df = pd.DataFrame(encoded_data, index=metadata_df.index)

    logger.info(f"   ✅ Created {len(feature_names)} one-hot encoded features")

    return encoded_df, feature_names


# ================================================================================
# 📊 ANALYSIS FUNCTIONS
# ================================================================================

def analyze_variance_explained(lemur_result, config):
    """Compute and visualize variance explained by each LEMUR component"""

    logger.info("📈 Analyzing variance explained...")

    embeddings = lemur_result['embedding']
    n_components = lemur_result['n_components']

    # Compute variance for each component
    component_variances = np.var(embeddings, axis=0)

    # Normalize to sum to 1
    variance_explained = component_variances / np.sum(component_variances)

    # Cumulative variance explained
    cumulative_variance = np.cumsum(variance_explained)

    logger.info(f"   ✅ Top 5 components explain {cumulative_variance[4]:.3f} of total variance")

    # Create visualization
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # Bar plot of variance explained
    components = np.arange(1, n_components + 1)
    bars = ax1.bar(components, variance_explained, alpha=0.7, color='steelblue')
    ax1.set_xlabel('LEMUR Component', fontsize=12)
    ax1.set_ylabel('Variance Explained', fontsize=12)
    ax1.set_title('Variance Explained by LEMUR Components', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)

    # Add value labels on bars for top 5 components
    for i, (bar, var_exp) in enumerate(zip(bars[:5], variance_explained[:5])):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + 0.001,
                f'{var_exp:.3f}', ha='center', va='bottom', fontsize=10)

    # Cumulative variance plot
    ax2.plot(components, cumulative_variance, 'o-', color='darkred', linewidth=2, markersize=4)
    ax2.set_xlabel('Number of Components', fontsize=12)
    ax2.set_ylabel('Cumulative Variance Explained', fontsize=12)
    ax2.set_title('Cumulative Variance Explained', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(0, 1)

    # Add horizontal lines at common thresholds
    for threshold in [0.5, 0.8, 0.9]:
        ax2.axhline(y=threshold, color='gray', linestyle='--', alpha=0.5)
        ax2.text(n_components * 0.7, threshold + 0.02, f'{threshold:.0%}',
                fontsize=10, alpha=0.7)

    plt.tight_layout()

    # Save plot
    variance_plot_file = Path(config.PLOTS_DIR) / "variance_explained_barplot.png"
    plt.savefig(variance_plot_file, dpi=config.DPI, bbox_inches='tight')
    plt.close()

    logger.info(f"   💾 Saved variance plot: {variance_plot_file}")

    # Return results
    variance_results = {
        'variance_explained': variance_explained,
        'cumulative_variance': cumulative_variance,
        'component_variances': component_variances,
        'top_5_cumulative': cumulative_variance[4]
    }

    return variance_results


def compute_embedding_covariate_correlations(lemur_result, metadata_df, encoded_features, config):
    """Compute correlations between LEMUR components and covariates with robust type handling"""

    logger.info("🔗 Computing embedding-covariate correlations...")

    embeddings = lemur_result['embedding']
    n_components = lemur_result['n_components']

    # Initialize correlation matrix
    correlation_matrix = np.zeros((n_components, len(encoded_features)))
    pvalue_matrix = np.zeros((n_components, len(encoded_features)))

    # Compute correlations with robust type checking
    for i in range(n_components):
        embedding_component = embeddings[:, i]

        # Check for valid embedding component
        if np.isnan(embedding_component).any() or np.isinf(embedding_component).any():
            logger.warning(f"   ⚠️ Invalid values in embedding component {i+1}, skipping")
            continue

        for j, feature in enumerate(encoded_features):
            try:
                # Ensure feature data is numeric and handle missing values
                feature_data = metadata_df[feature].values

                # Convert to numeric if needed and handle type issues
                if not pd.api.types.is_numeric_dtype(feature_data):
                    feature_data = pd.to_numeric(feature_data, errors='coerce')

                # Handle missing or infinite values
                if np.isnan(feature_data).any() or np.isinf(feature_data).any():
                    feature_data = np.nan_to_num(feature_data, nan=np.nanmedian(feature_data),
                                                posinf=np.nanmax(feature_data[np.isfinite(feature_data)]),
                                                neginf=np.nanmin(feature_data[np.isfinite(feature_data)]))

                # Check if there's variation in the feature
                if np.var(feature_data) == 0:
                    logger.debug(f"   No variation in feature {feature}, correlation set to 0")
                    correlation_matrix[i, j] = 0.0
                    pvalue_matrix[i, j] = 1.0
                    continue

                # Compute correlation with error handling
                corr, pval = spearmanr(embedding_component, feature_data)

                # Handle correlation computation results
                correlation_matrix[i, j] = corr if not np.isnan(corr) else 0.0
                pvalue_matrix[i, j] = pval if not np.isnan(pval) else 1.0

            except Exception as e:
                logger.warning(f"   ⚠️ Error computing correlation for {feature}: {e}")
                correlation_matrix[i, j] = 0.0
                pvalue_matrix[i, j] = 1.0

    # Create DataFrame for results
    component_names = [f'Component_{i+1}' for i in range(n_components)]
    corr_df = pd.DataFrame(
        correlation_matrix,
        index=component_names,
        columns=encoded_features
    )

    pval_df = pd.DataFrame(
        pvalue_matrix,
        index=component_names,
        columns=encoded_features
    )

    # Find significant correlations
    significant_mask = pval_df < config.SIGNIFICANCE_THRESHOLD
    n_significant = significant_mask.sum().sum()

    logger.info(f"   ✅ Computed {n_components} × {len(encoded_features)} correlations")
    logger.info(f"   📊 {n_significant} significant correlations (p < {config.SIGNIFICANCE_THRESHOLD})")

    # Create heatmap with error handling
    try:
        plt.figure(figsize=(max(12, len(encoded_features) * 0.5), max(8, n_components * 0.4)))

        # Create mask for non-significant correlations
        mask = ~significant_mask.values

        # Plot heatmap
        sns.heatmap(
            corr_df,
            annot=True,
            fmt='.3f',
            cmap='RdBu_r',
            center=0,
            mask=mask,
            cbar_kws={'label': 'Spearman Correlation'},
            square=False,
            linewidths=0.5
        )

        plt.title('LEMUR Component - Covariate Correlations\n(Only significant correlations shown)',
                  fontsize=14, fontweight='bold', pad=20)
        plt.xlabel('Covariates', fontsize=12)
        plt.ylabel('LEMUR Components', fontsize=12)
        plt.xticks(rotation=45, ha='right')
        plt.yticks(rotation=0)

        plt.tight_layout()

        # Save heatmap
        heatmap_file = Path(config.PLOTS_DIR) / "embedding_covariate_correlation_heatmap.png"
        plt.savefig(heatmap_file, dpi=config.DPI, bbox_inches='tight')
        plt.close()

        logger.info(f"   💾 Saved correlation heatmap: {heatmap_file}")

    except Exception as e:
        logger.error(f"   ❌ Error creating correlation heatmap: {e}")

    # Return results
    correlation_results = {
        'correlation_matrix': corr_df,
        'pvalue_matrix': pval_df,
        'significant_correlations': significant_mask,
        'n_significant': n_significant
    }

    return correlation_results


def visualize_latent_space_overlays(lemur_result, adata, metadata_columns, config):
    """Create latent space visualizations with covariate overlays"""

    logger.info("🎨 Creating latent space visualizations...")

    embeddings = lemur_result['embedding']

    # Prepare multiple coordinate systems for comprehensive visualization
    coordinate_systems = []

    # 1. Check for existing scVI UMAP
    if 'X_umap' in adata.obsm:
        logger.info("   📊 Found existing scVI UMAP coordinates")
        coordinate_systems.append({
            'coords': adata.obsm['X_umap'],
            'names': ['UMAP1', 'UMAP2'],
            'prefix': 'scvi_umap',
            'title': 'scVI UMAP',
            'point_size': 1
        })

    # 2. Check for LEMUR-based UMAP (if computed on LEMUR embeddings)
    lemur_embedding_key = lemur_result.get('embedding_key', 'X_lemur')
    if lemur_embedding_key in adata.obsm:
        # Compute UMAP on LEMUR embeddings if not already done
        try:
            logger.info(f"   🌊 Computing UMAP on LEMUR embeddings ({lemur_embedding_key})")
            adata_temp = adata.copy()
            sc.pp.neighbors(adata_temp, use_rep=lemur_embedding_key, key_added='lemur')
            sc.tl.umap(adata_temp, neighbors_key='lemur')

            coordinate_systems.append({
                'coords': adata_temp.obsm['X_umap'],
                'names': ['LEMUR-UMAP1', 'LEMUR-UMAP2'],
                'prefix': 'lemur_umap',
                'title': 'LEMUR UMAP',
                'point_size': 1
            })
            logger.info("   ✅ LEMUR UMAP computed successfully")
        except Exception as e:
            logger.warning(f"   ⚠️ Failed to compute LEMUR UMAP: {e}")

    # 3. Always include first 2 LEMUR components as fallback
    coordinate_systems.append({
        'coords': embeddings[:, :2],
        'names': ['LEMUR Component 1', 'LEMUR Component 2'],
        'prefix': 'lemur_2d',
        'title': 'LEMUR Components',
        'point_size': 3
    })

    logger.info(f"   📊 Will create visualizations for {len(coordinate_systems)} coordinate systems")

    # Create plots for each metadata column across all coordinate systems
    for coord_system in coordinate_systems:
        coords = coord_system['coords']
        coord_names = coord_system['names']
        plot_prefix = coord_system['prefix']
        system_title = coord_system['title']
        point_size = coord_system['point_size']

        logger.info(f"   🎨 Creating plots for {system_title}...")

        for col in metadata_columns:
            if col not in adata.obs.columns:
                logger.warning(f"   ⚠️ Skipping {col} - not found in metadata")
                continue

            # Get unique values and colors
            unique_vals = sorted(adata.obs[col].unique())
            n_unique = len(unique_vals)

            if n_unique > 20:
                logger.warning(f"   ⚠️ Skipping {col} - too many categories ({n_unique})")
                continue

            # Create color palette
            if n_unique <= 10:
                colors = sns.color_palette("tab10", n_unique)
            else:
                colors = sns.color_palette("husl", n_unique)

            # Create plot
            plt.figure(figsize=config.FIGURE_SIZE)

            # Plot each category
            for i, val in enumerate(unique_vals):
                mask = adata.obs[col] == val
                plt.scatter(
                    coords[mask, 0],
                    coords[mask, 1],
                    c=[colors[i]],
                    label=str(val),
                    alpha=0.6,
                    s=point_size
                )

            plt.xlabel(coord_names[0], fontsize=12)
            plt.ylabel(coord_names[1], fontsize=12)
            plt.title(f'{system_title} - Colored by {col}', fontsize=14, fontweight='bold')

            # Add legend
            if n_unique <= 10:
                plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
            else:
                plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, ncol=2)

            plt.grid(True, alpha=0.3)
            plt.tight_layout()

            # Save plot
            plot_file = Path(config.PLOTS_DIR) / f"{plot_prefix}_overlay_{col.lower()}.png"
            plt.savefig(plot_file, dpi=config.DPI, bbox_inches='tight')
            plt.close()

            logger.info(f"   💾 Saved {col} overlay ({system_title}): {plot_file}")

    # Create comprehensive comparison plots
    if len(metadata_columns) >= 3 and len(coordinate_systems) > 1:
        # Multi-system comparison for each metadata column
        for col in metadata_columns[:3]:  # Limit to first 3 columns
            if col not in adata.obs.columns:
                continue

            # Create subplot grid (systems × metadata)
            n_systems = len(coordinate_systems)
            fig, axes = plt.subplots(1, n_systems, figsize=(6 * n_systems, 6))
            if n_systems == 1:
                axes = [axes]

            unique_vals = sorted(adata.obs[col].unique())
            n_unique = len(unique_vals)

            if n_unique > 20:
                continue

            # Use consistent colors across all systems
            if n_unique <= 10:
                colors = sns.color_palette("tab10", n_unique)
            else:
                colors = sns.color_palette("husl", n_unique)

            for sys_idx, coord_system in enumerate(coordinate_systems):
                coords = coord_system['coords']
                coord_names = coord_system['names']
                system_title = coord_system['title']
                point_size = coord_system['point_size']

                ax = axes[sys_idx]

                for i, val in enumerate(unique_vals):
                    mask = adata.obs[col] == val
                    ax.scatter(
                        coords[mask, 0],
                        coords[mask, 1],
                        c=[colors[i]],
                        label=str(val),
                        alpha=0.6,
                        s=point_size
                    )

                ax.set_xlabel(coord_names[0], fontsize=10)
                ax.set_ylabel(coord_names[1], fontsize=10)
                ax.set_title(f'{system_title} - {col}', fontsize=12, fontweight='bold')
                ax.grid(True, alpha=0.3)

                if n_unique <= 6:
                    ax.legend(fontsize=8)

            plt.suptitle(f'Coordinate System Comparison - {col}', fontsize=16, fontweight='bold')
            plt.tight_layout()

            # Save comparison plot
            comparison_file = Path(config.PLOTS_DIR) / f"system_comparison_{col.lower()}.png"
            plt.savefig(comparison_file, dpi=config.DPI, bbox_inches='tight')
            plt.close()

            logger.info(f"   💾 Saved system comparison for {col}: {comparison_file}")

    # Create individual combined plots for each coordinate system
    for coord_system in coordinate_systems:
        if len(metadata_columns) >= 3:
            coords = coord_system['coords']
            coord_names = coord_system['names']
            plot_prefix = coord_system['prefix']
            system_title = coord_system['title']
            point_size = coord_system['point_size']

            fig, axes = plt.subplots(1, 3, figsize=(18, 6))

            for idx, col in enumerate(metadata_columns[:3]):
                if col not in adata.obs.columns:
                    continue

                ax = axes[idx]
                unique_vals = sorted(adata.obs[col].unique())
                n_unique = len(unique_vals)

                if n_unique <= 10:
                    colors = sns.color_palette("tab10", n_unique)
                else:
                    colors = sns.color_palette("husl", n_unique)

                for i, val in enumerate(unique_vals):
                    mask = adata.obs[col] == val
                    ax.scatter(
                        coords[mask, 0],
                        coords[mask, 1],
                        c=[colors[i]],
                        label=str(val),
                        alpha=0.6,
                        s=point_size
                    )

                ax.set_xlabel(coord_names[0], fontsize=10)
                ax.set_ylabel(coord_names[1], fontsize=10)
                ax.set_title(f'{col}', fontsize=12, fontweight='bold')
                ax.grid(True, alpha=0.3)

                if n_unique <= 6:
                    ax.legend(fontsize=8)

            plt.suptitle(f'{system_title} - All Covariates', fontsize=16, fontweight='bold')
            plt.tight_layout()

            # Save combined plot
            combined_file = Path(config.PLOTS_DIR) / f"{plot_prefix}_combined_overlays.png"
            plt.savefig(combined_file, dpi=config.DPI, bbox_inches='tight')
            plt.close()

            logger.info(f"   💾 Saved combined overlays ({system_title}): {combined_file}")

    # Return information about coordinate systems used
    coordinate_info = {
        'systems_used': [sys['title'] for sys in coordinate_systems],
        'total_systems': len(coordinate_systems),
        'has_scvi_umap': any('scVI UMAP' in sys['title'] for sys in coordinate_systems),
        'has_lemur_umap': any('LEMUR UMAP' in sys['title'] for sys in coordinate_systems),
        'has_lemur_components': any('LEMUR Components' in sys['title'] for sys in coordinate_systems)
    }

    return coordinate_info


def generate_summary_report(variance_results, correlation_results, lemur_result, viz_results, config):
    """Generate comprehensive summary report"""

    logger.info("📋 Generating summary report...")

    report_file = Path(config.OUTPUT_DIR) / "latent_space_analysis_report.txt"

    with open(report_file, 'w') as f:
        f.write("LEMUR LATENT SPACE ANALYSIS REPORT\n")
        f.write("=" * 50 + "\n\n")

        f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("Dataset: GSE225158 - OUD vs Control\n\n")

        # Embedding summary
        f.write("EMBEDDING SUMMARY\n")
        f.write("-" * 20 + "\n")
        f.write(f"Number of cells: {lemur_result['n_cells']:,}\n")
        f.write(f"Number of components: {lemur_result['n_components']}\n")
        f.write(f"Embedding key: {lemur_result['embedding_key']}\n\n")

        # Variance analysis
        f.write("VARIANCE EXPLAINED ANALYSIS\n")
        f.write("-" * 30 + "\n")
        f.write(f"Top 5 components explain: {variance_results['top_5_cumulative']:.3f} of total variance\n")
        f.write(f"Component 1 variance: {variance_results['variance_explained'][0]:.3f}\n")
        f.write(f"Component 2 variance: {variance_results['variance_explained'][1]:.3f}\n")
        f.write(f"Component 3 variance: {variance_results['variance_explained'][2]:.3f}\n\n")

        # Correlation analysis
        f.write("COVARIATE CORRELATION ANALYSIS\n")
        f.write("-" * 35 + "\n")
        f.write(f"Total correlations computed: {correlation_results['correlation_matrix'].size}\n")
        f.write(f"Significant correlations: {correlation_results['n_significant']}\n")
        f.write(f"Significance threshold: p < {config.SIGNIFICANCE_THRESHOLD}\n\n")

        # Visualization systems
        if viz_results:
            f.write("COORDINATE SYSTEMS ANALYZED\n")
            f.write("-" * 30 + "\n")
            f.write(f"Total systems: {viz_results['total_systems']}\n")
            f.write(f"Systems used: {', '.join(viz_results['systems_used'])}\n")
            if viz_results['has_scvi_umap']:
                f.write("✅ scVI UMAP available (from preprocessing)\n")
            if viz_results['has_lemur_umap']:
                f.write("✅ LEMUR UMAP computed (from LEMUR embeddings)\n")
            if viz_results['has_lemur_components']:
                f.write("✅ LEMUR components used (direct embedding space)\n")
            f.write("\n")

        # Top correlations
        f.write("TOP ABSOLUTE CORRELATIONS\n")
        f.write("-" * 28 + "\n")

        # Flatten and sort correlations
        corr_flat = correlation_results['correlation_matrix'].stack()
        pval_flat = correlation_results['pvalue_matrix'].stack()

        # Get significant correlations
        sig_mask = pval_flat < config.SIGNIFICANCE_THRESHOLD
        sig_corr = corr_flat[sig_mask]

        if len(sig_corr) > 0:
            top_corr = sig_corr.abs().sort_values(ascending=False).head(10)

            for (component, covariate), corr_val in top_corr.items():
                actual_corr = correlation_results['correlation_matrix'].loc[component, covariate]
                pval = correlation_results['pvalue_matrix'].loc[component, covariate]
                f.write(f"{component} - {covariate}: r={actual_corr:.3f} (p={pval:.3e})\n")
        else:
            f.write("No significant correlations found.\n")

        f.write("\n" + "=" * 50 + "\n")
        f.write("Analysis completed successfully.\n")
        f.write(f"Plots saved to: {config.PLOTS_DIR}\n")
        if viz_results and viz_results['total_systems'] > 1:
            f.write("\nVisualization files include:\n")
            f.write("- Individual overlay plots for each coordinate system\n")
            f.write("- System comparison plots showing differences\n")
            f.write("- Combined overlay plots for comprehensive view\n")

    logger.info(f"   💾 Summary report saved: {report_file}")


# ================================================================================
# 🏗️ MAIN ANALYSIS CLASS
# ================================================================================

class LatentSpaceAnalyzer:
    """Main class for LEMUR latent space analysis"""

    def __init__(self, config):
        self.config = config

    def analyze_latent_structure(self, lemur_result, adata):
        """
        Perform comprehensive latent space diagnostics for LEMUR embeddings

        Parameters:
        -----------
        lemur_result : dict
            Dictionary containing LEMUR results with 'embedding' key
        adata : AnnData
            AnnData object with metadata columns
        """

        logger.info("\n🌊 LEMUR LATENT SPACE ANALYSIS")
        logger.info("=" * 35)

        try:
            # 1. Prepare metadata
            metadata_df, available_columns = prepare_metadata_for_analysis(
                adata, self.config.REQUIRED_METADATA, self.config.CELLTYPE_CANDIDATES
            )

            # 2. One-hot encode metadata
            encoded_df, encoded_features = one_hot_encode_metadata(
                metadata_df, available_columns
            )

            # 3. Analyze variance explained
            logger.info("\n📊 Step 1: Variance Explained Analysis")
            variance_results = analyze_variance_explained(lemur_result, self.config)

            # 4. Compute embedding-covariate correlations
            logger.info("\n🔗 Step 2: Embedding-Covariate Correlations")
            correlation_results = compute_embedding_covariate_correlations(
                lemur_result, encoded_df, encoded_features, self.config
            )

            # 5. Create latent space visualizations
            logger.info("\n🎨 Step 3: Latent Space Visualizations")
            viz_results = visualize_latent_space_overlays(
                lemur_result, adata, available_columns, self.config
            )

            # 6. Generate summary report
            logger.info("\n📋 Step 4: Summary Report Generation")
            generate_summary_report(
                variance_results, correlation_results, lemur_result, viz_results, self.config
            )

            logger.info("\n✅ LATENT SPACE ANALYSIS COMPLETE!")
            logger.info(f"📁 Results saved to: {self.config.OUTPUT_DIR}")

            return {
                'variance_results': variance_results,
                'correlation_results': correlation_results,
                'visualization_results': viz_results,
                'metadata_columns': available_columns,
                'success': True
            }

        except Exception as e:
            logger.error(f"❌ Latent space analysis failed: {e}")
            import traceback
            logger.error(f"Full traceback: {traceback.format_exc()}")
            return {'success': False, 'error': str(e)}


# ================================================================================
# 🔍 DIAGNOSTIC FUNCTIONS
# ================================================================================

def check_data_availability(config):
    """Check if required data and embeddings are available"""

    print("\n🔍 CHECKING DATA AVAILABILITY")
    print("=" * 35)

    # Check for h5ad files - prioritize LEMUR-processed files
    possible_paths = [
        # LEMUR-processed files (highest priority)
        f"{config.BASE_DIR}/data/processed/snrna_scvi/GSE225158_lemur_processed.h5ad",
        f"{config.LEMUR_RESULTS_DIR}/lemur_processed_data.h5ad",
        f"{config.BASE_DIR}/results/snrna_scvi/lemur_comprehensive/lemur_processed_data.h5ad",
        # Fallback files
        f"{config.BASE_DIR}/data/processed/snrna_scvi/GSE225158_annotated_scvi.h5ad",
        f"{config.BASE_DIR}/data/processed/snrna_scvi/GSE225158_reprocessed_scvi.h5ad",
        f"{config.LEMUR_RESULTS_DIR}/processed_data.h5ad",
        f"{config.LEMUR_RESULTS_DIR}/../lemur_analysis/processed_data.h5ad"
    ]

    found_data = False
    for path in possible_paths:
        if Path(path).exists():
            print(f"✅ Found data file: {path}")
            try:
                adata = sc.read_h5ad(path)
                print(f"   📊 Data: {adata.n_obs:,} cells, {adata.n_vars:,} genes")

                # Check embeddings
                embeddings = list(adata.obsm.keys())
                print(f"   🔧 Available embeddings: {embeddings}")

                # Check for LEMUR embeddings
                lemur_embeddings = [key for key in embeddings if 'lemur' in key.lower()]
                if lemur_embeddings:
                    print(f"   ✅ LEMUR embeddings found: {lemur_embeddings}")
                    print(f"   🎯 File type: {'LEMUR-processed' if 'lemur' in str(path).lower() else 'Fallback'}")
                else:
                    print("   ❌ No LEMUR embeddings found")
                    print("   💡 Run 03b_LEMUR.py first to generate LEMUR embeddings")
                    print("   📂 Expected LEMUR files:")
                    print(f"     • {config.BASE_DIR}/data/processed/snrna_scvi/GSE225158_lemur_processed.h5ad")
                    print(f"     • {config.LEMUR_RESULTS_DIR}/lemur_processed_data.h5ad")

                # Check metadata
                available_meta = [col for col in config.REQUIRED_METADATA if col in adata.obs.columns]
                celltype_meta = [col for col in config.CELLTYPE_CANDIDATES if col in adata.obs.columns]

                print(f"   📋 Required metadata available: {available_meta}")
                print(f"   📋 Cell type columns available: {celltype_meta}")

                if len(available_meta) >= 2 and celltype_meta:
                    print("   ✅ Sufficient metadata for analysis")
                else:
                    print("   ⚠️  Limited metadata available")

                found_data = True
                break

            except Exception as e:
                print(f"   ❌ Error reading file: {e}")

    if not found_data:
        print("❌ No suitable data files found")
        print("💡 Expected locations:")
        for path in possible_paths:
            print(f"   • {path}")

    return found_data


# ================================================================================
# 🚀 MAIN EXECUTION
# ================================================================================

def main(check_only=False):
    """Main execution function"""

    print("🌊 LEMUR LATENT SPACE ANALYSIS")
    print("=" * 35)
    print("📊 Comprehensive Embedding Diagnostics")
    print("🎯 Variance Analysis + Covariate Correlations")
    print("🎨 Latent Space Visualizations")
    print("\n🔧 IMPLEMENTED FIXES:")
    print("   ✅ Type consistency in categorical → numerical conversions")
    print("   ✅ Robust correlation computation with error handling")
    print("   ✅ Enhanced missing value handling in metadata processing")
    print("   ✅ Proper .astype('category') before .cat.codes operations")
    print("\n⚠️  PREREQUISITES:")
    print("   1. Run 03b_LEMUR.py first to generate LEMUR embeddings")
    print("   2. Ensure LEMUR-processed h5ad file exists with embeddings")
    print("   3. Required metadata: Sex, level1, and cell type annotations")
    print("\n📁 Expected LEMUR files:")
    print("   • GSE225158_lemur_processed.h5ad (main data location)")
    print("   • lemur_processed_data.h5ad (results location)")

    # Initialize configuration
    config = LatentSpaceConfig()

    # Run diagnostic check
    data_available = check_data_availability(config)

    if check_only:
        return data_available

    if not data_available:
        print("\n❌ Cannot proceed - required data not found")
        print("💡 Run with --check flag to diagnose data availability")
        return

    config.setup_directories()

    try:
        # Load LEMUR results
        logger.info("\n📂 Loading Data...")
        lemur_result, adata = load_lemur_results(config.LEMUR_RESULTS_DIR)

        if lemur_result is None or adata is None:
            logger.error("❌ Failed to load required data.")
            logger.error("💡 Make sure you have run 03b_LEMUR.py first to generate LEMUR embeddings.")
            logger.error("💡 The script expects a LEMUR-processed h5ad file with embeddings.")
            logger.error("📁 Expected locations:")
            logger.error("   • /data/processed/snrna_scvi/GSE225158_lemur_processed.h5ad")
            logger.error("   • /results/snrna_scvi/lemur_comprehensive/lemur_processed_data.h5ad")
            return

        # Initialize analyzer
        analyzer = LatentSpaceAnalyzer(config)

        # Run analysis
        results = analyzer.analyze_latent_structure(lemur_result, adata)

        if results['success']:
            print("\n🎉 ANALYSIS COMPLETED SUCCESSFULLY!")
            print("=" * 40)
            print(f"📊 Analyzed {lemur_result['n_cells']:,} cells")
            print(f"🧬 {lemur_result['n_components']} LEMUR components")
            print(f"📁 Results: {config.OUTPUT_DIR}")
            print(f"🎨 Plots: {config.PLOTS_DIR}")

            # Summary of key findings
            if 'variance_results' in results:
                var_res = results['variance_results']
                print(f"   • Top 5 components: {var_res['top_5_cumulative']:.1%} variance")
                print(f"   • Component 1: {var_res['variance_explained'][0]:.1%} variance")

            if 'correlation_results' in results:
                corr_res = results['correlation_results']
                print(f"   • {corr_res['n_significant']} significant covariate correlations")

            if 'visualization_results' in results:
                viz_res = results['visualization_results']
                print(f"   • {viz_res['total_systems']} coordinate systems analyzed")
                print(f"   • Systems: {', '.join(viz_res['systems_used'])}")

            print(f"\n✨ Review plots in: {config.PLOTS_DIR}")

        else:
            print(f"\n❌ Analysis failed: {results['error']}")

    except Exception as e:
        logger.error(f"❌ Main execution failed: {e}")
        raise


if __name__ == "__main__":
    import sys

    # Check for diagnostic flag
    if len(sys.argv) > 1 and sys.argv[1] in ['--check', '-c', '--diagnostic']:
        print("🔍 RUNNING DIAGNOSTIC CHECK ONLY")
        main(check_only=True)
    else:
        main()
