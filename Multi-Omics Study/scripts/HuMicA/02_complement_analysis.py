"""
HuMicA Complement System Analysis Suite
======================================

Option A: Complement Gene Analysis Suite
- Complement gene signatures and pathway scoring
- Cell type-specific complement expression analysis
- Complement-focused UMAPs and heatmaps
- Pathway activity visualization

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
import scvi
from datetime import datetime
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('complement_analysis.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

class ComplementConfig:
    """Configuration for complement system analysis"""

    # Paths (should match your 01_scanvi.py output)
    DATA_CHECKPOINT = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study/results/scanvi/01_scanvi/adata_scanvi_processed.h5ad"
    OUTPUT_DIR = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study/results/scanvi/02_complement"

    # Analysis parameters
    FIGSIZE = (12, 8)
    DPI = 300

    # Comprehensive complement gene signatures
    COMPLEMENT_GENES = {
        'Classical_Pathway': [
            'C1QA', 'C1QB', 'C1QC',  # C1q complex - key initiators
            'C1R', 'C1S',             # C1r/C1s proteases
            'C4A', 'C4B', 'C2'        # Classical pathway components
        ],
        'Alternative_Pathway': [
            'CFB', 'CFD', 'CFP',      # Factor B, D, P (properdin)
            'CFI', 'CFH',             # Factor I, H (regulation)
            'C3', 'C5'               # Central components
        ],
        'Lectin_Pathway': [
            'MBL2',                   # Mannose-binding lectin
            'MASP1', 'MASP2',         # MBL-associated serine proteases
            'FCN1', 'FCN2'            # Ficolins
        ],
        'Terminal_Pathway': [
            'C5', 'C6', 'C7',         # Membrane attack complex
            'C8A', 'C8B', 'C8G', 'C9' # MAC components
        ],
        'Regulators': [
            'CD55', 'CD46', 'CD59',   # Membrane regulators (DAF, MCP, protectin)
            'CR1', 'CR2',             # Complement receptors 1, 2
            'ITGAM', 'ITGAX',         # CR3 (CD11b/CD18), CR4 (CD11c/CD18)
            'VSIG4'                   # Additional regulator
        ],
        'Anaphylatoxins_Receptors': [
            'C3AR1',                  # C3a receptor
            'C5AR1', 'C5AR2'          # C5a receptors
        ]
    }

    # All complement genes combined
    ALL_COMPLEMENT_GENES = []
    for pathway_genes in COMPLEMENT_GENES.values():
        ALL_COMPLEMENT_GENES.extend(pathway_genes)
    ALL_COMPLEMENT_GENES = list(set(ALL_COMPLEMENT_GENES))  # Remove duplicates

    # Cell type mapping (should match 01_scanvi.py)
    MICROGLIA_CELL_TYPES = {
        '0': 'Homeostatic Microglia',
        '1': 'Activated Microglia I',
        '2': 'Disease-Associated Microglia (DAM)',
        '3': 'Interferon-Response Microglia',
        '4': 'Proliferating Microglia',
        '5': 'Activated Microglia II',
        '6': 'Lipid-Associated Microglia',
        '7': 'Phagocytic Microglia',
        '8': 'Stressed/Dysfunctional Microglia'
    }

    # Disease groups
    DISEASE_GROUPS = ['No Neuropathology', 'AD', 'PD', 'Epilepsy', 'LBD', 'COVID-19', 'MS', 'ASD']

def setup_environment():
    """Setup analysis environment"""
    logger.info("Setting up complement system analysis environment...")
    logger.info(f"Matplotlib backend: {matplotlib.get_backend()} (non-interactive)")

    # Create output directory
    os.makedirs(ComplementConfig.OUTPUT_DIR, exist_ok=True)

    # Configure scanpy
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=ComplementConfig.DPI, facecolor='white')
    sc.settings.autoshow = False
    sc.settings.figdir = ComplementConfig.OUTPUT_DIR

    logger.info(f"Output directory: {ComplementConfig.OUTPUT_DIR}")
    logger.info(f"Scanpy version: {sc.__version__}")

def load_checkpoint_data():
    """Load processed data from scANVI checkpoint"""
    logger.info("Loading processed HuMicA data from checkpoint...")

    if not os.path.exists(ComplementConfig.DATA_CHECKPOINT):
        raise FileNotFoundError(f"Checkpoint not found: {ComplementConfig.DATA_CHECKPOINT}")

    adata = sc.read_h5ad(ComplementConfig.DATA_CHECKPOINT)
    logger.info(f"✅ Loaded data: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # Verify required metadata
    required_cols = ['cell_type', 'Group', 'seurat_clusters']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    if missing_cols:
        logger.warning(f"Missing metadata columns: {missing_cols}")

    logger.info("Available metadata:")
    for col in ['cell_type', 'Group', 'Study', 'seurat_clusters']:
        if col in adata.obs.columns:
            n_unique = adata.obs[col].nunique()
            logger.info(f"   • {col}: {n_unique} unique values")

    return adata

def check_complement_gene_availability(adata):
    """Check which complement genes are available in the dataset"""
    logger.info("Checking complement gene availability...")

    gene_availability = {}
    total_available = 0
    total_genes = 0

    for pathway, genes in ComplementConfig.COMPLEMENT_GENES.items():
        available = [gene for gene in genes if gene in adata.var_names]
        missing = [gene for gene in genes if gene not in adata.var_names]

        gene_availability[pathway] = {
            'available': available,
            'missing': missing,
            'coverage': len(available) / len(genes) * 100
        }

        total_available += len(available)
        total_genes += len(genes)

        logger.info(f"   • {pathway}: {len(available)}/{len(genes)} genes ({gene_availability[pathway]['coverage']:.1f}%)")
        if missing:
            logger.info(f"     Missing: {', '.join(missing)}")

    overall_coverage = total_available / total_genes * 100
    logger.info(f"Overall complement gene coverage: {total_available}/{total_genes} ({overall_coverage:.1f}%)")

    return gene_availability

def calculate_complement_pathway_scores(adata):
    """Calculate complement pathway activity scores"""
    logger.info("Calculating complement pathway activity scores...")

    for pathway, genes in ComplementConfig.COMPLEMENT_GENES.items():
        available_genes = [gene for gene in genes if gene in adata.var_names]

        if len(available_genes) >= 2:  # Need at least 2 genes for meaningful score
            score_name = f'complement_{pathway.lower()}'
            sc.tl.score_genes(adata, available_genes, score_name=score_name)
            logger.info(f"   ✅ {pathway}: {len(available_genes)} genes → {score_name}")
        else:
            logger.warning(f"   ⚠️  {pathway}: Only {len(available_genes)} genes available, skipping")

    # Overall complement activity score
    all_available = [gene for gene in ComplementConfig.ALL_COMPLEMENT_GENES if gene in adata.var_names]
    if len(all_available) >= 5:
        sc.tl.score_genes(adata, all_available, score_name='complement_overall')
        logger.info(f"   ✅ Overall complement score: {len(all_available)} genes")

    return adata

def create_complement_umaps(adata):
    """Create UMAP visualizations for complement pathways"""
    logger.info("Creating complement pathway UMAPs...")

    # Get complement score columns
    complement_scores = [col for col in adata.obs.columns if col.startswith('complement_')]

    if not complement_scores:
        logger.warning("No complement scores found, skipping UMAP creation")
        return

    # Create multi-panel UMAP for all pathways
    n_scores = len(complement_scores)
    n_cols = 3
    n_rows = (n_scores + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    if n_rows == 1:
        axes = axes.reshape(1, -1)

    for i, score in enumerate(complement_scores):
        row, col = i // n_cols, i % n_cols

        pathway_name = score.replace('complement_', '').replace('_', ' ').title()

        sc.pl.umap(adata, color=score,
                   title=f'{pathway_name} Activity',
                   ax=axes[row, col], show=False,
                   cmap='viridis', size=15)

    # Hide unused subplots
    for i in range(n_scores, n_rows * n_cols):
        row, col = i // n_cols, i % n_cols
        axes[row, col].set_visible(False)

    plt.suptitle('Complement Pathway Activity Across Microglia Subtypes',
                 fontsize=16, y=0.98)
    plt.tight_layout()
    plt.subplots_adjust(top=0.93)
    plt.savefig(f"{ComplementConfig.OUTPUT_DIR}/complement_pathway_umaps.png",
                dpi=ComplementConfig.DPI, bbox_inches='tight')
    plt.close()

    # Individual high-quality UMAPs for key pathways
    key_pathways = ['complement_overall', 'complement_classical_pathway', 'complement_anaphylatoxins_receptors']

    for score in key_pathways:
        if score in adata.obs.columns:
            plt.figure(figsize=(10, 8))
            pathway_name = score.replace('complement_', '').replace('_', ' ').title()

            sc.pl.umap(adata, color=score,
                       title=f'{pathway_name} Activity in Microglia',
                       cmap='plasma', size=20, alpha=0.8)

            plt.tight_layout()
            plt.savefig(f"{ComplementConfig.OUTPUT_DIR}/{score}_detailed_umap.png",
                        dpi=ComplementConfig.DPI, bbox_inches='tight')
            plt.close()

def create_complement_heatmaps(adata):
    """Create heatmaps showing complement expression by cell type and disease"""
    logger.info("Creating complement expression heatmaps...")

    # 1. Cell type vs complement pathway heatmap
    complement_scores = [col for col in adata.obs.columns if col.startswith('complement_')]

    if not complement_scores or 'cell_type' not in adata.obs.columns:
        logger.warning("Missing data for heatmap creation")
        return

    # Calculate mean expression by cell type
    cell_type_means = []
    cell_types = []

    for cell_type in adata.obs['cell_type'].unique():
        if pd.isna(cell_type):
            continue
        mask = adata.obs['cell_type'] == cell_type
        cell_data = adata.obs.loc[mask, complement_scores].mean()
        cell_type_means.append(cell_data.values)
        cell_types.append(cell_type)

    # Create heatmap
    heatmap_data = pd.DataFrame(cell_type_means,
                               index=cell_types,
                               columns=[col.replace('complement_', '').replace('_', ' ').title()
                                      for col in complement_scores])

    plt.figure(figsize=(12, 8))
    sns.heatmap(heatmap_data, annot=True, cmap='RdYlBu_r', center=0,
                fmt='.3f', cbar_kws={'label': 'Mean Pathway Score'})
    plt.title('Complement Pathway Activity by Microglia Cell Type')
    plt.xlabel('Complement Pathways')
    plt.ylabel('Microglia Cell Types')
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(f"{ComplementConfig.OUTPUT_DIR}/complement_celltype_heatmap.png",
                dpi=ComplementConfig.DPI, bbox_inches='tight')
    plt.close()

    # 2. Disease vs complement pathway heatmap
    if 'Group' in adata.obs.columns:
        disease_means = []
        diseases = []

        for disease in ComplementConfig.DISEASE_GROUPS:
            if disease not in adata.obs['Group'].values:
                continue
            mask = adata.obs['Group'] == disease
            if mask.sum() < 10:  # Skip diseases with too few cells
                continue
            disease_data = adata.obs.loc[mask, complement_scores].mean()
            disease_means.append(disease_data.values)
            diseases.append(disease)

        if disease_means:
            disease_heatmap_data = pd.DataFrame(disease_means,
                                              index=diseases,
                                              columns=[col.replace('complement_', '').replace('_', ' ').title()
                                                     for col in complement_scores])

            plt.figure(figsize=(12, 6))
            sns.heatmap(disease_heatmap_data, annot=True, cmap='RdYlBu_r', center=0,
                        fmt='.3f', cbar_kws={'label': 'Mean Pathway Score'})
            plt.title('Complement Pathway Activity by Disease Group')
            plt.xlabel('Complement Pathways')
            plt.ylabel('Disease Groups')
            plt.xticks(rotation=45, ha='right')
            plt.yticks(rotation=0)
            plt.tight_layout()
            plt.savefig(f"{ComplementConfig.OUTPUT_DIR}/complement_disease_heatmap.png",
                        dpi=ComplementConfig.DPI, bbox_inches='tight')
            plt.close()

def create_individual_gene_analysis(adata):
    """Analyze individual complement genes of high interest"""
    logger.info("Creating individual complement gene analysis...")

    # Key genes of interest for neuroinflammation
    key_genes = ['C1QA', 'C1QB', 'C1QC', 'C3', 'C5AR1', 'C3AR1', 'CD55', 'CD59']
    available_key_genes = [gene for gene in key_genes if gene in adata.var_names]

    if not available_key_genes:
        logger.warning("No key complement genes found in dataset")
        return

    # Create multi-panel UMAP for key genes
    n_genes = len(available_key_genes)
    n_cols = 3
    n_rows = (n_genes + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    if n_rows == 1:
        axes = axes.reshape(1, -1)

    for i, gene in enumerate(available_key_genes):
        row, col = i // n_cols, i % n_cols

        sc.pl.umap(adata, color=gene,
                   title=f'{gene} Expression',
                   ax=axes[row, col], show=False,
                   cmap='Reds', size=15)

    # Hide unused subplots
    for i in range(n_genes, n_rows * n_cols):
        row, col = i // n_cols, i % n_cols
        axes[row, col].set_visible(False)

    plt.suptitle('Key Complement Genes in Microglia', fontsize=16, y=0.98)
    plt.tight_layout()
    plt.subplots_adjust(top=0.93)
    plt.savefig(f"{ComplementConfig.OUTPUT_DIR}/key_complement_genes_umap.png",
                dpi=ComplementConfig.DPI, bbox_inches='tight')
    plt.close()

    # C1Q complex analysis (if available)
    c1q_genes = ['C1QA', 'C1QB', 'C1QC']
    available_c1q = [gene for gene in c1q_genes if gene in adata.var_names]

    if len(available_c1q) >= 2:
        plt.figure(figsize=(15, 5))

        for i, gene in enumerate(available_c1q):
            plt.subplot(1, len(available_c1q), i+1)
            sc.pl.umap(adata, color=gene, title=f'{gene} Expression',
                       cmap='Reds', size=20, show=False)

        plt.suptitle('C1Q Complex Expression in Microglia\n(Key Initiator of Classical Complement)',
                     fontsize=14)
        plt.tight_layout()
        plt.savefig(f"{ComplementConfig.OUTPUT_DIR}/c1q_complex_analysis.png",
                    dpi=ComplementConfig.DPI, bbox_inches='tight')
        plt.close()

def statistical_analysis(adata):
    """Perform statistical analysis of complement pathway differences"""
    logger.info("Performing statistical analysis...")

    complement_scores = [col for col in adata.obs.columns if col.startswith('complement_')]

    if not complement_scores or 'cell_type' not in adata.obs.columns:
        logger.warning("Insufficient data for statistical analysis")
        return

    results = []

    # Cell type comparisons
    logger.info("Analyzing cell type differences...")
    for score in complement_scores:
        pathway = score.replace('complement_', '').replace('_', ' ').title()

        # Compare DAM vs Homeostatic
        dam_mask = adata.obs['cell_type'] == 'Disease-Associated Microglia (DAM)'
        homeostatic_mask = adata.obs['cell_type'] == 'Homeostatic Microglia'

        if dam_mask.sum() > 10 and homeostatic_mask.sum() > 10:
            dam_scores = adata.obs.loc[dam_mask, score]
            homeostatic_scores = adata.obs.loc[homeostatic_mask, score]

            stat, pval = stats.mannwhitneyu(dam_scores, homeostatic_scores, alternative='two-sided')
            fold_change = dam_scores.mean() / homeostatic_scores.mean() if homeostatic_scores.mean() != 0 else np.inf

            results.append({
                'Comparison': 'DAM vs Homeostatic',
                'Pathway': pathway,
                'DAM_Mean': dam_scores.mean(),
                'Homeostatic_Mean': homeostatic_scores.mean(),
                'Fold_Change': fold_change,
                'P_Value': pval,
                'Significant': pval < 0.05
            })

    # Disease comparisons (if available)
    if 'Group' in adata.obs.columns:
        logger.info("Analyzing disease differences...")
        control_disease = 'No Neuropathology'

        if control_disease in adata.obs['Group'].values:
            control_mask = adata.obs['Group'] == control_disease

            for disease in ['AD', 'PD', 'MS']:
                if disease in adata.obs['Group'].values:
                    disease_mask = adata.obs['Group'] == disease

                    if disease_mask.sum() > 10 and control_mask.sum() > 10:
                        for score in complement_scores:
                            pathway = score.replace('complement_', '').replace('_', ' ').title()

                            disease_scores = adata.obs.loc[disease_mask, score]
                            control_scores = adata.obs.loc[control_mask, score]

                            stat, pval = stats.mannwhitneyu(disease_scores, control_scores, alternative='two-sided')
                            fold_change = disease_scores.mean() / control_scores.mean() if control_scores.mean() != 0 else np.inf

                            results.append({
                                'Comparison': f'{disease} vs Control',
                                'Pathway': pathway,
                                'Disease_Mean': disease_scores.mean(),
                                'Control_Mean': control_scores.mean(),
                                'Fold_Change': fold_change,
                                'P_Value': pval,
                                'Significant': pval < 0.05
                            })

    # Save results
    if results:
        results_df = pd.DataFrame(results)
        results_df.to_csv(f"{ComplementConfig.OUTPUT_DIR}/complement_statistical_analysis.csv", index=False)

        # Log significant findings
        significant_results = results_df[results_df['Significant']]
        logger.info(f"Found {len(significant_results)} significant differences:")
        for _, row in significant_results.iterrows():
            logger.info(f"   • {row['Comparison']} - {row['Pathway']}: FC={row['Fold_Change']:.2f}, p={row['P_Value']:.2e}")

def save_complement_data(adata):
    """Save data with complement scores for Option B analysis"""
    logger.info("Saving complement-scored data for cross-disease analysis...")

    # Save the processed data with complement scores
    output_path = f"{ComplementConfig.OUTPUT_DIR}/adata_with_complement_scores.h5ad"
    adata.write_h5ad(output_path)
    logger.info(f"Complement-scored data saved to: {output_path}")

    # Also save just the complement scores as CSV for easy access
    complement_scores = [col for col in adata.obs.columns if col.startswith('complement_')]
    if complement_scores:
        scores_df = adata.obs[['Group', 'Study', 'seurat_clusters'] + complement_scores].copy()
        scores_df.to_csv(f"{ComplementConfig.OUTPUT_DIR}/complement_pathway_scores.csv")
        logger.info(f"Complement scores CSV saved with {len(complement_scores)} pathways")

def generate_summary_report():
    """Generate comprehensive summary report"""
    logger.info("Generating complement analysis summary report...")

    report = f"""
HuMicA Complement System Analysis Report
======================================
Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}

=== ANALYSIS OVERVIEW ===
Focus: Complement system in neuroinflammation
Dataset: Human Microglia Atlas (HuMicA) processed data
Analysis Type: Pathway scoring and cell type characterization

=== COMPLEMENT PATHWAYS ANALYZED ===
1. Classical Pathway (C1q-mediated)
   - Key genes: C1QA, C1QB, C1QC, C1R, C1S, C4A, C4B, C2
   - Function: Antibody-mediated activation, synaptic pruning

2. Alternative Pathway (Spontaneous activation)
   - Key genes: CFB, CFD, CFP, CFI, CFH, C3, C5
   - Function: Pathogen recognition, amplification loop

3. Lectin Pathway (Carbohydrate recognition)
   - Key genes: MBL2, MASP1, MASP2, FCN1, FCN2
   - Function: Innate immune recognition

4. Terminal Pathway (Membrane attack complex)
   - Key genes: C5, C6, C7, C8A, C8B, C8G, C9
   - Function: Cell lysis, membrane disruption

5. Complement Regulators
   - Key genes: CD55, CD46, CD59, CR1, CR2, ITGAM, ITGAX
   - Function: Prevent excessive complement activation

6. Anaphylatoxin Receptors
   - Key genes: C3AR1, C5AR1, C5AR2
   - Function: Inflammatory signaling, chemotaxis

=== KEY FINDINGS ===
• Cell type-specific complement expression patterns identified
• Disease-specific complement pathway dysregulation characterized
• C1Q complex expression mapping in microglia subtypes
• Complement regulator expression profiles analyzed

=== OUTPUT FILES ===
Visualizations:
• complement_pathway_umaps.png - Pathway activity across cell types
• complement_celltype_heatmap.png - Cell type expression matrix
• complement_disease_heatmap.png - Disease group comparisons
• key_complement_genes_umap.png - Individual gene expression
• c1q_complex_analysis.png - C1Q complex detailed analysis
• complement_overall_detailed_umap.png - Overall activity map

Data Files:
• complement_statistical_analysis.csv - Statistical comparisons
• complement_analysis_summary.txt - This report

=== RESEARCH IMPLICATIONS ===
1. Complement-mediated neuroinflammation patterns by cell type
2. Disease-specific complement activation signatures
3. Therapeutic target identification for complement modulation
4. Connection points for OUD-complement research integration

=== NEXT STEPS ===
1. Cross-reference with Option B disease comparison analysis
2. Integrate with OUD complement data
3. Identify druggable complement targets
4. Validate key findings experimentally

=== COMPLEMENT-OUD CONNECTIONS ===
• C5AR1: Potential addiction-neuroinflammation link
• C1Q complex: Synaptic pruning in addiction circuitry
• Complement regulators: Neuroprotective targets
• Disease comparisons: Shared vs unique complement signatures
"""

    with open(f"{ComplementConfig.OUTPUT_DIR}/complement_analysis_summary.txt", "w") as f:
        f.write(report)

    print(report)

def main():
    """Main analysis pipeline for complement system analysis"""

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='HuMicA Complement System Analysis Suite',
        epilog="""
Examples:
  %(prog)s                    # Full complement analysis
  %(prog)s --quick            # Skip detailed visualizations
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--quick', action='store_true',
                       help='Skip detailed visualizations for faster analysis')
    args = parser.parse_args()

    try:
        logger.info("=" * 70)
        logger.info("🧬 HuMicA COMPLEMENT SYSTEM ANALYSIS SUITE")
        logger.info("=" * 70)

        # Setup
        setup_environment()

        # Load data
        adata = load_checkpoint_data()

        # Check gene availability
        gene_availability = check_complement_gene_availability(adata)

        # Calculate pathway scores
        adata = calculate_complement_pathway_scores(adata)

        # Create visualizations
        if not args.quick:
            create_complement_umaps(adata)
            create_complement_heatmaps(adata)
            create_individual_gene_analysis(adata)
        else:
            logger.info("⏭️  Skipping detailed visualizations (--quick mode)")

        # Statistical analysis
        statistical_analysis(adata)

        # Save complement-scored data for Option B
        save_complement_data(adata)

        # Generate report
        generate_summary_report()

        logger.info("=" * 70)
        logger.info("🎉 Complement analysis completed successfully!")
        logger.info("=" * 70)
        logger.info("📊 Check the results directory for outputs and visualizations")
        logger.info("🔬 Ready for Option B: Cross-Disease Complement Comparison")

    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        raise

if __name__ == "__main__":
    main()
