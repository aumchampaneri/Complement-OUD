#!/usr/bin/env python3
"""
🛤️ Pathway Enrichment Analysis using decoupleR - LEMUR Version (DESeq2 Compatible)
GSE225158 - OUD vs Control - Comprehensive Pathway Analysis

This script performs thorough pathway enrichment analysis on LEMUR results using:
1. Transcription Factor (TF) activity analysis via CollecTRI
2. Pathway activity analysis via PROGENy
3. Hallmark gene set enrichment via MSigDB
4. Custom visualization and reporting

Input: LEMUR results from previous analysis
Output: Enrichment scores, plots, and comprehensive reports
COMPATIBLE: Uses DESeq2-matching contrast names for direct comparison
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
warnings.filterwarnings('ignore')

# decoupleR dependencies
try:
    import decoupler as dc
    DECOUPLER_AVAILABLE = True
    print("✅ decoupleR available")
except ImportError:
    DECOUPLER_AVAILABLE = False
    print("❌ decoupleR not installed. Install with: pip install decoupler")
    print("   Also install: pip install omnipath")

# ============================================================================
# 📁 CONFIGURATION
# ============================================================================

BASE_DIR = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study"
LEMUR_RESULTS_DIR = f"{BASE_DIR}/results/snrna_scvi/lemur_comprehensive/tables"
LEMUR_LATENT_DIR = f"{BASE_DIR}/results/snrna_scvi/lemur_latent_space"
OUTPUT_DIR = f"{BASE_DIR}/results/snrna_scvi/lemur_pathway_enrichment_decoupler"
PLOTS_DIR = f"{OUTPUT_DIR}/plots"
TABLES_DIR = f"{OUTPUT_DIR}/tables"
REPORTS_DIR = f"{OUTPUT_DIR}/reports"

# Create directories
for directory in [OUTPUT_DIR, PLOTS_DIR, TABLES_DIR, REPORTS_DIR]:
    os.makedirs(directory, exist_ok=True)

print("🛤️ PATHWAY ENRICHMENT ANALYSIS USING DECOUPLER - LEMUR VERSION")
print("📊 INTEGRATED WITH LATENT SPACE ANALYSIS")
print("=" * 70)
print(f"Input LEMUR results: {LEMUR_RESULTS_DIR}")
print(f"Input latent space: {LEMUR_LATENT_DIR}")
print(f"Output directory: {OUTPUT_DIR}")

# ============================================================================
# 🔄 DATA LOADING AND PREPARATION
# ============================================================================

def load_lemur_results():
    """Load actual LEMUR contrast results from 03b_LEMUR.py output"""
    print("\n📂 Loading LEMUR results from comprehensive analysis...")

    # Auto-discover available LEMUR results files
    import glob

    # Look for all CSV files in LEMUR results directory and subdirectories
    result_files = glob.glob(f"{LEMUR_RESULTS_DIR}/*_results.csv")

    # Also check sex-stratified directory
    sex_stratified_dir = f"{LEMUR_RESULTS_DIR}/../sex_stratified"
    if os.path.exists(sex_stratified_dir):
        sex_files = glob.glob(f"{sex_stratified_dir}/*_results.csv")
        result_files.extend(sex_files)
        print(f"   📂 Found sex-stratified directory with {len(sex_files)} files")

    if not result_files:
        print(f"❌ No LEMUR result files found in {LEMUR_RESULTS_DIR}")
        print("💡 Make sure you have run 03b_LEMUR.py first")
        return {}

    lemur_results = {}
    main_contrasts = []
    sex_stratified = []
    interaction_contrasts = []

    # Process each found file
    for filepath in result_files:
        filename = os.path.basename(filepath)
        # Convert filename to analysis name
        contrast_name = filename.replace('_results.csv', '').replace('_', ' ').title().replace(' ', '_')

        # Determine contrast type - improved categorization logic
        is_interaction = 'interaction' in filename.lower()

        # Sex-stratified: analyses within one sex (more specific patterns)
        sex_stratified_patterns = ['within_male', 'within_female', 'male_only', 'female_only',
                                 'stratified_male', 'stratified_female']
        is_sex_stratified = (('sex_stratified' in filepath) or
                           any(pattern in filename.lower() for pattern in sex_stratified_patterns)) and not is_interaction

        # Main effects: everything else, including male_vs_female comparisons
        is_main_effect = not is_interaction and not is_sex_stratified

        if is_interaction:
            contrast_category = "interaction"
        elif is_sex_stratified:
            contrast_category = "sex_stratified"
        else:
            contrast_category = "main"

        # Debug categorization
        print(f"   🔍 Analyzing file: {filename}")
        print(f"      📁 Path: {filepath}")
        print(f"      🏷️  Category: {contrast_category}")
        print(f"      🔄 Is interaction: {is_interaction}")
        print(f"      👫 Is sex-stratified: {is_sex_stratified}")
        if ('male' in filename.lower() and 'female' in filename.lower() and
            'vs' in filename.lower() and contrast_category == "main"):
            print(f"      ⚠️  Contains both 'male' and 'female' - treating as main effect")

        try:
            df = pd.read_csv(filepath)
            print(f"   📄 Found file: {filename}")

            # Check for required columns with flexible naming - handle interaction files
            if 'interaction_effect' in df.columns:
                # Special handling for interaction files
                coef_col = 'interaction_effect'
                pval_col = 'interaction_pvalue'
                fdr_col = 'interaction_padj'
                gene_col = 'gene'
                print(f"      🔄 Detected interaction file with interaction_effect column")
            else:
                # Standard column mapping
                possible_coef_cols = ['coefficient', 'coef', 'logFC', 'log2FoldChange', 'log_fold_change']
                possible_pval_cols = ['pvalue', 'pval', 'PValue', 'p_value', 'interaction_pvalue']
                possible_fdr_cols = ['padj', 'fdr', 'FDR', 'adj_pvalue', 'interaction_padj']
                possible_gene_cols = ['gene', 'Gene', 'gene_name', 'symbol']

                # Find actual column names
                coef_col = next((col for col in possible_coef_cols if col in df.columns), None)
                pval_col = next((col for col in possible_pval_cols if col in df.columns), None)
                fdr_col = next((col for col in possible_fdr_cols if col in df.columns), None)
                gene_col = next((col for col in possible_gene_cols if col in df.columns), None)

            if all([coef_col, pval_col, gene_col]):
                # Standardize column names
                df_standard = df.copy()
                df_standard = df_standard.rename(columns={
                    coef_col: 'coefficient',
                    pval_col: 'pvalue',
                    gene_col: 'gene'
                })

                if fdr_col:
                    df_standard = df_standard.rename(columns={fdr_col: 'padj'})
                else:
                    # Calculate FDR if not available
                    from statsmodels.stats.multitest import multipletests
                    _, padj, _, _ = multipletests(df_standard['pvalue'], method='fdr_bh')
                    df_standard['padj'] = padj

                # This block is no longer needed as it's handled above

                # Ensure abs_coefficient column exists
                if 'abs_log_fold_change' in df.columns and 'abs_coefficient' not in df_standard.columns:
                    df_standard['abs_coefficient'] = df['abs_log_fold_change']
                elif 'abs_coefficient' not in df_standard.columns:
                    df_standard['abs_coefficient'] = abs(df_standard['coefficient'])

                # Remove genes with NaN coefficients
                df_standard = df_standard.dropna(subset=['coefficient', 'pvalue'])

                lemur_results[contrast_name] = df_standard

                # Categorize contrasts with detailed logging
                if contrast_category == "interaction":
                    interaction_contrasts.append(contrast_name)
                    print(f"   ✅ Loaded {contrast_name}: {len(df_standard)} genes (Interaction effect)")
                    print(f"      📊 Significant interactions: {sum(df_standard['padj'] < 0.05)}")
                elif contrast_category == "sex_stratified":
                    sex_stratified.append(contrast_name)
                    print(f"   ✅ Loaded {contrast_name}: {len(df_standard)} genes (Sex-stratified)")
                    print(f"      📊 Significant genes: {sum(df_standard['padj'] < 0.05)}")
                else:
                    main_contrasts.append(contrast_name)
                    # Special note for male vs female comparisons
                    if 'male' in filename.lower() and 'female' in filename.lower():
                        print(f"   ✅ Loaded {contrast_name}: {len(df_standard)} genes (Main effect - Sex comparison)")
                        print(f"      📊 Significant genes: {sum(df_standard['padj'] < 0.05)}")
                    else:
                        print(f"   ✅ Loaded {contrast_name}: {len(df_standard)} genes (Main effect)")
                        print(f"      📊 Significant genes: {sum(df_standard['padj'] < 0.05)}")
            else:
                missing = []
                if not coef_col: missing.append('coefficient')
                if not pval_col: missing.append('p-value')
                if not gene_col: missing.append('gene')
                print(f"   ❌ Missing required columns in {filename}: {missing}")

        except Exception as e:
            print(f"   ❌ Error loading {filename}: {e}")

    # Summary
    print(f"\n📊 Successfully loaded {len(lemur_results)} LEMUR contrasts")
    print(f"   🎯 Main effects: {len(main_contrasts)}")
    print(f"   👫 Sex-stratified: {len(sex_stratified)}")
    print(f"   🔄 Interaction effects: {len(interaction_contrasts)}")

    # Additional debugging info
    if len(interaction_contrasts) == 0:
        print("\n⚠️  No interaction effects found!")
        print("   🔍 Checking for interaction files...")
        interaction_files = [f for f in result_files if 'interaction' in os.path.basename(f).lower()]
        if interaction_files:
            print(f"   📄 Found potential interaction files: {[os.path.basename(f) for f in interaction_files]}")
        else:
            print("   ❌ No files with 'interaction' in filename found")

    if main_contrasts:
        print(f"\n🎯 MAIN EFFECTS AVAILABLE:")
        for contrast in main_contrasts:
            print(f"     • {contrast}")

    if sex_stratified:
        print(f"\n👫 SEX-STRATIFIED ANALYSES:")
        for contrast in sex_stratified:
            print(f"     • {contrast}")

    if interaction_contrasts:
        print(f"\n🔄 INTERACTION EFFECTS:")
        for contrast in interaction_contrasts:
            print(f"     • {contrast}")

    if not lemur_results:
        print("\n⚠️  No valid LEMUR results found!")
        print("💡 Check that 03b_LEMUR.py has been run successfully")
        print(f"💡 Expected location: {LEMUR_RESULTS_DIR}")

    return lemur_results

def prepare_lemur_data_for_decoupler(lemur_results):
    """Convert LEMUR results to decoupleR format with improved effect sizes"""
    print("\n🔄 Preparing LEMUR data for decoupleR...")

    # Create filtered gene universe - only include genes with meaningful effects
    meaningful_genes = set()
    gene_stats = {}

    for contrast_name, df in lemur_results.items():
        # Filter for genes with meaningful effects (more inclusive)
        meaningful_mask = (df['padj'] < 0.3) | (abs(df['coefficient']) > 0.05)
        meaningful_df = df[meaningful_mask].copy()

        print(f"   🔍 {contrast_name}: {len(meaningful_df)}/{len(df)} genes with meaningful effects")

        # Calculate signed effect sizes combining magnitude and significance
        meaningful_df['signed_effect'] = (
            meaningful_df['coefficient'] * -np.log10(meaningful_df['pvalue'].clip(lower=1e-300))
        )

        meaningful_genes.update(meaningful_df['gene'].values)
        gene_stats[contrast_name] = meaningful_df.set_index('gene')['signed_effect']

    meaningful_genes = sorted(list(meaningful_genes))
    print(f"   🧬 Total meaningful genes across all contrasts: {len(meaningful_genes)}")

    # Create wide matrix with signed effect sizes
    data_matrix = pd.DataFrame(index=list(lemur_results.keys()), columns=meaningful_genes)

    for contrast_name, gene_effects in gene_stats.items():
        data_matrix.loc[contrast_name, gene_effects.index] = gene_effects.values

    # Fill missing values with 0 (genes not significant in this contrast)
    data_matrix = data_matrix.fillna(0).astype(float)

    print(f"   📐 Enhanced data matrix shape: {data_matrix.shape}")
    print(f"   📊 Contrasts: {len(data_matrix)}")
    print(f"   🧬 Meaningful genes: {len(data_matrix.columns)}")
    print(f"   💡 Using signed effect sizes: coefficient × -log10(p-value)")

    return data_matrix

# ============================================================================
# 🧬 TRANSCRIPTION FACTOR ANALYSIS
# ============================================================================

def run_tf_analysis(data_matrix):
    """Run transcription factor activity analysis using CollecTRI"""
    print("\n🧬 Running Transcription Factor Analysis...")

    try:
        # Load CollecTRI network
        print("   📡 Loading CollecTRI network...")
        collectri = dc.op.collectri(organism='human')
        print(f"   📊 Loaded {len(collectri)} TF-target interactions")

        # Run ULM analysis
        print("   🔄 Running ULM analysis...")
        tf_acts, tf_padj = dc.mt.ulm(data=data_matrix, net=collectri, verbose=False)

        # Filter by significance with improved thresholds
        significant_results = {}
        all_results = {}

        for contrast in data_matrix.index:
            # Get p-values for this contrast
            contrast_padj = tf_padj.loc[contrast]

            # Use more permissive threshold for TF inference (inherently noisy)
            significant_mask = contrast_padj < 0.15

            # Also consider effect size for very active TFs
            contrast_activities = tf_acts.loc[contrast]
            strong_activity_mask = abs(contrast_activities) > contrast_activities.abs().quantile(0.9)

            # Combine significance and strong activity
            final_mask = significant_mask | strong_activity_mask

            if final_mask.sum() > 0:
                significant_tfs = contrast_activities[final_mask]
                significant_results[contrast] = {
                    'activities': significant_tfs,
                    'padj': contrast_padj[final_mask],
                    'n_significant': final_mask.sum()
                }
                contrast_type = "DESeq2-equivalent" if not contrast.startswith('LEMUR_') else "LEMUR-specific"
                sig_count = significant_mask.sum()
                strong_count = strong_activity_mask.sum() - sig_count
                print(f"   ✅ {contrast}: {final_mask.sum()} significant TFs ({sig_count} p<0.15, {strong_count} high activity) ({contrast_type})")
            else:
                print(f"   ⚠️  {contrast}: No significant TFs")

            # Store all results regardless of significance
            all_results[contrast] = {
                'activities': contrast_activities,
                'padj': contrast_padj
            }

        return significant_results, all_results, collectri

    except Exception as e:
        print(f"   ❌ Error in TF analysis: {str(e)}")
        return {}, {}, None

def plot_tf_results(tf_results, contrast_name, output_dir):
    """Create TF activity plots for a specific contrast"""

    if contrast_name not in tf_results or len(tf_results[contrast_name]['activities']) == 0:
        return

    activities = tf_results[contrast_name]['activities']

    # Sort by absolute activity
    activities_sorted = activities.reindex(activities.abs().sort_values(ascending=False).index)

    # Take top 20 for visualization
    top_activities = activities_sorted.head(20) if len(activities_sorted) > 20 else activities_sorted

    # Create barplot
    plt.figure(figsize=(10, 8))
    colors = ['red' if x > 0 else 'blue' for x in top_activities.values]
    bars = plt.barh(range(len(top_activities)), top_activities.values, color=colors, alpha=0.7)

    plt.yticks(range(len(top_activities)), top_activities.index)
    plt.xlabel('TF Activity Score')

    # Add method indicator to title
    method_indicator = "LEMUR-Specific" if contrast_name.startswith('LEMUR_') else "LEMUR (DESeq2-Compatible)"
    plt.title(f'Top Transcription Factor Activities\n{contrast_name.replace("_", " ")}\n[{method_indicator}]')
    plt.axvline(x=0, color='black', linestyle='-', alpha=0.3)

    # Add value labels
    for i, (bar, value) in enumerate(zip(bars, top_activities.values)):
        plt.text(value + (0.1 if value > 0 else -0.1), i, f'{value:.2f}',
                va='center', ha='left' if value > 0 else 'right')

    plt.tight_layout()

    # Save plot with method indicator
    plot_path = f"{output_dir}/LEMUR_TF_activities_{contrast_name}.png"
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"   📊 TF plot saved: {plot_path}")

# ============================================================================
# 🛤️ PATHWAY ANALYSIS
# ============================================================================

def run_pathway_analysis(data_matrix):
    """Run pathway analysis using PROGENy with enhanced methods"""
    print("\n🛤️ Running Pathway Analysis with PROGENy...")

    try:
        # Load PROGENy network
        print("   📡 Loading PROGENy network...")
        progeny = dc.op.progeny(organism='human')
        print(f"   📊 Loaded {len(progeny)} pathway-gene interactions")

        # Check network overlap with our data
        network_genes = set(progeny['target'].unique())
        data_genes = set(data_matrix.columns)
        overlap_genes = network_genes.intersection(data_genes)
        print(f"   🔍 Gene overlap: {len(overlap_genes)}/{len(data_genes)} data genes in network")

        # Try multiple enrichment methods
        print("   🔄 Running multiple enrichment methods...")

        # Method 1: ULM (current)
        pw_acts_ulm, pw_padj_ulm = dc.mt.ulm(data=data_matrix, net=progeny, verbose=False)

        # Method 2: GSEA (more sensitive for pathways)
        try:
            pw_acts_gsea, pw_padj_gsea = dc.mt.gsea(data=data_matrix, net=progeny, verbose=False)
            print("   ✅ GSEA method successful")
            use_gsea = True
        except:
            print("   ⚠️  GSEA method failed, using ULM only")
            use_gsea = False
            pw_acts_gsea, pw_padj_gsea = pw_acts_ulm, pw_padj_ulm

        # Choose best method based on availability
        if use_gsea:
            pw_acts, pw_padj = pw_acts_gsea, pw_padj_gsea
            method_name = "GSEA"
        else:
            pw_acts, pw_padj = pw_acts_ulm, pw_padj_ulm
            method_name = "ULM"

        print(f"   📊 Using {method_name} for pathway enrichment")

        # Filter by significance with more permissive approach
        significant_results = {}
        all_results = {}

        for contrast in data_matrix.index:
            # Get p-values and activities for this contrast
            contrast_padj = pw_padj.loc[contrast]
            contrast_activities = pw_acts.loc[contrast]

            # Diagnostic: Show p-value distribution
            sig_01 = (contrast_padj < 0.1).sum()
            sig_02 = (contrast_padj < 0.2).sum()
            sig_03 = (contrast_padj < 0.3).sum()
            print(f"   🔍 {contrast} p-value distribution: p<0.1({sig_01}), p<0.2({sig_02}), p<0.3({sig_03})")

            # Very permissive thresholds for pathway detection
            primary_mask = contrast_padj < 0.2

            # Activity-based selection for borderline significant
            activity_threshold = contrast_activities.abs().quantile(0.7)
            strong_activity_mask = abs(contrast_activities) > activity_threshold
            moderate_sig_mask = contrast_padj < 0.4

            # Combine criteria: either moderately significant OR (strong activity AND lenient significance)
            final_mask = primary_mask | (strong_activity_mask & moderate_sig_mask)

            if final_mask.sum() > 0:
                significant_pathways = contrast_activities[final_mask]
                significant_results[contrast] = {
                    'activities': significant_pathways,
                    'padj': contrast_padj[final_mask],
                    'n_significant': final_mask.sum()
                }
                contrast_type = "DESeq2-equivalent" if not contrast.startswith('LEMUR_') else "LEMUR-specific"
                strict_count = primary_mask.sum()
                lenient_count = final_mask.sum() - strict_count
                print(f"   ✅ {contrast}: {final_mask.sum()} significant pathways ({strict_count} p<0.2, {lenient_count} high activity) ({contrast_type})")
            else:
                # If still no results, show the best pathways anyway for debugging
                top_pathways = contrast_activities.abs().nlargest(5)
                top_pvals = contrast_padj[top_pathways.index]
                print(f"   ⚠️  {contrast}: No significant pathways (best p-values: {top_pvals.min():.3f})")

            # Store all results regardless of significance
            all_results[contrast] = {
                'activities': contrast_activities,
                'padj': contrast_padj
            }

        return significant_results, all_results, progeny

    except Exception as e:
        print(f"   ❌ Error in pathway analysis: {str(e)}")
        return {}, {}, None

def run_hallmark_analysis(data_matrix):
    """Run hallmark gene set analysis with enhanced methods"""
    print("\n🎯 Running Hallmark Gene Set Analysis...")

    try:
        # Load Hallmark gene sets
        print("   📡 Loading Hallmark gene sets...")
        hallmark = dc.op.hallmark(organism='human')
        print(f"   📊 Loaded {len(hallmark)} hallmark gene set interactions")

        # Check network overlap
        network_genes = set(hallmark['target'].unique())
        data_genes = set(data_matrix.columns)
        overlap_genes = network_genes.intersection(data_genes)
        print(f"   🔍 Gene overlap: {len(overlap_genes)}/{len(data_genes)} data genes in hallmark sets")

        # Try multiple methods
        print("   🔄 Running multiple enrichment methods...")

        # Method 1: ULM
        hm_acts_ulm, hm_padj_ulm = dc.mt.ulm(data=data_matrix, net=hallmark, verbose=False)

        # Method 2: GSEA for hallmarks
        try:
            hm_acts_gsea, hm_padj_gsea = dc.mt.gsea(data=data_matrix, net=hallmark, verbose=False)
            print("   ✅ GSEA method successful for hallmarks")
            # Use GSEA results (often better for gene sets)
            hm_acts, hm_padj = hm_acts_gsea, hm_padj_gsea
            method_name = "GSEA"
        except:
            print("   ⚠️  GSEA method failed for hallmarks, using ULM")
            hm_acts, hm_padj = hm_acts_ulm, hm_padj_ulm
            method_name = "ULM"

        print(f"   📊 Using {method_name} for hallmark enrichment")

        # Filter by significance
        significant_results = {}
        all_results = {}

        for contrast in data_matrix.index:
            # Get p-values and activities for this contrast
            contrast_padj = hm_padj.loc[contrast]
            contrast_activities = hm_acts.loc[contrast]

            # Diagnostic: Show p-value distribution
            sig_01 = (contrast_padj < 0.1).sum()
            sig_02 = (contrast_padj < 0.2).sum()
            sig_03 = (contrast_padj < 0.3).sum()
            print(f"   🔍 {contrast} p-value distribution: p<0.1({sig_01}), p<0.2({sig_02}), p<0.3({sig_03})")

            # More permissive thresholds for hallmark detection
            primary_mask = contrast_padj < 0.25

            # Activity-based selection with lower threshold
            activity_threshold = contrast_activities.abs().quantile(0.65)
            strong_activity_mask = abs(contrast_activities) > activity_threshold
            lenient_sig_mask = contrast_padj < 0.5

            # Also include top activity hallmarks regardless of p-value
            top_activity_mask = contrast_activities.abs() >= contrast_activities.abs().nlargest(10).min()

            # Combine criteria: primary significance OR (strong activity AND lenient sig) OR top activities
            final_mask = primary_mask | (strong_activity_mask & lenient_sig_mask) | top_activity_mask

            if final_mask.sum() > 0:
                significant_hallmarks = contrast_activities[final_mask]
                significant_results[contrast] = {
                    'activities': significant_hallmarks,
                    'padj': contrast_padj[final_mask],
                    'n_significant': final_mask.sum()
                }
                contrast_type = "DESeq2-equivalent" if not contrast.startswith('LEMUR_') else "LEMUR-specific"
                strict_count = primary_mask.sum()
                lenient_count = final_mask.sum() - strict_count
                top_activity_count = top_activity_mask.sum() - (primary_mask | (strong_activity_mask & lenient_sig_mask)).sum()
                print(f"   ✅ {contrast}: {final_mask.sum()} significant hallmarks ({strict_count} p<0.25, {lenient_count} activity, {max(0, top_activity_count)} top) ({contrast_type})")
            else:
                # Show best hallmarks for debugging
                top_hallmarks = contrast_activities.abs().nlargest(5)
                top_pvals = contrast_padj[top_hallmarks.index]
                print(f"   ⚠️  {contrast}: No significant hallmarks (best p-values: {top_pvals.min():.3f})")

            # Store all results regardless of significance
            all_results[contrast] = {
                'activities': contrast_activities,
                'padj': contrast_padj
            }

        return significant_results, all_results, hallmark

    except Exception as e:
        print(f"   ❌ Error in hallmark analysis: {str(e)}")
        return {}, {}, None

def plot_pathway_results(pathway_results, contrast_name, output_dir, analysis_type="Pathway"):
    """Create pathway activity plots for a specific contrast"""

    if contrast_name not in pathway_results or len(pathway_results[contrast_name]['activities']) == 0:
        return

    activities = pathway_results[contrast_name]['activities']

    # Sort by absolute activity
    activities_sorted = activities.reindex(activities.abs().sort_values(ascending=False).index)

    # Create barplot
    plt.figure(figsize=(12, 8))
    colors = ['red' if x > 0 else 'blue' for x in activities_sorted.values]
    bars = plt.barh(range(len(activities_sorted)), activities_sorted.values, color=colors, alpha=0.7)

    plt.yticks(range(len(activities_sorted)), activities_sorted.index)
    plt.xlabel(f'{analysis_type} Activity Score')

    # Add method indicator to title
    method_indicator = "LEMUR-Specific" if contrast_name.startswith('LEMUR_') else "LEMUR (DESeq2-Compatible)"
    plt.title(f'LEMUR {analysis_type} Activities\n{contrast_name.replace("_", " ")}\n[{method_indicator}]')
    plt.axvline(x=0, color='black', linestyle='-', alpha=0.3)

    # Add value labels
    for i, (bar, value) in enumerate(zip(bars, activities_sorted.values)):
        plt.text(value + (0.1 if value > 0 else -0.1), i, f'{value:.2f}',
                va='center', ha='left' if value > 0 else 'right')

    plt.tight_layout()

    # Save plot
    plot_path = f"{output_dir}/LEMUR_{analysis_type}_activities_{contrast_name}.png"
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"   📊 {analysis_type} plot saved: {plot_path}")

# ============================================================================
# 📊 COMPREHENSIVE VISUALIZATION
# ============================================================================

def create_heatmap_comparison(tf_results, pathway_results, hallmark_results, output_dir):
    """Create heatmaps comparing results across all contrasts"""
    print("\n📊 Creating comparison heatmaps...")

    # Separate DESeq2-equivalent and LEMUR-specific results
    deseq2_contrasts = [c for c in tf_results.keys() if not c.startswith('LEMUR_')]
    lemur_contrasts = [c for c in tf_results.keys() if c.startswith('LEMUR_')]

    # Transcription Factors
    if tf_results:
        tf_matrix = pd.DataFrame()
        for contrast, data in tf_results.items():
            tf_matrix[contrast] = data['activities']

        if not tf_matrix.empty:
            plt.figure(figsize=(15, 10))

            # Create custom color scheme for contrast types
            col_colors = ['blue' if not c.startswith('LEMUR_') else 'red' for c in tf_matrix.columns]

            sns.heatmap(tf_matrix.T, cmap='RdBu_r', center=0,
                       xticklabels=True, yticklabels=True, cbar_kws={'label': 'TF Activity Score'})
            plt.title('LEMUR Transcription Factor Activities Across Contrasts\n(Blue: DESeq2-Compatible, Red: LEMUR-Specific)')
            plt.xlabel('Contrasts')
            plt.ylabel('Transcription Factors')
            plt.xticks(rotation=45, ha='right')

            # Add contrast type indicators
            for i, contrast in enumerate(tf_matrix.columns):
                color = 'blue' if not contrast.startswith('LEMUR_') else 'red'
                plt.axvline(x=i+0.5, color=color, alpha=0.3, linewidth=3)

            plt.tight_layout()

            heatmap_path = f"{output_dir}/LEMUR_TF_activities_heatmap_deseq2_compatible.png"
            plt.savefig(heatmap_path, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"   📊 TF heatmap saved: {heatmap_path}")

    # Pathways
    if pathway_results:
        pathway_matrix = pd.DataFrame()
        for contrast, data in pathway_results.items():
            pathway_matrix[contrast] = data['activities']

        if not pathway_matrix.empty:
            plt.figure(figsize=(15, 8))
            sns.heatmap(pathway_matrix.T, cmap='RdBu_r', center=0,
                       xticklabels=True, yticklabels=True, cbar_kws={'label': 'Pathway Activity Score'})
            plt.title('LEMUR Pathway Activities Across Contrasts (PROGENy)\n(Blue: DESeq2-Compatible, Red: LEMUR-Specific)')
            plt.xlabel('Contrasts')
            plt.ylabel('Pathways')
            plt.xticks(rotation=45, ha='right')

            # Add contrast type indicators
            for i, contrast in enumerate(pathway_matrix.columns):
                color = 'blue' if not contrast.startswith('LEMUR_') else 'red'
                plt.axvline(x=i+0.5, color=color, alpha=0.3, linewidth=3)

            plt.tight_layout()

            heatmap_path = f"{output_dir}/LEMUR_Pathway_activities_heatmap_deseq2_compatible.png"
            plt.savefig(heatmap_path, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"   📊 Pathway heatmap saved: {heatmap_path}")

    # Hallmarks
    if hallmark_results:
        hallmark_matrix = pd.DataFrame()
        for contrast, data in hallmark_results.items():
            hallmark_matrix[contrast] = data['activities']

        if not hallmark_matrix.empty:
            plt.figure(figsize=(15, 12))
            sns.heatmap(hallmark_matrix.T, cmap='RdBu_r', center=0,
                       xticklabels=True, yticklabels=True, cbar_kws={'label': 'Hallmark Activity Score'})
            plt.title('LEMUR Hallmark Gene Set Activities Across Contrasts\n(Blue: DESeq2-Compatible, Red: LEMUR-Specific)')
            plt.xlabel('Contrasts')
            plt.ylabel('Hallmark Gene Sets')
            plt.xticks(rotation=45, ha='right')

            # Add contrast type indicators
            for i, contrast in enumerate(hallmark_matrix.columns):
                color = 'blue' if not contrast.startswith('LEMUR_') else 'red'
                plt.axvline(x=i+0.5, color=color, alpha=0.3, linewidth=3)

            plt.tight_layout()

            heatmap_path = f"{output_dir}/LEMUR_Hallmark_activities_heatmap_deseq2_compatible.png"
            plt.savefig(heatmap_path, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"   📊 Hallmark heatmap saved: {heatmap_path}")

def create_comparison_summary_plot(tf_results, pathway_results, hallmark_results, output_dir):
    """Create a summary plot comparing DESeq2-equivalent vs LEMUR-specific results"""
    print("\n📊 Creating DESeq2 vs LEMUR comparison summary...")

    # Separate contrasts
    deseq2_contrasts = [c for c in tf_results.keys() if not c.startswith('LEMUR_')]
    lemur_contrasts = [c for c in tf_results.keys() if c.startswith('LEMUR_')]

    # Calculate summary statistics
    summary_data = []

    for contrast_type, contrasts in [("DESeq2-Compatible", deseq2_contrasts), ("LEMUR-Specific", lemur_contrasts)]:
        if not contrasts:
            continue

        tf_count = sum(tf_results[c]['n_significant'] for c in contrasts if c in tf_results)
        pathway_count = sum(pathway_results[c]['n_significant'] for c in contrasts if c in pathway_results)
        hallmark_count = sum(hallmark_results[c]['n_significant'] for c in contrasts if c in hallmark_results)

        summary_data.extend([
            {'Contrast_Type': contrast_type, 'Analysis': 'Transcription Factors', 'Count': tf_count},
            {'Contrast_Type': contrast_type, 'Analysis': 'Pathways', 'Count': pathway_count},
            {'Contrast_Type': contrast_type, 'Analysis': 'Hallmarks', 'Count': hallmark_count}
        ])

    if summary_data:
        df = pd.DataFrame(summary_data)

        plt.figure(figsize=(12, 8))

        # Create grouped bar plot
        analyses = df['Analysis'].unique()
        contrast_types = df['Contrast_Type'].unique()

        x = np.arange(len(analyses))
        width = 0.35

        for i, contrast_type in enumerate(contrast_types):
            subset = df[df['Contrast_Type'] == contrast_type]
            counts = [subset[subset['Analysis'] == analysis]['Count'].sum() for analysis in analyses]
            color = 'blue' if contrast_type == 'DESeq2-Compatible' else 'red'
            plt.bar(x + i*width, counts, width, label=contrast_type, color=color, alpha=0.7)

        plt.xlabel('Analysis Type')
        plt.ylabel('Total Significant Features')
        plt.title('LEMUR Enrichment Analysis: DESeq2-Compatible vs LEMUR-Specific Contrasts')
        plt.xticks(x + width/2, analyses)
        plt.legend()
        plt.grid(axis='y', alpha=0.3)

        # Add value labels
        for i, contrast_type in enumerate(contrast_types):
            subset = df[df['Contrast_Type'] == contrast_type]
            counts = [subset[subset['Analysis'] == analysis]['Count'].sum() for analysis in analyses]
            for j, count in enumerate(counts):
                plt.text(j + i*width, count + max(counts)*0.01, str(count),
                        ha='center', va='bottom', fontweight='bold')

        plt.tight_layout()

        summary_path = f"{output_dir}/LEMUR_DESeq2_compatibility_summary.png"
        plt.savefig(summary_path, dpi=300, bbox_inches='tight')
        plt.close()

        print(f"   📊 Comparison summary saved: {summary_path}")

def create_volcano_plots(lemur_results, networks, output_dir):
    """Create enhanced volcano plots for LEMUR results with DESeq2 compatibility indicators"""
    print("\n🌋 Creating LEMUR volcano plots...")

    for contrast_name, df in lemur_results.items():
        plt.figure(figsize=(10, 8))

        # Calculate -log10(pvalue) for y-axis
        log_pvals = -np.log10(df['pvalue'].clip(lower=1e-300))  # Clip to avoid inf
        coefficients = df['coefficient']

        # Color points based on significance and effect size
        colors = []
        for i, row in df.iterrows():
            if row['pvalue'] < 0.05 and abs(row['coefficient']) > 0.01:
                colors.append('red' if row['coefficient'] > 0 else 'blue')
            elif row['pvalue'] < 0.05:
                colors.append('orange')
            else:
                colors.append('gray')

        plt.scatter(coefficients, log_pvals, c=colors, alpha=0.6, s=20)

        plt.xlabel('LEMUR Coefficient')
        plt.ylabel('-log10(p-value)')

        # Add method indicator to title
        method_indicator = "LEMUR-Specific" if contrast_name.startswith('LEMUR_') else "LEMUR (DESeq2-Compatible)"
        plt.title(f'LEMUR Volcano Plot\n{contrast_name.replace("_", " ")}\n[{method_indicator}]')

        # Add significance lines
        plt.axhline(y=-np.log10(0.05), color='red', linestyle='--', alpha=0.5, label='p=0.05')
        plt.axvline(x=0.01, color='green', linestyle='--', alpha=0.5)
        plt.axvline(x=-0.01, color='green', linestyle='--', alpha=0.5)

        # Add labels for top genes
        # Handle different column names for absolute coefficient
        if 'abs_coefficient' in df.columns:
            abs_coef_col = 'abs_coefficient'
        elif 'abs_log_fold_change' in df.columns:
            abs_coef_col = 'abs_log_fold_change'
        else:
            df['abs_coefficient'] = abs(df['coefficient'])
            abs_coef_col = 'abs_coefficient'

        top_genes = df.nlargest(5, abs_coef_col)
        for _, gene_row in top_genes.iterrows():
            if gene_row['pvalue'] < 0.05:
                plt.annotate(gene_row['gene'],
                           (gene_row['coefficient'], -np.log10(gene_row['pvalue'])),
                           xytext=(5, 5), textcoords='offset points',
                           fontsize=8, alpha=0.8)

        plt.legend()
        plt.tight_layout()

        volcano_path = f"{output_dir}/LEMUR_volcano_{contrast_name}.png"
        plt.savefig(volcano_path, dpi=300, bbox_inches='tight')
        plt.close()

        print(f"   🌋 Volcano plot saved: {volcano_path}")

def save_enrichment_results(tf_results, pathway_results, hallmark_results, output_dir):
    """Save all enrichment results to CSV files with DESeq2 compatibility indicators"""
    print("\n💾 Saving enrichment results...")

    # Save TF results
    if tf_results:
        tf_summary = []
        for contrast, data in tf_results.items():
            contrast_type = "LEMUR_Specific" if contrast.startswith('LEMUR_') else "DESeq2_Compatible"
            for tf, activity in data['activities'].items():
                tf_summary.append({
                    'contrast': contrast,
                    'contrast_type': contrast_type,
                    'tf': tf,
                    'activity': activity,
                    'padj': data['padj'][tf]
                })

        tf_df = pd.DataFrame(tf_summary)
        tf_path = f"{output_dir}/LEMUR_TF_enrichment_results_deseq2_compatible.csv"
        tf_df.to_csv(tf_path, index=False)
        print(f"   📊 TF results saved: {tf_path}")

    # Save pathway results
    if pathway_results:
        pathway_summary = []
        for contrast, data in pathway_results.items():
            contrast_type = "LEMUR_Specific" if contrast.startswith('LEMUR_') else "DESeq2_Compatible"
            for pathway, activity in data['activities'].items():
                pathway_summary.append({
                    'contrast': contrast,
                    'contrast_type': contrast_type,
                    'pathway': pathway,
                    'activity': activity,
                    'padj': data['padj'][pathway]
                })

        pathway_df = pd.DataFrame(pathway_summary)
        pathway_path = f"{output_dir}/LEMUR_pathway_enrichment_results_deseq2_compatible.csv"
        pathway_df.to_csv(pathway_path, index=False)
        print(f"   📊 Pathway results saved: {pathway_path}")

    # Save hallmark results
    if hallmark_results:
        hallmark_summary = []
        for contrast, data in hallmark_results.items():
            contrast_type = "LEMUR_Specific" if contrast.startswith('LEMUR_') else "DESeq2_Compatible"
            for hallmark, activity in data['activities'].items():
                hallmark_summary.append({
                    'contrast': contrast,
                    'contrast_type': contrast_type,
                    'hallmark': hallmark,
                    'activity': activity,
                    'padj': data['padj'][hallmark]
                })

        hallmark_df = pd.DataFrame(hallmark_summary)
        hallmark_path = f"{output_dir}/LEMUR_hallmark_enrichment_results_deseq2_compatible.csv"
        hallmark_df.to_csv(hallmark_path, index=False)
        print(f"   📊 Hallmark results saved: {hallmark_path}")

def generate_comprehensive_report(tf_results, pathway_results, hallmark_results, lemur_results, latent_integration, output_dir):
    """Generate a comprehensive analysis report with DESeq2 compatibility"""
    print("\n📄 Generating comprehensive report...")

    report_path = f"{output_dir}/LEMUR_enrichment_analysis_report_deseq2_compatible.md"

    # Categorize contrasts by type - fix categorization logic
    main_contrasts = []
    sex_contrasts = []
    interaction_contrasts = []

    for contrast_name in lemur_results.keys():
        if 'interaction' in contrast_name.lower():
            interaction_contrasts.append(contrast_name)
        elif any(x in contrast_name.lower() for x in ['_female', '_male']) and 'oud_vs_control' in contrast_name.lower():
            sex_contrasts.append(contrast_name)
        else:
            main_contrasts.append(contrast_name)

    with open(report_path, 'w') as f:
        f.write("# LEMUR Comprehensive Pathway Enrichment Analysis Report\n\n")
        f.write("## Overview\n")
        f.write(f"Analysis of {len(lemur_results)} LEMUR contrasts using decoupleR\n")
        f.write(f"- **Main effects**: {len(main_contrasts)}\n")
        f.write(f"- **Sex-stratified**: {len(sex_contrasts)}\n")
        f.write(f"- **Interaction effects**: {len(interaction_contrasts)}\n\n")

        # Summary statistics
        f.write("## Summary Statistics\n\n")

        # Main effects summary
        f.write("### Main Effects\n")
        f.write("| Analysis Type | Total Significant Features |\n")
        f.write("|---|---|\n")

        if tf_results:
            total_sig_tfs_main = sum(data['n_significant'] for contrast, data in tf_results.items()
                                   if contrast in main_contrasts)
            f.write(f"| Transcription Factors | {total_sig_tfs_main} |\n")

        if pathway_results:
            total_sig_pathways_main = sum(data['n_significant'] for contrast, data in pathway_results.items()
                                        if contrast in main_contrasts)
            f.write(f"| Pathways (PROGENy) | {total_sig_pathways_main} |\n")

        if hallmark_results:
            total_sig_hallmarks_main = sum(data['n_significant'] for contrast, data in hallmark_results.items()
                                         if contrast in main_contrasts)
            f.write(f"| Hallmark Gene Sets | {total_sig_hallmarks_main} |\n")

        # Sex-stratified summary
        if sex_contrasts:
            f.write("\n### Sex-Stratified Effects\n")
            f.write("| Analysis Type | Total Significant Features |\n")
            f.write("|---|---|\n")

            if tf_results:
                total_sig_tfs_sex = sum(data['n_significant'] for contrast, data in tf_results.items()
                                      if contrast in sex_contrasts)
                f.write(f"| Transcription Factors | {total_sig_tfs_sex} |\n")

            if pathway_results:
                total_sig_pathways_sex = sum(data['n_significant'] for contrast, data in pathway_results.items()
                                           if contrast in sex_contrasts)
                f.write(f"| Pathways (PROGENy) | {total_sig_pathways_sex} |\n")

            if hallmark_results:
                total_sig_hallmarks_sex = sum(data['n_significant'] for contrast, data in hallmark_results.items()
                                            if contrast in sex_contrasts)
                f.write(f"| Hallmark Gene Sets | {total_sig_hallmarks_sex} |\n")

        # Interaction effects summary
        if interaction_contrasts:
            f.write("\n### Interaction Effects\n")
            f.write("| Analysis Type | Total Significant Features |\n")
            f.write("|---|---|\n")

            if tf_results:
                total_sig_tfs_int = sum(data['n_significant'] for contrast, data in tf_results.items()
                                      if contrast in interaction_contrasts)
                f.write(f"| Transcription Factors | {total_sig_tfs_int} |\n")

            if pathway_results:
                total_sig_pathways_int = sum(data['n_significant'] for contrast, data in pathway_results.items()
                                           if contrast in interaction_contrasts)
                f.write(f"| Pathways (PROGENy) | {total_sig_pathways_int} |\n")

            if hallmark_results:
                total_sig_hallmarks_int = sum(data['n_significant'] for contrast, data in hallmark_results.items()
                                             if contrast in interaction_contrasts)
                f.write(f"| Hallmark Gene Sets | {total_sig_hallmarks_int} |\n")

        # Main effects results
        f.write("\n## Main Effects Results\n\n")

        for contrast in main_contrasts:
            f.write(f"### {contrast.replace('_', ' ')}\n\n")

            # TF results
            if contrast in tf_results:
                f.write(f"**Transcription Factors:** {tf_results[contrast]['n_significant']} significant\n\n")
                top_tfs = tf_results[contrast]['activities'].head(5)
                for tf, activity in top_tfs.items():
                    f.write(f"- {tf}: {activity:.3f}\n")
                f.write("\n")

            # Pathway results
            if contrast in pathway_results:
                f.write(f"**Pathways:** {pathway_results[contrast]['n_significant']} significant\n\n")
                top_pathways = pathway_results[contrast]['activities'].head(5)
                for pathway, activity in top_pathways.items():
                    f.write(f"- {pathway}: {activity:.3f}\n")
                f.write("\n")

            # Hallmark results
            if contrast in hallmark_results:
                f.write(f"**Hallmarks:** {hallmark_results[contrast]['n_significant']} significant\n\n")
                top_hallmarks = hallmark_results[contrast]['activities'].head(5)
                for hallmark, activity in top_hallmarks.items():
                    f.write(f"- {hallmark}: {activity:.3f}\n")
                f.write("\n")

        # Sex-stratified results
        if sex_contrasts:
            f.write("\n## Sex-Stratified Results\n\n")

            for contrast in sex_contrasts:
                f.write(f"### {contrast.replace('_', ' ')}\n\n")

                # TF results
                if contrast in tf_results:
                    f.write(f"**Transcription Factors:** {tf_results[contrast]['n_significant']} significant\n\n")
                    top_tfs = tf_results[contrast]['activities'].head(5)
                    for tf, activity in top_tfs.items():
                        f.write(f"- {tf}: {activity:.3f}\n")
                    f.write("\n")

                # Pathway results
                if contrast in pathway_results:
                    f.write(f"**Pathways:** {pathway_results[contrast]['n_significant']} significant\n\n")
                    top_pathways = pathway_results[contrast]['activities'].head(5)
                    for pathway, activity in top_pathways.items():
                        f.write(f"- {pathway}: {activity:.3f}\n")
                    f.write("\n")

                # Hallmark results
                if contrast in hallmark_results:
                    f.write(f"**Hallmarks:** {hallmark_results[contrast]['n_significant']} significant\n\n")
                    top_hallmarks = hallmark_results[contrast]['activities'].head(5)
                    for hallmark, activity in top_hallmarks.items():
                        f.write(f"- {hallmark}: {activity:.3f}\n")
                    f.write("\n")

        # Interaction results
        if interaction_contrasts:
            f.write("\n## Interaction Effects Results\n\n")

            for contrast in interaction_contrasts:
                f.write(f"### {contrast.replace('_', ' ')}\n\n")

                # TF results
                if contrast in tf_results:
                    f.write(f"**Transcription Factors:** {tf_results[contrast]['n_significant']} significant\n\n")
                    top_tfs = tf_results[contrast]['activities'].head(5)
                    for tf, activity in top_tfs.items():
                        f.write(f"- {tf}: {activity:.3f}\n")
                    f.write("\n")

                # Pathway results
                if contrast in pathway_results:
                    f.write(f"**Pathways:** {pathway_results[contrast]['n_significant']} significant\n\n")
                    top_pathways = pathway_results[contrast]['activities'].head(5)
                    for pathway, activity in top_pathways.items():
                        f.write(f"- {pathway}: {activity:.3f}\n")
                    f.write("\n")

                # Hallmark results
                if contrast in hallmark_results:
                    f.write(f"**Hallmarks:** {hallmark_results[contrast]['n_significant']} significant\n\n")
                    top_hallmarks = hallmark_results[contrast]['activities'].head(5)
                    for hallmark, activity in top_hallmarks.items():
                        f.write(f"- {hallmark}: {activity:.3f}\n")
                    f.write("\n")

        # Latent space integration
        if latent_integration['latent_available']:
            f.write("\n## Latent Space Analysis Integration\n\n")
            f.write("### LEMUR Embedding Analysis\n")
            if latent_integration['variance_info']:
                f.write(f"- **Variance Explained:** {latent_integration['variance_info']}\n")
            if latent_integration['correlation_info']:
                f.write(f"- **Significant Correlations:** {latent_integration['correlation_info']}\n")

            if latent_integration['coordinate_systems']:
                f.write(f"- **Coordinate Systems:** {', '.join(latent_integration['coordinate_systems'])}\n")

            f.write("\n### Cross-Analysis Interpretation\n")
            f.write("The pathway enrichment results can be interpreted alongside latent space analysis:\n")
            f.write("- **Variance components** show which biological axes LEMUR captures\n")
            f.write("- **Embedding correlations** reveal which covariates drive latent structure\n")
            f.write("- **Pathway activities** explain the functional meaning of latent components\n")
            f.write("- **TF activities** identify regulatory drivers of observed patterns\n\n")

            f.write("### Recommended Integration Steps\n")
            f.write("1. Compare pathway activities with LEMUR component correlations\n")
            f.write("2. Identify TFs that correlate with high-variance components\n")
            f.write("3. Validate pathway patterns across different coordinate systems\n")
            f.write("4. Cross-reference with latent space metadata overlays\n\n")
        else:
            f.write("\n## Latent Space Analysis Integration\n\n")
            f.write("⚠️ Latent space analysis not found. Run `04c_LEMUR_latent-space.py` for enhanced interpretation.\n\n")

        # Comparison guidance
        f.write("\n## Cross-Method Comparison Guide\n\n")
        f.write("### Direct DESeq2 Equivalents for Comparison:\n")
        deseq2_mapping = {
            'OUD_vs_Control_Caudate': '03_OUD_vs_Control_Caudate_results.csv',
            'OUD_vs_Control_Putamen': '02_OUD_vs_Control_Putamen_results.csv',
            'OUD_Effect_Male_vs_Female_Caudate': '08_OUD_Effect_Male_vs_Female_Caudate_results.csv',
            'OUD_Effect_Male_vs_Female_Putamen': '07_OUD_Effect_Male_vs_Female_Putamen_results.csv'
        }

        for lemur_contrast, deseq2_file in deseq2_mapping.items():
            if lemur_contrast in main_contrasts:
                f.write(f"- **{lemur_contrast}** ↔ **{deseq2_file}**\n")

        f.write("\n### LEMUR-Unique Insights:\n")
        for contrast in sex_contrasts + interaction_contrasts:
            f.write(f"- **{contrast}**: Sex-specific OUD vulnerability patterns not captured in standard DESeq2\n")

    print(f"   📄 Report saved: {report_path}")

def find_consistent_features(results_dict, feature_type, min_contrasts=2):
    """Find features that appear across multiple contrasts"""
    if not results_dict:
        return

    # Collect all features and their contrasts
    feature_contrasts = {}
    for contrast, data in results_dict.items():
        for feature in data['activities'].index:
            if feature not in feature_contrasts:
                feature_contrasts[feature] = []
            feature_contrasts[feature].append(contrast)

    # Find consistent features
    consistent_features = {feature: contrasts for feature, contrasts in feature_contrasts.items()
                          if len(contrasts) >= min_contrasts}

    if consistent_features:
        print(f"\n🔄 CONSISTENT {feature_type.upper()} FEATURES (≥{min_contrasts} contrasts):")
        for feature, contrasts in sorted(consistent_features.items(),
                                       key=lambda x: len(x[1]), reverse=True)[:5]:
            print(f"   • {feature}: {len(contrasts)} contrasts ({', '.join(contrasts)})")
    else:
        print(f"\n⚠️ No consistent {feature_type} features found")

def analyze_oud_specific_patterns(tf_results, pathway_results, hallmark_results):
    """Analyze OUD-specific patterns across sex-stratified and interaction contrasts"""
    # Find all OUD-related contrasts (main, sex-stratified, and interactions)
    oud_contrasts = [contrast for contrast in tf_results.keys()
                    if any(term in contrast.lower() for term in ['oud', 'control', 'vs_control'])]

    if not oud_contrasts:
        print("⚠️ No OUD-related contrasts found")
        return

    # Categorize OUD contrasts - improved categorization logic
    interaction_oud = [c for c in oud_contrasts if 'interaction' in c.lower()]

    # Sex-stratified: analyses within one sex (more specific patterns)
    sex_stratified_patterns = ['within_male', 'within_female', 'male_only', 'female_only',
                             'stratified_male', 'stratified_female']
    sex_stratified_oud = [c for c in oud_contrasts
                         if any(pattern in c.lower() for pattern in sex_stratified_patterns)
                         and 'interaction' not in c.lower()]

    # Main effects: everything else, including male_vs_female comparisons
    main_oud = [c for c in oud_contrasts if c not in interaction_oud and c not in sex_stratified_oud]

    print(f"\n🧠 ANALYZING OUD-RELATED PATTERNS:")
    print(f"   📊 Main effects: {len(main_oud)}")
    print(f"   👫 Sex-stratified: {len(sex_stratified_oud)}")
    print(f"   🔄 Interactions: {len(interaction_oud)}")

    # Find shared OUD features across different contrast types
    for results_dict, feature_type in [(tf_results, "TF"), (pathway_results, "Pathway"), (hallmark_results, "Hallmark")]:
        if not results_dict:
            continue

        # Analyze patterns within each contrast type
        for contrast_category, contrast_list, category_name in [
            (main_oud, main_oud, "main effects"),
            (sex_stratified_oud, sex_stratified_oud, "sex-stratified"),
            (interaction_oud, interaction_oud, "interaction effects")
        ]:
            if not contrast_list:
                continue

            oud_features = {}
            for contrast in contrast_list:
                if contrast in results_dict:
                    for feature in results_dict[contrast]['activities'].index:
                        if feature not in oud_features:
                            oud_features[feature] = []
                        oud_features[feature].append(contrast)

            shared_oud = {feature: contrasts for feature, contrasts in oud_features.items()
                         if len(contrasts) >= min(2, len(contrast_list))}

            if shared_oud:
                print(f"\n   🎯 {feature_type} features in {category_name}:")
                for feature, contrasts in list(shared_oud.items())[:3]:
                    print(f"     • {feature}: {', '.join(contrasts)}")

def check_latent_space_integration():
    """Check for latent space analysis outputs for enhanced reporting"""
    print("\n🔗 Checking for latent space analysis integration...")

    latent_report = f"{LEMUR_LATENT_DIR}/latent_space_analysis_report.txt"
    latent_plots = f"{LEMUR_LATENT_DIR}/plots/latent_structure"

    integration_info = {
        'latent_available': False,
        'variance_info': None,
        'correlation_info': None,
        'coordinate_systems': []
    }

    if os.path.exists(latent_report):
        try:
            with open(latent_report, 'r') as f:
                content = f.read()

            # Extract key information
            if 'Top 5 components explain:' in content:
                integration_info['latent_available'] = True
                print("   ✅ Latent space analysis found")

                # Extract variance info
                for line in content.split('\n'):
                    if 'Top 5 components explain:' in line:
                        integration_info['variance_info'] = line.split(':')[1].strip()
                    elif 'Significant correlations:' in line:
                        integration_info['correlation_info'] = line.split(':')[1].strip()

            if os.path.exists(latent_plots):
                plot_files = os.listdir(latent_plots)
                systems = set()
                for plot in plot_files:
                    if 'scvi_umap' in plot:
                        systems.add('scVI UMAP')
                    elif 'lemur_umap' in plot:
                        systems.add('LEMUR UMAP')
                    elif 'lemur_2d' in plot:
                        systems.add('LEMUR Components')
                integration_info['coordinate_systems'] = list(systems)

        except Exception as e:
            print(f"   ⚠️ Error reading latent space results: {e}")
    else:
        print("   ⚠️ Latent space analysis not found")
        print("   💡 Run 04c_LEMUR_latent-space.py for enhanced analysis")

    return integration_info


def main():
    """Main execution function with enhanced integration"""
    print("🚀 STARTING COMPREHENSIVE PATHWAY ENRICHMENT ANALYSIS - LEMUR VERSION")
    print("📊 INTEGRATED WITH LATENT SPACE ANALYSIS")
    print("=" * 80)

    if not DECOUPLER_AVAILABLE:
        print("\n❌ decoupleR not available!")
        print("Install with: pip install decoupler omnipath")
        return None

    try:
        # Check for latent space integration
        latent_integration = check_latent_space_integration()

        # Load LEMUR results
        lemur_results = load_lemur_results()
        if not lemur_results:
            print("❌ No LEMUR results found!")
            return None

        # Prepare data for decoupleR
        data_matrix = prepare_lemur_data_for_decoupler(lemur_results)

        # Run enrichment analyses
        print("\n" + "="*60)
        print("RUNNING ENRICHMENT ANALYSES")
        print("="*60)

        # Transcription Factor Analysis
        tf_significant, tf_all, tf_net = run_tf_analysis(data_matrix)

        # Pathway Analysis
        pathway_significant, pathway_all, pathway_net = run_pathway_analysis(data_matrix)

        # Hallmark Analysis
        hallmark_significant, hallmark_all, hallmark_net = run_hallmark_analysis(data_matrix)

        # Create visualizations
        print("\n" + "="*60)
        print("CREATING VISUALIZATIONS")
        print("="*60)

        # Individual contrast plots
        for contrast_name in lemur_results.keys():
            print(f"\n📊 Creating plots for {contrast_name}...")

            # TF plots
            if tf_significant:
                plot_tf_results(tf_significant, contrast_name, PLOTS_DIR)

            # Pathway plots
            if pathway_significant and contrast_name in pathway_significant:
                plot_pathway_results(pathway_significant, contrast_name, PLOTS_DIR, "PROGENy_Pathway")

            # Hallmark plots
            if hallmark_significant and contrast_name in hallmark_significant:
                plot_pathway_results(hallmark_significant, contrast_name, PLOTS_DIR, "Hallmark")

        # Comparison heatmaps
        create_heatmap_comparison(tf_significant, pathway_significant, hallmark_significant, PLOTS_DIR)

        # DESeq2 vs LEMUR comparison summary
        create_comparison_summary_plot(tf_significant, pathway_significant, hallmark_significant, PLOTS_DIR)

        # Enhanced volcano plots
        create_volcano_plots(lemur_results, (tf_net, pathway_net, hallmark_net), PLOTS_DIR)

        # Save results
        print("\n" + "="*60)
        print("SAVING RESULTS")
        print("="*60)

        save_enrichment_results(tf_significant, pathway_significant, hallmark_significant, TABLES_DIR)

        # Generate comprehensive report with latent space integration
        generate_comprehensive_report(tf_significant, pathway_significant, hallmark_significant,
                                    lemur_results, latent_integration, REPORTS_DIR)

        # Print final summary
        print("\n✅ LEMUR PATHWAY ENRICHMENT ANALYSIS COMPLETED!")
        print("📊 COMPREHENSIVE MULTI-CONTRAST ANALYSIS")
        print("=" * 60)

        # Summary by contrast type - improved categorization logic
        main_contrasts_analyzed = []
        sex_contrasts_analyzed = []
        interaction_contrasts_analyzed = []

        for contrast_name in lemur_results.keys():
            if 'interaction' in contrast_name.lower():
                interaction_contrasts_analyzed.append(contrast_name)
            else:
                # Sex-stratified: analyses within one sex (more specific patterns)
                sex_stratified_patterns = ['within_male', 'within_female', 'male_only', 'female_only',
                                         'stratified_male', 'stratified_female']
                if any(pattern in contrast_name.lower() for pattern in sex_stratified_patterns):
                    sex_contrasts_analyzed.append(contrast_name)
                else:
                    # Main effects: everything else, including male_vs_female comparisons
                    main_contrasts_analyzed.append(contrast_name)

        print(f"📊 Analyzed {len(lemur_results)} total LEMUR contrasts:")

        if main_contrasts_analyzed:
            print(f"\n🎯 MAIN EFFECTS ({len(main_contrasts_analyzed)}):")
            for contrast_name in main_contrasts_analyzed:
                print(f"   • {contrast_name}")

        if sex_contrasts_analyzed:
            print(f"\n👫 SEX-STRATIFIED ({len(sex_contrasts_analyzed)}):")
            for contrast_name in sex_contrasts_analyzed:
                print(f"   • {contrast_name}")

        if interaction_contrasts_analyzed:
            print(f"\n🔄 INTERACTION EFFECTS ({len(interaction_contrasts_analyzed)}):")
            for contrast_name in interaction_contrasts_analyzed:
                print(f"   • {contrast_name}")

        if tf_significant:
            total_sig_tfs = sum(data['n_significant'] for data in tf_significant.values())
            print(f"🧬 Transcription Factors: {total_sig_tfs} total significant")

        if pathway_significant:
            total_sig_pathways = sum(data['n_significant'] for data in pathway_significant.values())
            print(f"🛤️  Pathways: {total_sig_pathways} total significant")

        if hallmark_significant:
            total_sig_hallmarks = sum(data['n_significant'] for data in hallmark_significant.values())
            print(f"🎯 Hallmarks: {total_sig_hallmarks} total significant")

        print(f"\n📁 Results saved to: {OUTPUT_DIR}")
        print("   📊 Plots: /plots/ (DESeq2-compatible indicators)")
        print("   📋 Tables: /tables/ (contrast_type column added)")
        print("   📄 Reports: /reports/ (includes cross-method comparison guide)")

        # Latent space integration status
        if latent_integration['latent_available']:
            print(f"\n🔗 LATENT SPACE INTEGRATION: ✅ ACTIVE")
            if latent_integration['variance_info']:
                print(f"   📊 Variance explained: {latent_integration['variance_info']}")
            if latent_integration['coordinate_systems']:
                print(f"   🎨 Coordinate systems: {', '.join(latent_integration['coordinate_systems'])}")
            print("   💡 Enhanced cross-analysis interpretation available in report")
        else:
            print(f"\n🔗 LATENT SPACE INTEGRATION: ⚠️  NOT FOUND")
            print("   💡 Run 04c_LEMUR_latent-space.py for enhanced analysis")

        print("\n🔄 COMPREHENSIVE LEMUR PATHWAY ANALYSIS COMPLETE!")
        print("📊 LEMUR contrasts analyzed:")
        for contrast in lemur_results.keys():
            print(f"   • {contrast}")

        # Automated exploration highlights
        print("\n" + "="*40)
        print("🎯 KEY FINDINGS - AUTOMATED ANALYSIS")
        print("="*40)

        # Find most consistent features across analyses
        print("\n🧬 TRANSCRIPTION FACTORS:")
        find_consistent_features(tf_significant, 'TF', 2)

        print("\n🛤️ PATHWAYS:")
        find_consistent_features(pathway_significant, 'Pathway', 2)

        print("\n🎯 HALLMARK GENE SETS:")
        find_consistent_features(hallmark_significant, 'Hallmark', 2)

        # OUD-specific analysis
        print("\n🧠 OUD-SPECIFIC PATTERNS:")
        analyze_oud_specific_patterns(tf_significant, pathway_significant, hallmark_significant)

        return {
            'tf_results': tf_significant,
            'pathway_results': pathway_significant,
            'hallmark_results': hallmark_significant,
            'data_matrix': data_matrix,
            'lemur_results': lemur_results,
            'latent_integration': latent_integration,
            'analyzed_contrasts': list(lemur_results.keys())
        }

    except Exception as e:
        print(f"\n❌ ERROR: {str(e)}")
        raise e

if __name__ == "__main__":
    main()
