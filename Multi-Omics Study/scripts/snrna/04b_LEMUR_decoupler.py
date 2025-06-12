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
LEMUR_RESULTS_DIR = f"{BASE_DIR}/results/snrna_scvi/lemur_analysis/tables"
OUTPUT_DIR = f"{BASE_DIR}/results/snrna_scvi/lemur_pathway_enrichment_decoupler"
PLOTS_DIR = f"{OUTPUT_DIR}/plots"
TABLES_DIR = f"{OUTPUT_DIR}/tables"
REPORTS_DIR = f"{OUTPUT_DIR}/reports"

# Create directories
for directory in [OUTPUT_DIR, PLOTS_DIR, TABLES_DIR, REPORTS_DIR]:
    os.makedirs(directory, exist_ok=True)

print("🛤️ PATHWAY ENRICHMENT ANALYSIS USING DECOUPLER - LEMUR VERSION")
print("📊 DESEQ2-COMPATIBLE CONTRAST NAMING")
print("=" * 70)
print(f"Input LEMUR results: {LEMUR_RESULTS_DIR}")
print(f"Output directory: {OUTPUT_DIR}")

# ============================================================================
# 🔄 DATA LOADING AND PREPARATION
# ============================================================================

def load_lemur_results():
    """Load all LEMUR contrast results with DESeq2-compatible naming"""
    print("\n📂 Loading LEMUR results with DESeq2-compatible contrast names...")
    
    # Define contrast mapping: DESeq2_name -> LEMUR_filename
    # Only include contrasts that have DIRECT EQUIVALENTS between methods
    contrast_mapping = {
        # DIRECT EQUIVALENTS (4 contrasts)
        'OUD_vs_Control_Caudate': 'strategic_caudate_oud_effect_results.csv',
        'OUD_vs_Control_Putamen': 'strategic_putamen_oud_effect_results.csv', 
        'OUD_Effect_Male_vs_Female_Caudate': 'strategic_caudate_sex_interaction_results.csv',
        'OUD_Effect_Male_vs_Female_Putamen': 'strategic_putamen_sex_interaction_results.csv',
        
        # LEMUR-SPECIFIC CONTRASTS (marked clearly)
        'LEMUR_Sex_OUD_Vulnerability_Caudate': 'strategic_caudate_sex_oud_results.csv',
        'LEMUR_Sex_OUD_Vulnerability_Putamen': 'strategic_putamen_sex_oud_results.csv'
    }
    
    lemur_results = {}
    direct_equivalents = []
    lemur_specific = []
    
    for contrast_name, filename in contrast_mapping.items():
        filepath = f"{LEMUR_RESULTS_DIR}/{filename}"
        
        if os.path.exists(filepath):
            df = pd.read_csv(filepath)
            
            # Ensure required columns exist for LEMUR data
            required_cols = ['gene', 'coefficient', 'pval', 'fdr', 'z_score']
            if all(col in df.columns for col in required_cols):
                # Remove genes with NaN coefficients
                df = df.dropna(subset=['coefficient', 'z_score'])
                lemur_results[contrast_name] = df
                
                # Categorize contrasts
                if contrast_name.startswith('LEMUR_'):
                    lemur_specific.append(contrast_name)
                    print(f"   ✅ Loaded {contrast_name}: {len(df)} genes (LEMUR-specific)")
                else:
                    direct_equivalents.append(contrast_name)
                    print(f"   ✅ Loaded {contrast_name}: {len(df)} genes (DESeq2 equivalent)")
            else:
                missing_cols = [col for col in required_cols if col not in df.columns]
                print(f"   ❌ Missing required columns in {contrast_name}: {missing_cols}")
        else:
            print(f"   ❌ File not found: {filepath}")
    
    # Summary
    print(f"\n📊 Successfully loaded {len(lemur_results)} LEMUR contrasts")
    print(f"   🔄 Direct DESeq2 equivalents: {len(direct_equivalents)}")
    print(f"   🧠 LEMUR-specific contrasts: {len(lemur_specific)}")
    
    if direct_equivalents:
        print(f"\n🎯 DESeq2-EQUIVALENT CONTRASTS:")
        for contrast in direct_equivalents:
            print(f"     • {contrast}")
    
    if lemur_specific:
        print(f"\n🧬 LEMUR-SPECIFIC CONTRASTS:")
        for contrast in lemur_specific:
            print(f"     • {contrast}")
    
    # Note missing DESeq2 contrasts
    missing_deseq2 = [
        'OUD_vs_Control_Pooled',
        'OUD_vs_Control_Male', 
        'OUD_vs_Control_Female',
        'OUD_Effect_Male_vs_Female',
        'Control_Effect_Male_vs_Female',
        'Control_Effect_Male_vs_Female_Putamen'
    ]
    
    print(f"\n⚠️  DESeq2 CONTRASTS NOT AVAILABLE IN LEMUR:")
    for contrast in missing_deseq2:
        print(f"     • {contrast}")
    
    return lemur_results

def prepare_lemur_data_for_decoupler(lemur_results):
    """Convert LEMUR results to decoupleR format (wide matrix)"""
    print("\n🔄 Preparing LEMUR data for decoupleR...")
    
    # Create a comprehensive gene universe
    all_genes = set()
    for df in lemur_results.values():
        all_genes.update(df['gene'].values)
    all_genes = sorted(list(all_genes))
    
    # Create wide matrix with coefficients (equivalent to fold changes)
    data_matrix = pd.DataFrame(index=list(lemur_results.keys()), columns=all_genes)
    
    for contrast_name, df in lemur_results.items():
        # Use coefficients as the effect size measure
        gene_coeffs = df.set_index('gene')['coefficient']
        data_matrix.loc[contrast_name, gene_coeffs.index] = gene_coeffs.values
    
    # Fill missing values with 0
    data_matrix = data_matrix.fillna(0).astype(float)
    
    print(f"   📐 Data matrix shape: {data_matrix.shape}")
    print(f"   📊 Contrasts: {len(data_matrix)}")
    print(f"   🧬 Genes: {len(data_matrix.columns)}")
    
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
        
        # Filter by significance
        significant_results = {}
        all_results = {}
        
        for contrast in data_matrix.index:
            # Get p-values for this contrast
            contrast_padj = tf_padj.loc[contrast]
            significant_mask = contrast_padj < 0.05
            
            if significant_mask.sum() > 0:
                significant_tfs = tf_acts.loc[contrast, significant_mask]
                significant_results[contrast] = {
                    'activities': significant_tfs,
                    'padj': contrast_padj[significant_mask],
                    'n_significant': significant_mask.sum()
                }
                contrast_type = "DESeq2-equivalent" if not contrast.startswith('LEMUR_') else "LEMUR-specific"
                print(f"   ✅ {contrast}: {significant_mask.sum()} significant TFs ({contrast_type})")
            else:
                print(f"   ⚠️  {contrast}: No significant TFs")
            
            # Store all results regardless of significance
            all_results[contrast] = {
                'activities': tf_acts.loc[contrast],
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
    """Run pathway analysis using PROGENy"""
    print("\n🛤️ Running Pathway Analysis with PROGENy...")
    
    try:
        # Load PROGENy network
        print("   📡 Loading PROGENy network...")
        progeny = dc.op.progeny(organism='human')
        print(f"   📊 Loaded {len(progeny)} pathway-gene interactions")
        
        # Run ULM analysis
        print("   🔄 Running ULM analysis...")
        pw_acts, pw_padj = dc.mt.ulm(data=data_matrix, net=progeny, verbose=False)
        
        # Filter by significance
        significant_results = {}
        all_results = {}
        
        for contrast in data_matrix.index:
            # Get p-values for this contrast
            contrast_padj = pw_padj.loc[contrast]
            significant_mask = contrast_padj < 0.05
            
            if significant_mask.sum() > 0:
                significant_pathways = pw_acts.loc[contrast, significant_mask]
                significant_results[contrast] = {
                    'activities': significant_pathways,
                    'padj': contrast_padj[significant_mask],
                    'n_significant': significant_mask.sum()
                }
                contrast_type = "DESeq2-equivalent" if not contrast.startswith('LEMUR_') else "LEMUR-specific"
                print(f"   ✅ {contrast}: {significant_mask.sum()} significant pathways ({contrast_type})")
            else:
                print(f"   ⚠️  {contrast}: No significant pathways")
            
            # Store all results regardless of significance
            all_results[contrast] = {
                'activities': pw_acts.loc[contrast],
                'padj': contrast_padj
            }
        
        return significant_results, all_results, progeny
        
    except Exception as e:
        print(f"   ❌ Error in pathway analysis: {str(e)}")
        return {}, {}, None

def run_hallmark_analysis(data_matrix):
    """Run hallmark gene set analysis"""
    print("\n🎯 Running Hallmark Gene Set Analysis...")
    
    try:
        # Load Hallmark gene sets
        print("   📡 Loading Hallmark gene sets...")
        hallmark = dc.op.hallmark(organism='human')
        print(f"   📊 Loaded {len(hallmark)} hallmark gene set interactions")
        
        # Run ULM analysis
        print("   🔄 Running ULM analysis...")
        hm_acts, hm_padj = dc.mt.ulm(data=data_matrix, net=hallmark, verbose=False)
        
        # Filter by significance
        significant_results = {}
        all_results = {}
        
        for contrast in data_matrix.index:
            # Get p-values for this contrast
            contrast_padj = hm_padj.loc[contrast]
            significant_mask = contrast_padj < 0.05
            
            if significant_mask.sum() > 0:
                significant_hallmarks = hm_acts.loc[contrast, significant_mask]
                significant_results[contrast] = {
                    'activities': significant_hallmarks,
                    'padj': contrast_padj[significant_mask],
                    'n_significant': significant_mask.sum()
                }
                contrast_type = "DESeq2-equivalent" if not contrast.startswith('LEMUR_') else "LEMUR-specific"
                print(f"   ✅ {contrast}: {significant_mask.sum()} significant hallmarks ({contrast_type})")
            else:
                print(f"   ⚠️  {contrast}: No significant hallmarks")
            
            # Store all results regardless of significance
            all_results[contrast] = {
                'activities': hm_acts.loc[contrast],
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
        log_pvals = -np.log10(df['pval'].clip(lower=1e-300))  # Clip to avoid inf
        coefficients = df['coefficient']
        
        # Color points based on significance and effect size
        colors = []
        for i, row in df.iterrows():
            if row['pval'] < 0.05 and abs(row['coefficient']) > 0.01:
                colors.append('red' if row['coefficient'] > 0 else 'blue')
            elif row['pval'] < 0.05:
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
        top_genes = df.nlargest(5, 'abs_coefficient')
        for _, gene_row in top_genes.iterrows():
            if gene_row['pval'] < 0.05:
                plt.annotate(gene_row['gene'], 
                           (gene_row['coefficient'], -np.log10(gene_row['pval'])),
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

def generate_comprehensive_report(tf_results, pathway_results, hallmark_results, lemur_results, output_dir):
    """Generate a comprehensive analysis report with DESeq2 compatibility"""
    print("\n📄 Generating comprehensive report...")
    
    report_path = f"{output_dir}/LEMUR_enrichment_analysis_report_deseq2_compatible.md"
    
    # Separate contrasts
    deseq2_contrasts = [c for c in lemur_results.keys() if not c.startswith('LEMUR_')]
    lemur_contrasts = [c for c in lemur_results.keys() if c.startswith('LEMUR_')]
    
    with open(report_path, 'w') as f:
        f.write("# LEMUR Pathway Enrichment Analysis Report (DESeq2-Compatible)\n\n")
        f.write("## Overview\n")
        f.write(f"Analysis of {len(lemur_results)} LEMUR contrasts using decoupleR\n")
        f.write(f"- **DESeq2-Compatible contrasts**: {len(deseq2_contrasts)}\n")
        f.write(f"- **LEMUR-Specific contrasts**: {len(lemur_contrasts)}\n\n")
        
        # Summary statistics
        f.write("## Summary Statistics\n\n")
        f.write("### DESeq2-Compatible Contrasts\n")
        f.write("| Analysis Type | Total Significant Features |\n")
        f.write("|---|---|\n")
        
        if tf_results:
            total_sig_tfs = sum(data['n_significant'] for contrast, data in tf_results.items() 
                              if not contrast.startswith('LEMUR_'))
            f.write(f"| Transcription Factors | {total_sig_tfs} |\n")
        
        if pathway_results:
            total_sig_pathways = sum(data['n_significant'] for contrast, data in pathway_results.items() 
                                   if not contrast.startswith('LEMUR_'))
            f.write(f"| Pathways (PROGENy) | {total_sig_pathways} |\n")
        
        if hallmark_results:
            total_sig_hallmarks = sum(data['n_significant'] for contrast, data in hallmark_results.items() 
                                    if not contrast.startswith('LEMUR_'))
            f.write(f"| Hallmark Gene Sets | {total_sig_hallmarks} |\n")
        
        f.write("\n### LEMUR-Specific Contrasts\n")
        f.write("| Analysis Type | Total Significant Features |\n")
        f.write("|---|---|\n")
        
        if tf_results:
            total_sig_tfs_lemur = sum(data['n_significant'] for contrast, data in tf_results.items() 
                                    if contrast.startswith('LEMUR_'))
            f.write(f"| Transcription Factors | {total_sig_tfs_lemur} |\n")
        
        if pathway_results:
            total_sig_pathways_lemur = sum(data['n_significant'] for contrast, data in pathway_results.items() 
                                         if contrast.startswith('LEMUR_'))
            f.write(f"| Pathways (PROGENy) | {total_sig_pathways_lemur} |\n")
        
        if hallmark_results:
            total_sig_hallmarks_lemur = sum(data['n_significant'] for contrast, data in hallmark_results.items() 
                                          if contrast.startswith('LEMUR_'))
            f.write(f"| Hallmark Gene Sets | {total_sig_hallmarks_lemur} |\n")
        
        # DESeq2-Compatible results
        f.write("\n## DESeq2-Compatible Results\n\n")
        
        for contrast in deseq2_contrasts:
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
        
        # LEMUR-Specific results
        if lemur_contrasts:
            f.write("\n## LEMUR-Specific Results\n\n")
            
            for contrast in lemur_contrasts:
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
            if lemur_contrast in deseq2_contrasts:
                f.write(f"- **{lemur_contrast}** ↔ **{deseq2_file}**\n")
        
        f.write("\n### LEMUR-Unique Insights:\n")
        for contrast in lemur_contrasts:
            f.write(f"- **{contrast}**: Sex-specific OUD vulnerability patterns not captured in DESeq2\n")
    
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
    """Analyze OUD-specific patterns across brain regions"""
    oud_contrasts = [contrast for contrast in tf_results.keys() if 'OUD' in contrast]
    
    if not oud_contrasts:
        print("⚠️ No OUD-specific contrasts found")
        return
    
    print(f"\n🧠 ANALYZING {len(oud_contrasts)} OUD-SPECIFIC CONTRASTS:")
    
    # Find shared OUD features
    all_oud_features = {}
    
    for results_dict, feature_type in [(tf_results, "TF"), (pathway_results, "Pathway"), (hallmark_results, "Hallmark")]:
        if not results_dict:
            continue
            
        oud_features = {}
        for contrast in oud_contrasts:
            if contrast in results_dict:
                for feature in results_dict[contrast]['activities'].index:
                    if feature not in oud_features:
                        oud_features[feature] = []
                    oud_features[feature].append(contrast)
        
        shared_oud = {feature: contrasts for feature, contrasts in oud_features.items() 
                     if len(contrasts) >= 2}
        
        if shared_oud:
            print(f"\n   {feature_type} features shared across OUD contrasts:")
            for feature, contrasts in list(shared_oud.items())[:3]:
                print(f"     • {feature}: {', '.join(contrasts)}")

def main():
    """Main execution function with DESeq2 compatibility"""
    print("🚀 STARTING COMPREHENSIVE PATHWAY ENRICHMENT ANALYSIS - LEMUR VERSION")
    print("📊 DESEQ2-COMPATIBLE CONTRAST NAMING & ANALYSIS")
    print("=" * 80)
    
    if not DECOUPLER_AVAILABLE:
        print("\n❌ decoupleR not available!")
        print("Install with: pip install decoupler omnipath")
        return None
    
    try:
        # Load LEMUR results with DESeq2-compatible naming
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
            if pathway_significant:
                plot_pathway_results(pathway_significant, contrast_name, PLOTS_DIR, "PROGENy_Pathway")
            
            # Hallmark plots
            if hallmark_significant:
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
        
        # Generate comprehensive report
        generate_comprehensive_report(tf_significant, pathway_significant, hallmark_significant, 
                                    lemur_results, REPORTS_DIR)
        
        # Print final summary
        print("\n✅ LEMUR PATHWAY ENRICHMENT ANALYSIS COMPLETED!")
        print("📊 DESEQ2-COMPATIBLE VERSION")
        print("=" * 60)
        
        # Separate summary by contrast type
        deseq2_contrasts = [c for c in lemur_results.keys() if not c.startswith('LEMUR_')]
        lemur_contrasts = [c for c in lemur_results.keys() if c.startswith('LEMUR_')]
        
        print(f"📊 Analyzed {len(lemur_results)} total LEMUR contrasts:")
        print(f"   🔄 DESeq2-compatible: {len(deseq2_contrasts)}")
        print(f"   🧠 LEMUR-specific: {len(lemur_contrasts)}")
        
        if tf_significant:
            total_sig_tfs_deseq2 = sum(data['n_significant'] for c, data in tf_significant.items() if not c.startswith('LEMUR_'))
            total_sig_tfs_lemur = sum(data['n_significant'] for c, data in tf_significant.items() if c.startswith('LEMUR_'))
            print(f"🧬 Transcription Factors: {total_sig_tfs_deseq2} (DESeq2-compatible) + {total_sig_tfs_lemur} (LEMUR-specific)")
        
        if pathway_significant:
            total_sig_pathways_deseq2 = sum(data['n_significant'] for c, data in pathway_significant.items() if not c.startswith('LEMUR_'))
            total_sig_pathways_lemur = sum(data['n_significant'] for c, data in pathway_significant.items() if c.startswith('LEMUR_'))
            print(f"🛤️  Pathways: {total_sig_pathways_deseq2} (DESeq2-compatible) + {total_sig_pathways_lemur} (LEMUR-specific)")
        
        if hallmark_significant:
            total_sig_hallmarks_deseq2 = sum(data['n_significant'] for c, data in hallmark_significant.items() if not c.startswith('LEMUR_'))
            total_sig_hallmarks_lemur = sum(data['n_significant'] for c, data in hallmark_significant.items() if c.startswith('LEMUR_'))
            print(f"🎯 Hallmarks: {total_sig_hallmarks_deseq2} (DESeq2-compatible) + {total_sig_hallmarks_lemur} (LEMUR-specific)")
        
        print(f"\n📁 Results saved to: {OUTPUT_DIR}")
        print("   📊 Plots: /plots/ (DESeq2-compatible indicators)")
        print("   📋 Tables: /tables/ (contrast_type column added)") 
        print("   📄 Reports: /reports/ (includes cross-method comparison guide)")
        
        print("\n🔄 CROSS-METHOD COMPARISON READY!")
        print("📊 Direct DESeq2 equivalents available for:")
        for contrast in deseq2_contrasts:
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
            'deseq2_compatible_contrasts': deseq2_contrasts,
            'lemur_specific_contrasts': lemur_contrasts
        }
        
    except Exception as e:
        print(f"\n❌ ERROR: {str(e)}")
        raise e

if __name__ == "__main__":
    main()