#!/usr/bin/env python3
"""
🛤️ Pathway Enrichment Analysis using decoupleR
GSE225158 - OUD vs Control - Comprehensive Pathway Analysis

This script performs thorough pathway enrichment analysis on DESeq2 results using:
1. Transcription Factor (TF) activity analysis via CollecTRI
2. Pathway activity analysis via PROGENy
3. Hallmark gene set enrichment via MSigDB
4. Custom visualization and reporting

Input: DESeq2 results from 03a_pyDESeq2.py
Output: Enrichment scores, plots, and comprehensive reports
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
DESEQ2_RESULTS_DIR = f"{BASE_DIR}/results/snrna_scvi/pydeseq2_analysis/tables/Full_Dataset"
OUTPUT_DIR = f"{BASE_DIR}/results/snrna_scvi/pathway_enrichment_decoupler"
PLOTS_DIR = f"{OUTPUT_DIR}/plots"
TABLES_DIR = f"{OUTPUT_DIR}/tables"
REPORTS_DIR = f"{OUTPUT_DIR}/reports"

# Create directories
for directory in [OUTPUT_DIR, PLOTS_DIR, TABLES_DIR, REPORTS_DIR]:
    os.makedirs(directory, exist_ok=True)

print("🛤️ PATHWAY ENRICHMENT ANALYSIS USING DECOUPLER")
print("=" * 50)
print(f"Input DESeq2 results: {DESEQ2_RESULTS_DIR}")
print(f"Output directory: {OUTPUT_DIR}")

# ============================================================================
# 🔄 DATA LOADING AND PREPARATION
# ============================================================================

def load_deseq2_results():
    """Load all DESeq2 contrast results"""
    print("\n📂 Loading DESeq2 results...")
    
    # Define all contrast files
    contrast_files = {
        'OUD_vs_Control_Pooled': '01_OUD_vs_Control_Pooled_results.csv',
        'OUD_vs_Control_Putamen': '02_OUD_vs_Control_Putamen_results.csv',
        'OUD_vs_Control_Caudate': '03_OUD_vs_Control_Caudate_results.csv',
        'OUD_vs_Control_Male': '04_OUD_vs_Control_Male_results.csv',
        'OUD_vs_Control_Female': '05_OUD_vs_Control_Female_results.csv',
        'OUD_Effect_Male_vs_Female': '06_OUD_Effect_Male_vs_Female_results.csv',
        'OUD_Effect_Male_vs_Female_Putamen': '07_OUD_Effect_Male_vs_Female_Putamen_results.csv',
        'OUD_Effect_Male_vs_Female_Caudate': '08_OUD_Effect_Male_vs_Female_Caudate_results.csv',
        'Control_Effect_Male_vs_Female': '09_Control_Effect_Male_vs_Female_results.csv',
        'Control_Effect_Male_vs_Female_Putamen': '10_Control_Effect_Male_vs_Female_Putamen_results.csv'
    }
    
    deseq2_results = {}
    
    for contrast_name, filename in contrast_files.items():
        filepath = f"{DESEQ2_RESULTS_DIR}/{filename}"
        
        if os.path.exists(filepath):
            df = pd.read_csv(filepath)
            
            # Ensure required columns exist
            if all(col in df.columns for col in ['gene', 'stat', 'log2FoldChange', 'pvalue', 'padj']):
                # Remove genes with NaN statistics
                df = df.dropna(subset=['stat', 'log2FoldChange'])
                deseq2_results[contrast_name] = df
                print(f"   ✅ Loaded {contrast_name}: {len(df)} genes")
            else:
                print(f"   ❌ Missing required columns in {contrast_name}")
        else:
            print(f"   ❌ File not found: {filepath}")
    
    print(f"\n📊 Successfully loaded {len(deseq2_results)} contrasts")
    return deseq2_results

def prepare_data_for_decoupler(deseq2_results):
    """Convert DESeq2 results to decoupleR format (wide matrix)"""
    print("\n🔄 Preparing data for decoupleR...")
    
    # Create a comprehensive gene universe
    all_genes = set()
    for df in deseq2_results.values():
        all_genes.update(df['gene'].values)
    all_genes = sorted(list(all_genes))
    
    # Create wide matrix with t-statistics
    data_matrix = pd.DataFrame(index=list(deseq2_results.keys()), columns=all_genes)
    
    for contrast_name, df in deseq2_results.items():
        # Use t-statistics (stat column) as recommended by decoupleR
        gene_stats = df.set_index('gene')['stat']
        data_matrix.loc[contrast_name, gene_stats.index] = gene_stats.values
    
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
                print(f"   ✅ {contrast}: {significant_mask.sum()} significant TFs")
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
    plt.title(f'Top Transcription Factor Activities\n{contrast_name.replace("_", " ")}')
    plt.axvline(x=0, color='black', linestyle='-', alpha=0.3)
    
    # Add value labels
    for i, (bar, value) in enumerate(zip(bars, top_activities.values)):
        plt.text(value + (0.1 if value > 0 else -0.1), i, f'{value:.2f}', 
                va='center', ha='left' if value > 0 else 'right')
    
    plt.tight_layout()
    
    # Save plot
    plot_path = f"{output_dir}/TF_activities_{contrast_name}.png"
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
                print(f"   ✅ {contrast}: {significant_mask.sum()} significant pathways")
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
                print(f"   ✅ {contrast}: {significant_mask.sum()} significant hallmarks")
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
    plt.title(f'{analysis_type} Activities\n{contrast_name.replace("_", " ")}')
    plt.axvline(x=0, color='black', linestyle='-', alpha=0.3)
    
    # Add value labels
    for i, (bar, value) in enumerate(zip(bars, activities_sorted.values)):
        plt.text(value + (0.1 if value > 0 else -0.1), i, f'{value:.2f}', 
                va='center', ha='left' if value > 0 else 'right')
    
    plt.tight_layout()
    
    # Save plot
    plot_path = f"{output_dir}/{analysis_type}_activities_{contrast_name}.png"
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"   📊 {analysis_type} plot saved: {plot_path}")

# ============================================================================
# 📊 COMPREHENSIVE VISUALIZATION
# ============================================================================

def create_heatmap_comparison(tf_results, pathway_results, hallmark_results, output_dir):
    """Create heatmaps comparing results across all contrasts"""
    print("\n📊 Creating comparison heatmaps...")
    
    # Transcription Factors
    if tf_results:
        tf_matrix = pd.DataFrame()
        for contrast, data in tf_results.items():
            tf_matrix[contrast] = data['activities']
        
        if not tf_matrix.empty:
            plt.figure(figsize=(15, 10))
            sns.heatmap(tf_matrix.T, cmap='RdBu_r', center=0, 
                       xticklabels=True, yticklabels=True, cbar_kws={'label': 'TF Activity Score'})
            plt.title('Transcription Factor Activities Across Contrasts')
            plt.xlabel('Contrasts')
            plt.ylabel('Transcription Factors')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            
            heatmap_path = f"{output_dir}/TF_activities_heatmap.png"
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
            plt.title('Pathway Activities Across Contrasts (PROGENy)')
            plt.xlabel('Contrasts')
            plt.ylabel('Pathways')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            
            heatmap_path = f"{output_dir}/Pathway_activities_heatmap.png"
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
            plt.title('Hallmark Gene Set Activities Across Contrasts')
            plt.xlabel('Contrasts')
            plt.ylabel('Hallmark Gene Sets')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            
            heatmap_path = f"{output_dir}/Hallmark_activities_heatmap.png"
            plt.savefig(heatmap_path, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"   📊 Hallmark heatmap saved: {heatmap_path}")

def create_tutorial_plots(deseq2_results, tf_results, pathway_results, hallmark_results, networks, output_dir):
    """Create all plots from the decoupleR tutorial"""
    print("\n🎨 Creating tutorial-style plots...")
    
    tf_net, pathway_net, hallmark_net = networks
    
    for contrast_name, df in deseq2_results.items():
        if df.empty:
            continue
            
        print(f"\n   📊 Creating plots for {contrast_name}...")
        
        try:
            # 1. Basic volcano plot (tutorial style)
            if not df.empty:
                plt.figure(figsize=(8, 6))
                plt.scatter(df['log2FoldChange'], -np.log10(df['pvalue'].clip(lower=1e-300)), 
                           c='lightblue', alpha=0.6, s=15)
                plt.axhline(y=-np.log10(0.05), color='red', linestyle='--', alpha=0.7)
                plt.axvline(x=0, color='black', linestyle='-', alpha=0.3)
                plt.xlabel('Log2 Fold Change')
                plt.ylabel('-Log10 P-value')
                plt.title(f'Volcano Plot\n{contrast_name.replace("_", " ")}')
                
                volcano_basic_path = f"{output_dir}/volcano_basic_{contrast_name}.png"
                plt.savefig(volcano_basic_path, dpi=300, bbox_inches='tight')
                plt.close()
                print(f"      🌋 Basic volcano plot: {volcano_basic_path}")
            
            # 2. TF activity barplot (tutorial style)
            if contrast_name in tf_results and tf_results[contrast_name]['activities'].size > 0:
                activities = tf_results[contrast_name]['activities']
                top_tfs = activities.abs().nlargest(15)
                
                plt.figure(figsize=(6, 8))
                colors = ['red' if x > 0 else 'blue' for x in top_tfs.values]
                plt.barh(range(len(top_tfs)), top_tfs.values, color=colors, alpha=0.7)
                plt.yticks(range(len(top_tfs)), top_tfs.index)
                plt.xlabel('TF Activity Score')
                plt.title(f'Transcription Factor Activities\n{contrast_name.replace("_", " ")}')
                plt.axvline(x=0, color='black', linestyle='-', alpha=0.3)
                plt.tight_layout()
                
                tf_bar_path = f"{output_dir}/tf_barplot_{contrast_name}.png"
                plt.savefig(tf_bar_path, dpi=300, bbox_inches='tight')
                plt.close()
                print(f"      🧬 TF barplot: {tf_bar_path}")
            
            # 3. Pathway activity barplot (tutorial style)
            if contrast_name in pathway_results and pathway_results[contrast_name]['activities'].size > 0:
                pw_activities = pathway_results[contrast_name]['activities']
                
                plt.figure(figsize=(6, 6))
                colors = ['red' if x > 0 else 'blue' for x in pw_activities.values]
                plt.barh(range(len(pw_activities)), pw_activities.values, color=colors, alpha=0.7)
                plt.yticks(range(len(pw_activities)), pw_activities.index)
                plt.xlabel('Pathway Activity Score')
                plt.title(f'Pathway Activities (PROGENy)\n{contrast_name.replace("_", " ")}')
                plt.axvline(x=0, color='black', linestyle='-', alpha=0.3)
                plt.tight_layout()
                
                pw_bar_path = f"{output_dir}/pathway_barplot_{contrast_name}.png"
                plt.savefig(pw_bar_path, dpi=300, bbox_inches='tight')
                plt.close()
                print(f"      🛤️ Pathway barplot: {pw_bar_path}")
            
            # 4. Pathway dotplot (tutorial style)
            if contrast_name in pathway_results and pathway_results[contrast_name]['activities'].size > 0:
                pw_acts = pathway_results[contrast_name]['activities']
                pw_padj = pathway_results[contrast_name]['padj']
                
                # Create dataframe for dotplot
                plot_df = pd.DataFrame({
                    'pathway': pw_acts.index,
                    'score': pw_acts.values,
                    'pvalue': pw_padj.values,
                    'logpval': -np.log10(pw_padj.values.clip(lower=1e-10))
                })
                
                plt.figure(figsize=(8, 6))
                scatter = plt.scatter(plot_df['score'], plot_df['pathway'], 
                                    s=plot_df['logpval']*20, c=plot_df['score'], 
                                    cmap='RdBu_r', alpha=0.8, edgecolors='black', linewidth=0.5)
                plt.xlabel('Pathway Activity Score')
                plt.title(f'Pathway Activities (Dotplot)\n{contrast_name.replace("_", " ")}')
                plt.axvline(x=0, color='black', linestyle='-', alpha=0.3)
                plt.colorbar(scatter, label='Activity Score')
                plt.tight_layout()
                
                pw_dot_path = f"{output_dir}/pathway_dotplot_{contrast_name}.png"
                plt.savefig(pw_dot_path, dpi=300, bbox_inches='tight')
                plt.close()
                print(f"      🎯 Pathway dotplot: {pw_dot_path}")
            
            # 5. Hallmark barplot (tutorial style)
            if contrast_name in hallmark_results and hallmark_results[contrast_name]['activities'].size > 0:
                hm_activities = hallmark_results[contrast_name]['activities']
                top_hallmarks = hm_activities.abs().nlargest(20)
                
                plt.figure(figsize=(10, 8))
                colors = ['red' if x > 0 else 'blue' for x in top_hallmarks.values]
                plt.barh(range(len(top_hallmarks)), top_hallmarks.values, color=colors, alpha=0.7)
                plt.yticks(range(len(top_hallmarks)), 
                          [name.replace('_', ' ') for name in top_hallmarks.index])
                plt.xlabel('Hallmark Activity Score')
                plt.title(f'Hallmark Gene Set Activities\n{contrast_name.replace("_", " ")}')
                plt.axvline(x=0, color='black', linestyle='-', alpha=0.3)
                plt.tight_layout()
                
                hm_bar_path = f"{output_dir}/hallmark_barplot_{contrast_name}.png"
                plt.savefig(hm_bar_path, dpi=300, bbox_inches='tight')
                plt.close()
                print(f"      🎯 Hallmark barplot: {hm_bar_path}")
                
        except Exception as e:
            print(f"      ❌ Error creating plots for {contrast_name}: {str(e)}")


def create_volcano_plots(deseq2_results, networks, output_dir):
    """Create volcano plots with pathway overlays"""
    print("\n🌋 Creating enhanced volcano plots...")
    
    tf_net, pathway_net, hallmark_net = networks
    
    for contrast_name, df in deseq2_results.items():
        if tf_net is not None and not df.empty:
            try:
                # Enhanced volcano plot with pathway overlay
                plt.figure(figsize=(10, 8))
                
                # Base volcano plot
                scatter = plt.scatter(df['log2FoldChange'], -np.log10(df['pvalue'].clip(lower=1e-300)), 
                                    c='lightgray', alpha=0.6, s=20)
                
                # Highlight significant genes
                sig_mask = (df['padj'] < 0.05) & (df['log2FoldChange'].abs() > 0.5)
                if sig_mask.sum() > 0:
                    plt.scatter(df.loc[sig_mask, 'log2FoldChange'], 
                              -np.log10(df.loc[sig_mask, 'pvalue'].clip(lower=1e-300)),
                              c='red', alpha=0.8, s=30, label=f'Significant ({sig_mask.sum()})')
                
                plt.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
                plt.axvline(x=0, color='black', linestyle='-', alpha=0.3)
                plt.xlabel('Log2 Fold Change')
                plt.ylabel('-Log10 P-value')
                plt.title(f'Enhanced Volcano Plot: {contrast_name.replace("_", " ")}')
                plt.legend()
                
                volcano_path = f"{output_dir}/volcano_enhanced_{contrast_name}.png"
                plt.savefig(volcano_path, dpi=300, bbox_inches='tight')
                plt.close()
                
                print(f"   🌋 Enhanced volcano plot saved: {volcano_path}")
                
            except Exception as e:
                print(f"   ❌ Error creating volcano plot for {contrast_name}: {str(e)}")


def create_source_target_plots(deseq2_results, pathway_results, pathway_net, output_dir):
    """Create source-target plots for top pathways (tutorial style)"""
    print("\n🎯 Creating source-target plots...")
    
    if pathway_net is None:
        return
    
    for contrast_name, df in deseq2_results.items():
        if contrast_name not in pathway_results or df.empty:
            continue
            
        try:
            pw_activities = pathway_results[contrast_name]['activities']
            top_pathways = pw_activities.abs().nlargest(3)
            
            for pathway_name in top_pathways.index:
                # Get pathway targets
                pathway_targets = pathway_net[pathway_net['source'] == pathway_name]
                
                if len(pathway_targets) == 0:
                    continue
                
                # Merge with DESeq2 results
                merged_data = pathway_targets.merge(df, left_on='target', right_on='gene', how='inner')
                
                if len(merged_data) < 5:
                    continue
                
                # Create source-target plot
                plt.figure(figsize=(8, 6))
                plt.scatter(merged_data['weight'], merged_data['stat'], 
                           c=merged_data['stat'], cmap='RdBu_r', s=50, alpha=0.7)
                plt.xlabel('Pathway Weight')
                plt.ylabel('Gene t-statistic')
                plt.title(f'{pathway_name} Target Genes\n{contrast_name.replace("_", " ")}')
                plt.axhline(y=0, color='black', linestyle='-', alpha=0.3)
                plt.axvline(x=0, color='black', linestyle='-', alpha=0.3)
                
                # Add top gene labels
                top_genes = merged_data.nlargest(5, 'stat')
                for _, gene_row in top_genes.iterrows():
                    plt.annotate(gene_row['gene'], 
                               (gene_row['weight'], gene_row['stat']),
                               xytext=(5, 5), textcoords='offset points',
                               fontsize=8, alpha=0.8)
                
                plt.colorbar(label='t-statistic')
                plt.tight_layout()
                
                st_path = f"{output_dir}/source_target_{pathway_name}_{contrast_name}.png"
                plt.savefig(st_path, dpi=300, bbox_inches='tight')
                plt.close()
                
                print(f"   🎯 Source-target plot: {st_path}")
                
        except Exception as e:
            print(f"   ❌ Error creating source-target plots for {contrast_name}: {str(e)}")


def create_leading_edge_plots(deseq2_results, pathway_results, pathway_net, output_dir):
    """Create leading edge plots (tutorial style)"""
    print("\n📈 Creating leading edge plots...")
    
    if pathway_net is None:
        return
    
    for contrast_name, df in deseq2_results.items():
        if contrast_name not in pathway_results or df.empty:
            continue
            
        try:
            pw_activities = pathway_results[contrast_name]['activities']
            top_pathways = pw_activities.abs().nlargest(2)
            
            for pathway_name in top_pathways.index:
                # Get pathway targets
                pathway_targets = pathway_net[pathway_net['source'] == pathway_name]
                
                if len(pathway_targets) == 0:
                    continue
                
                # Merge with DESeq2 results and sort by t-statistic
                merged_data = pathway_targets.merge(df, left_on='target', right_on='gene', how='inner')
                merged_data = merged_data.sort_values('stat', ascending=False)
                
                if len(merged_data) < 10:
                    continue
                
                # Create leading edge plot
                plt.figure(figsize=(10, 6))
                
                # Plot running enrichment score
                enrichment_score = np.cumsum(merged_data['weight'] * np.sign(merged_data['stat']))
                plt.plot(range(len(enrichment_score)), enrichment_score, 'b-', linewidth=2)
                
                plt.xlabel('Gene Rank')
                plt.ylabel('Running Enrichment Score')
                plt.title(f'Leading Edge: {pathway_name}\n{contrast_name.replace("_", " ")}')
                plt.axhline(y=0, color='black', linestyle='-', alpha=0.3)
                
                # Mark leading edge
                max_idx = np.argmax(np.abs(enrichment_score))
                plt.axvline(x=max_idx, color='red', linestyle='--', alpha=0.7, 
                           label=f'Leading Edge ({max_idx+1} genes)')
                plt.legend()
                
                plt.tight_layout()
                
                le_path = f"{output_dir}/leading_edge_{pathway_name}_{contrast_name}.png"
                plt.savefig(le_path, dpi=300, bbox_inches='tight')
                plt.close()
                
                print(f"   📈 Leading edge plot: {le_path}")
                
        except Exception as e:
            print(f"   ❌ Error creating leading edge plots for {contrast_name}: {str(e)}")

# ============================================================================
# 💾 RESULTS EXPORT
# ============================================================================

def save_enrichment_results(tf_results, pathway_results, hallmark_results, output_dir):
    """Save enrichment results to files"""
    print("\n💾 Saving enrichment results...")
    
    # Save TF results
    if tf_results:
        tf_summary = []
        for contrast, data in tf_results.items():
            for tf, activity in data['activities'].items():
                tf_summary.append({
                    'contrast': contrast,
                    'tf': tf,
                    'activity_score': activity,
                    'padj': data['padj'][tf],
                    'analysis_type': 'Transcription_Factor'
                })
        
        if tf_summary:
            tf_df = pd.DataFrame(tf_summary)
            tf_path = f"{output_dir}/transcription_factor_results.csv"
            tf_df.to_csv(tf_path, index=False)
            print(f"   💾 TF results saved: {tf_path}")
    
    # Save pathway results
    if pathway_results:
        pathway_summary = []
        for contrast, data in pathway_results.items():
            for pathway, activity in data['activities'].items():
                pathway_summary.append({
                    'contrast': contrast,
                    'pathway': pathway,
                    'activity_score': activity,
                    'padj': data['padj'][pathway],
                    'analysis_type': 'PROGENy_Pathway'
                })
        
        if pathway_summary:
            pathway_df = pd.DataFrame(pathway_summary)
            pathway_path = f"{output_dir}/progeny_pathway_results.csv"
            pathway_df.to_csv(pathway_path, index=False)
            print(f"   💾 Pathway results saved: {pathway_path}")
    
    # Save hallmark results
    if hallmark_results:
        hallmark_summary = []
        for contrast, data in hallmark_results.items():
            for hallmark, activity in data['activities'].items():
                hallmark_summary.append({
                    'contrast': contrast,
                    'hallmark': hallmark,
                    'activity_score': activity,
                    'padj': data['padj'][hallmark],
                    'analysis_type': 'Hallmark_Gene_Set'
                })
        
        if hallmark_summary:
            hallmark_df = pd.DataFrame(hallmark_summary)
            hallmark_path = f"{output_dir}/hallmark_geneset_results.csv"
            hallmark_df.to_csv(hallmark_path, index=False)
            print(f"   💾 Hallmark results saved: {hallmark_path}")
    
    # Create comprehensive summary
    all_results = []
    for results_dict, analysis_type in [(tf_results, 'Transcription_Factor'), 
                                       (pathway_results, 'PROGENy_Pathway'), 
                                       (hallmark_results, 'Hallmark_Gene_Set')]:
        if results_dict:
            for contrast, data in results_dict.items():
                all_results.append({
                    'contrast': contrast,
                    'analysis_type': analysis_type,
                    'n_significant': data['n_significant'],
                    'top_positive': data['activities'].nlargest(1).index[0] if len(data['activities']) > 0 else 'None',
                    'top_positive_score': data['activities'].nlargest(1).iloc[0] if len(data['activities']) > 0 else 0,
                    'top_negative': data['activities'].nsmallest(1).index[0] if len(data['activities']) > 0 else 'None',
                    'top_negative_score': data['activities'].nsmallest(1).iloc[0] if len(data['activities']) > 0 else 0
                })
    
    if all_results:
        summary_df = pd.DataFrame(all_results)
        summary_path = f"{output_dir}/enrichment_analysis_summary.csv"
        summary_df.to_csv(summary_path, index=False)
        print(f"   💾 Summary saved: {summary_path}")

def generate_comprehensive_report(tf_results, pathway_results, hallmark_results, deseq2_results, output_dir):
    """Generate a comprehensive analysis report"""
    print("\n📄 Generating comprehensive report...")
    
    report_path = f"{output_dir}/pathway_enrichment_report.md"
    
    with open(report_path, 'w') as f:
        f.write("# Pathway Enrichment Analysis Report\n")
        f.write("## GSE225158 - OUD vs Control Analysis\n\n")
        
        f.write("### Analysis Overview\n")
        f.write(f"- **Total Contrasts Analyzed**: {len(deseq2_results)}\n")
        f.write("- **Enrichment Methods Used**:\n")
        f.write("  - Transcription Factor Activity (CollecTRI)\n")
        f.write("  - Pathway Activity (PROGENy)\n")
        f.write("  - Hallmark Gene Set Enrichment (MSigDB)\n\n")
        
        # TF Results Summary
        if tf_results:
            f.write("### Transcription Factor Analysis Summary\n")
            total_sig_tfs = sum(data['n_significant'] for data in tf_results.values())
            f.write(f"- **Total Significant TFs**: {total_sig_tfs}\n")
            f.write("- **Significant TFs by Contrast**:\n")
            for contrast, data in tf_results.items():
                f.write(f"  - {contrast}: {data['n_significant']} TFs\n")
            f.write("\n")
        
        # Pathway Results Summary
        if pathway_results:
            f.write("### Pathway Analysis Summary (PROGENy)\n")
            total_sig_pathways = sum(data['n_significant'] for data in pathway_results.values())
            f.write(f"- **Total Significant Pathways**: {total_sig_pathways}\n")
            f.write("- **Significant Pathways by Contrast**:\n")
            for contrast, data in pathway_results.items():
                f.write(f"  - {contrast}: {data['n_significant']} pathways\n")
            f.write("\n")
        
        # Hallmark Results Summary
        if hallmark_results:
            f.write("### Hallmark Gene Set Analysis Summary\n")
            total_sig_hallmarks = sum(data['n_significant'] for data in hallmark_results.values())
            f.write(f"- **Total Significant Hallmarks**: {total_sig_hallmarks}\n")
            f.write("- **Significant Hallmarks by Contrast**:\n")
            for contrast, data in hallmark_results.items():
                f.write(f"  - {contrast}: {data['n_significant']} hallmarks\n")
            f.write("\n")
        
        f.write("### Key Findings\n")
        f.write("#### Most Activated Across Contrasts\n")
        
        # Find most commonly activated features
        all_activations = {}
        for analysis_name, results_dict in [("TF", tf_results), ("Pathway", pathway_results), ("Hallmark", hallmark_results)]:
            if results_dict:
                feature_counts = {}
                for contrast, data in results_dict.items():
                    for feature, score in data['activities'].items():
                        if score > 0:  # Activated
                            feature_counts[feature] = feature_counts.get(feature, 0) + 1
                
                if feature_counts:
                    top_features = sorted(feature_counts.items(), key=lambda x: x[1], reverse=True)[:5]
                    f.write(f"**{analysis_name} Analysis:**\n")
                    for feature, count in top_features:
                        f.write(f"- {feature}: activated in {count} contrasts\n")
                    f.write("\n")
        
        f.write("### Files Generated\n")
        f.write("- Individual contrast plots in `/plots/`\n")
        f.write("- Comprehensive heatmaps comparing all contrasts\n")
        f.write("- Detailed results tables in `/tables/`\n")
        f.write("- Summary statistics and this report in `/reports/`\n\n")
        
        f.write("### Analysis Parameters\n")
        f.write("- **Significance Threshold**: p-adj < 0.05\n")
        f.write("- **Enrichment Method**: Univariate Linear Model (ULM)\n")
        f.write("- **Input Statistics**: t-statistics from DESeq2\n")
        f.write("- **Networks Used**: CollecTRI, PROGENy, MSigDB Hallmarks\n")
    
    print(f"   📄 Report saved: {report_path}")

# ============================================================================
# 🔍 INTEGRATED EXPLORATION FUNCTIONS
# ============================================================================

def explore_contrast_results(tf_results, pathway_results, hallmark_results, contrast_name):
    """Explore results for a specific contrast"""
    print(f"\n🔍 Exploring results for: {contrast_name}")
    print("=" * 50)
    
    results_dict = {
        'Transcription Factor': tf_results,
        'Pathway': pathway_results, 
        'Hallmark': hallmark_results
    }
    
    for analysis_type, results in results_dict.items():
        if not results or contrast_name not in results:
            continue
            
        data = results[contrast_name]
        activities = data['activities']
        
        if len(activities) > 0:
            print(f"\n{analysis_type.upper()} Analysis:")
            print(f"   Significant features: {len(activities)}")
            
            # Top activated
            top_activated = activities.nlargest(3)
            print("   Top activated:")
            for feature, score in top_activated.items():
                print(f"     • {feature}: {score:.3f}")
            
            # Top repressed
            top_repressed = activities.nsmallest(3)
            print("   Top repressed:")
            for feature, score in top_repressed.items():
                print(f"     • {feature}: {score:.3f}")

def find_consistent_features(results_dict, analysis_type='tf', threshold=3):
    """Find features that are consistently activated/repressed across contrasts"""
    print(f"\n🎯 Finding consistent {analysis_type} features...")
    
    if not results_dict:
        print(f"   ❌ No {analysis_type} results found")
        return None
    
    # Count activations and repressions
    feature_stats = []
    all_features = set()
    
    # Collect all features
    for contrast, data in results_dict.items():
        all_features.update(data['activities'].keys())
    
    for feature in all_features:
        activated_contrasts = []
        repressed_contrasts = []
        all_scores = []
        
        for contrast, data in results_dict.items():
            if feature in data['activities']:
                score = data['activities'][feature]
                all_scores.append(score)
                if score > 0:
                    activated_contrasts.append(contrast)
                else:
                    repressed_contrasts.append(contrast)
        
        if len(activated_contrasts) >= threshold or len(repressed_contrasts) >= threshold:
            feature_stats.append({
                'feature': feature,
                'n_activated': len(activated_contrasts),
                'n_repressed': len(repressed_contrasts),
                'activated_contrasts': activated_contrasts,
                'repressed_contrasts': repressed_contrasts,
                'mean_activity': np.mean(all_scores),
                'max_activity': max(all_scores),
                'min_activity': min(all_scores)
            })
    
    if feature_stats:
        consistent_df = pd.DataFrame(feature_stats)
        consistent_df = consistent_df.sort_values('mean_activity', key=abs, ascending=False)
        
        print(f"   ✅ Found {len(consistent_df)} consistent features")
        print("\n   Top consistent features:")
        for _, row in consistent_df.head(5).iterrows():
            direction = "activated" if row['mean_activity'] > 0 else "repressed"
            count = row['n_activated'] if row['mean_activity'] > 0 else row['n_repressed']
            print(f"     • {row['feature']}: {direction} in {count} contrasts (mean: {row['mean_activity']:.3f})")
        
        return consistent_df
    else:
        print(f"   ⚠️  No features found in ≥{threshold} contrasts")
        return None

def analyze_oud_specific_patterns(tf_results, pathway_results, hallmark_results):
    """Analyze patterns specific to OUD vs Control contrasts"""
    print("\n🧠 Analyzing OUD-specific patterns...")
    
    # Define OUD vs Control contrasts
    oud_contrasts = [
        'OUD_vs_Control_Pooled',
        'OUD_vs_Control_Putamen', 
        'OUD_vs_Control_Caudate',
        'OUD_vs_Control_Male',
        'OUD_vs_Control_Female'
    ]
    
    oud_patterns = {}
    
    for analysis_name, results_dict in [("Transcription Factor", tf_results), 
                                       ("Pathway", pathway_results), 
                                       ("Hallmark", hallmark_results)]:
        if not results_dict:
            continue
        
        # Find features significant in multiple OUD contrasts
        feature_counts = {}
        oud_data = {}
        
        for contrast in oud_contrasts:
            if contrast in results_dict:
                oud_data[contrast] = results_dict[contrast]
                for feature in results_dict[contrast]['activities'].keys():
                    feature_counts[feature] = feature_counts.get(feature, 0) + 1
        
        multi_contrast_features = [f for f, count in feature_counts.items() if count >= 2]
        
        if multi_contrast_features:
            print(f"\n{analysis_name} features in multiple OUD contrasts:")
            for feature in multi_contrast_features[:5]:
                contrasts_with_feature = []
                scores = []
                for contrast, data in oud_data.items():
                    if feature in data['activities']:
                        contrasts_with_feature.append(contrast)
                        scores.append(data['activities'][feature])
                
                mean_activity = np.mean(scores)
                print(f"   • {feature}: {len(contrasts_with_feature)} contrasts, mean activity: {mean_activity:.3f}")
            
            oud_patterns[analysis_name] = {
                'multi_contrast_features': multi_contrast_features,
                'top_features': multi_contrast_features[:10]
            }
    
    return oud_patterns

def create_feature_comparison_plot(results_dict, feature_name, analysis_type, output_dir):
    """Create a comparison plot for a specific feature across contrasts"""
    
    if not results_dict:
        print(f"   ❌ No {analysis_type} results found")
        return
    
    feature_data = []
    for contrast, data in results_dict.items():
        if feature_name in data['activities']:
            feature_data.append({
                'contrast': contrast,
                'activity_score': data['activities'][feature_name]
            })
    
    if not feature_data:
        print(f"   ❌ Feature '{feature_name}' not found in {analysis_type} results")
        return
    
    feature_df = pd.DataFrame(feature_data)
    feature_df = feature_df.sort_values('activity_score')
    
    # Create the plot
    plt.figure(figsize=(12, 6))
    
    colors = ['red' if x > 0 else 'blue' for x in feature_df['activity_score']]
    bars = plt.barh(range(len(feature_df)), feature_df['activity_score'], 
                   color=colors, alpha=0.7)
    
    plt.yticks(range(len(feature_df)), 
              [c.replace('_', ' ') for c in feature_df['contrast']])
    plt.xlabel('Activity Score')
    plt.title(f'{feature_name} Activity Across Contrasts\n({analysis_type.upper()} Analysis)')
    plt.axvline(x=0, color='black', linestyle='-', alpha=0.3)
    
    # Add value labels
    for i, (bar, value) in enumerate(zip(bars, feature_df['activity_score'])):
        plt.text(value + (0.05 if value > 0 else -0.05), i, f'{value:.2f}', 
                va='center', ha='left' if value > 0 else 'right', fontsize=9)
    
    plt.tight_layout()
    
    plot_path = f"{output_dir}/feature_comparison_{analysis_type}_{feature_name.replace('/', '_')}.png"
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"   📊 Plot saved: {plot_path}")

def interactive_exploration(tf_results, pathway_results, hallmark_results, deseq2_results):
    """Interactive exploration of pathway enrichment results"""
    print("\n🔍 INTERACTIVE EXPLORATION MODE")
    print("=" * 50)
    
    available_contrasts = list(deseq2_results.keys())
    
    while True:
        print("\n📋 Exploration Options:")
        print("1. Explore specific contrast")
        print("2. Find consistent TF features")
        print("3. Find consistent pathway features") 
        print("4. Find consistent hallmark features")
        print("5. Analyze OUD-specific patterns")
        print("6. Create feature comparison plot")
        print("7. Show dataset overview")
        print("0. Exit exploration")
        
        try:
            choice = input("\nSelect option (0-7): ").strip()
            
            if choice == "0":
                print("👋 Exiting exploration mode")
                break
                
            elif choice == "1":
                print("\nAvailable contrasts:")
                for i, contrast in enumerate(available_contrasts, 1):
                    print(f"{i}. {contrast}")
                
                try:
                    idx = int(input("Select contrast number: ")) - 1
                    if 0 <= idx < len(available_contrasts):
                        explore_contrast_results(tf_results, pathway_results, hallmark_results, 
                                               available_contrasts[idx])
                except (ValueError, IndexError):
                    print("Invalid selection")
            
            elif choice == "2":
                threshold = int(input("Minimum number of contrasts (default 3): ") or "3")
                find_consistent_features(tf_results, 'TF', threshold)
            
            elif choice == "3":
                threshold = int(input("Minimum number of contrasts (default 3): ") or "3")
                find_consistent_features(pathway_results, 'Pathway', threshold)
            
            elif choice == "4":
                threshold = int(input("Minimum number of contrasts (default 3): ") or "3")
                find_consistent_features(hallmark_results, 'Hallmark', threshold)
            
            elif choice == "5":
                analyze_oud_specific_patterns(tf_results, pathway_results, hallmark_results)
            
            elif choice == "6":
                print("Analysis types: 1=TF, 2=Pathway, 3=Hallmark")
                analysis_choice = input("Select analysis type (1-3): ").strip()
                
                results_map = {'1': ('TF', tf_results), '2': ('Pathway', pathway_results), 
                              '3': ('Hallmark', hallmark_results)}
                
                if analysis_choice in results_map:
                    analysis_type, results_dict = results_map[analysis_choice]
                    feature_name = input("Feature name: ").strip()
                    if feature_name:
                        create_feature_comparison_plot(results_dict, feature_name, analysis_type, PLOTS_DIR)
            
            elif choice == "7":
                print("\n📊 Dataset Overview:")
                for analysis_name, results_dict in [("Transcription Factor", tf_results), 
                                                   ("Pathway", pathway_results), 
                                                   ("Hallmark", hallmark_results)]:
                    if results_dict:
                        total_features = set()
                        total_significant = sum(data['n_significant'] for data in results_dict.values())
                        for data in results_dict.values():
                            total_features.update(data['activities'].keys())
                        
                        print(f"\n{analysis_name}:")
                        print(f"  Contrasts: {len(results_dict)}")
                        print(f"  Unique features: {len(total_features)}")
                        print(f"  Total significant results: {total_significant}")
            
            else:
                print("Invalid choice")
                
        except KeyboardInterrupt:
            print("\n👋 Exiting exploration mode")
            break
        except Exception as e:
            print(f"❌ Error: {e}")

# ============================================================================
# 🚀 MAIN EXECUTION
# ============================================================================

def main():
    """Main execution function"""
    print("🚀 STARTING COMPREHENSIVE PATHWAY ENRICHMENT ANALYSIS")
    print("=" * 60)
    
    if not DECOUPLER_AVAILABLE:
        print("\n❌ decoupleR not available!")
        print("Install with: pip install decoupler omnipath")
        return None
    
    try:
        # Load DESeq2 results
        deseq2_results = load_deseq2_results()
        if not deseq2_results:
            print("❌ No DESeq2 results found!")
            return None
        
        # Prepare data for decoupleR
        data_matrix = prepare_data_for_decoupler(deseq2_results)
        
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
        for contrast_name in deseq2_results.keys():
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
        
        # Tutorial-style plots (decoupleR style)
        create_tutorial_plots(deseq2_results, tf_significant, pathway_significant, hallmark_significant, 
                             (tf_net, pathway_net, hallmark_net), PLOTS_DIR)
        
        # Enhanced volcano plots
        create_volcano_plots(deseq2_results, (tf_net, pathway_net, hallmark_net), PLOTS_DIR)
        
        # Source-target plots
        create_source_target_plots(deseq2_results, pathway_significant, pathway_net, PLOTS_DIR)
        
        # Leading edge plots
        create_leading_edge_plots(deseq2_results, pathway_significant, pathway_net, PLOTS_DIR)
        
        # Save results
        print("\n" + "="*60)
        print("SAVING RESULTS")
        print("="*60)
        
        save_enrichment_results(tf_significant, pathway_significant, hallmark_significant, TABLES_DIR)
        
        # Generate comprehensive report
        generate_comprehensive_report(tf_significant, pathway_significant, hallmark_significant, 
                                    deseq2_results, REPORTS_DIR)
        
        # Print final summary
        print("\n✅ PATHWAY ENRICHMENT ANALYSIS COMPLETED!")
        print("=" * 50)
        print(f"📊 Analyzed {len(deseq2_results)} contrasts")
        
        if tf_significant:
            total_sig_tfs = sum(data['n_significant'] for data in tf_significant.values())
            print(f"🧬 Found {total_sig_tfs} significant transcription factors")
        
        if pathway_significant:
            total_sig_pathways = sum(data['n_significant'] for data in pathway_significant.values())
            print(f"🛤️  Found {total_sig_pathways} significant pathways")
        
        if hallmark_significant:
            total_sig_hallmarks = sum(data['n_significant'] for data in hallmark_significant.values())
            print(f"🎯 Found {total_sig_hallmarks} significant hallmark gene sets")
        
        print(f"\n📁 Results saved to: {OUTPUT_DIR}")
        print("   📊 Plots: /plots/ (including tutorial-style plots)")
        print("   📋 Tables: /tables/") 
        print("   📄 Reports: /reports/")
        print("\n🎨 Plot types generated:")
        print("   • Basic volcano plots (tutorial style)")
        print("   • TF activity barplots")
        print("   • Pathway activity barplots & dotplots")
        print("   • Hallmark gene set barplots")
        print("   • Source-target scatter plots")
        print("   • Leading edge plots")
        print("   • Comparison heatmaps")
        
        # Offer interactive exploration
        print("\n" + "="*60)
        print("ANALYSIS COMPLETE - EXPLORATION OPTIONS")
        print("="*60)
        
        while True:
            explore_choice = input("\n🔍 Would you like to explore results interactively? (y/n/auto): ").strip().lower()
            
            if explore_choice in ['y', 'yes']:
                interactive_exploration(tf_significant, pathway_significant, hallmark_significant, deseq2_results)
                break
            elif explore_choice in ['n', 'no']:
                print("👍 Skipping interactive exploration")
                break
            elif explore_choice in ['auto', 'a']:
                print("\n🤖 Running automated exploration...")
                
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
                
                # Show most interesting contrast
                if tf_significant:
                    most_active_contrast = max(tf_significant.keys(), 
                                             key=lambda x: tf_significant[x]['n_significant'])
                    print(f"\n🔥 MOST ACTIVE CONTRAST: {most_active_contrast}")
                    explore_contrast_results(tf_significant, pathway_significant, hallmark_significant, 
                                           most_active_contrast)
                
                break
            else:
                print("Please enter 'y' for yes, 'n' for no, or 'auto' for automated exploration")
        
        return {
            'tf_results': tf_significant,
            'pathway_results': pathway_significant,
            'hallmark_results': hallmark_significant,
            'data_matrix': data_matrix,
            'deseq2_results': deseq2_results
        }
        
    except Exception as e:
        print(f"\n❌ ERROR: {str(e)}")
        raise e

if __name__ == "__main__":
    main()