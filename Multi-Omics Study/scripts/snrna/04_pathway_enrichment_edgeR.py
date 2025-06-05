"""
🔬 Pathway Enrichment Analysis - edgeR Results
GSE225158 - Comprehensive pathway enrichment from edgeR differential expression

Strategy:
1. Load edgeR results from multiple contrasts
2. Perform Gene Ontology (GO) enrichment analysis
3. Run KEGG pathway enrichment
4. Analyze Reactome pathways
5. Focus on complement system and neuroinflammation
6. Generate comprehensive visualizations

Dependencies:
- Python packages: pandas, numpy, matplotlib, seaborn, requests
- R packages: clusterProfiler, enrichplot, ReactomePA, org.Hs.eg.db
- rpy2 for Python-R integration
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import requests
import json
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# R integration for enrichment analysis
try:
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri, numpy2ri
    from rpy2.robjects.packages import importr
    R_AVAILABLE = True
    print("✅ R integration available")
except ImportError:
    R_AVAILABLE = False
    print("❌ R integration not available - install rpy2 and R packages")
    exit(1)

# Paths
BASE_DIR = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study"
EDGER_DIR = f"{BASE_DIR}/results/snrna_scvi/differential_expression_edgeR"
OUTPUT_DIR = f"{BASE_DIR}/results/snrna_scvi/pathway_enrichment_edgeR"

# Organized output directories
BASIC_DIR = f"{OUTPUT_DIR}/01_basic_enrichment"
ADVANCED_DIR = f"{OUTPUT_DIR}/02_advanced_analysis"
SPECIALIZED_DIR = f"{OUTPUT_DIR}/03_specialized_analysis"
VISUALIZATION_DIR = f"{OUTPUT_DIR}/04_visualizations"
COMPARISON_DIR = f"{OUTPUT_DIR}/05_comparative_analysis"
SUMMARY_DIR = f"{OUTPUT_DIR}/06_summary_reports"

# Sub-directories for specific analyses
PLOTS_DIR = f"{VISUALIZATION_DIR}/plots"
GSEA_DIR = f"{ADVANCED_DIR}/gsea_analysis"
MSIGDB_DIR = f"{ADVANCED_DIR}/msigdb_enrichment"
COMPLEMENT_DIR = f"{SPECIALIZED_DIR}/complement_system"
NEUROINFLAM_DIR = f"{SPECIALIZED_DIR}/neuroinflammation"
NETWORK_DIR = f"{VISUALIZATION_DIR}/network_maps"
UPSET_DIR = f"{VISUALIZATION_DIR}/upset_plots"
LEADING_EDGE_DIR = f"{ADVANCED_DIR}/leading_edge"
CROSSTALK_DIR = f"{COMPARISON_DIR}/pathway_crosstalk"
MODULES_DIR = f"{COMPARISON_DIR}/functional_modules"

# Create all output directories
output_dirs = [
    OUTPUT_DIR, BASIC_DIR, ADVANCED_DIR, SPECIALIZED_DIR, 
    VISUALIZATION_DIR, COMPARISON_DIR, SUMMARY_DIR, PLOTS_DIR,
    GSEA_DIR, MSIGDB_DIR, COMPLEMENT_DIR, NEUROINFLAM_DIR,
    NETWORK_DIR, UPSET_DIR, LEADING_EDGE_DIR, CROSSTALK_DIR, MODULES_DIR
]

for dir_path in output_dirs:
    os.makedirs(dir_path, exist_ok=True)

# Analysis parameters
PVALUE_THRESHOLD = 0.05
LOGFC_THRESHOLD = 0.5
MIN_GENES_PATHWAY = 5
MAX_GENES_PATHWAY = 500

# Use existing significance columns from edgeR results
# Priority: significant_nominal > significant_fdr10 > significant_fdr05 > significant_strict
SIGNIFICANCE_COLUMN_PRIORITY = ['significant_nominal', 'significant_fdr10', 'significant_fdr05', 'significant_strict', 'significant']

def load_edger_results():
    """Load all edgeR differential expression results"""
    print("\n📁 LOADING EDGER RESULTS")
    print("=" * 50)
    
    results_dict = {}
    
    # Find all edgeR result files
    result_files = list(Path(EDGER_DIR).glob("edgeR_*_results.csv"))
    
    if not result_files:
        raise FileNotFoundError(f"No edgeR results found in {EDGER_DIR}")
    
    for result_file in result_files:
        contrast_name = result_file.stem.replace("edgeR_", "").replace("_results", "")
        
        try:
            df = pd.read_csv(result_file)
            
            # Use the best available significance column from edgeR results
            sig_col = None
            for col in SIGNIFICANCE_COLUMN_PRIORITY:
                if col in df.columns:
                    sig_col = col
                    break
            
            if sig_col is None:
                # Fallback: create significance column
                df['pathway_significant'] = (df['FDR'] < PVALUE_THRESHOLD) & (abs(df['logFC']) > LOGFC_THRESHOLD)
                sig_col = 'pathway_significant'
            else:
                # Use existing column but rename for clarity
                df['pathway_significant'] = df[sig_col]
                sig_col = 'pathway_significant'
            
            results_dict[contrast_name] = df
            
            n_sig = df[sig_col].sum()
            sig_type = sig_col.replace('pathway_significant', f'created (FDR<{PVALUE_THRESHOLD}, |logFC|>{LOGFC_THRESHOLD})')
            if sig_col == 'pathway_significant':
                if 'significant_nominal' in df.columns:
                    sig_type = 'significant_nominal (p<0.05)'
                elif 'significant_fdr10' in df.columns:
                    sig_type = 'significant_fdr10 (FDR<0.10)'
                elif 'significant_fdr05' in df.columns:
                    sig_type = 'significant_fdr05 (FDR<0.05)'
            
            print(f"   📊 {contrast_name}: {len(df)} genes, {n_sig} significant ({sig_type})")
            
        except Exception as e:
            print(f"   ❌ Failed to load {result_file}: {e}")
    
    print(f"   ✅ Loaded {len(results_dict)} contrasts")
    return results_dict

def prepare_gene_lists(results_dict):
    """Prepare gene lists for pathway enrichment"""
    print("\n🧬 PREPARING GENE LISTS")
    print("=" * 50)
    
    gene_lists = {}
    
    for contrast_name, df in results_dict.items():
        print(f"   📋 Processing {contrast_name}...")
        
        # All significant genes (using the pathway_significant column)
        sig_genes = df[df['pathway_significant']]['gene'].tolist()
        
        # Upregulated and downregulated genes from significant set
        sig_df = df[df['pathway_significant']]
        up_genes = sig_df[sig_df['logFC'] > 0]['gene'].tolist()
        down_genes = sig_df[sig_df['logFC'] < 0]['gene'].tolist()
        
        # All genes ranked by signed log p-value (for GSEA)
        df_ranked = df.copy()
        df_ranked['signed_logp'] = -np.log10(df_ranked['PValue']) * np.sign(df_ranked['logFC'])
        df_ranked = df_ranked.sort_values('signed_logp', ascending=False)
        ranked_genes = dict(zip(df_ranked['gene'], df_ranked['signed_logp']))
        
        gene_lists[contrast_name] = {
            'all_significant': sig_genes,
            'upregulated': up_genes,
            'downregulated': down_genes,
            'ranked_all': ranked_genes
        }
        
        print(f"     • Significant: {len(sig_genes)}")
        print(f"     • Upregulated: {len(up_genes)}")
        print(f"     • Downregulated: {len(down_genes)}")
    
    return gene_lists

def run_go_enrichment(gene_lists):
    """Run Gene Ontology enrichment analysis using R/clusterProfiler"""
    print("\n🎯 GENE ONTOLOGY ENRICHMENT")
    print("=" * 50)
    
    if not R_AVAILABLE:
        raise RuntimeError("R integration not available")
    
    # Initialize R packages
    with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter + numpy2ri.converter):
        try:
            ro.r('library(clusterProfiler)')
            ro.r('library(org.Hs.eg.db)')
            ro.r('library(enrichplot)')
            print("   ✅ R packages loaded successfully")
        except Exception as e:
            raise RuntimeError(f"Could not load R packages: {e}")
        
        go_results = {}
        
        for contrast_name, genes in gene_lists.items():
            print(f"\n   🔬 Analyzing {contrast_name}...")
            
            # Skip if too few significant genes
            if len(genes['all_significant']) < MIN_GENES_PATHWAY:
                print(f"     ⚠️  Too few significant genes ({len(genes['all_significant'])}) - need at least {MIN_GENES_PATHWAY}")
                continue
            
            try:
                # Convert gene symbols to Entrez IDs
                gene_symbols = ro.StrVector(genes['all_significant'])
                ro.globalenv['gene_symbols'] = gene_symbols
                
                ro.r('''
                # Convert symbols to Entrez IDs
                gene_entrez <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", 
                                  OrgDb = org.Hs.eg.db, drop = TRUE)
                entrez_genes <- gene_entrez$ENTREZID
                
                # Background genes (all detected genes)
                universe_symbols <- unique(c(gene_symbols))  # Placeholder - should be all detected genes
                universe_entrez <- bitr(universe_symbols, fromType = "SYMBOL", toType = "ENTREZID",
                                      OrgDb = org.Hs.eg.db, drop = TRUE)$ENTREZID
                ''')
                
                # GO Biological Process
                ro.r('''
                go_bp <- enrichGO(gene = entrez_genes,
                                universe = universe_entrez,
                                OrgDb = org.Hs.eg.db,
                                ont = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.2,
                                readable = TRUE,
                                minGSSize = 5,
                                maxGSSize = 500)
                ''')
                
                # GO Molecular Function
                ro.r('''
                go_mf <- enrichGO(gene = entrez_genes,
                                universe = universe_entrez,
                                OrgDb = org.Hs.eg.db,
                                ont = "MF",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.2,
                                readable = TRUE,
                                minGSSize = 5,
                                maxGSSize = 500)
                ''')
                
                # GO Cellular Component
                ro.r('''
                go_cc <- enrichGO(gene = entrez_genes,
                                universe = universe_entrez,
                                OrgDb = org.Hs.eg.db,
                                ont = "CC",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.2,
                                readable = TRUE,
                                minGSSize = 5,
                                maxGSSize = 500)
                ''')
                
                # Extract results
                go_bp_df = extract_enrichment_results('go_bp', 'GO_BP')
                go_mf_df = extract_enrichment_results('go_mf', 'GO_MF')
                go_cc_df = extract_enrichment_results('go_cc', 'GO_CC')
                
                go_results[contrast_name] = {
                    'GO_BP': go_bp_df,
                    'GO_MF': go_mf_df,
                    'GO_CC': go_cc_df
                }
                
                print(f"     ✅ GO BP: {len(go_bp_df)} terms")
                print(f"     ✅ GO MF: {len(go_mf_df)} terms")
                print(f"     ✅ GO CC: {len(go_cc_df)} terms")
                
            except Exception as e:
                print(f"     ❌ Failed GO enrichment for {contrast_name}: {e}")
        
        return go_results

def run_kegg_enrichment(gene_lists):
    """Run KEGG pathway enrichment analysis"""
    print("\n🛤️  KEGG PATHWAY ENRICHMENT")
    print("=" * 50)
    
    with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter + numpy2ri.converter):
        kegg_results = {}
        
        for contrast_name, genes in gene_lists.items():
            if len(genes['all_significant']) < MIN_GENES_PATHWAY:
                print(f"     ⚠️  Skipping {contrast_name}: too few genes ({len(genes['all_significant'])})")
                continue
                
            print(f"\n   🔬 Analyzing {contrast_name}...")
            
            try:
                gene_symbols = ro.StrVector(genes['all_significant'])
                ro.globalenv['gene_symbols'] = gene_symbols
                
                ro.r('''
                # Convert to Entrez IDs
                gene_entrez <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", 
                                  OrgDb = org.Hs.eg.db, drop = TRUE)
                entrez_genes <- gene_entrez$ENTREZID
                
                # KEGG enrichment
                kegg_result <- enrichKEGG(gene = entrez_genes,
                                        organism = "hsa",
                                        pvalueCutoff = 0.05,
                                        qvalueCutoff = 0.2,
                                        minGSSize = 5,
                                        maxGSSize = 500)
                ''')
                
                kegg_df = extract_enrichment_results('kegg_result', 'KEGG')
                kegg_results[contrast_name] = kegg_df
                
                print(f"     ✅ KEGG: {len(kegg_df)} pathways")
                
            except Exception as e:
                print(f"     ❌ Failed KEGG enrichment for {contrast_name}: {e}")
        
        return kegg_results

def run_reactome_enrichment(gene_lists):
    """Run Reactome pathway enrichment analysis"""
    print("\n⚛️  REACTOME PATHWAY ENRICHMENT")
    print("=" * 50)
    
    with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter + numpy2ri.converter):
        try:
            ro.r('library(ReactomePA)')
            print("   ✅ ReactomePA loaded")
        except Exception as e:
            print(f"   ❌ ReactomePA not available: {e}")
            return {}
        
        reactome_results = {}
        
        for contrast_name, genes in gene_lists.items():
            if len(genes['all_significant']) < MIN_GENES_PATHWAY:
                print(f"     ⚠️  Skipping {contrast_name}: too few genes ({len(genes['all_significant'])})")
                continue
                
            print(f"\n   🔬 Analyzing {contrast_name}...")
            
            try:
                gene_symbols = ro.StrVector(genes['all_significant'])
                ro.globalenv['gene_symbols'] = gene_symbols
                
                ro.r('''
                # Convert to Entrez IDs
                gene_entrez <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", 
                                  OrgDb = org.Hs.eg.db, drop = TRUE)
                entrez_genes <- gene_entrez$ENTREZID
                
                # Reactome enrichment
                reactome_result <- enrichPathway(gene = entrez_genes,
                                               organism = "human",
                                               pvalueCutoff = 0.05,
                                               qvalueCutoff = 0.2,
                                               minGSSize = 5,
                                               maxGSSize = 500,
                                               readable = TRUE)
                ''')
                
                reactome_df = extract_enrichment_results('reactome_result', 'Reactome')
                reactome_results[contrast_name] = reactome_df
                
                print(f"     ✅ Reactome: {len(reactome_df)} pathways")
                
            except Exception as e:
                print(f"     ❌ Failed Reactome enrichment for {contrast_name}: {e}")
        
        return reactome_results

def extract_enrichment_results(r_object_name, analysis_type):
    """Extract enrichment results from R object to pandas DataFrame"""
    try:
        ro.r(f'''
        if(nrow({r_object_name}@result) > 0) {{
            result_df <- {r_object_name}@result
            result_ids <- result_df$ID
            result_descriptions <- result_df$Description
            result_pvalues <- result_df$pvalue
            result_qvalues <- result_df$qvalue
            result_counts <- result_df$Count
            result_genes <- result_df$geneID
        }} else {{
            result_ids <- character(0)
            result_descriptions <- character(0)
            result_pvalues <- numeric(0)
            result_qvalues <- numeric(0)
            result_counts <- numeric(0)
            result_genes <- character(0)
        }}
        ''')
        
        # Extract components
        ids = list(ro.r('result_ids'))
        descriptions = list(ro.r('result_descriptions'))
        pvalues = list(ro.r('result_pvalues'))
        qvalues = list(ro.r('result_qvalues'))
        counts = list(ro.r('result_counts'))
        genes = list(ro.r('result_genes'))
        
        if len(ids) == 0:
            return pd.DataFrame()
        
        df = pd.DataFrame({
            'ID': ids,
            'Description': descriptions,
            'pvalue': pvalues,
            'qvalue': qvalues,
            'Count': counts,
            'geneID': genes,
            'analysis_type': analysis_type
        })
        
        return df
        
    except Exception as e:
        print(f"   ❌ Failed to extract {analysis_type} results: {e}")
        return pd.DataFrame()

def analyze_complement_pathways(enrichment_results):
    """Focus analysis on complement system pathways"""
    print("\n🧬 COMPLEMENT SYSTEM ANALYSIS")
    print("=" * 50)
    
    # Complement-related keywords
    complement_keywords = [
        'complement', 'classical complement', 'alternative complement', 'lectin complement',
        'membrane attack complex', 'MAC', 'C1q', 'C3', 'C4', 'C5', 'C9',
        'complement activation', 'complement cascade', 'complement regulation',
        'complement receptor', 'complement inhibitor', 'complement factor',
        'innate immune', 'humoral immune', 'immune complex'
    ]
    
    complement_results = {}
    
    for contrast_name, results in enrichment_results.items():
        complement_pathways = []
        
        for analysis_type, df in results.items():
            if df.empty:
                continue
                
            # Find complement-related pathways
            mask = df['Description'].str.lower().str.contains('|'.join(complement_keywords), na=False)
            complement_df = df[mask].copy()
            
            if not complement_df.empty:
                complement_df['contrast'] = contrast_name
                complement_pathways.append(complement_df)
        
        if complement_pathways:
            complement_results[contrast_name] = pd.concat(complement_pathways, ignore_index=True)
            n_pathways = len(complement_results[contrast_name])
            print(f"   🎯 {contrast_name}: {n_pathways} complement pathways")
    
    return complement_results

def create_enrichment_plots(enrichment_results):
    """Create comprehensive visualization plots"""
    print("\n📊 CREATING ENRICHMENT PLOTS")
    print("=" * 50)
    
    # 1. Top pathways heatmap
    create_pathway_heatmap(enrichment_results)
    
    # 2. Dot plots for each contrast
    create_dotplots(enrichment_results)
    
    # 3. Network plots
    create_network_plots(enrichment_results)
    
    # 4. Comparison plots
    create_comparison_plots(enrichment_results)
    
    print("   ✅ All plots created successfully")

def create_pathway_heatmap(enrichment_results):
    """Create heatmap of top enriched pathways across contrasts"""
    print("   📊 Creating pathway heatmap...")
    
    # Collect top pathways from all contrasts
    all_pathways = []
    
    for contrast_name, results in enrichment_results.items():
        for analysis_type, df in results.items():
            if df.empty:
                continue
            
            # Get top 10 pathways
            top_pathways = df.nsmallest(10, 'pvalue').copy()
            top_pathways['contrast'] = contrast_name
            top_pathways['analysis_type'] = analysis_type
            all_pathways.append(top_pathways)
    
    if not all_pathways:
        print("     ⚠️  No pathways to plot")
        return
    
    combined_df = pd.concat(all_pathways, ignore_index=True)
    
    # Create pivot table for heatmap
    heatmap_data = combined_df.pivot_table(
        index='Description',
        columns='contrast',
        values='pvalue',
        aggfunc='min'
    )
    
    # Transform p-values for better visualization
    heatmap_data = -np.log10(heatmap_data)
    
    # Plot heatmap
    plt.figure(figsize=(15, 10))
    sns.heatmap(
        heatmap_data,
        cmap='Reds',
        annot=False,
        cbar_kws={'label': '-log10(p-value)'},
        xticklabels=True,
        yticklabels=True
    )
    plt.title('Top Enriched Pathways Across Contrasts', fontsize=16, fontweight='bold')
    plt.xlabel('Contrast', fontsize=12)
    plt.ylabel('Pathway', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0, fontsize=8)
    plt.tight_layout()
    plt.savefig(f"{PLOTS_DIR}/pathway_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()

def create_dotplots(enrichment_results):
    """Create dot plots for each contrast"""
    print("   📊 Creating dot plots...")
    
    for contrast_name, results in enrichment_results.items():
        fig, axes = plt.subplots(1, 3, figsize=(20, 8))
        fig.suptitle(f'Pathway Enrichment: {contrast_name}', fontsize=16, fontweight='bold')
        
        analysis_types = ['GO_BP', 'KEGG', 'Reactome']
        
        for i, analysis_type in enumerate(analysis_types):
            ax = axes[i]
            
            if analysis_type in results and not results[analysis_type].empty:
                df = results[analysis_type].copy()
                
                # Get top 15 pathways
                df_plot = df.nsmallest(15, 'pvalue')
                
                if len(df_plot) > 0:
                    # Create dot plot
                    scatter = ax.scatter(
                        df_plot['Count'],
                        range(len(df_plot)),
                        c=-np.log10(df_plot['pvalue']),
                        s=100,
                        cmap='Reds',
                        alpha=0.7
                    )
                    
                    ax.set_yticks(range(len(df_plot)))
                    ax.set_yticklabels([desc[:50] + '...' if len(desc) > 50 else desc 
                                       for desc in df_plot['Description']], fontsize=8)
                    ax.set_xlabel('Gene Count', fontsize=10)
                    ax.set_title(f'{analysis_type}', fontsize=12, fontweight='bold')
                    ax.grid(True, alpha=0.3)
                    
                    # Add colorbar
                    cbar = plt.colorbar(scatter, ax=ax)
                    cbar.set_label('-log10(p-value)', fontsize=8)
                else:
                    ax.text(0.5, 0.5, 'No significant pathways', 
                           ha='center', va='center', transform=ax.transAxes)
                    ax.set_title(f'{analysis_type}', fontsize=12, fontweight='bold')
            else:
                ax.text(0.5, 0.5, 'No data available', 
                       ha='center', va='center', transform=ax.transAxes)
                ax.set_title(f'{analysis_type}', fontsize=12, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(f"{PLOTS_DIR}/dotplot_{contrast_name}.png", dpi=300, bbox_inches='tight')
        plt.close()

def create_network_plots(enrichment_results):
    """Create network plots showing pathway relationships"""
    print("   📊 Creating network plots...")
    
    # This is a placeholder for network analysis
    # Would require additional packages like networkx for full implementation
    
    for contrast_name, results in enrichment_results.items():
        # Simple pathway count plot as placeholder
        counts = []
        labels = []
        
        for analysis_type, df in results.items():
            if not df.empty:
                counts.append(len(df))
                labels.append(analysis_type)
        
        if counts:
            plt.figure(figsize=(8, 6))
            bars = plt.bar(labels, counts, color=['skyblue', 'lightcoral', 'lightgreen'])
            plt.title(f'Enriched Pathway Counts: {contrast_name}', fontsize=14, fontweight='bold')
            plt.xlabel('Analysis Type', fontsize=12)
            plt.ylabel('Number of Pathways', fontsize=12)
            
            # Add value labels on bars
            for bar, count in zip(bars, counts):
                plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                        str(count), ha='center', va='bottom', fontweight='bold')
            
            plt.tight_layout()
            plt.savefig(f"{PLOTS_DIR}/pathway_counts_{contrast_name}.png", dpi=300, bbox_inches='tight')
            plt.close()

def create_comparison_plots(enrichment_results):
    """Create comparison plots across contrasts"""
    print("   📊 Creating comparison plots...")
    
    # Pathway overlap analysis
    all_pathways = {}
    
    for contrast_name, results in enrichment_results.items():
        pathways = set()
        for analysis_type, df in results.items():
            if not df.empty:
                pathways.update(df['Description'].tolist())
        all_pathways[contrast_name] = pathways
    
    # Create overlap matrix
    contrasts = list(all_pathways.keys())
    overlap_matrix = pd.DataFrame(0, index=contrasts, columns=contrasts)
    
    for i, contrast1 in enumerate(contrasts):
        for j, contrast2 in enumerate(contrasts):
            if i <= j:
                overlap = len(all_pathways[contrast1] & all_pathways[contrast2])
                overlap_matrix.loc[contrast1, contrast2] = overlap
                overlap_matrix.loc[contrast2, contrast1] = overlap
    
    # Plot overlap heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(
        overlap_matrix,
        annot=True,
        fmt='d',
        cmap='Blues',
        square=True,
        cbar_kws={'label': 'Shared Pathways'}
    )
    plt.title('Pathway Overlap Between Contrasts', fontsize=14, fontweight='bold')
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(f"{PLOTS_DIR}/pathway_overlap.png", dpi=300, bbox_inches='tight')
    plt.close()

def save_enrichment_results(all_results, complement_results):
    """Save all enrichment analysis results"""
    print("\n💾 SAVING ENRICHMENT RESULTS")
    print("=" * 50)
    
    # Save individual contrast results
    for contrast_name, results in all_results.items():
        for analysis_type, df in results.items():
            if not df.empty:
                output_file = f"{OUTPUT_DIR}/{analysis_type}_{contrast_name}_enrichment.csv"
                df.to_csv(output_file, index=False)
                print(f"   📁 {analysis_type} - {contrast_name}: {output_file}")
    
    # Save combined results
    all_combined = []
    for contrast_name, results in all_results.items():
        for analysis_type, df in results.items():
            if not df.empty:
                df_copy = df.copy()
                df_copy['contrast'] = contrast_name
                df_copy['analysis_type'] = analysis_type
                all_combined.append(df_copy)
    
    if all_combined:
        combined_df = pd.concat(all_combined, ignore_index=True)
        combined_file = f"{OUTPUT_DIR}/all_pathway_enrichment_results.csv"
        combined_df.to_csv(combined_file, index=False)
        print(f"   📁 Combined results: {combined_file}")
    
    # Save complement-focused results
    if complement_results:
        complement_combined = pd.concat(complement_results.values(), ignore_index=True)
        complement_file = f"{OUTPUT_DIR}/complement_pathway_enrichment.csv"
        complement_combined.to_csv(complement_file, index=False)
        print(f"   📁 Complement pathways: {complement_file}")
    
    # Create summary report
    create_summary_report(all_results, complement_results)
    
    print("   ✅ All results saved successfully")

def create_summary_report(all_results, complement_results):
    """Create comprehensive summary report"""
    summary_lines = []
    summary_lines.append("PATHWAY ENRICHMENT ANALYSIS SUMMARY")
    summary_lines.append("=" * 50)
    summary_lines.append("")
    
    # Overall statistics
    total_pathways = 0
    total_contrasts = len(all_results)
    
    for contrast_name, results in all_results.items():
        contrast_pathways = 0
        summary_lines.append(f"CONTRAST: {contrast_name}")
        summary_lines.append("-" * 30)
        
        for analysis_type, df in results.items():
            n_pathways = len(df)
            contrast_pathways += n_pathways
            total_pathways += n_pathways
            summary_lines.append(f"  {analysis_type}: {n_pathways} pathways")
        
        summary_lines.append(f"  Total: {contrast_pathways} pathways")
        summary_lines.append("")
    
    summary_lines.append(f"OVERALL SUMMARY:")
    summary_lines.append(f"  Total contrasts analyzed: {total_contrasts}")
    summary_lines.append(f"  Total pathways enriched: {total_pathways}")
    
    if complement_results:
        total_complement = sum(len(df) for df in complement_results.values())
        summary_lines.append(f"  Complement pathways: {total_complement}")
    
    summary_lines.append("")
    summary_lines.append(f"Analysis completed at: {pd.Timestamp.now()}")
    
    # Save summary
    summary_file = f"{OUTPUT_DIR}/pathway_enrichment_summary.txt"
    with open(summary_file, 'w') as f:
        f.write('\n'.join(summary_lines))
    
    print(f"   📁 Summary report: {summary_file}")

def main():
    """Main pathway enrichment workflow"""
    print("🛤️  PATHWAY ENRICHMENT ANALYSIS - EDGER RESULTS")
    print("=" * 70)
    print("Comprehensive pathway enrichment from edgeR differential expression")
    print("=" * 70)
    
    try:
        # Load edgeR results
        results_dict = load_edger_results()
        
        # Prepare gene lists
        gene_lists = prepare_gene_lists(results_dict)
        
        # Run enrichment analyses
        print("\n🔬 RUNNING ENRICHMENT ANALYSES")
        print("=" * 50)
        
        go_results = run_go_enrichment(gene_lists)
        kegg_results = run_kegg_enrichment(gene_lists)
        reactome_results = run_reactome_enrichment(gene_lists)
        
        # Combine all results
        all_results = {}
        for contrast_name in gene_lists.keys():
            all_results[contrast_name] = {}
            
            if contrast_name in go_results:
                all_results[contrast_name].update(go_results[contrast_name])
            
            if contrast_name in kegg_results:
                all_results[contrast_name]['KEGG'] = kegg_results[contrast_name]
            
            if contrast_name in reactome_results:
                all_results[contrast_name]['Reactome'] = reactome_results[contrast_name]
        
        # Analyze complement pathways
        complement_results = analyze_complement_pathways(all_results)
        
        # Create visualizations
        create_enrichment_plots(all_results)
        
        # Save results
        save_enrichment_results(all_results, complement_results)
        
        print("\n" + "=" * 70)
        print("✅ PATHWAY ENRICHMENT ANALYSIS COMPLETE!")
        
        total_pathways = sum(
            len(df) for results in all_results.values() 
            for df in results.values() if not df.empty
        )
        
        print(f"📊 Total enriched pathways: {total_pathways}")
        print(f"🧬 Complement pathways: {sum(len(df) for df in complement_results.values())}")
        print(f"📁 Results saved to: {OUTPUT_DIR}/")
        print("🎯 Ready for biological interpretation!")
        print("=" * 70)
        
        return all_results, complement_results
        
    except Exception as e:
        print(f"❌ Pathway enrichment analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return None, None

if __name__ == "__main__":
    main()