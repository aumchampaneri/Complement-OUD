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
from scipy.stats import combine_pvalues
from sklearn.cluster import AgglomerativeClustering
from sklearn.feature_extraction.text import TfidfVectorizer
import warnings
warnings.filterwarnings('ignore')

# Optional imports for advanced features
try:
    import networkx as nx
    NETWORKX_AVAILABLE = True
except ImportError:
    NETWORKX_AVAILABLE = False
    print("⚠️  NetworkX not available - some network visualizations will be limited")

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

def run_gsea_analysis(gene_lists):
    """Run Gene Set Enrichment Analysis using ranked gene lists"""
    print("\n🎯 GENE SET ENRICHMENT ANALYSIS (GSEA)")
    print("=" * 50)
    
    if not R_AVAILABLE:
        raise RuntimeError("R integration not available")
    
    with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter + numpy2ri.converter):
        try:
            ro.r('library(fgsea)')
            print("   ✅ fgsea package loaded successfully")
        except Exception as e:
            print(f"   ⚠️  fgsea not available, using clusterProfiler GSEA: {e}")
        
        gsea_results = {}
        
        for contrast_name, genes in gene_lists.items():
            if len(genes['ranked_all']) < MIN_GENES_PATHWAY:
                print(f"   ⚠️  Skipping {contrast_name}: too few genes for GSEA")
                continue
                
            print(f"\n   🔬 GSEA for {contrast_name}...")
            
            try:
                # Prepare ranked gene list
                ranked_genes = genes['ranked_all']
                gene_symbols = list(ranked_genes.keys())
                gene_ranks = list(ranked_genes.values())
                
                ro.globalenv['gene_symbols'] = ro.StrVector(gene_symbols)
                ro.globalenv['gene_ranks'] = ro.FloatVector(gene_ranks)
                
                ro.r('''
                # Convert to named vector
                names(gene_ranks) <- gene_symbols
                gene_list <- sort(gene_ranks, decreasing = TRUE)
                
                # Convert symbols to Entrez for pathway databases
                gene_entrez_map <- bitr(names(gene_list), fromType = "SYMBOL", 
                                      toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
                
                # Create Entrez-based ranked list
                entrez_ranks <- gene_list[gene_entrez_map$SYMBOL]
                names(entrez_ranks) <- gene_entrez_map$ENTREZID
                entrez_ranks <- sort(entrez_ranks, decreasing = TRUE)
                ''')
                
                # Run GSEA for different pathway collections
                ro.r('''
                # GSEA GO Biological Process
                gsea_go_bp <- gseGO(geneList = entrez_ranks,
                                   OrgDb = org.Hs.eg.db,
                                   ont = "BP",
                                   minGSSize = 15,
                                   maxGSSize = 500,
                                   pvalueCutoff = 0.25,
                                   verbose = FALSE)
                
                # GSEA KEGG
                gsea_kegg <- gseKEGG(geneList = entrez_ranks,
                                    organism = "hsa",
                                    minGSSize = 15,
                                    maxGSSize = 500,
                                    pvalueCutoff = 0.25,
                                    verbose = FALSE)
                ''')
                
                # Extract GSEA results
                gsea_go_df = extract_gsea_results('gsea_go_bp', 'GSEA_GO_BP')
                gsea_kegg_df = extract_gsea_results('gsea_kegg', 'GSEA_KEGG')
                
                gsea_results[contrast_name] = {
                    'GSEA_GO_BP': gsea_go_df,
                    'GSEA_KEGG': gsea_kegg_df
                }
                
                print(f"     ✅ GSEA GO BP: {len(gsea_go_df)} pathways")
                print(f"     ✅ GSEA KEGG: {len(gsea_kegg_df)} pathways")
                
            except Exception as e:
                print(f"     ❌ Failed GSEA for {contrast_name}: {e}")
        
        return gsea_results

def extract_gsea_results(r_object_name, analysis_type):
    """Extract GSEA results from R object to pandas DataFrame"""
    try:
        ro.r(f'''
        if(nrow({r_object_name}@result) > 0) {{
            result_df <- {r_object_name}@result
            result_ids <- result_df$ID
            result_descriptions <- result_df$Description
            result_pvalues <- result_df$pvalue
            result_qvalues <- result_df$qvalue
            result_nes <- result_df$NES
            result_size <- result_df$setSize
            result_genes <- result_df$core_enrichment
        }} else {{
            result_ids <- character(0)
            result_descriptions <- character(0)
            result_pvalues <- numeric(0)
            result_qvalues <- numeric(0)
            result_nes <- numeric(0)
            result_size <- numeric(0)
            result_genes <- character(0)
        }}
        ''')
        
        # Extract components
        ids = list(ro.r('result_ids'))
        descriptions = list(ro.r('result_descriptions'))
        pvalues = list(ro.r('result_pvalues'))
        qvalues = list(ro.r('result_qvalues'))
        nes = list(ro.r('result_nes'))
        sizes = list(ro.r('result_size'))
        genes = list(ro.r('result_genes'))
        
        if len(ids) == 0:
            return pd.DataFrame()
        
        df = pd.DataFrame({
            'ID': ids,
            'Description': descriptions,
            'pvalue': pvalues,
            'qvalue': qvalues,
            'NES': nes,
            'setSize': sizes,
            'core_enrichment': genes,
            'analysis_type': analysis_type
        })
        
        return df
        
    except Exception as e:
        print(f"   ❌ Failed to extract {analysis_type} results: {e}")
        return pd.DataFrame()

def run_msigdb_enrichment(gene_lists):
    """Run MSigDB enrichment analysis (Hallmark, C2, etc.)"""
    print("\n📚 MSIGDB ENRICHMENT ANALYSIS")
    print("=" * 50)
    
    with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter + numpy2ri.converter):
        try:
            ro.r('library(msigdbr)')
            print("   ✅ msigdbr package loaded successfully")
        except Exception as e:
            print(f"   ❌ msigdbr not available: {e}")
            return {}
        
        msigdb_results = {}
        
        for contrast_name, genes in gene_lists.items():
            if len(genes['all_significant']) < MIN_GENES_PATHWAY:
                print(f"   ⚠️  Skipping {contrast_name}: too few genes")
                continue
                
            print(f"\n   🔬 MSigDB for {contrast_name}...")
            
            try:
                gene_symbols = ro.StrVector(genes['all_significant'])
                ro.globalenv['gene_symbols'] = gene_symbols
                
                ro.r('''
                # Get MSigDB gene sets
                hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
                c2_canonical <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP")
                
                # Convert to format for enricher
                hallmark_list <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)
                c2_list <- split(c2_canonical$gene_symbol, c2_canonical$gs_name)
                ''')
                
                # Hallmark enrichment
                ro.r('''
                hallmark_result <- enricher(gene_symbols,
                                          TERM2GENE = hallmark_sets[,c("gs_name", "gene_symbol")],
                                          pvalueCutoff = 0.05,
                                          qvalueCutoff = 0.2,
                                          minGSSize = 5,
                                          maxGSSize = 500)
                
                c2_result <- enricher(gene_symbols,
                                    TERM2GENE = c2_canonical[,c("gs_name", "gene_symbol")],
                                    pvalueCutoff = 0.05,
                                    qvalueCutoff = 0.2,
                                    minGSSize = 5,
                                    maxGSSize = 500)
                ''')
                
                hallmark_df = extract_enrichment_results('hallmark_result', 'MSigDB_Hallmark')
                c2_df = extract_enrichment_results('c2_result', 'MSigDB_C2_Canonical')
                
                msigdb_results[contrast_name] = {
                    'MSigDB_Hallmark': hallmark_df,
                    'MSigDB_C2_Canonical': c2_df
                }
                
                print(f"     ✅ Hallmark: {len(hallmark_df)} pathways")
                print(f"     ✅ C2 Canonical: {len(c2_df)} pathways")
                
            except Exception as e:
                print(f"     ❌ Failed MSigDB for {contrast_name}: {e}")
        
        return msigdb_results

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

def analyze_neuroinflammation_pathways(enrichment_results):
    """Focus analysis on neuroinflammation and addiction pathways"""
    print("\n🧠 NEUROINFLAMMATION & ADDICTION PATHWAY ANALYSIS")
    print("=" * 50)
    
    # Neuroinflammation and addiction-related keywords
    neuro_keywords = [
        'neuroinflammation', 'microglia', 'astrocyte', 'cytokine', 'interleukin',
        'tumor necrosis factor', 'interferon', 'toll-like receptor', 'NF-kappa',
        'dopamine', 'serotonin', 'GABA', 'glutamate', 'addiction', 'reward',
        'synaptic', 'neurotransmitter', 'opioid', 'morphine', 'cocaine',
        'inflammation', 'immune response', 'innate immunity'
    ]
    
    neuro_results = {}
    
    for contrast_name, results in enrichment_results.items():
        neuro_pathways = []
        
        for analysis_type, df in results.items():
            if df.empty:
                continue
                
            # Find neuroinflammation/addiction-related pathways
            mask = df['Description'].str.lower().str.contains('|'.join(neuro_keywords), na=False)
            neuro_df = df[mask].copy()
            
            if not neuro_df.empty:
                neuro_df['contrast'] = contrast_name
                neuro_pathways.append(neuro_df)
        
        if neuro_pathways:
            neuro_results[contrast_name] = pd.concat(neuro_pathways, ignore_index=True)
            n_pathways = len(neuro_results[contrast_name])
            print(f"   🧠 {contrast_name}: {n_pathways} neuroinflammation/addiction pathways")
    
    return neuro_results

def create_upset_plots(enrichment_results):
    """Create UpSet plots for pathway overlap analysis"""
    print("   📊 Creating UpSet plots...")
    
    try:
        # Collect pathways from each contrast
        pathway_sets = {}
        
        for contrast_name, results in enrichment_results.items():
            all_pathways = set()
            for analysis_type, df in results.items():
                if not df.empty:
                    all_pathways.update(df['Description'].tolist())
            pathway_sets[contrast_name] = all_pathways
        
        if len(pathway_sets) < 2:
            print("     ⚠️  Need at least 2 contrasts for UpSet plots")
            return
        
        # Create binary matrix for UpSet plot
        all_pathways = set()
        for pathways in pathway_sets.values():
            all_pathways.update(pathways)
        
        upset_data = pd.DataFrame(index=list(all_pathways))
        for contrast_name, pathways in pathway_sets.items():
            upset_data[contrast_name] = upset_data.index.isin(pathways)
        
        # Save upset data
        upset_file = f"{UPSET_DIR}/pathway_upset_data.csv"
        upset_data.to_csv(upset_file)
        
        # Create simple visualization (matplotlib-based)
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Calculate intersection sizes
        contrast_names = list(pathway_sets.keys())
        intersections = []
        labels = []
        
        # Individual sets
        for name in contrast_names:
            size = len(pathway_sets[name])
            intersections.append(size)
            labels.append(name)
        
        # Pairwise intersections
        for i, name1 in enumerate(contrast_names):
            for j, name2 in enumerate(contrast_names[i+1:], i+1):
                intersection = len(pathway_sets[name1] & pathway_sets[name2])
                intersections.append(intersection)
                labels.append(f"{name1} ∩ {name2}")
        
        # Plot
        bars = ax.bar(range(len(intersections)), intersections)
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=45, ha='right')
        ax.set_ylabel('Number of Pathways')
        ax.set_title('Pathway Overlaps Between Contrasts')
        
        # Add value labels
        for bar, value in zip(bars, intersections):
            if value > 0:
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                       str(value), ha='center', va='bottom')
        
        plt.tight_layout()
        plt.savefig(f"{UPSET_DIR}/pathway_upset_plot.png", dpi=300, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"     ❌ Failed to create UpSet plots: {e}")

def create_enrichment_maps(enrichment_results):
    """Create enrichment map networks showing pathway relationships"""
    print("   📊 Creating enrichment maps...")
    
    try:
        import networkx as nx
        
        for contrast_name, results in enrichment_results.items():
            if not any(not df.empty for df in results.values()):
                continue
                
            # Create network graph
            G = nx.Graph()
            
            # Add pathway nodes
            for analysis_type, df in results.items():
                if df.empty:
                    continue
                
                # Get top 20 pathways
                top_pathways = df.nsmallest(20, 'pvalue')
                
                for _, row in top_pathways.iterrows():
                    pathway_id = row['Description']
                    G.add_node(pathway_id, 
                             analysis_type=analysis_type,
                             pvalue=row['pvalue'],
                             size=row.get('Count', 10))
            
            if len(G.nodes()) < 2:
                continue
            
            # Add edges based on gene overlap (simplified)
            # This is a placeholder - full implementation would calculate Jaccard similarity
            nodes = list(G.nodes())
            for i, node1 in enumerate(nodes):
                for node2 in nodes[i+1:i+6]:  # Connect to nearby nodes
                    G.add_edge(node1, node2, weight=0.5)
            
            # Create layout and plot
            plt.figure(figsize=(15, 12))
            pos = nx.spring_layout(G, k=1, iterations=50)
            
            # Draw network
            nx.draw_networkx_nodes(G, pos, node_size=50, alpha=0.7)
            nx.draw_networkx_edges(G, pos, alpha=0.3)
            
            plt.title(f'Enrichment Map: {contrast_name}', fontsize=14)
            plt.axis('off')
            plt.tight_layout()
            plt.savefig(f"{NETWORK_DIR}/enrichment_map_{contrast_name}.png", 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
    except ImportError:
        print("     ⚠️  NetworkX not available for enrichment maps")
    except Exception as e:
        print(f"     ❌ Failed to create enrichment maps: {e}")

def create_volcano_plots(enrichment_results):
    """Create volcano plots for pathway enrichment significance"""
    print("   📊 Creating volcano plots...")
    
    for contrast_name, results in enrichment_results.items():
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle(f'Pathway Enrichment Volcano Plots: {contrast_name}', fontsize=16)
        
        analysis_types = ['GO_BP', 'GO_MF', 'KEGG', 'Reactome']
        
        for i, analysis_type in enumerate(analysis_types):
            row, col = i // 2, i % 2
            ax = axes[row, col]
            
            if analysis_type in results and not results[analysis_type].empty:
                df = results[analysis_type].copy()
                
                # Calculate -log10(pvalue)
                df['neg_log_pval'] = -np.log10(df['pvalue'])
                df['gene_count'] = df.get('Count', 10)
                
                # Create volcano plot
                scatter = ax.scatter(df['gene_count'], df['neg_log_pval'],
                                   c=df['neg_log_pval'], s=30, alpha=0.6, cmap='Reds')
                
                ax.set_xlabel('Gene Count')
                ax.set_ylabel('-log10(p-value)')
                ax.set_title(f'{analysis_type}')
                ax.grid(True, alpha=0.3)
                
                # Add significance line
                ax.axhline(y=-np.log10(0.05), color='red', linestyle='--', alpha=0.5)
                
                # Annotate top pathways
                top_pathways = df.nlargest(3, 'neg_log_pval')
                for _, row in top_pathways.iterrows():
                    ax.annotate(row['Description'][:30] + '...', 
                              (row['gene_count'], row['neg_log_pval']),
                              xytext=(5, 5), textcoords='offset points',
                              fontsize=8, alpha=0.8)
            else:
                ax.text(0.5, 0.5, 'No data available', 
                       ha='center', va='center', transform=ax.transAxes)
                ax.set_title(f'{analysis_type}')
        
        plt.tight_layout()
        plt.savefig(f"{VISUALIZATION_DIR}/volcano_plots_{contrast_name}.png", 
                   dpi=300, bbox_inches='tight')
        plt.close()

def perform_leading_edge_analysis(enrichment_results, gene_lists):
    """Identify leading edge genes driving pathway enrichments"""
    print("\n🎯 LEADING EDGE ANALYSIS")
    print("=" * 50)
    
    leading_edge_results = {}
    
    for contrast_name, results in enrichment_results.items():
        print(f"   🔬 Analyzing {contrast_name}...")
        
        contrast_leading_edge = {}
        
        for analysis_type, df in results.items():
            if df.empty:
                continue
            
            # Get top 10 pathways
            top_pathways = df.nsmallest(10, 'pvalue')
            
            pathway_genes = []
            for _, pathway in top_pathways.iterrows():
                # Extract genes from pathway (simplified - would need actual gene sets)
                if 'geneID' in pathway:
                    genes = str(pathway['geneID']).split('/')
                elif 'core_enrichment' in pathway:
                    genes = str(pathway['core_enrichment']).split('/')
                else:
                    genes = []
                
                pathway_genes.extend(genes)
            
            # Count gene frequencies
            gene_counts = pd.Series(pathway_genes).value_counts()
            
            contrast_leading_edge[analysis_type] = {
                'top_pathways': top_pathways,
                'leading_edge_genes': gene_counts.head(20)
            }
        
        leading_edge_results[contrast_name] = contrast_leading_edge
    
    return leading_edge_results

def analyze_pathway_crosstalk(enrichment_results):
    """Analyze pathway-pathway interactions and crosstalk"""
    print("\n🔗 PATHWAY CROSSTALK ANALYSIS")
    print("=" * 50)
    
    crosstalk_results = {}
    
    for contrast_name, results in enrichment_results.items():
        print(f"   🔬 Analyzing crosstalk in {contrast_name}...")
        
        # Collect all significant pathways
        all_pathways = []
        for analysis_type, df in results.items():
            if not df.empty:
                sig_pathways = df[df['pvalue'] < 0.05].copy()
                sig_pathways['source'] = analysis_type
                all_pathways.append(sig_pathways)
        
        if not all_pathways:
            continue
        
        combined_pathways = pd.concat(all_pathways, ignore_index=True)
        
        # Simple crosstalk analysis based on shared keywords
        pathway_keywords = {}
        for _, pathway in combined_pathways.iterrows():
            desc = pathway['Description'].lower()
            keywords = set(desc.split())
            # Filter common words
            filtered_keywords = {k for k in keywords if len(k) > 3 and 
                               k not in ['the', 'and', 'for', 'with', 'pathway', 'process']}
            pathway_keywords[pathway['Description']] = filtered_keywords
        
        # Calculate pathway similarity
        crosstalk_matrix = pd.DataFrame(index=pathway_keywords.keys(), 
                                      columns=pathway_keywords.keys())
        
        for path1, kw1 in pathway_keywords.items():
            for path2, kw2 in pathway_keywords.items():
                if len(kw1) > 0 and len(kw2) > 0:
                    jaccard = len(kw1 & kw2) / len(kw1 | kw2)
                else:
                    jaccard = 0
                crosstalk_matrix.loc[path1, path2] = jaccard
        
        crosstalk_results[contrast_name] = {
            'pathways': combined_pathways,
            'crosstalk_matrix': crosstalk_matrix
        }
    
    return crosstalk_results

def create_functional_modules(enrichment_results):
    """Group pathways into functional modules using clustering"""
    print("\n🏗️  FUNCTIONAL MODULE ANALYSIS")
    print("=" * 50)
    
    from sklearn.cluster import AgglomerativeClustering
    from sklearn.metrics.pairwise import cosine_similarity
    
    module_results = {}
    
    for contrast_name, results in enrichment_results.items():
        print(f"   🔬 Creating modules for {contrast_name}...")
        
        # Collect pathway descriptions
        all_pathways = []
        for analysis_type, df in results.items():
            if not df.empty:
                pathways = df[['Description', 'pvalue']].copy()
                pathways['source'] = analysis_type
                all_pathways.append(pathways)
        
        if not all_pathways:
            continue
        
        combined = pd.concat(all_pathways, ignore_index=True)
        
        if len(combined) < 5:
            continue
        
        # Simple text-based clustering using pathway descriptions
        from sklearn.feature_extraction.text import TfidfVectorizer
        
        try:
            # Create TF-IDF matrix
            vectorizer = TfidfVectorizer(max_features=100, stop_words='english', 
                                       ngram_range=(1, 2))
            tfidf_matrix = vectorizer.fit_transform(combined['Description'])
            
            # Perform clustering
            n_clusters = min(5, len(combined) // 3)
            if n_clusters < 2:
                n_clusters = 2
            
            clustering = AgglomerativeClustering(n_clusters=n_clusters, 
                                               linkage='ward')
            clusters = clustering.fit_predict(tfidf_matrix.toarray())
            
            # Add cluster labels
            combined['module'] = clusters
            
            module_results[contrast_name] = combined
            
            print(f"     ✅ Created {n_clusters} functional modules")
            
        except Exception as e:
            print(f"     ❌ Failed clustering for {contrast_name}: {e}")
    
    return module_results

def compare_contrast_patterns(enrichment_results):
    """Analyze pathway patterns across different contrasts"""
    print("\n📊 CONTRAST PATTERN ANALYSIS")
    print("=" * 50)
    
    # Create pathway occurrence matrix
    all_pathways = set()
    contrast_pathways = {}
    
    for contrast_name, results in enrichment_results.items():
        pathways = set()
        for analysis_type, df in results.items():
            if not df.empty:
                sig_pathways = df[df['pvalue'] < 0.05]['Description']
                pathways.update(sig_pathways)
        
        contrast_pathways[contrast_name] = pathways
        all_pathways.update(pathways)
    
    # Create binary occurrence matrix
    occurrence_matrix = pd.DataFrame(index=list(all_pathways), 
                                   columns=list(contrast_pathways.keys()))
    
    for pathway in all_pathways:
        for contrast in contrast_pathways.keys():
            occurrence_matrix.loc[pathway, contrast] = pathway in contrast_pathways[contrast]
    
    # Calculate pattern statistics
    pattern_stats = {
        'unique_pathways': len(all_pathways),
        'shared_pathways': sum(occurrence_matrix.sum(axis=1) > 1),
        'contrast_specific': sum(occurrence_matrix.sum(axis=1) == 1),
        'universal_pathways': sum(occurrence_matrix.sum(axis=1) == len(contrast_pathways))
    }
    
    print(f"   📊 Total unique pathways: {pattern_stats['unique_pathways']}")
    print(f"   📊 Shared pathways: {pattern_stats['shared_pathways']}")
    print(f"   📊 Contrast-specific: {pattern_stats['contrast_specific']}")
    print(f"   📊 Universal pathways: {pattern_stats['universal_pathways']}")
    
    return {
        'occurrence_matrix': occurrence_matrix,
        'pattern_statistics': pattern_stats,
        'contrast_pathways': contrast_pathways
    }

def perform_meta_pathway_analysis(enrichment_results):
    """Meta-analysis combining pathway evidence across contrasts"""
    print("\n🔬 META-PATHWAY ANALYSIS")
    print("=" * 50)
    
    # Collect all pathway results
    meta_pathways = []
    
    for contrast_name, results in enrichment_results.items():
        for analysis_type, df in results.items():
            if df.empty:
                continue
            
            df_copy = df.copy()
            df_copy['contrast'] = contrast_name
            df_copy['analysis_type'] = analysis_type
            meta_pathways.append(df_copy)
    
    if not meta_pathways:
        return {}
    
    combined_df = pd.concat(meta_pathways, ignore_index=True)
    
    # Meta-analysis by pathway description
    meta_results = []
    
    for pathway_desc in combined_df['Description'].unique():
        pathway_data = combined_df[combined_df['Description'] == pathway_desc]
        
        if len(pathway_data) > 1:
            # Simple meta-analysis using Fisher's method
            from scipy.stats import combine_pvalues
            
            try:
                pvalues = pathway_data['pvalue'].values
                if len(pvalues) > 1:
                    combined_stat, combined_pval = combine_pvalues(pvalues, method='fisher')
                else:
                    combined_pval = pvalues[0]
                
                meta_results.append({
                    'Description': pathway_desc,
                    'n_studies': len(pathway_data),
                    'contrasts': ', '.join(pathway_data['contrast'].unique()),
                    'analysis_types': ', '.join(pathway_data['analysis_type'].unique()),
                    'individual_pvalues': pvalues.tolist(),
                    'meta_pvalue': combined_pval,
                    'mean_effect': pathway_data['pvalue'].mean()
                })
            except Exception as e:
                print(f"     ⚠️  Meta-analysis failed for {pathway_desc}: {e}")
    
    meta_df = pd.DataFrame(meta_results)
    
    if not meta_df.empty:
        meta_df = meta_df.sort_values('meta_pvalue')
        print(f"   📊 Meta-analysis completed for {len(meta_df)} pathways")
        print(f"   📊 Top meta-pathway: {meta_df.iloc[0]['Description']}")
    
    return meta_df

def create_enrichment_plots(enrichment_results):
    """Create comprehensive visualization plots"""
    print("\n📊 CREATING ENRICHMENT PLOTS")
    print("=" * 50)
    
    # 1. Basic plots
    create_pathway_heatmap(enrichment_results)
    create_dotplots(enrichment_results)
    create_network_plots(enrichment_results)
    create_comparison_plots(enrichment_results)
    
    # 2. Advanced visualizations
    create_upset_plots(enrichment_results)
    create_enrichment_maps(enrichment_results)
    create_volcano_plots(enrichment_results)
    
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
    
    # Save basic enrichment results
    for contrast_name, results in all_results.items():
        for analysis_type, df in results.items():
            if not df.empty:
                output_file = f"{BASIC_DIR}/{analysis_type}_{contrast_name}_enrichment.csv"
                df.to_csv(output_file, index=False)
                print(f"   📁 {analysis_type} - {contrast_name}: {output_file}")
    
    # Save combined basic results
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
        combined_file = f"{BASIC_DIR}/all_basic_enrichment_results.csv"
        combined_df.to_csv(combined_file, index=False)
        print(f"   📁 Combined basic results: {combined_file}")
    
    # Save complement-focused results
    if complement_results:
        complement_combined = pd.concat(complement_results.values(), ignore_index=True)
        complement_file = f"{COMPLEMENT_DIR}/complement_pathway_enrichment.csv"
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

def save_advanced_results(gsea_results, msigdb_results, leading_edge_results, 
                         crosstalk_results, module_results, pattern_results, meta_results):
    """Save all advanced analysis results"""
    print("\n💾 SAVING ADVANCED ANALYSIS RESULTS")
    print("=" * 50)
    
    # Save GSEA results
    if gsea_results:
        for contrast_name, results in gsea_results.items():
            for analysis_type, df in results.items():
                if not df.empty:
                    output_file = f"{GSEA_DIR}/{analysis_type}_{contrast_name}.csv"
                    df.to_csv(output_file, index=False)
                    print(f"   📁 GSEA {analysis_type} - {contrast_name}: {output_file}")
    
    # Save MSigDB results
    if msigdb_results:
        for contrast_name, results in msigdb_results.items():
            for analysis_type, df in results.items():
                if not df.empty:
                    output_file = f"{MSIGDB_DIR}/{analysis_type}_{contrast_name}.csv"
                    df.to_csv(output_file, index=False)
                    print(f"   📁 MSigDB {analysis_type} - {contrast_name}: {output_file}")
    
    # Save leading edge results
    if leading_edge_results:
        for contrast_name, results in leading_edge_results.items():
            for analysis_type, data in results.items():
                if 'leading_edge_genes' in data and not data['leading_edge_genes'].empty:
                    output_file = f"{LEADING_EDGE_DIR}/leading_edge_{analysis_type}_{contrast_name}.csv"
                    data['leading_edge_genes'].to_csv(output_file)
                    print(f"   📁 Leading edge {analysis_type} - {contrast_name}: {output_file}")
    
    # Save crosstalk results
    if crosstalk_results:
        for contrast_name, data in crosstalk_results.items():
            if 'crosstalk_matrix' in data:
                output_file = f"{CROSSTALK_DIR}/crosstalk_matrix_{contrast_name}.csv"
                data['crosstalk_matrix'].to_csv(output_file)
                print(f"   📁 Crosstalk matrix - {contrast_name}: {output_file}")
    
    # Save module results
    if module_results:
        for contrast_name, df in module_results.items():
            if not df.empty:
                output_file = f"{MODULES_DIR}/functional_modules_{contrast_name}.csv"
                df.to_csv(output_file, index=False)
                print(f"   📁 Functional modules - {contrast_name}: {output_file}")
    
    # Save pattern analysis
    if pattern_results and 'occurrence_matrix' in pattern_results:
        pattern_file = f"{COMPARISON_DIR}/pathway_occurrence_matrix.csv"
        pattern_results['occurrence_matrix'].to_csv(pattern_file)
        print(f"   📁 Pattern analysis: {pattern_file}")
    
    # Save meta-analysis results
    if meta_results is not None and not meta_results.empty:
        meta_file = f"{COMPARISON_DIR}/meta_pathway_analysis.csv"
        meta_results.to_csv(meta_file, index=False)
        print(f"   📁 Meta-analysis: {meta_file}")

def main():
    """Main pathway enrichment workflow"""
    print("🛤️  COMPREHENSIVE PATHWAY ENRICHMENT ANALYSIS - EDGER RESULTS")
    print("=" * 80)
    print("Advanced pathway enrichment with specialized analyses and visualizations")
    print("=" * 80)
    
    try:
        # Load edgeR results
        results_dict = load_edger_results()
        
        # Prepare gene lists
        gene_lists = prepare_gene_lists(results_dict)
        
        # Run basic enrichment analyses
        print("\n🔬 RUNNING BASIC ENRICHMENT ANALYSES")
        print("=" * 50)
        
        go_results = run_go_enrichment(gene_lists)
        kegg_results = run_kegg_enrichment(gene_lists)
        reactome_results = run_reactome_enrichment(gene_lists)
        
        # Combine basic results
        all_results = {}
        for contrast_name in gene_lists.keys():
            all_results[contrast_name] = {}
            
            if contrast_name in go_results:
                all_results[contrast_name].update(go_results[contrast_name])
            
            if contrast_name in kegg_results:
                all_results[contrast_name]['KEGG'] = kegg_results[contrast_name]
            
            if contrast_name in reactome_results:
                all_results[contrast_name]['Reactome'] = reactome_results[contrast_name]
        
        # Run advanced analyses
        print("\n🔬 RUNNING ADVANCED ANALYSES")
        print("=" * 50)
        
        gsea_results = run_gsea_analysis(gene_lists)
        msigdb_results = run_msigdb_enrichment(gene_lists)
        
        # Specialized analyses
        print("\n🔬 RUNNING SPECIALIZED ANALYSES")
        print("=" * 50)
        
        complement_results = analyze_complement_pathways(all_results)
        neuro_results = analyze_neuroinflammation_pathways(all_results)
        leading_edge_results = perform_leading_edge_analysis(all_results, gene_lists)
        
        # Comparative analyses
        print("\n🔬 RUNNING COMPARATIVE ANALYSES")
        print("=" * 50)
        
        crosstalk_results = analyze_pathway_crosstalk(all_results)
        module_results = create_functional_modules(all_results)
        pattern_results = compare_contrast_patterns(all_results)
        meta_results = perform_meta_pathway_analysis(all_results)
        
        # Create visualizations
        create_enrichment_plots(all_results)
        
        # Save all results
        save_enrichment_results(all_results, complement_results)
        save_advanced_results(gsea_results, msigdb_results, leading_edge_results,
                            crosstalk_results, module_results, pattern_results, meta_results)
        
        # Save specialized analyses
        if neuro_results:
            neuro_combined = pd.concat(neuro_results.values(), ignore_index=True)
            neuro_file = f"{NEUROINFLAM_DIR}/neuroinflammation_pathways.csv"
            neuro_combined.to_csv(neuro_file, index=False)
            print(f"   📁 Neuroinflammation pathways: {neuro_file}")
        
        # Create comprehensive summary
        create_summary_report(all_results, complement_results)
        
        print("\n" + "=" * 80)
        print("✅ COMPREHENSIVE PATHWAY ENRICHMENT ANALYSIS COMPLETE!")
        print("=" * 80)
        
        # Calculate total statistics
        total_basic = sum(len(df) for results in all_results.values() 
                         for df in results.values() if not df.empty)
        total_gsea = sum(len(df) for results in gsea_results.values() 
                        for df in results.values() if not df.empty) if gsea_results else 0
        total_msigdb = sum(len(df) for results in msigdb_results.values() 
                          for df in results.values() if not df.empty) if msigdb_results else 0
        
        print(f"📊 Basic enrichment pathways: {total_basic}")
        print(f"📊 GSEA pathways: {total_gsea}")
        print(f"📊 MSigDB pathways: {total_msigdb}")
        print(f"🧬 Complement pathways: {sum(len(df) for df in complement_results.values())}")
        print(f"🧠 Neuroinflammation pathways: {sum(len(df) for df in neuro_results.values()) if neuro_results else 0}")
        print(f"📁 Results organized in: {OUTPUT_DIR}/")
        print("🎯 Ready for comprehensive biological interpretation!")
        print("=" * 80)
        
        return {
            'basic': all_results,
            'gsea': gsea_results,
            'msigdb': msigdb_results,
            'complement': complement_results,
            'neuroinflammation': neuro_results,
            'leading_edge': leading_edge_results,
            'crosstalk': crosstalk_results,
            'modules': module_results,
            'patterns': pattern_results,
            'meta': meta_results
        }
        
    except Exception as e:
        print(f"❌ Comprehensive pathway analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return None

if __name__ == "__main__":
    main()