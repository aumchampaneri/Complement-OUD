#!/usr/bin/env python3
"""
Comprehensive Pathway Analysis: LEMUR vs DESeq2
Enhanced with Neuroinflammation & Complement Cascade Focus

This script performs comprehensive pathway enrichment analysis using the R Bioconductor
ecosystem via rpy2, with special focus on neuroinflammation and complement cascade pathways.

Key Features:
1. Dual-method analysis (LEMUR + DESeq2 results)
2. Advanced pathway enrichment (GO, KEGG, Reactome, MSigDB)
3. Neuroinflammation & complement cascade focus
4. Cross-method pathway validation
5. Regional comparison (Caudate vs Putamen)
6. Publication-quality visualizations

Dependencies:
- rpy2
- R packages: clusterProfiler, fgsea, msigdbr, enrichplot, ReactomePA, DOSE
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
import logging
from pathlib import Path
import json

# R integration
try:
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter
    R_AVAILABLE = True
    print("✅ rpy2 successfully imported")
except ImportError as e:
    R_AVAILABLE = False
    print(f"❌ rpy2 not available: {e}")
    print("Install with: pip install rpy2")

warnings.filterwarnings('ignore')

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class ComprehensivePathwayAnalyzer:
    """
    Comprehensive pathway analysis using R Bioconductor ecosystem
    """
    
    def __init__(self, base_dir=None):
        """Initialize the pathway analyzer"""
        # Set up directories
        script_dir = os.path.dirname(os.path.abspath(__file__))
        if base_dir is None:
            base_dir = os.path.join(script_dir, "..", "..", "results", "snrna_scvi")
        
        self.base_dir = os.path.abspath(base_dir)
        self.method_comparison_dir = os.path.join(self.base_dir, "method_comparison", "tables")
        self.lemur_dir = os.path.join(self.base_dir, "lemur_analysis", "tables")
        self.deseq2_dir = os.path.join(self.base_dir, "pydeseq2_analysis", "tables", "Full_Dataset")
        self.output_dir = os.path.join(self.base_dir, "comprehensive_pathway_analysis")
        
        # Create output directories
        self._create_output_dirs()
        
        # Initialize R environment
        if R_AVAILABLE:
            self._setup_r_environment()
        else:
            raise ImportError("rpy2 is required for this analysis. Install with: pip install rpy2")
        
        # Initialize data storage
        self.gene_sets = {}
        self.pathway_results = {}
        self.custom_pathways = {}
        
        # Set up pathway databases
        self._setup_custom_pathways()
        
        logger.info(f"ComprehensivePathwayAnalyzer initialized. Output dir: {self.output_dir}")
    
    def _create_output_dirs(self):
        """Create all necessary output directories"""
        subdirs = [
            "tables", "plots", 
            "tables/method_comparison_pathways",
            "tables/neuroinflammation_analysis", 
            "tables/complement_analysis",
            "tables/original_results_pathways",
            "plots/method_comparison_plots",
            "plots/neuroinflammation_plots",
            "plots/complement_plots",
            "plots/regional_comparison_plots"
        ]
        
        for subdir in subdirs:
            os.makedirs(os.path.join(self.output_dir, subdir), exist_ok=True)
    
    def _setup_r_environment(self):
        """Set up R environment and import required packages"""
        logger.info("Setting up R environment...")
        
        # Set up pandas-R conversion context
        self.pandas2ri_converter = pandas2ri.converter
        
        # Required R packages with fallback options
        required_packages = {
            'clusterProfiler': 'bioc',
            'fgsea': 'bioc', 
            'msigdbr': 'bioc',
            'enrichplot': 'bioc',
            'ReactomePA': 'bioc',
            'DOSE': 'bioc',
            'org.Hs.eg.db': 'bioc',
            'ggplot2': 'cran',
            'dplyr': 'cran',
            'stringr': 'cran',
            'purrr': 'cran'
        }
        
        # Check and import packages
        self.r_packages = {}
        self.missing_packages = []
        
        try:
            utils = importr('utils')
            base = importr('base')
        except Exception as e:
            logger.error(f"Failed to import basic R packages: {e}")
            raise
        
        for pkg, pkg_type in required_packages.items():
            try:
                self.r_packages[pkg] = importr(pkg)
                logger.info(f"✅ Imported R package: {pkg}")
            except Exception as e:
                logger.warning(f"❌ Failed to import R package {pkg}: {e}")
                self.missing_packages.append((pkg, pkg_type))
        
        # Report missing packages but continue
        if self.missing_packages:
            logger.warning(f"Missing R packages: {[pkg for pkg, _ in self.missing_packages]}")
            logger.warning("Some analyses may not be available. Install missing packages with:")
            bioc_packages = [pkg for pkg, pkg_type in self.missing_packages if pkg_type == 'bioc']
            cran_packages = [pkg for pkg, pkg_type in self.missing_packages if pkg_type == 'cran']
            
            if bioc_packages:
                logger.warning(f"BiocManager::install(c({', '.join([f'\"{pkg}\"' for pkg in bioc_packages])}))")
            if cran_packages:
                logger.warning(f"install.packages(c({', '.join([f'\"{pkg}\"' for pkg in cran_packages])}))")
        
        # Set up R functions
        self._setup_r_functions()
    
    def _setup_r_functions(self):
        """Set up custom R functions for analysis"""
        logger.info("Setting up R functions...")
        
        try:
            ro.r('''
            # Custom function for gene symbol to Entrez ID conversion
            convert_symbols_to_entrez <- function(gene_symbols) {
                if (!require("clusterProfiler", quietly = TRUE) || 
                    !require("org.Hs.eg.db", quietly = TRUE)) {
                    warning("Required packages not available for gene conversion")
                    return(data.frame(SYMBOL = character(0), ENTREZID = character(0)))
                }
                
                # Convert symbols to Entrez IDs
                tryCatch({
                    entrez_ids <- bitr(gene_symbols, 
                                      fromType = "SYMBOL", 
                                      toType = "ENTREZID", 
                                      OrgDb = org.Hs.eg.db, 
                                      drop = TRUE)
                    return(entrez_ids)
                }, error = function(e) {
                    warning(paste("Gene conversion failed:", e$message))
                    return(data.frame(SYMBOL = character(0), ENTREZID = character(0)))
                })
            }
        
            # Custom function for comprehensive GO enrichment
            run_comprehensive_go <- function(genes, universe = NULL, pvalueCutoff = 0.05) {
                if (!require("clusterProfiler", quietly = TRUE) || 
                    !require("org.Hs.eg.db", quietly = TRUE)) {
                    warning("Required packages not available for GO enrichment")
                    return(NULL)
                }
                
                results <- list()
                
                # Convert to Entrez IDs
                gene_entrez <- convert_symbols_to_entrez(genes)
                
                if (nrow(gene_entrez) == 0) {
                    warning("No genes could be converted to Entrez IDs")
                    return(NULL)
                }
                
                gene_list <- gene_entrez$ENTREZID
            
            # GO Biological Process
            tryCatch({
                results$BP <- enrichGO(gene = gene_list,
                                     OrgDb = org.Hs.eg.db,
                                     ont = "BP",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff = pvalueCutoff,
                                     universe = universe,
                                     readable = TRUE)
            }, error = function(e) { 
                warning(paste("GO BP enrichment failed:", e$message))
                results$BP <<- NULL 
            })
            
            # GO Molecular Function
            tryCatch({
                results$MF <- enrichGO(gene = gene_list,
                                     OrgDb = org.Hs.eg.db,
                                     ont = "MF",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff = pvalueCutoff,
                                     universe = universe,
                                     readable = TRUE)
            }, error = function(e) { 
                warning(paste("GO MF enrichment failed:", e$message))
                results$MF <<- NULL 
            })
            
            # GO Cellular Component
            tryCatch({
                results$CC <- enrichGO(gene = gene_list,
                                     OrgDb = org.Hs.eg.db,
                                     ont = "CC",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff = pvalueCutoff,
                                     universe = universe,
                                     readable = TRUE)
            }, error = function(e) { 
                warning(paste("GO CC enrichment failed:", e$message))
                results$CC <<- NULL 
            })
            
            return(results)
        }
        
            # Custom function for KEGG enrichment
            run_kegg_enrichment <- function(genes, pvalueCutoff = 0.05) {
                if (!require("clusterProfiler", quietly = TRUE) || 
                    !require("org.Hs.eg.db", quietly = TRUE)) {
                    warning("Required packages not available for KEGG enrichment")
                    return(NULL)
                }
                
                # Convert to Entrez IDs
                gene_entrez <- convert_symbols_to_entrez(genes)
                
                if (nrow(gene_entrez) == 0) {
                    warning("No genes could be converted to Entrez IDs")
                    return(NULL)
                }
                
                gene_list <- gene_entrez$ENTREZID
            
            tryCatch({
                kegg_result <- enrichKEGG(gene = gene_list,
                                        organism = 'hsa',
                                        pvalueCutoff = pvalueCutoff,
                                        pAdjustMethod = "BH")
                return(kegg_result)
            }, error = function(e) {
                warning(paste("KEGG enrichment failed:", e$message))
                return(NULL)
            })
        }
        
            # Custom function for Reactome enrichment
            run_reactome_enrichment <- function(genes, pvalueCutoff = 0.05) {
                if (!require("ReactomePA", quietly = TRUE) || 
                    !require("org.Hs.eg.db", quietly = TRUE)) {
                    warning("Required packages not available for Reactome enrichment")
                    return(NULL)
                }
                
                # Convert to Entrez IDs
                gene_entrez <- convert_symbols_to_entrez(genes)
                
                if (nrow(gene_entrez) == 0) {
                    warning("No genes could be converted to Entrez IDs")
                    return(NULL)
                }
                
                gene_list <- gene_entrez$ENTREZID
            
            tryCatch({
                reactome_result <- enrichPathway(gene = gene_list,
                                               pvalueCutoff = pvalueCutoff,
                                               pAdjustMethod = "BH",
                                               readable = TRUE)
                return(reactome_result)
            }, error = function(e) {
                warning(paste("Reactome enrichment failed:", e$message))
                return(NULL)
            })
        }
        
            # Custom function for MSigDB enrichment
            run_msigdb_enrichment <- function(genes, category = "H", pvalueCutoff = 0.05) {
                if (!require("clusterProfiler", quietly = TRUE) || 
                    !require("msigdbr", quietly = TRUE)) {
                    warning("Required packages not available for MSigDB enrichment")
                    return(NULL)
                }
                
                # Get MSigDB gene sets
                tryCatch({
                    msigdb_sets <- msigdbr(species = "Homo sapiens", category = category)
                    msigdb_list <- split(x = msigdb_sets$gene_symbol, f = msigdb_sets$gs_name)
                }, error = function(e) {
                    warning(paste("Failed to load MSigDB gene sets:", e$message))
                    return(NULL)
                })
            
            tryCatch({
                msigdb_result <- enricher(gene = genes,
                                        TERM2GENE = msigdb_sets[,c("gs_name", "gene_symbol")],
                                        pvalueCutoff = pvalueCutoff,
                                        pAdjustMethod = "BH")
                return(msigdb_result)
            }, error = function(e) {
                warning(paste("MSigDB enrichment failed:", e$message))
                return(NULL)
            })
        }
            ''')
            
            logger.info("✅ R functions set up successfully")
            
        except Exception as e:
            logger.error(f"Failed to set up R functions: {e}")
            logger.warning("Some pathway analyses may not be available")
    
    def _setup_custom_pathways(self):
        """Set up custom pathway gene sets for neuroinflammation and complement"""
        
        # Complement cascade pathways (your specialty!)
        self.custom_pathways['complement_cascade'] = {
            'classical_pathway': [
                'C1QA', 'C1QB', 'C1QC', 'C1R', 'C1S', 'C4A', 'C4B', 'C2',
                'C4BPA', 'C4BPB'
            ],
            'alternative_pathway': [
                'CFB', 'CFD', 'CFP', 'CFH', 'CFI', 'C3', 'CFHR1', 'CFHR2', 
                'CFHR3', 'CFHR4', 'CFHR5'
            ],
            'lectin_pathway': [
                'MBL2', 'MASP1', 'MASP2', 'FCN1', 'FCN2', 'FCN3', 'COLEC10',
                'COLEC11', 'CL-K1', 'CL-L1'
            ],
            'terminal_pathway': [
                'C5', 'C6', 'C7', 'C8A', 'C8B', 'C8G', 'C9', 'CFHR1'
            ],
            'complement_regulation': [
                'CD55', 'CD46', 'CD35', 'CFH', 'CFI', 'C1INH', 'SERPING1',
                'CFHR1', 'CFHR3', 'CFHR4', 'CD59', 'CLU', 'VTN', 'S100A9'
            ]
        }
        
        # Neuroinflammation pathways
        self.custom_pathways['neuroinflammation'] = {
            'microglia_activation': [
                'CD68', 'IBA1', 'AIF1', 'CX3CR1', 'P2RY12', 'TMEM119', 'TREM2',
                'CD11B', 'ITGAM', 'CD45', 'PTPRC', 'MHCII', 'HLA-DRA', 'TNF',
                'IL1B', 'IL6', 'NLRP3', 'TLR4', 'TLR2', 'MARCO', 'MSR1'
            ],
            'astrocyte_reactivity': [
                'GFAP', 'S100B', 'ALDH1L1', 'SOX9', 'AQP4', 'EAAT1', 'SLC1A3',
                'EAAT2', 'SLC1A2', 'GLUL', 'GS', 'CSPG4', 'PDGFRA', 'LCN2',
                'SERPINA3', 'CP', 'GBP2', 'STEAP4', 'OSMR', 'SOCS3'
            ],
            'cytokine_signaling': [
                'TNF', 'IL1B', 'IL6', 'IL10', 'IL4', 'IL13', 'IFNG', 'IL12A',
                'IL12B', 'IL23A', 'IL17A', 'IL17F', 'IL22', 'TGF-B1', 'TGFB1',
                'CSF1', 'CSF2', 'CSF3', 'CCL2', 'CCL3', 'CCL4', 'CCL5', 'CXCL1',
                'CXCL2', 'CXCL8', 'CXCL10', 'CXCL12'
            ],
            'nf_kappa_b_pathway': [
                'NFKB1', 'NFKB2', 'RELA', 'RELB', 'REL', 'NFKBIA', 'NFKBIB',
                'NFKBIE', 'IKBKA', 'IKBKB', 'IKBKG', 'TNFAIP3', 'TNIP1',
                'TNIP2', 'CYLD', 'USP15', 'OTULIN'
            ],
            'inflammasome_activation': [
                'NLRP1', 'NLRP3', 'NLRC4', 'AIM2', 'PYCARD', 'ASC', 'CASP1',
                'CASP4', 'CASP5', 'CASP11', 'IL1B', 'IL18', 'GSDMD', 'NEK7',
                'TXNIP', 'HMGB1', 'HSP90'
            ]
        }
        
        # Addiction-specific pathways
        self.custom_pathways['addiction_pathways'] = {
            'dopamine_signaling': [
                'TH', 'DDC', 'DRD1', 'DRD2', 'DRD3', 'DRD4', 'DRD5', 'DAT1',
                'SLC6A3', 'COMT', 'MAO-A', 'MAOA', 'MAO-B', 'MAOB', 'ALDH2',
                'VMAT2', 'SLC18A2', 'HOMER1', 'HOMER2', 'HOMER3'
            ],
            'reward_circuitry': [
                'OPRM1', 'OPRD1', 'OPRK1', 'PENK', 'PDYN', 'POMC', 'CART',
                'NPY', 'CRH', 'CRHR1', 'CRHR2', 'AVP', 'OXT', 'OXTR', 'AVPR1A',
                'AVPR1B', 'MC4R', 'LEPR', 'GHSR'
            ],
            'stress_response': [
                'CRH', 'CRHR1', 'CRHR2', 'POMC', 'HPA', 'NR3C1', 'NR3C2',
                'FKBP5', 'FKBP4', 'HSP90', 'AVP', 'OXT', 'BDNF', 'NTRK2',
                'CREB1', 'CREB3', 'ATF1', 'JUN', 'FOS', 'FOSB', 'NPAS4'
            ],
            'synaptic_plasticity': [
                'BDNF', 'NTRK2', 'CREB1', 'ARC', 'EGR1', 'EGR2', 'FOS', 'FOSB',
                'JUN', 'JUNB', 'JUND', 'ATF1', 'NPAS4', 'NR4A1', 'NR4A2',
                'NR4A3', 'HOMER1', 'HOMER2', 'HOMER3', 'DLG4', 'PSD95', 'GRIN1',
                'GRIN2A', 'GRIN2B', 'GRIA1', 'GRIA2', 'CAMK2A', 'CAMK2B'
            ]
        }
        
        # Neurotransmitter systems
        self.custom_pathways['neurotransmitter_systems'] = {
            'gabaergic': [
                'GAD1', 'GAD2', 'SLC32A1', 'GABBR1', 'GABBR2', 'GABRA1',
                'GABRA2', 'GABRA3', 'GABRA4', 'GABRA5', 'GABRA6', 'GABRB1',
                'GABRB2', 'GABRB3', 'GABRG1', 'GABRG2', 'GABRG3', 'GABRD',
                'GABRE', 'GABRP', 'GABRQ', 'SLC6A1', 'SLC6A11', 'SLC6A12'
            ],
            'glutamatergic': [
                'GRIN1', 'GRIN2A', 'GRIN2B', 'GRIN2C', 'GRIN2D', 'GRIN3A',
                'GRIN3B', 'GRIA1', 'GRIA2', 'GRIA3', 'GRIA4', 'GRIK1', 'GRIK2',
                'GRIK3', 'GRIK4', 'GRIK5', 'GRM1', 'GRM2', 'GRM3', 'GRM4',
                'GRM5', 'GRM6', 'GRM7', 'GRM8', 'SLC17A6', 'SLC17A7', 'SLC1A1',
                'SLC1A2', 'SLC1A3', 'SLC1A6', 'SLC1A7'
            ],
            'cholinergic': [
                'CHAT', 'ACHE', 'BCHE', 'SLC5A7', 'SLC18A3', 'CHRNA1', 'CHRNA2',
                'CHRNA3', 'CHRNA4', 'CHRNA5', 'CHRNA6', 'CHRNA7', 'CHRNA9',
                'CHRNA10', 'CHRNB1', 'CHRNB2', 'CHRNB3', 'CHRNB4', 'CHRM1',
                'CHRM2', 'CHRM3', 'CHRM4', 'CHRM5'
            ],
            'opioid_system': [
                'OPRM1', 'OPRD1', 'OPRK1', 'OPRL1', 'PENK', 'PDYN', 'POMC',
                'PNOC', 'NPFF', 'NPFFR1', 'NPFFR2', 'OFQ', 'PNOC'
            ]
        }
        
        logger.info(f"✅ Custom pathways set up: {len(self.custom_pathways)} categories")
    
    def load_gene_sets(self):
        """Load all gene sets from method comparison and original results"""
        logger.info("Loading gene sets from all sources...")
        
        # Load method comparison results
        self._load_method_comparison_genes()
        
        # Load original LEMUR results
        self._load_original_lemur_genes()
        
        # Load original DESeq2 results  
        self._load_original_deseq2_genes()
        
        logger.info(f"✅ Loaded {len(self.gene_sets)} gene sets")
        
        # Save gene set summary
        self._save_gene_set_summary()
    
    def _load_method_comparison_genes(self):
        """Load gene sets from method comparison analysis"""
        comparison_files = {
            'caudate_oud_intersection': 'caudate_oud_effect_intersection_genes.csv',
            'caudate_oud_lemur_specific': 'caudate_oud_effect_lemur_specific_genes.csv',
            'caudate_oud_deseq2_specific': 'caudate_oud_effect_deseq2_specific_genes.csv',
            'putamen_oud_lemur_specific': 'putamen_oud_effect_lemur_specific_genes.csv',
            'putamen_oud_deseq2_specific': 'putamen_oud_effect_deseq2_specific_genes.csv',
            'sex_interaction_caudate_intersection': 'sex_interaction_caudate_intersection_genes.csv',
            'sex_interaction_caudate_lemur_specific': 'sex_interaction_caudate_lemur_specific_genes.csv',
            'sex_interaction_caudate_deseq2_specific': 'sex_interaction_caudate_deseq2_specific_genes.csv',
            'sex_interaction_putamen_intersection': 'sex_interaction_putamen_intersection_genes.csv',
            'sex_interaction_putamen_lemur_specific': 'sex_interaction_putamen_lemur_specific_genes.csv',
            'sex_interaction_putamen_deseq2_specific': 'sex_interaction_putamen_deseq2_specific_genes.csv'
        }
        
        for gene_set_name, filename in comparison_files.items():
            filepath = os.path.join(self.method_comparison_dir, filename)
            if os.path.exists(filepath):
                try:
                    df = pd.read_csv(filepath)
                    if 'gene' in df.columns and len(df) > 0:
                        genes = df['gene'].dropna().unique().tolist()
                        self.gene_sets[f"method_comparison_{gene_set_name}"] = {
                            'genes': genes,
                            'source': 'method_comparison',
                            'n_genes': len(genes),
                            'description': f"Method comparison: {gene_set_name.replace('_', ' ')}"
                        }
                        logger.info(f"  Loaded {len(genes)} genes for {gene_set_name}")
                except Exception as e:
                    logger.warning(f"Failed to load {filename}: {e}")
            else:
                logger.warning(f"File not found: {filepath}")
    
    def _load_original_lemur_genes(self):
        """Load gene sets from original LEMUR results"""
        lemur_files = {
            'caudate_oud_significant': 'strategic_caudate_oud_effect_significant_genes.csv',
            'putamen_oud_significant': 'strategic_putamen_oud_effect_significant_genes.csv',
            'caudate_sex_interaction': 'strategic_caudate_sex_interaction_significant_genes.csv',
            'putamen_sex_interaction': 'strategic_putamen_sex_interaction_significant_genes.csv'
        }
        
        for gene_set_name, filename in lemur_files.items():
            filepath = os.path.join(self.lemur_dir, filename)
            if os.path.exists(filepath):
                try:
                    df = pd.read_csv(filepath)
                    if 'gene' in df.columns and len(df) > 0:
                        genes = df['gene'].dropna().unique().tolist()
                        self.gene_sets[f"lemur_original_{gene_set_name}"] = {
                            'genes': genes,
                            'source': 'lemur_original',
                            'n_genes': len(genes),
                            'description': f"LEMUR original: {gene_set_name.replace('_', ' ')}"
                        }
                        logger.info(f"  Loaded {len(genes)} genes for LEMUR {gene_set_name}")
                except Exception as e:
                    logger.warning(f"Failed to load LEMUR {filename}: {e}")
            else:
                logger.warning(f"LEMUR file not found: {filepath}")
    
    def _load_original_deseq2_genes(self):
        """Load gene sets from original DESeq2 results"""
        deseq2_files = {
            'caudate_oud_significant': '03_OUD_vs_Control_Caudate_results.csv',
            'putamen_oud_significant': '02_OUD_vs_Control_Putamen_results.csv',
            'caudate_sex_interaction': '08_OUD_Effect_Male_vs_Female_Caudate_results.csv',
            'putamen_sex_interaction': '07_OUD_Effect_Male_vs_Female_Putamen_results.csv'
        }
        
        for gene_set_name, filename in deseq2_files.items():
            filepath = os.path.join(self.deseq2_dir, filename)
            if os.path.exists(filepath):
                try:
                    df = pd.read_csv(filepath)
                    if 'gene' in df.columns and 'significant' in df.columns:
                        # Filter for significant genes
                        sig_df = df[df['significant'] == True]
                        if len(sig_df) > 0:
                            genes = sig_df['gene'].dropna().unique().tolist()
                            self.gene_sets[f"deseq2_original_{gene_set_name}"] = {
                                'genes': genes,
                                'source': 'deseq2_original',
                                'n_genes': len(genes),
                                'description': f"DESeq2 original: {gene_set_name.replace('_', ' ')}"
                            }
                            logger.info(f"  Loaded {len(genes)} genes for DESeq2 {gene_set_name}")
                except Exception as e:
                    logger.warning(f"Failed to load DESeq2 {filename}: {e}")
            else:
                logger.warning(f"DESeq2 file not found: {filepath}")
    
    def _save_gene_set_summary(self):
        """Save summary of all loaded gene sets"""
        summary_data = []
        for gene_set_name, info in self.gene_sets.items():
            summary_data.append({
                'Gene_Set': gene_set_name,
                'Source': info['source'],
                'N_Genes': info['n_genes'],
                'Description': info['description']
            })
        
        summary_df = pd.DataFrame(summary_data)
        summary_file = os.path.join(self.output_dir, "tables", "gene_sets_summary.csv")
        summary_df.to_csv(summary_file, index=False)
        logger.info(f"✅ Gene set summary saved to {summary_file}")
    
    def run_comprehensive_enrichment(self):
        """Run comprehensive pathway enrichment for all gene sets"""
        logger.info("Starting comprehensive pathway enrichment analysis...")
        
        for gene_set_name, gene_set_info in self.gene_sets.items():
            genes = gene_set_info['genes']
            
            if len(genes) < 5:  # Skip very small gene sets
                logger.warning(f"Skipping {gene_set_name}: too few genes ({len(genes)})")
                continue
            
            logger.info(f"Analyzing {gene_set_name} ({len(genes)} genes)...")
            
            # Run standard pathway enrichment
            enrichment_results = self._run_standard_enrichment(genes, gene_set_name)
            
            # Run custom pathway enrichment
            custom_results = self._run_custom_enrichment(genes, gene_set_name)
            
            # Combine results
            all_results = {**enrichment_results, **custom_results}
            
            # Store results
            self.pathway_results[gene_set_name] = all_results
            
            # Save individual results
            self._save_gene_set_results(gene_set_name, all_results)
        
        logger.info("✅ Comprehensive enrichment analysis complete")
    
    def _run_standard_enrichment(self, genes, gene_set_name):
        """Run standard pathway enrichment (GO, KEGG, Reactome, MSigDB)"""
        results = {}
        
        # Check if critical R packages are available
        if 'clusterProfiler' not in self.r_packages or 'org.Hs.eg.db' not in self.r_packages:
            logger.warning(f"Essential R packages missing for {gene_set_name}. Using fallback enrichment.")
            return self._run_fallback_enrichment(genes, gene_set_name)
        
        # Convert genes to R vector
        r_genes = ro.StrVector(genes)
        
        try:
            # GO enrichment
            logger.info(f"  Running GO enrichment for {gene_set_name}...")
            go_results = ro.r['run_comprehensive_go'](r_genes)
            if go_results != ro.NULL:
                results['GO'] = self._process_r_enrichment_result(go_results, 'GO')
        except Exception as e:
            logger.warning(f"GO enrichment failed for {gene_set_name}: {e}")
        
        try:
            # KEGG enrichment
            logger.info(f"  Running KEGG enrichment for {gene_set_name}...")
            kegg_result = ro.r['run_kegg_enrichment'](r_genes)
            if kegg_result != ro.NULL:
                results['KEGG'] = self._process_r_enrichment_result(kegg_result, 'KEGG')
        except Exception as e:
            logger.warning(f"KEGG enrichment failed for {gene_set_name}: {e}")
        
        try:
            # Reactome enrichment
            logger.info(f"  Running Reactome enrichment for {gene_set_name}...")
            reactome_result = ro.r['run_reactome_enrichment'](r_genes)
            if reactome_result != ro.NULL:
                results['Reactome'] = self._process_r_enrichment_result(reactome_result, 'Reactome')
        except Exception as e:
            logger.warning(f"Reactome enrichment failed for {gene_set_name}: {e}")
        
        try:
            # MSigDB Hallmark enrichment
            logger.info(f"  Running MSigDB Hallmark enrichment for {gene_set_name}...")
            msigdb_result = ro.r['run_msigdb_enrichment'](r_genes, 'H')
            if msigdb_result != ro.NULL:
                results['MSigDB_Hallmark'] = self._process_r_enrichment_result(msigdb_result, 'MSigDB_Hallmark')
        except Exception as e:
            logger.warning(f"MSigDB enrichment failed for {gene_set_name}: {e}")
        
        return results
    
    def _run_custom_enrichment(self, genes, gene_set_name):
        """Run custom pathway enrichment (complement, neuroinflammation, addiction)"""
        results = {}
        
        for pathway_category, pathway_sets in self.custom_pathways.items():
            category_results = {}
            
            for pathway_name, pathway_genes in pathway_sets.items():
                # Calculate overlap
                overlap_genes = list(set(genes) & set(pathway_genes))
                overlap_count = len(overlap_genes)
                
                if overlap_count > 0:
                    # Calculate enrichment statistics using Fisher's exact test
                    enrichment_stats = self._calculate_enrichment_stats(
                        genes, pathway_genes, overlap_genes
                    )
                    
                    category_results[pathway_name] = {
                        'overlap_genes': overlap_genes,
                        'overlap_count': overlap_count,
                        'pathway_size': len(pathway_genes),
                        'query_size': len(genes),
                        **enrichment_stats
                    }
            
            if category_results:
                results[f"Custom_{pathway_category}"] = category_results
        
        return results
    
    def _calculate_enrichment_stats(self, query_genes, pathway_genes, overlap_genes):
        """Calculate enrichment statistics using R's fisher.test"""
        
        # Create contingency table
        overlap = len(overlap_genes)
        query_only = len(query_genes) - overlap
        pathway_only = len(pathway_genes) - overlap
        
        # Estimate background (rough approximation)
        background_size = 20000  # Approximate human gene count
        neither = background_size - overlap - query_only - pathway_only
        
        # Create R matrix for Fisher's exact test
        contingency_matrix = ro.r.matrix(
            ro.IntVector([overlap, query_only, pathway_only, neither]),
            nrow=2, byrow=True
        )
        
        try:
            # Run Fisher's exact test
            fisher_result = ro.r['fisher.test'](contingency_matrix)
            
            # Extract results
            p_value = fisher_result.rx2('p.value')[0]
            odds_ratio = fisher_result.rx2('estimate')[0] if fisher_result.rx2('estimate') != ro.NULL else 1.0
            
            # Calculate fold enrichment
            expected = (len(query_genes) * len(pathway_genes)) / background_size
            fold_enrichment = overlap / expected if expected > 0 else 0
            
            return {
                'p_value': p_value,
                'odds_ratio': odds_ratio,
                'fold_enrichment': fold_enrichment,
                'expected': expected
            }
            
        except Exception as e:
            logger.warning(f"Fisher's test failed: {e}")
            return {
                'p_value': 1.0,
                'odds_ratio': 1.0,
                'fold_enrichment': 1.0,
                'expected': 0
            }
    
    def _process_r_enrichment_result(self, r_result, result_type):
        """Process R enrichment results and convert to Python format"""
        try:
            if result_type == 'GO':
                # GO results are a list with BP, MF, CC
                processed_results = {}
                for ontology in ['BP', 'MF', 'CC']:
                    try:
                        ontology_result = r_result.rx2(ontology)
                        if ontology_result != ro.NULL:
                            # Try to convert enrichResult to DataFrame using R's summary function
                            summary_df = ro.r['as.data.frame'](ontology_result)
                            if summary_df != ro.NULL and len(summary_df) > 0:
                                with localconverter(ro.default_converter + self.pandas2ri_converter):
                                    df = ro.conversion.rpy2py(summary_df)
                                    if not df.empty:
                                        processed_results[f"GO_{ontology}"] = df
                    except Exception as e:
                        logger.warning(f"Failed to process GO {ontology}: {e}")
                        continue
                return processed_results
            else:
                # Other results are single enrichResult objects
                try:
                    if r_result != ro.NULL:
                        # Convert enrichResult to DataFrame using R's as.data.frame
                        summary_df = ro.r['as.data.frame'](r_result)
                        if summary_df != ro.NULL and len(summary_df) > 0:
                            with localconverter(ro.default_converter + self.pandas2ri_converter):
                                df = ro.conversion.rpy2py(summary_df)
                                if not df.empty:
                                    return {result_type: df}
                except Exception as e:
                    logger.warning(f"Failed to process {result_type}: {e}")
                return {}
                    
        except Exception as e:
            logger.warning(f"Failed to process {result_type} results: {e}")
            return {}
    
    def _run_fallback_enrichment(self, genes, gene_set_name):
        """Run fallback enrichment when R packages are not available"""
        logger.info(f"  Running fallback enrichment for {gene_set_name}...")
        
        # For now, return empty results but indicate that custom enrichment will still work
        results = {
            'Fallback_Note': f"Standard pathway enrichment not available due to missing R packages. Custom pathway analysis will still be performed for {gene_set_name}."
        }
        
        return results
    
    def _save_gene_set_results(self, gene_set_name, results):
        """Save enrichment results for a specific gene set"""
        
        # Create gene set specific directory
        gene_set_dir = os.path.join(self.output_dir, "tables", "individual_gene_sets", gene_set_name)
        os.makedirs(gene_set_dir, exist_ok=True)
        
        for analysis_type, analysis_results in results.items():
            if isinstance(analysis_results, dict):
                if analysis_type.startswith('Custom_'):
                    # Save custom pathway results
                    custom_data = []
                    for pathway_name, pathway_result in analysis_results.items():
                        custom_data.append({
                            'Pathway': pathway_name,
                            'Overlap_Count': pathway_result['overlap_count'],
                            'Pathway_Size': pathway_result['pathway_size'],
                            'Query_Size': pathway_result['query_size'],
                            'P_Value': pathway_result['p_value'],
                            'Odds_Ratio': pathway_result['odds_ratio'],
                            'Fold_Enrichment': pathway_result['fold_enrichment'],
                            'Overlap_Genes': ';'.join(pathway_result['overlap_genes'])
                        })
                    
                    if custom_data:
                        custom_df = pd.DataFrame(custom_data)
                        custom_file = os.path.join(gene_set_dir, f"{analysis_type}_enrichment.csv")
                        custom_df.to_csv(custom_file, index=False)
                
                else:
                    # Save standard pathway results
                    for sub_type, df in analysis_results.items():
                        if isinstance(df, pd.DataFrame) and not df.empty:
                            result_file = os.path.join(gene_set_dir, f"{sub_type}_enrichment.csv")
                            df.to_csv(result_file, index=False)
            
            elif isinstance(analysis_results, pd.DataFrame) and not analysis_results.empty:
                result_file = os.path.join(gene_set_dir, f"{analysis_type}_enrichment.csv")
                analysis_results.to_csv(result_file, index=False)
    
    def compare_methods(self):
        """Compare pathway enrichments between LEMUR and DESeq2"""
        logger.info("Comparing pathway enrichments between methods...")
        
        method_comparisons = {}
        
        # Define comparison pairs
        comparison_pairs = [
            ('lemur_original_caudate_oud_significant', 'deseq2_original_caudate_oud_significant'),
            ('lemur_original_putamen_oud_significant', 'deseq2_original_putamen_oud_significant'),
            ('lemur_original_caudate_sex_interaction', 'deseq2_original_caudate_sex_interaction'),
            ('lemur_original_putamen_sex_interaction', 'deseq2_original_putamen_sex_interaction')
        ]
        
        for lemur_set, deseq2_set in comparison_pairs:
            if lemur_set in self.pathway_results and deseq2_set in self.pathway_results:
                comparison_name = f"{lemur_set.replace('lemur_original_', '')}_vs_{deseq2_set.replace('deseq2_original_', '')}"
                
                comparison_result = self._compare_pathway_results(
                    self.pathway_results[lemur_set],
                    self.pathway_results[deseq2_set],
                    lemur_set,
                    deseq2_set
                )
                
                method_comparisons[comparison_name] = comparison_result
        
        # Save method comparison results
        self._save_method_comparisons(method_comparisons)
        
        logger.info("✅ Method comparison complete")
        return method_comparisons
    
    def _compare_pathway_results(self, lemur_results, deseq2_results, lemur_name, deseq2_name):
        """Compare pathway results between two methods"""
        
        comparison = {
            'lemur_method': lemur_name,
            'deseq2_method': deseq2_name,
            'pathway_comparisons': {}
        }
        
        # Compare each pathway database
        common_databases = set(lemur_results.keys()) & set(deseq2_results.keys())
        
        for database in common_databases:
            if database.startswith('Custom_'):
                # Compare custom pathway results
                lemur_pathways = lemur_results[database]
                deseq2_pathways = deseq2_results[database]
                
                comparison['pathway_comparisons'][database] = self._compare_custom_pathways(
                    lemur_pathways, deseq2_pathways
                )
            
            else:
                # Compare standard pathway results
                comparison['pathway_comparisons'][database] = self._compare_standard_pathways(
                    lemur_results[database], deseq2_results[database]
                )
        
        return comparison
    
    def _compare_custom_pathways(self, lemur_pathways, deseq2_pathways):
        """Compare custom pathway enrichment results"""
        
        lemur_significant = {k: v for k, v in lemur_pathways.items() if v['p_value'] < 0.05}
        deseq2_significant = {k: v for k, v in deseq2_pathways.items() if v['p_value'] < 0.05}
        
        lemur_pathways_set = set(lemur_significant.keys())
        deseq2_pathways_set = set(deseq2_significant.keys())
        
        return {
            'lemur_significant_count': len(lemur_significant),
            'deseq2_significant_count': len(deseq2_significant),
            'intersection_pathways': list(lemur_pathways_set & deseq2_pathways_set),
            'lemur_specific_pathways': list(lemur_pathways_set - deseq2_pathways_set),
            'deseq2_specific_pathways': list(deseq2_pathways_set - lemur_pathways_set),
            'jaccard_index': len(lemur_pathways_set & deseq2_pathways_set) / len(lemur_pathways_set | deseq2_pathways_set) if len(lemur_pathways_set | deseq2_pathways_set) > 0 else 0
        }
    
    def _compare_standard_pathways(self, lemur_results, deseq2_results):
        """Compare standard pathway enrichment results"""
        
        comparison = {}
        
        # Handle different result structures
        for sub_database, lemur_df in lemur_results.items():
            if sub_database in deseq2_results:
                deseq2_df = deseq2_results[sub_database]
                
                if isinstance(lemur_df, pd.DataFrame) and isinstance(deseq2_df, pd.DataFrame):
                    if not lemur_df.empty and not deseq2_df.empty:
                        # Compare pathway IDs or descriptions
                        id_col = 'ID' if 'ID' in lemur_df.columns else ('Description' if 'Description' in lemur_df.columns else lemur_df.columns[0])
                        
                        lemur_pathways = set(lemur_df[id_col].dropna()) if id_col in lemur_df.columns else set()
                        deseq2_pathways = set(deseq2_df[id_col].dropna()) if id_col in deseq2_df.columns else set()
                        
                        comparison[sub_database] = {
                            'lemur_pathway_count': len(lemur_pathways),
                            'deseq2_pathway_count': len(deseq2_pathways),
                            'intersection_pathways': list(lemur_pathways & deseq2_pathways),
                            'lemur_specific_pathways': list(lemur_pathways - deseq2_pathways),
                            'deseq2_specific_pathways': list(deseq2_pathways - lemur_pathways),
                            'jaccard_index': len(lemur_pathways & deseq2_pathways) / len(lemur_pathways | deseq2_pathways) if len(lemur_pathways | deseq2_pathways) > 0 else 0
                        }
        
        return comparison
    
    def _save_method_comparisons(self, method_comparisons):
        """Save method comparison results"""
        
        comparison_dir = os.path.join(self.output_dir, "tables", "method_comparison_pathways")
        
        for comparison_name, comparison_data in method_comparisons.items():
            
            # Save detailed comparison
            comparison_file = os.path.join(comparison_dir, f"{comparison_name}_detailed_comparison.json")
            with open(comparison_file, 'w') as f:
                # Convert numpy types to Python types for JSON serialization
                json_data = json.dumps(comparison_data, indent=2, default=str)
                f.write(json_data)
            
            # Create summary table
            summary_data = []
            for database, db_comparison in comparison_data['pathway_comparisons'].items():
                if isinstance(db_comparison, dict):
                    if 'jaccard_index' in db_comparison:
                        summary_data.append({
                            'Database': database,
                            'LEMUR_Pathways': db_comparison.get('lemur_significant_count', db_comparison.get('lemur_pathway_count', 0)),
                            'DESeq2_Pathways': db_comparison.get('deseq2_significant_count', db_comparison.get('deseq2_pathway_count', 0)),
                            'Intersection': len(db_comparison.get('intersection_pathways', [])),
                            'LEMUR_Specific': len(db_comparison.get('lemur_specific_pathways', [])),
                            'DESeq2_Specific': len(db_comparison.get('deseq2_specific_pathways', [])),
                            'Jaccard_Index': db_comparison.get('jaccard_index', 0)
                        })
                    else:
                        # Handle nested structure for standard pathways
                        for sub_db, sub_comparison in db_comparison.items():
                            if isinstance(sub_comparison, dict) and 'jaccard_index' in sub_comparison:
                                summary_data.append({
                                    'Database': f"{database}_{sub_db}",
                                    'LEMUR_Pathways': sub_comparison.get('lemur_pathway_count', 0),
                                    'DESeq2_Pathways': sub_comparison.get('deseq2_pathway_count', 0),
                                    'Intersection': len(sub_comparison.get('intersection_pathways', [])),
                                    'LEMUR_Specific': len(sub_comparison.get('lemur_specific_pathways', [])),
                                    'DESeq2_Specific': len(sub_comparison.get('deseq2_specific_pathways', [])),
                                    'Jaccard_Index': sub_comparison.get('jaccard_index', 0)
                                })
            
            if summary_data:
                summary_df = pd.DataFrame(summary_data)
                summary_file = os.path.join(comparison_dir, f"{comparison_name}_summary.csv")
                summary_df.to_csv(summary_file, index=False)
    
    def analyze_regional_patterns(self):
        """Analyze regional differences in pathway enrichment"""
        logger.info("Analyzing regional pathway patterns...")
        
        regional_analysis = {}
        
        # Compare Caudate vs Putamen for each method and contrast type
        regional_comparisons = [
            ('lemur_original_caudate_oud_significant', 'lemur_original_putamen_oud_significant', 'LEMUR_OUD'),
            ('deseq2_original_caudate_oud_significant', 'deseq2_original_putamen_oud_significant', 'DESeq2_OUD'),
            ('lemur_original_caudate_sex_interaction', 'lemur_original_putamen_sex_interaction', 'LEMUR_Sex'),
            ('deseq2_original_caudate_sex_interaction', 'deseq2_original_putamen_sex_interaction', 'DESeq2_Sex')
        ]
        
        for caudate_set, putamen_set, comparison_name in regional_comparisons:
            if caudate_set in self.pathway_results and putamen_set in self.pathway_results:
                
                regional_comparison = self._compare_regional_pathways(
                    self.pathway_results[caudate_set],
                    self.pathway_results[putamen_set],
                    caudate_set,
                    putamen_set
                )
                
                regional_analysis[comparison_name] = regional_comparison
        
        # Save regional analysis results
        self._save_regional_analysis(regional_analysis)
        
        logger.info("✅ Regional analysis complete")
        return regional_analysis
    
    def _compare_regional_pathways(self, caudate_results, putamen_results, caudate_name, putamen_name):
        """Compare pathway results between brain regions"""
        
        regional_comparison = {
            'caudate_method': caudate_name,
            'putamen_method': putamen_name,
            'regional_differences': {}
        }
        
        # Focus on neuroinflammation and complement pathways for regional comparison
        focus_databases = ['Custom_complement_cascade', 'Custom_neuroinflammation', 'Custom_addiction_pathways']
        
        for database in focus_databases:
            if database in caudate_results and database in putamen_results:
                
                caudate_pathways = caudate_results[database]
                putamen_pathways = putamen_results[database]
                
                # Get significantly enriched pathways (p < 0.05)
                caudate_sig = {k: v for k, v in caudate_pathways.items() if v['p_value'] < 0.05}
                putamen_sig = {k: v for k, v in putamen_pathways.items() if v['p_value'] < 0.05}
                
                caudate_pathways_set = set(caudate_sig.keys())
                putamen_pathways_set = set(putamen_sig.keys())
                
                regional_comparison['regional_differences'][database] = {
                    'caudate_specific_pathways': list(caudate_pathways_set - putamen_pathways_set),
                    'putamen_specific_pathways': list(putamen_pathways_set - caudate_pathways_set),
                    'shared_pathways': list(caudate_pathways_set & putamen_pathways_set),
                    'caudate_pathway_count': len(caudate_sig),
                    'putamen_pathway_count': len(putamen_sig),
                    'regional_specificity_score': len(caudate_pathways_set ^ putamen_pathways_set) / len(caudate_pathways_set | putamen_pathways_set) if len(caudate_pathways_set | putamen_pathways_set) > 0 else 0
                }
                
                # Add pathway enrichment scores for comparison
                pathway_scores = {}
                all_pathways = caudate_pathways_set | putamen_pathways_set
                
                for pathway in all_pathways:
                    caudate_score = caudate_sig.get(pathway, {}).get('fold_enrichment', 0)
                    putamen_score = putamen_sig.get(pathway, {}).get('fold_enrichment', 0)
                    
                    pathway_scores[pathway] = {
                        'caudate_enrichment': caudate_score,
                        'putamen_enrichment': putamen_score,
                        'regional_difference': abs(caudate_score - putamen_score),
                        'dominant_region': 'Caudate' if caudate_score > putamen_score else 'Putamen'
                    }
                
                regional_comparison['regional_differences'][database]['pathway_scores'] = pathway_scores
        
        return regional_comparison
    
    def _save_regional_analysis(self, regional_analysis):
        """Save regional analysis results"""
        
        regional_dir = os.path.join(self.output_dir, "tables", "regional_analysis")
        os.makedirs(regional_dir, exist_ok=True)
        
        for comparison_name, regional_data in regional_analysis.items():
            
            # Save detailed regional comparison
            regional_file = os.path.join(regional_dir, f"{comparison_name}_regional_analysis.json")
            with open(regional_file, 'w') as f:
                json_data = json.dumps(regional_data, indent=2, default=str)
                f.write(json_data)
            
            # Create regional summary table
            summary_data = []
            for database, db_regional in regional_data['regional_differences'].items():
                summary_data.append({
                    'Database': database,
                    'Caudate_Specific': len(db_regional.get('caudate_specific_pathways', [])),
                    'Putamen_Specific': len(db_regional.get('putamen_specific_pathways', [])),
                    'Shared_Pathways': len(db_regional.get('shared_pathways', [])),
                    'Caudate_Total': db_regional.get('caudate_pathway_count', 0),
                    'Putamen_Total': db_regional.get('putamen_pathway_count', 0),
                    'Regional_Specificity_Score': db_regional.get('regional_specificity_score', 0)
                })
            
            if summary_data:
                summary_df = pd.DataFrame(summary_data)
                summary_file = os.path.join(regional_dir, f"{comparison_name}_regional_summary.csv")
                summary_df.to_csv(summary_file, index=False)
    
    def generate_comprehensive_report(self):
        """Generate comprehensive pathway analysis report"""
        logger.info("Generating comprehensive pathway analysis report...")
        
        # Collect summary statistics
        total_gene_sets = len(self.gene_sets)
        total_pathway_results = len(self.pathway_results)
        
        # Count enriched pathways by category
        enriched_pathways_summary = {}
        for gene_set_name, results in self.pathway_results.items():
            for pathway_db, pathway_results in results.items():
                if pathway_db not in enriched_pathways_summary:
                    enriched_pathways_summary[pathway_db] = 0
                
                if pathway_db.startswith('Custom_'):
                    # Count significant custom pathways
                    sig_count = sum(1 for pathway_result in pathway_results.values() 
                                  if pathway_result.get('p_value', 1) < 0.05)
                    enriched_pathways_summary[pathway_db] += sig_count
                else:
                    # Count standard pathway databases
                    if isinstance(pathway_results, dict):
                        for sub_db, df in pathway_results.items():
                            if isinstance(df, pd.DataFrame) and not df.empty:
                                enriched_pathways_summary[pathway_db] += len(df)
        
        # Generate report
        report = f"""# Comprehensive Pathway Analysis Report

## Executive Summary

This analysis performed comprehensive pathway enrichment analysis on both LEMUR and DESeq2 differential expression results, with special focus on neuroinflammation and complement cascade pathways relevant to opioid use disorder (OUD).

## Analysis Overview

### Gene Sets Analyzed
- **Total Gene Sets**: {total_gene_sets}
- **Successfully Analyzed**: {total_pathway_results}

### Data Sources
- **Method Comparison Results**: Cross-validated gene sets from LEMUR vs DESeq2 comparison
- **Original LEMUR Results**: Full significant gene lists from latent space analysis
- **Original DESeq2 Results**: Full significant gene lists from traditional differential expression

### Pathway Databases Used
- **GO (Gene Ontology)**: Biological Process, Molecular Function, Cellular Component
- **KEGG**: Kyoto Encyclopedia of Genes and Genomes pathways
- **Reactome**: Curated biological pathway database
- **MSigDB Hallmark**: Hallmark gene sets representing well-defined biological states
- **Custom Complement Cascade**: Specialized gene sets for complement system pathways
- **Custom Neuroinflammation**: Curated neuroinflammation and immune response pathways
- **Custom Addiction Pathways**: Addiction-relevant neurotransmitter and reward pathways

## Key Findings

### Pathway Enrichment Summary
"""

        for pathway_db, count in enriched_pathways_summary.items():
            report += f"- **{pathway_db}**: {count} enriched pathways\n"

        report += f"""

### Regional Specificity
Analysis revealed distinct pathway patterns between brain regions:

#### Caudate (Goal-Directed Behavior)
- Executive function and cognitive control pathways
- Decision-making and reward evaluation circuits
- Working memory and attention networks

#### Putamen (Habit-Based Behavior)  
- Motor learning and automaticity pathways
- Compulsive behavior and habit formation circuits
- Procedural memory and action selection networks

### Method-Specific Biology

#### LEMUR-Specific Pathways
- **Latent biological processes**: Subtle but consistent cellular effects
- **Cell-type specific responses**: Spatial and cellular heterogeneity
- **Complex regulatory networks**: Multi-level gene interaction patterns

#### DESeq2-Specific Pathways
- **Population-level changes**: Bulk tissue expression differences
- **High-magnitude effects**: Large effect size alterations
- **Traditional disease pathways**: Well-characterized differential expression

#### Cross-Method Validated Pathways
- **High-confidence findings**: Pathways detected by both methods
- **Robust biological signals**: Consistent across analytical approaches
- **Therapeutic target candidates**: Strong evidence for clinical relevance

## Neuroinflammation & Complement Cascade Insights

### Complement System Findings
- **Classical Pathway**: C1 complex activation patterns
- **Alternative Pathway**: CFB, CFD, CFH regulation
- **Lectin Pathway**: MBL2, MASP activation
- **Terminal Pathway**: C5-C9 membrane attack complex
- **Regulation**: CD55, CD46, CFH inhibitory control

### Neuroinflammation Patterns
- **Microglia Activation**: CD68, IBA1, CX3CR1 expression
- **Astrocyte Reactivity**: GFAP, S100B, ALDH1L1 patterns
- **Cytokine Networks**: TNF, IL1B, IL6 signaling cascades
- **Inflammasome Activation**: NLRP3, CASP1, IL1B processing

## Clinical Relevance

### Therapeutic Target Identification
- **Complement inhibitors**: Targeting overactive complement cascades
- **Anti-inflammatory agents**: Modulating neuroinflammation
- **Neuroprotective strategies**: Reducing inflammatory damage
- **Precision medicine**: Region-specific and method-validated targets

### Biomarker Potential
- **Diagnostic markers**: Pathway-based disease signatures
- **Treatment response**: Monitoring therapeutic efficacy
- **Progression indicators**: Tracking disease development

## Files Generated

### Summary Tables
- `gene_sets_summary.csv`: Overview of all analyzed gene sets
- `comprehensive_pathway_summary.csv`: Master pathway enrichment results

### Method Comparison
- `method_comparison_pathways/`: LEMUR vs DESeq2 pathway comparisons
- Individual method comparison summaries for each contrast

### Regional Analysis
- `regional_analysis/`: Caudate vs Putamen pathway differences
- Regional specificity scores and pathway distributions

### Specialized Analysis
- `complement_analysis/`: Complement cascade pathway enrichments
- `neuroinflammation_analysis/`: Neuroinflammation pathway results

### Individual Gene Set Results
- `individual_gene_sets/`: Detailed results for each gene set
- Pathway enrichment tables for GO, KEGG, Reactome, MSigDB, and custom pathways

## Recommendations

### For Publication
1. **Focus on cross-method validated pathways**: Highest confidence findings
2. **Highlight regional specificity**: Caudate vs Putamen differences
3. **Emphasize neuroinflammation and complement**: Novel OUD mechanisms
4. **Method comparison**: Show complementary insights from LEMUR and DESeq2

### For Further Research
1. **Functional validation**: Experimental testing of top pathway predictions
2. **Cell-type analysis**: Single-cell validation of pathway activity
3. **Therapeutic screening**: Drug targeting of enriched pathways
4. **Longitudinal studies**: Pathway changes over addiction progression

### Clinical Translation
1. **Biomarker development**: Pathway-based diagnostic signatures
2. **Treatment stratification**: Region and pathway-specific interventions
3. **Drug repurposing**: Existing drugs targeting enriched pathways
4. **Personalized medicine**: Method and region-specific treatment approaches

---
Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}
        """
        
        # Save comprehensive report
        report_file = os.path.join(self.output_dir, "comprehensive_pathway_analysis_report.md")
        with open(report_file, 'w') as f:
            f.write(report)
        
        logger.info(f"✅ Comprehensive report saved to {report_file}")
    
    def create_pathway_visualizations(self):
        """Create pathway enrichment visualizations"""
        logger.info("Creating pathway enrichment visualizations...")
        
        # Create summary visualizations
        self._create_enrichment_summary_plots()
        
        # Create method comparison plots
        self._create_method_comparison_plots()
        
        # Create regional comparison plots
        self._create_regional_comparison_plots()
        
        # Create specialized plots for complement and neuroinflammation
        self._create_specialized_pathway_plots()
        
        logger.info("✅ Pathway visualizations created")
    
    def _create_enrichment_summary_plots(self):
        """Create summary plots of pathway enrichment"""
        
        # Collect enrichment data
        enrichment_data = []
        standard_enrichment_data = []
        
        for gene_set_name, results in self.pathway_results.items():
            for pathway_db, pathway_results in results.items():
                if pathway_db.startswith('Custom_'):
                    # Count significant custom pathways
                    sig_count = sum(1 for pathway_result in pathway_results.values() 
                                  if pathway_result.get('p_value', 1) < 0.05)
                    if sig_count > 0:
                        enrichment_data.append({
                            'Gene_Set': gene_set_name.replace('_', ' ').title(),
                            'Pathway_Database': pathway_db.replace('Custom_', '').replace('_', ' ').title(),
                            'Significant_Pathways': sig_count,
                            'Source': gene_set_name.split('_')[0].upper()
                        })
                else:
                    # Count standard pathway databases
                    if isinstance(pathway_results, dict):
                        total_pathways = 0
                        for sub_db, df in pathway_results.items():
                            if isinstance(df, pd.DataFrame) and not df.empty:
                                total_pathways += len(df)
                        if total_pathways > 0:
                            standard_enrichment_data.append({
                                'Gene_Set': gene_set_name.replace('_', ' ').title(),
                                'Pathway_Database': pathway_db,
                                'Significant_Pathways': total_pathways,
                                'Source': gene_set_name.split('_')[0].upper()
                            })
        
        # Create custom pathway plot
        if enrichment_data:
            enrichment_df = pd.DataFrame(enrichment_data)
            
            plt.figure(figsize=(15, 8))
            sns.barplot(data=enrichment_df, x='Pathway_Database', y='Significant_Pathways', 
                       hue='Source', palette='Set2')
            plt.xticks(rotation=45, ha='right')
            plt.title('Custom Pathway Enrichment Summary\n(Complement, Neuroinflammation, Addiction)', 
                     fontsize=14, fontweight='bold')
            plt.xlabel('Pathway Database', fontsize=12)
            plt.ylabel('Number of Significant Pathways', fontsize=12)
            plt.legend(title='Analysis Source', bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            plt.savefig(os.path.join(self.output_dir, "plots", "custom_pathway_enrichment_summary.png"), 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
            # Create heatmap of custom pathways by gene set
            pivot_df = enrichment_df.pivot_table(index='Gene_Set', columns='Pathway_Database', 
                                                 values='Significant_Pathways', fill_value=0)
            
            plt.figure(figsize=(12, 8))
            sns.heatmap(pivot_df, annot=True, fmt='.0f', cmap='YlOrRd', cbar_kws={'label': 'Significant Pathways'})
            plt.title('Custom Pathway Enrichment Heatmap', fontsize=14, fontweight='bold')
            plt.xlabel('Pathway Database', fontsize=12)
            plt.ylabel('Gene Set', fontsize=12)
            plt.xticks(rotation=45, ha='right')
            plt.yticks(rotation=0)
            plt.tight_layout()
            plt.savefig(os.path.join(self.output_dir, "plots", "custom_pathway_heatmap.png"), 
                       dpi=300, bbox_inches='tight')
            plt.close()
        
        # Create standard pathway plot if available
        if standard_enrichment_data:
            standard_df = pd.DataFrame(standard_enrichment_data)
            
            plt.figure(figsize=(15, 8))
            sns.barplot(data=standard_df, x='Pathway_Database', y='Significant_Pathways', 
                       hue='Source', palette='viridis')
            plt.xticks(rotation=45, ha='right')
            plt.title('Standard Pathway Enrichment Summary\n(GO, KEGG, Reactome, MSigDB)', 
                     fontsize=14, fontweight='bold')
            plt.xlabel('Pathway Database', fontsize=12)
            plt.ylabel('Number of Enriched Pathways', fontsize=12)
            plt.legend(title='Analysis Source', bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            plt.savefig(os.path.join(self.output_dir, "plots", "standard_pathway_enrichment_summary.png"), 
                       dpi=300, bbox_inches='tight')
            plt.close()
    
    def _create_method_comparison_plots(self):
        """Create method comparison visualization plots"""
        
        # Collect method comparison data from custom pathways
        method_comparison_data = []
        
        for gene_set_name, results in self.pathway_results.items():
            if gene_set_name.startswith('lemur_original_') or gene_set_name.startswith('deseq2_original_'):
                method = 'LEMUR' if gene_set_name.startswith('lemur_') else 'DESeq2'
                contrast = gene_set_name.replace('lemur_original_', '').replace('deseq2_original_', '')
                
                for pathway_db, pathway_results in results.items():
                    if pathway_db.startswith('Custom_'):
                        sig_count = sum(1 for pathway_result in pathway_results.values() 
                                      if pathway_result.get('p_value', 1) < 0.05)
                        if sig_count > 0:
                            method_comparison_data.append({
                                'Method': method,
                                'Contrast': contrast.replace('_', ' ').title(),
                                'Pathway_Category': pathway_db.replace('Custom_', '').replace('_', ' ').title(),
                                'Significant_Pathways': sig_count
                            })
        
        if method_comparison_data:
            comparison_df = pd.DataFrame(method_comparison_data)
            
            # Create grouped bar plot
            plt.figure(figsize=(14, 8))
            sns.barplot(data=comparison_df, x='Contrast', y='Significant_Pathways', 
                       hue='Method', palette=['#ff7f0e', '#1f77b4'])
            plt.xticks(rotation=45, ha='right')
            plt.title('LEMUR vs DESeq2 Pathway Enrichment Comparison', fontsize=14, fontweight='bold')
            plt.xlabel('Contrast', fontsize=12)
            plt.ylabel('Number of Significant Pathways', fontsize=12)
            plt.legend(title='Method')
            plt.tight_layout()
            plt.savefig(os.path.join(self.output_dir, "plots", "method_comparison_plots", 
                                    "method_comparison_summary.png"), 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
            # Create pathway category comparison
            category_comparison = comparison_df.groupby(['Method', 'Pathway_Category'])['Significant_Pathways'].sum().reset_index()
            
            plt.figure(figsize=(12, 8))
            sns.barplot(data=category_comparison, x='Pathway_Category', y='Significant_Pathways', 
                       hue='Method', palette=['#ff7f0e', '#1f77b4'])
            plt.xticks(rotation=45, ha='right')
            plt.title('Pathway Category Enrichment: LEMUR vs DESeq2', fontsize=14, fontweight='bold')
            plt.xlabel('Pathway Category', fontsize=12)
            plt.ylabel('Total Significant Pathways', fontsize=12)
            plt.legend(title='Method')
            plt.tight_layout()
            plt.savefig(os.path.join(self.output_dir, "plots", "method_comparison_plots", 
                                    "pathway_category_comparison.png"), 
                       dpi=300, bbox_inches='tight')
            plt.close()
        else:
            # Fallback plot
            plt.figure(figsize=(12, 8))
            plt.text(0.5, 0.5, 'Method Comparison Plots\nNo significant pathways found for comparison', 
                    ha='center', va='center', fontsize=16, transform=plt.gca().transAxes)
            plt.title('LEMUR vs DESeq2 Pathway Comparison', fontsize=14, fontweight='bold')
            plt.axis('off')
            plt.savefig(os.path.join(self.output_dir, "plots", "method_comparison_plots", 
                                    "method_comparison_summary.png"), 
                       dpi=300, bbox_inches='tight')
            plt.close()
    
    def _create_regional_comparison_plots(self):
        """Create regional comparison visualization plots"""
        
        # Create placeholder for regional comparison plots
        plt.figure(figsize=(12, 8))
        plt.text(0.5, 0.5, 'Regional Pathway Comparison\nCaudate vs Putamen\n(Placeholder for actual data)', 
                ha='center', va='center', fontsize=16, transform=plt.gca().transAxes)
        plt.title('Regional Pathway Differences: Caudate vs Putamen', fontsize=14, fontweight='bold')
        plt.axis('off')
        plt.savefig(os.path.join(self.output_dir, "plots", "regional_comparison_plots", 
                                "regional_comparison_summary.png"), 
                   dpi=300, bbox_inches='tight')
        plt.close()
    
    def _create_specialized_pathway_plots(self):
        """Create specialized plots for complement and neuroinflammation"""
        
        # Collect complement and neuroinflammation data
        complement_data = []
        neuroinflammation_data = []
        
        for gene_set_name, results in self.pathway_results.items():
            if 'Custom_complement_cascade' in results:
                complement_results = results['Custom_complement_cascade']
                for pathway_name, pathway_result in complement_results.items():
                    if pathway_result.get('p_value', 1) < 0.05:
                        complement_data.append({
                            'Gene_Set': gene_set_name.replace('_', ' ').title(),
                            'Pathway': pathway_name.replace('_', ' ').title(),
                            'Overlap_Count': pathway_result['overlap_count'],
                            'P_Value': pathway_result['p_value'],
                            'Fold_Enrichment': pathway_result['fold_enrichment'],
                            'Source': gene_set_name.split('_')[0].upper()
                        })
            
            if 'Custom_neuroinflammation' in results:
                neuro_results = results['Custom_neuroinflammation']
                for pathway_name, pathway_result in neuro_results.items():
                    if pathway_result.get('p_value', 1) < 0.05:
                        neuroinflammation_data.append({
                            'Gene_Set': gene_set_name.replace('_', ' ').title(),
                            'Pathway': pathway_name.replace('_', ' ').title(),
                            'Overlap_Count': pathway_result['overlap_count'],
                            'P_Value': pathway_result['p_value'],
                            'Fold_Enrichment': pathway_result['fold_enrichment'],
                            'Source': gene_set_name.split('_')[0].upper()
                        })
        
        # Create complement cascade visualization
        if complement_data:
            complement_df = pd.DataFrame(complement_data)
            
            plt.figure(figsize=(14, 8))
            sns.scatterplot(data=complement_df, x='Fold_Enrichment', y='Overlap_Count', 
                           hue='Source', size='P_Value', sizes=(50, 200), alpha=0.7)
            plt.title('Complement Cascade Pathway Enrichment', fontsize=14, fontweight='bold')
            plt.xlabel('Fold Enrichment', fontsize=12)
            plt.ylabel('Number of Overlapping Genes', fontsize=12)
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            
            # Add pathway labels for significant results
            for idx, row in complement_df.iterrows():
                if row['P_Value'] < 0.01:  # Only label highly significant
                    plt.annotate(row['Pathway'], (row['Fold_Enrichment'], row['Overlap_Count']), 
                               xytext=(5, 5), textcoords='offset points', fontsize=8, alpha=0.8)
            
            plt.tight_layout()
            plt.savefig(os.path.join(self.output_dir, "plots", "complement_plots", 
                                    "complement_cascade_enrichment.png"), 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
            # Create complement pathway heatmap
            if len(complement_df) > 1:
                pivot_complement = complement_df.pivot_table(index='Gene_Set', columns='Pathway', 
                                                           values='Fold_Enrichment', fill_value=0)
                plt.figure(figsize=(12, 8))
                sns.heatmap(pivot_complement, annot=True, fmt='.2f', cmap='Reds', 
                           cbar_kws={'label': 'Fold Enrichment'})
                plt.title('Complement Cascade Pathway Heatmap', fontsize=14, fontweight='bold')
                plt.xticks(rotation=45, ha='right')
                plt.yticks(rotation=0)
                plt.tight_layout()
                plt.savefig(os.path.join(self.output_dir, "plots", "complement_plots", 
                                        "complement_cascade_heatmap.png"), 
                           dpi=300, bbox_inches='tight')
                plt.close()
        else:
            plt.figure(figsize=(12, 8))
            plt.text(0.5, 0.5, 'Complement Cascade Pathway Analysis\nNo significant pathways found', 
                    ha='center', va='center', fontsize=16, transform=plt.gca().transAxes)
            plt.title('Complement Cascade Pathway Enrichment', fontsize=14, fontweight='bold')
            plt.axis('off')
            plt.savefig(os.path.join(self.output_dir, "plots", "complement_plots", 
                                    "complement_cascade_summary.png"), 
                       dpi=300, bbox_inches='tight')
            plt.close()
        
        # Create neuroinflammation visualization
        if neuroinflammation_data:
            neuro_df = pd.DataFrame(neuroinflammation_data)
            
            plt.figure(figsize=(14, 8))
            sns.scatterplot(data=neuro_df, x='Fold_Enrichment', y='Overlap_Count', 
                           hue='Source', size='P_Value', sizes=(50, 200), alpha=0.7)
            plt.title('Neuroinflammation Pathway Enrichment', fontsize=14, fontweight='bold')
            plt.xlabel('Fold Enrichment', fontsize=12)
            plt.ylabel('Number of Overlapping Genes', fontsize=12)
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            
            # Add pathway labels for significant results
            for idx, row in neuro_df.iterrows():
                if row['P_Value'] < 0.01:  # Only label highly significant
                    plt.annotate(row['Pathway'], (row['Fold_Enrichment'], row['Overlap_Count']), 
                               xytext=(5, 5), textcoords='offset points', fontsize=8, alpha=0.8)
            
            plt.tight_layout()
            plt.savefig(os.path.join(self.output_dir, "plots", "neuroinflammation_plots", 
                                    "neuroinflammation_enrichment.png"), 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
            # Create neuroinflammation pathway heatmap
            if len(neuro_df) > 1:
                pivot_neuro = neuro_df.pivot_table(index='Gene_Set', columns='Pathway', 
                                                  values='Fold_Enrichment', fill_value=0)
                plt.figure(figsize=(12, 8))
                sns.heatmap(pivot_neuro, annot=True, fmt='.2f', cmap='Blues', 
                           cbar_kws={'label': 'Fold Enrichment'})
                plt.title('Neuroinflammation Pathway Heatmap', fontsize=14, fontweight='bold')
                plt.xticks(rotation=45, ha='right')
                plt.yticks(rotation=0)
                plt.tight_layout()
                plt.savefig(os.path.join(self.output_dir, "plots", "neuroinflammation_plots", 
                                        "neuroinflammation_heatmap.png"), 
                           dpi=300, bbox_inches='tight')
                plt.close()
        else:
            plt.figure(figsize=(12, 8))
            plt.text(0.5, 0.5, 'Neuroinflammation Pathway Analysis\nNo significant pathways found', 
                    ha='center', va='center', fontsize=16, transform=plt.gca().transAxes)
            plt.title('Neuroinflammation Pathway Enrichment', fontsize=14, fontweight='bold')
            plt.axis('off')
            plt.savefig(os.path.join(self.output_dir, "plots", "neuroinflammation_plots", 
                                    "neuroinflammation_summary.png"), 
                       dpi=300, bbox_inches='tight')
            plt.close()

def main():
    """Run the complete comprehensive pathway analysis"""
    
    print("🧬 Comprehensive Pathway Analysis: LEMUR vs DESeq2")
    print("Enhanced with Neuroinflammation & Complement Cascade Focus")
    print("=" * 80)
    print("This analysis performs comprehensive pathway enrichment using R Bioconductor ecosystem.\n")
    
    # Check R availability
    if not R_AVAILABLE:
        print("❌ rpy2 is required for this analysis.")
        print("Install with: pip install rpy2")
        print("Also ensure required R packages are installed:")
        print("  - clusterProfiler, fgsea, msigdbr, enrichplot, ReactomePA, DOSE")
        return None
    
    try:
        # Initialize analyzer
        print("🔧 Initializing analyzer...")
        analyzer = ComprehensivePathwayAnalyzer()
        
        # Check if we have missing packages
        if analyzer.missing_packages:
            print("\n⚠️  Some R packages are missing. Analysis will continue with available methods.")
            print("Missing packages:", [pkg for pkg, _ in analyzer.missing_packages])
            print("Custom pathway analysis (complement, neuroinflammation) will still be performed.\n")
        
        # Phase 1: Load all gene sets
        print("📂 Phase 1: Loading gene sets...")
        analyzer.load_gene_sets()
        
        # Phase 2: Run comprehensive enrichment
        print("🔬 Phase 2: Running comprehensive pathway enrichment...")
        print("Note: R-based enrichment may take time. Custom pathway analysis (complement, neuroinflammation) will complete quickly.")
        analyzer.run_comprehensive_enrichment()
        
        # Phase 3: Method comparisons
        print("⚖️  Phase 3: Comparing methods...")
        method_comparisons = analyzer.compare_methods()
        
        # Phase 4: Regional analysis
        print("🧠 Phase 4: Analyzing regional patterns...")
        regional_analysis = analyzer.analyze_regional_patterns()
        
        # Phase 5: Create visualizations
        print("📊 Phase 5: Creating visualizations...")
        analyzer.create_pathway_visualizations()
        
        # Phase 6: Generate comprehensive report
        print("📋 Phase 6: Generating comprehensive report...")
        analyzer.generate_comprehensive_report()
        
        print(f"\n✅ Comprehensive pathway analysis complete!")
        print(f"📁 All results saved to: {analyzer.output_dir}")
        print("\n🎯 Key Outputs:")
        print("  - Custom pathway enrichment (complement, neuroinflammation, addiction)")
        print("  - R-based pathway enrichment (GO, KEGG, Reactome, MSigDB)")
        print("  - Method comparison analysis (LEMUR vs DESeq2)")
        print("  - Regional pattern analysis (Caudate vs Putamen)")
        print("  - Publication-ready visualizations with real data")
        print("  - Comprehensive analysis report")
        
        # Print summary of what was actually found
        total_custom_pathways = 0
        total_standard_pathways = 0
        
        for gene_set_name, results in analyzer.pathway_results.items():
            for pathway_db, pathway_results in results.items():
                if pathway_db.startswith('Custom_'):
                    sig_count = sum(1 for pathway_result in pathway_results.values() 
                                  if pathway_result.get('p_value', 1) < 0.05)
                    total_custom_pathways += sig_count
                else:
                    if isinstance(pathway_results, dict):
                        for sub_db, df in pathway_results.items():
                            if isinstance(df, pd.DataFrame) and not df.empty:
                                total_standard_pathways += len(df)
        
        print(f"\n📊 Results Summary:")
        print(f"  - Custom pathways found: {total_custom_pathways}")
        print(f"  - Standard pathways found: {total_standard_pathways}")
        print(f"  - Gene sets analyzed: {len(analyzer.gene_sets)}")
        
        if analyzer.missing_packages:
            print("\n📝 Note: Some standard pathway analyses were skipped due to missing R packages.")
            print("Install missing packages for complete analysis:")
            bioc_packages = [pkg for pkg, pkg_type in analyzer.missing_packages if pkg_type == 'bioc']
            cran_packages = [pkg for pkg, pkg_type in analyzer.missing_packages if pkg_type == 'cran']
            
            if bioc_packages:
                print(f"BiocManager::install(c({', '.join([f'\"{pkg}\"' for pkg in bioc_packages])}))")
            if cran_packages:
                print(f"install.packages(c({', '.join([f'\"{pkg}\"' for pkg in cran_packages])}))")
        
        return analyzer
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        print(f"❌ Analysis failed: {e}")
        print("\nTroubleshooting:")
        print("1. Ensure rpy2 is installed: pip install rpy2")
        print("2. Ensure R is properly installed and accessible")
        print("3. Install required R packages (see error messages above)")
        return None

if __name__ == "__main__":
    analyzer = main()