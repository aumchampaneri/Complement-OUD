# ========================================================================
# GSE118918 Pathway Enrichment Analysis: Complement and Morphine Response
# ========================================================================
# 
# STUDY OVERVIEW:
# Dataset: GSE118918 - Nucleus Accumbens RNA-seq (Mock vs Morphine)
# Analysis: Comprehensive pathway enrichment with focus on complement system
# Integration: Morphine response, addiction pathways, and complement cascade
# 
# SCIENTIFIC RATIONALE:
# This analysis follows established best practices for pathway enrichment:
# - Yu et al. (2012) clusterProfiler: universal enrichment tool
# - Subramanian et al. (2005) GSEA: Gene Set Enrichment Analysis
# - Liberzon et al. (2011) MSigDB: molecular signatures database
# - Hänzelmann et al. (2013) GSVA: gene set variation analysis
# 
# PATHWAY DATABASES:
# - Gene Ontology (GO): Biological processes, molecular functions, cellular components
# - KEGG: Kyoto Encyclopedia of Genes and Genomes
# - Reactome: Curated biological pathways
# - MSigDB: Molecular Signatures Database
# - Custom complement pathway definitions
# 
# COMPLEMENT SYSTEM FOCUS:
# - Classical, alternative, and lectin pathways
# - Complement regulators and receptors
# - Integration with neuroinflammation and addiction pathways
# - Morphine-induced complement activation
# ========================================================================

# Set reproducibility parameters
set.seed(42)
options(stringsAsFactors = FALSE)

# Record analysis start time
analysis_start_time <- Sys.time()
cat("Pathway Enrichment Analysis started at:", as.character(analysis_start_time), "\n")

# ========================================================================
# SECTION 1: LIBRARY LOADING AND ENVIRONMENT SETUP
# ========================================================================

cat("\n=== LOADING REQUIRED LIBRARIES ===\n")

# Core enrichment analysis packages
required_packages <- c(
  "clusterProfiler",   # GO and KEGG enrichment
  "enrichplot",        # Enrichment visualization
  "DOSE",             # Disease ontology
  "ReactomePA",       # Reactome pathway analysis
  "msigdbr",          # MSigDB access
  "GSVA",             # Gene set variation analysis
  "GSEABase",         # Gene set utilities
  "org.Mm.eg.db",     # Mouse organism database
  "AnnotationDbi",    # Annotation utilities
  "edgeR",            # For DE results loading
  "limma",            # For statistical methods
  "dplyr",            # Data manipulation
  "ggplot2",          # Plotting
  "RColorBrewer",     # Color schemes
  "ComplexHeatmap",   # Advanced heatmaps
  "circlize",         # Circular plots
  "VennDiagram",      # Venn diagrams
  "UpSetR"            # UpSet plots
)

# Enhanced optional packages
optional_packages <- c(
  "pathview",         # KEGG pathway visualization
  "Rgraphviz",        # Graph visualization
  "igraph",           # Network analysis
  "networkD3",        # Interactive networks
  "plotly",           # Interactive plots
  "DT",               # Interactive tables
  "rmarkdown",        # Report generation
  "knitr"             # Document processing
)

# Load required packages
loaded_packages <- character()
for(pkg in required_packages) {
  if(require(pkg, character.only = TRUE, quietly = TRUE)) {
    loaded_packages <- c(loaded_packages, pkg)
  } else {
    warning(paste("Required package", pkg, "not available. Some analyses may be limited."))
  }
}

# Load optional packages
loaded_optional <- character()
for(pkg in optional_packages) {
  if(require(pkg, character.only = TRUE, quietly = TRUE)) {
    loaded_optional <- c(loaded_optional, pkg)
  }
}

cat("Required packages loaded:", length(loaded_packages), "/", length(required_packages), "\n")
cat("Optional packages available:", paste(loaded_optional, collapse = ", "), "\n")

# ========================================================================
# SECTION 2: DATA LOADING AND SETUP
# ========================================================================

cat("\n=== LOADING DIFFERENTIAL EXPRESSION RESULTS ===\n")

# Define paths
base_dir <- "/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus/GSE118918"
setwd(base_dir)

# Create output directory structure
output_structure <- list(
  main = "Outputs/02_Pathway_Enrichment",
  data = "Outputs/02_Pathway_Enrichment/Data",
  plots = "Outputs/02_Pathway_Enrichment/Plots",
  reports = "Outputs/02_Pathway_Enrichment/Reports",
  tables = "Outputs/02_Pathway_Enrichment/Tables",
  go_analysis = "Outputs/02_Pathway_Enrichment/GO_Analysis",
  kegg_analysis = "Outputs/02_Pathway_Enrichment/KEGG_Analysis",
  complement = "Outputs/02_Pathway_Enrichment/Complement_Analysis",
  gsea = "Outputs/02_Pathway_Enrichment/GSEA",
  networks = "Outputs/02_Pathway_Enrichment/Networks",
  custom_pathways = "Outputs/02_Pathway_Enrichment/Custom_Pathways"
)

# Create directories
for(dir_path in output_structure) {
  dir.create(dir_path, showWarnings = FALSE, recursive = TRUE)
}

# Load differential expression results (assuming from script 01)
if(file.exists("Outputs/01_Differential_Expression/Tables/differential_expression_results.csv")) {
  de_results <- read.csv("Outputs/01_Differential_Expression/Tables/differential_expression_results.csv")
  cat("Loaded DE results with", nrow(de_results), "genes\n")
} else {
  stop("Differential expression results not found. Please run 01_Differential_Expression.R first.")
}

# Load processed data for GSEA
dge_normalized <- readRDS("Outputs/01_Processing_QC/Data/dge_normalized_final.rds")
final_metadata <- readRDS("Outputs/01_Processing_QC/Data/sample_metadata_final.rds")
logcpm_final <- readRDS("Outputs/01_Processing_QC/Data/logcpm_normalized_final.rds")

# ========================================================================
# SECTION 3: GENE ID CONVERSION AND ANNOTATION
# ========================================================================

cat("\n=== GENE ID CONVERSION AND ANNOTATION ===\n")

# Function to convert gene symbols to Entrez IDs
convert_gene_ids <- function(gene_symbols) {
  cat("Converting gene symbols to Entrez IDs...\n")
  
  # Convert using org.Mm.eg.db
  entrez_ids <- mapIds(org.Mm.eg.db, 
                      keys = gene_symbols,
                      column = "ENTREZID",
                      keytype = "SYMBOL",
                      multiVals = "first")
  
  # Create conversion table
  conversion_table <- data.frame(
    symbol = gene_symbols,
    entrez = entrez_ids,
    stringsAsFactors = FALSE
  )
  
  # Remove NAs
  conversion_table <- conversion_table[!is.na(conversion_table$entrez), ]
  
  cat("Successfully converted", nrow(conversion_table), "out of", length(gene_symbols), "gene symbols\n")
  
  return(conversion_table)
}

# Convert DE results gene IDs
gene_conversion <- convert_gene_ids(de_results$gene_id)

# Add Entrez IDs to DE results
de_results_annotated <- de_results %>%
  left_join(gene_conversion, by = c("gene_id" = "symbol")) %>%
  filter(!is.na(entrez))

cat("DE results after ID conversion:", nrow(de_results_annotated), "genes\n")

# Save annotated results
write.csv(de_results_annotated, 
          file.path(output_structure$data, "de_results_with_entrez.csv"),
          row.names = FALSE)

# ========================================================================
# SECTION 4: COMPLEMENT PATHWAY DEFINITION
# ========================================================================

cat("\n=== DEFINING COMPLEMENT PATHWAY GENE SETS ===\n")

# Define comprehensive complement pathway gene sets
define_complement_pathways <- function() {
  
  complement_pathways <- list(
    # Classical pathway
    "Complement_Classical" = c(
      "C1qa", "C1qb", "C1qc", "C1r", "C1s", "C2", "C3", "C4a", "C4b"
    ),
    
    # Alternative pathway
    "Complement_Alternative" = c(
      "C3", "Cfb", "Cfd", "Cfh", "Cfi", "Cfp", "C5", "C6", "C7", "C8a", "C8b", "C8g", "C9"
    ),
    
    # Lectin pathway
    "Complement_Lectin" = c(
      "Mbl1", "Mbl2", "Masp1", "Masp2", "C2", "C3", "C4a", "C4b"
    ),
    
    # Terminal pathway (MAC formation)
    "Complement_Terminal" = c(
      "C5", "C6", "C7", "C8a", "C8b", "C8g", "C9"
    ),
    
    # Complement regulators
    "Complement_Regulators" = c(
      "Cd55", "Cd46", "Cd35", "Cr1", "Cr2", "Crry", "Cfh", "Cfi", "C4bp"
    ),
    
    # Complement receptors
    "Complement_Receptors" = c(
      "C3ar1", "C5ar1", "Cr1", "Cr2", "Itgam", "Itgax", "Itgb2", "Cd35"
    ),
    
    # Complete complement system
    "Complement_Complete" = c(
      "C1qa", "C1qb", "C1qc", "C1r", "C1s", "C2", "C3", "C4a", "C4b",
      "Cfb", "Cfd", "Cfh", "Cfi", "Cfp", "C5", "C6", "C7", "C8a", "C8b", "C8g", "C9",
      "Mbl1", "Mbl2", "Masp1", "Masp2", "Cd55", "Cd46", "Cd35", "Cr1", "Cr2", "Crry",
      "C3ar1", "C5ar1", "Itgam", "Itgax", "Itgb2", "C4bp"
    )
  )
  
  return(complement_pathways)
}

# Get complement pathways
complement_pathways <- define_complement_pathways()

# Convert to Entrez IDs
complement_pathways_entrez <- lapply(complement_pathways, function(genes) {
  converted <- gene_conversion$entrez[match(genes, gene_conversion$symbol)]
  converted[!is.na(converted)]
})

# Remove empty pathways
complement_pathways_entrez <- complement_pathways_entrez[sapply(complement_pathways_entrez, length) > 0]

cat("Defined", length(complement_pathways_entrez), "complement pathway gene sets\n")
for(i in 1:length(complement_pathways_entrez)) {
  cat("-", names(complement_pathways_entrez)[i], ":", length(complement_pathways_entrez[[i]]), "genes\n")
}

# ========================================================================
# SECTION 5: GENE ONTOLOGY (GO) ENRICHMENT ANALYSIS
# ========================================================================

cat("\n=== GENE ONTOLOGY ENRICHMENT ANALYSIS ===\n")

# Function to perform comprehensive GO analysis
perform_go_analysis <- function(de_genes, universe_genes, output_dir) {
  cat("Performing GO enrichment analysis...\n")
  
  go_results <- list()
  
  # Define significance thresholds
  sig_genes <- de_genes$entrez[de_genes$adj.P.Val < 0.05]
  up_genes <- de_genes$entrez[de_genes$adj.P.Val < 0.05 & de_genes$logFC > 0]
  down_genes <- de_genes$entrez[de_genes$adj.P.Val < 0.05 & de_genes$logFC < 0]
  
  # GO Biological Process
  if(length(sig_genes) > 0) {
    go_bp <- enrichGO(gene = sig_genes,
                      universe = universe_genes,
                      OrgDb = org.Mm.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.1,
                      readable = TRUE)
    go_results$BP_all <- go_bp
  }
  
  # GO Molecular Function
  if(length(sig_genes) > 0) {
    go_mf <- enrichGO(gene = sig_genes,
                      universe = universe_genes,
                      OrgDb = org.Mm.eg.db,
                      ont = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.1,
                      readable = TRUE)
    go_results$MF_all <- go_mf
  }
  
  # GO Cellular Component
  if(length(sig_genes) > 0) {
    go_cc <- enrichGO(gene = sig_genes,
                      universe = universe_genes,
                      OrgDb = org.Mm.eg.db,
                      ont = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.1,
                      readable = TRUE)
    go_results$CC_all <- go_cc
  }
  
  # Upregulated genes GO analysis
  if(length(up_genes) > 5) {
    go_bp_up <- enrichGO(gene = up_genes,
                         universe = universe_genes,
                         OrgDb = org.Mm.eg.db,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.1,
                         readable = TRUE)
    go_results$BP_up <- go_bp_up
  }
  
  # Downregulated genes GO analysis
  if(length(down_genes) > 5) {
    go_bp_down <- enrichGO(gene = down_genes,
                           universe = universe_genes,
                           OrgDb = org.Mm.eg.db,
                           ont = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.1,
                           readable = TRUE)
    go_results$BP_down <- go_bp_down
  }
  
  return(go_results)
}

# Perform GO analysis
universe_genes <- de_results_annotated$entrez
go_results <- perform_go_analysis(de_results_annotated, universe_genes, output_structure$go_analysis)

# Save GO results
for(analysis_name in names(go_results)) {
  if(!is.null(go_results[[analysis_name]]) && nrow(go_results[[analysis_name]]) > 0) {
    write.csv(as.data.frame(go_results[[analysis_name]]),
              file.path(output_structure$go_analysis, paste0("GO_", analysis_name, "_results.csv")),
              row.names = FALSE)
    cat("Saved", analysis_name, "results:", nrow(go_results[[analysis_name]]), "terms\n")
  }
}

# ========================================================================
# SECTION 6: KEGG PATHWAY ENRICHMENT ANALYSIS
# ========================================================================

cat("\n=== KEGG PATHWAY ENRICHMENT ANALYSIS ===\n")

# Function to perform KEGG analysis
perform_kegg_analysis <- function(de_genes, universe_genes, output_dir) {
  cat("Performing KEGG pathway analysis...\n")
  
  kegg_results <- list()
  
  # All significant genes
  sig_genes <- de_genes$entrez[de_genes$adj.P.Val < 0.05]
  up_genes <- de_genes$entrez[de_genes$adj.P.Val < 0.05 & de_genes$logFC > 0]
  down_genes <- de_genes$entrez[de_genes$adj.P.Val < 0.05 & de_genes$logFC < 0]
  
  # KEGG pathway enrichment - all significant genes
  if(length(sig_genes) > 0) {
    kegg_all <- enrichKEGG(gene = sig_genes,
                           universe = universe_genes,
                           organism = "mmu",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.1)
    if(!is.null(kegg_all) && nrow(kegg_all) > 0) {
      kegg_all <- setReadable(kegg_all, org.Mm.eg.db, keyType = "ENTREZID")
      kegg_results$all <- kegg_all
    }
  }
  
  # KEGG pathway enrichment - upregulated genes
  if(length(up_genes) > 5) {
    kegg_up <- enrichKEGG(gene = up_genes,
                          universe = universe_genes,
                          organism = "mmu",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.1)
    if(!is.null(kegg_up) && nrow(kegg_up) > 0) {
      kegg_up <- setReadable(kegg_up, org.Mm.eg.db, keyType = "ENTREZID")
      kegg_results$up <- kegg_up
    }
  }
  
  # KEGG pathway enrichment - downregulated genes
  if(length(down_genes) > 5) {
    kegg_down <- enrichKEGG(gene = down_genes,
                            universe = universe_genes,
                            organism = "mmu",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.1)
    if(!is.null(kegg_down) && nrow(kegg_down) > 0) {
      kegg_down <- setReadable(kegg_down, org.Mm.eg.db, keyType = "ENTREZID")
      kegg_results$down <- kegg_down
    }
  }
  
  return(kegg_results)
}

# Perform KEGG analysis
kegg_results <- perform_kegg_analysis(de_results_annotated, universe_genes, output_structure$kegg_analysis)

# Save KEGG results
for(analysis_name in names(kegg_results)) {
  if(!is.null(kegg_results[[analysis_name]]) && nrow(kegg_results[[analysis_name]]) > 0) {
    write.csv(as.data.frame(kegg_results[[analysis_name]]),
              file.path(output_structure$kegg_analysis, paste0("KEGG_", analysis_name, "_results.csv")),
              row.names = FALSE)
    cat("Saved KEGG", analysis_name, "results:", nrow(kegg_results[[analysis_name]]), "pathways\n")
  }
}

# ========================================================================
# SECTION 7: CUSTOM COMPLEMENT PATHWAY ANALYSIS
# ========================================================================

cat("\n=== CUSTOM COMPLEMENT PATHWAY ANALYSIS ===\n")

# Function to test complement pathway enrichment
test_complement_enrichment <- function(de_genes, complement_pathways, output_dir) {
  cat("Testing complement pathway enrichment...\n")
  
  complement_results <- list()
  
  # Get significant genes
  sig_genes <- de_genes$entrez[de_genes$adj.P.Val < 0.05]
  all_genes <- de_genes$entrez
  
  # Test each complement pathway
  for(pathway_name in names(complement_pathways)) {
    pathway_genes <- complement_pathways[[pathway_name]]
    
    # Fisher's exact test
    genes_in_pathway <- sum(sig_genes %in% pathway_genes)
    genes_not_in_pathway <- length(sig_genes) - genes_in_pathway
    background_in_pathway <- sum(all_genes %in% pathway_genes) - genes_in_pathway
    background_not_in_pathway <- length(all_genes) - sum(all_genes %in% pathway_genes) - genes_not_in_pathway
    
    # Create contingency table
    cont_table <- matrix(c(genes_in_pathway, genes_not_in_pathway,
                          background_in_pathway, background_not_in_pathway),
                        nrow = 2,
                        dimnames = list(c("In_pathway", "Not_in_pathway"),
                                       c("Significant", "Background")))
    
    # Perform Fisher's exact test
    fisher_result <- fisher.test(cont_table, alternative = "greater")
    
    # Calculate enrichment metrics
    overlap_genes <- intersect(sig_genes, pathway_genes)
    overlap_symbols <- de_genes$gene_id[de_genes$entrez %in% overlap_genes]
    
    complement_results[[pathway_name]] <- list(
      pathway = pathway_name,
      total_genes = length(pathway_genes),
      significant_genes = genes_in_pathway,
      overlap_genes = overlap_symbols,
      p_value = fisher_result$p.value,
      odds_ratio = fisher_result$estimate,
      enrichment_score = (genes_in_pathway / length(sig_genes)) / (length(pathway_genes) / length(all_genes))
    )
  }
  
  # Convert to data frame
  complement_df <- data.frame(
    pathway = names(complement_results),
    total_genes = sapply(complement_results, function(x) x$total_genes),
    significant_genes = sapply(complement_results, function(x) x$significant_genes),
    p_value = sapply(complement_results, function(x) x$p_value),
    odds_ratio = sapply(complement_results, function(x) x$odds_ratio),
    enrichment_score = sapply(complement_results, function(x) x$enrichment_score),
    stringsAsFactors = FALSE
  )
  
  # Adjust p-values
  complement_df$adj_p_value <- p.adjust(complement_df$p_value, method = "BH")
  
  # Add overlap genes
  complement_df$overlap_genes <- sapply(complement_results, function(x) paste(x$overlap_genes, collapse = ";"))
  
  # Sort by p-value
  complement_df <- complement_df[order(complement_df$p_value), ]
  
  return(list(results = complement_df, detailed = complement_results))
}

# Test complement pathway enrichment
complement_enrichment <- test_complement_enrichment(de_results_annotated, complement_pathways_entrez, output_structure$complement)

# Save complement results
write.csv(complement_enrichment$results,
          file.path(output_structure$complement, "complement_pathway_enrichment.csv"),
          row.names = FALSE)

cat("Complement pathway enrichment results:\n")
print(complement_enrichment$results[complement_enrichment$results$adj_p_value < 0.05, ])

# ========================================================================
# SECTION 8: GENE SET ENRICHMENT ANALYSIS (GSEA)
# ========================================================================

cat("\n=== GENE SET ENRICHMENT ANALYSIS (GSEA) ===\n")

# Function to perform GSEA
perform_gsea_analysis <- function(de_genes, output_dir) {
  cat("Performing GSEA analysis...\n")
  
  # Create ranked gene list
  gene_list <- de_genes$logFC
  names(gene_list) <- de_genes$entrez
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  gsea_results <- list()
  
  # GSEA GO Biological Process
  tryCatch({
    gsea_go_bp <- gseGO(geneList = gene_list,
                        OrgDb = org.Mm.eg.db,
                        ont = "BP",
                        keyType = "ENTREZID",
                        minGSSize = 15,
                        maxGSSize = 500,
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH")
    if(!is.null(gsea_go_bp) && nrow(gsea_go_bp) > 0) {
      gsea_results$GO_BP <- gsea_go_bp
    }
  }, error = function(e) {
    cat("GSEA GO BP failed:", e$message, "\n")
  })
  
  # GSEA KEGG
  tryCatch({
    gsea_kegg <- gseKEGG(geneList = gene_list,
                         organism = "mmu",
                         keyType = "kegg",
                         minGSSize = 15,
                         maxGSSize = 500,
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "BH")
    if(!is.null(gsea_kegg) && nrow(gsea_kegg) > 0) {
      gsea_kegg <- setReadable(gsea_kegg, org.Mm.eg.db, keyType = "ENTREZID")
      gsea_results$KEGG <- gsea_kegg
    }
  }, error = function(e) {
    cat("GSEA KEGG failed:", e$message, "\n")
  })
  
  return(list(results = gsea_results, gene_list = gene_list))
}

# Perform GSEA
gsea_analysis <- perform_gsea_analysis(de_results_annotated, output_structure$gsea)

# Save GSEA results
for(analysis_name in names(gsea_analysis$results)) {
  if(!is.null(gsea_analysis$results[[analysis_name]]) && nrow(gsea_analysis$results[[analysis_name]]) > 0) {
    write.csv(as.data.frame(gsea_analysis$results[[analysis_name]]),
              file.path(output_structure$gsea, paste0("GSEA_", analysis_name, "_results.csv")),
              row.names = FALSE)
    cat("Saved GSEA", analysis_name, "results:", nrow(gsea_analysis$results[[analysis_name]]), "pathways\n")
  }
}

# ========================================================================
# SECTION 9: PUBLICATION-QUALITY VISUALIZATIONS
# ========================================================================

cat("\n=== CREATING PUBLICATION-QUALITY PLOTS ===\n")

# Function to create comprehensive pathway visualization
create_pathway_plots <- function(go_results, kegg_results, complement_results, gsea_results, output_dir) {
  cat("Creating pathway enrichment visualizations...\n")
  
  # GO enrichment plots
  if(!is.null(go_results$BP_all) && nrow(go_results$BP_all) > 0) {
    # GO dotplot
    png(file.path(output_dir, "Figure_3A_GO_BP_dotplot.png"), 
        width = 3000, height = 2400, res = 300)
    p1 <- dotplot(go_results$BP_all, showCategory = 20) +
      ggtitle("GO Biological Process Enrichment") +
      theme(text = element_text(size = 12))
    print(p1)
    dev.off()
    
    # GO enrichment map
    if(nrow(go_results$BP_all) > 5) {
      tryCatch({
        go_simplified <- simplify(go_results$BP_all)
        png(file.path(output_dir, "Figure_3B_GO_BP_enrichmap.png"), 
            width = 3000, height = 2400, res = 300)
        p2 <- emapplot(pairwise_termsim(go_simplified), showCategory = 15) +
          ggtitle("GO BP Enrichment Network")
        print(p2)
        dev.off()
      }, error = function(e) {
        cat("GO enrichment map failed:", e$message, "\n")
      })
    }
  }
  
  # KEGG pathway plots
  if(!is.null(kegg_results$all) && nrow(kegg_results$all) > 0) {
    # KEGG dotplot
    png(file.path(output_dir, "Figure_3C_KEGG_dotplot.png"), 
        width = 3000, height = 2400, res = 300)
    p3 <- dotplot(kegg_results$all, showCategory = 15) +
      ggtitle("KEGG Pathway Enrichment") +
      theme(text = element_text(size = 12))
    print(p3)
    dev.off()
  }
  
  # Complement pathway visualization
  comp_sig <- complement_results[complement_results$adj_p_value < 0.05, ]
  if(nrow(comp_sig) > 0) {
    png(file.path(output_dir, "Figure_3D_Complement_enrichment.png"), 
        width = 2400, height = 1800, res = 300)
    
    comp_plot_data <- comp_sig
    comp_plot_data$neg_log_p <- -log10(comp_plot_data$adj_p_value)
    
    p4 <- ggplot(comp_plot_data, aes(x = reorder(pathway, neg_log_p), y = neg_log_p)) +
      geom_col(fill = "steelblue", alpha = 0.7) +
      geom_text(aes(label = significant_genes), hjust = -0.1, size = 3) +
      coord_flip() +
      labs(title = "Complement Pathway Enrichment",
           x = "Pathway",
           y = "-log10(Adjusted P-value)",
           subtitle = "Numbers indicate significant genes in pathway") +
      theme_minimal() +
      theme(text = element_text(size = 12))
    print(p4)
    dev.off()
  }
  
  # GSEA plots
  if(!is.null(gsea_results$GO_BP) && nrow(gsea_results$GO_BP) > 0) {
    # GSEA dotplot
    png(file.path(output_dir, "Figure_3E_GSEA_GO_dotplot.png"), 
        width = 3000, height = 2400, res = 300)
    p5 <- dotplot(gsea_results$GO_BP, showCategory = 15, split = ".sign") +
      facet_grid(.~.sign) +
      ggtitle("GSEA GO Biological Process") +
      theme(text = element_text(size = 12))
    print(p5)
    dev.off()
  }
  
  # Combined summary plot
  create_pathway_summary_plot(go_results, kegg_results, complement_results, output_dir)
}

# Function to create pathway summary plot
create_pathway_summary_plot <- function(go_results, kegg_results, complement_results, output_dir) {
  
  # Collect enrichment counts
  enrichment_summary <- data.frame(
    Database = character(),
    Category = character(),
    Significant_Terms = integer(),
    stringsAsFactors = FALSE
  )
  
  # GO results
  if(!is.null(go_results$BP_all)) {
    enrichment_summary <- rbind(enrichment_summary,
                               data.frame(Database = "GO", Category = "Biological Process", 
                                         Significant_Terms = sum(go_results$BP_all$p.adjust < 0.05)))
  }
  if(!is.null(go_results$MF_all)) {
    enrichment_summary <- rbind(enrichment_summary,
                               data.frame(Database = "GO", Category = "Molecular Function", 
                                         Significant_Terms = sum(go_results$MF_all$p.adjust < 0.05)))
  }
  if(!is.null(go_results$CC_all)) {
    enrichment_summary <- rbind(enrichment_summary,
                               data.frame(Database = "GO", Category = "Cellular Component", 
                                         Significant_Terms = sum(go_results$CC_all$p.adjust < 0.05)))
  }
  
  # KEGG results
  if(!is.null(kegg_results$all)) {
    enrichment_summary <- rbind(enrichment_summary,
                               data.frame(Database = "KEGG", Category = "Pathways", 
                                         Significant_Terms = sum(kegg_results$all$p.adjust < 0.05)))
  }
  
  # Complement results
  comp_sig_count <- sum(complement_results$adj_p_value < 0.05)
  enrichment_summary <- rbind(enrichment_summary,
                             data.frame(Database = "Custom", Category = "Complement", 
                                       Significant_Terms = comp_sig_count))
  
  # Create summary plot
  if(nrow(enrichment_summary) > 0) {
    png(file.path(output_dir, "Figure_3F_Pathway_Summary.png"), 
        width = 2400, height = 1800, res = 300)
    
    p_summary <- ggplot(enrichment_summary, aes(x = paste(Database, Category, sep = "\n"), 
                                                y = Significant_Terms, fill = Database)) +
      geom_col(alpha = 0.7) +
      geom_text(aes(label = Significant_Terms), vjust = -0.3, size = 4) +
      labs(title = "Pathway Enrichment Summary",
           x = "Database/Category",
           y = "Number of Significant Terms",
           fill = "Database") +
      theme_minimal() +
      theme(text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_fill_brewer(type = "qual", palette = "Set2")
    
    print(p_summary)
    dev.off()
  }
}

# Create all pathway plots
create_pathway_plots(go_results, kegg_results, complement_enrichment$results, gsea_analysis$results, output_structure$plots)

# ========================================================================
# SECTION 10: COMPREHENSIVE REPORT GENERATION
# ========================================================================

cat("\n=== GENERATING COMPREHENSIVE PATHWAY REPORT ===\n")

# Function to create pathway analysis report
create_pathway_report <- function(output_file, go_results, kegg_results, complement_results, gsea_results) {
  
  sink(output_file)
  
  cat("=================================================================\n")
  cat("GSE118918 PATHWAY ENRICHMENT ANALYSIS REPORT\n")
  cat("=================================================================\n\n")
  
  cat("ANALYSIS OVERVIEW\n")
  cat("-----------------\n")
  cat("Dataset: GSE118918 (Nucleus Accumbens, Mock vs Morphine)\n")
  cat("Analysis Date:", as.character(Sys.time()), "\n")
  cat("Focus: Complement pathway and morphine response mechanisms\n\n")
  
  cat("GENE ONTOLOGY ENRICHMENT\n")
  cat("------------------------\n")
  if(!is.null(go_results$BP_all)) {
    bp_sig <- sum(go_results$BP_all$p.adjust < 0.05)
    cat("Biological Process terms (FDR < 0.05):", bp_sig, "\n")
    if(bp_sig > 0) {
      top_bp <- head(go_results$BP_all[order(go_results$BP_all$p.adjust), ], 5)
      cat("Top 5 BP terms:\n")
      for(i in 1:min(5, nrow(top_bp))) {
        cat("  ", i, ". ", top_bp$Description[i], " (p.adj = ", 
            sprintf("%.2e", top_bp$p.adjust[i]), ")\n", sep = "")
      }
    }
  }
  
  if(!is.null(go_results$MF_all)) {
    mf_sig <- sum(go_results$MF_all$p.adjust < 0.05)
    cat("Molecular Function terms (FDR < 0.05):", mf_sig, "\n")
  }
  
  if(!is.null(go_results$CC_all)) {
    cc_sig <- sum(go_results$CC_all$p.adjust < 0.05)
    cat("Cellular Component terms (FDR < 0.05):", cc_sig, "\n")
  }
  cat("\n")
  
  cat("KEGG PATHWAY ENRICHMENT\n")
  cat("-----------------------\n")
  if(!is.null(kegg_results$all)) {
    kegg_sig <- sum(kegg_results$all$p.adjust < 0.05)
    cat("Significant KEGG pathways (FDR < 0.05):", kegg_sig, "\n")
    if(kegg_sig > 0) {
      top_kegg <- head(kegg_results$all[order(kegg_results$all$p.adjust), ], 5)
      cat("Top 5 KEGG pathways:\n")
      for(i in 1:min(5, nrow(top_kegg))) {
        cat("  ", i, ". ", top_kegg$Description[i], " (p.adj = ", 
            sprintf("%.2e", top_kegg$p.adjust[i]), ")\n", sep = "")
      }
    }
  } else {
    cat("No significant KEGG pathways found\n")
  }
  cat("\n")
  
  cat("COMPLEMENT PATHWAY ANALYSIS\n")
  cat("---------------------------\n")
  comp_sig <- complement_results[complement_results$adj_p_value < 0.05, ]
  cat("Significant complement pathways (FDR < 0.05):", nrow(comp_sig), "\n")
  if(nrow(comp_sig) > 0) {
    cat("Enriched complement pathways:\n")
    for(i in 1:nrow(comp_sig)) {
      cat("  ", i, ". ", comp_sig$pathway[i], 
          " (", comp_sig$significant_genes[i], " genes, p.adj = ", 
          sprintf("%.2e", comp_sig$adj_p_value[i]), ")\n", sep = "")
    }
    cat("\nGenes contributing to complement enrichment:\n")
    all_comp_genes <- unique(unlist(strsplit(comp_sig$overlap_genes, ";")))
    all_comp_genes <- all_comp_genes[all_comp_genes != ""]
    cat(paste(all_comp_genes, collapse = ", "), "\n")
  } else {
    cat("No significant complement pathway enrichment detected\n")
  }
  cat("\n")
  
  cat("GSEA ANALYSIS\n")
  cat("-------------\n")
  if(!is.null(gsea_results$GO_BP)) {
    gsea_sig <- sum(gsea_results$GO_BP$p.adjust < 0.05)
    cat("Significant GSEA GO BP terms (FDR < 0.05):", gsea_sig, "\n")
  }
  if(!is.null(gsea_results$KEGG)) {
    gsea_kegg_sig <- sum(gsea_results$KEGG$p.adjust < 0.05)
    cat("Significant GSEA KEGG pathways (FDR < 0.05):", gsea_kegg_sig, "\n")
  }
  cat("\n")
  
  cat("KEY FINDINGS\n")
  cat("------------\n")
  cat("1. Morphine treatment effects on complement system:\n")
  if(nrow(comp_sig) > 0) {
    cat("   - Significant complement pathway dysregulation detected\n")
    cat("   - Most affected pathway:", comp_sig$pathway[1], "\n")
  } else {
    cat("   - No significant complement pathway changes\n")
  }
  
  cat("2. Overall pathway landscape:\n")
  total_pathways <- 0
  if(!is.null(go_results$BP_all)) total_pathways <- total_pathways + sum(go_results$BP_all$p.adjust < 0.05)
  if(!is.null(kegg_results$all)) total_pathways <- total_pathways + sum(kegg_results$all$p.adjust < 0.05)
  total_pathways <- total_pathways + nrow(comp_sig)
  cat("   - Total significant pathways/terms:", total_pathways, "\n")
  
  cat("\n3. Therapeutic implications:\n")
  if(nrow(comp_sig) > 0) {
    cat("   - Complement system represents potential therapeutic target\n")
    cat("   - Neuroinflammatory pathways may be involved in morphine response\n")
  }
  
  cat("\nRECOMMENDations\n")
  cat("---------------\n")
  cat("1. Further investigate complement genes with highest fold changes\n")
  cat("2. Validate complement pathway activity with functional assays\n")
  cat("3. Consider complement inhibitors as potential therapeutics\n")
  cat("4. Integrate with other addiction-related datasets\n")
  
  cat("\n=================================================================\n")
  cat("PATHWAY ANALYSIS COMPLETED\n")
  cat("=================================================================\n")
  
  sink()
}

# Generate comprehensive report
create_pathway_report(
  file.path(output_structure$reports, "GSE118918_Pathway_Analysis_Report.txt"),
  go_results, kegg_results, complement_enrichment$results, gsea_analysis$results
)

# ========================================================================
# SECTION 11: DATA EXPORT AND SESSION INFO
# ========================================================================

cat("\n=== SAVING RESULTS AND SESSION INFO ===\n")

# Save all analysis objects
analysis_objects <- list(
  go_results = go_results,
  kegg_results = kegg_results,
  complement_enrichment = complement_enrichment,
  gsea_analysis = gsea_analysis,
  gene_conversion = gene_conversion,
  complement_pathways = complement_pathways_entrez
)

saveRDS(analysis_objects, file.path(output_structure$data, "pathway_enrichment_results.rds"))

# Save session information
session_info <- sessionInfo()
writeLines(capture.output(print(session_info)), 
           file.path(output_structure$reports, "session_info.txt"))

# Create analysis summary
analysis_summary <- data.frame(
  analysis_type = c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Complement", "GSEA_GO", "GSEA_KEGG"),
  total_terms = c(
    ifelse(!is.null(go_results$BP_all), nrow(go_results$BP_all), 0),
    ifelse(!is.null(go_results$MF_all), nrow(go_results$MF_all), 0),
    ifelse(!is.null(go_results$CC_all), nrow(go_results$CC_all), 0),
    ifelse(!is.null(kegg_results$all), nrow(kegg_results$all), 0),
    nrow(complement_enrichment$results),
    ifelse(!is.null(gsea_analysis$results$GO_BP), nrow(gsea_analysis$results$GO_BP), 0),
    ifelse(!is.null(gsea_analysis$results$KEGG), nrow(gsea_analysis$results$KEGG), 0)
  ),
  significant_terms = c(
    ifelse(!is.null(go_results$BP_all), sum(go_results$BP_all$p.adjust < 0.05), 0),
    ifelse(!is.null(go_results$MF_all), sum(go_results$MF_all$p.adjust < 0.05), 0),
    ifelse(!is.null(go_results$CC_all), sum(go_results$CC_all$p.adjust < 0.05), 0),
    ifelse(!is.null(kegg_results$all), sum(kegg_results$all$p.adjust < 0.05), 0),
    sum(complement_enrichment$results$adj_p_value < 0.05),
    ifelse(!is.null(gsea_analysis$results$GO_BP), sum(gsea_analysis$results$GO_BP$p.adjust < 0.05), 0),
    ifelse(!is.null(gsea_analysis$results$KEGG), sum(gsea_analysis$results$KEGG$p.adjust < 0.05), 0)
  )
)

write.csv(analysis_summary, 
          file.path(output_structure$tables, "pathway_analysis_summary.csv"),
          row.names = FALSE)

cat("\n=== PATHWAY ENRICHMENT ANALYSIS COMPLETED ===\n")
cat("Output directory:", output_structure$main, "\n")
cat("Key files generated:\n")
cat("- GO enrichment results in:", output_structure$go_analysis, "\n")
cat("- KEGG pathway results in:", output_structure$kegg_analysis, "\n")
cat("- Complement analysis in:", output_structure$complement, "\n")
cat("- GSEA results in:", output_structure$gsea, "\n")
cat("- Publication figures in:", output_structure$plots, "\n")
cat("- Comprehensive report:", file.path(output_structure$reports, "GSE118918_Pathway_Analysis_Report.txt"), "\n")

cat("\nNext recommended step: Cross-dataset integration analysis\n")
