# Comprehensive Pathway Enrichment Analysis for Morphine OUD Study
# ==============================================================

# Set reproducibility seed
set.seed(42)

# Load required libraries
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(DOSE)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(RColorBrewer)
  library(pheatmap)
  library(gridExtra)
})

# ----------------------
# SETUP AND DIRECTORIES
# ----------------------

# Base directory and paths
base_dir <- "/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus/GSE289002"
de_dir <- file.path(base_dir, "Outputs/Differential_Expression/Tables")
results_dir <- file.path(base_dir, "Outputs/Differential_Expression/Results")

# Create pathway analysis output directories
pathway_dirs <- list(
  base = file.path(base_dir, "Outputs/Pathway_Analysis"),
  results = file.path(base_dir, "Outputs/Pathway_Analysis/Results"),
  plots = file.path(base_dir, "Outputs/Pathway_Analysis/Plots"),
  tables = file.path(base_dir, "Outputs/Pathway_Analysis/Tables"),
  comparisons = file.path(base_dir, "Outputs/Pathway_Analysis/Comparisons")
)

# Create directories
lapply(pathway_dirs, function(x) dir.create(x, showWarnings = FALSE, recursive = TRUE))

cat("=== COMPREHENSIVE PATHWAY ENRICHMENT ANALYSIS ===\n")
cat("Analysis Date:", as.character(Sys.time()), "\n")
cat("Output Directory:", pathway_dirs$base, "\n\n")

# ----------------------
# 1. LOAD AND PREPARE GENE LISTS
# ----------------------

cat("=== LOADING DIFFERENTIAL EXPRESSION RESULTS ===\n")

# Function to load and prepare gene lists
load_gene_list <- function(file_path, list_name) {
  if (file.exists(file_path)) {
    genes <- read.csv(file_path)
    cat("Loaded", list_name, ":", nrow(genes), "genes\n")
    return(genes$gene_id)
  } else {
    cat("Warning:", file_path, "not found\n")
    return(NULL)
  }
}

# Load all gene lists
gene_lists <- list(
  # Main treatment comparisons
  acute_vs_control = load_gene_list(file.path(de_dir, "Acute_vs_Control_significant_genes.csv"), "Acute vs Control"),
  short_vs_control = load_gene_list(file.path(de_dir, "Short_vs_Control_significant_genes.csv"), "Short vs Control"),
  chronic_vs_control = load_gene_list(file.path(de_dir, "Chronic_vs_Control_significant_genes.csv"), "Chronic vs Control"),
  
  # Time course comparisons
  short_vs_acute = load_gene_list(file.path(de_dir, "Short_vs_Acute_significant_genes.csv"), "Short vs Acute"),
  chronic_vs_acute = load_gene_list(file.path(de_dir, "Chronic_vs_Acute_significant_genes.csv"), "Chronic vs Acute"),
  chronic_vs_short = load_gene_list(file.path(de_dir, "Chronic_vs_Short_significant_genes.csv"), "Chronic vs Short"),
  
  # Special gene sets
  acute_specific = load_gene_list(file.path(results_dir, "acute_specific_genes.csv"), "Acute-specific"),
  chronic_specific = load_gene_list(file.path(results_dir, "chronic_specific_genes.csv"), "Chronic-specific"),
  shared_acute_chronic = load_gene_list(file.path(results_dir, "acute_chronic_shared_genes.csv"), "Shared Acute/Chronic")
)

# Remove NULL entries
gene_lists <- gene_lists[!sapply(gene_lists, is.null)]

cat("\nSuccessfully loaded", length(gene_lists), "gene lists\n")

# Check the format of gene IDs
cat("\nChecking gene ID format:\n")
sample_genes <- head(gene_lists[[1]], 5)
cat("Sample genes:", paste(sample_genes, collapse = ", "), "\n")

# Determine if these are Ensembl IDs or symbols
is_ensembl <- all(grepl("^ENSMUSG", sample_genes))
cat("Gene ID format:", ifelse(is_ensembl, "Ensembl IDs", "Gene Symbols"), "\n\n")

# ----------------------
# 2. GENE ID CONVERSION (ENSEMBL TO SYMBOLS TO ENTREZ)
# ----------------------

cat("=== CONVERTING GENE IDs FOR PATHWAY ANALYSIS ===\n")

# Function to convert Ensembl IDs to gene symbols and then to Entrez IDs
convert_ensembl_to_entrez <- function(ensembl_ids, list_name) {
  cat("Converting", list_name, "- Input:", length(ensembl_ids), "Ensembl IDs\n")
  
  if (length(ensembl_ids) == 0) {
    return(NULL)
  }
  
  # Step 1: Convert Ensembl to Symbols
  tryCatch({
    ensembl_to_symbol <- bitr(ensembl_ids, 
                             fromType = "ENSEMBL", 
                             toType = "SYMBOL", 
                             OrgDb = org.Mm.eg.db)
    
    cat("  - Ensembl to Symbol: ", nrow(ensembl_to_symbol), " genes converted\n", sep = "")
    
    if (nrow(ensembl_to_symbol) == 0) {
      cat("  - No valid conversions found\n")
      return(NULL)
    }
    
    # Step 2: Convert Symbols to Entrez IDs
    symbol_to_entrez <- bitr(ensembl_to_symbol$SYMBOL, 
                            fromType = "SYMBOL", 
                            toType = "ENTREZID", 
                            OrgDb = org.Mm.eg.db)
    
    cat("  - Symbol to Entrez: ", nrow(symbol_to_entrez), " genes converted\n", sep = "")
    cat("  - Final conversion rate: ", round(nrow(symbol_to_entrez)/length(ensembl_ids)*100, 1), "%\n", sep = "")
    
    # Also save the gene symbol mappings for reference
    final_mapping <- merge(ensembl_to_symbol, symbol_to_entrez, by = "SYMBOL")
    write.csv(final_mapping, 
              file.path(pathway_dirs$tables, paste0(list_name, "_gene_id_mapping.csv")), 
              row.names = FALSE)
    
    return(list(
      entrez_ids = symbol_to_entrez$ENTREZID,
      symbols = symbol_to_entrez$SYMBOL,
      mapping = final_mapping
    ))
    
  }, error = function(e) {
    cat("  - Conversion failed:", e$message, "\n")
    return(NULL)
  })
}

# Function to convert gene symbols to Entrez IDs (if already symbols)
convert_symbols_to_entrez <- function(gene_symbols, list_name) {
  cat("Converting", list_name, "- Input:", length(gene_symbols), "gene symbols\n")
  
  if (length(gene_symbols) == 0) {
    return(NULL)
  }
  
  tryCatch({
    converted <- bitr(gene_symbols, 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Mm.eg.db)
    
    cat("  - Successfully converted:", nrow(converted), "genes\n")
    cat("  - Conversion rate:", round(nrow(converted)/length(gene_symbols)*100, 1), "%\n")
    
    return(list(
      entrez_ids = converted$ENTREZID,
      symbols = converted$SYMBOL,
      mapping = converted
    ))
    
  }, error = function(e) {
    cat("  - Conversion failed:", e$message, "\n")
    return(NULL)
  })
}

# Convert all gene lists based on their format
converted_gene_lists <- list()

for (list_name in names(gene_lists)) {
  if (is_ensembl) {
    converted_gene_lists[[list_name]] <- convert_ensembl_to_entrez(gene_lists[[list_name]], list_name)
  } else {
    converted_gene_lists[[list_name]] <- convert_symbols_to_entrez(gene_lists[[list_name]], list_name)
  }
}

# Remove NULL entries and extract just the Entrez IDs
entrez_lists <- list()
for (list_name in names(converted_gene_lists)) {
  if (!is.null(converted_gene_lists[[list_name]]) && 
      length(converted_gene_lists[[list_name]]$entrez_ids) >= 5) {
    entrez_lists[[list_name]] <- converted_gene_lists[[list_name]]$entrez_ids
  }
}

cat("\nConversion completed for", length(entrez_lists), "gene lists with >=5 genes\n\n")

# ----------------------
# 3. PATHWAY ENRICHMENT ANALYSIS
# ----------------------

cat("=== PERFORMING PATHWAY ENRICHMENT ANALYSIS ===\n")

# Function to perform comprehensive enrichment analysis
perform_enrichment <- function(gene_entrez, analysis_name, min_genes = 5) {
  
  if (length(gene_entrez) < min_genes) {
    cat("Skipping", analysis_name, "- too few genes (", length(gene_entrez), ")\n")
    return(NULL)
  }
  
  cat("Analyzing", analysis_name, "with", length(gene_entrez), "genes\n")
  
  results <- list()
  
  # GO Biological Process
  tryCatch({
    results$go_bp <- enrichGO(gene = gene_entrez,
                              OrgDb = org.Mm.eg.db,
                              ont = "BP",
                              pAdjustMethod = "BH",
                              qvalueCutoff = 0.05,
                              readable = TRUE)
    
    if (!is.null(results$go_bp) && nrow(results$go_bp@result) > 0) {
      cat("  - GO BP: ", nrow(results$go_bp@result), " terms\n", sep = "")
    } else {
      cat("  - GO BP: No significant terms\n")
      results$go_bp <- NULL
    }
  }, error = function(e) {
    cat("  - GO BP: Failed -", e$message, "\n")
    results$go_bp <<- NULL
  })
  
  # GO Molecular Function
  tryCatch({
    results$go_mf <- enrichGO(gene = gene_entrez,
                              OrgDb = org.Mm.eg.db,
                              ont = "MF",
                              pAdjustMethod = "BH",
                              qvalueCutoff = 0.05,
                              readable = TRUE)
    
    if (!is.null(results$go_mf) && nrow(results$go_mf@result) > 0) {
      cat("  - GO MF: ", nrow(results$go_mf@result), " terms\n", sep = "")
    } else {
      cat("  - GO MF: No significant terms\n")
      results$go_mf <- NULL
    }
  }, error = function(e) {
    cat("  - GO MF: Failed -", e$message, "\n")
    results$go_mf <<- NULL
  })
  
  # GO Cellular Component
  tryCatch({
    results$go_cc <- enrichGO(gene = gene_entrez,
                              OrgDb = org.Mm.eg.db,
                              ont = "CC",
                              pAdjustMethod = "BH",
                              qvalueCutoff = 0.05,
                              readable = TRUE)
    
    if (!is.null(results$go_cc) && nrow(results$go_cc@result) > 0) {
      cat("  - GO CC: ", nrow(results$go_cc@result), " terms\n", sep = "")
    } else {
      cat("  - GO CC: No significant terms\n")
      results$go_cc <- NULL
    }
  }, error = function(e) {
    cat("  - GO CC: Failed -", e$message, "\n")
    results$go_cc <<- NULL
  })
  
  # KEGG Pathways
  tryCatch({
    results$kegg <- enrichKEGG(gene = gene_entrez,
                               organism = 'mmu',
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.05)
    
    if (!is.null(results$kegg) && nrow(results$kegg@result) > 0) {
      cat("  - KEGG: ", nrow(results$kegg@result), " pathways\n", sep = "")
    } else {
      cat("  - KEGG: No significant pathways\n")
      results$kegg <- NULL
    }
  }, error = function(e) {
    cat("  - KEGG: Failed -", e$message, "\n")
    results$kegg <<- NULL
  })
  
  # Disease Ontology
  tryCatch({
    results$do <- enrichDO(gene = gene_entrez,
                           ont = "DO",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05,
                           readable = TRUE)
    
    if (!is.null(results$do) && nrow(results$do@result) > 0) {
      cat("  - DO: ", nrow(results$do@result), " terms\n", sep = "")
    } else {
      cat("  - DO: No significant terms\n")
      results$do <- NULL
    }
  }, error = function(e) {
    cat("  - DO: Failed -", e$message, "\n")
    results$do <<- NULL
  })
  
  return(results)
}

# Perform enrichment for all gene lists
enrichment_results <- list()
for (list_name in names(entrez_lists)) {
  enrichment_results[[list_name]] <- perform_enrichment(entrez_lists[[list_name]], list_name)
}

# Remove NULL results
enrichment_results <- enrichment_results[!sapply(enrichment_results, is.null)]

cat("\nEnrichment analysis completed for", length(enrichment_results), "gene lists\n\n")

# ----------------------
# 4. SAVE ENRICHMENT RESULTS
# ----------------------

cat("=== SAVING ENRICHMENT RESULTS ===\n")

# Function to save enrichment results
save_enrichment_results <- function(results, analysis_name) {
  if (is.null(results)) return()
  
  analysis_dir <- file.path(pathway_dirs$tables, analysis_name)
  dir.create(analysis_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (db_name in names(results)) {
    if (!is.null(results[[db_name]]) && nrow(results[[db_name]]@result) > 0) {
      
      # Convert to data frame and save
      result_df <- as.data.frame(results[[db_name]])
      write.csv(result_df, 
                file.path(analysis_dir, paste0(db_name, "_enrichment.csv")), 
                row.names = FALSE)
      
      cat("Saved", db_name, "results for", analysis_name, "(", nrow(result_df), "terms )\n")
    }
  }
}

# Save all results
for (analysis_name in names(enrichment_results)) {
  save_enrichment_results(enrichment_results[[analysis_name]], analysis_name)
}

# ----------------------
# 5. CREATE VISUALIZATIONS
# ----------------------

cat("\n=== CREATING PATHWAY VISUALIZATIONS ===\n")

# Function to create enrichment plots
create_enrichment_plots <- function(results, analysis_name, top_n = 20) {
  if (is.null(results)) return()
  
  plot_dir <- file.path(pathway_dirs$plots, analysis_name)
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (db_name in names(results)) {
    if (!is.null(results[[db_name]]) && nrow(results[[db_name]]@result) > 0) {
      
      # Dot plot
      tryCatch({
        p_dot <- dotplot(results[[db_name]], showCategory = min(top_n, nrow(results[[db_name]]@result))) +
          ggtitle(paste(analysis_name, "-", db_name, "Enrichment")) +
          theme(plot.title = element_text(hjust = 0.5, size = 12))
        
        ggsave(file.path(plot_dir, paste0(db_name, "_dotplot.png")), 
               p_dot, width = 12, height = 8, dpi = 300)
        
        cat("Created dotplot for", analysis_name, "-", db_name, "\n")
      }, error = function(e) {
        cat("Failed to create dotplot for", analysis_name, "-", db_name, ":", e$message, "\n")
      })
      
      # Bar plot
      tryCatch({
        p_bar <- barplot(results[[db_name]], showCategory = min(top_n, nrow(results[[db_name]]@result))) +
          ggtitle(paste(analysis_name, "-", db_name, "Enrichment")) +
          theme(plot.title = element_text(hjust = 0.5, size = 12))
        
        ggsave(file.path(plot_dir, paste0(db_name, "_barplot.png")), 
               p_bar, width = 12, height = 8, dpi = 300)
        
        cat("Created barplot for", analysis_name, "-", db_name, "\n")
      }, error = function(e) {
        cat("Failed to create barplot for", analysis_name, "-", db_name, ":", e$message, "\n")
      })
    }
  }
}

# Create plots for all analyses
for (analysis_name in names(enrichment_results)) {
  create_enrichment_plots(enrichment_results[[analysis_name]], analysis_name)
}

# ----------------------
# 6. ADDICTION-SPECIFIC ANALYSIS
# ----------------------

cat("\n=== ADDICTION-SPECIFIC PATHWAY ANALYSIS ===\n")

# Function to search for addiction-related terms
find_addiction_terms <- function(enrichment_result, search_terms) {
  if (is.null(enrichment_result) || nrow(enrichment_result@result) == 0) {
    return(NULL)
  }
  
  result_df <- as.data.frame(enrichment_result)
  
  # Search in Description column
  addiction_mask <- grepl(paste(search_terms, collapse = "|"), 
                         result_df$Description, ignore.case = TRUE)
  
  if (any(addiction_mask)) {
    return(result_df[addiction_mask, ])
  } else {
    return(NULL)
  }
}

# Addiction-related search terms
addiction_terms <- c("dopamine", "serotonin", "neurotransmitter", "synap", "neural", 
                     "behavior", "reward", "addiction", "drug", "opioid", "morphine",
                     "pleasure", "motivat", "learning", "memory", "stress", "anxiety",
                     "calcium", "signal", "receptor", "binding")

# Search for addiction-related pathways in all results
addiction_pathways <- list()

for (analysis_name in names(enrichment_results)) {
  cat("Searching for addiction-related pathways in", analysis_name, "\n")
  
  analysis_addiction <- list()
  
  for (db_name in names(enrichment_results[[analysis_name]])) {
    addiction_results <- find_addiction_terms(enrichment_results[[analysis_name]][[db_name]], 
                                            addiction_terms)
    
    if (!is.null(addiction_results)) {
      analysis_addiction[[db_name]] <- addiction_results
      cat("  - Found", nrow(addiction_results), "addiction-related", db_name, "terms\n")
    }
  }
  
  if (length(analysis_addiction) > 0) {
    addiction_pathways[[analysis_name]] <- analysis_addiction
  }
}

# Save addiction-specific results
if (length(addiction_pathways) > 0) {
  addiction_dir <- file.path(pathway_dirs$results, "Addiction_Specific")
  dir.create(addiction_dir, showWarnings = FALSE, recursive = TRUE)

  for (analysis_name in names(addiction_pathways)) {
    for (db_name in names(addiction_pathways[[analysis_name]])) {
      write.csv(addiction_pathways[[analysis_name]][[db_name]],
                file.path(addiction_dir, paste0(analysis_name, "_", db_name, "_addiction_terms.csv")),
                row.names = FALSE)
    }
  }
}

# ----------------------
# 7. SUMMARY REPORT
# ----------------------

cat("\n=== GENERATING PATHWAY ANALYSIS SUMMARY ===\n")

# Create comprehensive summary
create_pathway_summary <- function(results_list, output_file) {
  
  sink(output_file)
  
  cat("================================================================\n")
  cat("COMPREHENSIVE PATHWAY ENRICHMENT ANALYSIS REPORT\n")
  cat("================================================================\n\n")
  
  cat("Analysis Date:", as.character(Sys.time()), "\n")
  cat("Gene Lists Analyzed:", length(results_list), "\n")
  cat("Gene ID Format:", ifelse(is_ensembl, "Ensembl IDs (converted)", "Gene Symbols"), "\n\n")
  
  if (length(results_list) == 0) {
    cat("No enrichment results found. This could be due to:\n")
    cat("1. Gene ID conversion issues\n")
    cat("2. Too few genes in each list\n")
    cat("3. No significant pathway enrichment\n\n")
    cat("Check the gene ID mapping files in the Tables directory.\n")
  } else {
    
    # Summary table
    summary_data <- data.frame(
      Analysis = character(),
      GO_BP = integer(),
      GO_MF = integer(),
      GO_CC = integer(),
      KEGG = integer(),
      DO = integer(),
      stringsAsFactors = FALSE
    )
    
    for (analysis_name in names(results_list)) {
      if (!is.null(results_list[[analysis_name]])) {
        
        row_data <- data.frame(
          Analysis = analysis_name,
          GO_BP = ifelse(is.null(results_list[[analysis_name]]$go_bp), 0, 
                        nrow(results_list[[analysis_name]]$go_bp@result)),
          GO_MF = ifelse(is.null(results_list[[analysis_name]]$go_mf), 0, 
                        nrow(results_list[[analysis_name]]$go_mf@result)),
          GO_CC = ifelse(is.null(results_list[[analysis_name]]$go_cc), 0, 
                        nrow(results_list[[analysis_name]]$go_cc@result)),
          KEGG = ifelse(is.null(results_list[[analysis_name]]$kegg), 0, 
                       nrow(results_list[[analysis_name]]$kegg@result)),
          DO = ifelse(is.null(results_list[[analysis_name]]$do), 0, 
                     nrow(results_list[[analysis_name]]$do@result)),
          stringsAsFactors = FALSE
        )
        
        summary_data <- rbind(summary_data, row_data)
      }
    }
    
    cat("ENRICHMENT SUMMARY:\n")
    print(summary_data)
    cat("\n")
    
    # Top pathways for each analysis
    cat("TOP PATHWAYS BY ANALYSIS:\n\n")
    
    for (analysis_name in names(results_list)) {
      if (!is.null(results_list[[analysis_name]])) {
        
        cat("=== ", toupper(analysis_name), " ===\n", sep = "")
        
        # Show top GO BP terms
        if (!is.null(results_list[[analysis_name]]$go_bp) && 
            nrow(results_list[[analysis_name]]$go_bp@result) > 0) {
          cat("\nTop GO Biological Processes:\n")
          top_go <- head(as.data.frame(results_list[[analysis_name]]$go_bp), 5)
          for (i in 1:nrow(top_go)) {
            cat(i, ". ", top_go$Description[i], " (p=", 
                format(top_go$pvalue[i], digits = 3), ")\n", sep = "")
          }
        }
        
        # Show top KEGG pathways
        if (!is.null(results_list[[analysis_name]]$kegg) && 
            nrow(results_list[[analysis_name]]$kegg@result) > 0) {
          cat("\nTop KEGG Pathways:\n")
          top_kegg <- head(as.data.frame(results_list[[analysis_name]]$kegg), 5)
          for (i in 1:nrow(top_kegg)) {
            cat(i, ". ", top_kegg$Description[i], " (p=", 
                format(top_kegg$pvalue[i], digits = 3), ")\n", sep = "")
          }
        }
        
        cat("\n")
      }
    }
    
    # Addiction-specific findings
    if (length(addiction_pathways) > 0) {
      cat("ADDICTION-RELATED PATHWAYS:\n\n")
      
      for (analysis_name in names(addiction_pathways)) {
        cat("=== ", toupper(analysis_name), " ===\n", sep = "")
        
        for (db_name in names(addiction_pathways[[analysis_name]])) {
          terms <- addiction_pathways[[analysis_name]][[db_name]]
          cat(db_name, ": ", nrow(terms), " addiction-related terms\n", sep = "")
          
          if (nrow(terms) > 0) {
            top_terms <- head(terms, 3)
            for (i in 1:nrow(top_terms)) {
              cat("  ", i, ". ", top_terms$Description[i], "\n", sep = "")
            }
          }
        }
        cat("\n")
      }
    }
  }
  
  cat("================================================================\n")
  cat("ANALYSIS COMPLETE\n")
  cat("================================================================\n")
  
  sink()
}

# Generate summary report
create_pathway_summary(enrichment_results, 
                      file.path(pathway_dirs$results, "pathway_enrichment_summary.txt"))

# Save summary table
if (length(enrichment_results) > 0) {
  summary_table <- data.frame(
    Analysis = names(enrichment_results),
    Total_Pathways = sapply(enrichment_results, function(x) {
      sum(sapply(x, function(y) ifelse(is.null(y), 0, nrow(y@result))))
    }),
    GO_Terms = sapply(enrichment_results, function(x) {
      sum(sapply(x[c("go_bp", "go_mf", "go_cc")], function(y) ifelse(is.null(y), 0, nrow(y@result))))
    }),
    KEGG_Pathways = sapply(enrichment_results, function(x) {
      ifelse(is.null(x$kegg), 0, nrow(x$kegg@result))
    }),
    stringsAsFactors = FALSE
  )
  
  write.csv(summary_table, file.path(pathway_dirs$tables, "enrichment_summary_table.csv"), row.names = FALSE)
}

# ----------------------
# 8. FINAL STATUS
# ----------------------

cat("\n")
cat("================================================================\n")
cat("PATHWAY ENRICHMENT ANALYSIS COMPLETE\n")
cat("================================================================\n")
cat("Analysis completed at:", as.character(Sys.time()), "\n\n")

cat("RESULTS GENERATED:\n")
cat("- Gene ID mapping tables for all lists\n")
cat("- Enrichment tables for", length(enrichment_results), "gene sets\n")
cat("- Visualization plots for all significant pathways\n")
if (length(addiction_pathways) > 0) {
  cat("- Addiction-specific pathway identification\n")
}
cat("- Comprehensive summary report\n\n")

cat("OUTPUT DIRECTORIES:\n")
cat("- Results:", pathway_dirs$results, "\n")
cat("- Tables:", pathway_dirs$tables, "\n")
cat("- Plots:", pathway_dirs$plots, "\n\n")

if (length(enrichment_results) > 0) {
  cat("KEY FILES TO EXAMINE:\n")
  cat("1. pathway_enrichment_summary.txt - Overall summary\n")
  cat("2. enrichment_summary_table.csv - Quantitative summary\n")
  cat("3. Gene ID mapping files in Tables/ - Check conversion success\n")
  cat("4. Individual analysis folders in Tables/\n")
  cat("5. Visualization plots in Plots/\n")
  if (length(addiction_pathways) > 0) {
    cat("6. Addiction-specific results in Results/Addiction_Specific/\n")
  }
} else {
  cat("NOTE: No pathway enrichment results found.\n")
  cat("Check gene ID mapping files for conversion issues.\n")
}

cat("\nREADY FOR:\n")
cat("- Biological interpretation of pathway results\n")
cat("- Literature comparison of identified pathways\n")
cat("- Target identification for validation\n")
cat("- Network analysis of pathway interactions\n")

cat("================================================================\n")