# ==============================================================================
# Comprehensive Neuroinflammatory Analysis
# ==============================================================================
# Purpose: Cross-dataset pathway enrichment analysis with Excel export
# Features: Multi-method enrichment, comprehensive Excel output, visualizations
# Datasets: GSE174409 (Bulk RNA-seq) + GSE225158 (snRNA-seq pseudobulk)
# Methods: Paired limma, Mixed effects, DESeq2 across GO, KEGG, Reactome, Hallmark
# ==============================================================================

# Configuration - ADD MISSING VARIABLES
BASE_DIR <- "/Users/aumchampaneri/Complement-OUD/Translational_Study/Homo_sapiens"
RESULTS_DIR <- file.path(BASE_DIR, "Results")
OUTPUTS_BASE <- file.path(BASE_DIR, "Outputs")
NEUROINFLAMM_DIR <- file.path(RESULTS_DIR, "Neuroinflammatory_Analysis")

# Load required libraries
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(msigdbr)
  library(dplyr)
  library(ggplot2)
  library(enrichplot)
  library(openxlsx)
  library(stringr)
})

# ==============================================================================
# SETUP AND CONFIGURATION
# ==============================================================================

#' Create organized output directories for neuroinflammatory analysis
create_output_directories <- function() {
  cat("=== Creating Organized Output Directories ===\n")
  
  # Main directories
  dirs_to_create <- c(
    NEUROINFLAMM_DIR,
    file.path(NEUROINFLAMM_DIR, "Excel_Exports"),
    file.path(NEUROINFLAMM_DIR, "Summary_Figures"),
    file.path(NEUROINFLAMM_DIR, "Summary_Figures", "Dataset_Overviews"),
    file.path(NEUROINFLAMM_DIR, "Individual_Plots"),
    file.path(NEUROINFLAMM_DIR, "Enrichment_CSV_Results")
  )
  
  # Create directories
  for (dir in dirs_to_create) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      cat("✓ Created:", dir, "\n")
    }
  }
  
  cat("✓ Output directories ready\n")
}

# ==============================================================================
# DATA LOADING FUNCTIONS
# ==============================================================================

#' Load both datasets for comprehensive analysis
#' @return List containing both datasets with their results and expression data
load_both_datasets <- function() {
  cat("=== Loading Both Datasets for Comprehensive Analysis ===\n")
  
  # Define file paths
  gse174409_file <- file.path(RESULTS_DIR, "GSE174409", "GSE174409_region_analysis_results.rds")
  gse225158_file <- file.path(RESULTS_DIR, "GSE225158", "GSE225158_region_analysis_results.rds")
  
  # Check if files exist
  if (!file.exists(gse174409_file)) {
    stop("GSE174409 results not found: ", gse174409_file, 
         "\nPlease run 00_GSE174409_Paired-Analysis.R first")
  }
  
  if (!file.exists(gse225158_file)) {
    stop("GSE225158 results not found: ", gse225158_file,
         "\nPlease run 00_GSE225158_Python-Loading.py and 00_GSE225158_Paired-Pseudobulk.R first")
  }
  
  # Load datasets
  gse174409_data <- readRDS(gse174409_file)
  gse225158_data <- readRDS(gse225158_file)
  
  cat("✓ GSE174409 loaded:", length(gse174409_data), "components\n")
  cat("✓ GSE225158 loaded:", length(gse225158_data), "components\n")
  
  # Validate that both datasets have the required structure
  required_methods <- c("paired_limma", "mixed_effects", "deseq2")
  
  # Check GSE174409 structure
  missing_174409 <- setdiff(required_methods, names(gse174409_data))
  if (length(missing_174409) > 0) {
    stop("GSE174409 missing methods: ", paste(missing_174409, collapse = ", "))
  }
  
  # Check GSE225158 structure
  missing_225158 <- setdiff(required_methods, names(gse225158_data))
  if (length(missing_225158) > 0) {
    stop("GSE225158 missing methods: ", paste(missing_225158, collapse = ", "))
  }
  
  # Verify expression data is available (critical for comprehensive analysis)
  if (!"expression_data" %in% names(gse174409_data)) {
    cat("⚠ GSE174409 missing expression_data - some analyses may be limited\n")
  }
  
  if (!"expression_data" %in% names(gse225158_data)) {
    cat("⚠ GSE225158 missing expression_data - some analyses may be limited\n")
  }
  
  # Return both datasets
  datasets <- list(
    GSE174409 = gse174409_data,
    GSE225158 = gse225158_data
  )
  
  cat("✓ Both datasets successfully loaded and validated\n")
  return(datasets)
}

#' Extract significant genes from analysis results
#' @param results_obj Analysis results object from limma or DESeq2
#' @param method Method name (for column selection)
#' @param fdr_threshold FDR threshold for significance
#' @return Vector of significant gene names
extract_significant_genes <- function(results_obj, method = "limma", fdr_threshold = 0.05) {
  results_df <- results_obj$results
  
  # Determine p-value column based on method
  if (method %in% c("paired_limma", "mixed_effects")) {
    pval_col <- "adj.P.Val"
  } else if (method == "deseq2") {
    pval_col <- "padj"
  } else {
    # Try to detect automatically
    if ("adj.P.Val" %in% colnames(results_df)) {
      pval_col <- "adj.P.Val"
    } else if ("padj" %in% colnames(results_df)) {
      pval_col <- "padj"
    } else {
      stop("Cannot determine p-value column for method: ", method)
    }
  }
  
  # Extract significant genes
  sig_mask <- results_df[[pval_col]] < fdr_threshold & !is.na(results_df[[pval_col]])
  sig_genes <- rownames(results_df)[sig_mask]
  
  return(sig_genes)
}

# ==============================================================================
# ENRICHMENT ANALYSIS FUNCTIONS
# ==============================================================================

#' Run enrichment analysis for both datasets
#' @param datasets List containing both GSE174409 and GSE225158 data
#' @return Comprehensive results object with all enrichment analyses
run_enrichment_analysis <- function(datasets) {
  cat("\n=== Running Enrichment Analysis ===\n")
  
  # Initialize results structure
  enrichment_results <- list()
  
  # Analysis methods to process
  methods <- c("paired_limma", "mixed_effects", "deseq2")
  
  # Pathway databases
  databases <- c("GO_BP", "KEGG", "Reactome", "Hallmark")
  
  # Process each dataset
  for (dataset_name in names(datasets)) {
    cat("Processing", dataset_name, "...\n")
    
    dataset <- datasets[[dataset_name]]
    enrichment_results[[dataset_name]] <- list()
    
    # Process each method
    for (method in methods) {
      cat("  Running", method, "enrichment...\n")
      
      # Extract significant genes
      sig_genes <- extract_significant_genes(dataset[[method]], method)
      
      if (length(sig_genes) < 10) {
        cat("    ⚠ Too few significant genes (", length(sig_genes), ") for", method, "\n")
        next
      }
      
      enrichment_results[[dataset_name]][[method]] <- list()
      
      # Run enrichment for each database
      for (db in databases) {
        tryCatch({
          enrichment_result <- run_pathway_enrichment(sig_genes, database = db, dataset_name = dataset_name, method = method)
          enrichment_results[[dataset_name]][[method]][[db]] <- enrichment_result
        }, error = function(e) {
          cat("    ⚠ Enrichment failed for", db, ":", e$message, "\n")
          enrichment_results[[dataset_name]][[method]][[db]] <- NULL
        })
      }
    }
  }
  
  cat("✓ Enrichment analysis complete\n")
  
  # Create comprehensive results object
  comprehensive_results <- list(
    enrichment_results = enrichment_results,
    datasets = datasets,
    analysis_date = Sys.Date(),
    methods_analyzed = methods,
    databases_analyzed = databases
  )
  
  return(comprehensive_results)
}

#' Run pathway enrichment analysis for a set of genes
#' @param genes Vector of gene symbols
#' @param database Pathway database ("GO_BP", "KEGG", "Reactome", "Hallmark")
#' @param dataset_name Dataset name for logging
#' @param method Method name for logging
#' @return enrichResult object or NULL if failed
run_pathway_enrichment <- function(genes, database = "GO_BP", dataset_name = "", method = "") {
  
  # Convert gene symbols to Entrez IDs for most databases
  if (database %in% c("GO_BP", "KEGG", "Reactome")) {
    gene_mapping <- clusterProfiler::bitr(genes, 
                                          fromType = "SYMBOL", 
                                          toType = "ENTREZID", 
                                          OrgDb = org.Hs.eg.db)
    
    if (nrow(gene_mapping) == 0) {
      stop("No genes could be mapped to Entrez IDs")
    }
    
    entrez_genes <- gene_mapping$ENTREZID
    cat(sprintf("    Mapped %d genes to Entrez IDs\n", length(entrez_genes)))
  } else {
    # For databases that work with gene symbols
    entrez_genes <- genes
  }
  
  # Run enrichment based on database
  if (database == "GO_BP") {
    result <- clusterProfiler::enrichGO(
      gene = entrez_genes,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    )
  } else if (database == "KEGG") {
    result <- clusterProfiler::enrichKEGG(
      gene = entrez_genes,
      organism = "hsa",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
  } else if (database == "Reactome") {
    result <- ReactomePA::enrichPathway(
      gene = entrez_genes,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    )
  } else if (database == "Hallmark") {
    # FIXED: Use 'collection' instead of deprecated 'category'
    hallmark_sets <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
    hallmark_list <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)
    
    result <- clusterProfiler::enricher(
      gene = genes,  # Use symbols for Hallmark
      TERM2GENE = hallmark_sets[, c("gs_name", "gene_symbol")],
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
  } else {
    stop("Unsupported database: ", database)
  }
  
  # Log results
  if (!is.null(result) && nrow(result@result) > 0) {
    total_terms <- nrow(result@result)
    sig_terms <- sum(result@result$p.adjust < 0.05)
    cat(sprintf("    %s: %d enriched terms\n", database, total_terms))
  } else {
    cat(sprintf("    %s: No enrichment found\n", database))
  }
  
  return(result)
}

# ==============================================================================
# VISUALIZATION FUNCTIONS
# ==============================================================================

#' Create comparative visualizations across datasets and methods
#' @param comprehensive_results Results from run_enrichment_analysis
create_comparative_visualizations <- function(comprehensive_results) {
  cat("\n=== Creating Comparative Visualizations ===\n")
  
  # Create dataset-specific plots
  create_dataset_specific_plots(comprehensive_results)
  
  # Create cross-dataset summary plots
  create_summary_plots(comprehensive_results)
  
  cat("✓ Comparative visualizations complete\n")
}

#' Create individual plots for each dataset-method-database combination
create_dataset_specific_plots <- function(comprehensive_results) {
  cat("\n=== Creating Dataset-Specific Enrichment Plots ===\n")
  
  for (dataset in names(comprehensive_results$enrichment_results)) {
    for (method in names(comprehensive_results$enrichment_results[[dataset]])) {
      for (db in names(comprehensive_results$enrichment_results[[dataset]][[method]])) {
        
        result_obj <- comprehensive_results$enrichment_results[[dataset]][[method]][[db]]
        
        if (!is.null(result_obj) && nrow(result_obj@result) > 0) {
          
          # Check if there are significant results
          sig_results <- result_obj@result[result_obj@result$p.adjust < 0.05, ]
          
          if (nrow(sig_results) > 0) {
            tryCatch({
              # Create dotplot
              p <- enrichplot::dotplot(result_obj, showCategory = 15) +
                ggtitle(paste(dataset, method, db, "Enrichment")) +
                theme_minimal()
              
              # Save plot
              plot_dir <- file.path(NEUROINFLAMM_DIR, "Individual_Plots", dataset, method)
              if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
              
              plot_file <- file.path(plot_dir, paste0(db, "_dotplot.png"))
              ggplot2::ggsave(plot_file, p, width = 10, height = 8, dpi = 300)
              
              cat("✓", db, "plots saved for", dataset, method, "\n")
              
            }, error = function(e) {
              cat("⚠ Plot creation failed for", dataset, method, db, ":", e$message, "\n")
            })
          } else {
            cat("⚠ No significant", db, "terms for", dataset, method, "\n")
          }
        }
      }
    }
  }
}

#' Create summary plots across all analyses
create_summary_plots <- function(comprehensive_results) {
  cat("Creating enrichment summary plots...\n")
  
  # Create summary data
  summary_data <- data.frame()
  
  for (dataset in names(comprehensive_results$enrichment_results)) {
    for (method in names(comprehensive_results$enrichment_results[[dataset]])) {
      for (db in names(comprehensive_results$enrichment_results[[dataset]][[method]])) {
        
        result_obj <- comprehensive_results$enrichment_results[[dataset]][[method]][[db]]
        
        if (!is.null(result_obj) && nrow(result_obj@result) > 0) {
          total_terms <- nrow(result_obj@result)
          sig_terms <- sum(result_obj@result$p.adjust < 0.05, na.rm = TRUE)
          
          summary_data <- rbind(summary_data, data.frame(
            Dataset = dataset,
            Method = method,
            Database = db,
            Total_Terms = total_terms,
            Significant_Terms = sig_terms,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  if (nrow(summary_data) > 0) {
    # Create summary plot
    p_summary <- ggplot(summary_data, aes(x = Database, y = Significant_Terms, fill = Dataset)) +
      geom_col(position = "dodge", alpha = 0.8) +
      facet_wrap(~ Method, scales = "free_y") +
      scale_fill_manual(values = c("GSE174409" = "#FF6B6B", "GSE225158" = "#4ECDC4")) +
      labs(title = "Pathway Enrichment Summary",
           subtitle = "Significant pathways (FDR < 0.05) by dataset, method, and database",
           x = "Pathway Database",
           y = "Number of Significant Pathways") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Save summary plot
    summary_dir <- file.path(NEUROINFLAMM_DIR, "Summary_Figures", "Dataset_Overviews")
    if (!dir.exists(summary_dir)) dir.create(summary_dir, recursive = TRUE)
    
    summary_file <- file.path(summary_dir, "enrichment_summary.png")
    ggplot2::ggsave(summary_file, p_summary, width = 14, height = 10, dpi = 300)
    
    cat("✓ Summary plot saved:", summary_file, "\n")
  }
}

# ==============================================================================
# EXCEL EXPORT FUNCTIONS (ENHANCED)
# ==============================================================================

#' Export comprehensive enrichment results to Excel with multiple sheets (FIXED)
export_enrichment_excel <- function(comprehensive_results, output_dir) {
  cat("\n=== Exporting Enrichment Results to Excel ===\n")
  
  # Create Excel file
  wb <- createWorkbook()
  excel_file <- file.path(output_dir, "Comprehensive_Pathway_Enrichment_Results.xlsx")
  
  # Track all data for summary sheets
  all_significant_pathways <- list()  # FIXED: Use list instead of data.frame
  summary_data <- data.frame()
  
  # Process each dataset-method-database combination
  for (dataset in names(comprehensive_results$enrichment_results)) {
    for (method in names(comprehensive_results$enrichment_results[[dataset]])) {
      for (db in names(comprehensive_results$enrichment_results[[dataset]][[method]])) {
        
        result_obj <- comprehensive_results$enrichment_results[[dataset]][[method]][[db]]
        
        if (!is.null(result_obj) && nrow(result_obj@result) > 0) {
          
          # Get enrichment data and STANDARDIZE columns
          enrichment_df <- standardize_enrichment_columns(result_obj@result, db)
          enrichment_df$Dataset <- dataset
          enrichment_df$Method <- method
          enrichment_df$Database <- db
          
          # Create sheet name (Excel has 31 character limit)
          sheet_name <- paste0(dataset, "_", method, "_", substring(db, 1, 10))
          if (nchar(sheet_name) > 31) {
            sheet_name <- substring(sheet_name, 1, 31)
          }
          
          # Add sheet
          addWorksheet(wb, sheet_name)
          writeData(wb, sheet_name, enrichment_df)
          
          # Track summary stats
          total_pathways <- nrow(enrichment_df)
          sig_pathways <- sum(enrichment_df$p.adjust < 0.05, na.rm = TRUE)
          
          summary_data <- rbind(summary_data, data.frame(
            Dataset = dataset,
            Method = method,
            Database = db,
            Total_Pathways = total_pathways,
            Significant_Pathways = sig_pathways,
            Sheet_Name = sheet_name,
            stringsAsFactors = FALSE
          ))
          
          cat(sprintf("✓ Processed: %s (%d pathways, %d significant)\n", 
                     sheet_name, total_pathways, sig_pathways))
          
          # Add significant pathways to master list (FIXED)
          if (sig_pathways > 0) {
            sig_data <- enrichment_df[enrichment_df$p.adjust < 0.05, ]
            # Store each significant result separately to avoid column mismatch
            all_significant_pathways[[paste(dataset, method, db, sep = "_")]] <- sig_data
          }
        }
      }
    }
  }
  
  # Add summary sheet
  addWorksheet(wb, "Summary", tabColour = "red")
  writeData(wb, "Summary", summary_data)
  cat("✓ Summary sheet created\n")
  
  # FIXED: Combine significant pathways safely
  if (length(all_significant_pathways) > 0) {
    # Find common columns across all significant results
    all_columns <- unique(unlist(lapply(all_significant_pathways, colnames)))
    
    # Standardize all data frames to have the same columns
    standardized_pathways <- lapply(all_significant_pathways, function(df) {
      missing_cols <- setdiff(all_columns, colnames(df))
      for (col in missing_cols) {
        df[[col]] <- NA  # Add missing columns as NA
      }
      return(df[, all_columns])  # Reorder columns consistently
    })
    
    # Now rbind should work
    combined_significant <- do.call(rbind, standardized_pathways)
    
    # Add all significant pathways sheet
    addWorksheet(wb, "All_Significant_Pathways", tabColour = "green")
    writeData(wb, "All_Significant_Pathways", combined_significant)
    cat("✓ All significant pathways sheet:", nrow(combined_significant), "pathways\n")
    
    # Add top 100 pathways sheet
    if (nrow(combined_significant) > 0) {
      top_pathways <- combined_significant %>%
        arrange(p.adjust) %>%
        head(100)
      
      addWorksheet(wb, "Top_100_Pathways", tabColour = "blue")
      writeData(wb, "Top_100_Pathways", top_pathways)
      cat("✓ Top 100 pathways sheet created\n")
      
      # Create database-specific sheets
      for (db in unique(combined_significant$Database)) {
        db_data <- combined_significant[combined_significant$Database == db, ]
        if (nrow(db_data) > 0) {
          sheet_name <- paste0("Significant_", db)
          addWorksheet(wb, sheet_name)
          writeData(wb, sheet_name, db_data)
          cat(sprintf("✓ %s sheet: %d pathways\n", sheet_name, nrow(db_data)))
        }
      }
      
      # Create complement-specific sheet
      complement_keywords <- c("complement", "C1Q", "C3", "C5", "classical complement")
      complement_pathways <- combined_significant[
        grepl(paste(complement_keywords, collapse = "|"), 
              combined_significant$Description, ignore.case = TRUE), ]
      
      if (nrow(complement_pathways) > 0) {
        addWorksheet(wb, "Complement_Pathways", tabColour = "orange")
        writeData(wb, "Complement_Pathways", complement_pathways)
        cat(sprintf("✓ Complement pathways sheet: %d pathways\n", nrow(complement_pathways)))
      }
    }
  }
  
  # Save Excel file
  saveWorkbook(wb, excel_file, overwrite = TRUE)
  
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("📊 COMPREHENSIVE EXCEL EXPORT COMPLETE!\n")
  cat("📁 File saved:", excel_file, "\n")
  cat("📋 Sheets created:\n")
  for (i in 1:length(names(wb))) {
    cat(sprintf("  %d. %s\n", i, names(wb)[i]))
  }
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  return(excel_file)
}

#' Standardize enrichment result columns for consistent processing (ENHANCED)
#' @param enrichment_df Enrichment results data frame
#' @param database Database name for context
#' @return Standardized data frame with consistent columns
standardize_enrichment_columns <- function(enrichment_df, database) {
  # Define core columns that should always be present
  core_cols <- c("ID", "Description", "pvalue", "p.adjust")
  
  # Define optional columns that may vary by database
  optional_cols <- c("qvalue", "GeneRatio", "BgRatio", "Count", "geneID", "FoldEnrichment")
  
  # Create base standardized data frame with core columns
  standardized_df <- data.frame(
    ID = ifelse("ID" %in% colnames(enrichment_df), 
                enrichment_df$ID, 
                paste0(database, "_", seq_len(nrow(enrichment_df)))),
    Description = enrichment_df$Description,
    pvalue = enrichment_df$pvalue,
    p.adjust = enrichment_df$p.adjust,
    stringsAsFactors = FALSE
  )
  
  # Add optional columns if they exist, otherwise fill with appropriate defaults
  for (col in optional_cols) {
    if (col %in% colnames(enrichment_df)) {
      standardized_df[[col]] <- enrichment_df[[col]]
    } else {
      # Provide appropriate defaults based on column type
      if (col == "qvalue") {
        standardized_df[[col]] <- enrichment_df$p.adjust  # Use p.adjust as fallback
      } else if (col == "Count" && "GeneRatio" %in% colnames(enrichment_df)) {
        # Extract count from GeneRatio if Count is missing
        gene_ratios <- as.character(enrichment_df$GeneRatio)
        counts <- sapply(strsplit(gene_ratios, "/"), function(x) {
          if (length(x) >= 1) as.numeric(x[1]) else NA
        })
        standardized_df[[col]] <- counts
      } else if (col %in% c("GeneRatio", "BgRatio", "geneID")) {
        standardized_df[[col]] <- ""  # Empty string for character columns
      } else {
        standardized_df[[col]] <- NA  # NA for numeric columns
      }
    }
  }
  
  return(standardized_df)
}

# ==============================================================================
# MAIN PIPELINE
# ==============================================================================

#' Run comprehensive neuroinflammatory analysis with Excel export (FIXED)
run_comprehensive_neuroinflammatory_analysis <- function() {
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("COMPREHENSIVE NEUROINFLAMMATORY ANALYSIS\n")
  cat("Cross-dataset pathway enrichment with Excel export\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  # Create organized output directories
  create_output_directories()
  
  # Create analysis results directory (separate from figures)
  if (!dir.exists(NEUROINFLAMM_DIR)) dir.create(NEUROINFLAMM_DIR, recursive = TRUE)
  
  tryCatch({
    # Load both datasets
    datasets <- load_both_datasets()
    
    # Run enrichment analysis
    comprehensive_results <- run_enrichment_analysis(datasets)
    
    # Create visualizations
    create_comparative_visualizations(comprehensive_results)
    
    # Export to Excel file
    excel_file <- export_enrichment_excel(comprehensive_results, NEUROINFLAMM_DIR)
    
    # FIXED: Export CSV with proper error handling
    csv_dir <- tryCatch({
      export_enrichment_csvs(comprehensive_results, NEUROINFLAMM_DIR)
    }, error = function(e) {
      cat("⚠ CSV export failed:", e$message, "\n")
      cat("Excel export succeeded, CSV export skipped\n")
      NULL
    })
    
    # Save comprehensive results
    results_file <- file.path(NEUROINFLAMM_DIR, "comprehensive_neuroinflammatory_analysis_enhanced.rds")
    saveRDS(comprehensive_results, results_file)
    
    cat("\n", paste(rep("=", 70), collapse = ""), "\n")
    cat("SUCCESS: Comprehensive analysis complete!\n")
    cat("✓ Enrichment results saved as RDS:", results_file, "\n")
    cat("✓ All pathway results exported to Excel:", excel_file, "\n")
    
    if (!is.null(csv_dir)) {
      cat("✓ CSV files also available:", csv_dir, "\n")
    } else {
      cat("⚠ CSV export skipped due to column mismatch - Excel export complete\n")
    }
    
    cat("✓ Ready for dashboard creation\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")
    
    return(list(
      results = comprehensive_results,
      excel_file = excel_file,
      csv_directory = csv_dir
    ))
    
  }, error = function(e) {
    cat("\n❌ ANALYSIS ERROR:", e$message, "\n")
    stop(e)
  })
}

#' Export all enrichment results to CSV files for easy access
export_enrichment_csvs <- function(comprehensive_results, output_dir) {
  cat("\n=== Exporting Enrichment Results to CSV ===\n")
  
  # Create CSV output directory
  csv_dir <- file.path(output_dir, "Enrichment_CSV_Results")
  if (!dir.exists(csv_dir)) dir.create(csv_dir, recursive = TRUE)
  
  # Initialize summary data frame
  all_enrichment_summary <- data.frame()
  
  # Process each dataset and method
  for (dataset in names(comprehensive_results$enrichment_results)) {
    for (method in names(comprehensive_results$enrichment_results[[dataset]])) {
      for (db in names(comprehensive_results$enrichment_results[[dataset]][[method]])) {
        
        result_obj <- comprehensive_results$enrichment_results[[dataset]][[method]][[db]]
        
        if (!is.null(result_obj) && nrow(result_obj@result) > 0) {
          # Standardize columns before processing
          enrichment_df <- standardize_enrichment_columns(result_obj@result, db)
          
          # Add metadata columns
          enrichment_df$Dataset <- dataset
          enrichment_df$Method <- method
          enrichment_df$Database <- db
          enrichment_df$Analysis_ID <- paste(dataset, method, db, sep = "_")
          
          # Save individual CSV file
          csv_filename <- paste0(dataset, "_", method, "_", db, "_enrichment.csv")
          csv_filepath <- file.path(csv_dir, csv_filename)
          write.csv(enrichment_df, csv_filepath, row.names = FALSE)
          
          cat("✓ Saved:", csv_filename, "\n")
          
          # Add to master summary (significant pathways only)
          if (any(enrichment_df$p.adjust < 0.05, na.rm = TRUE)) {
            sig_pathways <- enrichment_df[enrichment_df$p.adjust < 0.05 & !is.na(enrichment_df$p.adjust), ]
            
            if (nrow(sig_pathways) > 0) {
              # Ensure consistent column structure before rbind
              if (nrow(all_enrichment_summary) == 0) {
                all_enrichment_summary <- sig_pathways
              } else {
                # Find common columns and use only those
                common_cols <- intersect(colnames(all_enrichment_summary), colnames(sig_pathways))
                all_enrichment_summary <- rbind(
                  all_enrichment_summary[, common_cols, drop = FALSE],
                  sig_pathways[, common_cols, drop = FALSE]
                )
              }
            }
          }
        }
      }
    }
  }
  
  # Save master summary of all significant pathways
  if (nrow(all_enrichment_summary) > 0) {
    master_file <- file.path(csv_dir, "ALL_Significant_Pathways_Master.csv")
    write.csv(all_enrichment_summary, master_file, row.names = FALSE)
    cat("✓ Master summary saved:", basename(master_file), "\n")
  }
  
  cat("📁 All CSV files saved to:", csv_dir, "\n")
  return(csv_dir)
}

# ==============================================================================
# EXECUTION
# ==============================================================================

if (!exists("SOURCED")) {
  comprehensive_analysis_results <- run_comprehensive_neuroinflammatory_analysis()
}
