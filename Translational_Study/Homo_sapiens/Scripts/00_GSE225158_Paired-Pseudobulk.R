# ==============================================================================
# GSE225158 Paired Pseudobulk Analysis
# ==============================================================================
# Purpose: Paired differential expression analysis on pseudobulk data from snRNA-seq
# Dataset: GSE225158 - Striatal single-nucleus RNA-seq (Caudate + Putamen)
# Design:  10 paired subjects (OUD vs CTL) × 2 brain regions
# Methods: Paired limma-voom, Mixed effects, DESeq2 (identical to GSE174409)
# Input:   CSV files created by 00_GSE225158_Python-Loading.py
# ==============================================================================

# Configuration
OUTPUT_DIR <- "/Users/aumchampaneri/Complement-OUD/Translational_Study/Homo_sapiens/Results/GSE225158"
COUNTS_FILE <- file.path(OUTPUT_DIR, "GSE225158_pseudobulk_counts.csv")
METADATA_FILE <- file.path(OUTPUT_DIR, "GSE225158_pseudobulk_metadata.csv")

# NEW: Organized figure outputs
OUTPUTS_BASE <- "/Users/aumchampaneri/Complement-OUD/Translational_Study/Homo_sapiens/Outputs"
GSE225158_OUTPUTS <- file.path(OUTPUTS_BASE, "GSE225158_Analysis")

# Load required libraries
suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)  # Add this
  library(org.Hs.eg.db)     # Add this
  library(tibble)           # Add this
})

# ==============================================================================
# DATA LOADING AND VALIDATION
# ==============================================================================

#' Load and validate GSE225158 pseudobulk data
#' @return List containing counts matrix and metadata
load_pseudobulk_data <- function() {
  cat("=== Loading GSE225158 Pseudobulk Data ===\n")
  
  # Check input files exist
  if (!file.exists(COUNTS_FILE)) {
    stop("Counts file not found: ", COUNTS_FILE, 
         "\nPlease run 00_GSE225158_Python-Loading.py first")
  }
  if (!file.exists(METADATA_FILE)) {
    stop("Metadata file not found: ", METADATA_FILE,
         "\nPlease run 00_GSE225158_Python-Loading.py first")
  }
  
  # Load data
  counts <- read.csv(COUNTS_FILE, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  metadata <- read.csv(METADATA_FILE, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  
  # Quick verification of gene names
  cat("Gene names loaded correctly:", !all(grepl("^[0-9]+$", rownames(counts)[1:10])), "\n")
  cat("Sample gene names:", paste(rownames(counts)[1:5], collapse = ", "), "\n")
  
  # Fix R's automatic "X" prefix for numeric column names
  colnames(counts) <- gsub("^X", "", colnames(counts))
  
  # Align samples
  common_samples <- intersect(colnames(counts), rownames(metadata))
  counts <- counts[, common_samples]
  metadata <- metadata[common_samples, ]
  
  # Convert to factors
  metadata$Region <- factor(metadata$Region)
  metadata$Sex <- factor(metadata$Sex) 
  metadata$OUD_Status <- factor(metadata$OUD_Status)
  metadata$Subject_ID <- factor(metadata$Subject_ID)
  
  # Validate paired design
  paired_subjects <- metadata %>%
    group_by(Subject_ID) %>%
    summarise(n_regions = n_distinct(Region), .groups = 'drop') %>%
    filter(n_regions == 2) %>%
    pull(Subject_ID)
  
  if (length(paired_subjects) < 3) {
    stop("Insufficient paired subjects found: ", length(paired_subjects), 
         "\nNeed at least 3 for analysis")
  }
  
  # Filter to paired subjects and sort
  metadata <- metadata %>%
    filter(Subject_ID %in% paired_subjects) %>%
    arrange(Subject_ID, Region)
  
  counts <- counts[, rownames(metadata)]
  
  # Print summary
  cat("Loaded:", nrow(counts), "genes ×", ncol(counts), "samples\n")
  cat("Paired subjects:", length(paired_subjects), "\n")
  cat("Design:", paste(names(table(metadata$Region)), collapse = " vs "), "\n")
  cat("Groups:", paste(names(table(metadata$OUD_Status)), collapse = " vs "), "\n")
  
  # CRITICAL: Verify gene identifier format
  sample_genes <- rownames(counts)[1:10]
  cat("Sample gene identifiers:", paste(sample_genes, collapse = ", "), "\n")
  
  # Check if gene names look like Ensembl IDs
  if (any(grepl("^ENSG", sample_genes))) {
    cat("⚠ Detected Ensembl IDs - converting to gene symbols for consistency\n")
    
    tryCatch({
      # Convert Ensembl to gene symbols
      gene_mapping <- clusterProfiler::bitr(rownames(counts), 
                                            fromType = "ENSEMBL", 
                                            toType = "SYMBOL", 
                                            OrgDb = org.Hs.eg.db)
      
      cat("Mapped", nrow(gene_mapping), "out of", nrow(counts), "genes\n")
      
      # Keep only mapped genes
      counts_mapped <- counts[gene_mapping$ENSEMBL, ]
      
      # Handle duplicates using base R approach
      if (any(duplicated(gene_mapping$SYMBOL))) {
        cat("Handling", sum(duplicated(gene_mapping$SYMBOL)), "duplicate gene symbols\n")
        
        # Simple approach: aggregate using base R
        ensembl_to_symbol <- setNames(gene_mapping$SYMBOL, gene_mapping$ENSEMBL)
        gene_symbols <- ensembl_to_symbol[rownames(counts_mapped)]
        
        # Create data frame for aggregation
        temp_df <- data.frame(
          gene_symbol = gene_symbols,
          counts_mapped,
          stringsAsFactors = FALSE
        )
        
        # Aggregate by gene symbol
        counts_agg <- aggregate(temp_df[, -1], by = list(gene_symbol = temp_df$gene_symbol), FUN = sum)
        
        # Convert back to matrix
        rownames(counts_agg) <- counts_agg$gene_symbol
        counts_agg$gene_symbol <- NULL
        counts <- as.matrix(counts_agg)
        
      } else {
        # No duplicates, simple assignment
        rownames(counts_mapped) <- gene_mapping$SYMBOL
        counts <- as.matrix(counts_mapped)
      }
      
      cat("✓ Gene symbol conversion complete:", nrow(counts), "genes retained\n")
      
    }, error = function(e) {
      cat("Warning: Gene conversion failed:", e$message, "\n")
    })
  } else {
    cat("✓ Gene identifiers appear to be symbols already\n")
  }
  
  return(list(counts = counts, metadata = metadata))
}

# ==============================================================================
# STATISTICAL ANALYSIS METHODS
# ==============================================================================

#' Method 1: Paired limma-voom analysis with duplicate correlation
paired_limma_analysis <- function(counts, metadata) {
  cat("\n=== Method 1: Paired Limma-Voom ===\n")
  
  # Create DGEList and filter low-expressed genes
  dge <- DGEList(counts = counts)
  keep <- filterByExpr(dge, group = metadata$Region)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  cat("Genes after filtering:", nrow(dge), "\n")
  
  # Design matrix
  design <- model.matrix(~ Region + Sex + OUD_Status, data = metadata)
  
  # Voom transformation with duplicate correlation
  v <- voom(dge, design, plot = FALSE)
  corfit <- duplicateCorrelation(v, design, block = metadata$Subject_ID)
  cat("Intra-subject correlation:", round(corfit$consensus, 3), "\n")
  
  # Fit linear model
  fit <- lmFit(v, design, block = metadata$Subject_ID, correlation = corfit$consensus)
  fit <- eBayes(fit)
  
  # Extract region effect results
  region_coef <- grep("Region", colnames(design), value = TRUE)[1]
  results <- topTable(fit, coef = region_coef, number = Inf, sort.by = "P")
  results$Method <- "Paired_Limma"
  
  return(list(results = results, fit = fit, correlation = corfit$consensus))
}

#' Method 2: Mixed effects analysis with quality weights
mixed_effects_analysis <- function(counts, metadata) {
  cat("\n=== Method 2: Mixed Effects ===\n")
  
  # Create DGEList and filter
  dge <- DGEList(counts = counts)
  keep <- filterByExpr(dge, group = metadata$Region)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  cat("Genes after filtering:", nrow(dge), "\n")
  
  # Design matrix
  design <- model.matrix(~ Region + Sex + OUD_Status, data = metadata)
  
  # Voom with quality weights (iterative)
  v <- voomWithQualityWeights(dge, design, plot = FALSE)
  corfit <- duplicateCorrelation(v, design, block = metadata$Subject_ID)
  
  v <- voomWithQualityWeights(dge, design, plot = FALSE, 
                              block = metadata$Subject_ID, 
                              correlation = corfit$consensus)
  corfit <- duplicateCorrelation(v, design, block = metadata$Subject_ID)
  cat("Intra-subject correlation:", round(corfit$consensus, 3), "\n")
  
  # Fit model
  fit <- lmFit(v, design, block = metadata$Subject_ID, correlation = corfit$consensus)
  fit <- eBayes(fit)
  
  # Extract results
  region_coef <- grep("Region", colnames(design), value = TRUE)[1]
  results <- topTable(fit, coef = region_coef, number = Inf, sort.by = "P")
  results$Method <- "Mixed_Effects"
  
  return(list(results = results, fit = fit, correlation = corfit$consensus))
}

#' Method 3: DESeq2 analysis
deseq2_analysis <- function(counts, metadata) {
  cat("\n=== Method 3: DESeq2 ===\n")
  
  # Prepare integer count matrix
  counts_int <- round(as.matrix(counts))
  design_formula <- ~ Region + Sex + OUD_Status
  
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = counts_int,
    colData = metadata,
    design = design_formula
  )
  
  # Filter low counts
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  cat("Genes after filtering:", nrow(dds), "\n")
  
  # Run DESeq2
  dds <- DESeq(dds, quiet = TRUE)
  
  # Extract region comparison results
  region_levels <- levels(metadata$Region)
  results_obj <- results(dds, contrast = c("Region", region_levels[2], region_levels[1]))
  results_df <- as.data.frame(results_obj)
  results_df$Method <- "DESeq2"
  
  return(list(results = results_df, dds = dds))
}

# ==============================================================================
# RESULTS OUTPUT AND SUMMARY
# ==============================================================================

#' Save analysis results and create summary
save_results <- function(results_list, metadata) {
  cat("\n=== Saving Results ===\n")
  
  # Create output directory
  if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
  
  # Save main results object (CRITICAL for integration)
  results_file <- file.path(OUTPUT_DIR, "GSE225158_region_analysis_results.rds")
  saveRDS(results_list, results_file)
  cat("✓ Main results saved:", results_file, "\n")
  
  # Save individual method results as CSV
  for (method in names(results_list)) {
    method_file <- file.path(OUTPUT_DIR, paste0("GSE225158_", method, "_results.csv"))
    write.csv(results_list[[method]]$results, method_file, row.names = TRUE)
    cat("✓", method, "results:", method_file, "\n")
  }
  
  # Create method comparison summary
  summary_df <- data.frame(
    Method = names(results_list),
    N_Genes_Tested = sapply(results_list, function(x) nrow(x$results)),
    N_Significant_FDR05 = sapply(results_list, function(x) {
      res <- x$results
      pval_col <- ifelse("adj.P.Val" %in% colnames(res), "adj.P.Val", "padj")
      sum(res[[pval_col]] < 0.05, na.rm = TRUE)
    }),
    Intra_Subject_Correlation = sapply(results_list, function(x) {
      if ("correlation" %in% names(x)) round(x$correlation, 3) else NA
    }),
    N_Subjects = nrow(metadata) / 2,
    N_Samples = nrow(metadata)
  )
  
  summary_file <- file.path(OUTPUT_DIR, "GSE225158_method_comparison.csv")
  write.csv(summary_df, summary_file, row.names = FALSE)
  cat("✓ Summary saved:", summary_file, "\n")
  
  # Verify the RDS file was created properly
  if (file.exists(results_file)) {
    cat("✓ RDS file confirmed for integration analysis\n")
  } else {
    warning("RDS file not created - integration analysis may fail")
  }
  
  return(summary_df)
}

# ==============================================================================
# MAIN ANALYSIS PIPELINE
# ==============================================================================

#' Run complete GSE225158 paired pseudobulk analysis
#' @return List of analysis results from all three methods
run_gse225158_analysis <- function() {
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("GSE225158 PAIRED PSEUDOBULK ANALYSIS\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  tryCatch({
    # Load and validate data
    data <- load_pseudobulk_data()
    
    # Run three analysis methods (same as GSE174409)
    results_list <- list(
      paired_limma = paired_limma_analysis(data$counts, data$metadata),
      mixed_effects = mixed_effects_analysis(data$counts, data$metadata),
      deseq2 = deseq2_analysis(data$counts, data$metadata)
    )
    
    # Print analysis summary
    cat("\n=== Analysis Summary ===\n")
    for (method in names(results_list)) {
      res <- results_list[[method]]$results
      pval_col <- ifelse("adj.P.Val" %in% colnames(res), "adj.P.Val", "padj")
      sig_genes <- sum(res[[pval_col]] < 0.05, na.rm = TRUE)
      cat(sprintf("%-15s: %4d significant genes (FDR < 0.05)\n", method, sig_genes))
    }
    
    # Save results
    summary_df <- save_results(results_list, data$metadata)
    
    cat("\n", paste(rep("=", 70), collapse = ""), "\n")
    cat("SUCCESS: Analysis complete! Results saved to:", OUTPUT_DIR, "\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")
    
    return(results_list)
    
  }, error = function(e) {
    cat("\nERROR:", e$message, "\n")
    stop(e)
  })
}

# ==============================================================================
# EXECUTION
# ==============================================================================

# Run analysis when script is sourced
if (!exists("SOURCED")) {
  results <- run_gse225158_analysis()
}

# ==============================================================================
# ENHANCED DATA PREPARATION FOR COMPREHENSIVE ANALYSIS
# ==============================================================================

#' Prepare and save expression matrices for downstream analysis
prepare_expression_matrices <- function(counts, metadata) {
  cat("\n=== Preparing Expression Matrices for Comprehensive Analysis ===\n")
  
  # Create DGEList for normalization
  dge <- DGEList(counts = counts)
  keep <- filterByExpr(dge, group = metadata$Region)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  
  cat("Genes after filtering:", nrow(dge), "\n")
  
  # Log2-CPM transformation for expression matrix
  log_cpm <- cpm(dge, log = TRUE, prior.count = 1)
  
  # Voom-transformed values (for compatibility with limma downstream)
  design <- model.matrix(~ Region + Sex + OUD_Status, data = metadata)
  v <- voom(dge, design, plot = FALSE)
  voom_expression <- v$E
  
  # Raw counts (filtered)
  filtered_counts <- dge$counts
  
  expression_data <- list(
    log_cpm = log_cpm,
    voom_expression = voom_expression,
    filtered_counts = filtered_counts,
    metadata = metadata,
    design_matrix = design,
    normalization_factors = dge$samples$norm.factors,
    lib_sizes = dge$samples$lib.size,
    technology = "snRNA_pseudobulk"
  )
  
  cat("✓ Expression matrices prepared:\n")
  cat("  - Log2-CPM:", dim(log_cpm)[1], "genes ×", dim(log_cpm)[2], "samples\n")
  cat("  - Voom-transformed:", dim(voom_expression)[1], "genes ×", dim(voom_expression)[2], "samples\n")
  cat("  - Filtered counts:", dim(filtered_counts)[1], "genes ×", dim(filtered_counts)[2], "samples\n")
  
  return(expression_data)
}

# ==============================================================================
# ENHANCED RESULTS OUTPUT
# ==============================================================================

#' Save analysis results with organized outputs
save_enhanced_results <- function(results_list, metadata, expression_data) {
  cat("\n=== Saving Enhanced Results with Expression Data ===\n")
  
  # Create output directories
  if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
  if (!dir.exists(GSE225158_OUTPUTS)) dir.create(GSE225158_OUTPUTS, recursive = TRUE)
  
  # Combine results with expression data
  enhanced_results <- list(
    # Statistical results
    paired_limma = results_list$paired_limma,
    mixed_effects = results_list$mixed_effects,
    deseq2 = results_list$deseq2,
    
    # Expression matrices (CRITICAL for comprehensive analysis)
    expression_data = expression_data,
    
    # Metadata
    sample_metadata = metadata,
    
    # Analysis info
    analysis_date = Sys.Date(),
    dataset = "GSE225158",
    design = "paired_pseudobulk",
    note = "snRNA-seq pseudobulk with expression matrices for comprehensive analysis"
  )
  
  # Save main results object (ENHANCED)
  results_file <- file.path(OUTPUT_DIR, "GSE225158_region_analysis_results.rds")
  saveRDS(enhanced_results, results_file)
  cat("✓ Enhanced results saved:", results_file, "\n")
  
  # Save expression matrices separately for easy access
  expression_file <- file.path(OUTPUT_DIR, "GSE225158_expression_matrices.rds")
  saveRDS(expression_data, expression_file)
  cat("✓ Expression matrices saved:", expression_file, "\n")
  
  # Save individual method results as CSV (existing code)
  for (method in names(results_list)) {
    method_file <- file.path(OUTPUT_DIR, paste0("GSE225158_", method, "_results.csv"))
    write.csv(results_list[[method]]$results, method_file, row.names = TRUE)
    cat("✓", method, "results:", method_file, "\n")
  }
  
  # Create enhanced method comparison summary
  summary_df <- data.frame(
    Method = names(results_list),
    N_Genes_Tested = sapply(results_list, function(x) nrow(x$results)),
    N_Significant_FDR05 = sapply(results_list, function(x) {
      res <- x$results
      pval_col <- ifelse("adj.P.Val" %in% colnames(res), "adj.P.Val", "padj")
      sum(res[[pval_col]] < 0.05, na.rm = TRUE)
    }),
    Intra_Subject_Correlation = sapply(results_list, function(x) {
      if ("correlation" %in% names(x)) round(x$correlation, 3) else NA
    }),
    N_Subjects = nrow(metadata) / 2,
    N_Samples = nrow(metadata),
    Technology = "snRNA_pseudobulk",
    Expression_Matrices_Available = TRUE
  )
  
  summary_file <- file.path(OUTPUT_DIR, "GSE225158_method_comparison.csv")
  write.csv(summary_df, summary_file, row.names = FALSE)
  cat("✓ Enhanced summary saved:", summary_file, "\n")
  
  # Create a summary visualization in organized outputs
  summary_plot_dir <- file.path(GSE225158_OUTPUTS, "Summary_Plots")
  if (!dir.exists(summary_plot_dir)) dir.create(summary_plot_dir, recursive = TRUE)
  
  # Method comparison plot
  p_summary <- ggplot(summary_df, aes(x = Method, y = N_Significant_FDR05)) +
    geom_col(fill = "darkorange", alpha = 0.7) +
    geom_text(aes(label = N_Significant_FDR05), vjust = -0.5) +
    labs(title = "GSE225158: Significant Genes by Method",
         subtitle = "Caudate vs Putamen comparison (snRNA-seq pseudobulk)",
         x = "Analysis Method", 
         y = "Number of Significant Genes (FDR < 0.05)") +
    theme_minimal()
  
  ggsave(file.path(summary_plot_dir, "GSE225158_method_comparison.png"),
         p_summary, width = 10, height = 6, dpi = 300)
  
  cat("✓ Summary plot saved to organized outputs\n")
  cat("📊 Figures location:", GSE225158_OUTPUTS, "\n")
  cat("📁 Data location:", OUTPUT_DIR, "\n")
  
  return(summary_df)
}

# ==============================================================================
# ENHANCED MAIN ANALYSIS PIPELINE (FIXED)
# ==============================================================================

#' Run enhanced GSE225158 analysis with expression matrix preparation
run_gse225158_enhanced_analysis <- function() {
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("GSE225158 ENHANCED PAIRED PSEUDOBULK ANALYSIS\n")
  cat("Expression matrices for comprehensive neuroinflammatory analysis\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  tryCatch({
    # Load and validate data
    data <- load_pseudobulk_data()
    
    # Prepare expression matrices for comprehensive analysis
    expression_data <- prepare_expression_matrices(data$counts, data$metadata)
    
    # Run three analysis methods (same as before)
    results_list <- list(
      paired_limma = paired_limma_analysis(data$counts, data$metadata),
      mixed_effects = mixed_effects_analysis(data$counts, data$metadata),
      deseq2 = deseq2_analysis(data$counts, data$metadata)
    )
    
    # Print analysis summary
    cat("\n=== Analysis Summary ===\n")
    for (method in names(results_list)) {
      res <- results_list[[method]]$results
      pval_col <- ifelse("adj.P.Val" %in% colnames(res), "adj.P.Val", "padj")
      sig_genes <- sum(res[[pval_col]] < 0.05, na.rm = TRUE)
      cat(sprintf("%-15s: %4d significant genes (FDR < 0.05)\n", method, sig_genes))
    }
    
    # Save enhanced results with expression matrices
    summary_df <- save_enhanced_results(results_list, data$metadata, expression_data)
    
    cat("\n", paste(rep("=", 70), collapse = ""), "\n")
    cat("SUCCESS: Enhanced analysis complete!\n")
    cat("✓ Statistical results saved\n")
    cat("✓ Expression matrices saved for comprehensive analysis\n")
    cat("✓ Ready for neuroinflammatory pathway analysis\n")
    cat("Results saved to:", OUTPUT_DIR, "\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")
    
    return(list(
      results = results_list,
      expression_data = expression_data,
      metadata = data$metadata
    ))
    
  }, error = function(e) {
    cat("\nERROR:", e$message, "\n")
    stop(e)
  })
}

# ==============================================================================
# EXECUTION (UPDATED)
# ==============================================================================

# Run enhanced analysis when script is sourced
if (!exists("SOURCED")) {
  results <- run_gse225158_enhanced_analysis()
}