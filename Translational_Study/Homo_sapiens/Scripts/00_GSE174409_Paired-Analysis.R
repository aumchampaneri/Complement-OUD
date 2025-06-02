# ==============================================================================
# GSE174409 Paired Brain Region Analysis
# ==============================================================================
# Purpose: Perform paired differential expression analysis comparing DLPFC vs NAc
#          brain regions in OUD patients using bulk RNA-seq data
# Dataset: GSE174409 - 40 subjects with paired samples (DLPFC and NAc regions)
# Methods: Paired limma-voom, Mixed effects, and DESeq2 approaches
# Enhanced: Now includes expression matrices for comprehensive neuroinflammatory analysis
# ==============================================================================

# File paths
DATA_DIR <- "/Users/aumchampaneri/Complement-OUD/Translational_Study/Homo_sapiens/Data/GSE174409"
OUTPUT_DIR <- "/Users/aumchampaneri/Complement-OUD/Translational_Study/Homo_sapiens/Results/GSE174409"
COUNTS_FILE <- file.path(DATA_DIR, "GSE174409_raw_counts_02102020.csv")
METADATA_FILE <- file.path(DATA_DIR, "GSE174409_series_matrix.txt")

# NEW: Organized figure outputs
OUTPUTS_BASE <- "/Users/aumchampaneri/Complement-OUD/Translational_Study/Homo_sapiens/Outputs"
GSE174409_OUTPUTS <- file.path(OUTPUTS_BASE, "GSE174409_Analysis")

# Load required libraries
suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(tibble)  # Add tibble for proper rownames_to_column
})

# ==============================================================================
# DATA LOADING AND PREPROCESSING
# ==============================================================================

#' Load and prepare GSE174409 bulk RNA-seq data
#' @param counts_file Path to counts CSV file
#' @param metadata_file Path to series matrix file
#' @return List with counts matrix and metadata
load_and_prepare_bulk <- function(counts_file, metadata_file) {
  cat("=== Loading GSE174409 Data ===\n")
  
  # Load count matrix
  counts <- read.csv(counts_file, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  cat("Count matrix dimensions:", dim(counts), "\n")
  
  # CRITICAL: Convert Ensembl IDs to gene symbols for consistency with GSE225158
  cat("Converting Ensembl IDs to gene symbols...\n")
  ensembl_ids <- rownames(counts)
  cat("Sample Ensembl IDs:", paste(head(ensembl_ids), collapse = ", "), "\n")
  
  # Map Ensembl IDs to gene symbols
  tryCatch({
    gene_mapping <- clusterProfiler::bitr(ensembl_ids, 
                                          fromType = "ENSEMBL", 
                                          toType = "SYMBOL", 
                                          OrgDb = org.Hs.eg.db)
    
    cat("Mapped", nrow(gene_mapping), "out of", length(ensembl_ids), "Ensembl IDs\n")
    
    # Keep only genes with successful mapping
    counts_mapped <- counts[gene_mapping$ENSEMBL, ]
    
    # Handle duplicate gene symbols BEFORE any operations
    if (any(duplicated(gene_mapping$SYMBOL))) {
      cat("Handling", sum(duplicated(gene_mapping$SYMBOL)), "duplicate gene symbols...\n")
      
      # Simple approach: aggregate using base R
      # Create mapping between Ensembl and symbols
      ensembl_to_symbol <- setNames(gene_mapping$SYMBOL, gene_mapping$ENSEMBL)
      
      # Add gene symbols to the matrix
      gene_symbols <- ensembl_to_symbol[rownames(counts_mapped)]
      
      # Create a data frame for aggregation
      temp_df <- data.frame(
        gene_symbol = gene_symbols,
        counts_mapped,
        stringsAsFactors = FALSE
      )
      
      # Aggregate by gene symbol using base R
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
    
    cat("Final count matrix after gene symbol conversion:", dim(counts), "\n")
    
  }, error = function(e) {
    cat("Warning: Gene symbol conversion failed:", e$message, "\n")
    cat("Proceeding with original Ensembl IDs\n")
  })
  
  # Parse metadata from series matrix file
  lines <- readLines(metadata_file)
  
  # Extract sample information
  sample_id_line <- grep("!Sample_geo_accession", lines, value = TRUE)
  sample_ids <- unlist(strsplit(gsub("!Sample_geo_accession\\s+", "", sample_id_line), "\\t"))
  sample_ids <- sample_ids[sample_ids != ""]
  
  title_line <- grep("!Sample_title", lines, value = TRUE)
  sample_titles <- unlist(strsplit(gsub("!Sample_title\\s+", "", title_line), "\\t"))
  sample_titles <- sample_titles[sample_titles != ""]
  
  # Create initial metadata
  metadata <- data.frame(
    Sample_ID = sample_ids,
    Title = sample_titles,
    stringsAsFactors = FALSE
  )
  
  # Parse characteristics
  char_lines <- grep("!Sample_characteristics_ch1", lines, value = TRUE)
  for(i in seq_along(char_lines)) {
    char_values <- unlist(strsplit(gsub("!Sample_characteristics_ch1\\s+", "", char_lines[i]), "\\t"))
    char_values <- char_values[char_values != ""]
    
    if(length(char_values) > 0 && grepl(":", char_values[1])) {
      var_name <- gsub("([^:]+):.*", "\\1", char_values[1])
      var_values <- gsub(".*?:\\s*", "", char_values)
      
      # Clean variable name
      var_name <- gsub("[^A-Za-z0-9_]", "_", var_name)
      var_name <- gsub("_+", "_", var_name)
      var_name <- gsub("^_|_$", "", var_name)
      
      if(length(var_values) == length(sample_ids)) {
        metadata[[var_name]] <- var_values
      }
    }
  }
  
  # Match samples between counts and metadata using cleaned titles
  sample_titles_clean <- gsub('^"|"$', '', sample_titles)
  metadata$Sample_ID <- sample_titles_clean
  counts <- counts[, sample_titles_clean]
  
  # Clean metadata - remove quotes
  for(col in colnames(metadata)) {
    if(is.character(metadata[[col]])) {
      metadata[[col]] <- gsub('^"|"$', '', metadata[[col]])
    }
  }
  
  # Create standardized variables
  metadata$Subject_ID <- gsub("^(HN|HBA)", "", metadata$Title)  # Extract subject number
  metadata$Brain_Region <- ifelse(grepl("^HN", metadata$Title), "DLPFC", "NAc")  # Infer brain region
  metadata$OUD_Status <- metadata$diagnosis
  metadata$Sex <- factor(metadata$Sex)
  metadata$Region <- factor(metadata$Brain_Region)
  metadata$OUD_Status <- factor(metadata$OUD_Status)
  
  # Filter for paired subjects only
  paired_subjects <- metadata %>%
    group_by(Subject_ID) %>%
    filter(n_distinct(Brain_Region) == 2) %>%
    pull(Subject_ID) %>%
    unique()
  
  metadata <- metadata %>%
    filter(Subject_ID %in% paired_subjects) %>%
    arrange(Subject_ID, Brain_Region)
  
  counts <- counts[, metadata$Sample_ID]
  
  cat("Paired subjects found:", length(paired_subjects), "\n")
  cat("Final dimensions - Genes:", nrow(counts), "| Samples:", ncol(counts), "\n")
  cat("Brain regions:", table(metadata$Region), "\n")
  cat("OUD status:", table(metadata$OUD_Status), "\n\n")
  
  return(list(counts = counts, metadata = metadata))
}

# ==============================================================================
# ENHANCED EXPRESSION MATRIX PREPARATION
# ==============================================================================

#' Prepare expression matrices for comprehensive analysis
prepare_expression_matrices <- function(counts, metadata) {
  cat("\n=== Preparing Expression Matrices for Comprehensive Analysis ===\n")
  
  # Create DGEList for normalization
  dge <- DGEList(counts = counts)
  keep <- filterByExpr(dge, group = metadata$Region)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  
  cat("Genes after filtering:", nrow(dge), "\n")
  
  # Log2-CPM transformation
  log_cpm <- cpm(dge, log = TRUE, prior.count = 1)
  
  # Voom-transformed values
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
    technology = "bulk_RNA_seq"
  )
  
  cat("✓ Expression matrices prepared:\n")
  cat("  - Log2-CPM:", dim(log_cpm)[1], "genes ×", dim(log_cpm)[2], "samples\n")
  cat("  - Voom-transformed:", dim(voom_expression)[1], "genes ×", dim(voom_expression)[2], "samples\n")
  cat("  - Filtered counts:", dim(filtered_counts)[1], "genes ×", dim(filtered_counts)[2], "samples\n")
  
  return(expression_data)
}

# ==============================================================================
# ANALYSIS METHODS
# ==============================================================================

#' Method 1: Paired limma-voom analysis
#' @param counts Count matrix
#' @param metadata Sample metadata
#' @return Analysis results
paired_limma_analysis <- function(counts, metadata) {
  cat("=== Paired Limma-Voom Analysis ===\n")
  
  # Create DGEList and filter
  dge <- DGEList(counts = counts)
  keep <- filterByExpr(dge, group = metadata$Region)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  cat("Genes after filtering:", nrow(dge), "\n")
  
  # Design matrix
  design <- model.matrix(~ Region + Sex + OUD_Status, data = metadata)
  
  # Voom transformation with paired correlation
  v <- voom(dge, design, plot = FALSE)
  corfit <- duplicateCorrelation(v, design, block = metadata$Subject_ID)
  cat("Intra-subject correlation:", round(corfit$consensus, 3), "\n")
  
  # Fit model
  fit <- lmFit(v, design, block = metadata$Subject_ID, correlation = corfit$consensus)
  fit <- eBayes(fit)
  
  # Extract results
  results <- topTable(fit, coef = "RegionNAc", number = Inf, sort.by = "P")
  results$Method <- "Paired_Limma"
  
  return(list(results = results, fit = fit, correlation = corfit$consensus))
}

#' Method 2: Mixed effects limma analysis
#' @param counts Count matrix
#' @param metadata Sample metadata
#' @return Analysis results
mixed_effects_analysis <- function(counts, metadata) {
  cat("=== Mixed Effects Analysis ===\n")
  
  # Create DGEList and filter
  dge <- DGEList(counts = counts)
  keep <- filterByExpr(dge, group = metadata$Region)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  cat("Genes after filtering:", nrow(dge), "\n")
  
  # Design matrix
  design <- model.matrix(~ Region + Sex + OUD_Status, data = metadata)
  
  # Voom with quality weights and correlation
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
  results <- topTable(fit, coef = "RegionNAc", number = Inf, sort.by = "P")
  results$Method <- "Mixed_Effects"
  
  return(list(results = results, fit = fit, correlation = corfit$consensus))
}

#' Method 3: DESeq2 analysis
#' @param counts Count matrix
#' @param metadata Sample metadata
#' @return Analysis results
deseq2_analysis <- function(counts, metadata) {
  cat("=== DESeq2 Analysis ===\n")
  
  # Create DESeq2 object
  counts_int <- round(as.matrix(counts))
  design_formula <- ~ Region + Sex + OUD_Status
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts_int,
    colData = metadata,
    design = design_formula
  )
  
  # Filter and run DESeq2
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  cat("Genes after filtering:", nrow(dds), "\n")
  
  dds <- DESeq(dds, quiet = TRUE)
  
  # Extract results
  results_obj <- results(dds, contrast = c("Region", "NAc", "DLPFC"))
  results_df <- as.data.frame(results_obj)
  results_df$Method <- "DESeq2"
  
  return(list(results = results_df, dds = dds))
}

# ==============================================================================
# ENHANCED RESULTS SAVING
# ==============================================================================

#' Save enhanced analysis results with organized outputs
save_enhanced_results <- function(results_list, metadata, expression_data) {
  cat("\n=== Saving Enhanced Results with Expression Data ===\n")
  
  # Create output directories
  if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
  if (!dir.exists(GSE174409_OUTPUTS)) dir.create(GSE174409_OUTPUTS, recursive = TRUE)
  
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
    dataset = "GSE174409",
    design = "paired_region",
    note = "Bulk RNA-seq with expression matrices for comprehensive analysis"
  )
  
  # Save enhanced results object
  results_file <- file.path(OUTPUT_DIR, "GSE174409_region_analysis_results.rds")
  saveRDS(enhanced_results, results_file)
  cat("✓ Enhanced results saved:", results_file, "\n")
  
  # Save expression matrices separately
  expression_file <- file.path(OUTPUT_DIR, "GSE174409_expression_matrices.rds")
  saveRDS(expression_data, expression_file)
  cat("✓ Expression matrices saved:", expression_file, "\n")
  
  # Save individual method results as CSV
  for (method in names(results_list)) {
    method_file <- file.path(OUTPUT_DIR, paste0("GSE174409_", method, "_results.csv"))
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
    N_Subjects = length(unique(metadata$Subject_ID)),
    N_Samples = nrow(metadata),
    Technology = "bulk_RNA_seq",
    Expression_Matrices_Available = TRUE
  )
  
  summary_file <- file.path(OUTPUT_DIR, "GSE174409_method_comparison.csv")
  write.csv(summary_df, summary_file, row.names = FALSE)
  cat("✓ Enhanced summary saved:", summary_file, "\n")
  
  # Create a summary visualization in organized outputs
  summary_plot_dir <- file.path(GSE174409_OUTPUTS, "Summary_Plots")
  if (!dir.exists(summary_plot_dir)) dir.create(summary_plot_dir, recursive = TRUE)
  
  # Method comparison plot
  p_summary <- ggplot(summary_df, aes(x = Method, y = N_Significant_FDR05)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_text(aes(label = N_Significant_FDR05), vjust = -0.5) +
    labs(title = "GSE174409: Significant Genes by Method",
         subtitle = "DLPFC vs NAc comparison",
         x = "Analysis Method",
         y = "Number of Significant Genes (FDR < 0.05)") +
    theme_minimal()
  
  ggsave(file.path(summary_plot_dir, "GSE174409_method_comparison.png"),
         p_summary, width = 10, height = 6, dpi = 300)
  
  cat("✓ Summary plot saved to organized outputs\n")
  cat("📊 Figures location:", GSE174409_OUTPUTS, "\n")
  cat("📁 Data location:", OUTPUT_DIR, "\n")
  
  return(summary_df)
}

# ==============================================================================
# MAIN ANALYSIS FUNCTION (ENHANCED)
# ==============================================================================

#' Run complete GSE174409 paired analysis with expression matrices
#' @return List of results from all methods plus expression data
run_gse174409_analysis <- function() {
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("GSE174409 ENHANCED PAIRED ANALYSIS\n")
  cat("Bulk RNA-seq with expression matrices for comprehensive analysis\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  tryCatch({
    # Load and prepare data
    bulk_data <- load_and_prepare_bulk(COUNTS_FILE, METADATA_FILE)
    
    # Prepare expression matrices for comprehensive analysis
    expression_data <- prepare_expression_matrices(bulk_data$counts, bulk_data$metadata)
    
    # Run analyses
    results_list <- list(
      paired_limma = paired_limma_analysis(bulk_data$counts, bulk_data$metadata),
      mixed_effects = mixed_effects_analysis(bulk_data$counts, bulk_data$metadata),
      deseq2 = deseq2_analysis(bulk_data$counts, bulk_data$metadata)
    )
    
    # Summary statistics
    cat("\n=== Analysis Summary ===\n")
    for(method in names(results_list)) {
      res <- results_list[[method]]$results
      pval_col <- ifelse("adj.P.Val" %in% colnames(res), "adj.P.Val", "padj")
      sig_genes <- sum(res[[pval_col]] < 0.05, na.rm = TRUE)
      cat(sprintf("%-15s: %4d significant genes (FDR < 0.05)\n", method, sig_genes))
    }
    
    # Save enhanced results
    summary_df <- save_enhanced_results(results_list, bulk_data$metadata, expression_data)
    
    cat("\n", paste(rep("=", 70), collapse = ""), "\n")
    cat("SUCCESS: Enhanced GSE174409 analysis complete!\n")
    cat("✓ Statistical results saved\n")
    cat("✓ Expression matrices saved for comprehensive analysis\n")
    cat("✓ Ready for neuroinflammatory pathway analysis\n")
    cat("Results saved to:", OUTPUT_DIR, "\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")
    
    return(list(
      results = results_list,
      expression_data = expression_data,
      metadata = bulk_data$metadata
    ))
    
  }, error = function(e) {
    cat("\nERROR:", e$message, "\n")
    stop(e)
  })
}

# ==============================================================================
# EXECUTION
# ==============================================================================

# Run analysis if script is sourced directly
if(!exists("SOURCED")) {
  results <- run_gse174409_analysis()
}