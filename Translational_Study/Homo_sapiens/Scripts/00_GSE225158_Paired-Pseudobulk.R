# ==============================================================================
# GSE225158 Paired Pseudobulk Analysis
# ==============================================================================
# Purpose: Perform paired differential expression analysis on pseudobulk data
#          created from snRNA-seq data (mirrors GSE174409 workflow exactly)
# Dataset: GSE225158 - Pseudobulk data from single-nucleus RNA-seq
# Methods: Paired limma-voom, Mixed effects, and DESeq2 approaches (same as GSE174409)
# ==============================================================================

# File paths
DATA_DIR <- "/Users/aumchampaneri/Complement-OUD/Translational_Study/Homo_sapiens/Data/GSE225158"
OUTPUT_DIR <- "/Users/aumchampaneri/Complement-OUD/Translational_Study/Homo_sapiens/Results/GSE225158"
COUNTS_FILE <- file.path(OUTPUT_DIR, "GSE225158_pseudobulk_counts.csv")
METADATA_FILE <- file.path(OUTPUT_DIR, "GSE225158_pseudobulk_metadata.csv")

# Load required libraries
library(limma)
library(edgeR)
library(DESeq2)
library(dplyr)

# ==============================================================================
# DATA LOADING AND PREPROCESSING
# ==============================================================================

#' Load and prepare GSE225158 pseudobulk data (mirrors GSE174409 exactly)
#' @param counts_file Path to pseudobulk counts CSV file
#' @param metadata_file Path to pseudobulk metadata CSV file
#' @return List with counts matrix and metadata
load_and_prepare_pseudobulk <- function(counts_file, metadata_file) {
  cat("=== Loading GSE225158 Pseudobulk Data ===\n")
  
  # Load count matrix (same as GSE174409)
  counts <- read.csv(counts_file, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  cat("Count matrix dimensions:", dim(counts), "\n")
  
  # Load metadata (same as GSE174409)
  metadata <- read.csv(metadata_file, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  cat("Metadata dimensions:", dim(metadata), "\n")
  
  # Ensure sample order matches
  common_samples <- intersect(colnames(counts), rownames(metadata))
  counts <- counts[, common_samples]
  metadata <- metadata[common_samples, ]
  
  # Convert to factors (same as GSE174409)
  metadata$Region <- factor(metadata$Region)
  metadata$Sex <- factor(metadata$Sex) 
  metadata$OUD_Status <- factor(metadata$OUD_Status)
  
  # Filter for paired subjects only (same logic as GSE174409)
  paired_subjects <- metadata %>%
    group_by(Subject_ID) %>%
    filter(n_distinct(Region) == 2) %>%
    pull(Subject_ID) %>%
    unique()
  
  metadata <- metadata %>%
    filter(Subject_ID %in% paired_subjects) %>%
    arrange(Subject_ID, Region)
  
  counts <- counts[, metadata$Sample_ID]
  
  cat("Paired subjects found:", length(paired_subjects), "\n")
  cat("Final dimensions - Genes:", nrow(counts), "| Samples:", ncol(counts), "\n")
  cat("Brain regions:", table(metadata$Region), "\n")
  cat("OUD status:", table(metadata$OUD_Status), "\n\n")
  
  return(list(counts = counts, metadata = metadata))
}

# ==============================================================================
# ANALYSIS METHODS (IDENTICAL TO GSE174409)
# ==============================================================================

#' Method 1: Paired limma-voom analysis (identical to GSE174409)
paired_limma_analysis <- function(counts, metadata) {
  cat("=== Paired Limma-Voom Analysis ===\n")
  
  dge <- DGEList(counts = counts)
  keep <- filterByExpr(dge, group = metadata$Region)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  cat("Genes after filtering:", nrow(dge), "\n")
  
  design <- model.matrix(~ Region + Sex + OUD_Status, data = metadata)
  
  v <- voom(dge, design, plot = FALSE)
  corfit <- duplicateCorrelation(v, design, block = metadata$Subject_ID)
  cat("Intra-subject correlation:", round(corfit$consensus, 3), "\n")
  
  fit <- lmFit(v, design, block = metadata$Subject_ID, correlation = corfit$consensus)
  fit <- eBayes(fit)
  
  region_coef <- grep("Region", colnames(design), value = TRUE)[1]
  results <- topTable(fit, coef = region_coef, number = Inf, sort.by = "P")
  results$Method <- "Paired_Limma"
  
  return(list(results = results, fit = fit, correlation = corfit$consensus))
}

#' Method 2: Mixed effects analysis (identical to GSE174409)
mixed_effects_analysis <- function(counts, metadata) {
  cat("=== Mixed Effects Analysis ===\n")
  
  dge <- DGEList(counts = counts)
  keep <- filterByExpr(dge, group = metadata$Region)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  cat("Genes after filtering:", nrow(dge), "\n")
  
  design <- model.matrix(~ Region + Sex + OUD_Status, data = metadata)
  
  v <- voomWithQualityWeights(dge, design, plot = FALSE)
  corfit <- duplicateCorrelation(v, design, block = metadata$Subject_ID)
  
  v <- voomWithQualityWeights(dge, design, plot = FALSE, 
                              block = metadata$Subject_ID, 
                              correlation = corfit$consensus)
  corfit <- duplicateCorrelation(v, design, block = metadata$Subject_ID)
  cat("Intra-subject correlation:", round(corfit$consensus, 3), "\n")
  
  fit <- lmFit(v, design, block = metadata$Subject_ID, correlation = corfit$consensus)
  fit <- eBayes(fit)
  
  region_coef <- grep("Region", colnames(design), value = TRUE)[1]
  results <- topTable(fit, coef = region_coef, number = Inf, sort.by = "P")
  results$Method <- "Mixed_Effects"
  
  return(list(results = results, fit = fit, correlation = corfit$consensus))
}

#' Method 3: DESeq2 analysis (identical to GSE174409)
deseq2_analysis <- function(counts, metadata) {
  cat("=== DESeq2 Analysis ===\n")
  
  counts_int <- round(as.matrix(counts))
  design_formula <- ~ Region + Sex + OUD_Status
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts_int,
    colData = metadata,
    design = design_formula
  )
  
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  cat("Genes after filtering:", nrow(dds), "\n")
  
  dds <- DESeq(dds, quiet = TRUE)
  
  region_levels <- levels(metadata$Region)
  results_obj <- results(dds, contrast = c("Region", region_levels[2], region_levels[1]))
  results_df <- as.data.frame(results_obj)
  results_df$Method <- "DESeq2"
  
  return(list(results = results_df, dds = dds))
}

# ==============================================================================
# MAIN ANALYSIS FUNCTION (IDENTICAL TO GSE174409)
# ==============================================================================

#' Run complete GSE225158 pseudobulk paired analysis
run_gse225158_analysis <- function() {
  # Load and prepare data
  pseudobulk_data <- load_and_prepare_pseudobulk(COUNTS_FILE, METADATA_FILE)
  
  # Run analyses (same three methods as GSE174409)
  results_list <- list(
    paired_limma = paired_limma_analysis(pseudobulk_data$counts, pseudobulk_data$metadata),
    mixed_effects = mixed_effects_analysis(pseudobulk_data$counts, pseudobulk_data$metadata),
    deseq2 = deseq2_analysis(pseudobulk_data$counts, pseudobulk_data$metadata)
  )
  
  # Summary statistics (identical to GSE174409)
  cat("=== Analysis Summary ===\n")
  for(method in names(results_list)) {
    res <- results_list[[method]]$results
    pval_col <- ifelse("adj.P.Val" %in% colnames(res), "adj.P.Val", "padj")
    sig_genes <- sum(res[[pval_col]] < 0.05, na.rm = TRUE)
    cat(sprintf("%-15s - Significant genes (FDR < 0.05): %d\n", method, sig_genes))
  }
  
  # Save results (identical structure to GSE174409)
  if(!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
  
  saveRDS(results_list, file.path(OUTPUT_DIR, "GSE225158_region_analysis_results.rds"))
  
  for(method in names(results_list)) {
    write.csv(results_list[[method]]$results, 
              file.path(OUTPUT_DIR, paste0("GSE225158_", method, "_results.csv")), 
              row.names = TRUE)
  }
  
  # Create comparison summary (identical to GSE174409)
  summary_df <- data.frame(
    Method = names(results_list),
    N_Genes_Tested = sapply(results_list, function(x) nrow(x$results)),
    N_Significant_FDR05 = sapply(results_list, function(x) {
      res <- x$results
      pval_col <- ifelse("adj.P.Val" %in% colnames(res), "adj.P.Val", "padj")
      sum(res[[pval_col]] < 0.05, na.rm = TRUE)
    }),
    Correlation = sapply(results_list, function(x) {
      if("correlation" %in% names(x)) round(x$correlation, 3) else NA
    })
  )
  
  write.csv(summary_df, file.path(OUTPUT_DIR, "GSE225158_method_comparison.csv"), row.names = FALSE)
  
  cat("\n=== Results Saved ===\n")
  cat("Main results:", file.path(OUTPUT_DIR, "GSE225158_region_analysis_results.rds"), "\n")
  
  # Print structure of saved results (identical to GSE174409)
  cat("\n=== Results Structure ===\n")
  cat("The saved .rds file contains:\n")
  str(results_list, max.level = 2, give.attr = FALSE)
  
  return(results_list)
}

# ==============================================================================
# EXECUTION
# ==============================================================================

# Run analysis if script is sourced directly
if(!exists("SOURCED")) {
  results <- run_gse225158_analysis()
}