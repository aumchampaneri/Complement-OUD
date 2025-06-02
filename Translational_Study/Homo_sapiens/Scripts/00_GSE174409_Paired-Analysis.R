# ==============================================================================
# GSE174409 Paired Brain Region Analysis
# ==============================================================================
# Purpose: Perform paired differential expression analysis comparing DLPFC vs NAc
#          brain regions in OUD patients using bulk RNA-seq data
# Dataset: GSE174409 - 40 subjects with paired samples (DLPFC and NAc regions)
# Methods: Paired limma-voom, Mixed effects, and DESeq2 approaches
# ==============================================================================

# File paths
DATA_DIR <- "/Users/aumchampaneri/Complement-OUD/Translational_Study/Homo_sapiens/Data/GSE174409"
OUTPUT_DIR <- "/Users/aumchampaneri/Complement-OUD/Translational_Study/Homo_sapiens/Results/GSE174409"
COUNTS_FILE <- file.path(DATA_DIR, "GSE174409_raw_counts_02102020.csv")
METADATA_FILE <- file.path(DATA_DIR, "GSE174409_series_matrix.txt")

# Load required libraries
library(limma)
library(edgeR)
library(DESeq2)
library(dplyr)

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
# MAIN ANALYSIS FUNCTION
# ==============================================================================

#' Run complete GSE174409 paired analysis
#' @return List of results from all methods
run_gse174409_analysis <- function() {
  # Load and prepare data
  bulk_data <- load_and_prepare_bulk(COUNTS_FILE, METADATA_FILE)
  
  # Run analyses
  results_list <- list(
    paired_limma = paired_limma_analysis(bulk_data$counts, bulk_data$metadata),
    mixed_effects = mixed_effects_analysis(bulk_data$counts, bulk_data$metadata),
    deseq2 = deseq2_analysis(bulk_data$counts, bulk_data$metadata)
  )
  
  # Summary statistics
  cat("=== Analysis Summary ===\n")
  for(method in names(results_list)) {
    res <- results_list[[method]]$results
    pval_col <- ifelse("adj.P.Val" %in% colnames(res), "adj.P.Val", "padj")
    sig_genes <- sum(res[[pval_col]] < 0.05, na.rm = TRUE)
    cat(sprintf("%-15s - Significant genes (FDR < 0.05): %d\n", method, sig_genes))
  }
  
  # Save results
  if(!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
  
  saveRDS(results_list, file.path(OUTPUT_DIR, "GSE174409_region_analysis_results.rds"))
  
  for(method in names(results_list)) {
    write.csv(results_list[[method]]$results, 
              file.path(OUTPUT_DIR, paste0("GSE174409_", method, "_results.csv")), 
              row.names = TRUE)
  }
  
  # Create comparison summary
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
  
  write.csv(summary_df, file.path(OUTPUT_DIR, "GSE174409_method_comparison.csv"), row.names = FALSE)
  
  cat("\n=== Results Saved ===\n")
  cat("Main results:", file.path(OUTPUT_DIR, "GSE174409_region_analysis_results.rds"), "\n")
  
  # Print structure of saved results
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
  results <- run_gse174409_analysis()
}