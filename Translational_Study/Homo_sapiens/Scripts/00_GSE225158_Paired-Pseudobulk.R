# ==============================================================================
# GSE225158 Paired Pseudobulk Analysis
# ==============================================================================
# Purpose: Perform paired differential expression analysis on snRNA-seq data
#          by creating pseudobulk samples and comparing brain regions
# Dataset: GSE225158 - Single-nucleus RNA-seq data aggregated to pseudobulk
# Methods: Paired limma-voom with duplicateCorrelation
# ==============================================================================

# File paths
DATA_DIR <- "/Users/aumchampaneri/Complement-OUD/Translational_Study/Homo_sapiens/Data/GSE225158"
OUTPUT_DIR <- "/Users/aumchampaneri/Complement-OUD/Translational_Study/Homo_sapiens/Results"

# Load required libraries
library(Seurat)
library(limma)
library(edgeR)
library(dplyr)
library(tidyr)

# ==============================================================================
# DATA LOADING AND PREPROCESSING
# ==============================================================================

# Load and prepare snRNA-seq data - aggregate across all cell types
load_and_prepare_snrnaseq <- function(seurat_object) {
  # Create pseudobulk by aggregating counts per subject, region, and disease status
  # This sums across all cell types to create bulk-like data
  pseudobulk_counts <- AggregateExpression(
    seurat_object,
    group.by = c("Subject_ID", "Region", "Sex", "OUD_Status"),
    assays = "RNA",
    slot = "counts",
    return.seurat = FALSE
  )
  
  return(pseudobulk_counts$RNA)  # Return just the count matrix
}

# Model region effect with paired design
model_region_effect_pseudobulk <- function(counts_matrix) {
  # Create sample metadata from column names
  sample_info <- data.frame(
    Sample_ID = colnames(counts_matrix),
    stringsAsFactors = FALSE
  )
  
  # Parse sample info (assuming format: SubjectID_Region_Sex_OUDStatus)
  sample_info <- sample_info %>%
    separate(Sample_ID, 
             into = c("Subject_ID", "Region", "Sex", "OUD_Status"), 
             sep = "_", 
             remove = FALSE)
  
  # Filter for paired samples (subjects with both brain regions)
  paired_subjects <- sample_info %>%
    group_by(Subject_ID) %>%
    filter(n_distinct(Region) == 2) %>%
    pull(Subject_ID) %>%
    unique()
  
  cat("Found", length(paired_subjects), "subjects with paired samples\n")
  
  if(length(paired_subjects) < 3) {
    stop("Insufficient paired subjects. Need at least 3 for analysis.")
  }
  
  sample_info_paired <- sample_info %>%
    filter(Subject_ID %in% paired_subjects) %>%
    arrange(Subject_ID, Region)
  
  counts_paired <- counts_matrix[, sample_info_paired$Sample_ID]
  
  # Create DGEList object
  dge <- DGEList(counts = counts_paired)
  
  # Filter low expressed genes
  keep <- filterByExpr(dge, group = sample_info_paired$Region)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  cat("After filtering:", nrow(dge), "genes retained\n")
  
  # Normalize
  dge <- calcNormFactors(dge)
  
  # Create design matrix for paired analysis
  # Make sure factor levels are set correctly
  sample_info_paired$Region <- factor(sample_info_paired$Region)
  sample_info_paired$Sex <- factor(sample_info_paired$Sex)
  sample_info_paired$OUD_Status <- factor(sample_info_paired$OUD_Status)
  
  # Design matrix including region, sex, and disease status
  design <- model.matrix(~ Region + Sex + OUD_Status, data = sample_info_paired)
  
  # Account for subject pairing using duplicateCorrelation
  v <- voom(dge, design, plot = TRUE)
  
  # Calculate correlation for paired subjects
  corfit <- duplicateCorrelation(v, design, block = sample_info_paired$Subject_ID)
  
  cat("Intra-subject correlation:", round(corfit$consensus, 3), "\n")
  
  # Fit model with correlation
  fit <- lmFit(v, design, block = sample_info_paired$Subject_ID, 
               correlation = corfit$consensus)
  fit <- eBayes(fit)
  
  # Extract results for each coefficient
  results_list <- list()
  
  # Region effect
  region_coef <- grep("Region", colnames(design), value = TRUE)
  if(length(region_coef) > 0) {
    results_list$Region <- topTable(fit, coef = region_coef[1], 
                                   number = Inf, sort.by = "P")
    results_list$Region$Comparison <- "Region_Effect"
  }
  
  # Sex effect
  sex_coef <- grep("Sex", colnames(design), value = TRUE)
  if(length(sex_coef) > 0) {
    results_list$Sex <- topTable(fit, coef = sex_coef[1], 
                                number = Inf, sort.by = "P")
    results_list$Sex$Comparison <- "Sex_Effect"
  }
  
  # OUD Status effect
  oud_coef <- grep("OUD_Status", colnames(design), value = TRUE)
  if(length(oud_coef) > 0) {
    results_list$OUD_Status <- topTable(fit, coef = oud_coef[1], 
                                       number = Inf, sort.by = "P")
    results_list$OUD_Status$Comparison <- "OUD_Effect"
  }
  
  return(list(
    results = results_list,
    fit = fit,
    design = design,
    correlation = corfit$consensus,
    sample_info = sample_info_paired
  ))
}

# Load the Seurat object
load_seurat_data <- function() {
  h5_files <- list.files(DATA_DIR, pattern = "\\.h5seurat$", full.names = TRUE)
  
  if(length(h5_files) == 0) {
    stop("No .h5seurat files found in the data directory")
  }
  
  cat("Loading Seurat object from:", h5_files[1], "\n")
  seurat_obj <- LoadH5Seurat(h5_files[1])
  
  # Check metadata structure
  cat("Metadata columns:", colnames(seurat_obj@meta.data), "\n")
  cat("Number of cells:", ncol(seurat_obj), "\n")
  cat("Number of genes:", nrow(seurat_obj), "\n")
  
  # Check for required columns
  required_cols <- c("Region", "Subject_ID", "Sex", "OUD_Status")
  missing_cols <- setdiff(required_cols, colnames(seurat_obj@meta.data))
  
  if(length(missing_cols) > 0) {
    cat("Missing required columns:", paste(missing_cols, collapse = ", "), "\n")
    cat("Please ensure metadata contains: Region, Subject_ID, Sex, OUD_Status\n")
  }
  
  return(seurat_obj)
}

# Main analysis function
run_pseudobulk_analysis <- function() {
  # Load data
  cat("=== Loading GSE225158 Data ===\n")
  seurat_obj <- load_seurat_data()
  
  # Create pseudobulk data
  cat("\n=== Creating Pseudobulk Data ===\n")
  pseudobulk_counts <- load_and_prepare_snrnaseq(seurat_obj)
  
  cat("Pseudobulk matrix dimensions:", dim(pseudobulk_counts), "\n")
  cat("Sample names:", head(colnames(pseudobulk_counts)), "\n")
  
  # Run paired analysis
  cat("\n=== Running Paired Analysis ===\n")
  results <- model_region_effect_pseudobulk(pseudobulk_counts)
  
  # Print summary
  cat("\n=== Analysis Summary ===\n")
  for(comparison in names(results$results)) {
    sig_genes <- sum(results$results[[comparison]]$adj.P.Val < 0.05, na.rm = TRUE)
    cat(comparison, "- Significant genes (FDR < 0.05):", sig_genes, "\n")
  }
  
  # Save results
  if(!dir.exists(OUTPUT_DIR)) {
    dir.create(OUTPUT_DIR, recursive = TRUE)
  }
  
  saveRDS(results, file.path(OUTPUT_DIR, "GSE225158_pseudobulk_paired_results.rds"))
  
  # Save individual result tables
  for(comparison in names(results$results)) {
    write.csv(results$results[[comparison]], 
              file.path(OUTPUT_DIR, paste0("GSE225158_", comparison, "_results.csv")), 
              row.names = TRUE)
  }
  
  # Create summary table
  summary_df <- data.frame(
    Comparison = names(results$results),
    N_Genes_Tested = sapply(results$results, nrow),
    N_Significant_FDR05 = sapply(results$results, function(x) sum(x$adj.P.Val < 0.05, na.rm = TRUE)),
    N_Subjects = nrow(results$sample_info) / 2,  # Paired design
    N_Samples = nrow(results$sample_info),
    Correlation = round(results$correlation, 3)
  )
  
  write.csv(summary_df, file.path(OUTPUT_DIR, "GSE225158_pseudobulk_summary.csv"), row.names = FALSE)
  
  cat("\n=== Results Saved ===\n")
  cat("Main results:", file.path(OUTPUT_DIR, "GSE225158_pseudobulk_paired_results.rds"), "\n")
  cat("Summary table:", file.path(OUTPUT_DIR, "GSE225158_pseudobulk_summary.csv"), "\n")
  
  # Print structure of saved results
  cat("\n=== Results Structure ===\n")
  cat("The saved .rds file contains:\n")
  str(results, max.level = 2, give.attr = FALSE)
  
  return(results)
}

# ==============================================================================
# EXECUTION
# ==============================================================================

# Run the analysis
if(!exists("SOURCED")) {
  results <- run_pseudobulk_analysis()
}