# ==============================================================================
# Pseudobulk Differential Expression Analysis using DESeq2
# Aggregates cells by sample and condition for proper statistical modeling
# Accounts for biological replicates and reduces single-cell noise
# ==============================================================================

# Load required libraries
library(Seurat)
library(DESeq2)
library(dplyr)
library(tidyr)
library(openxlsx)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# Set working directory and paths
base_dir <- "/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus/GSE207128"
data_dir <- file.path(base_dir, "Outputs/04_Annotated_Data")
output_dir <- file.path(base_dir, "Outputs/05_Analysis_Results/Pseudobulk_Analysis")

# Create output directories
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "Tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "Plots"), recursive = TRUE, showWarnings = FALSE)

cat("Starting pseudobulk differential expression analysis...\n")
cat("Date:", Sys.Date(), "\n")
cat("Time:", Sys.time(), "\n\n")

# ==============================================================================
# Load Data
# ==============================================================================

cat("Loading annotated Seurat object...\n")
annotated_files <- list.files(data_dir, pattern = "*.rds", full.names = TRUE)

# Load filtered annotated file
filtered_file <- grep("filtered_annotated", annotated_files, value = TRUE)
if(length(filtered_file) > 0) {
  seurat_obj <- readRDS(filtered_file[1])
  cat("Loaded:", basename(filtered_file[1]), "\n")
} else if(length(annotated_files) > 0) {
  seurat_obj <- readRDS(annotated_files[1])
  cat("Loaded:", basename(annotated_files[1]), "\n")
} else {
  stop("No RDS files found in the annotated data directory")
}

# Use predicted cell types
seurat_obj$cell_type_final <- seurat_obj$predicted_celltype

# Data summary
cat("\nData summary:\n")
cat("- Total cells:", ncol(seurat_obj), "\n")
cat("- Total genes:", nrow(seurat_obj), "\n")
cat("- Cell types:", length(unique(seurat_obj$cell_type_final)), "\n")
cat("- Conditions:", paste(unique(seurat_obj$condition), collapse = ", "), "\n")
cat("- Samples:", paste(unique(seurat_obj$sample_id), collapse = ", "), "\n")

# Check sample-condition mapping
sample_info <- seurat_obj@meta.data %>%
  select(sample_id, condition) %>%
  distinct() %>%
  arrange(condition, sample_id)

cat("\nSample-condition mapping:\n")
print(sample_info)

# ==============================================================================
# Pseudobulk Aggregation Function
# ==============================================================================

create_pseudobulk <- function(seurat_obj, cell_type) {
  
  cat(paste("Creating pseudobulk for", cell_type, "...\n"))
  
  # Subset to specific cell type
  subset_obj <- subset(seurat_obj, subset = cell_type_final == cell_type)
  
  # Check if enough cells
  min_cells_per_sample <- 10
  cells_per_sample <- table(subset_obj$sample_id)
  
  if(any(cells_per_sample < min_cells_per_sample)) {
    low_samples <- names(cells_per_sample)[cells_per_sample < min_cells_per_sample]
    cat(paste("  Warning: Low cell counts in samples:", paste(low_samples, collapse = ", "), "\n"))
  }
  
  # Get raw count matrix
  counts <- GetAssayData(subset_obj, layer = "counts")
  
  # Create sample-level metadata
  meta <- subset_obj@meta.data %>%
    select(sample_id, condition, cell_type_final) %>%
    group_by(sample_id, condition) %>%
    summarise(
      n_cells = n(),
      cell_type = first(cell_type_final),
      .groups = 'drop'
    )
  
  # Aggregate counts by sample
  sample_ids <- unique(subset_obj$sample_id)
  pseudobulk_counts <- matrix(0, nrow = nrow(counts), ncol = length(sample_ids))
  rownames(pseudobulk_counts) <- rownames(counts)
  colnames(pseudobulk_counts) <- sample_ids
  
  for(sample in sample_ids) {
    sample_cells <- subset_obj$sample_id == sample
    if(sum(sample_cells) > 0) {
      pseudobulk_counts[, sample] <- Matrix::rowSums(counts[, sample_cells, drop = FALSE])
    }
  }
  
  # Filter out samples with very low cell counts
  keep_samples <- meta$n_cells >= min_cells_per_sample
  pseudobulk_counts <- pseudobulk_counts[, meta$sample_id[keep_samples], drop = FALSE]
  meta <- meta[keep_samples, ]
  
  cat(paste("  Pseudobulk created:", ncol(pseudobulk_counts), "samples,", nrow(pseudobulk_counts), "genes\n"))
  cat(paste("  Cell counts per sample:", paste(meta$n_cells, collapse = ", "), "\n"))
  
  return(list(
    counts = pseudobulk_counts,
    meta = meta,
    cell_type = cell_type
  ))
}

# ==============================================================================
# DESeq2 Analysis Function
# ==============================================================================

run_deseq2_analysis <- function(pseudobulk_data, comparison) {
  
  cat(paste("  Running DESeq2 for", pseudobulk_data$cell_type, ":", comparison, "\n"))
  
  counts <- pseudobulk_data$counts
  meta <- pseudobulk_data$meta
  
  # Set up comparison groups
  if(comparison == "Dep_vs_Naive") {
    keep_conditions <- c("Dependent", "Naive")
    group1 <- "Dependent"
    group2 <- "Naive"
  } else if(comparison == "With_vs_Dep") {
    keep_conditions <- c("Withdrawal", "Dependent")
    group1 <- "Withdrawal"
    group2 <- "Dependent"
  } else if(comparison == "With_vs_Naive") {
    keep_conditions <- c("Withdrawal", "Naive")
    group1 <- "Withdrawal"
    group2 <- "Naive"
  } else {
    stop("Invalid comparison")
  }
  
  # Filter samples for this comparison
  keep_samples <- meta$condition %in% keep_conditions
  if(sum(keep_samples) < 4) {  # Need at least 2 samples per group
    cat(paste("    Skipping - insufficient samples for comparison\n"))
    return(NULL)
  }
  
  counts_filt <- counts[, keep_samples, drop = FALSE]
  meta_filt <- meta[keep_samples, ]
  
  # Check if we have samples in both groups
  group_counts <- table(meta_filt$condition)
  if(any(group_counts < 2)) {
    cat(paste("    Skipping - need at least 2 samples per group\n"))
    return(NULL)
  }
  
  # Create DESeq2 object
  meta_filt$condition <- factor(meta_filt$condition, levels = c(group2, group1))
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts_filt,
    colData = meta_filt,
    design = ~ condition
  )
  
  # Filter lowly expressed genes
  keep_genes <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep_genes, ]
  
  cat(paste("    Testing", nrow(dds), "genes after filtering\n"))
  
  # Run DESeq2
  dds <- DESeq(dds, quiet = TRUE)
  
  # Extract results
  res <- results(dds, contrast = c("condition", group1, group2))
  
  # Convert to data frame and add gene names
  results_df <- as.data.frame(res)
  results_df$gene <- rownames(results_df)
  
  # Add significance categories
  results_df$significance <- "Not Significant"
  results_df$significance[results_df$padj < 0.05 & !is.na(results_df$padj)] <- "Significant"
  results_df$significance[results_df$padj < 0.01 & !is.na(results_df$padj)] <- "Highly Significant"
  
  # Add effect size categories
  results_df$effect_size <- "Small"
  results_df$effect_size[abs(results_df$log2FoldChange) > 0.5 & !is.na(results_df$log2FoldChange)] <- "Medium"
  results_df$effect_size[abs(results_df$log2FoldChange) > 1.0 & !is.na(results_df$log2FoldChange)] <- "Large"
  
  # Add regulation direction
  results_df$regulation <- "Not Significant"
  results_df$regulation[results_df$padj < 0.05 & results_df$log2FoldChange > 0 & !is.na(results_df$padj)] <- "Upregulated"
  results_df$regulation[results_df$padj < 0.05 & results_df$log2FoldChange < 0 & !is.na(results_df$padj)] <- "Downregulated"
  
  # Add metadata
  results_df$cell_type <- pseudobulk_data$cell_type
  results_df$comparison <- comparison
  results_df$group1 <- group1
  results_df$group2 <- group2
  results_df$total_samples <- ncol(counts_filt)
  results_df$group1_samples <- sum(meta_filt$condition == group1)
  results_df$group2_samples <- sum(meta_filt$condition == group2)
  results_df$genes_tested <- nrow(results_df)
  results_df$analysis_date <- Sys.Date()
  
  # Reorder columns
  results_df <- results_df[, c("gene", "log2FoldChange", "baseMean", "stat", "pvalue", "padj",
                              "significance", "effect_size", "regulation",
                              "cell_type", "comparison", "group1", "group2",
                              "total_samples", "group1_samples", "group2_samples", 
                              "genes_tested", "analysis_date")]
  
  # Sort by adjusted p-value
  results_df <- results_df[order(results_df$padj, na.last = TRUE), ]
  
  # Remove rows with NA padj
  results_df <- results_df[!is.na(results_df$padj), ]
  
  cat(paste("    Completed:", nrow(results_df), "genes tested,", 
            sum(results_df$padj < 0.05), "significant (padj < 0.05)\n"))
  
  return(list(
    results = results_df,
    dds = dds
  ))
}

# ==============================================================================
# Run Comprehensive Pseudobulk Analysis
# ==============================================================================

# Define analysis parameters
cell_types <- unique(seurat_obj$cell_type_final)
cell_types <- cell_types[!is.na(cell_types)]
comparisons <- c("Dep_vs_Naive", "With_vs_Dep", "With_vs_Naive")

cat("\nRunning pseudobulk analysis for:\n")
cat("Cell types:", paste(cell_types, collapse = ", "), "\n")
cat("Comparisons:", paste(comparisons, collapse = ", "), "\n\n")

# Initialize results storage
all_pseudobulk_results <- list()
all_dds_objects <- list()
summary_stats <- data.frame()

# Create pseudobulk data for each cell type
pseudobulk_data <- list()
for(cell_type in cell_types) {
  pseudobulk_data[[cell_type]] <- create_pseudobulk(seurat_obj, cell_type)
}

# Run DESeq2 analysis for each cell type and comparison
for(cell_type in cell_types) {
  for(comparison in comparisons) {
    
    analysis_name <- paste(cell_type, comparison, sep = "_")
    
    # Run DESeq2 analysis
    deseq_result <- run_deseq2_analysis(pseudobulk_data[[cell_type]], comparison)
    
    if(!is.null(deseq_result)) {
      # Store results
      all_pseudobulk_results[[analysis_name]] <- deseq_result$results
      all_dds_objects[[analysis_name]] <- deseq_result$dds
      
      # Calculate summary statistics
      results_df <- deseq_result$results
      summary_stats <- rbind(summary_stats, data.frame(
        cell_type = cell_type,
        comparison = comparison,
        total_genes = nrow(results_df),
        significant_genes = sum(results_df$padj < 0.05),
        highly_significant = sum(results_df$padj < 0.01),
        upregulated = sum(results_df$regulation == "Upregulated"),
        downregulated = sum(results_df$regulation == "Downregulated"),
        large_effect = sum(results_df$padj < 0.05 & abs(results_df$log2FoldChange) > 1),
        total_samples = unique(results_df$total_samples),
        group1_samples = unique(results_df$group1_samples),
        group2_samples = unique(results_df$group2_samples),
        analysis_date = Sys.Date(),
        stringsAsFactors = FALSE
      ))
    }
  }
}

cat("\n", rep("=", 60), "\n")
cat("PSEUDOBULK ANALYSIS COMPLETE\n")
cat(rep("=", 60), "\n")

# Display summary
cat("Summary of successful analyses:\n")
print(summary_stats)

# ==============================================================================
# Save Results
# ==============================================================================

cat("\nSaving results...\n")

# Save R objects
saveRDS(all_pseudobulk_results, file.path(output_dir, "pseudobulk_all_results.rds"))
saveRDS(all_dds_objects, file.path(output_dir, "pseudobulk_dds_objects.rds"))
saveRDS(summary_stats, file.path(output_dir, "pseudobulk_summary_stats.rds"))
saveRDS(pseudobulk_data, file.path(output_dir, "pseudobulk_count_data.rds"))

# Save summary table
write.csv(summary_stats, file.path(output_dir, "Tables", "pseudobulk_summary.csv"), row.names = FALSE)

# Save individual results as CSV
for(analysis_name in names(all_pseudobulk_results)) {
  write.csv(all_pseudobulk_results[[analysis_name]], 
            file.path(output_dir, "Tables", paste0(analysis_name, "_pseudobulk.csv")), 
            row.names = FALSE)
}

# Create comprehensive Excel workbook
excel_data <- all_pseudobulk_results
excel_data[["Summary"]] <- summary_stats
write.xlsx(excel_data, file.path(output_dir, "Tables", "pseudobulk_comprehensive_results.xlsx"))

# ==============================================================================
# Generate Quick Report
# ==============================================================================

cat("\n", rep("=", 60), "\n")
cat("PSEUDOBULK ANALYSIS REPORT\n")
cat(rep("=", 60), "\n")

# Overall statistics
successful_analyses <- length(all_pseudobulk_results)

cat("OVERVIEW:\n")
cat("- Planned analyses:", length(cell_types) * length(comparisons), "\n")
cat("- Successful analyses:", successful_analyses, "\n")
cat("- Failed analyses (insufficient samples):", (length(cell_types) * length(comparisons)) - successful_analyses, "\n\n")

# Results by comparison
cat("RESULTS BY COMPARISON:\n")
for(comp in comparisons) {
  comp_data <- summary_stats[summary_stats$comparison == comp, ]
  if(nrow(comp_data) > 0) {
    cat(paste(comp, ":\n"))
    cat(paste("  - Cell types analyzed:", nrow(comp_data), "\n"))
    cat(paste("  - Total significant genes:", sum(comp_data$significant_genes), "\n"))
    cat(paste("  - Avg significant per cell type:", round(mean(comp_data$significant_genes), 1), "\n"))
  }
}

# Top cell types by significant genes
cat("\nTOP CELL TYPES (by significant genes):\n")
if(nrow(summary_stats) > 0) {
  celltype_totals <- summary_stats %>%
    group_by(cell_type) %>%
    summarise(total_significant = sum(significant_genes), 
              analyses = n(), .groups = 'drop') %>%
    arrange(desc(total_significant))
  
  print(head(celltype_totals, 10))
}

# Top genes across all analyses
cat("\n", rep("=", 60), "\n")
cat("TOP GENES ACROSS ALL ANALYSES\n")
cat(rep("=", 60), "\n")

if(length(all_pseudobulk_results) > 0) {
  all_combined <- do.call(rbind, all_pseudobulk_results)
  top_genes <- all_combined %>%
    filter(padj < 0.01) %>%
    group_by(gene) %>%
    summarise(
      times_significant = n(),
      avg_logFC = mean(log2FoldChange),
      min_padj = min(padj),
      .groups = 'drop'
    ) %>%
    arrange(desc(times_significant), min_padj)
  
  cat("Most frequently significant genes (padj < 0.01):\n")
  print(head(top_genes, 15))
}

# Files created
cat("\nFILES CREATED:\n")
cat("- R objects: pseudobulk_all_results.rds, pseudobulk_dds_objects.rds, pseudobulk_summary_stats.rds\n")
cat("- Excel workbook: pseudobulk_comprehensive_results.xlsx\n")
cat("- Summary CSV: pseudobulk_summary.csv\n")
cat("- Individual analysis CSVs: [cell_type]_[comparison]_pseudobulk.csv\n")
cat("- Output directory:", output_dir, "\n")

cat("\nPseudobulk analysis complete! Ready for visualization and pathway analysis.\n")
cat("Script completed successfully!\n")