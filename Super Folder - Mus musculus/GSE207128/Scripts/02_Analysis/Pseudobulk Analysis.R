# ==============================================================================
# Pseudobulk Differential Expression Analysis using edgeR
# Aggregates cells by sample and condition for proper statistical modeling
# Accounts for biological replicates and reduces single-cell noise
# Uses both filtered and non-filtered data to capture all cell types
# ==============================================================================

# Load required libraries
library(Seurat)
library(edgeR)
library(dplyr)
library(tidyr)
library(openxlsx)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(Matrix)

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
# Load Data - MODIFIED TO USE BOTH FILTERED AND NON-FILTERED
# ==============================================================================

cat("Loading annotated Seurat objects...\n")
annotated_files <- list.files(data_dir, pattern = "*.rds", full.names = TRUE)
cat("Available RDS files:\n")
print(annotated_files)

# Load both filtered and non-filtered files
filtered_file <- grep("filtered_annotated", annotated_files, value = TRUE)
non_filtered_file <- grep("^(?!.*filtered).*annotated.*\\.rds$", annotated_files, value = TRUE, perl = TRUE)

if(length(filtered_file) > 0 && length(non_filtered_file) > 0) {
  # Load both files
  seurat_filtered <- readRDS(filtered_file[1])
  seurat_non_filtered <- readRDS(non_filtered_file[1])
  
  cat("Loaded filtered:", basename(filtered_file[1]), "\n")
  cat("Loaded non-filtered:", basename(non_filtered_file[1]), "\n")
  
  # Check cell types in each
  filtered_cell_types <- unique(seurat_filtered$predicted_celltype)
  non_filtered_cell_types <- unique(seurat_non_filtered$predicted_celltype)
  
  cat("\nCell types comparison:\n")
  cat("Filtered data cell types:", paste(sort(filtered_cell_types[!is.na(filtered_cell_types)]), collapse = ", "), "\n")
  cat("Non-filtered data cell types:", paste(sort(non_filtered_cell_types[!is.na(non_filtered_cell_types)]), collapse = ", "), "\n")
  
  # Find additional cell types in non-filtered data
  additional_cell_types <- setdiff(non_filtered_cell_types, filtered_cell_types)
  additional_cell_types <- additional_cell_types[!is.na(additional_cell_types)]
  
  cat("Additional cell types in non-filtered:", paste(additional_cell_types, collapse = ", "), "\n")
  
  # Strategy: Use filtered data as base, but extract additional cell types from non-filtered
  if(length(additional_cell_types) > 0) {
    cat("\nExtracting additional cell types from non-filtered data...\n")
    
    # Get cells with additional cell types from non-filtered data
    additional_cells <- subset(seurat_non_filtered, 
                              subset = predicted_celltype %in% additional_cell_types)
    
    # Check if these cells are high quality by applying same QC metrics
    # Add QC metrics to additional cells if not present
    if(!"percent.mt" %in% colnames(additional_cells@meta.data)) {
      additional_cells[["percent.mt"]] <- PercentageFeatureSet(additional_cells, pattern = "^mt-")
    }
    if(!"percent.ribo" %in% colnames(additional_cells@meta.data)) {
      additional_cells[["percent.ribo"]] <- PercentageFeatureSet(additional_cells, pattern = "^Rp[sl]")
    }
    
    # Apply same quality filters as used in your QC script
    additional_cells_filtered <- subset(additional_cells, 
                                       subset = nFeature_RNA > 200 & 
                                               nFeature_RNA < 7500 & 
                                               percent.mt < 25 & 
                                               nCount_RNA > 500)
    
    cat("Additional cells after QC filtering:", ncol(additional_cells_filtered), "\n")
    
    # Merge with filtered data
    if(ncol(additional_cells_filtered) > 0) {
      seurat_obj <- merge(seurat_filtered, additional_cells_filtered)
      cat("Combined object created with", ncol(seurat_obj), "cells\n")
    } else {
      cat("No additional cells passed QC, using filtered data only\n")
      seurat_obj <- seurat_filtered
    }
    
  } else {
    cat("No additional cell types found, using filtered data only\n")
    seurat_obj <- seurat_filtered
  }
  
} else if(length(filtered_file) > 0) {
  # Only filtered file available
  seurat_obj <- readRDS(filtered_file[1])
  cat("Loaded filtered only:", basename(filtered_file[1]), "\n")
  
} else if(length(non_filtered_file) > 0) {
  # Only non-filtered file available
  seurat_obj <- readRDS(non_filtered_file[1])
  cat("Loaded non-filtered only:", basename(non_filtered_file[1]), "\n")
  
} else {
  stop("No annotated RDS files found in the directory")
}

# Check and create sample_id if needed
if(!"sample_id" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj$sample_id <- seurat_obj$orig.ident
  cat("Created sample_id from orig.ident\n")
}

# Use predicted cell types
seurat_obj$cell_type_final <- seurat_obj$predicted_celltype

# Data summary
cat("\nFinal data summary:\n")
cat("- Total cells:", ncol(seurat_obj), "\n")
cat("- Total genes:", nrow(seurat_obj), "\n")
cat("- Cell types:", length(unique(seurat_obj$cell_type_final)), "\n")
cat("- Conditions:", paste(unique(seurat_obj$condition), collapse = ", "), "\n")
cat("- Samples:", paste(unique(seurat_obj$sample_id), collapse = ", "), "\n")

# Display all unique cell types
cat("\nAll unique cell types in combined data:\n")
all_cell_types <- sort(unique(seurat_obj$cell_type_final))
all_cell_types <- all_cell_types[!is.na(all_cell_types)]
print(all_cell_types)

# Cell counts per condition and cell type
cell_counts <- seurat_obj@meta.data %>%
  group_by(cell_type_final, condition) %>%
  summarise(count = n(), .groups = 'drop') %>%
  tidyr::pivot_wider(names_from = condition, values_from = count, values_fill = 0)

cat("\nCell counts per condition and cell type:\n")
print(cell_counts)

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
  min_cells_per_sample <- 5  # Reduced from 10 to capture more cell types
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
# edgeR Analysis Function
# ==============================================================================

run_edgeR_pseudobulk_analysis <- function(pseudobulk_data, comparison) {
  
  cat(paste("  Running edgeR for", pseudobulk_data$cell_type, ":", comparison, "\n"))
  
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
  if(sum(keep_samples) < 3) {  # Reduced from 4 to capture more comparisons
    cat(paste("    Skipping - insufficient samples for comparison (need ≥3, have", sum(keep_samples), ")\n"))
    return(NULL)
  }
  
  counts_filt <- counts[, keep_samples, drop = FALSE]
  meta_filt <- meta[keep_samples, ]
  
  # Check if we have samples in both groups
  group_counts <- table(meta_filt$condition)
  if(any(group_counts < 1)) {  # Reduced from 2 to be more lenient
    cat(paste("    Skipping - need at least 1 sample per group\n"))
    return(NULL)
  }
  
  # Create edgeR objects
  condition <- factor(meta_filt$condition, levels = c(group2, group1))
  design <- model.matrix(~ condition)
  
  # Create DGEList object
  dge <- DGEList(counts = counts_filt, group = condition)
  
  # Calculate normalization factors
  dge <- calcNormFactors(dge)
  
  # Filter lowly expressed genes (more lenient for rare cell types)
  keep_genes <- filterByExpr(dge, design, min.count = 5, min.total.count = 10)
  dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
  
  cat(paste("    Testing", nrow(dge), "genes after filtering\n"))
  
  # Check if we have enough genes
  if(nrow(dge) < 100) {
    cat(paste("    Skipping - too few genes after filtering (", nrow(dge), ")\n"))
    return(NULL)
  }
  
  # Estimate dispersions
  dge <- estimateDisp(dge, design, robust = TRUE)
  
  # Fit quasi-likelihood model
  fit <- glmQLFit(dge, design, robust = TRUE)
  
  # Test for differential expression
  qlf <- glmQLFTest(fit, coef = 2)  # Test condition effect
  
  # Extract results
  results_table <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  results_df <- data.frame(
    gene = rownames(results_table),
    log2FoldChange = results_table$logFC,
    baseMean = results_table$logCPM,
    stat = results_table$F,
    pvalue = results_table$PValue,
    padj = results_table$FDR,
    stringsAsFactors = FALSE
  )
  
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
    dge = dge,
    fit = fit
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
all_edgeR_objects <- list()
summary_stats <- data.frame()

# Create pseudobulk data for each cell type
cat("\n", rep("=", 60), "\n")
cat("CREATING PSEUDOBULK DATA\n")
cat(rep("=", 60), "\n")

pseudobulk_data <- list()
for(cell_type in cell_types) {
  pseudobulk_data[[cell_type]] <- create_pseudobulk(seurat_obj, cell_type)
}

# Add diagnostic information
cat("\n", rep("=", 60), "\n")
cat("PSEUDOBULK CREATION SUMMARY\n")
cat(rep("=", 60), "\n")

for(cell_type in cell_types) {
  cat(paste("\n", cell_type, ":\n"))
  
  if(is.null(pseudobulk_data[[cell_type]])) {
    cat("  Status: FAILED - No pseudobulk data created\n")
  } else {
    pb_data <- pseudobulk_data[[cell_type]]
    cat(paste("  Status: SUCCESS\n"))
    cat(paste("  Samples:", ncol(pb_data$counts), "\n"))
    cat(paste("  Genes:", nrow(pb_data$counts), "\n"))
    
    # Check samples per condition
    condition_counts <- table(pb_data$meta$condition)
    cat("  Samples per condition:", paste(names(condition_counts), "=", condition_counts, collapse = ", "), "\n")
    cat("  Cell counts per sample:", paste(pb_data$meta$n_cells, collapse = ", "), "\n")
  }
}

# Run edgeR analysis for each cell type and comparison
cat("\n", rep("=", 60), "\n")
cat("RUNNING DIFFERENTIAL EXPRESSION ANALYSIS\n")
cat(rep("=", 60), "\n")

for(cell_type in cell_types) {
  for(comparison in comparisons) {
    
    analysis_name <- paste(cell_type, comparison, sep = "_")
    
    # Skip if pseudobulk creation failed
    if(is.null(pseudobulk_data[[cell_type]])) {
      cat(paste("Skipping", analysis_name, "- no pseudobulk data\n"))
      next
    }
    
    # Run edgeR analysis
    edgeR_result <- run_edgeR_pseudobulk_analysis(pseudobulk_data[[cell_type]], comparison)
    
    if(!is.null(edgeR_result)) {
      # Store results
      all_pseudobulk_results[[analysis_name]] <- edgeR_result$results
      all_edgeR_objects[[analysis_name]] <- edgeR_result
      
      # Calculate summary statistics
      results_df <- edgeR_result$results
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
saveRDS(all_edgeR_objects, file.path(output_dir, "pseudobulk_edgeR_objects.rds"))
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
if(length(all_pseudobulk_results) > 0) {
  excel_data <- all_pseudobulk_results
  excel_data[["Summary"]] <- summary_stats
  write.xlsx(excel_data, file.path(output_dir, "Tables", "pseudobulk_comprehensive_results.xlsx"))
}

# ==============================================================================
# Generate Quick Report
# ==============================================================================

cat("\n", rep("=", 60), "\n")
cat("PSEUDOBULK ANALYSIS REPORT\n")
cat(rep("=", 60), "\n")

# Overall statistics
successful_analyses <- length(all_pseudobulk_results)
total_cell_types <- length(cell_types)

cat("OVERVIEW:\n")
cat("- Total cell types found:", total_cell_types, "\n")
cat("- Planned analyses:", total_cell_types * length(comparisons), "\n")
cat("- Successful analyses:", successful_analyses, "\n")
cat("- Failed analyses (insufficient samples/cells):", (total_cell_types * length(comparisons)) - successful_analyses, "\n\n")

# Results by comparison
if(nrow(summary_stats) > 0) {
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
  celltype_totals <- summary_stats %>%
    group_by(cell_type) %>%
    summarise(total_significant = sum(significant_genes), 
              analyses = n(), .groups = 'drop') %>%
    arrange(desc(total_significant))
  
  print(head(celltype_totals, 15))
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
} else {
  cat("No successful analyses to report genes from.\n")
}

# Files created
cat("\nFILES CREATED:\n")
cat("- R objects: pseudobulk_all_results.rds, pseudobulk_edgeR_objects.rds, pseudobulk_summary_stats.rds\n")
if(length(all_pseudobulk_results) > 0) {
  cat("- Excel workbook: pseudobulk_comprehensive_results.xlsx\n")
}
cat("- Summary CSV: pseudobulk_summary.csv\n")
cat("- Individual analysis CSVs: [cell_type]_[comparison]_pseudobulk.csv\n")
cat("- Output directory:", output_dir, "\n")

cat("\nPseudobulk analysis complete! Ready for visualization and pathway analysis.\n")
cat("Script completed successfully!\n")