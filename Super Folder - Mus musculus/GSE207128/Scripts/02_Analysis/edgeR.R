# ==============================================================================
# Differential Expression Analysis using edgeR
# General genome-wide analysis for all contrasts
# Following original study methodology (edgeRQLFDetRate)
# Data is already batch-corrected via Harmony integration
# ==============================================================================

# Load required libraries
library(Seurat)
library(edgeR)
library(dplyr)
library(openxlsx)
library(tidyr)

# Set working directory and paths
base_dir <- "/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus/GSE207128"
data_dir <- file.path(base_dir, "Outputs/04_Annotated_Data")
output_dir <- file.path(base_dir, "Outputs/05_Analysis_Results/DEG_Analysis")

# Create output directories
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "Tables"), recursive = TRUE, showWarnings = FALSE)

cat("Starting edgeR differential expression analysis...\n")
cat("Date:", Sys.Date(), "\n")
cat("Time:", Sys.time(), "\n\n")

# ==============================================================================
# Load Data
# ==============================================================================

cat("Loading annotated Seurat object...\n")
annotated_files <- list.files(data_dir, pattern = "*.rds", full.names = TRUE)
cat("Available RDS files in annotated data directory:\n")
print(annotated_files)

# Load filtered annotated file (has cell type predictions and quality filtering)
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

# Data summary
cat("\nData summary:\n")
cat("- Total cells:", ncol(seurat_obj), "\n")
cat("- Total genes:", nrow(seurat_obj), "\n")
cat("- Cell types:", length(unique(seurat_obj$predicted_celltype)), "\n")
cat("- Conditions:", paste(unique(seurat_obj$condition), collapse = ", "), "\n")

# Use predicted cell types from annotation pipeline
seurat_obj$cell_type_final <- seurat_obj$predicted_celltype
cat("Using predicted_celltype from annotation pipeline\n")

# Display unique cell types
cat("\nUnique cell types:\n")
print(sort(unique(seurat_obj$cell_type_final)))

# Cell counts per condition and cell type
cell_counts <- seurat_obj@meta.data %>%
  group_by(cell_type_final, condition) %>%
  summarise(count = n(), .groups = 'drop') %>%
  tidyr::pivot_wider(names_from = condition, values_from = count, values_fill = 0)

cat("\nCell counts per condition and cell type:\n")
print(cell_counts)
cat("\n")

# ==============================================================================
# edgeR Analysis Function - FIXED VERSION
# ==============================================================================

run_edgeR_analysis <- function(seurat_obj, cell_type, comparison, min_cells = 50) {
  
  cat(paste("Analyzing", cell_type, ":", comparison, "\n"))
  
  # Subset to specific cell type
  subset_obj <- subset(seurat_obj, subset = cell_type_final == cell_type)
  
  # Check if enough cells
  if(ncol(subset_obj) < min_cells) {
    cat(paste("  Skipping - insufficient cells (<", min_cells, ")\n"))
    return(NULL)
  }
  
  # Set up comparison groups
  if(comparison == "Dep_vs_Naive") {
    subset_obj <- subset(subset_obj, subset = condition %in% c("Dependent", "Naive"))
    group1 <- "Dependent"
    group2 <- "Naive"
  } else if(comparison == "With_vs_Dep") {
    subset_obj <- subset(subset_obj, subset = condition %in% c("Withdrawal", "Dependent"))
    group1 <- "Withdrawal"
    group2 <- "Dependent"
  } else if(comparison == "With_vs_Naive") {
    subset_obj <- subset(subset_obj, subset = condition %in% c("Withdrawal", "Naive"))
    group1 <- "Withdrawal"
    group2 <- "Naive"
  } else {
    stop("Invalid comparison. Use 'Dep_vs_Naive', 'With_vs_Dep', or 'With_vs_Naive'")
  }
  
  # Check if both groups have sufficient cells
  group_counts <- table(subset_obj$condition)
  cat(paste("  Cell counts -", group1, ":", group_counts[group1], 
            group2, ":", group_counts[group2], "\n"))
  
  if(any(group_counts < 10)) {
    cat(paste("  Skipping - insufficient cells in one group (<10)\n"))
    return(NULL)
  }
  
  # Get count matrix - USE ORIGINAL COUNTS
  counts <- GetAssayData(subset_obj, layer = "counts")
  
  # Create experimental design - SIMPLE DESIGN (data already batch-corrected)
  condition <- factor(subset_obj$condition, levels = c(group2, group1))
  design <- model.matrix(~ condition)
  cat(paste("  Using simple design (data already batch-corrected via Harmony)\n"))
  
  # Print design matrix column names for debugging
  cat(paste("  Design matrix columns:", paste(colnames(design), collapse = ", "), "\n"))
  
  # Create DGEList object
  dge <- DGEList(counts = counts, group = condition)
  
  # Calculate normalization factors
  dge <- calcNormFactors(dge)
  
  # Filter lowly expressed genes
  keep <- filterByExpr(dge, design, min.count = 5, min.total.count = 15)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  cat(paste("  Testing", nrow(dge), "genes after filtering\n"))
  
  # Estimate dispersions
  dge <- estimateDisp(dge, design, robust = TRUE)
  
  # Fit quasi-likelihood model
  fit <- glmQLFit(dge, design, robust = TRUE)
  
  # Test for differential expression - FIXED: USE COEFFICIENT NUMBER
  qlf <- glmQLFTest(fit, coef = 2)  # Test second coefficient (condition effect)
  
  # Extract results
  results <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  results$gene <- rownames(results)
  
  # Add significance categories
  results$significance <- "Not Significant"
  results$significance[results$FDR < 0.05] <- "Significant"
  results$significance[results$FDR < 0.01] <- "Highly Significant"
  
  # Add effect size categories
  results$effect_size <- "Small"
  results$effect_size[abs(results$logFC) > 0.5] <- "Medium" 
  results$effect_size[abs(results$logFC) > 1.0] <- "Large"
  
  # Add regulation direction
  results$regulation <- "Not Significant"
  results$regulation[results$FDR < 0.05 & results$logFC > 0] <- "Upregulated"
  results$regulation[results$FDR < 0.05 & results$logFC < 0] <- "Downregulated"
  
  # Add metadata
  results$cell_type <- cell_type
  results$comparison <- comparison
  results$group1 <- group1
  results$group2 <- group2
  results$total_cells <- ncol(subset_obj)
  results$group1_cells <- sum(subset_obj$condition == group1)
  results$group2_cells <- sum(subset_obj$condition == group2)
  results$genes_tested <- nrow(results)
  results$analysis_date <- Sys.Date()
  
  # Reorder columns
  results <- results[, c("gene", "logFC", "logCPM", "F", "PValue", "FDR", 
                        "significance", "effect_size", "regulation",
                        "cell_type", "comparison", "group1", "group2", 
                        "total_cells", "group1_cells", "group2_cells", 
                        "genes_tested", "analysis_date")]
  
  # Sort by FDR
  results <- results[order(results$FDR), ]
  
  cat(paste("  Completed:", nrow(results), "genes tested,", 
            sum(results$FDR < 0.05), "significant (FDR < 0.05)\n"))
  
  return(results)
}

# ==============================================================================
# Run Comprehensive Analysis
# ==============================================================================

# Define analysis parameters
cell_types <- unique(seurat_obj$cell_type_final)
cell_types <- cell_types[!is.na(cell_types)]
comparisons <- c("Dep_vs_Naive", "With_vs_Dep", "With_vs_Naive")

cat("Running analysis for:\n")
cat("Cell types:", paste(cell_types, collapse = ", "), "\n")
cat("Comparisons:", paste(comparisons, collapse = ", "), "\n\n")

# Initialize results storage
all_results <- list()
summary_stats <- data.frame()

# Run analysis for each cell type and comparison
for(cell_type in cell_types) {
  for(comparison in comparisons) {
    
    analysis_name <- paste(cell_type, comparison, sep = "_")
    
    # Run edgeR analysis
    results <- run_edgeR_analysis(seurat_obj, cell_type, comparison)
    
    if(!is.null(results)) {
      # Store results
      all_results[[analysis_name]] <- results
      
      # Calculate summary statistics
      summary_stats <- rbind(summary_stats, data.frame(
        cell_type = cell_type,
        comparison = comparison,
        total_genes = nrow(results),
        significant_genes = sum(results$FDR < 0.05),
        highly_significant = sum(results$FDR < 0.01),
        upregulated = sum(results$regulation == "Upregulated"),
        downregulated = sum(results$regulation == "Downregulated"),
        large_effect = sum(results$FDR < 0.05 & abs(results$logFC) > 1),
        total_cells = unique(results$total_cells),
        group1_cells = unique(results$group1_cells),
        group2_cells = unique(results$group2_cells),
        analysis_date = Sys.Date(),
        stringsAsFactors = FALSE
      ))
    }
  }
}

cat("\n" , rep("=", 60), "\n")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 60), "\n")

# Display summary
cat("Summary of successful analyses:\n")
print(summary_stats)

# ==============================================================================
# Save Results
# ==============================================================================

cat("\nSaving results...\n")

# Save R objects
saveRDS(all_results, file.path(output_dir, "edgeR_all_results.rds"))
saveRDS(summary_stats, file.path(output_dir, "edgeR_summary_stats.rds"))

# Save summary table
write.csv(summary_stats, file.path(output_dir, "Tables", "edgeR_summary.csv"), row.names = FALSE)

# Save individual results as CSV
for(analysis_name in names(all_results)) {
  write.csv(all_results[[analysis_name]], 
            file.path(output_dir, "Tables", paste0(analysis_name, "_edgeR.csv")), 
            row.names = FALSE)
}

# Create comprehensive Excel workbook
excel_data <- all_results
excel_data[["Summary"]] <- summary_stats
write.xlsx(excel_data, file.path(output_dir, "Tables", "edgeR_comprehensive_results.xlsx"))

# ==============================================================================
# Generate Quick Report
# ==============================================================================

cat("\n" , rep("=", 60), "\n")
cat("QUICK ANALYSIS REPORT\n")
cat(rep("=", 60), "\n")

# Overall statistics
successful_analyses <- length(all_results)

cat("OVERVIEW:\n")
cat("- Planned analyses:", length(cell_types) * length(comparisons), "\n")
cat("- Successful analyses:", successful_analyses, "\n")
cat("- Failed analyses (insufficient cells):", (length(cell_types) * length(comparisons)) - successful_analyses, "\n\n")

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
celltype_totals <- summary_stats %>%
  group_by(cell_type) %>%
  summarise(total_significant = sum(significant_genes), 
            analyses = n(), .groups = 'drop') %>%
  arrange(desc(total_significant))

print(head(celltype_totals, 10))

# Top genes across all analyses
cat("\n" , rep("=", 60), "\n")
cat("TOP GENES ACROSS ALL ANALYSES\n")
cat(rep("=", 60), "\n")

if(length(all_results) > 0) {
  all_combined <- do.call(rbind, all_results)
  top_genes <- all_combined %>%
    filter(FDR < 0.01) %>%
    group_by(gene) %>%
    summarise(
      times_significant = n(),
      avg_logFC = mean(logFC),
      min_FDR = min(FDR),
      .groups = 'drop'
    ) %>%
    arrange(desc(times_significant), min_FDR)
  
  cat("Most frequently significant genes (FDR < 0.01):\n")
  print(head(top_genes, 15))
}

# Files created
cat("\nFILES CREATED:\n")
cat("- R objects: edgeR_all_results.rds, edgeR_summary_stats.rds\n")
cat("- Excel workbook: edgeR_comprehensive_results.xlsx\n")
cat("- Summary CSV: edgeR_summary.csv\n")
cat("- Individual analysis CSVs: [cell_type]_[comparison]_edgeR.csv\n")
cat("- Output directory:", output_dir, "\n")

cat("\nAnalysis complete! Results ready for downstream pathway analysis.\n")
cat("Script completed successfully!\n")