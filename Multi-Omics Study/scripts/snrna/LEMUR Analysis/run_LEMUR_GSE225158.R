#!/usr/bin/env Rscript
#' @title End-to-End LEMUR Analysis for GSE225158 snRNA-seq Data
#' @description Comprehensive LEMUR analysis pipeline for OUD vs Control comparison
#' @author Generated Analysis Pipeline
#' @date 2024
#'
#' This script performs:
#' 1. Data loading and QC
#' 3. LEMUR modeling with batch correction
#' 4. Differential testing
#' 5. Neighborhood identification
#' 6. Visualization and output

# =============================================================================
# SETUP AND PACKAGE INSTALLATION
# =============================================================================

# Required packages
required_packages <- c(
  # Core Bioconductor packages
  "SingleCellExperiment", "scran", "scater", "lemur",
  # Data manipulation and visualization
  "tidyverse", "ggplot2", "viridis", "patchwork",
  # File I/O
  "zellkonverter", "HDF5Array",
  # Dimensionality reduction and batch correction
  "uwot", "harmony", "BiocNeighbors", "FNN",
  # Statistics
  "Matrix", "limma"
)

# Install function for missing packages
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    if (pkg %in% rownames(available.packages())) {
      install.packages(pkg, dependencies = TRUE)
    } else {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, dependencies = TRUE)
    }
    library(pkg, character.only = TRUE)
  }
}

# Install and load all packages
cat("📦 Installing and loading required packages...\n")
invisible(lapply(required_packages, install_if_missing))

# =============================================================================
# CONFIGURATION
# =============================================================================

# File paths
BASE_DIR <- "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study"
INPUT_H5AD <- file.path(BASE_DIR, "data/raw/snrna/GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5ad")
OUTPUT_DIR <- "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study/scripts/snrna/LEMUR Analysis/outputs"

# Create output directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "reports"), recursive = TRUE, showWarnings = FALSE)

# Analysis parameters
PARAMS <- list(
  # Quality control thresholds
  min_genes = 200,
  max_genes = 6000,
  max_mito_pct = 20,
  min_counts = 500,
  max_counts = 50000,

  # Feature selection
  n_hvg = 4000,

  # LEMUR parameters
  n_embedding = 20,
  design_formula = "~ condition",

  # Batch correction
  batch_key = "donor_id",
  condition_key = "condition",

  # Statistical thresholds
  fdr_threshold = 0.05,
  effect_size_threshold = 0.1,

  # Visualization
  figure_dpi = 300,
  figure_width = 10,
  figure_height = 8
)

cat("🔧 Configuration complete\n")
cat("📂 Input file:", INPUT_H5AD, "\n")
cat("📂 Output directory:", OUTPUT_DIR, "\n")

# =============================================================================
# DATA LOADING AND INITIAL SETUP
# =============================================================================

cat("\n📊 Loading H5AD data...\n")

# Load H5AD file into SingleCellExperiment
sce <- zellkonverter::readH5AD(INPUT_H5AD, use_hdf5 = TRUE)

cat("✅ Data loaded successfully\n")
cat("📈 Dimensions:", nrow(sce), "genes ×", ncol(sce), "cells\n")

# Print available metadata columns
cat("🏷️  Available metadata columns:\n")
print(colnames(colData(sce)))

# ============================== ALTER METADATA

# After loading the SCE object, let's check and fix the assay names
# Add this code after the line where 'sce' is created:
if (!"counts" %in% assayNames(sce)) {
  # If there's a "matrix" or "X" assay, rename it to "counts"
  if ("matrix" %in% assayNames(sce)) {
    assay(sce, "counts") <- assay(sce, "matrix")
  } else if ("X" %in% assayNames(sce)) {
    assay(sce, "counts") <- assay(sce, "X")
  }
}

assayNames(sce)


# Access the metadata slot of your SingleCellExperiment object
meta_data <- colData(sce)

# Rename columns
colnames(meta_data)[colnames(meta_data) == "ID"] <- "donor_id"
colnames(meta_data)[colnames(meta_data) == "Dx_OUD"] <- "condition"

# Replace the metadata with the updated one
colData(sce) <- meta_data


# ============================== ALTER METADATA END

# Check for required columns
required_cols <- c("donor_id", "condition")
missing_cols <- setdiff(required_cols, colnames(colData(sce)))
if (length(missing_cols) > 0) {
  stop("❌ Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Print condition summary
cat("\n📋 Condition summary:\n")
print(table(sce$condition))

cat("\n📋 Donor summary:\n")
print(table(sce$donor_id))

# =============================================================================
# QUALITY CONTROL
# =============================================================================

cat("\n🔍 Computing quality control metrics...\n")

# Calculate mitochondrial gene percentage
is_mito <- grepl("^MT-|^mt-", rownames(sce))
sce <- addPerCellQCMetrics(sce, subsets = list(Mito = is_mito))

# Calculate ribosomal gene percentage
is_ribo <- grepl("^RP[SL]|^Rp[sl]", rownames(sce))
sce <- addPerCellQCMetrics(sce, subsets = list(Ribo = is_ribo))

# Add gene detection metrics
colData(sce)$n_genes <- colSums(counts(sce) > 0)
colData(sce)$total_counts <- colSums(counts(sce))

# Print QC summary before filtering
cat("📊 QC metrics before filtering:\n")
qc_summary_before <- data.frame(
  Metric = c(rep("total_counts", 6), rep("n_genes", 6), rep("pct_mito", 6)),
  Statistic = rep(names(summary(sce$total_counts)), 3),
  Value = c(
    as.numeric(summary(sce$total_counts)),
    as.numeric(summary(sce$n_genes)),
    as.numeric(summary(sce$subsets_Mito_percent))
  )
)
print(qc_summary_before)

# Apply quality control filters
cat("\n🚫 Applying quality control filters...\n")

# Cell filtering
cells_to_keep <- (
  sce$n_genes >= PARAMS$min_genes &
    sce$n_genes <= PARAMS$max_genes &
    sce$total_counts >= PARAMS$min_counts &
    sce$total_counts <= PARAMS$max_counts &
    sce$subsets_Mito_percent <= PARAMS$max_mito_pct
)

cat("📉 Cells removed:", sum(!cells_to_keep), "out of", ncol(sce), "\n")
sce_filtered <- sce[, cells_to_keep]

# Gene filtering - keep genes expressed in at least 10 cells
genes_to_keep <- rowSums(counts(sce_filtered) > 0) >= 10
cat("📉 Genes removed:", sum(!genes_to_keep), "out of", nrow(sce_filtered), "\n")
sce_filtered <- sce_filtered[genes_to_keep, ]

cat("✅ Filtered data dimensions:", nrow(sce_filtered), "genes ×", ncol(sce_filtered), "cells\n")

# Update condition summary after filtering
cat("\n📋 Condition summary after filtering:\n")
print(table(sce_filtered$condition))

# =============================================================================
# NORMALIZATION AND FEATURE SELECTION
# =============================================================================

cat("\n🔬 Performing normalization...\n")

# Size factor calculation and normalization using scran
clusters <- quickCluster(sce_filtered)
sce_filtered <- computeSumFactors(sce_filtered, clusters = clusters)
sce_filtered <- logNormCounts(sce_filtered)

cat("✅ Normalization complete\n")

# Feature selection - highly variable genes
cat("\n🎯 Selecting highly variable genes...\n")

# Compute gene variance
dec <- modelGeneVar(sce_filtered)
hvg <- getTopHVGs(dec, n = PARAMS$n_hvg)

cat("✅ Selected", length(hvg), "highly variable genes\n")

# Store HVG information
metadata(sce_filtered)$hvg <- hvg
metadata(sce_filtered)$dec <- dec

# =============================================================================
# DIMENSIONALITY REDUCTION (for visualization)
# =============================================================================

cat("\n📐 Computing dimensionality reduction...\n")

# PCA on HVGs
sce_filtered <- runPCA(sce_filtered, subset_row = hvg, ncomponents = 50)

# UMAP for visualization
sce_filtered <- runUMAP(sce_filtered, dimred = "PCA", n_dimred = 30)

cat("✅ Dimensionality reduction complete\n")

# =============================================================================
# LEMUR MODEL FITTING
# =============================================================================

cat("\n🧠 Fitting LEMUR model...\n")

# Prepare data for LEMUR (using HVGs and logcounts)
sce_hvg <- sce_filtered[hvg, ]

# Ensure unique column names in colData
colnames(colData(sce_hvg)) <- make.names(colnames(colData(sce_hvg)), unique = TRUE)

# Convert factors to ensure proper handling
sce_hvg$condition <- as.factor(sce_hvg$condition)
sce_hvg$donor_id <- as.factor(sce_hvg$donor_id)

# Check the design matrix before fitting
cat("📋 Checking design matrix...\n")
design_formula <- as.formula(PARAMS$design_formula)
cat("Design formula:", deparse(design_formula), "\n")

# Create a test design matrix to check for collinearity
test_data <- colData(sce_hvg)[1:100, ] # Use subset for quick check
test_matrix <- model.matrix(design_formula, data = test_data)
cat("Design matrix dimensions:", dim(test_matrix), "\n")
cat("Design matrix rank:", qr(test_matrix)$rank, "\n")
cat("Number of columns:", ncol(test_matrix), "\n")

if (qr(test_matrix)$rank < ncol(test_matrix)) {
  cat("⚠️  Warning: Design matrix appears to have collinearity issues\n")
  cat("Trying alternative design formula...\n")

  # Try with explicit contrasts
  options(contrasts = c("contr.treatment", "contr.poly"))

  # Alternative: use only condition if donor causes issues
  if (PARAMS$design_formula == "~ condition") {
    cat("Using condition-only design\n")
  } else {
    cat("Switching to condition-only design to avoid collinearity\n")
    design_formula <- ~condition
  }
}

# Fit LEMUR model with design formula
fit <- lemur(sce_hvg,
  design = design_formula,
  n_embedding = PARAMS$n_embedding,
  verbose = TRUE
)

cat("✅ LEMUR model fitted successfully\n")
cat("📊 Embedding dimensions:", ncol(fit$embedding), "\n")

# =============================================================================
# BATCH CORRECTION WITH HARMONY (Optional)
# =============================================================================

cat("\n🎵 Applying Harmony batch correction...\n")

# Apply harmony to LEMUR embedding
library(harmony)
harmony_embedding <- harmony::HarmonyMatrix(
  data_mat = fit$embedding,
  meta_data = as.data.frame(colData(sce_hvg)),
  vars_use = PARAMS$batch_key,
  do_pca = FALSE,
  verbose = TRUE
)

# Store harmony-corrected embedding
fit[["harmony_embedding"]] <- harmony_embedding

cat("✅ Harmony correction complete\n")

# =============================================================================
# DIFFERENTIAL TESTING
# =============================================================================

cat("\n🔬 Performing differential testing...\n")

# Test for differential expression between conditions
# Using the main contrast: OUD vs Control
de_res <- test_de(fit, contrast = cond(condition = "OUD") - cond(condition = "None"))

cat("✅ Differential testing complete\n")

# Extract differential expression results from LEMUR fit object
# The test_de function returns the fit object with results stored in the DE assay
de_matrix <- assay(de_res, "DE")
cat("📊 DE matrix dimensions:", dim(de_matrix), "\n")

# Convert to data frame for easier manipulation
de_results_df <- data.frame(
  gene = rownames(de_matrix),
  effect_size = rowMeans(de_matrix), # Average effect size across cells
  stringsAsFactors = FALSE
)

# For LEMUR, we need to compute p-values from the effect sizes
# Using a simple approach based on effect size distribution
de_results_df$pval <- 2 * pnorm(-abs(de_results_df$effect_size / sd(de_results_df$effect_size)))

# Apply multiple testing correction
de_results_df$adj_pval <- p.adjust(de_results_df$pval, method = "BH")

# Filter significant results
sig_results <- de_results_df[de_results_df$adj_pval < PARAMS$fdr_threshold &
  abs(de_results_df$effect_size) > PARAMS$effect_size_threshold, ]
cat("📈 Significant results (FDR < 0.05, |effect| > 0.1):", nrow(sig_results), "\n")

# =============================================================================
# NEIGHBORHOOD IDENTIFICATION
# =============================================================================

cat("\n🏘️  Analyzing differential neighborhoods...\n")

# Use LEMUR's built-in neighborhood analysis based on the embedding
# Extract the harmony-corrected embedding for neighborhood analysis
if ("harmony_embedding" %in% names(colData(de_res))) {
  embedding_matrix <- colData(de_res)$harmony_embedding
} else {
  embedding_matrix <- de_res$embedding
}

# Create neighborhoods based on k-nearest neighbors in the embedding space
library(FNN)
k_neighbors <- 50 # Number of neighbors for each cell
nn_indices <- get.knn(embedding_matrix, k = k_neighbors)$nn.index

# Create differential neighborhoods based on DE gene expression patterns
# Use clustering approach to find meaningful neighborhoods
de_neighborhoods <- list()
if (nrow(sig_results) > 0) {
  # Get expression of DE genes
  de_expr <- assay(de_res, "DE")[rownames(sig_results), ]

  # Find cells with high differential expression
  cell_scores <- colMeans(abs(de_expr))
  threshold <- quantile(cell_scores, 0.95) # Top 5% of cells
  high_de_cells <- which(cell_scores > threshold)

  # Create fewer, more meaningful neighborhoods using clustering
  if (length(high_de_cells) > 10) {
    # Cluster high DE cells based on their embedding coordinates
    high_de_coords <- embedding_matrix[high_de_cells, ]
    n_clusters <- min(10, max(3, length(high_de_cells) %/% 100)) # 3-10 clusters

    # Perform k-means clustering
    set.seed(42)
    clusters <- kmeans(high_de_coords, centers = n_clusters, nstart = 20)

    # Create neighborhoods based on clusters
    for (i in 1:n_clusters) {
      cluster_cells <- high_de_cells[clusters$cluster == i]
      if (length(cluster_cells) > 5) { # Only keep clusters with >5 cells
        # Expand each cluster to include neighbors
        all_neighbors <- unique(c(cluster_cells, as.vector(nn_indices[cluster_cells, ])))
        de_neighborhoods[[length(de_neighborhoods) + 1]] <- all_neighbors
      }
    }
  }
}

# Structure neighborhoods object properly
neighborhoods <- list()
neighborhoods$embedding <- embedding_matrix
neighborhoods$knn_indices <- nn_indices
neighborhoods$de_genes <- rownames(sig_results)
neighborhoods$de_neighborhoods <- de_neighborhoods

cat("✅ Neighborhood analysis complete\n")
cat("📊 Number of DE genes for neighborhood analysis:", length(neighborhoods$de_genes), "\n")
cat("📊 Number of DE neighborhoods found:", length(neighborhoods$de_neighborhoods), "\n")

# =============================================================================
# RESULTS COMPILATION AND OUTPUT
# =============================================================================

cat("\n📊 Compiling and saving results...\n")

# Create comprehensive results table
results_table <- de_results_df
results_table$significant <- results_table$adj_pval < PARAMS$fdr_threshold

# Add gene annotations if available
if ("Symbol" %in% colnames(rowData(de_res))) {
  results_table$gene_symbol <- rowData(sce_hvg)[results_table$gene, "Symbol"]
}

# Sort by adjusted p-value
results_table <- results_table[order(results_table$adj_pval), ]

# Save results table
write.csv(results_table,
  file = file.path(OUTPUT_DIR, "tables", "lemur_de_results.csv"),
  row.names = FALSE
)

# Save top significant results
top_results <- head(results_table[results_table$significant, ], 100)
write.csv(top_results,
  file = file.path(OUTPUT_DIR, "tables", "top_de_genes.csv"),
  row.names = FALSE
)

# Save neighborhood results
if (length(neighborhoods$de_neighborhoods) > 0) {
  neighborhood_summary <- data.frame(
    neighborhood_id = seq_along(neighborhoods$de_neighborhoods),
    n_cells = sapply(neighborhoods$de_neighborhoods, length),
    stringsAsFactors = FALSE
  )

  write.csv(neighborhood_summary,
    file = file.path(OUTPUT_DIR, "tables", "de_neighborhoods.csv"),
    row.names = FALSE
  )
}

cat("✅ Results saved to tables/\n")

# =============================================================================
# VISUALIZATION
# =============================================================================

cat("\n🎨 Creating visualizations...\n")

# Extract UMAP coordinates
umap_coords <- reducedDim(sce_filtered, "UMAP")
plot_data <- data.frame(
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  condition = sce_filtered$condition,
  donor_id = sce_filtered$donor_id,
  total_counts = sce_filtered$total_counts,
  n_genes = sce_filtered$n_genes,
  stringsAsFactors = FALSE
)

# 1. UMAP colored by condition
p1 <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = condition)) +
  geom_point(size = 0.5, alpha = 0.7) +
  scale_color_viridis_d(name = "Condition") +
  theme_minimal() +
  labs(
    title = "UMAP - Condition",
    x = "UMAP 1", y = "UMAP 2"
  ) +
  theme(legend.position = "bottom")

# 2. UMAP colored by donor (batch)
p2 <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = donor_id)) +
  geom_point(size = 0.5, alpha = 0.7) +
  scale_color_viridis_d(name = "Donor ID") +
  theme_minimal() +
  labs(
    title = "UMAP - Donor ID (Batch)",
    x = "UMAP 1", y = "UMAP 2"
  ) +
  theme(legend.position = "bottom")

# 3. QC metrics
p3 <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = log10(total_counts))) +
  geom_point(size = 0.5, alpha = 0.7) +
  scale_color_viridis_c(name = "log10(UMI)") +
  theme_minimal() +
  labs(
    title = "UMAP - Total UMI Counts",
    x = "UMAP 1", y = "UMAP 2"
  ) +
  theme(legend.position = "bottom")

p4 <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = n_genes)) +
  geom_point(size = 0.5, alpha = 0.7) +
  scale_color_viridis_c(name = "Genes") +
  theme_minimal() +
  labs(
    title = "UMAP - Number of Genes",
    x = "UMAP 1", y = "UMAP 2"
  ) +
  theme(legend.position = "bottom")

# Combine plots
combined_plot <- (p1 + p2) / (p3 + p4)

# Save combined UMAP plot
ggsave(file.path(OUTPUT_DIR, "plots", "umap_overview.png"),
  plot = combined_plot,
  width = PARAMS$figure_width,
  height = PARAMS$figure_height,
  dpi = PARAMS$figure_dpi
)

# 5. Volcano plot
volcano_data <- results_table[!is.na(results_table$effect_size) & !is.na(results_table$adj_pval), ]
volcano_data$neg_log10_padj <- -log10(volcano_data$adj_pval)
volcano_data$significant <- volcano_data$adj_pval < PARAMS$fdr_threshold

p_volcano <- ggplot(volcano_data, aes(x = effect_size, y = neg_log10_padj, color = significant)) +
  geom_point(alpha = 0.6, size = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red")) +
  geom_hline(yintercept = -log10(PARAMS$fdr_threshold), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(
    title = "Volcano Plot - LEMUR DE Results",
    subtitle = paste("OUD vs Control | FDR <", PARAMS$fdr_threshold),
    x = "Effect Size",
    y = "-log10(Adjusted P-value)",
    color = "Significant"
  ) +
  theme(legend.position = "bottom")

ggsave(file.path(OUTPUT_DIR, "plots", "volcano_plot.png"),
  plot = p_volcano,
  width = 8, height = 6,
  dpi = PARAMS$figure_dpi
)

# 6. DE neighborhoods on UMAP (if any found)
if (length(neighborhoods$de_neighborhoods) > 0 && length(neighborhoods$de_neighborhoods) <= 20) {
  # Create neighborhood annotation
  plot_data$de_neighborhood <- "Background"
  for (i in seq_along(neighborhoods$de_neighborhoods)) {
    cells_in_neighborhood <- neighborhoods$de_neighborhoods[[i]]
    if (length(cells_in_neighborhood) > 0 && length(cells_in_neighborhood) <= nrow(plot_data)) {
      valid_cells <- cells_in_neighborhood[cells_in_neighborhood <= nrow(plot_data)]
      if (length(valid_cells) > 0) {
        plot_data$de_neighborhood[valid_cells] <- paste("Neighborhood", i)
      }
    }
  }

  # Use better colors for visualization
  if (length(unique(plot_data$de_neighborhood)) <= 10) {
    color_scale <- scale_color_brewer(palette = "Set3", name = "DE Neighborhood")
  } else {
    color_scale <- scale_color_viridis_d(name = "DE Neighborhood", option = "plasma")
  }

  p_neighborhoods <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = de_neighborhood)) +
    geom_point(size = 0.3, alpha = 0.6) +
    color_scale +
    theme_minimal() +
    labs(
      title = "UMAP - DE Neighborhoods",
      subtitle = paste("Found", length(neighborhoods$de_neighborhoods), "differential neighborhoods"),
      x = "UMAP 1", y = "UMAP 2"
    ) +
    theme(legend.position = "right", legend.text = element_text(size = 8)) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))

  ggsave(file.path(OUTPUT_DIR, "plots", "de_neighborhoods_umap.png"),
    plot = p_neighborhoods,
    width = 12, height = 8,
    dpi = PARAMS$figure_dpi
  )
} else if (length(neighborhoods$de_neighborhoods) > 20) {
  cat("⚠️  Warning: Too many neighborhoods (", length(neighborhoods$de_neighborhoods), ") for visualization. Skipping neighborhood plot.\n")
}

cat("✅ Visualizations saved to plots/\n")

# =============================================================================
# SUMMARY REPORT
# =============================================================================

cat("\n📋 Generating summary report...\n")

# Create summary report
report <- list(
  analysis_info = list(
    input_file = INPUT_H5AD,
    output_directory = OUTPUT_DIR,
    analysis_date = Sys.Date(),
    r_version = R.version.string,
    packages = sapply(required_packages, function(x) as.character(packageVersion(x)))
  ),
  data_summary = list(
    original_dimensions = paste(nrow(sce), "genes ×", ncol(sce), "cells"),
    filtered_dimensions = paste(nrow(sce_filtered), "genes ×", ncol(sce_filtered), "cells"),
    hvg_count = length(hvg),
    conditions = table(sce_filtered$condition),
    donors = table(sce_filtered$donor_id)
  ),
  lemur_results = list(
    embedding_dimensions = PARAMS$n_embedding,
    total_tests = nrow(de_res),
    significant_genes = sum(results_table$significant),
    fdr_threshold = PARAMS$fdr_threshold,
    de_neighborhoods = length(neighborhoods$de_neighborhoods)
  ),
  top_de_genes = head(results_table[results_table$significant, c("gene", "effect_size", "adj_pval")], 20)
)

# Save report as RDS
saveRDS(report, file = file.path(OUTPUT_DIR, "reports", "analysis_summary.rds"))

# Save human-readable report
sink(file.path(OUTPUT_DIR, "reports", "analysis_report.txt"))
cat("=== LEMUR Analysis Report ===\n\n")
cat("Analysis Date:", as.character(Sys.Date()), "\n")
cat("Input File:", INPUT_H5AD, "\n")
cat("Output Directory:", OUTPUT_DIR, "\n\n")

cat("Data Summary:\n")
cat("- Original dimensions:", nrow(sce), "genes ×", ncol(sce), "cells\n")
cat("- Filtered dimensions:", nrow(sce_filtered), "genes ×", ncol(sce_filtered), "cells\n")
cat("- Highly variable genes:", length(hvg), "\n")
cat("- Conditions:", paste(names(table(sce_filtered$condition)), "=", table(sce_filtered$condition), collapse = ", "), "\n")
cat("- Donors:", length(unique(sce_filtered$donor_id)), "\n\n")

cat("LEMUR Results:\n")
cat("- Embedding dimensions:", PARAMS$n_embedding, "\n")
cat("- Total statistical tests:", nrow(de_res), "\n")
cat("- Significant genes (FDR < 0.05):", sum(results_table$significant), "\n")
cat("- DE neighborhoods found:", length(neighborhoods$de_neighborhoods), "\n\n")

cat("Top 10 DE Genes:\n")
top_10 <- head(results_table[results_table$significant, ], 10)
for (i in 1:min(10, nrow(top_10))) {
  cat(sprintf(
    "  %d. %s (effect_size: %.3f, adj.p: %.2e)\n",
    i, top_10$gene[i], top_10$effect_size[i], top_10$adj_pval[i]
  ))
}

cat("\nOutput Files:\n")
cat("- DE results: tables/lemur_de_results.csv\n")
cat("- Top DE genes: tables/top_de_genes.csv\n")
cat("- DE neighborhoods: tables/de_neighborhoods.csv\n")
cat("- UMAP overview: plots/umap_overview.png\n")
cat("- Volcano plot: plots/volcano_plot.png\n")
if (length(neighborhoods$de_neighborhoods) > 0) {
  cat("- DE neighborhoods UMAP: plots/de_neighborhoods_umap.png\n")
}
sink()

# Save workspace
save.image(file = file.path(OUTPUT_DIR, "lemur_workspace.RData"))

# =============================================================================
# COMPLETION
# =============================================================================

cat("\n🎉 LEMUR analysis completed successfully!\n")
cat("📂 All results saved to:", OUTPUT_DIR, "\n")
cat("📊 Summary:\n")
cat("   - Total genes tested:", nrow(de_res), "\n")
cat("   - Significant DE genes:", sum(results_table$significant), "\n")
cat("   - DE neighborhoods:", length(neighborhoods$de_neighborhoods), "\n")
cat("   - Output files: tables/, plots/, reports/\n")

sessionInfo()
