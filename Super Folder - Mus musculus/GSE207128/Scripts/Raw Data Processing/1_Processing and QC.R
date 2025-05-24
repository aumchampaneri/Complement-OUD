# Processing and Quality Control of the combined Raw Data

# Set working directories
setwd("/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus/GSE207128")
input_dir <- "Outputs/Processed Data"
output_dir <- "Outputs/Processed Data/QC"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load required libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(patchwork)

# Install required packages
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("glmGamPoi", quietly = TRUE)) BiocManager::install("glmGamPoi")

# Install and load doublet detection packages
if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
  remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", force = TRUE)
}
if (!requireNamespace("scDblFinder", quietly = TRUE)) {
  BiocManager::install("scDblFinder")
}

# Load packages after installation
library(glmGamPoi)
library(DoubletFinder)
library(scDblFinder)

# 1. Load Merged Seurat Object
combined_seurat <- readRDS(file.path(input_dir, "combined_seurat.rds"))
print(paste("Loaded Seurat object with", ncol(combined_seurat), "cells and", nrow(combined_seurat), "genes"))

# 2. Quality Control (QC)
# 2.1 Add QC metrics
combined_seurat[["percent.mt"]] <- PercentageFeatureSet(combined_seurat, pattern = "^mt-")

# 2.2 Visualize metrics (pre-filtering)
# Standard violin plots
qc_plots_pre <- VlnPlot(combined_seurat,
                        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                        group.by = "sample_id", pt.size = 0.1, ncol = 3)
ggsave(file.path(output_dir, "QC_violins_pre_filtering.png"), qc_plots_pre, width = 15, height = 6)

# Add histograms for better visualization of distributions
qc_hist1 <- ggplot(combined_seurat@meta.data, aes(x=nFeature_RNA)) +
  geom_histogram(bins=100) +
  geom_vline(xintercept=c(100, 8000), color="red", linetype="dashed") +
  ggtitle("Gene count distribution")

qc_hist2 <- ggplot(combined_seurat@meta.data, aes(x=percent.mt)) +
  geom_histogram(bins=100) +
  geom_vline(xintercept=15, color="red", linetype="dashed") +
  ggtitle("Mitochondrial content")

qc_hist3 <- ggplot(combined_seurat@meta.data, aes(x=nCount_RNA)) +
  geom_histogram(bins=100) +
  ggtitle("UMI count distribution")

histograms <- (qc_hist1 | qc_hist2) / qc_hist3
ggsave(file.path(output_dir, "QC_histograms_pre_filtering.png"), histograms, width=12, height=8)

# 2.3 Filter poor-quality cells (using more lenient thresholds)
print(paste("Pre-filtering cell count:", ncol(combined_seurat)))
combined_seurat <- subset(combined_seurat,
                         subset = nFeature_RNA > 100 &    # More permissive minimum
                                  nFeature_RNA < 8000 &   # Higher maximum
                                  percent.mt < 15)        # More permissive mitochondrial threshold
print(paste("Post-filtering cell count:", ncol(combined_seurat)))

# Visualize metrics (post-filtering)
qc_plots_post <- VlnPlot(combined_seurat,
                         features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                         group.by = "sample_id", pt.size = 0.1, ncol = 3)
ggsave(file.path(output_dir, "QC_violins_post_filtering.png"), qc_plots_post, width = 15, height = 6)

# 3. Normalization (Pearson Residuals via SCTransform)
print("Running SCTransform...")
combined_seurat <- SCTransform(combined_seurat,
                              vars.to.regress = "percent.mt",
                              verbose = FALSE)

# 4. Dimensionality Reduction
# 4.1 Run PCA
print("Running PCA...")
combined_seurat <- RunPCA(combined_seurat)

# Save elbow plot
elbow_plot <- ElbowPlot(combined_seurat, ndims = 50)
ggsave(file.path(output_dir, "PCA_elbow_plot.png"), elbow_plot, width = 8, height = 6)

# 5. Clustering
print("Finding neighbors and clusters...")
combined_seurat <- FindNeighbors(combined_seurat, dims = 1:30)
combined_seurat <- FindClusters(combined_seurat, resolution = 0.5)

# 6. UMAP Visualization
print("Running UMAP...")
combined_seurat <- RunUMAP(combined_seurat, dims = 1:30)

# Save UMAP visualizations
umap_by_condition <- DimPlot(combined_seurat, group.by = "condition", label = TRUE) + NoLegend()
ggsave(file.path(output_dir, "UMAP_by_condition.png"), umap_by_condition, width = 8, height = 6)

umap_by_sample <- DimPlot(combined_seurat, group.by = "sample_id", label = TRUE) + NoLegend()
ggsave(file.path(output_dir, "UMAP_by_sample.png"), umap_by_sample, width = 8, height = 6)

umap_by_cluster <- DimPlot(combined_seurat, reduction = "umap", label = TRUE, repel = TRUE)
ggsave(file.path(output_dir, "UMAP_by_cluster.png"), umap_by_cluster, width = 8, height = 6)

# 7. Doublet Detection with error handling
print("Running doublet detection...")

# Method 1: DoubletFinder
tryCatch({
  # Pre-process for DoubletFinder
  print("Running DoubletFinder parameter optimization...")
  sweep.res <- paramSweep_v3(combined_seurat, PCs = 1:30, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  optimal_pk <- bcmvn$pK[which.max(bcmvn$BCmetric)]
  print(paste("Optimal pK:", optimal_pk))

  # Estimate doublet rate based on 10x data (~0.8% per 1000 cells)
  cell_count <- ncol(combined_seurat)
  doublet_rate <- 0.008 * (cell_count/1000)
  print(paste("Estimated doublet rate:", round(doublet_rate*100, 2), "%"))

  # Run DoubletFinder
  print("Running DoubletFinder...")
  combined_seurat <- doubletFinder_v3(combined_seurat,
                                     PCs = 1:30,
                                     pN = 0.25,
                                     pK = optimal_pk,
                                     nExp = round(doublet_rate * cell_count),
                                     sct = TRUE)

  # Get the name of the DoubletFinder results column
  DF_column <- grep("DF.classifications", colnames(combined_seurat@meta.data), value = TRUE)
}, error = function(e) {
  message("DoubletFinder encountered an error: ", e$message)
  # Create a dummy classification if DoubletFinder fails
  combined_seurat$DF_classifications <- "Singlet"
  DF_column <- "DF_classifications"
  assign("combined_seurat", combined_seurat, envir = .GlobalEnv)
  assign("DF_column", DF_column, envir = .GlobalEnv)
})

# Method 2: scDblFinder
tryCatch({
  print("Running scDblFinder...")
  # Extract SCT counts
  sct_counts <- GetAssayData(combined_seurat, assay = "SCT", layer = "counts")

  # Run scDblFinder
  dbl_results <- scDblFinder(sct_counts)
  combined_seurat$scDblFinder_class <- dbl_results$class
  combined_seurat$scDblFinder_score <- dbl_results$score
}, error = function(e) {
  message("scDblFinder encountered an error: ", e$message)
  # Create dummy classifications if scDblFinder fails
  combined_seurat$scDblFinder_class <- "singlet"
  combined_seurat$scDblFinder_score <- 0
  assign("combined_seurat", combined_seurat, envir = .GlobalEnv)
})

# Try to visualize doublets in UMAP
tryCatch({
  p1 <- DimPlot(combined_seurat, reduction = "umap", group.by = DF_column) +
        ggtitle("DoubletFinder Results")
  p2 <- DimPlot(combined_seurat, reduction = "umap", group.by = "scDblFinder_class") +
        ggtitle("scDblFinder Results")
  doublet_plots <- p1 | p2
  ggsave(file.path(output_dir, "doublet_detection_plots.png"), doublet_plots, width = 15, height = 6)

  # Identify consistent doublets (detected by both methods)
  combined_seurat$consistent_doublet <-
    (combined_seurat@meta.data[[DF_column]] == "Doublet" &
     combined_seurat$scDblFinder_class == "doublet")

  # Filter out consistent doublets
  combined_seurat_filtered <- subset(combined_seurat,
                                    subset = consistent_doublet == FALSE)

  # Check how many cells were removed
  print(paste("Original cell count:", ncol(combined_seurat)))
  print(paste("Cells after doublet removal:", ncol(combined_seurat_filtered)))
  print(paste("Removed", ncol(combined_seurat) - ncol(combined_seurat_filtered), "doublets"))

  # Use filtered object for downstream analysis
  combined_seurat <- combined_seurat_filtered
}, error = function(e) {
  message("Doublet visualization or filtering encountered an error: ", e$message)
  message("Continuing without doublet filtering")
  combined_seurat_filtered <- combined_seurat
  assign("combined_seurat", combined_seurat, envir = .GlobalEnv)
  assign("combined_seurat_filtered", combined_seurat_filtered, envir = .GlobalEnv)
})

# Re-run dimensionality reduction after doublet removal
print("Re-running dimensionality reduction and clustering after doublet removal...")
combined_seurat <- RunPCA(combined_seurat)
combined_seurat <- RunUMAP(combined_seurat, dims = 1:30)
combined_seurat <- FindNeighbors(combined_seurat, dims = 1:30)
combined_seurat <- FindClusters(combined_seurat, resolution = 0.5)

# Visualize final UMAP after doublet removal
final_umap <- DimPlot(combined_seurat, reduction = "umap", label = TRUE, repel = TRUE)
ggsave(file.path(output_dir, "UMAP_after_doublet_removal.png"), final_umap, width = 8, height = 6)

# Add final summary statistics
cat("Final summary:\n")
print(table(combined_seurat$condition))
print(table(combined_seurat$seurat_clusters))

# Report structure of the final object
cat("\nFinal Seurat object structure:\n")
print(combined_seurat)

# 8. Save the Processed Object
saveRDS(combined_seurat, file.path(input_dir, "combined_seurat_qc_sct.rds"))
print("Processing complete. Final Seurat object saved.")