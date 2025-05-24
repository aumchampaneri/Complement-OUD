# =============================================================================
# Enhanced Single-Cell RNA-seq Processing and Quality Control Pipeline
# =============================================================================
# 
# Purpose: Comprehensive processing optimized for large datasets
#
# Date: May 24, 2025
# Version: 4.0 (Streamlined - No Doublet Detection)
#
# Key Features:
# - Skip memory-intensive doublet detection
# - Streamlined processing for large datasets
# - All essential QC and analysis steps preserved
# =============================================================================

# Set seed for reproducibility
set.seed(42)

# =============================================================================
# 1. SETUP AND CONFIGURATION
# =============================================================================

# Set working directories
setwd("/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus/GSE207128")
input_dir <- "Outputs/Processed Data"
output_dir <- "Outputs/Processed Data/QC"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Create subdirectories for organized output
dir.create(file.path(output_dir, "plots"), showWarnings = FALSE)
dir.create(file.path(output_dir, "statistics"), showWarnings = FALSE)
dir.create(file.path(output_dir, "thresholds"), showWarnings = FALSE)

# Define analysis parameters
QC_PARAMS <- list(
  # Basic QC thresholds
  min_features_base = 200,
  max_features_base = 8000,
  min_counts_base = 500,
  max_counts_base = 40000,
  max_mt_percent_base = 15,
  
  # Analysis parameters
  n_variable_features = 2000,
  pca_dims = 1:30,
  clustering_resolutions = seq(0.1, 1.2, by = 0.1),
  default_resolution = 0.5,
  mad_threshold = 3  # MAD multiplier for adaptive thresholds
)

# Package installation function
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing package:", pkg))
    if(pkg == "clustree") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install("clustree")
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
  }
}

# Core packages
core_packages <- c("Seurat", "dplyr", "Matrix", "ggplot2", "patchwork", 
                   "glmGamPoi", "cluster", "jsonlite")

# Optional packages
optional_packages <- c("clustree", "gridExtra")

# Install and load core packages
for(pkg in core_packages) {
  install_if_missing(pkg)
  library(pkg, character.only = TRUE)
}

# Try to load optional packages
clustree_available <- FALSE
for(pkg in optional_packages) {
  tryCatch({
    install_if_missing(pkg)
    library(pkg, character.only = TRUE)
    if(pkg == "clustree") clustree_available <- TRUE
  }, error = function(e) {
    message(paste("Optional package", pkg, "not available:", e$message))
  })
}

# Save session information
session_info <- capture.output(sessionInfo())
writeLines(session_info, file.path(output_dir, "session_info.txt"))

# Save analysis parameters
saveRDS(QC_PARAMS, file.path(output_dir, "analysis_parameters.rds"))
writeLines(jsonlite::toJSON(QC_PARAMS, pretty = TRUE), 
           file.path(output_dir, "analysis_parameters.json"))

print("=== Enhanced scRNA-seq Processing Pipeline (Streamlined) ===")
print(paste("Analysis started at:", Sys.time()))
print(paste("Random seed set to:", 42))

# =============================================================================
# 2. DATA LOADING AND INITIAL ASSESSMENT
# =============================================================================

print("Loading merged Seurat object...")
combined_seurat <- readRDS(file.path(input_dir, "combined_seurat.rds"))

# Initial data summary
initial_stats <- list(
  total_cells = ncol(combined_seurat),
  total_genes = nrow(combined_seurat),
  samples = unique(combined_seurat$sample_id),
  conditions = unique(combined_seurat$condition)
)

print(paste("Loaded Seurat object with", initial_stats$total_cells, "cells and", 
            initial_stats$total_genes, "genes"))
print(paste("Samples:", paste(initial_stats$samples, collapse = ", ")))
print(paste("Conditions:", paste(initial_stats$conditions, collapse = ", ")))

# =============================================================================
# 3. QC METRICS CALCULATION
# =============================================================================

print("Calculating QC metrics...")

# Standard QC metrics
combined_seurat[["percent.mt"]] <- PercentageFeatureSet(combined_seurat, pattern = "^mt-")
combined_seurat[["percent.ribo"]] <- PercentageFeatureSet(combined_seurat, pattern = "^Rp[sl]")
combined_seurat[["log10GenesPerUMI"]] <- log10(combined_seurat$nFeature_RNA) / log10(combined_seurat$nCount_RNA)
combined_seurat[["gene_umi_ratio"]] <- combined_seurat$nFeature_RNA / combined_seurat$nCount_RNA

# =============================================================================
# 4. ADAPTIVE QC THRESHOLD CALCULATION
# =============================================================================

print("Calculating adaptive QC thresholds...")

# Function to calculate sample-specific thresholds using MAD
calculate_adaptive_thresholds <- function(seurat_obj, metric, n_mad = 3, min_val = NULL, max_val = NULL) {
  per_sample_stats <- seurat_obj@meta.data %>%
    group_by(sample_id) %>%
    summarise(
      n_cells = n(),
      median_val = median(!!sym(metric), na.rm = TRUE),
      mad_val = mad(!!sym(metric), na.rm = TRUE),
      q25 = quantile(!!sym(metric), 0.25, na.rm = TRUE),
      q75 = quantile(!!sym(metric), 0.75, na.rm = TRUE),
      iqr = q75 - q25,
      lower_thresh_mad = median_val - n_mad * mad_val,
      upper_thresh_mad = median_val + n_mad * mad_val,
      lower_thresh_iqr = q25 - 1.5 * iqr,
      upper_thresh_iqr = q75 + 1.5 * iqr,
      .groups = 'drop'
    )
  
  # Apply global constraints if provided
  if(!is.null(min_val)) {
    per_sample_stats$lower_thresh_mad <- pmax(per_sample_stats$lower_thresh_mad, min_val)
    per_sample_stats$lower_thresh_iqr <- pmax(per_sample_stats$lower_thresh_iqr, min_val)
  }
  if(!is.null(max_val)) {
    per_sample_stats$upper_thresh_mad <- pmin(per_sample_stats$upper_thresh_mad, max_val)
    per_sample_stats$upper_thresh_iqr <- pmin(per_sample_stats$upper_thresh_iqr, max_val)
  }
  
  return(per_sample_stats)
}

# Calculate adaptive thresholds for key metrics
feature_thresholds <- calculate_adaptive_thresholds(combined_seurat, "nFeature_RNA", 
                                                   n_mad = QC_PARAMS$mad_threshold,
                                                   min_val = QC_PARAMS$min_features_base)

count_thresholds <- calculate_adaptive_thresholds(combined_seurat, "nCount_RNA",
                                                 n_mad = QC_PARAMS$mad_threshold,
                                                 min_val = QC_PARAMS$min_counts_base)

mt_thresholds <- calculate_adaptive_thresholds(combined_seurat, "percent.mt",
                                              n_mad = QC_PARAMS$mad_threshold,
                                              max_val = QC_PARAMS$max_mt_percent_base)

# Save threshold information
write.csv(feature_thresholds, file.path(output_dir, "thresholds", "adaptive_feature_thresholds.csv"), row.names = FALSE)
write.csv(count_thresholds, file.path(output_dir, "thresholds", "adaptive_count_thresholds.csv"), row.names = FALSE)
write.csv(mt_thresholds, file.path(output_dir, "thresholds", "adaptive_mt_thresholds.csv"), row.names = FALSE)

print("Adaptive thresholds calculated and saved.")

# =============================================================================
# 5. PRE-FILTERING VISUALIZATION
# =============================================================================

print("Creating pre-filtering QC visualizations...")

# Violin plots
qc_plots_pre <- VlnPlot(combined_seurat,
                        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
                        group.by = "sample_id", pt.size = 0.05, ncol = 2, raster = TRUE)
ggsave(file.path(output_dir, "plots", "QC_violins_pre_filtering.png"), 
       qc_plots_pre, width = 16, height = 10, dpi = 300)

# Enhanced histograms with thresholds
create_qc_histograms <- function(seurat_obj, thresholds_list) {
  
  # Feature histogram
  p1 <- ggplot(seurat_obj@meta.data, aes(x = nFeature_RNA)) +
    geom_histogram(bins = 100, alpha = 0.7, fill = "skyblue") +
    geom_vline(xintercept = c(QC_PARAMS$min_features_base, QC_PARAMS$max_features_base), 
               color = "red", linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = median(thresholds_list$feature$lower_thresh_mad, na.rm = TRUE),
               color = "blue", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = median(thresholds_list$feature$upper_thresh_mad, na.rm = TRUE),
               color = "blue", linetype = "dotted", linewidth = 1) +
    labs(title = "Gene Count Distribution",
         subtitle = "Red: Fixed thresholds, Blue: Adaptive thresholds",
         x = "Number of Features", y = "Frequency") +
    theme_minimal()
  
  # MT content histogram
  p2 <- ggplot(seurat_obj@meta.data, aes(x = percent.mt)) +
    geom_histogram(bins = 100, alpha = 0.7, fill = "lightcoral") +
    geom_vline(xintercept = QC_PARAMS$max_mt_percent_base, 
               color = "red", linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = median(thresholds_list$mt$upper_thresh_mad, na.rm = TRUE),
               color = "blue", linetype = "dotted", linewidth = 1) +
    labs(title = "Mitochondrial Content",
         subtitle = "Red: Fixed threshold, Blue: Adaptive threshold",
         x = "Mitochondrial %", y = "Frequency") +
    theme_minimal()
  
  # UMI count histogram
  p3 <- ggplot(seurat_obj@meta.data, aes(x = nCount_RNA)) +
    geom_histogram(bins = 100, alpha = 0.7, fill = "lightgreen") +
    geom_vline(xintercept = c(QC_PARAMS$min_counts_base, QC_PARAMS$max_counts_base),
               color = "red", linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = median(thresholds_list$count$lower_thresh_mad, na.rm = TRUE),
               color = "blue", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = median(thresholds_list$count$upper_thresh_mad, na.rm = TRUE),
               color = "blue", linetype = "dotted", linewidth = 1) +
    labs(title = "UMI Count Distribution",
         subtitle = "Red: Fixed thresholds, Blue: Adaptive thresholds",
         x = "UMI Count", y = "Frequency") +
    theme_minimal()
  
  # Complexity histogram
  p4 <- ggplot(seurat_obj@meta.data, aes(x = log10GenesPerUMI)) +
    geom_histogram(bins = 100, alpha = 0.7, fill = "plum") +
    labs(title = "Transcriptional Complexity",
         x = "log10(Genes/UMI)", y = "Frequency") +
    theme_minimal()
  
  return(list(p1, p2, p3, p4))
}

threshold_list <- list(feature = feature_thresholds, count = count_thresholds, mt = mt_thresholds)
hist_plots <- create_qc_histograms(combined_seurat, threshold_list)

# Combine histograms
histograms_combined <- (hist_plots[[1]] | hist_plots[[2]]) / (hist_plots[[3]] | hist_plots[[4]])
ggsave(file.path(output_dir, "plots", "QC_histograms_pre_filtering.png"), 
       histograms_combined, width = 16, height = 12, dpi = 300)

# =============================================================================
# 6. CELL FILTERING
# =============================================================================

print("Applying cell filtering...")

# Store pre-filtering counts
pre_filter_count <- ncol(combined_seurat)
print(paste("Pre-filtering cell count:", pre_filter_count))

# Store original for comparison
original_seurat <- combined_seurat

# Apply filtering
combined_seurat_filtered <- subset(combined_seurat,
                                  subset = nFeature_RNA >= QC_PARAMS$min_features_base &
                                          nFeature_RNA <= QC_PARAMS$max_features_base &
                                          nCount_RNA >= QC_PARAMS$min_counts_base &
                                          nCount_RNA <= QC_PARAMS$max_counts_base &
                                          percent.mt <= QC_PARAMS$max_mt_percent_base)

post_filter_count <- ncol(combined_seurat_filtered)
cells_removed <- pre_filter_count - post_filter_count
removal_percentage <- round((cells_removed / pre_filter_count) * 100, 2)

print(paste("Post-filtering cell count:", post_filter_count))
print(paste("Cells removed:", cells_removed, "(", removal_percentage, "%)"))

# Update object
combined_seurat <- combined_seurat_filtered

# Create filtering summary
sample_filtering_summary <- combined_seurat@meta.data %>%
  dplyr::count(sample_id, name = "cells_retained") %>%
  left_join(
    original_seurat@meta.data %>% 
      dplyr::count(sample_id, name = "total_cells"),
    by = "sample_id"
  ) %>%
  mutate(
    cells_removed = total_cells - cells_retained,
    removal_percentage = round((cells_removed / total_cells) * 100, 2)
  )

write.csv(sample_filtering_summary, 
          file.path(output_dir, "statistics", "filtering_summary_by_sample.csv"), 
          row.names = FALSE)

# =============================================================================
# 7. POST-FILTERING VISUALIZATION
# =============================================================================

print("Creating post-filtering QC visualizations...")

# Post-filtering violin plots
qc_plots_post <- VlnPlot(combined_seurat,
                         features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
                         group.by = "sample_id", pt.size = 0.05, ncol = 2, raster = TRUE)
ggsave(file.path(output_dir, "plots", "QC_violins_post_filtering.png"), 
       qc_plots_post, width = 16, height = 10, dpi = 300)

# Comparison plots (pre vs post filtering)
original_seurat_metrics <- original_seurat@meta.data
original_seurat_metrics$percent.mt <- PercentageFeatureSet(original_seurat, pattern = "^mt-")

comparison_data <- rbind(
  data.frame(
    nFeature_RNA = original_seurat_metrics$nFeature_RNA,
    nCount_RNA = original_seurat_metrics$nCount_RNA,
    percent.mt = original_seurat_metrics$percent.mt,
    stage = "Pre-filtering"
  ),
  data.frame(
    nFeature_RNA = combined_seurat@meta.data$nFeature_RNA,
    nCount_RNA = combined_seurat@meta.data$nCount_RNA,
    percent.mt = combined_seurat@meta.data$percent.mt,
    stage = "Post-filtering"
  )
)

comparison_plots <- list()
for(metric in c("nFeature_RNA", "nCount_RNA", "percent.mt")) {
  comparison_plots[[metric]] <- ggplot(comparison_data, aes(x = .data[[metric]], fill = stage)) +
    geom_density(alpha = 0.6) +
    labs(title = paste("Distribution of", metric), x = metric, y = "Density") +
    theme_minimal() +
    scale_fill_manual(values = c("Pre-filtering" = "lightcoral", "Post-filtering" = "lightblue"))
}

comparison_combined <- wrap_plots(comparison_plots, ncol = 1)
ggsave(file.path(output_dir, "plots", "filtering_comparison.png"), 
       comparison_combined, width = 10, height = 12, dpi = 300)

# =============================================================================
# 8. NORMALIZATION AND FEATURE SELECTION
# =============================================================================

print("Performing normalization and feature selection...")

# Standard log-normalization
combined_seurat <- NormalizeData(combined_seurat, 
                                normalization.method = "LogNormalize", 
                                scale.factor = 10000,
                                verbose = FALSE)

# Find variable features
combined_seurat <- FindVariableFeatures(combined_seurat, 
                                       selection.method = "vst", 
                                       nfeatures = QC_PARAMS$n_variable_features,
                                       verbose = FALSE)

# Visualize variable features
variable_features_plot <- VariableFeaturePlot(combined_seurat) +
  labs(title = paste("Top", QC_PARAMS$n_variable_features, "Variable Features"))
ggsave(file.path(output_dir, "plots", "variable_features.png"), 
       variable_features_plot, width = 12, height = 8, dpi = 300)

# Scale data (only variable features for memory efficiency)
print("Scaling data...")
combined_seurat <- ScaleData(combined_seurat, 
                            features = VariableFeatures(combined_seurat),
                            vars.to.regress = c("percent.mt"),
                            verbose = FALSE)

# =============================================================================
# 9. DIMENSIONALITY REDUCTION
# =============================================================================

print("Performing PCA...")
combined_seurat <- RunPCA(combined_seurat, 
                         features = VariableFeatures(object = combined_seurat),
                         verbose = FALSE)

# PCA visualization
pca_plots <- list(
  ElbowPlot(combined_seurat, ndims = 50) + 
    labs(title = "PCA Elbow Plot") + 
    geom_vline(xintercept = max(QC_PARAMS$pca_dims), color = "red", linetype = "dashed"),
  
  DimPlot(combined_seurat, reduction = "pca", group.by = "sample_id", raster = TRUE) + 
    labs(title = "PCA by Sample"),
  
  DimPlot(combined_seurat, reduction = "pca", group.by = "condition", raster = TRUE) + 
    labs(title = "PCA by Condition")
)

pca_combined <- wrap_plots(pca_plots, ncol = 2)
ggsave(file.path(output_dir, "plots", "PCA_analysis.png"), 
       pca_combined, width = 16, height = 10, dpi = 300)

# =============================================================================
# 10. CLUSTERING OPTIMIZATION
# =============================================================================

print("Optimizing clustering resolution...")

# Test multiple clustering resolutions
combined_seurat <- FindNeighbors(combined_seurat, dims = QC_PARAMS$pca_dims, verbose = FALSE)

# Test resolutions
for(res in QC_PARAMS$clustering_resolutions) {
  combined_seurat <- FindClusters(combined_seurat, resolution = res, verbose = FALSE)
}

# Create clustree plot if available
if(clustree_available) {
  tryCatch({
    clustree_plot <- clustree(combined_seurat, prefix = "RNA_snn_res.")
    ggsave(file.path(output_dir, "plots", "clustering_resolution_tree.png"), 
           clustree_plot, width = 12, height = 10, dpi = 300)
  }, error = function(e) {
    message("Could not create clustree plot: ", e$message)
  })
}

# Calculate silhouette scores (using sample for large datasets)
pca_embeddings <- Embeddings(combined_seurat, "pca")[, QC_PARAMS$pca_dims]

silhouette_scores <- sapply(QC_PARAMS$clustering_resolutions, function(res) {
  clusters <- combined_seurat@meta.data[[paste0("RNA_snn_res.", res)]]
  if(length(unique(clusters)) > 1 && length(unique(clusters)) < nrow(pca_embeddings)) {
    tryCatch({
      # Use sample for large datasets
      sample_size <- min(5000, nrow(pca_embeddings))
      sample_idx <- sample(nrow(pca_embeddings), sample_size)
      dist_matrix <- dist(pca_embeddings[sample_idx, ])
      clusters_sample <- clusters[sample_idx]
      
      sil <- cluster::silhouette(as.numeric(as.factor(clusters_sample)), dist_matrix)
      return(mean(sil[,3]))
    }, error = function(e) return(NA))
  } else {
    return(NA)
  }
})

# Save resolution analysis
res_analysis <- data.frame(
  resolution = QC_PARAMS$clustering_resolutions, 
  silhouette_score = silhouette_scores,
  n_clusters = sapply(QC_PARAMS$clustering_resolutions, function(res) {
    length(unique(combined_seurat@meta.data[[paste0("RNA_snn_res.", res)]]))
  })
)

write.csv(res_analysis, file.path(output_dir, "statistics", "resolution_analysis.csv"), row.names = FALSE)

# Plot resolution analysis
res_plot <- ggplot(res_analysis, aes(x = resolution)) +
  geom_line(aes(y = silhouette_score), color = "blue", linewidth = 1, na.rm = TRUE) +
  geom_point(aes(y = silhouette_score), color = "blue", size = 2, na.rm = TRUE) +
  geom_line(aes(y = n_clusters/max(n_clusters, na.rm = TRUE)), color = "red", linewidth = 1) +
  geom_point(aes(y = n_clusters/max(n_clusters, na.rm = TRUE)), color = "red", size = 2) +
  scale_y_continuous(
    name = "Silhouette Score",
    sec.axis = sec_axis(~.*max(res_analysis$n_clusters, na.rm = TRUE), name = "Number of Clusters")
  ) +
  labs(title = "Clustering Resolution Optimization",
       x = "Resolution",
       subtitle = "Blue: Silhouette Score, Red: Number of Clusters") +
  theme_minimal()

ggsave(file.path(output_dir, "plots", "resolution_optimization.png"), 
       res_plot, width = 10, height = 6, dpi = 300)

# Select optimal resolution
valid_scores <- !is.na(res_analysis$silhouette_score)
if(any(valid_scores)) {
  optimal_resolution <- res_analysis$resolution[which.max(res_analysis$silhouette_score[valid_scores])]
} else {
  optimal_resolution <- QC_PARAMS$default_resolution
}

print(paste("Optimal resolution selected:", optimal_resolution))

# Set optimal clustering
combined_seurat$seurat_clusters <- combined_seurat@meta.data[[paste0("RNA_snn_res.", optimal_resolution)]]
Idents(combined_seurat) <- combined_seurat$seurat_clusters

# =============================================================================
# 11. UMAP VISUALIZATION
# =============================================================================

print("Generating UMAP embedding...")
combined_seurat <- RunUMAP(combined_seurat, dims = QC_PARAMS$pca_dims, verbose = FALSE)

# Create UMAP visualizations
umap_plots <- list(
  DimPlot(combined_seurat, group.by = "condition", label = TRUE, repel = TRUE, raster = TRUE) + 
    labs(title = "UMAP by Condition") + NoLegend(),
  
  DimPlot(combined_seurat, group.by = "sample_id", label = TRUE, repel = TRUE, raster = TRUE) + 
    labs(title = "UMAP by Sample") + NoLegend(),
  
  DimPlot(combined_seurat, reduction = "umap", label = TRUE, repel = TRUE, raster = TRUE) + 
    labs(title = paste("UMAP by Cluster (res =", optimal_resolution, ")")),
  
  FeaturePlot(combined_seurat, features = "nFeature_RNA", raster = TRUE) + 
    labs(title = "Gene Count per Cell"),
  
  FeaturePlot(combined_seurat, features = "nCount_RNA", raster = TRUE) + 
    labs(title = "UMI Count per Cell"),
  
  FeaturePlot(combined_seurat, features = "percent.mt", raster = TRUE) + 
    labs(title = "Mitochondrial %")
)

# Save individual UMAP plots
plot_names <- c("condition", "sample", "cluster", "gene_count", "umi_count", "mt_percent")
for(i in 1:length(umap_plots)) {
  ggsave(file.path(output_dir, "plots", paste0("UMAP_by_", plot_names[i], ".png")), 
         umap_plots[[i]], width = 10, height = 8, dpi = 300)
}

# Combined UMAP plot
umap_combined <- wrap_plots(umap_plots[1:3], ncol = 2)
ggsave(file.path(output_dir, "plots", "UMAP_combined_overview.png"), 
       umap_combined, width = 16, height = 12, dpi = 300)

# =============================================================================
# 12. CELL CYCLE SCORING
# =============================================================================

print("Running cell cycle scoring with mouse-specific genes...")

tryCatch({
  # Mouse cell cycle genes
  mouse_s_genes <- c("Mcm5", "Pcna", "Tyms", "Fen1", "Mcm7", "Mcm4", "Rrm1", "Ung", 
                     "Gins2", "Mcm6", "Cdca7", "Dtl", "Prim1", "Uhrf1", "Cenpu", 
                     "Hells", "Rfc2", "Rpa2", "Nasp", "Rad51ap1", "Gmnn", "Wdr76", 
                     "Slbp", "Ccne2", "Ubr7", "Pold3", "Msh2", "Atad2", "Rad51", 
                     "Rrm2", "Cdc45", "Cdc6", "Exo1", "Tipin", "Dscc1", "Blm", 
                     "Casp8ap2", "Usp1", "Clspn", "Pola1", "Chaf1b", "Mrpl36", "E2f8")
  
  mouse_g2m_genes <- c("Hmgb2", "Cdk1", "Nusap1", "Ube2c", "Birc5", "Tpx2", "Top2a", 
                       "Ndc80", "Cks2", "Nuf2", "Cks1b", "Mki67", "Tmpo", "Cenpf", 
                       "Tacc3", "Fam64a", "Smc4", "Ccnb2", "Ckap2l", "Ckap2", "Aurkb", 
                       "Bub1", "Kif11", "Anp32e", "Tubb4b", "Gtse1", "Kif20b", "Hjurp", 
                       "Cdca3", "Hn1", "Cdc20", "Ttk", "Cdc25c", "Kif2c", "Rangap1", 
                       "Ncapd2", "Dlgap5", "Cdca2", "Cdca8", "Ect2", "Kif23", "Hmmr", 
                       "Aurka", "Psrc1", "Anln", "Lbr", "Ckap5", "Cenpe", "Ctcf", 
                       "Nek2", "G2e3", "Gas2l3", "Cbx5", "Cenpa")
  
  # Filter available genes
  available_s_genes <- mouse_s_genes[mouse_s_genes %in% rownames(combined_seurat)]
  available_g2m_genes <- mouse_g2m_genes[mouse_g2m_genes %in% rownames(combined_seurat)]
  
  print(paste("Using", length(available_s_genes), "S phase genes and", 
              length(available_g2m_genes), "G2M phase genes"))
  
  if(length(available_s_genes) >= 5 && length(available_g2m_genes) >= 5) {
    combined_seurat <- CellCycleScoring(combined_seurat, 
                                       s.features = available_s_genes,
                                       g2m.features = available_g2m_genes,
                                       set.ident = FALSE)
    
    # Visualize cell cycle effects
    cc_plot <- DimPlot(combined_seurat, group.by = "Phase", raster = TRUE) +
      labs(title = "Cell Cycle Phase Distribution")
    ggsave(file.path(output_dir, "plots", "cell_cycle_phases.png"), 
           cc_plot, width = 10, height = 8, dpi = 300)
    
    # Cell cycle summary
    cc_summary <- table(combined_seurat$Phase, combined_seurat$condition)
    write.csv(as.data.frame.matrix(cc_summary), 
              file.path(output_dir, "statistics", "cell_cycle_summary.csv"))
    
    print("Cell cycle scoring completed successfully")
    
  } else {
    print("Insufficient cell cycle genes found in dataset")
    combined_seurat$Phase <- "Unknown"
    combined_seurat$S.Score <- 0
    combined_seurat$G2M.Score <- 0
  }
  
}, error = function(e) {
  message("Cell cycle scoring failed: ", e$message)
  combined_seurat$Phase <- "Unknown"
  combined_seurat$S.Score <- 0
  combined_seurat$G2M.Score <- 0
})

# =============================================================================
# 13. MARKER GENE ANALYSIS
# =============================================================================

print("Creating marker gene visualizations...")

# Known mouse brain markers
known_markers <- c("Cd68", "Cx3cr1", "Tmem119", "Gfap", "Mbp", "Snap25", "Slc17a7", "Gad1")
available_markers <- known_markers[known_markers %in% rownames(combined_seurat)]

if(length(available_markers) > 0) {
  print(paste("Creating plots for:", paste(available_markers, collapse = ", ")))
  
  # Create individual marker plots for memory efficiency
  for(marker in available_markers) {
    marker_plot <- FeaturePlot(combined_seurat, features = marker, raster = TRUE) +
      labs(title = paste("Expression of", marker))
    
    ggsave(file.path(output_dir, "plots", paste0("marker_", marker, ".png")), 
           marker_plot, width = 8, height = 6, dpi = 300)
  }
  
  print("Marker expression plots created")
} else {
  print("No known marker genes found in dataset")
}

# =============================================================================
# 14. FINAL PROCESSING AND VALIDATION
# =============================================================================

print("Performing final processing steps...")

# Final UMAP
final_umap <- DimPlot(combined_seurat, reduction = "umap", label = TRUE, repel = TRUE, raster = TRUE) +
  labs(title = "Final UMAP - Complete Processing")
ggsave(file.path(output_dir, "plots", "UMAP_final.png"), 
       final_umap, width = 10, height = 8, dpi = 300)

# =============================================================================
# 15. FINAL METRICS AND SAVE
# =============================================================================

print("Generating final metrics...")

# Calculate final metrics
final_metrics <- list(
  processing_date = Sys.Date(),
  processing_time = Sys.time(),
  initial_cells = initial_stats$total_cells,
  final_cells = ncol(combined_seurat),
  cells_removed = initial_stats$total_cells - ncol(combined_seurat),
  removal_percentage = round((initial_stats$total_cells - ncol(combined_seurat)) / initial_stats$total_cells * 100, 2),
  final_genes = nrow(combined_seurat),
  final_clusters = length(unique(combined_seurat$seurat_clusters)),
  optimal_resolution = optimal_resolution,
  median_genes_per_cell = median(combined_seurat$nFeature_RNA),
  median_umis_per_cell = median(combined_seurat$nCount_RNA),
  median_mt_percent = median(combined_seurat$percent.mt),
  cells_per_condition = as.list(table(combined_seurat$condition)),
  cells_per_sample = as.list(table(combined_seurat$sample_id)),
  doublet_detection_performed = FALSE,
  streamlined_processing = TRUE,
  available_markers = available_markers
)

# Save final metrics
saveRDS(final_metrics, file.path(output_dir, "statistics", "final_processing_metrics.rds"))

# Print summary
cat("\n=== FINAL PROCESSING SUMMARY ===\n")
cat("Initial cells:", final_metrics$initial_cells, "\n")
cat("Final cells:", final_metrics$final_cells, "\n")
cat("Cells removed (QC filtering only):", final_metrics$cells_removed, "(", final_metrics$removal_percentage, "%)\n")
cat("Final clusters:", final_metrics$final_clusters, "\n")
cat("Optimal resolution:", final_metrics$optimal_resolution, "\n")
cat("Median genes per cell:", final_metrics$median_genes_per_cell, "\n")
cat("Median UMIs per cell:", final_metrics$median_umis_per_cell, "\n")

print("Final cell distribution by condition:")
print(table(combined_seurat$condition))

print("Final cell distribution by sample:")
print(table(combined_seurat$sample_id))

print("Final cluster distribution:")
print(table(combined_seurat$seurat_clusters))

# =============================================================================
# 16. SAVE PROCESSED DATA
# =============================================================================

print("Saving processed data...")

# Save the final processed Seurat object
saveRDS(combined_seurat, file.path(input_dir, "processed_seurat.rds"))

# Create processing report
create_processing_report <- function(output_dir, final_metrics, QC_PARAMS) {
  report_lines <- c(
    "# Single-Cell RNA-seq Processing Report - STREAMLINED",
    paste("Generated on:", Sys.Date()),
    "",
    "## Processing Parameters",
    paste("- Random seed:", 42),
    paste("- Minimum features:", QC_PARAMS$min_features_base),
    paste("- Maximum features:", QC_PARAMS$max_features_base),
    paste("- Minimum counts:", QC_PARAMS$min_counts_base),
    paste("- Maximum counts:", QC_PARAMS$max_counts_base),
    paste("- Maximum mitochondrial %:", QC_PARAMS$max_mt_percent_base),
    paste("- Variable features:", QC_PARAMS$n_variable_features),
    paste("- PCA dimensions:", paste(range(QC_PARAMS$pca_dims), collapse = "-")),
    paste("- Optimal clustering resolution:", final_metrics$optimal_resolution),
    "",
    "## Processing Features",
    "- ✅ Comprehensive QC analysis",
    "- ✅ Adaptive threshold calculation",
    "- ✅ Variable feature scaling only",
    "- ✅ Resolution optimization",
    "- ✅ Cell cycle scoring",
    "- ❌ Doublet detection (skipped for performance)",
    "",
    "## Results Summary",
    paste("- Initial cells:", final_metrics$initial_cells),
    paste("- Final cells:", final_metrics$final_cells),
    paste("- Cells removed (QC only):", final_metrics$cells_removed, "(", final_metrics$removal_percentage, "%)"),
    paste("- Final clusters:", final_metrics$final_clusters),
    paste("- Median genes per cell:", final_metrics$median_genes_per_cell),
    paste("- Median UMIs per cell:", final_metrics$median_umis_per_cell),
    "",
    "## Output Files",
    "### Main Results",
    "- processed_seurat.rds - Final processed Seurat object",
    "",
    "### Plots",
    "- QC_violins_pre_filtering.png - Pre-filtering QC metrics",
    "- QC_violins_post_filtering.png - Post-filtering QC metrics", 
    "- QC_histograms_pre_filtering.png - QC histograms with thresholds",
    "- filtering_comparison.png - Pre vs post filtering comparison",
    "- variable_features.png - Highly variable features",
    "- PCA_analysis.png - PCA visualization",
    "- clustering_resolution_tree.png - Resolution optimization tree",
    "- resolution_optimization.png - Resolution analysis plot",
    "- UMAP_combined_overview.png - UMAP overview",
    "- cell_cycle_phases.png - Cell cycle analysis",
    "- marker_*.png - Individual marker expression plots",
    "- UMAP_final.png - Final processed UMAP",
    "",
    "### Statistics",
    "- filtering_summary_by_sample.csv - Filtering statistics per sample",
    "- resolution_analysis.csv - Clustering resolution analysis",
    "- cell_cycle_summary.csv - Cell cycle phase distribution",
    "- final_processing_metrics.rds - Complete processing metrics"
  )
  
  writeLines(report_lines, file.path(output_dir, "Processing_Report_Streamlined.md"))
}

create_processing_report(output_dir, final_metrics, QC_PARAMS)

# Final message
end_time <- Sys.time()
cat("\n=== PROCESSING COMPLETE ===\n")
cat("Processed Seurat object saved as: processed_seurat.rds\n")
cat("Outputs saved in:", output_dir, "\n")
cat("Processing report: Processing_Report_Streamlined.md\n")
cat("Analysis completed at:", as.character(end_time), "\n")

print("Streamlined scRNA-seq processing pipeline completed successfully!")

# Final garbage collection
gc(verbose = TRUE, reset = TRUE)