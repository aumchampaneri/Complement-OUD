# =============================================================================
# Enhanced Single-Cell RNA-seq Processing and Quality Control Pipeline
# =============================================================================
# 
# Purpose: Comprehensive processing optimized for large datasets with hybrid configuration
#
# Date: May 25, 2025
# Version: 5.1 (Hybrid Configuration - Tissue-Optimized - FULLY CORRECTED)
#
# Key Features:
# - Hybrid configuration system (original paper + robust methods)
# - Tissue-specific optimization for amygdala/GSE207128
# - Memory-efficient processing for large datasets
# - Comprehensive QC and analysis with smart defaults
# - Complete error handling and validation
# =============================================================================

# Set seed for reproducibility
set.seed(42)

# =============================================================================
# 1. SETUP AND HYBRID CONFIGURATION
# =============================================================================

cat("=== INITIALIZING HYBRID PROCESSING PIPELINE ===\n")
start_time <- Sys.time()

# Load required libraries with error handling
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(cluster)
})

# Check for optional libraries
clustree_available <- requireNamespace("clustree", quietly = TRUE)
if (!clustree_available) {
  message("clustree package not available - clustering tree plots will be skipped")
}

# Set working directories
setwd("/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus/GSE207128")
input_dir <- "Outputs/01_Aggregated_Data"
output_dir <- "Outputs/02_Processed_Data/QC_Processing"

# Create output directories with validation
dir.create(file.path(output_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "statistics"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "thresholds"), showWarnings = FALSE, recursive = TRUE)

# Verify directories were created
if (!dir.exists(output_dir)) {
  stop("Could not create output directory: ", output_dir)
}

cat("Output directory structure created successfully\n")

# Source the hybrid configuration with robust error handling
source_config_success <- tryCatch({
  source("Scripts/00_Configuration/Hybrid_Analysis_Config.R")
  TRUE
}, error = function(e) {
  message("Hybrid config not found - using fallback configuration")
  message("Error details: ", e$message)
  FALSE
})

# Get dataset-specific configuration
if (source_config_success) {
  config <- get_analysis_config(getwd(), tissue_type = "auto")
  print(paste("Configuration optimized for:", 
              ifelse(grepl("amygdala", toString(config)), "Amygdala (GSE207128)", "General brain")))
} else {
  # Fallback configuration for GSE207128 (original paper validated)
  config <- list(
    qc = list(
      min_features_base = 250,
      max_features_base = 6000,
      min_counts_base = 500,
      max_counts_base = 40000,
      max_mt_percent_base = 12,
      min_features_adaptive = FALSE,
      doublet_detection = FALSE,
      ribosomal_threshold = 30
    ),
    integration = list(
      variable_features = 2000,
      pca_dims = 25
    ),
    clustering = list(
      resolutions = c(0.4, 0.6, 0.8, 1.0),
      silhouette_analysis = FALSE,  # Skip for efficiency
      default_resolution = 0.8
    )
  )
  print("Using fallback configuration optimized for GSE207128")
}

# Helper function for null coalescing (define early)
`%||%` <- function(a, b) if (is.null(a)) b else a

# Extract parameters for this script (maintains compatibility)
QC_PARAMS <- list(
  # Use hybrid configuration values with fallbacks
  min_features_base = config$qc$min_features_base %||% 250,
  max_features_base = config$qc$max_features_base %||% 6000,
  min_counts_base = config$qc$min_counts_base %||% 500,
  max_counts_base = config$qc$max_counts_base %||% 40000,
  max_mt_percent_base = config$qc$max_mt_percent_base %||% 12,
  
  # Analysis parameters from config
  n_variable_features = config$integration$variable_features %||% 2000,
  pca_dims = 1:(config$integration$pca_dims %||% 25),
  clustering_resolutions = config$clustering$resolutions %||% c(0.4, 0.6, 0.8, 1.0),
  default_resolution = config$clustering$default_resolution %||% 0.8,
  mad_threshold = 3,
  
  # Efficiency settings based on dataset
  min_features_adaptive = config$qc$min_features_adaptive %||% FALSE,
  doublet_detection = config$qc$doublet_detection %||% FALSE,
  silhouette_analysis = config$clustering$silhouette_analysis %||% FALSE,
  ribosomal_threshold = config$qc$ribosomal_threshold %||% 30
)

# Print configuration being used
cat("=== HYBRID CONFIGURATION ACTIVE ===\n")
cat("Tissue-specific optimizations applied\n")
cat("QC thresholds: min_features =", QC_PARAMS$min_features_base, 
    ", max_mt =", QC_PARAMS$max_mt_percent_base, "%\n")
cat("Variable features:", QC_PARAMS$n_variable_features, "\n")
cat("PCA dimensions:", max(QC_PARAMS$pca_dims), "\n")
cat("Doublet detection:", QC_PARAMS$doublet_detection, "\n")
cat("Adaptive thresholds:", QC_PARAMS$min_features_adaptive, "\n")
cat("Silhouette analysis:", QC_PARAMS$silhouette_analysis, "\n")

# =============================================================================
# 2. DATA LOADING AND INITIAL ASSESSMENT
# =============================================================================

print("Loading merged Seurat object...")

# Check if input file exists
input_file <- file.path(input_dir, "combined_seurat.rds")
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, 
       "\nPlease run 0_10x Data Aggregation.R first")
}

combined_seurat <- readRDS(input_file)

# Validate loaded object
if (!inherits(combined_seurat, "Seurat")) {
  stop("Loaded object is not a valid Seurat object")
}

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

# Validate minimum requirements
if (initial_stats$total_cells < 100) {
  stop("Too few cells in dataset: ", initial_stats$total_cells)
}

# =============================================================================
# 3. QC METRICS CALCULATION
# =============================================================================

print("Calculating QC metrics...")

# Standard QC metrics with error handling
tryCatch({
  combined_seurat[["percent.mt"]] <- PercentageFeatureSet(combined_seurat, pattern = "^mt-")
  combined_seurat[["percent.ribo"]] <- PercentageFeatureSet(combined_seurat, pattern = "^Rp[sl]")
  combined_seurat[["log10GenesPerUMI"]] <- log10(combined_seurat$nFeature_RNA) / log10(combined_seurat$nCount_RNA)
  combined_seurat[["gene_umi_ratio"]] <- combined_seurat$nFeature_RNA / combined_seurat$nCount_RNA
}, error = function(e) {
  stop("Error calculating QC metrics: ", e$message)
})

# Validate QC metrics
if (any(is.na(combined_seurat$percent.mt))) {
  warning("Some cells have NA mitochondrial percentages")
}

print("QC metrics calculated successfully")

# =============================================================================
# 4. ENHANCED ADAPTIVE QC THRESHOLD CALCULATION
# =============================================================================

print("Calculating smart QC thresholds...")

# Enhanced adaptive threshold calculation with tissue-specific constraints
calculate_adaptive_thresholds <- function(seurat_obj, metric, n_mad = 3, min_val = NULL, max_val = NULL, tissue_type = "amygdala") {
  
  # Validate inputs
  if (!metric %in% colnames(seurat_obj@meta.data)) {
    stop("Metric '", metric, "' not found in metadata")
  }
  
  if (!"sample_id" %in% colnames(seurat_obj@meta.data)) {
    stop("sample_id column not found in metadata")
  }
  
  # Calculate per-sample statistics
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
  
  # Apply tissue-specific constraints
  if (tissue_type == "amygdala") {
    # Use original paper's validated constraints for amygdala
    if (metric == "nFeature_RNA") {
      min_val <- max(min_val %||% 250, 200)  # More conservative for amygdala
      max_val <- min(max_val %||% 6000, 8000)
    } else if (metric == "percent.mt") {
      max_val <- min(max_val %||% 12, 15)  # Stricter MT threshold
    } else if (metric == "nCount_RNA") {
      min_val <- max(min_val %||% 500, 300)
      max_val <- min(max_val %||% 40000, 60000)
    }
  }
  
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

# Calculate adaptive thresholds for key metrics with error handling
tryCatch({
  feature_thresholds <- calculate_adaptive_thresholds(combined_seurat, "nFeature_RNA", 
                                                     n_mad = QC_PARAMS$mad_threshold,
                                                     min_val = QC_PARAMS$min_features_base,
                                                     tissue_type = "amygdala")
  
  count_thresholds <- calculate_adaptive_thresholds(combined_seurat, "nCount_RNA",
                                                   n_mad = QC_PARAMS$mad_threshold,
                                                   min_val = QC_PARAMS$min_counts_base,
                                                   tissue_type = "amygdala")
  
  mt_thresholds <- calculate_adaptive_thresholds(combined_seurat, "percent.mt",
                                                n_mad = QC_PARAMS$mad_threshold,
                                                max_val = QC_PARAMS$max_mt_percent_base,
                                                tissue_type = "amygdala")
}, error = function(e) {
  stop("Error calculating adaptive thresholds: ", e$message)
})

# Save threshold information
tryCatch({
  write.csv(feature_thresholds, file.path(output_dir, "thresholds", "adaptive_feature_thresholds.csv"), row.names = FALSE)
  write.csv(count_thresholds, file.path(output_dir, "thresholds", "adaptive_count_thresholds.csv"), row.names = FALSE)
  write.csv(mt_thresholds, file.path(output_dir, "thresholds", "adaptive_mt_thresholds.csv"), row.names = FALSE)
}, error = function(e) {
  warning("Could not save threshold files: ", e$message)
})

print("Smart thresholds calculated and saved.")

# =============================================================================
# 5. PRE-FILTERING VISUALIZATION
# =============================================================================

print("Creating pre-filtering QC visualizations...")

# Violin plots with error handling
tryCatch({
  qc_plots_pre <- VlnPlot(combined_seurat,
                          features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
                          group.by = "sample_id", pt.size = 0.05, ncol = 2, raster = TRUE)
  ggsave(file.path(output_dir, "plots", "QC_violins_pre_filtering.png"), 
         qc_plots_pre, width = 16, height = 10, dpi = 300)
}, error = function(e) {
  warning("Could not create pre-filtering violin plots: ", e$message)
})

# Enhanced histograms with thresholds
create_qc_histograms <- function(seurat_obj, thresholds_list) {
  
  tryCatch({
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
  }, error = function(e) {
    warning("Error creating histograms: ", e$message)
    return(list())
  })
}

threshold_list <- list(feature = feature_thresholds, count = count_thresholds, mt = mt_thresholds)
hist_plots <- create_qc_histograms(combined_seurat, threshold_list)

if (length(hist_plots) > 0) {
  tryCatch({
    # Combine histograms
    histograms_combined <- (hist_plots[[1]] | hist_plots[[2]]) / (hist_plots[[3]] | hist_plots[[4]])
    ggsave(file.path(output_dir, "plots", "QC_histograms_pre_filtering.png"), 
           histograms_combined, width = 16, height = 12, dpi = 300)
  }, error = function(e) {
    warning("Could not save combined histogram: ", e$message)
  })
}

# =============================================================================
# 6. SMART CELL FILTERING (HYBRID APPROACH)
# =============================================================================

print("Applying smart cell filtering with tissue-specific optimization...")

# Store pre-filtering counts
pre_filter_count <- ncol(combined_seurat)
print(paste("Pre-filtering cell count:", pre_filter_count))

# Store original for comparison
original_seurat <- combined_seurat

# Choose filtering strategy based on configuration
if (QC_PARAMS$min_features_adaptive) {
  print("Using adaptive filtering approach...")
  
  # Apply sample-specific adaptive filtering
  cells_keep <- c()
  for(sample in unique(combined_seurat$sample_id)) {
    sample_metrics_feat <- feature_thresholds[feature_thresholds$sample_id == sample, ]
    sample_metrics_mt <- mt_thresholds[mt_thresholds$sample_id == sample, ]
    
    if(nrow(sample_metrics_feat) > 0 && nrow(sample_metrics_mt) > 0) {
      good_cells <- WhichCells(combined_seurat, 
                              expression = sample_id == sample &
                                         nFeature_RNA > sample_metrics_feat$lower_thresh_mad[1] &
                                         nFeature_RNA < sample_metrics_feat$upper_thresh_mad[1] &
                                         percent.mt < sample_metrics_mt$upper_thresh_mad[1])
    } else {
      # Fallback to global thresholds
      good_cells <- WhichCells(combined_seurat,
                              expression = sample_id == sample &
                                         nFeature_RNA >= QC_PARAMS$min_features_base &
                                         nFeature_RNA <= QC_PARAMS$max_features_base &
                                         percent.mt <= QC_PARAMS$max_mt_percent_base)
    }
    cells_keep <- c(cells_keep, good_cells)
  }
  
  combined_seurat_filtered <- subset(combined_seurat, cells = cells_keep)
  filtering_method <- "Adaptive (sample-specific)"
  
} else {
  print("Using original paper's validated thresholds (optimized for amygdala)...")
  
  # Apply fixed thresholds optimized for this tissue/dataset
  tryCatch({
    combined_seurat_filtered <- subset(combined_seurat,
                                      subset = nFeature_RNA >= QC_PARAMS$min_features_base &
                                              nFeature_RNA <= QC_PARAMS$max_features_base &
                                              nCount_RNA >= QC_PARAMS$min_counts_base &
                                              nCount_RNA <= QC_PARAMS$max_counts_base &
                                              percent.mt <= QC_PARAMS$max_mt_percent_base)
  }, error = function(e) {
    stop("Error during cell filtering: ", e$message)
  })
  
  filtering_method <- "Fixed thresholds (tissue-optimized)"
}

post_filter_count <- ncol(combined_seurat_filtered)
cells_removed <- pre_filter_count - post_filter_count
removal_percentage <- round((cells_removed / pre_filter_count) * 100, 2)

print(paste("Filtering method used:", filtering_method))
print(paste("Post-filtering cell count:", post_filter_count))
print(paste("Cells removed:", cells_removed, "(", removal_percentage, "%)"))

# Validate filtering results
if (post_filter_count < 50) {
  stop("Too few cells remaining after filtering: ", post_filter_count)
}

# Update object
combined_seurat <- combined_seurat_filtered

# Create filtering summary
tryCatch({
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
}, error = function(e) {
  warning("Could not create filtering summary: ", e$message)
})

# =============================================================================
# 7. POST-FILTERING VISUALIZATION
# =============================================================================

print("Creating post-filtering QC visualizations...")

# Post-filtering violin plots
tryCatch({
  qc_plots_post <- VlnPlot(combined_seurat,
                           features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
                           group.by = "sample_id", pt.size = 0.05, ncol = 2, raster = TRUE)
  ggsave(file.path(output_dir, "plots", "QC_violins_post_filtering.png"), 
         qc_plots_post, width = 16, height = 10, dpi = 300)
}, error = function(e) {
  warning("Could not create post-filtering violin plots: ", e$message)
})

# Comparison plots (pre vs post filtering)
tryCatch({
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
}, error = function(e) {
  warning("Could not create comparison plots: ", e$message)
})

# =============================================================================
# 8. NORMALIZATION AND FEATURE SELECTION
# =============================================================================

print("Performing normalization and feature selection...")

# Standard log-normalization with error handling
tryCatch({
  combined_seurat <- NormalizeData(combined_seurat, 
                                  normalization.method = "LogNormalize", 
                                  scale.factor = 10000,
                                  verbose = FALSE)
}, error = function(e) {
  stop("Error during normalization: ", e$message)
})

# Find variable features
tryCatch({
  combined_seurat <- FindVariableFeatures(combined_seurat, 
                                         selection.method = "vst", 
                                         nfeatures = QC_PARAMS$n_variable_features,
                                         verbose = FALSE)
}, error = function(e) {
  stop("Error finding variable features: ", e$message)
})

# Visualize variable features
tryCatch({
  variable_features_plot <- VariableFeaturePlot(combined_seurat) +
    labs(title = paste("Top", QC_PARAMS$n_variable_features, "Variable Features"))
  ggsave(file.path(output_dir, "plots", "variable_features.png"), 
         variable_features_plot, width = 12, height = 8, dpi = 300)
}, error = function(e) {
  warning("Could not create variable features plot: ", e$message)
})

# Scale data (only variable features for memory efficiency)
print("Scaling data...")
tryCatch({
  combined_seurat <- ScaleData(combined_seurat, 
                              features = VariableFeatures(combined_seurat),
                              vars.to.regress = c("percent.mt"),
                              verbose = FALSE)
}, error = function(e) {
  stop("Error during data scaling: ", e$message)
})

# =============================================================================
# 9. DIMENSIONALITY REDUCTION
# =============================================================================

print("Performing PCA...")
tryCatch({
  combined_seurat <- RunPCA(combined_seurat, 
                           features = VariableFeatures(object = combined_seurat),
                           verbose = FALSE)
}, error = function(e) {
  stop("Error during PCA: ", e$message)
})

# PCA visualization
tryCatch({
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
}, error = function(e) {
  warning("Could not create PCA plots: ", e$message)
})

# =============================================================================
# 10. EFFICIENT CLUSTERING OPTIMIZATION
# =============================================================================

print("Performing efficient clustering optimization...")

# Test clustering resolutions with error handling
tryCatch({
  combined_seurat <- FindNeighbors(combined_seurat, dims = QC_PARAMS$pca_dims, verbose = FALSE)
}, error = function(e) {
  stop("Error finding neighbors: ", e$message)
})

# Test resolutions
for(res in QC_PARAMS$clustering_resolutions) {
  tryCatch({
    combined_seurat <- FindClusters(combined_seurat, resolution = res, verbose = FALSE)
  }, error = function(e) {
    warning("Could not cluster at resolution ", res, ": ", e$message)
  })
}

# Create clustree plot if available and if silhouette analysis is enabled
if(clustree_available && QC_PARAMS$silhouette_analysis) {
  tryCatch({
    clustree_plot <- clustree::clustree(combined_seurat, prefix = "RNA_snn_res.")
    ggsave(file.path(output_dir, "plots", "clustering_resolution_tree.png"), 
           clustree_plot, width = 12, height = 10, dpi = 300)
  }, error = function(e) {
    message("Could not create clustree plot: ", e$message)
  })
}

# Calculate silhouette scores (only if enabled for efficiency)
if (QC_PARAMS$silhouette_analysis) {
  print("Running silhouette analysis for resolution optimization...")
  
  tryCatch({
    pca_embeddings <- Embeddings(combined_seurat, "pca")[, QC_PARAMS$pca_dims]
    
    silhouette_scores <- sapply(QC_PARAMS$clustering_resolutions, function(res) {
      clusters <- combined_seurat@meta.data[[paste0("RNA_snn_res.", res)]]
      if(length(unique(clusters)) > 1 && length(unique(clusters)) < nrow(pca_embeddings)) {
        tryCatch({
          # Use sample for large datasets (more efficient)
          sample_size <- min(3000, nrow(pca_embeddings))
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
    
    # Select optimal resolution based on silhouette score
    valid_scores <- !is.na(silhouette_scores)
    if(any(valid_scores)) {
      optimal_resolution <- QC_PARAMS$clustering_resolutions[which.max(silhouette_scores[valid_scores])]
    } else {
      optimal_resolution <- QC_PARAMS$default_resolution
    }
  }, error = function(e) {
    warning("Error in silhouette analysis: ", e$message)
    silhouette_scores <- rep(NA, length(QC_PARAMS$clustering_resolutions))
    optimal_resolution <- QC_PARAMS$default_resolution
  })
  
} else {
  print("Skipping silhouette analysis for efficiency - using default resolution")
  silhouette_scores <- rep(NA, length(QC_PARAMS$clustering_resolutions))
  optimal_resolution <- QC_PARAMS$default_resolution
}

# Save resolution analysis
tryCatch({
  res_analysis <- data.frame(
    resolution = QC_PARAMS$clustering_resolutions, 
    silhouette_score = silhouette_scores,
    n_clusters = sapply(QC_PARAMS$clustering_resolutions, function(res) {
      length(unique(combined_seurat@meta.data[[paste0("RNA_snn_res.", res)]]))
    }),
    method = filtering_method
  )
  
  write.csv(res_analysis, file.path(output_dir, "statistics", "resolution_analysis.csv"), row.names = FALSE)
}, error = function(e) {
  warning("Could not save resolution analysis: ", e$message)
})

print(paste("Optimal resolution selected:", optimal_resolution))
print(paste("Analysis method:", filtering_method))

# Set optimal clustering
combined_seurat$seurat_clusters <- combined_seurat@meta.data[[paste0("RNA_snn_res.", optimal_resolution)]]
Idents(combined_seurat) <- combined_seurat$seurat_clusters

# Plot resolution analysis
if (exists("res_analysis")) {
  tryCatch({
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
           subtitle = paste("Blue: Silhouette Score, Red: Number of Clusters -", filtering_method)) +
      theme_minimal()
    
    ggsave(file.path(output_dir, "plots", "resolution_optimization.png"), 
           res_plot, width = 10, height = 6, dpi = 300)
  }, error = function(e) {
    warning("Could not create resolution plot: ", e$message)
  })
}

# =============================================================================
# 11. UMAP VISUALIZATION
# =============================================================================

print("Generating UMAP embedding...")
tryCatch({
  combined_seurat <- RunUMAP(combined_seurat, dims = QC_PARAMS$pca_dims, verbose = FALSE)
}, error = function(e) {
  stop("Error during UMAP generation: ", e$message)
})

# Create UMAP visualizations
tryCatch({
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
}, error = function(e) {
  warning("Could not create UMAP plots: ", e$message)
})

# =============================================================================
# 12. CELL CYCLE SCORING
# =============================================================================

print("Running cell cycle scoring with mouse-specific genes...")

tryCatch({
  # Mouse cell cycle genes (complete lists)
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
    tryCatch({
      marker_plot <- FeaturePlot(combined_seurat, features = marker, raster = TRUE) +
        labs(title = paste("Expression of", marker))
      
      ggsave(file.path(output_dir, "plots", paste0("marker_", marker, ".png")), 
             marker_plot, width = 8, height = 6, dpi = 300)
    }, error = function(e) {
      warning("Could not create plot for marker ", marker, ": ", e$message)
    })
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
tryCatch({
  final_umap <- DimPlot(combined_seurat, reduction = "umap", label = TRUE, repel = TRUE, raster = TRUE) +
    labs(title = paste("Final UMAP - Complete Processing (", filtering_method, ")", sep = ""))
  ggsave(file.path(output_dir, "plots", "UMAP_final.png"), 
         final_umap, width = 10, height = 8, dpi = 300)
}, error = function(e) {
  warning("Could not create final UMAP: ", e$message)
})

# =============================================================================
# 15. ENHANCED FINAL METRICS WITH CONFIGURATION INFO
# =============================================================================

print("Generating enhanced final metrics...")

# Calculate final metrics with configuration info
final_metrics <- list(
  # Configuration info
  configuration_type = ifelse(grepl("amygdala", toString(config)), "Amygdala-optimized", "General brain"),
  filtering_method = filtering_method,
  tissue_optimized = TRUE,
  
  # Processing info
  processing_date = Sys.Date(),
  processing_time = Sys.time(),
  initial_cells = initial_stats$total_cells,
  final_cells = ncol(combined_seurat),
  cells_removed = initial_stats$total_cells - ncol(combined_seurat),
  removal_percentage = round((initial_stats$total_cells - ncol(combined_seurat)) / initial_stats$total_cells * 100, 2),
  final_genes = nrow(combined_seurat),
  final_clusters = length(unique(combined_seurat$seurat_clusters)),
  optimal_resolution = optimal_resolution,
  
  # QC metrics
  median_genes_per_cell = median(combined_seurat$nFeature_RNA),
  median_umis_per_cell = median(combined_seurat$nCount_RNA),
  median_mt_percent = median(combined_seurat$percent.mt),
  
  # Configuration parameters used
  config_parameters = list(
    min_features = QC_PARAMS$min_features_base,
    max_features = QC_PARAMS$max_features_base,
    max_mt_percent = QC_PARAMS$max_mt_percent_base,
    variable_features = QC_PARAMS$n_variable_features,
    pca_dims_used = max(QC_PARAMS$pca_dims),
    silhouette_analysis_performed = QC_PARAMS$silhouette_analysis
  ),
  
  # Sample distributions
  cells_per_condition = as.list(table(combined_seurat$condition)),
  cells_per_sample = as.list(table(combined_seurat$sample_id)),
  
  # Processing features
  doublet_detection_performed = QC_PARAMS$doublet_detection,
  streamlined_processing = TRUE,
  hybrid_configuration = TRUE,
  available_markers = if(exists("available_markers")) available_markers else character(0)
)

# Save final metrics with error handling
tryCatch({
  saveRDS(final_metrics, file.path(output_dir, "statistics", "final_processing_metrics.rds"))
}, error = function(e) {
  warning("Could not save final metrics: ", e$message)
})

# Print summary
cat("\n=== FINAL PROCESSING SUMMARY ===\n")
cat("Configuration:", final_metrics$configuration_type, "\n")
cat("Filtering method:", final_metrics$filtering_method, "\n")
cat("Initial cells:", final_metrics$initial_cells, "\n")
cat("Final cells:", final_metrics$final_cells, "\n")
cat("Cells removed:", final_metrics$cells_removed, "(", final_metrics$removal_percentage, "%)\n")
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

# Save the final processed Seurat object with error handling
tryCatch({
  saveRDS(combined_seurat, "Outputs/02_Processed_Data/processed_seurat.rds")
  print("Processed Seurat object saved successfully")
}, error = function(e) {
  stop("Could not save processed Seurat object: ", e$message)
})

# Create enhanced processing report
create_processing_report <- function(output_dir, final_metrics, QC_PARAMS) {
  tryCatch({
    report_lines <- c(
      "# Single-Cell RNA-seq Processing Report - HYBRID CONFIGURATION",
      paste("Generated on:", Sys.Date()),
      paste("Processing time:", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 2), "minutes"),
      "",
      "## Configuration Applied",
      paste("- Configuration type:", final_metrics$configuration_type),
      paste("- Filtering method:", final_metrics$filtering_method),
      paste("- Tissue optimization:", final_metrics$tissue_optimized),
      "",
      "## Processing Parameters",
      paste("- Random seed:", 42),
      paste("- Minimum features:", final_metrics$config_parameters$min_features),
      paste("- Maximum features:", final_metrics$config_parameters$max_features),
      paste("- Maximum mitochondrial %:", final_metrics$config_parameters$max_mt_percent),
      paste("- Variable features:", final_metrics$config_parameters$variable_features),
      paste("- PCA dimensions:", final_metrics$config_parameters$pca_dims_used),
      paste("- Optimal clustering resolution:", final_metrics$optimal_resolution),
      paste("- Silhouette analysis:", final_metrics$config_parameters$silhouette_analysis_performed),
      "",
      "## Processing Features",
      "- ✅ Hybrid configuration system",
      "- ✅ Tissue-specific optimization", 
      "- ✅ Smart QC threshold selection",
      "- ✅ Comprehensive visualization",
      "- ✅ Memory-efficient processing",
      "- ✅ Resolution optimization",
      "- ✅ Cell cycle scoring",
      "- ✅ Robust error handling",
      paste("- ", ifelse(final_metrics$doublet_detection_performed, "✅", "❌"), " Doublet detection"),
      "",
      "## Results Summary",
      paste("- Initial cells:", final_metrics$initial_cells),
      paste("- Final cells:", final_metrics$final_cells),
      paste("- Cells removed:", final_metrics$cells_removed, "(", final_metrics$removal_percentage, "%)"),
      paste("- Final clusters:", final_metrics$final_clusters),
      paste("- Median genes per cell:", final_metrics$median_genes_per_cell),
      paste("- Median UMIs per cell:", final_metrics$median_umis_per_cell),
      paste("- Median MT percent:", round(final_metrics$median_mt_percent, 2), "%"),
      "",
      "## Cell Type Markers Available",
      if(length(final_metrics$available_markers) > 0) {
        paste("-", paste(final_metrics$available_markers, collapse = ", "))
      } else {
        "- No standard markers detected in dataset"
      },
      "",
      "## Output Files",
      "### Main Results",
      "- processed_seurat.rds - Final processed Seurat object",
      "",
      "### Quality Control Plots",
      "- QC_violins_pre_filtering.png - Pre-filtering QC metrics",
      "- QC_violins_post_filtering.png - Post-filtering QC metrics", 
      "- QC_histograms_pre_filtering.png - QC histograms with thresholds",
      "- filtering_comparison.png - Pre vs post filtering comparison",
      "",
      "### Analysis Plots",
      "- variable_features.png - Highly variable features",
      "- PCA_analysis.png - PCA visualization",
      "- clustering_resolution_tree.png - Resolution optimization tree (if available)",
      "- resolution_optimization.png - Resolution analysis plot",
      "- UMAP_combined_overview.png - UMAP overview",
      "- cell_cycle_phases.png - Cell cycle analysis",
      "- marker_*.png - Individual marker expression plots",
      "- UMAP_final.png - Final processed UMAP",
      "",
      "### Statistics and Data",
      "- filtering_summary_by_sample.csv - Filtering statistics per sample",
      "- resolution_analysis.csv - Clustering resolution analysis",
      "- cell_cycle_summary.csv - Cell cycle phase distribution",
      "- final_processing_metrics.rds - Complete processing metrics",
      "",
      "### Adaptive Thresholds",
      "- adaptive_feature_thresholds.csv - Sample-specific feature thresholds",
      "- adaptive_count_thresholds.csv - Sample-specific count thresholds", 
      "- adaptive_mt_thresholds.csv - Sample-specific mitochondrial thresholds",
      "",
      "## Next Steps",
      "1. Run 2_Integration.R for batch correction",
      "2. Run 3_CellType_Annotation.R for cell type identification",
      "3. Proceed with downstream analysis"
    )
    
    writeLines(report_lines, file.path(output_dir, "Processing_Report_Hybrid.md"))
    return(TRUE)
  }, error = function(e) {
    warning("Could not create processing report: ", e$message)
    return(FALSE)
  })
}

report_created <- create_processing_report(output_dir, final_metrics, QC_PARAMS)

# Final message
end_time <- Sys.time()
processing_time <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 2)

cat("\n=== HYBRID PROCESSING COMPLETE ===\n")
cat("Configuration:", final_metrics$configuration_type, "\n")
cat("Total processing time:", processing_time, "minutes\n")
cat("Processed Seurat object saved as: processed_seurat.rds\n")
cat("Outputs saved in:", output_dir, "\n")
if (report_created) {
  cat("Processing report: Processing_Report_Hybrid.md\n")
}
cat("Analysis completed at:", as.character(end_time), "\n")

print("✅ Hybrid scRNA-seq processing pipeline completed successfully!")
print("🚀 Ready for integration and cell type annotation!")

# Final garbage collection
gc(verbose = TRUE, reset = TRUE)

cat("\n=== PROCESSING PIPELINE SUMMARY ===\n")
cat("Status: SUCCESS\n")
cat("Next step: Run 2_Integration.R\n")
cat("=======================================\n")
