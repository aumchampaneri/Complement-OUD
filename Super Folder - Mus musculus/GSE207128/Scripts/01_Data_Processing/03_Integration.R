# =============================================================================
# Enhanced Single-Cell RNA-seq Integration Pipeline
# =============================================================================
# 
# Purpose: Robust integration with Harmony or fallback methods
#
# Date: May 25, 2025
# Version: 5.1 (Hybrid Configuration - Tissue-Optimized - FULLY REWORKED)
#
# Key Features:
# - Hybrid integration approach (Harmony + fallbacks)
# - Tissue-specific optimization for amygdala/GSE207128
# - Automatic method selection based on availability
# - Comprehensive error handling and validation
# - Memory-efficient processing for large datasets
# - Detailed integration metrics and visualizations
# =============================================================================

# Set seed for reproducibility
set.seed(42)

# =============================================================================
# 1. SETUP AND CONFIGURATION
# =============================================================================

cat("=== INITIALIZING HYBRID INTEGRATION PIPELINE ===\n")
start_time <- Sys.time()

# Helper function for null coalescing
`%||%` <- function(a, b) if (is.null(a)) b else a

# Load configuration with robust error handling
source_config_success <- tryCatch({
  source("/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus/GSE207128/Scripts/00_Configuration/Hybrid_Analysis_Config.R")
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
  # Fallback configuration for GSE207128
  config <- list(
    integration = list(
      pca_dims = 25,
      integration_dims = 20,
      variable_features = 2000
    ),
    clustering = list(
      resolutions = c(0.4, 0.6, 0.8, 1.0),
      default_resolution = 0.8
    ),
    dataset_info = list(
      type = "GSE207128",
      tissue = "amygdala",
      optimization = "fallback_configuration"
    )
  )
  print("Using fallback configuration optimized for GSE207128")
}

# Load required libraries with error handling
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# Check for integration methods availability
harmony_available <- FALSE
fastmnn_available <- FALSE

# Test Harmony availability
tryCatch({
  library(harmony)
  harmony_available <- TRUE
  cat("✓ Harmony is available - using Harmony integration\n")
}, error = function(e) {
  cat("✗ Harmony not available - will use batch regression\n")
  harmony_available <- FALSE
})

# Test batchelor/fastMNN availability (optional)
tryCatch({
  library(batchelor)
  fastmnn_available <- TRUE
  cat("✓ FastMNN is available as backup option\n")
}, error = function(e) {
  cat("✗ FastMNN not available\n")
  fastmnn_available <- FALSE
})

# Set working directories with validation
setwd("/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus/GSE207128")
input_path <- "Outputs/02_Processed_Data/processed_seurat.rds"
output_dir <- "Outputs/03_Integrated_Data/Integration"

# Create output directories
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "metrics"), showWarnings = FALSE, recursive = TRUE)

# Verify directories were created
if (!dir.exists(output_dir)) {
  stop("Could not create output directory: ", output_dir)
}

cat("Output directory structure created successfully\n")

# =============================================================================
# 2. DATA LOADING AND VALIDATION
# =============================================================================

print("=== STEP 2: DATA LOADING AND VALIDATION ===")

# Check if input file exists
if (!file.exists(input_path)) {
  stop("Processed Seurat object not found: ", input_path, 
       "\nPlease run 1_Processing and QC.R first")
}

# Load processed data with validation
seurat_obj <- readRDS(input_path)

# Validate loaded object
if (!inherits(seurat_obj, "Seurat")) {
  stop("Loaded object is not a valid Seurat object")
}

print(paste("✓ Loaded Seurat object with", ncol(seurat_obj), "cells and", nrow(seurat_obj), "genes"))

# Handle Seurat v5 objects (join layers if needed)
DefaultAssay(seurat_obj) <- "RNA"
if (inherits(seurat_obj[["RNA"]], "Assay5")) {
  print("Detected Seurat v5 object. Joining layers...")
  seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
}

# Validate required metadata columns
required_cols <- c("sample_id", "condition")
missing_cols <- setdiff(required_cols, colnames(seurat_obj@meta.data))

if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
}

# Check batch variables and sample distribution
available_batch_vars <- c("sample_id", "condition", "orig.ident")[
  c("sample_id", "condition", "orig.ident") %in% colnames(seurat_obj@meta.data)]

print("Available batch variables:")
print(available_batch_vars)

print("Sample distribution:")
sample_counts <- table(seurat_obj$sample_id)
print(sample_counts)

print("Condition distribution:")
condition_counts <- table(seurat_obj$condition)
print(condition_counts)

# Validate minimum requirements for integration
if (length(unique(seurat_obj$sample_id)) < 2) {
  warning("Only one sample detected - integration may not be necessary")
}

if (ncol(seurat_obj) < 100) {
  stop("Too few cells for reliable integration: ", ncol(seurat_obj))
}

# =============================================================================
# 3. PREPROCESSING VALIDATION AND PREPARATION
# =============================================================================

print("=== STEP 3: PREPROCESSING VALIDATION ===")

# Check if data is already normalized
max_data_value <- max(GetAssayData(seurat_obj, assay = "RNA", layer = "data"))
if (max_data_value > 10) {
  print("Normalizing data...")
  seurat_obj <- NormalizeData(seurat_obj, 
                             normalization.method = "LogNormalize",
                             scale.factor = 10000,
                             verbose = FALSE)
} else {
  print("✓ Data appears to be already normalized")
}

# Find variable features if not present or update count
current_var_features <- length(VariableFeatures(seurat_obj))
target_var_features <- config$integration$variable_features %||% 2000

if (current_var_features == 0 || current_var_features != target_var_features) {
  print(paste("Finding", target_var_features, "variable features..."))
  seurat_obj <- FindVariableFeatures(seurat_obj, 
                                    selection.method = "vst", 
                                    nfeatures = target_var_features,
                                    verbose = FALSE)
} else {
  print(paste("✓ Using existing", current_var_features, "variable features"))
}

# Scale data if not already scaled (for PCA)
if (!"scale.data" %in% slotNames(seurat_obj[["RNA"]]) || 
    length(seurat_obj[["RNA"]]@scale.data) == 0) {
  print("Scaling data for PCA...")
  
  # Use batch regression as fallback if no other integration method available
  vars_to_regress <- if (!harmony_available && !fastmnn_available) {
    c("percent.mt", "sample_id")  # Include sample_id for batch regression
  } else {
    c("percent.mt")  # Standard regression for mitochondrial effects
  }
  
  seurat_obj <- ScaleData(seurat_obj, 
                         features = VariableFeatures(seurat_obj),
                         vars.to.regress = vars_to_regress,
                         verbose = FALSE)
  print(paste("✓ Data scaled with regression of:", paste(vars_to_regress, collapse = ", ")))
} else {
  print("✓ Data appears to be already scaled")
}

# Run PCA if not present or update dimensions
target_pca_dims <- config$integration$pca_dims %||% 50

if (!"pca" %in% Reductions(seurat_obj) || 
    ncol(Embeddings(seurat_obj, "pca")) < target_pca_dims) {
  print(paste("Running PCA with", target_pca_dims, "dimensions..."))
  seurat_obj <- RunPCA(seurat_obj, 
                      features = VariableFeatures(seurat_obj),
                      npcs = target_pca_dims,
                      verbose = FALSE)
} else {
  print("✓ PCA already present")
}

# Determine optimal number of PCs to use
print("Determining optimal number of PCs...")
pct <- seurat_obj[["pca"]]@stdev / sum(seurat_obj[["pca"]]@stdev) * 100
cumu <- cumsum(pct)

# Find where cumulative variance > 90% and individual PC contributes < 5%
co1 <- which(cumu > 90 & pct < 5)[1]

# Find where consecutive PCs difference < 0.1%
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1

# Use the minimum of the two cutoffs, with sensible bounds
pcs_to_use <- min(co1, co2, target_pca_dims, na.rm = TRUE)

if (is.na(pcs_to_use) || pcs_to_use < 10) {
  pcs_to_use <- min(30, target_pca_dims)
  print("Using default number of PCs due to calculation issues")
}

print(paste("✓ Using", pcs_to_use, "principal components for integration"))

# Create PCA plots for quality assessment
tryCatch({
  elbow_plot <- ElbowPlot(seurat_obj, ndims = min(50, target_pca_dims)) +
    geom_vline(xintercept = pcs_to_use, color = "red", linetype = "dashed") +
    labs(title = "PCA Elbow Plot", subtitle = paste("Using", pcs_to_use, "PCs"))
  
  pca_batch_plot <- DimPlot(seurat_obj, reduction = "pca", group.by = "sample_id", raster = TRUE) +
    labs(title = "PCA by Sample (Pre-Integration)")
  
  pca_plots <- elbow_plot | pca_batch_plot
  ggsave(file.path(output_dir, "plots", "PCA_pre_integration.png"), 
         pca_plots, width = 16, height = 8, dpi = 300)
}, error = function(e) {
  warning("Could not create PCA plots: ", e$message)
})

# =============================================================================
# 4. INTEGRATION METHOD SELECTION AND EXECUTION
# =============================================================================

print("=== STEP 4: INTEGRATION METHOD SELECTION ===")

# Determine integration strategy
if (harmony_available) {
  integration_method <- "Harmony"
  print("✓ Using Harmony for integration")
} else if (fastmnn_available) {
  integration_method <- "FastMNN"
  print("✓ Using FastMNN for integration")
} else {
  integration_method <- "Batch Regression"
  print("✓ Using Batch Regression (already applied during scaling)")
}

# Execute integration based on available method
if (integration_method == "Harmony") {
  
  print("Running Harmony integration with tissue-specific parameters...")
  
  # Apply tissue-specific parameters
  if (grepl("amygdala", toString(config$dataset_info)) || grepl("GSE207128", toString(config$dataset_info))) {
    print("Applying amygdala-optimized Harmony parameters...")
    
    # Conservative parameters for amygdala tissue (original paper validated)
    seurat_integrated <- RunHarmony(seurat_obj, 
                                    group.by.vars = "sample_id",
                                    reduction = "pca",
                                    dims.use = 1:pcs_to_use,
                                    theta = 2,        # Moderate correction strength
                                    lambda = 1,       # Default diversity penalty
                                    sigma = 0.1,      # Width of soft kmeans clusters
                                    nclust = 50,      # Number of clusters in model
                                    tau = 0,          # Protection against overclustering small datasets
                                    block.size = 0.05, # What fraction of cells to update during clustering
                                    max.iter.harmony = 10,   # Maximum harmony iterations
                                    max.iter.cluster = 200,  # Maximum clustering iterations  
                                    epsilon.cluster = 1e-5,  # Convergence tolerance for clustering
                                    epsilon.harmony = 1e-4,  # Convergence tolerance for harmony
                                    plot_convergence = FALSE, # Don't plot for automated pipeline
                                    verbose = TRUE,
                                    reduction.save = "harmony")
    
    # Integration dimensions optimized for amygdala
    integration_dims_to_use <- min(pcs_to_use, config$integration$integration_dims %||% 20)
    
  } else {
    print("Applying standard Harmony parameters for general brain tissue...")
    
    # Standard parameters for other brain regions
    seurat_integrated <- RunHarmony(seurat_obj, 
                                    group.by.vars = "sample_id",
                                    reduction = "pca",
                                    dims.use = 1:pcs_to_use,
                                    theta = 2,
                                    max.iter.harmony = 20,
                                    plot_convergence = FALSE,
                                    verbose = TRUE,
                                    reduction.save = "harmony")
    
    integration_dims_to_use <- min(pcs_to_use, 30)
  }
  
  reduction_to_use <- "harmony"
  print(paste("✓ Harmony integration completed using", integration_dims_to_use, "dimensions"))

} else if (integration_method == "FastMNN") {
  
  print("Running FastMNN integration...")
  
  # Split by sample for FastMNN
  sample_list <- SplitObject(seurat_obj, split.by = "sample_id")
  
  # Run FastMNN
  tryCatch({
    # This is a simplified FastMNN approach - you may need to adjust
    mnn_result <- batchelor::fastMNN(sample_list, 
                                     d = pcs_to_use,
                                     k = 20,
                                     subset.row = VariableFeatures(seurat_obj))
    
    # Convert back to Seurat (this is complex and may need adjustment)
    seurat_integrated <- seurat_obj
    seurat_integrated[["mnn"]] <- CreateDimReducObject(embeddings = reducedDim(mnn_result),
                                                       key = "MNN_",
                                                       assay = "RNA")
    reduction_to_use <- "mnn"
    integration_dims_to_use <- min(pcs_to_use, 30)
    
    print("✓ FastMNN integration completed")
    }, error = function(e) {
      warning("FastMNN failed, falling back to batch regression: ", e$message)
      seurat_integrated <<- seurat_obj  # Use <<- for parent scope
      reduction_to_use <<- "pca"        # Use <<- for parent scope
      integration_dims_to_use <<- pcs_to_use  # Use <<- for parent scope
      integration_method <<- "Batch Regression (FastMNN fallback)"  # Use <<- for parent scope
    })

} else {
  # Batch regression already applied during scaling
  seurat_integrated <- seurat_obj
  reduction_to_use <- "pca"
  integration_dims_to_use <- pcs_to_use
  print("✓ Using batch regression (applied during data scaling)")
}

print(paste("Final integration method:", integration_method))
print(paste("Using reduction:", reduction_to_use))
print(paste("Using dimensions:", integration_dims_to_use))

# =============================================================================
# 5. POST-INTEGRATION CLUSTERING AND UMAP
# =============================================================================

print("=== STEP 5: POST-INTEGRATION CLUSTERING ===")

# Find neighbors using the integrated reduction
print("Finding neighbors on integrated data...")
seurat_integrated <- FindNeighbors(seurat_integrated, 
                                  reduction = reduction_to_use,
                                  dims = 1:integration_dims_to_use,
                                  verbose = FALSE)

# Test multiple clustering resolutions
print("Testing multiple clustering resolutions...")
clustering_resolutions <- config$clustering$resolutions %||% c(0.4, 0.6, 0.8, 1.0)

for(res in clustering_resolutions) {
  tryCatch({
    seurat_integrated <- FindClusters(seurat_integrated, 
                                     resolution = res, 
                                     verbose = FALSE)
    print(paste("✓ Clustering at resolution", res, "- found", 
                length(unique(seurat_integrated@meta.data[[paste0("RNA_snn_res.", res)]])), "clusters"))
  }, error = function(e) {
    warning("Could not cluster at resolution ", res, ": ", e$message)
  })
}

# Set default clustering
default_resolution <- config$clustering$default_resolution %||% 0.8
default_cluster_col <- paste0("RNA_snn_res.", default_resolution)

if (default_cluster_col %in% colnames(seurat_integrated@meta.data)) {
  seurat_integrated$seurat_clusters <- seurat_integrated@meta.data[[default_cluster_col]]
  Idents(seurat_integrated) <- seurat_integrated$seurat_clusters
  print(paste("✓ Set default clustering to resolution", default_resolution))
} else {
  # Use the last successful resolution
  available_res <- grep("RNA_snn_res", colnames(seurat_integrated@meta.data), value = TRUE)
  if (length(available_res) > 0) {
    seurat_integrated$seurat_clusters <- seurat_integrated@meta.data[[available_res[length(available_res)]]]
    Idents(seurat_integrated) <- seurat_integrated$seurat_clusters
    print(paste("✓ Using fallback clustering:", available_res[length(available_res)]))
  } else {
    stop("No successful clustering found")
  }
}

# Generate UMAP on integrated data
print("Generating UMAP on integrated data...")

# Set UMAP name based on integration method
umap_name <- switch(integration_method,
                   "Harmony" = "umap_harmony",
                   "FastMNN" = "umap_fastmnn", 
                   "umap_integrated")

if (startsWith(integration_method, "Batch Regression")) {
  umap_name <- "umap"
}

tryCatch({
  seurat_integrated <- RunUMAP(seurat_integrated, 
                              reduction = reduction_to_use,
                              dims = 1:integration_dims_to_use,
                              reduction.name = umap_name,
                              verbose = FALSE)
  print(paste("✓ UMAP generated:", umap_name))
}, error = function(e) {
  stop("Error generating UMAP: ", e$message)
})

# =============================================================================
# 6. INTEGRATION QUALITY ASSESSMENT
# =============================================================================

print("=== STEP 6: INTEGRATION QUALITY ASSESSMENT ===")

# Calculate integration metrics
print("Calculating integration quality metrics...")

# Sample mixing in clusters
sample_mixing <- table(seurat_integrated$seurat_clusters, seurat_integrated$sample_id)
condition_mixing <- table(seurat_integrated$seurat_clusters, seurat_integrated$condition)

# Calculate mixing entropy (higher = better mixing)
calculate_entropy <- function(count_matrix) {
  apply(count_matrix, 1, function(x) {
    p <- x / sum(x)
    p <- p[p > 0]  # Remove zeros to avoid log(0)
    -sum(p * log2(p))
  })
}

sample_entropy <- calculate_entropy(sample_mixing)
condition_entropy <- calculate_entropy(condition_mixing)

# Integration metrics summary
integration_metrics <- list(
  method = integration_method,
  dimensions_used = integration_dims_to_use,
  total_cells = ncol(seurat_integrated),
  n_clusters = length(unique(seurat_integrated$seurat_clusters)),
  n_samples = length(unique(seurat_integrated$sample_id)),
  n_conditions = length(unique(seurat_integrated$condition)),
  mean_sample_entropy = mean(sample_entropy, na.rm = TRUE),
  mean_condition_entropy = mean(condition_entropy, na.rm = TRUE),
  clustering_resolution = default_resolution
)

print("Integration metrics:")
print(integration_metrics)

# Save mixing matrices
write.csv(sample_mixing, file.path(output_dir, "metrics", "sample_mixing_by_cluster.csv"))
write.csv(condition_mixing, file.path(output_dir, "metrics", "condition_mixing_by_cluster.csv"))

# =============================================================================
# 7. COMPREHENSIVE VISUALIZATIONS
# =============================================================================

print("=== STEP 7: CREATING VISUALIZATIONS ===")

# Create pre-integration UMAP for comparison (if not using batch regression)
if (integration_method != "Batch Regression" && !startsWith(integration_method, "Batch Regression")) {
  print("Creating pre-integration UMAP for comparison...")
  
  tryCatch({
    seurat_pre <- seurat_obj
    seurat_pre <- RunUMAP(seurat_pre, 
                         reduction = "pca",
                         dims = 1:pcs_to_use,
                         reduction.name = "umap_pre_integration",
                         verbose = FALSE)
    
    # Before/after integration comparison
    p1 <- DimPlot(seurat_pre, reduction = "umap_pre_integration", 
                  group.by = "sample_id", raster = TRUE) +
      labs(title = "Before Integration", subtitle = "PCA → UMAP") +
      theme(legend.position = "bottom")
    
    p2 <- DimPlot(seurat_integrated, reduction = umap_name, 
                  group.by = "sample_id", raster = TRUE) +
      labs(title = paste("After", integration_method), 
           subtitle = paste(reduction_to_use, "→ UMAP")) +
      theme(legend.position = "bottom")
    
    integration_comparison <- p1 / p2
    
    # Clean up
    rm(seurat_pre)
    
  }, error = function(e) {
    warning("Could not create pre-integration comparison: ", e$message)
    
    # Fallback: just show post-integration
    p1 <- DimPlot(seurat_integrated, reduction = umap_name, 
                  group.by = "sample_id", raster = TRUE) +
      labs(title = "Integrated Samples")
    
    p2 <- DimPlot(seurat_integrated, reduction = umap_name, 
                  group.by = "condition", raster = TRUE) +
      labs(title = "Integrated Conditions")
    
    integration_comparison <- p1 | p2
  })
  
} else {
  # For batch regression, show sample and condition distribution
  p1 <- DimPlot(seurat_integrated, reduction = umap_name, 
                group.by = "sample_id", raster = TRUE) +
    labs(title = "Batch-Corrected Samples")
  
  p2 <- DimPlot(seurat_integrated, reduction = umap_name, 
                group.by = "condition", raster = TRUE) +
    labs(title = "Batch-Corrected Conditions") 
  
  integration_comparison <- p1 | p2
}

# Save integration comparison
ggsave(file.path(output_dir, "plots", "integration_comparison.png"), 
       integration_comparison, width = 12, height = 10, dpi = 300)

# Post-integration overview plots
print("Creating post-integration overview plots...")

tryCatch({
  # Main overview plots
  p3 <- DimPlot(seurat_integrated, reduction = umap_name, 
                group.by = "condition", raster = TRUE) +
    labs(title = "Conditions After Integration")
  
  p4 <- DimPlot(seurat_integrated, reduction = umap_name, 
                group.by = "seurat_clusters", 
                label = TRUE, repel = TRUE, raster = TRUE) +
    labs(title = paste("Clusters (res =", default_resolution, ")")) +
    NoLegend()
  
  p5 <- DimPlot(seurat_integrated, reduction = umap_name, 
                group.by = "sample_id", raster = TRUE) +
    labs(title = "Samples After Integration")
  
  # QC metrics overlay
  p6 <- FeaturePlot(seurat_integrated, reduction = umap_name,
                    features = "nFeature_RNA", raster = TRUE) +
    labs(title = "Genes per Cell")
  
  # Combined overview
  overview_plots <- (p3 | p4) / (p5 | p6)
  ggsave(file.path(output_dir, "plots", "integration_overview.png"),
         overview_plots, width = 16, height = 12, dpi = 300)
  
  # Individual high-quality plots
  ggsave(file.path(output_dir, "plots", "UMAP_by_condition.png"), p3, width = 10, height = 8, dpi = 300)
  ggsave(file.path(output_dir, "plots", "UMAP_by_clusters.png"), p4, width = 10, height = 8, dpi = 300)
  ggsave(file.path(output_dir, "plots", "UMAP_by_sample.png"), p5, width = 10, height = 8, dpi = 300)
  
}, error = function(e) {
  warning("Could not create some visualization plots: ", e$message)
})

# Clustering resolution comparison (if multiple resolutions tested)
if (length(clustering_resolutions) > 1) {
  tryCatch({
    res_plots <- list()
    for(i in 1:min(4, length(clustering_resolutions))) {
      res <- clustering_resolutions[i]
      res_col <- paste0("RNA_snn_res.", res)
      if (res_col %in% colnames(seurat_integrated@meta.data)) {
        res_plots[[i]] <- DimPlot(seurat_integrated, reduction = umap_name,
                                 group.by = res_col, label = TRUE, raster = TRUE) +
          labs(title = paste("Resolution", res)) + NoLegend()
      }
    }
    
    if (length(res_plots) > 0) {
      res_comparison <- wrap_plots(res_plots, ncol = 2)
      ggsave(file.path(output_dir, "plots", "clustering_resolutions_comparison.png"),
             res_comparison, width = 12, height = 10, dpi = 300)
    }
  }, error = function(e) {
    warning("Could not create resolution comparison: ", e$message)
  })
}

# =============================================================================
# 8. SAVE RESULTS AND CREATE SUMMARY
# =============================================================================

print("=== STEP 8: SAVING RESULTS ===")

# Save integrated object with descriptive name
output_filename <- switch(integration_method,
                         "Harmony" = "harmony_integrated_seurat.rds",
                         "FastMNN" = "fastmnn_integrated_seurat.rds",
                         "seurat_integrated.rds")

tryCatch({
  saveRDS(seurat_integrated, file.path(output_dir, output_filename))
  print(paste("✓ Integrated object saved:", output_filename))
}, error = function(e) {
  stop("Could not save integrated object: ", e$message)
})

# Save integration metrics
saveRDS(integration_metrics, file.path(output_dir, "metrics", "integration_metrics.rds"))

# Create detailed integration summary
integration_summary <- list(
  # Method and parameters
  integration_method = integration_method,
  reduction_used = reduction_to_use,
  dimensions_used = integration_dims_to_use,
  pcs_calculated = pcs_to_use,
  
  # Data characteristics
  input_cells = ncol(seurat_obj),
  output_cells = ncol(seurat_integrated),
  input_genes = nrow(seurat_obj),
  variable_features = length(VariableFeatures(seurat_integrated)),
  
  # Clustering results
  clustering_resolutions_tested = clustering_resolutions,
  default_resolution = default_resolution,
  final_clusters = length(unique(seurat_integrated$seurat_clusters)),
  
  # Sample and condition info
  samples = unique(seurat_integrated$sample_id),
  conditions = unique(seurat_integrated$condition),
  cells_per_sample = as.list(table(seurat_integrated$sample_id)),
  cells_per_condition = as.list(table(seurat_integrated$condition)),
  
  # Integration quality
  mean_sample_mixing_entropy = integration_metrics$mean_sample_entropy,
  mean_condition_mixing_entropy = integration_metrics$mean_condition_entropy,
  
  # Processing info
  processing_date = Sys.Date(),
  processing_time = as.character(Sys.time()),
  configuration_used = config$dataset_info$optimization %||% "fallback"
)

saveRDS(integration_summary, file.path(output_dir, "integration_summary.rds"))

# Create integration report
create_integration_report <- function(output_dir, summary, method) {
  tryCatch({
    report_lines <- c(
      "# Single-Cell RNA-seq Integration Report",
      paste("Generated on:", Sys.Date()),
      paste("Processing time:", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 2), "minutes"),
      "",
      "## Integration Method",
      paste("- Method used:", method),
      paste("- Reduction:", summary$reduction_used),
      paste("- Dimensions used:", summary$dimensions_used),
      paste("- Configuration:", summary$configuration_used),
      "",
      "## Dataset Information", 
      paste("- Input cells:", summary$input_cells),
      paste("- Output cells:", summary$output_cells),
      paste("- Variable features:", summary$variable_features),
      paste("- Samples:", length(summary$samples)),
      paste("- Conditions:", length(summary$conditions)),
      "",
      "## Clustering Results",
      paste("- Resolutions tested:", paste(summary$clustering_resolutions_tested, collapse = ", ")),
      paste("- Default resolution:", summary$default_resolution),
      paste("- Final clusters:", summary$final_clusters),
      "",
      "## Integration Quality",
      paste("- Sample mixing entropy:", round(summary$mean_sample_mixing_entropy, 3)),
      paste("- Condition mixing entropy:", round(summary$mean_condition_mixing_entropy, 3)),
      "",
      "## Output Files",
      "### Main Results",
      paste("- ", output_filename, " - Integrated Seurat object"),
      "",
      "### Visualizations",
      "- integration_comparison.png - Before/after integration",
      "- integration_overview.png - Post-integration overview",
      "- UMAP_by_condition.png - Conditions on UMAP",
      "- UMAP_by_clusters.png - Clusters on UMAP", 
      "- UMAP_by_sample.png - Samples on UMAP",
      "",
      "### Metrics and Data",
      "- integration_summary.rds - Complete integration summary",
      "- integration_metrics.rds - Quality metrics",
      "- sample_mixing_by_cluster.csv - Sample distribution per cluster",
      "- condition_mixing_by_cluster.csv - Condition distribution per cluster",
      "",
      "## Next Steps",
      "1. Run 3_CellType_Annotation.R for cell type identification",
      "2. Perform downstream analysis on integrated data",
      "3. Validate integration quality with known markers"
    )
    
    writeLines(report_lines, file.path(output_dir, "Integration_Report.md"))
    return(TRUE)
  }, error = function(e) {
    warning("Could not create integration report: ", e$message)
    return(FALSE)
  })
}

report_created <- create_integration_report(output_dir, integration_summary, integration_method)

# Final summary
end_time <- Sys.time()
processing_time <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 2)

cat("\n=== INTEGRATION COMPLETED SUCCESSFULLY ===\n")
cat("Method used:", integration_method, "\n")
cat("Total processing time:", processing_time, "minutes\n")
cat("Input cells:", integration_summary$input_cells, "\n")
cat("Output cells:", integration_summary$output_cells, "\n")
cat("Final clusters:", integration_summary$final_clusters, "\n")
cat("Integration quality (sample mixing):", round(integration_summary$mean_sample_mixing_entropy, 3), "\n")
cat("Integrated object saved as:", output_filename, "\n")
cat("All outputs saved in:", output_dir, "\n")
if (report_created) {
  cat("Integration report: Integration_Report.md\n")
}

print("✅ Hybrid integration pipeline completed successfully!")
print("🚀 Ready for cell type annotation!")

# Clean up
gc(verbose = FALSE, reset = TRUE)

cat("\n=== INTEGRATION PIPELINE SUMMARY ===\n")
cat("Status: SUCCESS\n")
cat("Method:", integration_method, "\n") 
cat("Next step: Run 3_CellType_Annotation.R\n")
cat("==========================================\n")
