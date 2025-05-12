# GSE207128 Analysis - Harmony Integration
# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(harmony)

# Set seed for reproducibility
set.seed(2025)

# Set paths
qc_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE207128/qc_results"
output_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE207128/harmony_results"

# Detect Seurat version - with proper conversion
tryCatch({
  seurat_version <- as.character(packageVersion("Seurat"))
  cat("Using Seurat version:", seurat_version, "\n")
  is_seurat_v5 <- numeric_version(seurat_version) >= numeric_version("5.0.0")
}, error = function(e) {
  cat("Could not determine Seurat version:", e$message, "\n")
  is_seurat_v5 <- FALSE  # Default to Seurat v4 behavior
})

# Check if output directory is writable
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(output_dir)) {
  stop("Cannot create directory: ", output_dir)
}

# Test write permission with a small file
test_file <- file.path(output_dir, "write_test.txt")
tryCatch({
  writeLines("test", test_file)
  file.remove(test_file)
}, error = function(e) {
  stop("Cannot write to directory: ", output_dir, " (", e$message, ")")
})

# Check if merged file exists to save processing time
merged_file_path <- file.path(output_dir, "seurat_merged.rds")

if (file.exists(merged_file_path)) {
  cat("Loading existing merged Seurat object...\n")
  seurat_merged <- readRDS(merged_file_path)

  # Check if we need to fix sample identifiers for Harmony
  if (!"sample" %in% colnames(seurat_merged@meta.data)) {
    cat("Creating 'sample' column from 'orig.ident'...\n")
    seurat_merged$sample <- seurat_merged$orig.ident
  }

  # Check if PCA exists, and run it if not
  if (!"pca" %in% Reductions(seurat_merged)) {
    cat("PCA reduction not found. Running preprocessing...\n")
    DefaultAssay(seurat_merged) <- "RNA"

    # Check if variable features are defined
    if (length(VariableFeatures(seurat_merged)) == 0) {
      cat("Finding variable features...\n")
      seurat_merged <- FindVariableFeatures(seurat_merged, selection.method = "vst", nfeatures = 2000)
      gc()
    }

    # Scale data using the appropriate method for the Seurat version
    cat("Scaling data...\n")
    seurat_merged <- ScaleData(seurat_merged, verbose = FALSE)
    gc()

    cat("Running PCA...\n")
    seurat_merged <- RunPCA(seurat_merged, npcs = 30, verbose = FALSE)
    gc()

    # Save the updated merged object
    saveRDS(seurat_merged, merged_file_path)

    # Save elbow plot
    p_elbow <- ElbowPlot(seurat_merged)
    ggsave(file.path(output_dir, "merged_elbow_plot.png"), p_elbow, width = 8, height = 5)
  }
} else {
  # Get list of all QC processed RDS files
  rds_files <- list.files(qc_dir, pattern = ".*_QC_processed.rds$", full.names = TRUE)

  if (length(rds_files) == 0) {
    stop("No QC-processed RDS files found in: ", qc_dir)
  }

  sample_ids <- gsub("_QC_processed.rds$", "", basename(rds_files))
  cat("Found", length(rds_files), "QC-processed Seurat objects\n")

  # Load QC-processed objects
  seurat_list <- list()
  for (i in 1:length(rds_files)) {
    sample_id <- sample_ids[i]
    cat("Loading", sample_id, "\n")

    tryCatch({
      seurat_obj <- readRDS(rds_files[i])
      seurat_list[[sample_id]] <- seurat_obj
    }, error = function(e) {
      cat("Error loading", sample_id, ":", e$message, "\n")
    })
  }

  # Check if we have valid objects to merge
  if (length(seurat_list) == 0) {
    stop("No valid Seurat objects loaded")
  }

  cat("Performing dataset integration using Harmony...\n")

  # Handle single sample case
  if (length(seurat_list) == 1) {
    cat("Only one sample found - no merging needed\n")
    seurat_merged <- seurat_list[[1]]
    seurat_merged$sample <- names(seurat_list)[1]
  } else {
    # Merge all Seurat objects
    seurat_merged <- merge(seurat_list[[1]],
                          y = seurat_list[2:length(seurat_list)],
                          add.cell.ids = names(seurat_list),
                          project = "GSE207128")
    # Create sample column for Harmony
    seurat_merged$sample <- seurat_merged$orig.ident
  }

  # Free memory
  rm(seurat_list)
  gc()

  # Standard preprocessing on merged object
  DefaultAssay(seurat_merged) <- "RNA"
  seurat_merged <- NormalizeData(seurat_merged)
  gc()

  seurat_merged <- FindVariableFeatures(seurat_merged, selection.method = "vst", nfeatures = 2000)
  gc()

  # Scale data
  seurat_merged <- ScaleData(seurat_merged, verbose = FALSE)
  gc()

  seurat_merged <- RunPCA(seurat_merged, npcs = 30, verbose = FALSE)
  gc()

  # Save the merged object for future runs
  saveRDS(seurat_merged, merged_file_path)

  # Save elbow plot
  p_elbow <- ElbowPlot(seurat_merged)
  ggsave(file.path(output_dir, "merged_elbow_plot.png"), p_elbow, width = 8, height = 5)
}

# Print available reductions and metadata for debugging
cat("Available reductions:", paste(Reductions(seurat_merged), collapse = ", "), "\n")
cat("Available metadata columns:", paste(colnames(seurat_merged@meta.data), collapse = ", "), "\n")

# Only run Harmony if we have multiple samples
sample_levels <- unique(seurat_merged$sample)
if (length(sample_levels) <= 1) {
  cat("Only one sample level detected:", sample_levels, "\n")
  cat("Skipping Harmony integration - proceeding with PCA-based analysis\n")
  seurat_harmony <- seurat_merged
  reduction_to_use <- "pca"
} else {
  # Use the direct harmony function
  cat("Running Harmony integration on", length(sample_levels), "samples...\n")

  # Extract PCA embeddings
  pca_embeddings <- Embeddings(seurat_merged, reduction = "pca")

  # Ensure metadata rows match PCA embedding rows
  if (!all(rownames(pca_embeddings) == rownames(seurat_merged@meta.data))) {
    cat("Warning: PCA embeddings and metadata rows don't match. Reordering metadata.\n")
    meta_data <- seurat_merged@meta.data[rownames(pca_embeddings), ]
  } else {
    meta_data <- seurat_merged@meta.data
  }

  # Run harmony directly with enhanced error checking
  harmony_embeddings <- tryCatch({
    harmony::HarmonyMatrix(
      data_mat = pca_embeddings,
      meta_data = meta_data,
      vars_use = "sample",
      do_pca = FALSE
    )
  }, error = function(e) {
    stop("Harmony integration failed: ", e$message, "\n",
         "Check that 'sample' column exists in metadata and matches PCA dimensions.")
  })

  # Create a new dimensional reduction with harmony embeddings
  seurat_harmony <- seurat_merged
  seurat_harmony[["harmony"]] <- CreateDimReducObject(
    embeddings = harmony_embeddings,
    key = "harmony_",
    assay = DefaultAssay(seurat_merged)
  )

  gc()
  reduction_to_use <- "harmony"
}

# Run dimension reduction and clustering
cat("Running UMAP and clustering using", reduction_to_use, "reduction...\n")

seurat_harmony <- RunUMAP(seurat_harmony, reduction = reduction_to_use, dims = 1:30)
gc()

seurat_harmony <- FindNeighbors(seurat_harmony, reduction = reduction_to_use, dims = 1:30)
seurat_harmony <- FindClusters(seurat_harmony, resolution = 0.5)
gc()

# Create UMAP visualizations
p_clusters <- DimPlot(seurat_harmony, reduction = "umap", group.by = "seurat_clusters")
ggsave(file.path(output_dir, "harmony_umap_clusters.png"), p_clusters, width = 10, height = 8)

p_samples <- DimPlot(seurat_harmony, reduction = "umap", group.by = "sample")
ggsave(file.path(output_dir, "harmony_umap_samples.png"), p_samples, width = 10, height = 8)

# Find marker genes for each cluster
cat("Finding marker genes...\n")
DefaultAssay(seurat_harmony) <- "RNA"
harmony_markers <- FindAllMarkers(seurat_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(harmony_markers, file.path(output_dir, "harmony_cluster_markers.csv"), row.names = FALSE)

# Save top 10 markers per cluster (using version-aware approach)
if (packageVersion("dplyr") >= numeric_version("1.0.0")) {
  top10 <- harmony_markers %>% group_by(cluster) %>% slice_max(n = 10, order_by = avg_log2FC)
} else {
  top10 <- harmony_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
}
write.csv(top10, file.path(output_dir, "harmony_top10_markers_per_cluster.csv"), row.names = FALSE)

# Save integrated object
saveRDS(seurat_harmony, file.path(output_dir, "seurat_harmony_integrated.rds"))

# Save session info for reproducibility
writeLines(capture.output(sessionInfo()), file.path(output_dir, "harmony_sessionInfo.txt"))

cat("Analysis complete! Results saved to", output_dir, "\n")