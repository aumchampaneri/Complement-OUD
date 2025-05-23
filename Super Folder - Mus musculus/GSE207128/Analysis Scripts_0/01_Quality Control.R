# Analysis Outputs_0 Analysis - QC and Integration
# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# Set working directory to processed data
setwd("/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE207128/processed_data")

# Create output directory for QC results
qc_dir <- "../qc_results"
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

# Get list of all RDS files
rds_files <- list.files(pattern = ".*_seurat_object.rds$")
sample_ids <- gsub("_seurat_object.rds$", "", rds_files)

cat("Found", length(rds_files), "Seurat objects\n")

# Create a summary data frame for QC stats
qc_stats <- data.frame(
  sample_id = character(),
  total_cells = integer(),
  retained_cells = integer(),
  retention_rate = numeric(),
  stringsAsFactors = FALSE
)

# Load and QC each object
seurat_list <- list()
for (i in 1:length(rds_files)) {
  sample_id <- sample_ids[i]
  cat("Processing", sample_id, "\n")

  # Load object with error handling
  tryCatch({
    seurat_obj <- readRDS(rds_files[i])

    # Calculate mitochondrial percentage
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

    # Add log-transformed metrics
    seurat_obj$log10_nCount_RNA <- log10(seurat_obj$nCount_RNA + 1)
    seurat_obj$log10_nFeature_RNA <- log10(seurat_obj$nFeature_RNA + 1)

    # QC plots before filtering
    p1 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",
                                          "log10_nFeature_RNA", "log10_nCount_RNA"),
                 ncol = 3)
    ggsave(file.path(qc_dir, paste0(sample_id, "_pre_QC.png")), p1, width = 15, height = 10)

    # Count cells before filtering
    total_cells <- ncol(seurat_obj)

    # Filter cells
    seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 &
                                            nFeature_RNA < 5000 &
                                            percent.mt < 15)

    # Count cells after filtering
    retained_cells <- ncol(seurat_obj)
    retention_rate <- 100 * retained_cells / total_cells

    # Log filtering statistics
    cat(sprintf("Retained %d of %d cells (%.1f%%) for %s\n",
                retained_cells, total_cells, retention_rate, sample_id))

    # Add to QC stats dataframe
    qc_stats <- rbind(qc_stats, data.frame(
      sample_id = sample_id,
      total_cells = total_cells,
      retained_cells = retained_cells,
      retention_rate = retention_rate
    ))

    # QC plots after filtering
    p2 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",
                                          "log10_nFeature_RNA", "log10_nCount_RNA"),
                 ncol = 3)
    ggsave(file.path(qc_dir, paste0(sample_id, "_post_QC.png")), p2, width = 15, height = 10)

    # Normalize and find variable features
    seurat_obj <- NormalizeData(seurat_obj)
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

    # Add sample metadata
    seurat_obj$sample <- sample_id
    seurat_obj$orig.ident <- sample_id

    # Store in list
    seurat_list[[sample_id]] <- seurat_obj

    # Save processed object
    saveRDS(seurat_obj, file.path(qc_dir, paste0(sample_id, "_QC_processed.rds")))

    # Clean up to save memory
    rm(seurat_obj, p1, p2)
    gc()

  }, error = function(e) {
    cat("Error processing", sample_id, ":", e$message, "\n")
  })
}

# Save QC stats
write.csv(qc_stats, file.path(qc_dir, "qc_filtering_stats.csv"), row.names = FALSE)



##### Here on down failed on computer due to memory issues -> use 02_Harmony.R instead

# Integrate datasets (using Seurat v4 approach)
cat("Performing dataset integration...\n")

# Select features for integration
features <- SelectIntegrationFeatures(object.list = seurat_list)

# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = seurat_list,
                                 anchor.features = features,
                                 normalization.method = "LogNormalize")

# Create integrated object
seurat_integrated <- IntegrateData(anchorset = anchors)

# Switch to integrated assay for downstream analysis
DefaultAssay(seurat_integrated) <- "integrated"

# Run dimension reduction
seurat_integrated <- ScaleData(seurat_integrated)
seurat_integrated <- RunPCA(seurat_integrated)

# Create and save elbow plot for PC selection
p_elbow <- ElbowPlot(seurat_integrated)
ggsave(file.path(qc_dir, "integrated_elbow_plot.png"), p_elbow, width = 8, height = 5)

# Continue with dimension reduction
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:30)

# Clustering
seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:30)
seurat_integrated <- FindClusters(seurat_integrated, resolution = 0.5)

# Create UMAP visualization
p3 <- DimPlot(seurat_integrated, reduction = "umap", group.by = "seurat_clusters")
ggsave(file.path(qc_dir, "integrated_umap_clusters.png"), p3, width = 10, height = 8)

p4 <- DimPlot(seurat_integrated, reduction = "umap", group.by = "sample")
ggsave(file.path(qc_dir, "integrated_umap_samples.png"), p4, width = 10, height = 8)

# Find marker genes for each cluster
DefaultAssay(seurat_integrated) <- "RNA"
markers <- FindAllMarkers(seurat_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, file.path(qc_dir, "cluster_markers.csv"), row.names = FALSE)

# Save top 10 markers per cluster
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10, file.path(qc_dir, "top10_markers_per_cluster.csv"), row.names = FALSE)

# Save integrated object
saveRDS(seurat_integrated, file.path(qc_dir, "seurat_integrated.rds"))

# Save session info for reproducibility
writeLines(capture.output(sessionInfo()), file.path(qc_dir, "sessionInfo.txt"))

cat("Analysis complete! Integrated object saved.\n")