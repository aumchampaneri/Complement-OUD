# GSE207128 Analysis - Harmony Integration
# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(harmony)

# Set paths
qc_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE207128/qc_results"
output_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE207128/harmony_results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Get list of all QC processed RDS files
rds_files <- list.files(qc_dir, pattern = ".*_QC_processed.rds$", full.names = TRUE)
sample_ids <- gsub("_QC_processed.rds$", "", basename(rds_files))

cat("Found", length(rds_files), "QC-processed Seurat objects\n")

# Load QC-processed objects
seurat_list <- list()
for (i in 1:length(rds_files)) {
  sample_id <- sample_ids[i]
  cat("Loading", sample_id, "\n")

  # Load object with error handling
  tryCatch({
    seurat_obj <- readRDS(rds_files[i])
    seurat_list[[sample_id]] <- seurat_obj
  }, error = function(e) {
    cat("Error loading", sample_id, ":", e$message, "\n")
  })
}

# Perform dataset integration using Harmony
cat("Performing dataset integration using Harmony...\n")

# Merge all Seurat objects
seurat_merged <- merge(seurat_list[[1]],
                      y = seurat_list[2:length(seurat_list)],
                      add.cell.ids = names(seurat_list),
                      project = "GSE207128")

# Save the merged object before Harmony
saveRDS(seurat_merged, file.path(output_dir, "seurat_merged.rds"))

# Standard preprocessing on merged object
DefaultAssay(seurat_merged) <- "RNA"
seurat_merged <- NormalizeData(seurat_merged)
seurat_merged <- FindVariableFeatures(seurat_merged, selection.method = "vst", nfeatures = 2000)
seurat_merged <- ScaleData(seurat_merged, verbose = FALSE)
seurat_merged <- RunPCA(seurat_merged, npcs = 30, verbose = FALSE)

# Save elbow plot
p_elbow <- ElbowPlot(seurat_merged)
ggsave(file.path(output_dir, "merged_elbow_plot.png"), p_elbow, width = 8, height = 5)

# Run Harmony for batch correction
cat("Running Harmony integration...\n")
seurat_harmony <- RunHarmony(seurat_merged,
                           group.by.vars = "sample",
                           reduction = "pca",
                           project.dim = FALSE)

# Run dimension reduction using Harmony embeddings
cat("Running UMAP and clustering...\n")
seurat_harmony <- RunUMAP(seurat_harmony, reduction = "harmony", dims = 1:30)
seurat_harmony <- FindNeighbors(seurat_harmony, reduction = "harmony", dims = 1:30)
seurat_harmony <- FindClusters(seurat_harmony, resolution = 0.5)

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

# Save top 10 markers per cluster
top10 <- harmony_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10, file.path(output_dir, "harmony_top10_markers_per_cluster.csv"), row.names = FALSE)

# Save integrated object
saveRDS(seurat_harmony, file.path(output_dir, "seurat_harmony_integrated.rds"))

# Save session info for reproducibility
writeLines(capture.output(sessionInfo()), file.path(output_dir, "harmony_sessionInfo.txt"))

cat("Harmony integration complete! Results saved to", output_dir, "\n")