# 03a_Manual Annotation.R - Fixed version

# Load required packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(biomaRt)

# Define paths
harmony_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE207128/harmony_results"
output_dir <- file.path(harmony_dir, "cell_annotation")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define brain-specific cell markers for mouse
brain_markers <- list(
  "Neurons" = c("Rbfox3", "Snap25", "Syt1", "Tubb3"),
  "Excitatory_neurons" = c("Slc17a7", "Camk2a", "Nrgn", "Satb2"),
  "Inhibitory_neurons" = c("Gad1", "Gad2", "Slc32a1", "Dlx1"),
  "Astrocytes" = c("Gfap", "Aldh1l1", "Aqp4", "Slc1a3"),
  "Oligodendrocytes" = c("Mbp", "Mag", "Mog", "Plp1"),
  "OPC" = c("Pdgfra", "Cspg4", "Sox10", "Olig1"),
  "Microglia" = c("Cx3cr1", "P2ry12", "Tmem119", "Csf1r", "C1qa"),
  "Endothelial" = c("Cldn5", "Pecam1", "Flt1", "Tek"),
  "Pericytes" = c("Pdgfrb", "Anpep", "Rgs5", "Notch3")
)

# Load Harmony integrated Seurat object
seurat_harmony <- readRDS(file.path(harmony_dir, "seurat_harmony_integrated.rds"))

# Check if gene IDs are in Ensembl format
gene_ids <- head(rownames(seurat_harmony), 5)
cat("Gene ID format in dataset:", gene_ids, "\n")

if(any(grepl("^ENSMUS", gene_ids))) {
  cat("Dataset uses Ensembl IDs. Converting markers from gene symbols...\n")
  # Set up biomaRt connection
  mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

  # Convert brain_markers list to use Ensembl IDs
  all_markers <- unique(unlist(brain_markers))
  mapping <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters = "external_gene_name",
    values = all_markers,
    mart = mart
  )

  brain_markers_ensembl <- lapply(brain_markers, function(genes) {
    ensembl_ids <- mapping$ensembl_gene_id[mapping$external_gene_name %in% genes]
    cat("Found", length(ensembl_ids), "Ensembl IDs for markers\n")
    return(ensembl_ids)
  })

  brain_markers <- brain_markers_ensembl
}

# Add module scores for each cell type
cat("Adding module scores for cell types...\n")
for(cell_type in names(brain_markers)) {
  markers <- brain_markers[[cell_type]]
  if(length(markers) > 0) {
    seurat_harmony <- AddModuleScore(
      seurat_harmony,
      features = list(markers),
      name = cell_type,
      seed = 123
    )
    # Rename score column
    old_col <- paste0(cell_type, "1")
    new_col <- paste0(cell_type, "_score")
    seurat_harmony[[new_col]] <- seurat_harmony[[old_col]]
    seurat_harmony[[old_col]] <- NULL
  }
}

# Get score columns
score_cols <- grep("_score$", colnames(seurat_harmony@meta.data), value=TRUE)
cat("Created", length(score_cols), "cell type scores\n")

# Assign cell types based on highest score
scores_df <- as.matrix(seurat_harmony@meta.data[, score_cols])
top_indices <- max.col(scores_df, ties.method = "first")
cell_types <- gsub("_score$", "", score_cols)
seurat_harmony$predicted_celltype <- cell_types[top_indices]

# Calculate score gap (confidence)
get_score_gap <- function(row) {
  sorted <- sort(row, decreasing = TRUE)
  if(length(sorted) > 1) return(sorted[1] - sorted[2])
  return(1)
}
seurat_harmony$score_gap <- apply(scores_df, 1, get_score_gap)

# Mark uncertain cells
threshold <- 0.1
seurat_harmony$predicted_celltype[seurat_harmony$score_gap < threshold] <- "Uncertain"

# Create summary tables
cluster_type_table <- table(seurat_harmony$seurat_clusters, seurat_harmony$predicted_celltype)
cluster_type_pct <- prop.table(cluster_type_table, margin = 1) * 100
write.csv(cluster_type_table, file.path(output_dir, "cluster_celltype_counts.csv"))
write.csv(cluster_type_pct, file.path(output_dir, "cluster_celltype_percent.csv"))

# Visualizations - Key marker genes
cat("Creating marker visualizations...\n")

# Check which markers are present in the dataset
neuron_markers <- c("Snap25", "Gad1", "Slc17a7")
glial_markers <- c("Aldh1l1", "Gfap", "Csf1r", "C1qa", "Mbp")

# Convert to Ensembl IDs if needed
if(any(grepl("^ENSMUS", gene_ids))) {
  neuron_mapping <- mapping[mapping$external_gene_name %in% neuron_markers, ]
  glial_mapping <- mapping[mapping$external_gene_name %in% glial_markers, ]

  neuron_markers <- neuron_mapping$ensembl_gene_id
  glial_markers <- glial_mapping$ensembl_gene_id
}

# Filter to markers present in dataset
present_neuron_markers <- intersect(neuron_markers, rownames(seurat_harmony))
present_glial_markers <- intersect(glial_markers, rownames(seurat_harmony))

# Create feature plots
if(length(present_neuron_markers) > 0) {
  p1 <- FeaturePlot(seurat_harmony, features = present_neuron_markers,
                  reduction = "umap", ncol = length(present_neuron_markers))
  ggsave(file.path(output_dir, "neuron_markers.png"), p1, width = min(15, 5*length(present_neuron_markers)), height = 5)
}

if(length(present_glial_markers) > 0) {
  p2 <- FeaturePlot(seurat_harmony, features = present_glial_markers,
                  reduction = "umap", ncol = min(3, length(present_glial_markers)))
  ggsave(file.path(output_dir, "glial_markers.png"), p2, width = min(15, 5*length(present_glial_markers)),
         height = ceiling(length(present_glial_markers)/3) * 5)
}

# Create dot plot summary - Note: RotatedAxis is from Seurat, not ggplot2
all_markers <- c(present_neuron_markers, present_glial_markers)
if(length(all_markers) > 0) {
  tryCatch({
    p3 <- DotPlot(seurat_harmony, features = all_markers, group.by = "predicted_celltype") +
         Seurat::RotatedAxis()
    ggsave(file.path(output_dir, "marker_dotplot.png"), p3, width = min(15, length(all_markers)), height = 8)
  }, error = function(e) {
    cat("Error creating dotplot:", e$message, "\n")
  })
}

# Generate suggested cluster renaming
cat("Creating cluster to cell type mapping...\n")
dominant_type <- apply(cluster_type_pct, 1, function(row) {
  if(max(row) < 50) return("Mixed")
  return(names(which.max(row)))
})

# Create a mapping dataframe
cluster_mapping <- data.frame(
  Cluster = rownames(cluster_type_pct),
  DominantType = dominant_type,
  PercentDominant = apply(cluster_type_pct, 1, max)
)
write.csv(cluster_mapping, file.path(output_dir, "cluster_to_celltype_mapping.csv"))

# Create manual renaming function with proper handling of duplicates
rename_script <- "# Function to rename clusters based on dominant cell types\n"
rename_script <- paste0(rename_script,
"rename_clusters <- function(seurat_obj) {
  # Load cluster mapping
  mapping <- read.csv('", file.path(output_dir, "cluster_to_celltype_mapping.csv"), "')

  # Create named vector for renaming
  new_ids <- mapping$DominantType
  names(new_ids) <- mapping$Cluster

  # Rename identities
  Idents(seurat_obj) <- 'seurat_clusters'
  seurat_obj <- RenameIdents(seurat_obj, new_ids)

  # Save to metadata
  seurat_obj$cluster_celltype <- Idents(seurat_obj)

  return(seurat_obj)
}
")

writeLines(rename_script, file.path(output_dir, "rename_clusters.R"))

# Create cell type UMAP visualization
p4 <- DimPlot(seurat_harmony, reduction = "umap", group.by = "predicted_celltype",
            label = TRUE, repel = TRUE) + ggtitle("Predicted Cell Types")
ggsave(file.path(output_dir, "celltype_umap.png"), p4, width = 12, height = 10)

# Save the annotated object
cat("Saving annotated Seurat object...\n")
saveRDS(seurat_harmony, file.path(output_dir, "seurat_annotated.rds"))

# Manual cluster renaming
cat("Applying cluster renaming based on dominant cell types...\n")
source(file.path(output_dir, "rename_clusters.R"))
seurat_harmony <- rename_clusters(seurat_harmony)

# Visualize renamed clusters
p5 <- DimPlot(seurat_harmony, reduction = "umap", group.by = "cluster_celltype",
            label = TRUE, repel = TRUE) + ggtitle("Renamed Clusters")
ggsave(file.path(output_dir, "renamed_clusters_umap.png"), p5, width = 12, height = 10)

# Save final object with renamed clusters
saveRDS(seurat_harmony, file.path(output_dir, "seurat_annotated_final.rds"))

cat("Manual annotation complete! Results saved to:", output_dir, "\n")