# Script to filter out cardiomyocytes from SingleR annotations

# Load required libraries
library(Seurat)
library(ggplot2)

# Load the SingleR annotated object
seurat_obj <- readRDS("/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE207128/harmony_results/singleR_annotation/seurat_dual_annotated.rds")

# Print original annotation structure
cat("Original SingleR annotation distribution:\n")
print(table(seurat_obj$singler_label))

# Create a new column to preserve original labels
seurat_obj$singler_label_original <- seurat_obj$singler_label

# Count cells to be removed
cells_to_remove <- sum(seurat_obj$singler_label == "Cardiomyocytes")
cat("Removing", cells_to_remove, "Cardiomyocytes cells from the dataset\n")

# Filter out the cardiomyocytes
seurat_obj_filtered <- subset(seurat_obj, subset = singler_label != "Cardiomyocytes")

# Print new cell count
cat("\nCells before filtering:", ncol(seurat_obj), "\n")
cat("Cells after filtering:", ncol(seurat_obj_filtered), "\n")

# Print new annotation distribution
cat("\nUpdated SingleR annotation distribution after removal:\n")
print(table(seurat_obj_filtered$singler_label))

# Check if gene names are present in the data
cat("\nChecking data structure for gene names...\n")
gene_ids <- head(rownames(seurat_obj_filtered), 10)
print("First 10 gene IDs:")
print(gene_ids)

# Check if these are Ensembl IDs or gene symbols
if(any(grepl("^ENSMUS", gene_ids))) {
  cat("Data contains Ensembl IDs (e.g., ENSMUSG...)\n")
} else {
  cat("Data contains gene symbols (not Ensembl IDs)\n")
}

# Fix directory path
output_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE207128/harmony_results/singleR_annotation"

# Save the filtered object
saveRDS(seurat_obj_filtered, file.path(output_dir, "seurat_filtered.rds"))

# Visualize before and after
p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "singler_label",
            label = TRUE, repel = TRUE) + ggtitle("Original SingleR Labels")

p2 <- DimPlot(seurat_obj_filtered, reduction = "umap", group.by = "singler_label",
            label = TRUE, repel = TRUE) + ggtitle("After Removing Cardiomyocytes")

# Save plots
plot_dir <- file.path(output_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE)
ggsave(file.path(plot_dir, "original_labels.png"), p1, width = 12, height = 10)
ggsave(file.path(plot_dir, "filtered_labels.png"), p2, width = 12, height = 10)

cat("Filtering complete! Updated Seurat object saved.\n")