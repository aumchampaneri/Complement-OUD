library(Seurat)
library(ggplot2)
library(pheatmap)
library(patchwork)

# Install missing package with specified repository
if (!requireNamespace("networkD3", quietly = TRUE)) {
  install.packages("networkD3", repos = "https://cloud.r-project.org")
  library(networkD3)
}

# Load object with both annotations
seurat_dual_annotated <- readRDS("/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE207128/harmony_results/annotation_validation/seurat_dual_annotated.rds")

# Create output directory
merged_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE207128/harmony_results/merged_annotations"
dir.create(merged_dir, showWarnings = FALSE, recursive = TRUE)

# Print the first few cell type labels from each method to examine
cat("First few marker-based annotations:", head(unique(seurat_dual_annotated$predicted_celltype)), "\n")
cat("First few SingleR annotations:", head(unique(seurat_dual_annotated$singler_label)), "\n")

# Create consensus table
consensus_table <- table(
  Marker = seurat_dual_annotated$predicted_celltype,
  SingleR = seurat_dual_annotated$singler_label
)

# Save the raw correspondence table for manual review
write.csv(consensus_table, file.path(merged_dir, "annotation_raw_correspondence.csv"))

# Calculate agreement score
agreement_score <- sum(diag(consensus_table)) / sum(consensus_table)
cat("Overall agreement between methods:", round(agreement_score * 100, 2), "%\n")

# Create simple direct match variable
seurat_dual_annotated$direct_match <- seurat_dual_annotated$predicted_celltype == seurat_dual_annotated$singler_label

# Create integrated annotation (with priority to marker-based)
seurat_dual_annotated$integrated_celltype <- seurat_dual_annotated$predicted_celltype

# Calculate stats on direct agreement
celltype_confidence <- table(seurat_dual_annotated$predicted_celltype,
                            seurat_dual_annotated$direct_match)
confidence_stats <- data.frame(
  CellType = rownames(celltype_confidence),
  Total = rowSums(celltype_confidence),
  Matching = celltype_confidence[, "TRUE"],
  Agreement = celltype_confidence[, "TRUE"] / rowSums(celltype_confidence)
)
confidence_stats <- confidence_stats[order(-confidence_stats$Agreement), ]
write.csv(confidence_stats, file.path(merged_dir, "celltype_confidence_stats.csv"))

# Create confidence heatmap
pheatmap(consensus_table,
         filename = file.path(merged_dir, "annotation_correspondence_heatmap.png"),
         main = "Cell Type Correspondence Between Methods",
         fontsize_row = 8,
         fontsize_col = 8)

# Create individual plots with better dimensions
p1 <- DimPlot(seurat_dual_annotated, reduction = "umap", group.by = "predicted_celltype",
             label = TRUE, repel = TRUE) + ggtitle("Marker-based")
ggsave(file.path(merged_dir, "1_marker_based_annotation.png"), p1, width = 12, height = 10)

p2 <- DimPlot(seurat_dual_annotated, reduction = "umap", group.by = "singler_label",
             label = TRUE, repel = TRUE) + ggtitle("SingleR")
ggsave(file.path(merged_dir, "2_singler_annotation.png"), p2, width = 12, height = 10)

p3 <- DimPlot(seurat_dual_annotated, reduction = "umap", group.by = "integrated_celltype",
             label = TRUE, repel = TRUE) + ggtitle("Integrated Annotation")
ggsave(file.path(merged_dir, "3_integrated_annotation.png"), p3, width = 12, height = 10)

p4 <- DimPlot(seurat_dual_annotated, reduction = "umap", group.by = "direct_match",
             cols = c("red", "blue")) + ggtitle("Annotation Agreement")
ggsave(file.path(merged_dir, "4_annotation_agreement.png"), p4, width = 10, height = 8)

# Save the final annotated object
seurat_dual_annotated$annotation_version <- "v1.0_marker+singler"
saveRDS(seurat_dual_annotated, file.path(merged_dir, "seurat_final_annotated.rds"))

cat("Integrated annotation complete! Results saved to", merged_dir, "\n")
cat("NOTE: With only 1.45% agreement, manually review the cell type labels from both methods.\n")
cat("You may need to create a custom mapping between the two annotation schemes.\n")