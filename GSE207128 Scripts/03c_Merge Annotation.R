# Create script for merging annotation results
# 05_Merge_Annotations.R

library(Seurat)
library(ggplot2)
library(pheatmap)

# Load object with both annotations
seurat_dual_annotated <- readRDS("/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE207128/harmony_results/annotation_validation/seurat_dual_annotated.rds")

# Create output directory
merged_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE207128/harmony_results/merged_annotations"
dir.create(merged_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Create consensus annotation
consensus_table <- table(
  Marker = seurat_dual_annotated$predicted_celltype,
  SingleR = seurat_dual_annotated$singler_label
)

# Calculate agreement score
agreement_score <- sum(diag(consensus_table)) / sum(consensus_table)
cat("Overall agreement between methods:", round(agreement_score * 100, 2), "%\n")

# Create an annotation confidence score
seurat_dual_annotated$annotation_match <- seurat_dual_annotated$predicted_celltype == seurat_dual_annotated$singler_label

# Create integrated annotation (with confidence label)
seurat_dual_annotated$integrated_celltype <- ifelse(
  seurat_dual_annotated$annotation_match,
  seurat_dual_annotated$predicted_celltype,  # Use consistent name where methods agree
  paste0(seurat_dual_annotated$predicted_celltype, "/", seurat_dual_annotated$singler_label)  # Show both where they differ
)

# Sankey diagram of annotation mapping
# (Optional, requires networkD3 package)
if(require(networkD3)) {
  # Create data for Sankey diagram
  links <- as.data.frame(matrix(0, nrow=nrow(consensus_table), ncol=3))
  colnames(links) <- c("source", "target", "value")

  counter <- 0
  for(i in 1:nrow(consensus_table)) {
    for(j in 1:ncol(consensus_table)) {
      if(consensus_table[i,j] > 0) {
        counter <- counter + 1
        links[counter,] <- c(i-1, j-1+nrow(consensus_table), consensus_table[i,j])
      }
    }
  }
  links <- links[1:counter,]
  links$source <- as.numeric(links$source)
  links$target <- as.numeric(links$target)
  links$value <- as.numeric(links$value)

  # Create nodes
  nodes <- data.frame(
    name = c(rownames(consensus_table), colnames(consensus_table))
  )

  # Create Sankey diagram
  sankey <- sankeyNetwork(
    Links = links, Nodes = nodes,
    Source = "source", Target = "target",
    Value = "value", NodeID = "name",
    sinksRight = FALSE, width = 800, height = 600
  )

  # Save as HTML
  saveNetwork(sankey, file.path(merged_dir, "annotation_comparison_sankey.html"))
}

# 2. Create integrated annotation plots
p1 <- DimPlot(seurat_dual_annotated, reduction = "umap", group.by = "predicted_celltype",
             label = TRUE, repel = TRUE) + ggtitle("Marker-based")
p2 <- DimPlot(seurat_dual_annotated, reduction = "umap", group.by = "singler_label",
             label = TRUE, repel = TRUE) + ggtitle("SingleR")
p3 <- DimPlot(seurat_dual_annotated, reduction = "umap", group.by = "integrated_celltype",
             label = TRUE, repel = TRUE) + ggtitle("Integrated Annotation")
p4 <- DimPlot(seurat_dual_annotated, reduction = "umap", group.by = "annotation_match",
             cols = c("red", "blue")) + ggtitle("Annotation Agreement")

combined <- (p1 + p2) / (p3 + p4)
ggsave(file.path(merged_dir, "integrated_annotation.png"), combined, width = 16, height = 12)

# 3. Calculate stats on annotation agreement
celltype_confidence <- table(seurat_dual_annotated$predicted_celltype,
                            seurat_dual_annotated$annotation_match)
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
         display_numbers = TRUE,
         filename = file.path(merged_dir, "annotation_correspondence_heatmap.png"),
         main = "Cell Type Correspondence Between Methods")

# Save the final annotated object
seurat_dual_annotated$annotation_version <- "v1.0_marker+singler"
saveRDS(seurat_dual_annotated, file.path(merged_dir, "seurat_final_annotated.rds"))

cat("Integrated annotation complete! Results saved to", merged_dir, "\n")