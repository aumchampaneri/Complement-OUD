# ==============================================================================
# Complement Gene Expression Heatmap
# GSE207128 - Complement-OUD Analysis
# ==============================================================================

# Load required libraries
library(Seurat)
library(ComplexHeatmap)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(grid)

# Set up paths
base_dir <- "/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus/GSE207128"
data_dir <- file.path(base_dir, "Outputs/04_Annotated_Data")
output_dir <- file.path(base_dir, "Outputs/06_Visualizations")

# Create output directory
dir.create(file.path(output_dir, "Complement_Analysis"), recursive = TRUE, showWarnings = FALSE)

cat("Loading annotated Seurat object...\n")
seurat_obj <- readRDS(file.path(data_dir, "filtered_annotated_seurat.rds"))

# ==============================================================================
# Define Comprehensive Complement Gene Set
# ==============================================================================

complement_genes <- c(
  # Classical pathway
  "C1qa", "C1qb", "C1qc", "C1r", "C1s", "C4a", "C4b", "C2",
  # Alternative pathway  
  "C3", "Cfb", "Cfd", "C5", "C6", "C7", "C8a", "C8b", "C8g", "C9",
  # Lectin pathway
  "Mbl1", "Mbl2", "Masp1", "Masp2", "Fcn1", "Fcn2", "Fcn3",
  # Regulators
  "Cfh", "Cfi", "C4bp", "Cd55", "Cd46", "Cd35", "Cr1", "Cr2",
  # Receptors
  "C3ar1", "C5ar1", "C5ar2", "Itgam", "Itgax", "Itgb2",
  # Additional complement components
  "Cfp", "Clu", "Vtn", "Serping1", "Cfhr1", "Cfhr2"
)

# Check which complement genes are available in the dataset
available_complement <- complement_genes[complement_genes %in% rownames(seurat_obj)]
cat("Available complement genes:", length(available_complement), "out of", length(complement_genes), "\n")
cat("Missing genes:", paste(setdiff(complement_genes, available_complement), collapse = ", "), "\n")

# ==============================================================================
# Create Pseudobulk Expression Matrix
# ==============================================================================

cat("Creating pseudobulk expression matrix...\n")

# Get metadata
metadata <- seurat_obj@meta.data

# Create pseudobulk by averaging expression within each cell type-condition combination
create_pseudobulk <- function(seurat_obj, group_vars) {
  # Normalize data if not already done
  if(!"data" %in% names(seurat_obj@assays$RNA@layers)) {
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  }
  
  # Get expression matrix
  expr_matrix <- GetAssayData(seurat_obj, slot = "data")
  
  # Create grouping variable
  metadata$group <- paste(metadata[[group_vars[1]]], metadata[[group_vars[2]]], sep = "_")
  
  # Calculate mean expression for each group
  pseudobulk <- sapply(unique(metadata$group), function(group) {
    cells_in_group <- rownames(metadata)[metadata$group == group]
    if(length(cells_in_group) > 1) {
      rowMeans(expr_matrix[, cells_in_group, drop = FALSE])
    } else if(length(cells_in_group) == 1) {
      expr_matrix[, cells_in_group]
    } else {
      rep(0, nrow(expr_matrix))
    }
  })
  
  return(pseudobulk)
}

# Create pseudobulk matrix
pseudobulk_matrix <- create_pseudobulk(seurat_obj, c("predicted_celltype", "condition"))

# Filter for complement genes
complement_matrix <- pseudobulk_matrix[available_complement, , drop = FALSE]

# Remove groups with all zeros
complement_matrix <- complement_matrix[, colSums(complement_matrix) > 0, drop = FALSE]

cat("Final matrix dimensions:", nrow(complement_matrix), "genes x", ncol(complement_matrix), "cell type-condition combinations\n")

# ==============================================================================
# Prepare Column Annotations
# ==============================================================================

cat("Preparing column annotations...\n")

# Extract cell type and condition from column names
col_info <- data.frame(
  sample = colnames(complement_matrix),
  stringsAsFactors = FALSE
)

col_info$cell_type <- sapply(strsplit(col_info$sample, "_"), function(x) {
  paste(x[1:(length(x)-1)], collapse = "_")
})

col_info$condition <- sapply(strsplit(col_info$sample, "_"), function(x) {
  x[length(x)]
})

# Define colors
condition_colors <- c("Naive" = "#2E8B57", "Dependent" = "#CD853F", "Withdrawal" = "#B22222")
cell_type_colors <- rainbow(length(unique(col_info$cell_type)))
names(cell_type_colors) <- unique(col_info$cell_type)

# Create column annotation
col_annotation <- HeatmapAnnotation(
  Condition = col_info$condition,
  `Cell Type` = col_info$cell_type,
  col = list(
    Condition = condition_colors,
    `Cell Type` = cell_type_colors
  ),
  annotation_height = unit(c(0.5, 1), "cm"),
  annotation_name_gp = gpar(fontsize = 10),
  annotation_legend_param = list(
    Condition = list(title = "Condition", title_gp = gpar(fontsize = 12, fontface = "bold")),
    `Cell Type` = list(title = "Cell Type", title_gp = gpar(fontsize = 12, fontface = "bold"))
  )
)

# ==============================================================================
# Prepare Row Annotations (Complement Pathways)
# ==============================================================================

cat("Preparing row annotations...\n")

# Define complement pathway categories
complement_pathways <- list(
  "Classical" = c("C1qa", "C1qb", "C1qc", "C1r", "C1s", "C4a", "C4b", "C2"),
  "Alternative" = c("C3", "Cfb", "Cfd", "C5", "C6", "C7", "C8a", "C8b", "C8g", "C9"),
  "Lectin" = c("Mbl1", "Mbl2", "Masp1", "Masp2", "Fcn1", "Fcn2", "Fcn3"),
  "Regulators" = c("Cfh", "Cfi", "C4bp", "Cd55", "Cd46", "Cd35", "Cr1", "Cr2"),
  "Receptors" = c("C3ar1", "C5ar1", "C5ar2", "Itgam", "Itgax", "Itgb2"),
  "Other" = c("Cfp", "Clu", "Vtn", "Serping1", "Cfhr1", "Cfhr2")
)

# Create pathway annotation for available genes
gene_pathways <- rep("Other", nrow(complement_matrix))
names(gene_pathways) <- rownames(complement_matrix)

for(pathway in names(complement_pathways)) {
  pathway_genes <- intersect(complement_pathways[[pathway]], rownames(complement_matrix))
  gene_pathways[pathway_genes] <- pathway
}

# Define pathway colors
pathway_colors <- c(
  "Classical" = "#FF6B6B",
  "Alternative" = "#4ECDC4", 
  "Lectin" = "#45B7D1",
  "Regulators" = "#96CEB4",
  "Receptors" = "#FFEAA7",
  "Other" = "#DDA0DD"
)

# Create row annotation
row_annotation <- rowAnnotation(
  Pathway = gene_pathways,
  col = list(Pathway = pathway_colors),
  width = unit(0.5, "cm"),
  annotation_name_gp = gpar(fontsize = 10),
  annotation_legend_param = list(
    Pathway = list(title = "Complement Pathway", title_gp = gpar(fontsize = 12, fontface = "bold"))
  )
)

# ==============================================================================
# Create Heatmap
# ==============================================================================

cat("Creating complement gene heatmap...\n")

# Scale data for better visualization (z-score by row)
scaled_matrix <- t(scale(t(complement_matrix)))
scaled_matrix[is.na(scaled_matrix)] <- 0

# Define color scale
col_fun <- colorRamp2(
  c(-2, -1, 0, 1, 2), 
  c("#053061", "#2166AC", "#F7F7F7", "#D6604D", "#67001F")
)

# Create main heatmap
main_heatmap <- Heatmap(
  scaled_matrix,
  name = "Expression\n(Z-score)",
  
  # Colors
  col = col_fun,
  
  # Column settings
  top_annotation = col_annotation,
  column_title = "Cell Type - Condition Combinations",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 45,
  
  # Row settings
  left_annotation = row_annotation,
  row_title = "Complement Genes",
  row_title_gp = gpar(fontsize = 14, fontface = "bold"),
  row_names_gp = gpar(fontsize = 9),
  
  # Clustering
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  
  # Size
  width = unit(0.3 * ncol(scaled_matrix), "cm"),
  height = unit(0.4 * nrow(scaled_matrix), "cm"),
  
  # Legend
  heatmap_legend_param = list(
    title = "Expression\n(Z-score)",
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 10),
    legend_height = unit(4, "cm")
  )
)

# ==============================================================================
# Save Heatmap
# ==============================================================================

cat("Saving heatmap...\n")

# Calculate appropriate figure size
fig_width <- max(12, 0.3 * ncol(scaled_matrix) + 6)
fig_height <- max(8, 0.4 * nrow(scaled_matrix) + 4)

# Save as PNG
png(file.path(output_dir, "Complement_Analysis/complement_heatmap.png"), 
    width = fig_width, height = fig_height, units = "in", res = 300)
draw(main_heatmap, 
     heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()

# Save as PDF
pdf(file.path(output_dir, "Complement_Analysis/complement_heatmap.pdf"), 
    width = fig_width, height = fig_height)
draw(main_heatmap,
     heatmap_legend_side = "right", 
     annotation_legend_side = "right")
dev.off()

# ==============================================================================
# Create Summary Statistics
# ==============================================================================

cat("Creating summary statistics...\n")

# Cell type with highest complement expression per condition
summary_stats <- col_info %>%
  mutate(total_expression = colSums(complement_matrix[, sample])) %>%
  group_by(condition) %>%
  arrange(desc(total_expression)) %>%
  slice_head(n = 3) %>%
  select(condition, cell_type, total_expression)

# Most highly expressed complement genes overall
gene_summary <- data.frame(
  gene = rownames(complement_matrix),
  pathway = gene_pathways,
  mean_expression = rowMeans(complement_matrix),
  max_expression = rowMaxs(complement_matrix),
  stringsAsFactors = FALSE
) %>%
  arrange(desc(mean_expression))

# Save summaries
write.csv(summary_stats, 
          file.path(output_dir, "Complement_Analysis/top_complement_celltype_by_condition.csv"), 
          row.names = FALSE)

write.csv(gene_summary, 
          file.path(output_dir, "Complement_Analysis/complement_gene_summary.csv"), 
          row.names = FALSE)

# Save expression matrix
write.csv(complement_matrix, 
          file.path(output_dir, "Complement_Analysis/complement_pseudobulk_matrix.csv"))

# ==============================================================================
# Create Additional Visualizations
# ==============================================================================

cat("Creating additional complement visualizations...\n")

# Pathway-specific heatmaps
for(pathway in names(complement_pathways)) {
  pathway_genes <- intersect(complement_pathways[[pathway]], rownames(complement_matrix))
  
  if(length(pathway_genes) >= 2) {
    pathway_matrix <- scaled_matrix[pathway_genes, , drop = FALSE]
    
    pathway_heatmap <- Heatmap(
      pathway_matrix,
      name = paste(pathway, "\nExpression\n(Z-score)"),
      col = col_fun,
      top_annotation = col_annotation,
      column_title = paste(pathway, "Pathway Genes"),
      column_title_gp = gpar(fontsize = 14, fontface = "bold"),
      column_names_gp = gpar(fontsize = 8),
      column_names_rot = 45,
      row_names_gp = gpar(fontsize = 10),
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      width = unit(0.3 * ncol(pathway_matrix), "cm"),
      height = unit(0.5 * nrow(pathway_matrix), "cm")
    )
    
    # Save pathway-specific heatmap
    png(file.path(output_dir, "Complement_Analysis", paste0("complement_", tolower(pathway), "_pathway.png")), 
        width = max(10, 0.3 * ncol(pathway_matrix) + 4), 
        height = max(6, 0.5 * nrow(pathway_matrix) + 3), 
        units = "in", res = 300)
    draw(pathway_heatmap)
    dev.off()
  }
}

# ==============================================================================
# Completion Message
# ==============================================================================

cat("\n", rep("=", 60), "\n")
cat("COMPLEMENT HEATMAP ANALYSIS COMPLETED!\n")
cat(rep("=", 60), "\n")

cat("\nFiles created in:", file.path(output_dir, "Complement_Analysis"), "\n")
cat("- complement_heatmap.png/.pdf: Main complement gene heatmap\n")
cat("- complement_[pathway]_pathway.png: Pathway-specific heatmaps\n")
cat("- complement_pseudobulk_matrix.csv: Expression data matrix\n")
cat("- top_complement_celltype_by_condition.csv: Summary statistics\n")
cat("- complement_gene_summary.csv: Gene expression summary\n")

cat("\nKey features:\n")
cat("- Multi-layered column annotations (Condition + Cell Type)\n")
cat("- Row annotations by complement pathway\n")
cat("- Z-score normalized expression values\n")
cat("- Hierarchical clustering of genes and samples\n")
cat("- Color-coded pathway categories\n")

cat("\nComplement genes analyzed:", length(available_complement), "\n")
cat("Cell type-condition combinations:", ncol(complement_matrix), "\n")

cat("\nHeatmap shows complement gene expression patterns across all conditions and cell types!\n")