# ==============================================================================
# UMAP and Data Structure Visualization Script
# GSE207128 - Complement-OUD Analysis
# ==============================================================================

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(patchwork)
library(gridExtra)

# Set up paths
base_dir <- "/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus/GSE207128"
data_dir <- file.path(base_dir, "Outputs/04_Annotated_Data")
output_dir <- file.path(base_dir, "Outputs/06_Visualizations")

# Create output directories
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "UMAP_Plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "QC_Plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "Cell_Type_Plots"), recursive = TRUE, showWarnings = FALSE)

cat("Loading annotated Seurat object...\n")
seurat_obj <- readRDS(file.path(data_dir, "filtered_annotated_seurat.rds"))

cat("Data loaded successfully!\n")
cat("Cells:", ncol(seurat_obj), "\n")
cat("Genes:", nrow(seurat_obj), "\n")
cat("Conditions:", paste(unique(seurat_obj$condition), collapse = ", "), "\n")
cat("Cell types:", length(unique(seurat_obj$predicted_celltype)), "\n")

# ==============================================================================
# 1. BASIC UMAP PLOTS
# ==============================================================================

cat("\nCreating basic UMAP plots...\n")

# Define custom color palettes
condition_colors <- c("Naive" = "#2E8B57", "Dependent" = "#CD853F", "Withdrawal" = "#B22222")
sample_colors <- brewer.pal(9, "Set1")

# UMAP by condition
p1 <- DimPlot(seurat_obj, 
              reduction = "umap", 
              group.by = "condition",
              cols = condition_colors,
              pt.size = 0.5) +
  ggtitle("UMAP by Condition") +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# UMAP by sample
p2 <- DimPlot(seurat_obj, 
              reduction = "umap", 
              group.by = "orig.ident",
              cols = sample_colors,
              pt.size = 0.5) +
  ggtitle("UMAP by Sample") +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# UMAP by cell type
p3 <- DimPlot(seurat_obj, 
              reduction = "umap", 
              group.by = "predicted_celltype",
              label = TRUE,
              label.size = 3,
              pt.size = 0.5) +
  ggtitle("UMAP by Cell Type") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# UMAP by clusters
p4 <- DimPlot(seurat_obj, 
              reduction = "umap", 
              group.by = "seurat_clusters",
              label = TRUE,
              label.size = 3,
              pt.size = 0.5) +
  ggtitle("UMAP by Seurat Clusters") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# Save individual plots
ggsave(file.path(output_dir, "UMAP_Plots/UMAP_by_condition.png"), p1, width = 8, height = 6, dpi = 300)
ggsave(file.path(output_dir, "UMAP_Plots/UMAP_by_sample.png"), p2, width = 10, height = 6, dpi = 300)
ggsave(file.path(output_dir, "UMAP_Plots/UMAP_by_celltype.png"), p3, width = 12, height = 8, dpi = 300)
ggsave(file.path(output_dir, "UMAP_Plots/UMAP_by_clusters.png"), p4, width = 10, height = 8, dpi = 300)

# Combined panel
combined_basic <- (p1 | p2) / (p3 | p4)
ggsave(file.path(output_dir, "UMAP_Plots/UMAP_combined_basic.png"), combined_basic, width = 16, height = 12, dpi = 300)

# ==============================================================================
# 2. SPLIT UMAP PLOTS
# ==============================================================================

cat("Creating split UMAP plots...\n")

# Split by condition
p5 <- DimPlot(seurat_obj, 
              reduction = "umap", 
              group.by = "predicted_celltype",
              split.by = "condition",
              label = TRUE,
              label.size = 2.5,
              pt.size = 0.3,
              ncol = 3) +
  ggtitle("Cell Types Across Conditions") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

ggsave(file.path(output_dir, "UMAP_Plots/UMAP_celltype_split_condition.png"), p5, width = 18, height = 6, dpi = 300)

# Cell type by condition (colored)
p6 <- DimPlot(seurat_obj, 
              reduction = "umap", 
              group.by = "condition",
              split.by = "condition",
              cols = condition_colors,
              pt.size = 0.5,
              ncol = 3) +
  ggtitle("Conditions (Split View)") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

ggsave(file.path(output_dir, "UMAP_Plots/UMAP_condition_split.png"), p6, width = 18, height = 6, dpi = 300)

# ==============================================================================
# 3. QUALITY CONTROL VISUALIZATIONS
# ==============================================================================

cat("Creating QC visualization plots...\n")

# QC metrics on UMAP
qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")

qc_plots <- list()
for(feature in qc_features) {
  qc_plots[[feature]] <- FeaturePlot(seurat_obj, 
                                     features = feature,
                                     reduction = "umap",
                                     pt.size = 0.5) +
    scale_colour_viridis_c() +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
}

# Combine QC plots
qc_combined <- wrap_plots(qc_plots, ncol = 2)
ggsave(file.path(output_dir, "QC_Plots/QC_metrics_UMAP.png"), qc_combined, width = 12, height = 10, dpi = 300)

# QC violin plots by condition
qc_violin_plots <- list()
for(feature in qc_features) {
  qc_violin_plots[[feature]] <- VlnPlot(seurat_obj, 
                                        features = feature,
                                        group.by = "condition",
                                        cols = condition_colors,
                                        pt.size = 0) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.title.x = element_blank())
}

qc_violin_combined <- wrap_plots(qc_violin_plots, ncol = 2)
ggsave(file.path(output_dir, "QC_Plots/QC_violins_by_condition.png"), qc_violin_combined, width = 12, height = 10, dpi = 300)

# ==============================================================================
# 4. CELL TYPE DISTRIBUTION ANALYSIS
# ==============================================================================

cat("Creating cell type distribution plots...\n")

# Cell type proportions by condition
cell_props <- seurat_obj@meta.data %>%
  group_by(condition, predicted_celltype) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(condition) %>%
  mutate(proportion = count / sum(count) * 100)

# Stacked bar plot
p7 <- ggplot(cell_props, aes(x = condition, y = proportion, fill = predicted_celltype)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d() +
  theme_minimal() +
  labs(title = "Cell Type Proportions by Condition",
       x = "Condition",
       y = "Proportion (%)",
       fill = "Cell Type") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "Cell_Type_Plots/celltype_proportions_stacked.png"), p7, width = 10, height = 8, dpi = 300)

# Grouped bar plot
p8 <- ggplot(cell_props, aes(x = predicted_celltype, y = proportion, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = condition_colors) +
  theme_minimal() +
  labs(title = "Cell Type Proportions by Condition",
       x = "Cell Type",
       y = "Proportion (%)",
       fill = "Condition") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "Cell_Type_Plots/celltype_proportions_grouped.png"), p8, width = 12, height = 8, dpi = 300)

# Cell count table
cell_counts <- seurat_obj@meta.data %>%
  group_by(condition, predicted_celltype) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = condition, values_from = count, values_fill = 0)

write.csv(cell_counts, file.path(output_dir, "Cell_Type_Plots/cell_counts_by_condition.csv"), row.names = FALSE)

# ==============================================================================
# 5. SAMPLE INTEGRATION QUALITY
# ==============================================================================

cat("Creating integration quality plots...\n")

# Sample mixing on UMAP
p9 <- DimPlot(seurat_obj, 
              reduction = "umap", 
              group.by = "orig.ident",
              pt.size = 0.3) +
  ggtitle("Sample Integration Quality") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave(file.path(output_dir, "UMAP_Plots/sample_integration_quality.png"), p9, width = 10, height = 8, dpi = 300)

# Sample distribution by cell type
sample_celltype <- seurat_obj@meta.data %>%
  group_by(orig.ident, predicted_celltype) %>%
  summarise(count = n(), .groups = 'drop')

p10 <- ggplot(sample_celltype, aes(x = predicted_celltype, y = count, fill = orig.ident)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  labs(title = "Sample Contribution by Cell Type",
       x = "Cell Type",
       y = "Cell Count",
       fill = "Sample") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "Cell_Type_Plots/sample_contribution_celltype.png"), p10, width = 12, height = 8, dpi = 300)

# ==============================================================================
# 6. BRAIN-SPECIFIC MARKER VISUALIZATION
# ==============================================================================

cat("Creating brain marker visualization plots...\n")

# Define key brain markers
brain_markers <- c(
  "Neurod6",    # Excitatory neurons
  "Gad1",       # Inhibitory neurons  
  "Aif1",       # Microglia
  "Gfap",       # Astrocytes
  "Olig1",      # Oligodendrocytes
  "Pdgfra",     # OPCs
  "Cldn5",      # Endothelial
  "Dcn"         # Fibroblasts
)

# Check which markers are available
available_markers <- brain_markers[brain_markers %in% rownames(seurat_obj)]
cat("Available brain markers:", paste(available_markers, collapse = ", "), "\n")

if(length(available_markers) > 0) {
  # Feature plots for brain markers
  marker_plots <- list()
  for(marker in available_markers) {
    marker_plots[[marker]] <- FeaturePlot(seurat_obj, 
                                          features = marker,
                                          reduction = "umap",
                                          pt.size = 0.5) +
      scale_colour_viridis_c() +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  # Save individual marker plots
  for(marker in names(marker_plots)) {
    ggsave(file.path(output_dir, "Cell_Type_Plots", paste0("marker_", marker, ".png")), 
           marker_plots[[marker]], width = 8, height = 6, dpi = 300)
  }
  
  # Combined marker plot
  if(length(marker_plots) >= 4) {
    marker_combined <- wrap_plots(marker_plots[1:min(8, length(marker_plots))], ncol = 4)
    ggsave(file.path(output_dir, "Cell_Type_Plots/brain_markers_combined.png"), 
           marker_combined, width = 16, height = 8, dpi = 300)
  }
}

# ==============================================================================
# 7. SUMMARY STATISTICS
# ==============================================================================

cat("Creating summary statistics...\n")

# Create comprehensive summary
summary_stats <- list(
  total_cells = ncol(seurat_obj),
  total_genes = nrow(seurat_obj),
  conditions = unique(seurat_obj$condition),
  n_conditions = length(unique(seurat_obj$condition)),
  samples = unique(seurat_obj$orig.ident),
  n_samples = length(unique(seurat_obj$orig.ident)),
  cell_types = unique(seurat_obj$predicted_celltype),
  n_cell_types = length(unique(seurat_obj$predicted_celltype)),
  clusters = unique(seurat_obj$seurat_clusters),
  n_clusters = length(unique(seurat_obj$seurat_clusters))
)

# Cell counts summary
cells_by_condition <- table(seurat_obj$condition)
cells_by_sample <- table(seurat_obj$orig.ident)
cells_by_celltype <- table(seurat_obj$predicted_celltype)

# Save summary
sink(file.path(output_dir, "data_summary.txt"))
cat("GSE207128 Data Summary\n")
cat("======================\n\n")
cat("Total cells:", summary_stats$total_cells, "\n")
cat("Total genes:", summary_stats$total_genes, "\n")
cat("Number of conditions:", summary_stats$n_conditions, "\n")
cat("Number of samples:", summary_stats$n_samples, "\n")
cat("Number of cell types:", summary_stats$n_cell_types, "\n")
cat("Number of clusters:", summary_stats$n_clusters, "\n\n")

cat("Cells by condition:\n")
print(cells_by_condition)
cat("\nCells by sample:\n")
print(cells_by_sample)
cat("\nCells by cell type:\n")
print(cells_by_celltype)
sink()

# ==============================================================================
# COMPLETION MESSAGE
# ==============================================================================

cat("\n", rep("=", 60), "\n")
cat("VISUALIZATION SCRIPT COMPLETED SUCCESSFULLY!\n")
cat(rep("=", 60), "\n")

cat("\nFiles created in:", output_dir, "\n")
cat("- UMAP_Plots/: Basic and split UMAP visualizations\n")
cat("- QC_Plots/: Quality control metrics\n") 
cat("- Cell_Type_Plots/: Cell type distributions and markers\n")
cat("- data_summary.txt: Comprehensive data summary\n")

cat("\nKey visualizations created:\n")
cat("- Basic UMAPs (condition, sample, cell type, clusters)\n")
cat("- Split UMAPs by condition\n")
cat("- QC metrics on UMAP and violin plots\n")
cat("- Cell type proportion analysis\n")
cat("- Sample integration quality assessment\n")
cat("- Brain-specific marker expression\n")

cat("\nNext steps: Review visualizations and proceed with pathway analysis!\n")