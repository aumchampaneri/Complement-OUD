# ==============================================================================
# Complete DEG Visualization Suite - Phase 1
# GSE207128 - Complement-OUD Analysis
# Implements: Summary, Volcanos, Heatmaps, Complement Analysis, Success Metrics
# ==============================================================================

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(patchwork)
library(viridis)
library(readr)
library(tibble)
library(ggrepel)

# Set up paths
base_dir <- "/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus/GSE207128"
deg_dir <- file.path(base_dir, "Outputs/05_Analysis_Results/DEG_Analysis")
data_dir <- file.path(base_dir, "Outputs/04_Annotated_Data")
output_dir <- file.path(base_dir, "Outputs/06_Visualizations/DEG_Complete_Suite")

# Create output directories
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "01_Summary_Plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "02_Volcano_Plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "03_Heatmaps"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "04_Complement_Analysis"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "05_Success_Metrics"), recursive = TRUE, showWarnings = FALSE)

cat("Loading differential expression results...\n")

# Load results - handle both possible file structures
tryCatch({
  all_results <- readRDS(file.path(deg_dir, "edgeR_all_results.rds"))
  deg_summary <- readRDS(file.path(deg_dir, "edgeR_summary_stats.rds"))
}, error = function(e) {
  stop("Could not load DEG results. Please ensure edgeR analysis has been completed first.")
})

# Also load Seurat object for additional context
seurat_obj <- readRDS(file.path(data_dir, "filtered_annotated_seurat.rds"))

cat("Results loaded successfully!\n")
cat("Total analyses:", length(all_results), "\n")

# ==============================================================================
# 1. SUMMARY OVERVIEW PLOTS
# ==============================================================================

cat("\n=== Creating Summary Overview Plots ===\n")

# Color schemes
comparison_colors <- c(
  "Dep_vs_Naive" = "#FF6B6B", 
  "With_vs_Dep" = "#4ECDC4", 
  "With_vs_Naive" = "#45B7D1",
  "vs_Naive" = "#FF6B6B",
  "vs_Dep" = "#4ECDC4",
  "Naive" = "#45B7D1"
)
regulation_colors <- c("Upregulated" = "#FF6B6B", "Downregulated" = "#4ECDC4", "Not Significant" = "grey70")

# 1A. DEG counts per cell type and comparison
p1a <- ggplot(deg_summary, aes(x = reorder(cell_type, significant_genes), y = significant_genes, fill = comparison)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = comparison_colors) +
  coord_flip() +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  labs(title = "Significant DEGs per Cell Type and Comparison", 
       x = "Cell Type", y = "Number of Significant Genes", fill = "Comparison")

ggsave(file.path(output_dir, "01_Summary_Plots/deg_counts_by_celltype.png"), p1a, width = 12, height = 10, dpi = 300)

# 1B. Up/Down regulation summary
deg_summary_regulation <- deg_summary %>%
  select(cell_type, comparison, upregulated, downregulated) %>%
  pivot_longer(cols = c(upregulated, downregulated), 
               names_to = "direction", values_to = "count") %>%
  mutate(direction = case_when(
    direction == "upregulated" ~ "Upregulated",
    direction == "downregulated" ~ "Downregulated"
  ))

p1b <- ggplot(deg_summary_regulation, aes(x = reorder(cell_type, count), y = count)) +
  geom_bar(stat = "identity", aes(fill = direction), position = "stack", alpha = 0.8) +
  scale_fill_manual(values = c("Upregulated" = "#FF6B6B", "Downregulated" = "#4ECDC4")) +
  coord_flip() +
  facet_wrap(~comparison, ncol = 3) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 8),
    strip.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  labs(title = "Up vs Down Regulated Genes by Cell Type and Comparison", 
       x = "Cell Type", y = "Number of Genes", fill = "Regulation")

ggsave(file.path(output_dir, "01_Summary_Plots/regulation_direction_summary.png"), p1b, width = 16, height = 10, dpi = 300)

# 1C. Total genes tested vs significant genes
p1c <- ggplot(deg_summary, aes(x = total_genes, y = significant_genes, color = comparison)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.5) +
  scale_color_manual(values = comparison_colors) +
  theme_minimal() +
  labs(title = "Genes Tested vs Significant Genes Found",
       x = "Total Genes Tested", y = "Significant Genes", color = "Comparison") +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "01_Summary_Plots/genes_tested_vs_significant.png"), p1c, width = 10, height = 8, dpi = 300)

cat("✓ Summary plots completed\n")

# ==============================================================================
# 2. ENHANCED VOLCANO PLOTS
# ==============================================================================

cat("\n=== Creating Enhanced Volcano Plots ===\n")

# Function to create enhanced volcano plot
create_enhanced_volcano <- function(results_data, title, highlight_genes = NULL) {
  
  # Prepare data
  results_data$neg_log10_FDR <- -log10(results_data$FDR)
  results_data$neg_log10_FDR[is.infinite(results_data$neg_log10_FDR)] <- max(results_data$neg_log10_FDR[!is.infinite(results_data$neg_log10_FDR)], na.rm = TRUE)
  
  # Add significance categories
  results_data$category <- "Not Significant"
  results_data$category[results_data$FDR < 0.05 & results_data$logFC > 0.5] <- "Upregulated"
  results_data$category[results_data$FDR < 0.05 & results_data$logFC < -0.5] <- "Downregulated"
  results_data$category[results_data$FDR < 0.05 & abs(results_data$logFC) <= 0.5] <- "Significant (low FC)"
  
  # Color scheme
  colors <- c("Not Significant" = "grey70", 
              "Significant (low FC)" = "lightblue",
              "Upregulated" = "#FF6B6B", 
              "Downregulated" = "#4ECDC4")
  
  # Create base plot
  p <- ggplot(results_data, aes(x = logFC, y = neg_log10_FDR)) +
    geom_point(aes(color = category), alpha = 0.6, size = 1) +
    scale_color_manual(values = colors) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5, color = "red") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", alpha = 0.5, color = "blue") +
    theme_minimal() +
    labs(title = title, x = "log2 Fold Change", y = "-log10(FDR)", color = "Significance") +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
    )
  
  # Highlight specific genes if provided
  if(!is.null(highlight_genes)) {
    highlight_data <- results_data[results_data$gene %in% highlight_genes, ]
    if(nrow(highlight_data) > 0) {
      p <- p + geom_point(data = highlight_data, color = "black", size = 2, shape = 21, fill = "yellow", alpha = 0.8)
    }
  }
  
  return(p)
}

# Combine all results for comprehensive analysis
all_combined <- do.call(rbind, all_results)
all_combined$analysis_id <- paste(all_combined$cell_type, all_combined$comparison, sep = "_")

# Get unique comparisons
comparisons <- unique(all_combined$comparison)

# 2A. Create volcano plots for each comparison (all cell types combined)
volcano_plots_comparison <- list()
for(comp in comparisons) {
  comp_data <- all_combined[all_combined$comparison == comp, ]
  volcano_plots_comparison[[comp]] <- create_enhanced_volcano(comp_data, paste("All Cell Types:", comp))
  
  # Save individual comparison volcano
  ggsave(file.path(output_dir, "02_Volcano_Plots", paste0("volcano_", comp, "_all_celltypes.png")), 
         volcano_plots_comparison[[comp]], width = 10, height = 8, dpi = 300)
}

# 2B. Combined comparison volcanos
combined_volcano_comp <- wrap_plots(volcano_plots_comparison, ncol = 3)
ggsave(file.path(output_dir, "02_Volcano_Plots/volcano_all_comparisons_combined.png"), 
       combined_volcano_comp, width = 24, height = 8, dpi = 300)

# 2C. Create cell type-specific volcanos for most interesting cell types
top_cell_types <- deg_summary %>%
  group_by(cell_type) %>%
  summarise(total_sig = sum(significant_genes), .groups = 'drop') %>%
  arrange(desc(total_sig)) %>%
  slice_head(n = 6) %>%
  pull(cell_type)

volcano_plots_celltype <- list()
for(ct in top_cell_types) {
  ct_data <- all_combined[all_combined$cell_type == ct, ]
  volcano_plots_celltype[[ct]] <- create_enhanced_volcano(ct_data, paste("Cell Type:", ct))
  
  # Save individual cell type volcano
  ggsave(file.path(output_dir, "02_Volcano_Plots", paste0("volcano_", gsub("[^A-Za-z0-9]", "_", ct), ".png")), 
         volcano_plots_celltype[[ct]], width = 10, height = 8, dpi = 300)
}

# Combined cell type volcanos
combined_volcano_ct <- wrap_plots(volcano_plots_celltype, ncol = 3)
ggsave(file.path(output_dir, "02_Volcano_Plots/volcano_top_celltypes_combined.png"), 
       combined_volcano_ct, width = 18, height = 12, dpi = 300)

cat("✓ Volcano plots completed\n")

# ==============================================================================
# 3. TOP GENES HEATMAPS
# ==============================================================================

cat("\n=== Creating Top Genes Heatmaps ===\n")

# 3A. Get top significant genes across all analyses
top_genes_global <- all_combined %>%
  filter(FDR < 0.01) %>%
  group_by(gene) %>%
  summarise(
    times_significant = n(),
    avg_logFC = mean(logFC),
    min_FDR = min(FDR),
    max_abs_logFC = max(abs(logFC)),
    .groups = 'drop'
  ) %>%
  arrange(desc(times_significant), min_FDR) %>%
  slice_head(n = 50)

cat("Top genes identified:", nrow(top_genes_global), "\n")

# Create logFC matrix for heatmap
heatmap_data <- all_combined %>%
  filter(gene %in% top_genes_global$gene) %>%
  select(gene, cell_type, comparison, logFC) %>%
  mutate(analysis = paste(cell_type, comparison, sep = "_")) %>%
  select(gene, analysis, logFC) %>%
  pivot_wider(names_from = analysis, values_from = logFC, values_fill = 0) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# Remove columns/analyses with all zeros
heatmap_data <- heatmap_data[, colSums(abs(heatmap_data)) > 0]

# 3B. Create comprehensive heatmap
col_fun <- colorRamp2(c(-2, -1, 0, 1, 2), c("#4ECDC4", "#B8E6E6", "white", "#FFB3B3", "#FF6B6B"))

# Prepare column annotations
col_info <- data.frame(
  analysis = colnames(heatmap_data),
  stringsAsFactors = FALSE
)
col_info$cell_type <- sapply(strsplit(col_info$analysis, "_(?=[^_]*_[^_]*$)", perl = TRUE), `[`, 1)
col_info$comparison <- sapply(strsplit(col_info$analysis, "_"), function(x) paste(tail(x, 2), collapse = "_"))

# Create cell type colors properly
unique_cell_types <- unique(col_info$cell_type)
cell_type_colors <- rainbow(length(unique_cell_types))
names(cell_type_colors) <- unique_cell_types

# Create annotations
col_annotation <- HeatmapAnnotation(
  Comparison = col_info$comparison,
  `Cell Type` = col_info$cell_type,
  col = list(
    Comparison = comparison_colors,
    `Cell Type` = cell_type_colors
  ),
  annotation_name_gp = gpar(fontsize = 10)
)

# Create main heatmap
ht_main <- Heatmap(
  heatmap_data,
  name = "log2FC",
  col = col_fun,
  top_annotation = col_annotation,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  column_title = "Top 50 Significant Genes Across All Analyses",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_param = list(
    title = "log2FC",
    title_gp = gpar(fontsize = 12, fontface = "bold")
  )
)

# Save heatmap
png(file.path(output_dir, "03_Heatmaps/top50_genes_heatmap.png"), width = 16, height = 12, units = "in", res = 300)
draw(ht_main)
dev.off()

# 3C. Cell type-specific heatmap (top genes per cell type)
ct_specific_genes <- all_combined %>%
  filter(FDR < 0.05) %>%
  group_by(cell_type) %>%
  arrange(FDR) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  pull(gene) %>%
  unique()

if(length(ct_specific_genes) > 5) {
  ct_heatmap_data <- all_combined %>%
    filter(gene %in% ct_specific_genes) %>%
    select(gene, cell_type, comparison, logFC) %>%
    mutate(analysis = paste(cell_type, comparison, sep = "_")) %>%
    select(gene, analysis, logFC) %>%
    pivot_wider(names_from = analysis, values_from = logFC, values_fill = 0) %>%
    column_to_rownames("gene") %>%
    as.matrix()
  
  ht_ct <- Heatmap(
    ct_heatmap_data,
    name = "log2FC",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_names_rot = 45,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    column_title = "Top Cell Type-Specific Genes",
    column_title_gp = gpar(fontsize = 14, fontface = "bold")
  )
  
  png(file.path(output_dir, "03_Heatmaps/celltype_specific_genes_heatmap.png"), width = 14, height = 10, units = "in", res = 300)
  draw(ht_ct)
  dev.off()
}

cat("✓ Heatmap analysis completed\n")

# ==============================================================================
# 4. COMPLEMENT-SPECIFIC ANALYSIS
# ==============================================================================

cat("\n=== Creating Complement-Specific Analysis ===\n")

# Define comprehensive complement genes
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
  # Other components
  "Cfp", "Clu", "Vtn", "Serping1"
)

# 4A. Extract complement results
complement_results <- all_combined %>%
  filter(gene %in% complement_genes) %>%
  arrange(FDR)

cat("Complement genes found in results:", length(unique(complement_results$gene)), "\n")

# 4B. Complement-specific volcano plots
for(comp in comparisons) {
  comp_data <- all_combined[all_combined$comparison == comp, ]
  p_comp <- create_enhanced_volcano(comp_data, paste("Complement Genes -", comp), 
                                   highlight_genes = complement_genes)
  
  ggsave(file.path(output_dir, "04_Complement_Analysis", paste0("volcano_complement_", comp, ".png")), 
         p_comp, width = 10, height = 8, dpi = 300)
}

# 4C. Complement gene summary table
complement_summary <- complement_results %>%
  group_by(gene) %>%
  summarise(
    significant_analyses = sum(FDR < 0.05),
    total_analyses = n(),
    avg_logFC = mean(logFC),
    min_FDR = min(FDR),
    max_abs_logFC = max(abs(logFC)),
    .groups = 'drop'
  ) %>%
  arrange(desc(significant_analyses), min_FDR)

write.csv(complement_summary, file.path(output_dir, "04_Complement_Analysis/complement_gene_summary.csv"), row.names = FALSE)

# 4D. Complement pathway heatmap
if(nrow(complement_results) > 0) {
  
  # Define pathway categories
  complement_pathways <- list(
    "Classical" = c("C1qa", "C1qb", "C1qc", "C1r", "C1s", "C4a", "C4b", "C2"),
    "Alternative" = c("C3", "Cfb", "Cfd", "C5", "C6", "C7", "C8a", "C8b", "C8g", "C9"),
    "Lectin" = c("Mbl1", "Mbl2", "Masp1", "Masp2", "Fcn1", "Fcn2", "Fcn3"),
    "Regulators" = c("Cfh", "Cfi", "C4bp", "Cd55", "Cd46", "Cd35", "Cr1", "Cr2"),
    "Receptors" = c("C3ar1", "C5ar1", "C5ar2", "Itgam", "Itgax", "Itgb2"),
    "Other" = c("Cfp", "Clu", "Vtn", "Serping1")
  )
  
  complement_heatmap_data <- complement_results %>%
    select(gene, cell_type, comparison, logFC) %>%
    mutate(analysis = paste(cell_type, comparison, sep = "_")) %>%
    select(gene, analysis, logFC) %>%
    pivot_wider(names_from = analysis, values_from = logFC, values_fill = 0) %>%
    column_to_rownames("gene") %>%
    as.matrix()
  
  # Create pathway annotations
  gene_pathways <- rep("Other", nrow(complement_heatmap_data))
  names(gene_pathways) <- rownames(complement_heatmap_data)
  
  for(pathway in names(complement_pathways)) {
    pathway_genes <- intersect(complement_pathways[[pathway]], rownames(complement_heatmap_data))
    gene_pathways[pathway_genes] <- pathway
  }
  
  pathway_colors <- c(
    "Classical" = "#FF6B6B",
    "Alternative" = "#4ECDC4", 
    "Lectin" = "#45B7D1",
    "Regulators" = "#96CEB4",
    "Receptors" = "#FFEAA7",
    "Other" = "#DDA0DD"
  )
  
  row_annotation <- rowAnnotation(
    Pathway = gene_pathways,
    col = list(Pathway = pathway_colors),
    annotation_name_gp = gpar(fontsize = 10)
  )
  
  ht_complement <- Heatmap(
    complement_heatmap_data,
    name = "log2FC",
    col = col_fun,
    left_annotation = row_annotation,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_names_rot = 45,
    row_names_gp = gpar(fontsize = 9),
    column_names_gp = gpar(fontsize = 8),
    column_title = "Complement Genes Expression Changes",
    column_title_gp = gpar(fontsize = 14, fontface = "bold")
  )
  
  png(file.path(output_dir, "04_Complement_Analysis/complement_genes_heatmap.png"), 
      width = 14, height = 10, units = "in", res = 300)
  draw(ht_complement)
  dev.off()
}

cat("✓ Complement analysis completed\n")

# ==============================================================================
# 5. SUCCESS METRICS AND METHOD VALIDATION
# ==============================================================================

cat("\n=== Creating Success Metrics Analysis ===\n")

# 5A. Analysis success rate by cell type
success_data <- deg_summary %>%
  group_by(cell_type) %>%
  summarise(
    total_analyses = n(),
    successful_analyses = sum(significant_genes > 0),
    success_rate = successful_analyses / total_analyses,
    avg_genes_tested = mean(total_genes),
    avg_significant = mean(significant_genes),
    .groups = 'drop'
  ) %>%
  arrange(desc(success_rate))

p5a <- ggplot(success_data, aes(x = reorder(cell_type, success_rate), y = success_rate)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  geom_text(aes(label = paste0(round(success_rate*100, 1), "%")), 
            hjust = -0.1, size = 3) +
  coord_flip() +
  theme_minimal() +
  ylim(0, 1.1) +
  labs(title = "Analysis Success Rate by Cell Type", 
       x = "Cell Type", y = "Success Rate") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave(file.path(output_dir, "05_Success_Metrics/analysis_success_rate.png"), p5a, width = 10, height = 8, dpi = 300)

# 5B. Genes tested vs success rate
p5b <- ggplot(success_data, aes(x = avg_genes_tested, y = success_rate)) +
  geom_point(size = 3, alpha = 0.7, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3) +
  ggrepel::geom_text_repel(aes(label = cell_type), size = 3) +
  theme_minimal() +
  labs(title = "Genes Tested vs Success Rate",
       x = "Average Genes Tested", y = "Success Rate") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave(file.path(output_dir, "05_Success_Metrics/genes_tested_vs_success.png"), p5b, width = 10, height = 8, dpi = 300)

# 5C. Effect size distribution
effect_size_data <- all_combined %>%
  filter(FDR < 0.05) %>%
  mutate(effect_category = case_when(
    abs(logFC) < 0.5 ~ "Small",
    abs(logFC) < 1.0 ~ "Medium", 
    abs(logFC) < 2.0 ~ "Large",
    TRUE ~ "Very Large"
  ))

p5c <- ggplot(effect_size_data, aes(x = effect_category, fill = comparison)) +
  geom_bar(position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = comparison_colors) +
  theme_minimal() +
  labs(title = "Effect Size Distribution of Significant Genes",
       x = "Effect Size Category", y = "Number of Genes", fill = "Comparison") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "bottom")

ggsave(file.path(output_dir, "05_Success_Metrics/effect_size_distribution.png"), p5c, width = 10, height = 8, dpi = 300)

# 5D. Cell type performance matrix
performance_matrix <- deg_summary %>%
  select(cell_type, comparison, significant_genes) %>%
  pivot_wider(names_from = comparison, values_from = significant_genes, values_fill = 0) %>%
  column_to_rownames("cell_type") %>%
  as.matrix()

ht_performance <- Heatmap(
  performance_matrix,
  name = "Significant\nGenes",
  col = viridis::viridis(100),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_rot = 0,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.0f", performance_matrix[i, j]), x, y, gp = gpar(fontsize = 8, col = "white"))
  },
  column_title = "Cell Type Performance Across Comparisons",
  column_title_gp = gpar(fontsize = 14, fontface = "bold")
)

png(file.path(output_dir, "05_Success_Metrics/celltype_performance_matrix.png"), 
    width = 10, height = 8, units = "in", res = 300)
draw(ht_performance)
dev.off()

cat("✓ Success metrics completed\n")

# ==============================================================================
# FINAL SUMMARY REPORT
# ==============================================================================

# Create comprehensive summary
total_significant <- sum(deg_summary$significant_genes)
total_analyses <- nrow(deg_summary)
successful_analyses <- sum(deg_summary$significant_genes > 0)

# Generate summary report
sink(file.path(output_dir, "Complete_DEG_Analysis_Summary.txt"))
cat("=== COMPLETE DEG ANALYSIS SUMMARY ===\n")
cat("Analysis Date:", Sys.Date(), "\n")
cat("Total Analyses Performed:", total_analyses, "\n")
cat("Successful Analyses:", successful_analyses, "\n")
cat("Success Rate:", round(successful_analyses/total_analyses*100, 1), "%\n")
cat("Total Significant Genes Found:", total_significant, "\n")
cat("Average Significant Genes per Analysis:", round(total_significant/total_analyses, 1), "\n\n")

cat("=== TOP PERFORMING CELL TYPES ===\n")
top_performers <- success_data %>% arrange(desc(avg_significant)) %>% slice_head(n = 5)
for(i in 1:nrow(top_performers)) {
  cat(paste0(i, ". ", top_performers$cell_type[i], ": ", 
             round(top_performers$avg_significant[i], 1), " avg significant genes\n"))
}

cat("\n=== COMPLEMENT GENES SUMMARY ===\n")
if(nrow(complement_results) > 0) {
  sig_complement <- sum(complement_results$FDR < 0.05)
  cat("Complement genes tested:", length(unique(complement_results$gene)), "\n")
  cat("Significant complement findings:", sig_complement, "\n")
  if(sig_complement > 0) {
    top_complement <- complement_results %>% 
      filter(FDR < 0.05) %>% 
      arrange(FDR) %>% 
      slice_head(n = 5)
    cat("Top significant complement genes:\n")
    for(i in 1:nrow(top_complement)) {
      cat(paste0("  ", top_complement$gene[i], " (", top_complement$cell_type[i], 
                 ", ", top_complement$comparison[i], "): FDR = ", 
                 format(top_complement$FDR[i], scientific = TRUE, digits = 3), "\n"))
    }
  }
} else {
  cat("No complement genes found in results\n")
}

cat("\n=== FILES CREATED ===\n")
cat("01_Summary_Plots: Overview and regulation summaries\n")
cat("02_Volcano_Plots: Enhanced volcano plots with gene highlighting\n") 
cat("03_Heatmaps: Top genes and cell type-specific heatmaps\n")
cat("04_Complement_Analysis: Complement-focused visualizations\n")
cat("05_Success_Metrics: Method validation and performance analysis\n")
sink()

# ==============================================================================
# COMPLETION MESSAGE
# ==============================================================================

cat("\n", rep("=", 80), "\n")
cat("COMPLETE DEG VISUALIZATION SUITE FINISHED SUCCESSFULLY!\n")
cat(rep("=", 80), "\n")

cat("\nOutput directory:", output_dir, "\n")
cat("Total visualizations created: ~25 plots and heatmaps\n")
cat("Summary report: Complete_DEG_Analysis_Summary.txt\n")

cat("\n🎯 KEY FINDINGS:\n")
cat("• Total significant DEGs:", total_significant, "\n")
cat("• Analysis success rate:", round(successful_analyses/total_analyses*100, 1), "%\n")
cat("• Complement genes analyzed:", length(unique(complement_results$gene)), "\n")
cat("• Top performing cell type:", success_data$cell_type[1], "\n")

cat("\n📊 VISUALIZATION CATEGORIES COMPLETED:\n")
cat("✓ 1. Summary Overview Plots - DEG counts and regulation patterns\n")
cat("✓ 2. Enhanced Volcano Plots - Significance visualization with highlighting\n")
cat("✓ 3. Top Genes Heatmaps - Expression patterns across analyses\n")
cat("✓ 4. Complement Analysis - Pathway-specific focus\n")
cat("✓ 5. Success Metrics - Method validation and performance\n")

cat("\n🚀 READY FOR PHASE 2: Cross-method validation with pseudobulk results!\n")