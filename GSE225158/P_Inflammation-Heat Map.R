# R

library(ComplexHeatmap)
library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(circlize)
library(viridis)

selected_pathway <- NULL # e.g., "TNF signaling pathway" or set via commandArgs

# Possible pathway options (from kegg_inflammatory_pathways.csv):
# "TNF signaling pathway"
# "IL-17 signaling pathway"
# "NF-kappa B signaling pathway"
# "NOD-like receptor signaling pathway"
# "Toll-like receptor signaling pathway"
# "Chemokine signaling pathway"
# "Cytokine-cytokine receptor interaction"
# "Jak-STAT signaling pathway"
# "MAPK signaling pathway"
# (Add or update this list based on the actual unique values in your CSV)

# 1. Load gene list
gene_list <- read_csv('GSE225158/KEGG outputs/kegg_unique_genes.csv', col_names = "gene")
selected_genes <- gene_list$gene

# Load pathway mapping
pathway_map <- read_csv('GSE225158/KEGG outputs/kegg_inflammatory_pathways.csv') # columns: gene, pathway

# Subset by pathway if argument is provided
if (!is.null(selected_pathway)) {
  pathway_genes <- pathway_map %>%
    filter(pathway == selected_pathway) %>%
    pull(gene)
  selected_genes <- intersect(selected_genes, pathway_genes)
  cat("Subsetting to pathway:", selected_pathway, "with", length(selected_genes), "genes\n")
}

# 2. Load expression matrix
expr <- read_csv('GSE225158/counts.csv') %>% column_to_rownames('gene')
cat("Genes in expression matrix:", nrow(expr), "\n")
expr_subset <- expr[rownames(expr) %in% selected_genes, ]
cat("Genes after subsetting:", nrow(expr_subset), "\n")

# 3. Load sample metadata
meta <- read_csv('GSE225158/meta.csv')
meta <- meta %>%
  filter(`...1` %in% colnames(expr_subset))
cat("Samples in meta after filtering:", nrow(meta), "\n")

# 4. Color palettes for annotation
diagnosis_col <- setNames(c('#a6cee3', '#b2df8a'), c('OUD', 'None'))
sex_col <- setNames(c('black', 'grey'), c('F', 'M'))
celltypes <- unique(meta$celltype3)
palette_func <- colorRampPalette(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"))
celltype_colors <- setNames(palette_func(length(celltypes)), celltypes)

# 5. Z-score normalization (row-wise)
use_log_transform <- TRUE

expr_subset <- expr_subset[, meta$`...1`, drop = FALSE]
cat("expr_subset dim after column filter:", dim(expr_subset), "\n")
if (use_log_transform) {
  expr_scaled <- t(scale(t(log2(as.matrix(expr_subset) + 1))))
} else {
  expr_scaled <- t(scale(t(as.matrix(expr_subset))))
}
cat("expr_scaled dim:", dim(expr_scaled), "\n")

# 6. Prepare column annotation and grouping
anno_df <- meta %>%
  arrange(match(`...1`, colnames(expr_scaled))) %>%
  select(Sex, Dx_OUD, celltype3)

# Grouping factor: first by Sex, then by Dx_OUD
column_group <- paste(anno_df$Sex, anno_df$Dx_OUD, sep = "_")
# Ensure correct order of levels: all F_None, F_OUD, M_None, M_OUD, etc.
sex_levels <- unique(anno_df$Sex)
dx_levels <- unique(anno_df$Dx_OUD)
group_levels <- as.vector(outer(sex_levels, dx_levels, paste, sep = "_"))
column_group <- factor(column_group, levels = group_levels)

column_ha <- HeatmapAnnotation(
  Sex = anno_df$Sex,
  Dx = anno_df$Dx_OUD,
  CellType = anno_df$celltype3,
  col = list(
    Sex = sex_col,
    Dx = diagnosis_col,
    CellType = celltype_colors
  )
)

# 7. Remove rows with all NA/NaN/Inf
expr_scaled <- expr_scaled[apply(expr_scaled, 1, function(x) all(is.finite(x))), , drop = FALSE]
cat("expr_scaled dim after filtering:", dim(expr_scaled), "\n")

# 8. Plot heatmap with column_split
# R

# Calculate dynamic sizes
heatmap_height <- unit(nrow(expr_scaled) * 0.3, "cm")
heatmap_width <- unit(ncol(expr_scaled) * 0.3, "cm")

pdf('GSE225158/figures/heatmap_output.pdf',
    width = as.numeric(heatmap_width, "cm") + 2,
    height = as.numeric(heatmap_height, "cm") + 2)

Heatmap(
  expr_scaled,
  name = "Z-score",
  top_annotation = column_ha,
  col = colorRamp2(seq(-4, 4, length.out = 100), colorRampPalette(c("#50c878", "white", "#FF69B4"))(100)),
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  column_split = column_group,
  heatmap_width = heatmap_width,
  heatmap_height = heatmap_height,
  row_names_gp = gpar(fontsize = 6) # adjust as needed
)
dev.off()