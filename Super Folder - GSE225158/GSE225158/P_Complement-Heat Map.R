library(ComplexHeatmap)
library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(circlize)
library(viridis)

# 1. Load gene list
gene_list <- read_csv('GSE225158/KEGG outputs/kegg_complement_unique_genes.csv')
selected_genes <- gene_list$gene

# 2. Load expression matrix
expr <- read_csv('GSE225158/counts.csv') %>% column_to_rownames('gene')
expr_subset <- expr[rownames(expr) %in% selected_genes, ]

# 3. Load sample metadata
meta <- read_csv('GSE225158/meta.csv')

# Filter to include only columns in expression data
meta <- meta %>%
  filter(`...1` %in% colnames(expr_subset))

# 4. Simplify cell types
meta$celltype_simplified <- recode(meta$celltype3,
  'Astrocytes' = 'Astrocytes',
  'D1-Matrix' = 'Medium Spiny Neurons',
  'D1-Striosome' = 'Medium Spiny Neurons',
  'D1/D2-Hybrid' = 'Medium Spiny Neurons',
  'D2-Matrix' = 'Medium Spiny Neurons',
  'D2-Striosome' = 'Medium Spiny Neurons',
  'Int-CCK' = 'Interneurons',
  'Int-PTHLH' = 'Interneurons',
  'Int-SST' = 'Interneurons',
  'Int-TH' = 'Interneurons',
  'Oligos' = 'Oligodendrocytes',
  'Oligos_Pre' = 'Oligodendrocytes',
  'Endothelial' = 'Endothelial',
  'Microglia' = 'Microglia',
  'Mural' = 'Mural',
  'OPCs' = 'OPCs'
)

# Now filter based on the simplified cell types
glial_celltypes <- c("Oligodendrocytes", "Astrocytes", "OPCs", "Microglia")
meta <- meta %>%
  filter(celltype_simplified %in% glial_celltypes)

cat("Samples in meta after filtering for glial cells:", nrow(meta), "\n")

# 5. Color palettes
celltypes_simple <- unique(meta$celltype_simplified)
palette_func_simple <- colorRampPalette(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"))
celltype_colors_simple <- setNames(palette_func_simple(length(celltypes_simple)), celltypes_simple)
diagnosis_col <- setNames(c('#a6cee3', '#b2df8a'), c('OUD', 'None'))
sex_col <- setNames(c('black', 'grey'), c('F', 'M'))

# 6. Z-score normalization (row-wise) with optional log transform
use_log_transform <- FALSE  # Set to TRUE for log2 transform before scaling

expr_subset <- expr_subset[, meta$`...1`]
if (use_log_transform) {
  expr_scaled <- t(scale(t(log2(as.matrix(expr_subset) + 1))))
} else {
  expr_scaled <- t(scale(t(as.matrix(expr_subset))))
}

# 7. Aggregate by group and cell type
meta$Group <- paste(meta$Sex, meta$Dx_OUD, sep = "_")
agg_expr <- meta %>%
  mutate(CellType = celltype_simplified) %>%
  select(Sample = `...1`, Group, CellType) %>%
  left_join(
    as.data.frame(t(expr_scaled)) %>% rownames_to_column("Sample"),
    by = "Sample"
  ) %>%
  group_by(Group, CellType) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")

# 8. Reshape to matrix: rows=genes, columns=Group_CellType
agg_expr_long <- agg_expr %>%
  pivot_longer(-c(Group, CellType), names_to = "Gene", values_to = "Value") %>%
  unite("Group_CellType", Group, CellType, sep = "_") %>%
  pivot_wider(names_from = Group_CellType, values_from = Value)

agg_expr_mat <- as.matrix(agg_expr_long[,-1])
rownames(agg_expr_mat) <- agg_expr_long$Gene

# 9. Count frequencies for proportional barplot
freqs <- meta %>%
  count(Group, celltype_simplified) %>%
  unite("Group_CellType", Group, celltype_simplified, sep = "_")
col_freqs <- freqs$n[match(colnames(agg_expr_mat), freqs$Group_CellType)]
col_props <- col_freqs / sum(col_freqs)

# 10. Prepare annotation for columns
col_anno <- strsplit(colnames(agg_expr_mat), "_")
anno_df <- data.frame(
  Group = sapply(col_anno, `[`, 1),
  Dx = sapply(col_anno, `[`, 2),
  CellType = sapply(col_anno, `[`, 3)
)

column_ha <- HeatmapAnnotation(
  Proportion = anno_barplot(col_props, border = FALSE, gp = gpar(fill = "grey")),
  Group = anno_df$Group,
  Dx = anno_df$Dx,
  CellType = anno_df$CellType,
  col = list(
    Group = sex_col,
    Dx = diagnosis_col,
    CellType = celltype_colors_simple
  )
)

# 11. Plot heatmap
# Save to PDF with custom size
pdf('GSE225158/figures/heatmap_output.pdf', width = 15, height = 15)

Heatmap(
  agg_expr_mat,
  name = "Z-score",
  top_annotation = column_ha,
  col = colorRamp2(seq(-4, 4, length.out = 100), colorRampPalette(c("#50c878", "white", "#FF69B4"))(100)),
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  heatmap_width = unit(14, "cm"),
  heatmap_height = unit(32, "cm")
)

dev.off()