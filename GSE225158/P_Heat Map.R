library(ComplexHeatmap)
library(readr)
library(dplyr)
library(tibble)
library(circlize)

# 1. Load gene list
gene_list <- read_csv('GSE225158/KEGG outputs/kegg_complement_unique_genes.csv')
selected_genes <- gene_list$gene

# 2. Load expression matrix
expr <- read_csv('GSE225158/counts.csv') %>% column_to_rownames('gene')
expr_subset <- expr[rownames(expr) %in% selected_genes, ]

# 3. Load sample metadata
meta <- read_csv('GSE225158/meta.csv')
meta <- meta %>% filter(ID %in% colnames(expr_subset)) %>% arrange(match(ID, colnames(expr_subset)))

# 4. Generate color palette for cell types (works for any number)
celltypes <- unique(meta$celltype3)
n_celltypes <- length(celltypes)
palette_func <- colorRampPalette(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"))
celltype_colors <- setNames(palette_func(n_celltypes), celltypes)

# 5. Other color palettes
diagnosis_col <- setNames(c('#a6cee3', '#b2df8a'), c('OUD', 'None'))
sex_col <- setNames(c('black', 'yellow'), c('F', 'M'))

# 6. Column annotation
column_ha <- HeatmapAnnotation(
  CellType = meta$celltype3,
  Diagnosis = meta$Dx_OUD,
  Sex = meta$Sex,
  col = list(CellType = celltype_colors, Diagnosis = diagnosis_col, Sex = sex_col)
)

# 7. Plot heatmap
Heatmap(
  as.matrix(expr_subset),
  name = "Expression",
  top_annotation = column_ha,
  col = colorRamp2(c(min(expr_subset, na.rm=TRUE), 0, max(expr_subset, na.rm=TRUE)), c("blue", "white", "red")),
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_columns = TRUE,
  cluster_rows = TRUE
)