# R

library(ComplexHeatmap)
library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(circlize)
library(viridis)
library(stringr)
library(org.Mm.eg.db)
library(AnnotationDbi)

cat("==== 0. Load gene symbol <-> Ensembl mapping ====\n")
gene_map <- read_csv('GSE289002/ensembl_to_symbol.csv', show_col_types = FALSE)
ensembl_to_symbol <- setNames(gene_map$gene_symbol, gene_map$ensembl_id)
symbol_to_ensembl <- setNames(gene_map$ensembl_id, toupper(gene_map$gene_symbol))

cat("==== 1. Load complement cascade gene list ====\n")
comp_genes <- read_csv('GSE289002/kegg_outputs/GSE207128_mouse_complement_cascade_genes.csv', show_col_types = FALSE)$gene
comp_genes <- toupper(comp_genes) # ensure case-insensitive matching
ensembl_ids <- symbol_to_ensembl[comp_genes]
ensembl_ids <- ensembl_ids[!is.na(ensembl_ids)]
cat("Number of complement genes mapped to Ensembl IDs:", length(ensembl_ids), "\n")

cat("==== 2. Load expression matrix ====\n")
expr <- read_csv('GSE289002/GSE289002_mouse_raw_counts.csv', show_col_types = FALSE) %>% column_to_rownames('...1')
cat("expr dim:", dim(expr), "\n")

cat("==== 3. Subset expression matrix to complement genes ====\n")
expr_subset <- expr[rownames(expr) %in% ensembl_ids, , drop = FALSE]
cat("expr_subset dim:", dim(expr_subset), "\n")
if (nrow(expr_subset) == 0) stop("No complement genes found in expression matrix.")

cat("==== 4. Load sample metadata ====\n")
meta <- read_csv('GSE289002/mouse_metadata.csv', show_col_types = FALSE)
meta <- meta %>% dplyr::filter(title %in% colnames(expr_subset))
expr_subset <- expr_subset[, meta$title, drop = FALSE]

cat("==== 5. Z-score normalization (row-wise) ====\n")
expr_scaled <- t(scale(t(as.matrix(expr_subset))))

cat("==== 6. Aggregate by sex and treatment ====\n")
meta$Group <- paste(meta$sex, meta$treatment, sep = "_")
agg_expr <- meta %>%
  dplyr::select(title, Group) %>%
  left_join(
    as.data.frame(t(expr_scaled)) %>% tibble::rownames_to_column("title"),
    by = "title"
  ) %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(across(where(is.numeric), mean), .groups = "drop")

agg_expr_long <- agg_expr %>%
  tidyr::pivot_longer(-Group, names_to = "Gene", values_to = "Value") %>%
  tidyr::pivot_wider(names_from = Group, values_from = Value)
agg_expr_mat <- as.matrix(agg_expr_long[,-1])
rownames(agg_expr_mat) <- agg_expr_long$Gene

# Set y-axis to gene symbols
gene_symbols <- ensembl_to_symbol[rownames(agg_expr_mat)]
gene_symbols[is.na(gene_symbols) | gene_symbols == ""] <- rownames(agg_expr_mat)[is.na(gene_symbols) | gene_symbols == ""]
rownames(agg_expr_mat) <- gene_symbols

cat("==== 7. Prepare annotation for columns ====\n")
col_anno <- strsplit(colnames(agg_expr_mat), "_")
anno_df <- data.frame(
  Sex = sapply(col_anno, `[`, 1),
  Treatment = sapply(col_anno, `[`, 2)
)
sex_col <- setNames(c('black', 'grey'), unique(anno_df$Sex))
treatment_col <- setNames(viridis(length(unique(anno_df$Treatment))), unique(anno_df$Treatment))
column_ha <- HeatmapAnnotation(
  Sex = anno_df$Sex,
  Treatment = anno_df$Treatment,
  col = list(
    Sex = sex_col,
    Treatment = treatment_col
  )
)

cat("==== 8. Plot heatmap ====\n")
dir.create('GSE289002/figures', showWarnings = FALSE, recursive = TRUE)
pdf('GSE289002/figures/heatmap_output.pdf', width = 5, height = 10)
Heatmap(
  agg_expr_mat,
  name = "Z-score",
  top_annotation = column_ha,
  col = colorRamp2(seq(-1, 1, length.out = 100), colorRampPalette(c("#50c878", "white", "#FF69B4"))(100)),
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_columns = TRUE,
  cluster_rows = TRUE
)
dev.off()
cat("==== DONE ====\n")