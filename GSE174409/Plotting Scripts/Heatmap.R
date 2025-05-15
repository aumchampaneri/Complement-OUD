library(ComplexHeatmap)
library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(circlize)
library(viridis)
library(stringr)
library(org.Hs.eg.db)
library(AnnotationDbi)

cat("==== 0. Load gene symbol <-> Ensembl mapping ====\n")
gene_map <- read_csv('GSE174409/NeuroinflammationResults/ensembl_to_symbol_mapping.csv', show_col_types = FALSE)
ensembl_to_symbol <- setNames(gene_map$Symbol, gene_map$EnsemblID)
symbol_to_ensembl <- setNames(gene_map$EnsemblID, toupper(gene_map$Symbol))

cat("==== 1. Load complement cascade gene list ====\n")
comp_genes <- read_csv('/Users/aumchampaneri/PycharmProjects/Complement-OUD/Super Folder - GSE225158/GSE225158/KEGG outputs/kegg_complement_unique_genes.csv', show_col_types = FALSE)$gene
comp_genes <- toupper(comp_genes) # ensure case-insensitive matching
ensembl_ids <- symbol_to_ensembl[comp_genes]
ensembl_ids <- ensembl_ids[!is.na(ensembl_ids)]
cat("Number of complement genes mapped to Ensembl IDs:", length(ensembl_ids), "\n")

cat("==== 2. Load expression matrix ====\n")
expr <- readRDS('/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE174409/QC/logcpm_filtered_normalized.rds')
cat("expr dim:", dim(expr), "\n")

cat("==== 3. Subset expression matrix to complement genes ====\n")
expr_subset <- expr[rownames(expr) %in% ensembl_ids, , drop = FALSE]
cat("expr_subset dim:", dim(expr_subset), "\n")
if (nrow(expr_subset) == 0) stop("No complement genes found in expression matrix.")

cat("==== 4. Load sample metadata ====\n")
meta <- readRDS('/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE174409/QC/metadata.rds')

# Print available columns in metadata to help identify the right one
cat("Available columns in metadata:", paste(colnames(meta), collapse=", "), "\n")

# Look for potential sample identifier columns
potential_id_cols <- c("title", "sample_name", "sample_id", "sample", "id", "geo_accession", "run", "sample_accession")
found_cols <- intersect(potential_id_cols, colnames(meta))

if (length(found_cols) > 0) {
    # Use the first matching column
    id_col <- found_cols[1]
    cat("Using column '", id_col, "' as sample identifier\n", sep="")
    meta$title <- meta[[id_col]]
} else {
    # If no common ID columns found, try to identify by matching with expression data
    cat("Checking which columns might contain sample identifiers...\n")

    # Check character/factor columns
    char_cols <- sapply(meta, function(x) is.character(x) || is.factor(x))
    if (any(char_cols)) {
        char_col_names <- names(char_cols)[char_cols]

        # Check which column values match with column names in expression data
        match_counts <- sapply(char_col_names, function(col) {
            sum(as.character(meta[[col]]) %in% colnames(expr_subset))
        })

        if (any(match_counts > 0)) {
            best_col <- char_col_names[which.max(match_counts)]
            cat("Using column '", best_col, "' as sample identifier (",
                max(match_counts), " matches)\n", sep="")
            meta$title <- as.character(meta[[best_col]])
        } else {
            stop("Could not identify a suitable sample identifier column.")
        }
    } else {
        stop("No character or factor columns found that could be sample identifiers.")
    }
}

# Ensure 'title' is a character vector
meta$title <- as.character(meta$title)

# Filter metadata to match expression matrix columns
meta <- meta %>% dplyr::filter(title %in% colnames(expr_subset))
if (nrow(meta) == 0) {
    stop("No matches found between metadata and expression data.")
}
expr_subset <- expr_subset[, meta$title, drop = FALSE]
cat("Matched", nrow(meta), "samples between metadata and expression data\n")

cat("==== 5. Z-score normalization (row-wise) ====\n")
expr_scaled <- t(scale(t(as.matrix(expr_subset))))

cat("==== 6. Aggregate by sex, diagnosis, and region ====\n")
meta$Group <- paste(meta$sex, meta$diagnosis, meta$region, sep = "_")
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
  Diagnosis = sapply(col_anno, `[`, 2),
  Region = sapply(col_anno, `[`, 3)
)

# Better color assignment for categorical variables
sex_col <- setNames(c('black', 'grey'), unique(anno_df$Sex))

# Distinct colors for diagnosis
diagnosis_values <- unique(anno_df$Diagnosis)
diagnosis_col <- setNames(
  c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33")[1:length(diagnosis_values)],
  diagnosis_values
)

# Different distinct colors for region
region_values <- unique(anno_df$Region)
region_col <- setNames(
  c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462")[1:length(region_values)],
  region_values
)

column_ha <- HeatmapAnnotation(
  Sex = anno_df$Sex,
  Diagnosis = anno_df$Diagnosis,
  Region = anno_df$Region,
  col = list(
    Sex = sex_col,
    Diagnosis = diagnosis_col,
    Region = region_col
  )
)

cat("==== 8. Plot heatmap ====\n")
dir.create('GSE174409/figures', showWarnings = FALSE, recursive = TRUE)
pdf('GSE174409/figures/heatmap_output.pdf', width = 5, height = 10)
Heatmap(
  agg_expr_mat,
  name = "Z-score",
  top_annotation = column_ha,
  col = colorRamp2(seq(-1, 1, length.out = 100), colorRampPalette(c("#4575B4", "white", "#D73027"))(100)),
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_columns = FALSE,
  cluster_rows = TRUE
)
dev.off()
cat("==== DONE ====\n")