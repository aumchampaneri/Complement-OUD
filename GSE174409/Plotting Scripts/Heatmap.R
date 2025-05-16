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
# Better color assignment for categorical variables with more contrasting colors
library(RColorBrewer)

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
cat("Number of complement genes that couldn't be mapped:", sum(is.na(symbol_to_ensembl[comp_genes])), "\n")
if(sum(is.na(symbol_to_ensembl[comp_genes])) > 0) {
  cat("Unmapped genes:", paste(comp_genes[is.na(symbol_to_ensembl[comp_genes])], collapse=", "), "\n")
}

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
unmapped_count <- sum(is.na(gene_symbols) | gene_symbols == "")
if (unmapped_count > 0) {
  cat("Warning:", unmapped_count, "genes couldn't be mapped to gene symbols, using Ensembl IDs instead\n")
}
gene_symbols[is.na(gene_symbols) | gene_symbols == ""] <- rownames(agg_expr_mat)[is.na(gene_symbols) | gene_symbols == ""]
rownames(agg_expr_mat) <- gene_symbols

# Optionally order columns
col_order <- c()
# First group by region
regions <- unique(sapply(strsplit(colnames(agg_expr_mat), "_"), `[`, 3))
diagnoses <- unique(sapply(strsplit(colnames(agg_expr_mat), "_"), `[`, 2))
sexes <- unique(sapply(strsplit(colnames(agg_expr_mat), "_"), `[`, 1))

# Create custom column order if needed
# (Uncomment and modify this section if you want a specific order)
# for (region in regions) {
#   for (sex in sexes) {
#     for (diag in diagnoses) {
#       col_name <- paste(sex, diag, region, sep="_")
#       if (col_name %in% colnames(agg_expr_mat)) {
#         col_order <- c(col_order, col_name)
#       }
#     }
#   }
# }
# if (length(col_order) > 0) {
#   agg_expr_mat <- agg_expr_mat[, col_order]
# }

cat("==== 7. Prepare annotation for columns ====\n")
col_anno <- strsplit(colnames(agg_expr_mat), "_")
anno_df <- data.frame(
  Sex = sapply(col_anno, `[`, 1),
  Diagnosis = sapply(col_anno, `[`, 2),
  Region = sapply(col_anno, `[`, 3)
)

# Better color assignment for categorical variables
sex_col <- setNames(c('black', 'grey50'), unique(anno_df$Sex))

# Distinct colors for diagnosis using a colorblind-friendly palette
diagnosis_values <- unique(anno_df$Diagnosis)
diagnosis_col <- setNames(
  brewer.pal(min(9, max(3, length(diagnosis_values))), "Set1")[1:length(diagnosis_values)],
  diagnosis_values
)

# Different distinct colors for region
region_values <- unique(anno_df$Region)
region_col <- setNames(
  brewer.pal(min(8, max(3, length(region_values))), "Dark2")[1:length(region_values)],
  region_values
)

# If more than 8 regions, use a combination approach
if (length(region_values) > 8) {
  region_col <- setNames(
    colorRampPalette(brewer.pal(8, "Paired"))(length(region_values)),
    region_values
  )
}

# Update the legend parameters for vertical list-style format
column_ha <- HeatmapAnnotation(
  Sex = anno_df$Sex,
  Diagnosis = anno_df$Diagnosis,
  Region = anno_df$Region,
  col = list(
    Sex = sex_col,
    Diagnosis = diagnosis_col,
    Region = region_col
  ),
  annotation_legend_param = list(
    Sex = list(
      title = "Sex",
      direction = "vertical",
      title_position = "topcenter",
      legend_width = unit(3, "cm")
    ),
    Diagnosis = list(
      title = "Diagnosis",
      direction = "vertical",
      title_position = "topcenter"
    ),
    Region = list(
      title = "Brain Region",
      direction = "vertical",
      title_position = "topcenter"
    )
  )
)

cat("==== 8. Plot heatmap ====\n")
# Calculate actual z-score range
zmin <- min(agg_expr_mat, na.rm = TRUE)
zmax <- max(agg_expr_mat, na.rm = TRUE)
cat("Z-score range: [", round(zmin, 2), ", ", round(zmax, 2), "]\n", sep="")

# Use wider range for better color representation
zrange <- max(abs(c(zmin, zmax)), na.rm = TRUE)
if (zrange < 1) zrange <- 1  # Minimum range of 1 for visibility

# Group rows if annotation is available (optional)
# This would group genes by function if you had a function annotation
# row_split <- NULL  # Replace with actual grouping if available

dir.create('GSE174409/figures', showWarnings = FALSE, recursive = TRUE)
pdf('GSE174409/figures/heatmap_output.pdf', width = 7, height = 10)

# Create the main heatmap
ht <- Heatmap(
  agg_expr_mat,
  name = "Z-score",
  top_annotation = column_ha,
  col = colorRamp2(c(-zrange, 0, zrange), c("#4575B4", "white", "#D73027")),
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_columns = TRUE,   # Enable column clustering for biological insight
  cluster_rows = TRUE,
  # row_split = row_split,  # Uncomment if you have row grouping
  row_names_gp = gpar(fontsize = 8),
  heatmap_legend_param = list(
    title = "Z-score",
    at = c(-zrange, -zrange/2, 0, zrange/2, zrange),
    labels = c(sprintf("%.1f", -zrange), sprintf("%.1f", -zrange/2), "0",
               sprintf("%.1f", zrange/2), sprintf("%.1f", zrange))
  )
)

# Draw the heatmap with a title
draw(ht, padding = unit(c(2, 2, 10, 2), "mm"),
     column_title = "Complement Cascade Gene Expression",
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     heatmap_legend_side = "right",
     annotation_legend_side = "right")

dev.off()

# Create an additional heatmap with viridis color scheme for colorblind accessibility
pdf('GSE174409/figures/heatmap_output_viridis.pdf', width = 7, height = 10)
ht_viridis <- Heatmap(
  agg_expr_mat,
  name = "Z-score",
  top_annotation = column_ha,
  col = colorRamp2(seq(-zrange, zrange, length.out = 100), viridis(100, option="C")),
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  row_names_gp = gpar(fontsize = 8),
  heatmap_legend_param = list(
    title = "Z-score",
    at = c(-zrange, -zrange/2, 0, zrange/2, zrange),
    labels = c(sprintf("%.1f", -zrange), sprintf("%.1f", -zrange/2), "0",
               sprintf("%.1f", zrange/2), sprintf("%.1f", zrange))
  )
)

draw(ht_viridis, padding = unit(c(2, 2, 10, 2), "mm"),
     column_title = "Complement Cascade Gene Expression (Viridis)",
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     heatmap_legend_side = "right",
     annotation_legend_side = "right")

dev.off()

cat("==== DONE ====\n")