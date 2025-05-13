library(ComplexHeatmap)
library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(circlize)
library(viridis)

# Set base directory for all files
base_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD"

# Load gene list
gene_list <- read_csv(file.path(base_dir, 'GSE225158/KEGG outputs/kegg_complement_unique_genes.csv'))
selected_genes <- gene_list$gene
cat("Complement genes list (first 10):", head(selected_genes, 10), "\n")

# Create output folder
heatmap_output_dir <- file.path(base_dir, 'GSE207128/heatmap_outputs')
dir.create(heatmap_output_dir, showWarnings = FALSE, recursive = TRUE)

# Helper function to convert gene IDs
process_gene_symbols <- function(expr, gene_mapping) {
  if(all(grepl("^ENSMUSG", rownames(expr)))) {
    gene_dict <- setNames(gene_mapping$gene_symbol, gene_mapping$ensembl_id)
    gene_symbols <- gene_dict[rownames(expr)]
    missing_idx <- is.na(gene_symbols)
    if(any(missing_idx)) gene_symbols[missing_idx] <- rownames(expr)[missing_idx]
    rownames(expr) <- make.unique(gene_symbols, sep="_")
  }
  return(expr)
}

# Helper function to find complement genes
find_complement_genes <- function(expr, selected_genes) {
  normalized_gene_list <- toupper(selected_genes)
  original_rownames <- rownames(expr)
  normalized_rownames <- toupper(original_rownames)
  clean_rownames <- sub("_\\d+$", "", normalized_rownames)

  direct_matches <- original_rownames[normalized_rownames %in% normalized_gene_list]
  if(length(direct_matches) == 0) {
    suffix_matches <- original_rownames[clean_rownames %in% normalized_gene_list]
    if(length(suffix_matches) == 0) {
      pattern <- paste0("^(", paste(normalized_gene_list, collapse="|"), ")")
      partial_matches <- original_rownames[grepl(pattern, normalized_rownames, ignore.case=TRUE)]

      if(length(partial_matches) == 0) {
        cat("No complement genes found. Using all genes.\n")
        return(expr)
      } else {
        cat("Found", length(partial_matches), "complement genes using partial matching\n")
        return(expr[partial_matches, ])
      }
    } else {
      cat("Found", length(suffix_matches), "complement genes after removing suffixes\n")
      return(expr[suffix_matches, ])
    }
  } else {
    cat("Found", length(direct_matches), "complement genes with direct matching\n")
    cat("First few matched genes:", head(direct_matches), "\n")
    return(expr[direct_matches, ])
  }
}

# Create heatmap function
create_heatmap <- function(expr_mat, column_ha, output_file, title) {
  # Check if matrix has data
  if(nrow(expr_mat) == 0 || ncol(expr_mat) == 0) {
    cat("ERROR: Empty matrix for", output_file, "\n")
    return(FALSE)
  }

  # Check for NA or infinite values
  if(any(is.na(expr_mat)) || any(is.infinite(as.matrix(expr_mat)))) {
    cat("Warning: Matrix contains NA or infinite values. Replacing with zeros.\n")
    expr_mat[is.na(expr_mat) | is.infinite(expr_mat)] <- 0
  }

  cat("Creating", output_file, "with dimensions:", nrow(expr_mat), "x", ncol(expr_mat), "\n")

  # Create PDF
  pdf_file <- file.path(heatmap_output_dir, output_file)
  pdf(pdf_file, width = 12, height = 12)

  # Create heatmap
  ht <- Heatmap(
    expr_mat,
    name = "Z-score",
    top_annotation = column_ha,
    col = colorRamp2(seq(-2, 2, length.out = 100), colorRampPalette(c("#50c878", "white", "#FF69B4"))(100)),
    show_row_names = TRUE,
    show_column_names = FALSE,
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_title = title,
    column_title_gp = gpar(fontsize = 12)
  )

  # Draw heatmap
  tryCatch({
    draw(ht)
    invisible(dev.off())
    cat("Heatmap saved to:", pdf_file, "\n")
    if(file.exists(pdf_file) && file.size(pdf_file) > 1000) {
      cat("File created successfully with size:", file.size(pdf_file), "bytes\n")
      return(TRUE)
    } else {
      cat("WARNING: File may be empty or too small:", file.size(pdf_file), "bytes\n")
      return(FALSE)
    }
  }, error = function(e) {
    dev.off()
    cat("ERROR creating heatmap:", e$message, "\n")
    return(FALSE)
  })
}

# Function for cell type processing
process_cell_group <- function(meta_filtered, expr_subset, group_name) {
  if(nrow(meta_filtered) == 0) {
    cat("No", group_name, "cell types found in the data\n")
    return()
  }

  cat("Processing", group_name, "cells. Found", nrow(meta_filtered), "samples\n")

  expr_filtered <- expr_subset[, meta_filtered$group]
  cat("Expression matrix dimensions:", nrow(expr_filtered), "genes x", ncol(expr_filtered), "samples\n")

  # Z-score normalization
  expr_scaled <- t(scale(t(as.matrix(expr_filtered))))
  cat("After scaling, NA values:", sum(is.na(expr_scaled)),
      "Inf values:", sum(is.infinite(expr_scaled)), "\n")

  # Replace NA/Inf values if any exist
  if(any(is.na(expr_scaled)) || any(is.infinite(expr_scaled))) {
    expr_scaled[is.na(expr_scaled) | is.infinite(expr_scaled)] <- 0
  }

  # Create aggregated expression matrix
  agg_expr <- meta_filtered %>%
    select(Sample = group, Group, CellType = celltype_simplified) %>%
    left_join(
      as.data.frame(t(expr_scaled)) %>% rownames_to_column("Sample"),
      by = "Sample"
    ) %>%
    group_by(Group, CellType) %>%
    summarise(across(where(is.numeric), mean), .groups = "drop")

  # Reshape to matrix
  agg_expr_long <- agg_expr %>%
    pivot_longer(-c(Group, CellType), names_to = "Gene", values_to = "Value") %>%
    unite("Group_CellType", Group, CellType, sep = "_") %>%
    pivot_wider(names_from = Group_CellType, values_from = Value)

  agg_expr_mat <- as.matrix(agg_expr_long[,-1])
  rownames(agg_expr_mat) <- agg_expr_long$Gene

  cat("Final aggregated matrix dimensions:", nrow(agg_expr_mat), "x", ncol(agg_expr_mat), "\n")

  # Count frequencies
  freqs <- meta_filtered %>%
    count(Group, celltype_simplified) %>%
    unite("Group_CellType", Group, celltype_simplified, sep = "_")
  col_freqs <- freqs$n[match(colnames(agg_expr_mat), freqs$Group_CellType)]
  col_props <- col_freqs / sum(col_freqs)

  # Prepare column annotation
  col_anno <- strsplit(colnames(agg_expr_mat), "_")
  anno_df <- data.frame(
    Condition = sapply(col_anno, `[`, 1),
    CellType = sapply(col_anno, `[`, 2)
  )

  celltypes_filtered <- unique(meta_filtered$celltype_simplified)
  palette_func <- colorRampPalette(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"))
  celltype_colors <- setNames(palette_func(length(celltypes_filtered)), celltypes_filtered)

  column_ha <- HeatmapAnnotation(
    Proportion = anno_barplot(col_props, border = FALSE, gp = gpar(fill = "grey")),
    Condition = anno_df$Condition,
    CellType = anno_df$CellType,
    col = list(
      Condition = diagnosis_col,
      CellType = celltype_colors
    )
  )

  # Create heatmap
  create_heatmap(
    agg_expr_mat,
    column_ha,
    paste0('complement_heatmap_', tolower(group_name), '.pdf'),
    paste(group_name, "Cells")
  )
}

# Load gene mapping
gene_mapping_path <- file.path(base_dir, "GSE207128/pseudobulk_results/GSE207128_pseudobulk_gene_mapping.csv")
gene_mapping <- read_csv(gene_mapping_path, col_names = FALSE)
colnames(gene_mapping) <- c("ensembl_id", "gene_symbol")
gene_mapping <- gene_mapping %>%
  filter(!is.na(gene_symbol) & gene_symbol != "") %>%
  distinct(ensembl_id, .keep_all = TRUE)

# Define cell type groups
glial_cells <- c("Astrocytes", "Microglia", "Oligodendrocytes")
immune_cells <- c("Macrophages", "T cells", "Monocytes", "B cells", "Granulocytes", "NK cells")

# Color palettes
diagnosis_col <- setNames(c('#a6cee3', '#b2df8a'), c('OpioidDependence', 'Normal'))

# Load the data
cat("\n========== Processing celltype_by_treatment dataset ==========\n")
pseudobulk_path <- file.path(base_dir, "GSE207128/pseudobulk_results/GSE207128_pseudobulk_celltype_by_treatment.rds")
pseudobulk_data <- readRDS(pseudobulk_path)
expr <- as.data.frame(pseudobulk_data)
cat("Original dimensions:", nrow(expr), "genes x", ncol(expr), "samples\n")

# Process gene symbols
expr <- process_gene_symbols(expr, gene_mapping)
expr_subset <- find_complement_genes(expr, selected_genes)
cat("After filtering for complement genes:", nrow(expr_subset), "genes x", ncol(expr_subset), "samples\n")

# Create metadata
meta <- data.frame(group = colnames(expr_subset)) %>%
  mutate(
    celltype = gsub("^(.+)_(.+)$", "\\1", group),
    Dx_OUD = gsub("^(.+)_(.+)$", "\\2", group),
    celltype_simplified = celltype,
    Group = Dx_OUD
  )

cat("Cell types found:", paste(unique(meta$celltype), collapse=", "), "\n")

# Create cell type specific heatmaps
process_cell_group(meta, expr_subset, "All")
process_cell_group(meta %>% filter(celltype %in% glial_cells), expr_subset, "Glial")
process_cell_group(meta %>% filter(celltype %in% immune_cells), expr_subset, "Immune")

cat("\nAll processing complete. Check files in:", heatmap_output_dir, "\n")