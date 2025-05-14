# GSE207128 - Pseudobulk Heatmap of Complement Genes

library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(RColorBrewer)
library(tibble)

#----- 1. Set up paths and directories -----#
base_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD"
heatmap_output_dir <- file.path(base_dir, "GSE207128/figures/heatmaps")
dir.create(heatmap_output_dir, showWarnings = FALSE, recursive = TRUE)

# Set working directory to save PDFs
setwd(heatmap_output_dir)

#----- 2. Define complement genes -----#
# List of complement-related genes to visualize
selected_genes <- c(
  "C1QA", "C1QB", "C1QC", "C1R", "C1S",
  "C2", "C3", "C4A", "C4B", "C5", "C6", "C7", "C8A", "C8B", "C8G", "C9",
  "CFB", "CFD", "CFH", "CFI", "CFP",
  "CR1", "CR2", "C3AR1", "C5AR1", "C5AR2",
  "ITGAM", "ITGAX", "VSIG4",
  "MASP1", "MASP2", "FCN1", "FCN2", "FCN3",
  "CD55", "CD46", "CD59", "SERPING1", "ADIPOQ"
)

#----- 3. Define utility functions -----#
# Function to process gene symbols in expression data
process_gene_symbols <- function(expr_data, gene_mapping) {
  cat("Sample rownames:", head(rownames(expr_data)), "\n")

  # Check if rownames are already gene symbols
  if(all(grepl("^[A-Za-z]", rownames(expr_data)))) {
    cat("Expression data already using gene symbols\n")
    return(expr_data)
  }

  # Map Ensembl IDs to gene symbols
  new_rownames <- gene_mapping$gene_symbol[match(rownames(expr_data), gene_mapping$ensembl_id)]

  # Filter out rows with missing or duplicate gene symbols
  has_symbol <- !is.na(new_rownames) & new_rownames != ""
  expr_data <- expr_data[has_symbol, ]
  new_rownames <- new_rownames[has_symbol]

  # Handle duplicate gene symbols by keeping the one with highest mean expression
  if(any(duplicated(new_rownames))) {
    cat("Found", sum(duplicated(new_rownames)), "duplicate gene symbols. Keeping highest expressed per gene.\n")
    expr_data <- expr_data %>%
      as.data.frame() %>%
      mutate(gene_symbol = new_rownames,
             mean_expr = rowMeans(select(., -gene_symbol), na.rm = TRUE)) %>%
      group_by(gene_symbol) %>%
      slice_max(order_by = mean_expr, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      select(-mean_expr) %>%
      column_to_rownames("gene_symbol")
  } else {
    rownames(expr_data) <- new_rownames
  }

  cat("Expression data: processed", nrow(expr_data), "genes with symbols\n")
  return(expr_data)
}

# Function to filter for complement genes with case-insensitive matching
find_complement_genes <- function(expr_data, gene_list) {
  if(nrow(expr_data) == 0) {
    cat("Expression matrix is empty, cannot find complement genes\n")
    return(expr_data)
  }

  # Try case-insensitive matching
  genes_in_data <- rownames(expr_data)
  gene_list_upper <- toupper(gene_list)
  genes_in_data_upper <- toupper(genes_in_data)

  idx <- which(genes_in_data_upper %in% gene_list_upper)
  genes_present <- genes_in_data[idx]

  cat("Found", length(genes_present), "complement genes out of", length(gene_list), "desired genes\n")

  if(length(genes_present) == 0) {
    # Try substring matching in case of prefix/suffix differences
    cat("No exact matches found. Looking for partial matches...\n")
    partial_matches <- character(0)
    for(gene in gene_list) {
      pattern <- gene
      matches <- grep(pattern, genes_in_data, ignore.case = TRUE, value = TRUE)
      if(length(matches) > 0) {
        partial_matches <- c(partial_matches, matches)
      }
    }
    genes_present <- unique(partial_matches)
    cat("Found", length(genes_present), "genes with partial matches\n")
  }

  if(length(genes_present) == 0) {
    # Fall back to top expressed genes if no complement genes found
    cat("No complement genes found. Using top 20 expressed genes instead.\n")
    mean_expr <- rowMeans(expr_data)
    genes_present <- names(sort(mean_expr, decreasing = TRUE))[1:min(20, nrow(expr_data))]
  } else {
    cat("Missing genes:", paste(setdiff(gene_list, genes_present), collapse=", "), "\n")
  }

  # Return subset of expression data with only selected genes
  return(expr_data[genes_present, , drop=FALSE])
}

# Function to create and save heatmap
create_heatmap <- function(expr_mat, column_ha, pdf_file, title) {
  if(nrow(expr_mat) == 0) {
    cat("No genes to plot in heatmap:", title, "\n")
    return(FALSE)
  }

  cat("Creating heatmap:", title, "\n")
  pdf(pdf_file, width = 10, height = max(8, nrow(expr_mat) * 0.25))

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
    return(TRUE)
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

  if(nrow(expr_subset) == 0) {
    cat("No genes to process for", group_name, "cells\n")
    return()
  }

  cat("Processing", group_name, "cells. Found", nrow(meta_filtered), "samples\n")

  # Make sure all sample names in metadata exist in expression data
  valid_samples <- intersect(meta_filtered$sample_name, colnames(expr_subset))
  if(length(valid_samples) == 0) {
    cat("No matching samples found between metadata and expression data for", group_name, "cells\n")
    return()
  }

  meta_filtered <- meta_filtered %>% filter(sample_name %in% valid_samples)

  expr_filtered <- expr_subset[, meta_filtered$sample_name]
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
    select(Sample = sample_name, Group = derived_condition, CellType = derived_celltype) %>%
    left_join(
      as.data.frame(t(expr_scaled)) %>% rownames_to_column("Sample"),
      by = "Sample"
    ) %>%
    group_by(Group, CellType) %>%
    summarise(across(where(is.numeric), mean), .groups = "drop")

  # Get raw counts for each group/celltype combination
  raw_counts <- meta_filtered %>%
    group_by(derived_condition, derived_celltype) %>%
    summarise(Count = n(), .groups = "drop") %>%
    unite("Group_CellType", derived_condition, derived_celltype, sep = "_")

  # Reshape the aggregated expression to matrix format
  agg_expr_long <- agg_expr %>%
    pivot_longer(cols = -c(Group, CellType), names_to = "Gene", values_to = "Value") %>%
    unite("Group_CellType", Group, CellType, sep = "_") %>%
    pivot_wider(names_from = Group_CellType, values_from = Value)

  agg_expr_mat <- as.matrix(agg_expr_long[,-1])
  rownames(agg_expr_mat) <- agg_expr_long$Gene

  # Set count values to display in the barplot
  column_counts <- numeric(ncol(agg_expr_mat))
  for(i in 1:ncol(agg_expr_mat)) {
    idx <- which(raw_counts$Group_CellType == colnames(agg_expr_mat)[i])
    if(length(idx) > 0) {
      column_counts[i] <- raw_counts$Count[idx]
    }
  }

  cat("Column counts:", paste(column_counts, collapse=", "), "\n")

  # Prepare column annotation
  col_anno <- strsplit(colnames(agg_expr_mat), "_")
  anno_df <- data.frame(
    Condition = sapply(col_anno, `[`, 1),
    CellType = sapply(col_anno, `[`, 2)
  )

  # Create color palettes
  celltypes_filtered <- unique(meta_filtered$derived_celltype)
  palette_func <- colorRampPalette(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"))
  celltype_colors <- setNames(palette_func(length(celltypes_filtered)), celltypes_filtered)

  conditions_filtered <- unique(meta_filtered$derived_condition)
  condition_colors <- setNames(
    c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c")[1:length(conditions_filtered)],
    conditions_filtered
  )

  column_ha <- HeatmapAnnotation(
    df = data.frame(
      Condition = anno_df$Condition,
      CellType = anno_df$CellType
    ),
    "Sample\nCount" = anno_barplot(
      column_counts,
      gp = gpar(fill = "gray50", col = NA),
      border = FALSE,
      axis_param = list(
        gp = gpar(fontsize = 8),
        side = "left"
      ),
      height = unit(2, "cm")
    ),
    col = list(
      Condition = condition_colors,
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

#----- 4. Load the pseudobulk data and metadata -----#
cat("\n========== Loading pseudobulk data ==========\n")
pseudobulk_path <- file.path(base_dir, "GSE207128/pseudobulk_results/GSE207128_pseudobulk_counts_for_deseq2.rds")
metadata_path <- file.path(base_dir, "GSE207128/pseudobulk_results/GSE207128_pseudobulk_metadata_for_deseq2.rds")

# Load the counts and metadata
expr <- readRDS(pseudobulk_path)
if(is.list(expr) && !is.data.frame(expr)) {
  # If the RDS file contains a list object (like a DESeq2 result)
  cat("Expression data is a list object. Extracting counts...\n")
  if("vst" %in% names(expr)) {
    expr <- assay(expr$vst)
  } else if("rlog" %in% names(expr)) {
    expr <- assay(expr$rlog)
  } else if("normalized_counts" %in% names(expr)) {
    expr <- expr$normalized_counts
  } else if("counts" %in% names(expr)) {
    expr <- expr$counts
  } else {
    expr <- as.matrix(expr[[1]])
  }
}
expr <- as.data.frame(expr)

metadata <- readRDS(metadata_path)
cat("Original dimensions:", nrow(expr), "genes x", ncol(expr), "samples\n")
cat("First few column names:", paste(head(colnames(expr)), collapse=", "), "\n")
cat("First few metadata sample names:", paste(head(metadata$sample_name), collapse=", "), "\n")

# Extract condition and cell type from sample name
metadata <- metadata %>%
  mutate(
    # Extract condition (Normal/Disease) from first letter of sample name
    derived_condition = case_when(
      grepl("^GSM\\d+_N", sample_name) ~ "Normal",
      grepl("^GSM\\d+_D", sample_name) ~ "OpioidDependence",
      TRUE ~ "Unknown"
    ),
    # Extract cell type without the prefix
    derived_celltype = sub("^GSM\\d+_[ND]\\d+_", "", sample_name)
  )

# Load gene mapping
gene_mapping_path <- file.path(base_dir, "GSE207128/pseudobulk_results/GSE207128_pseudobulk_gene_mapping.csv")
if(file.exists(gene_mapping_path)) {
  gene_mapping <- read_csv(gene_mapping_path, col_names = c("ensembl_id", "gene_symbol"))
  gene_mapping <- gene_mapping %>%
    filter(!is.na(gene_symbol) & gene_symbol != "") %>%
    distinct(ensembl_id, .keep_all = TRUE)
} else {
  cat("Gene mapping file not found:", gene_mapping_path, "\n")
  cat("Creating dummy mapping from rownames\n")
  gene_mapping <- data.frame(
    ensembl_id = rownames(expr),
    gene_symbol = rownames(expr)
  )
}

# Process gene symbols
expr <- process_gene_symbols(expr, gene_mapping)
expr_subset <- find_complement_genes(expr, selected_genes)
cat("After filtering for genes:", nrow(expr_subset), "genes x", ncol(expr_subset), "samples\n")

# Define cell type groups
glial_cells <- c("Astrocytes", "Microglia", "Oligodendrocytes")
immune_cells <- c("Macrophages", "T cells", "Monocytes", "B cells", "Granulocytes", "NK cells")

# Check if cell types in metadata match our defined groups
actual_cell_types <- unique(metadata$derived_celltype)
cat("Cell types in metadata:", paste(head(actual_cell_types, 10), collapse=", "), "...\n")

# Adjust cell type lists based on actual data
glial_cells <- intersect(glial_cells, actual_cell_types)
immune_cells <- intersect(immune_cells, actual_cell_types)

cat("Glial cells found:", paste(glial_cells, collapse=", "), "\n")
cat("Immune cells found:", paste(immune_cells, collapse=", "), "\n")

#----- 5. Process data and create heatmaps -----#
cat("\n========== Creating heatmaps ==========\n")

# Create cell type specific heatmaps
process_cell_group(metadata, expr_subset, "All")
if(length(glial_cells) > 0) {
  process_cell_group(metadata %>% filter(derived_celltype %in% glial_cells), expr_subset, "Glial")
}
if(length(immune_cells) > 0) {
  process_cell_group(metadata %>% filter(derived_celltype %in% immune_cells), expr_subset, "Immune")
}

cat("\nAll processing complete. Check files in:", heatmap_output_dir, "\n")