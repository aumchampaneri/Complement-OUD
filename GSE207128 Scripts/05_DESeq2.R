# DESeq2 Analysis for GSE207128 Pseudobulk Data

library(DESeq2)
library(dplyr)
library(tibble)
library(readr)
library(ggplot2)
library(pheatmap)

# Set base directory
base_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD"
output_dir <- file.path(base_dir, "GSE207128/DESeq2_outputs")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load pseudobulk data (using the celltype_by_treatment version)
pseudobulk_path <- file.path(base_dir, "GSE207128/pseudobulk_results/GSE207128_pseudobulk_celltype_by_treatment.rds")
pseudobulk_data <- readRDS(pseudobulk_path)
counts <- as.matrix(pseudobulk_data)

# Add a small pseudocount to avoid zeros (helps with size factor estimation)
counts <- counts + 1

# Create metadata from column names
sample_names <- colnames(counts)
meta <- data.frame(
  sample_id = sample_names,
  celltype = gsub("^(.+)_(.+)$", "\\1", sample_names),
  condition = gsub("^(.+)_(.+)$", "\\2", sample_names),
  row.names = sample_names
)

# More stringent pre-filtering to remove genes with too many zeros
min_samples <- ceiling(ncol(counts) * 0.25)
keep <- rowSums(counts >= 10) >= min_samples
counts_filtered <- counts[keep, ]
cat("Keeping", sum(keep), "out of", nrow(counts), "genes after filtering\n")

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_filtered),
  colData = meta,
  design = ~ condition
)

# Use alternative size factor estimation for sparse data
dds <- estimateSizeFactors(dds, type="poscounts")

# Run DESeq2
dds <- DESeq(dds)

# Get results for OpioidDependence vs Normal
res_all <- results(dds, contrast = c("condition", "OpioidDependence", "Normal"))
res_all_df <- as.data.frame(res_all) %>%
  rownames_to_column("gene_id")

# Save results
write.csv(res_all_df, file.path(output_dir, "deseq2_all_cells_OUD_vs_Normal.csv"), row.names = FALSE)

# Run DESeq2 for each cell type separately
cell_types <- unique(meta$celltype)
cell_type_results <- list()

for (cell in cell_types) {
  cat("Processing", cell, "\n")

  # Subset data for this cell type
  cell_samples <- rownames(meta)[meta$celltype == cell]
  if (length(cell_samples) < 2) {
    cat("Skipping", cell, "- not enough samples\n")
    next
  }

  # Check if we have samples from both conditions
  conditions_present <- unique(meta$condition[meta$celltype == cell])
  if (length(conditions_present) < 2) {
    cat("Skipping", cell, "- only one condition present\n")
    next
  }

  # Subset the already filtered counts
  cell_counts <- counts_filtered[, cell_samples]
  cell_meta <- meta[cell_samples, ]

  # Create DESeq2 object
  cell_dds <- DESeqDataSetFromMatrix(
    countData = round(cell_counts),
    colData = cell_meta,
    design = ~ condition
  )

  # Use alternative size factor estimation
  cell_dds <- estimateSizeFactors(cell_dds, type="poscounts")

  # Run DESeq2
  cell_dds <- tryCatch({
    DESeq(cell_dds)
  }, error = function(e) {
    cat("Error in DESeq for", cell, ":", e$message, "\n")
    return(NULL)
  })

  if (is.null(cell_dds)) next

  # Get results
  cell_res <- results(cell_dds, contrast = c("condition", "OpioidDependence", "Normal"))
  cell_res_df <- as.data.frame(cell_res) %>%
    rownames_to_column("gene_id")

  # Save results
  write.csv(cell_res_df, file.path(output_dir, paste0("deseq2_", cell, "_OUD_vs_Normal.csv")), row.names = FALSE)

  # Store for summary
  cell_type_results[[cell]] <- cell_res_df
}

# Create a summary of DE genes across cell types
create_volcano_plot <- function(res_df, title, file_name) {
  # Handle NA values in padj
  res_df <- res_df %>%
    mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
    mutate(sig = ifelse(padj < 0.05,
                        ifelse(log2FoldChange > 0.5, "Up",
                               ifelse(log2FoldChange < -0.5, "Down", "NS")),
                        "NS"))

  # Create volcano plot
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = sig)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("Down" = "blue", "NS" = "gray", "Up" = "red")) +
    theme_minimal() +
    labs(title = title,
         x = "Log2 Fold Change",
         y = "-Log10 P-value") +
    theme(legend.title = element_blank())

  # Save plot
  ggsave(file.path(output_dir, file_name), p, width = 8, height = 6)
}

# Create volcano plot for all cells
create_volcano_plot(res_all_df, "OUD vs Normal - All Cells", "volcano_all_cells.png")

# Create volcano plots for each cell type
for (cell in names(cell_type_results)) {
  create_volcano_plot(
    cell_type_results[[cell]],
    paste0("OUD vs Normal - ", cell),
    paste0("volcano_", cell, ".png")
  )
}

# Examine complement pathway genes specifically
if (file.exists(file.path(base_dir, 'GSE225158/KEGG outputs/kegg_complement_unique_genes.csv'))) {
  complement_genes <- read_csv(file.path(base_dir, 'GSE225158/KEGG outputs/kegg_complement_unique_genes.csv'))$gene

  # Extract results for complement genes from all-cell analysis
  complement_results <- res_all_df %>%
    filter(toupper(gene_id) %in% toupper(complement_genes))

  # Save complement-specific results
  write.csv(complement_results, file.path(output_dir, "deseq2_complement_genes_all_cells.csv"), row.names = FALSE)

  # Create heatmap of log2FC for complement genes across cell types
  if (length(cell_type_results) > 0) {
    # Create a matrix of log2FC values
    log2fc_matrix <- matrix(NA, nrow = length(complement_genes), ncol = length(cell_type_results))
    rownames(log2fc_matrix) <- complement_genes
    colnames(log2fc_matrix) <- names(cell_type_results)

    # Create p-value significance matrix
    pval_matrix <- matrix("", nrow = length(complement_genes), ncol = length(cell_type_results))
    rownames(pval_matrix) <- complement_genes
    colnames(pval_matrix) <- names(cell_type_results)

    for (i in seq_along(cell_type_results)) {
      cell <- names(cell_type_results)[i]
      cell_df <- cell_type_results[[cell]]

      for (gene in complement_genes) {
        # Case-insensitive matching
        gene_row <- which(toupper(cell_df$gene_id) == toupper(gene))
        if (length(gene_row) > 0) {
          log2fc_matrix[gene, cell] <- cell_df$log2FoldChange[gene_row[1]]

          # Add significance marker
          if (!is.na(cell_df$padj[gene_row[1]]) && cell_df$padj[gene_row[1]] < 0.05) {
            pval_matrix[gene, cell] <- "*"
          }
        }
      }
    }

    # Handle NAs for visualization
    log2fc_matrix_viz <- log2fc_matrix
    log2fc_matrix_viz[is.na(log2fc_matrix_viz)] <- 0

    # Create heatmap
    pdf(file.path(output_dir, "complement_genes_log2fc_heatmap.pdf"), width = 10, height = 12)
    pheatmap(log2fc_matrix_viz,
             display_numbers = pval_matrix,
             main = "Log2FC of Complement Genes Across Cell Types",
             color = colorRampPalette(c("blue", "white", "red"))(100),
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             fontsize_row = 10,
             fontsize_col = 10,
             angle_col = 45)
    dev.off()
  }
}

cat("DESeq2 analysis complete. Results saved in:", output_dir, "\n")