# DESeq2 Analysis for GSE207128 using gene-wise dispersion estimates

library(DESeq2)
library(dplyr)
library(readr)

# Set base directory
base_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD"
output_dir <- file.path(base_dir, "GSE207128/DESeq2_outputs")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load pseudobulk data and ensure integer counts
pseudobulk_path <- file.path(base_dir, "GSE207128/pseudobulk_results/GSE207128_pseudobulk_celltype_by_treatment.rds")
raw_counts <- readRDS(pseudobulk_path)
counts <- round(as.matrix(raw_counts))  # Round to integers for DESeq2

# Create metadata from column names
sample_names <- colnames(counts)
meta <- data.frame(
  row.names = sample_names,
  celltype = gsub("^(.+)_(.+)$", "\\1", sample_names),
  condition = gsub("^(.+)_(.+)$", "\\2", sample_names)
)

# Print dataset info
cat("Total samples:", nrow(meta), "\n")
cat("Cell types:", paste(unique(meta$celltype), collapse=", "), "\n")
cat("Conditions:", paste(unique(meta$condition), collapse=", "), "\n")

# Set up DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta,
  design = ~ celltype + condition
)

# Handle dispersion estimates separately as recommended in error message
dds <- estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds)
dispersions(dds) <- mcols(dds)$dispGeneEst
dds <- nbinomWaldTest(dds)

# Get overall results
res_overall <- results(dds, contrast=c("condition", "OpioidDependence", "Normal"))
res_df_overall <- as.data.frame(res_overall)
res_df_overall$gene <- rownames(res_df_overall)

# Save overall results
write.csv(res_df_overall, file.path(output_dir, "deseq2_results_Overall_OUD_vs_Normal.csv"), row.names=FALSE)
cat("Overall_OUD_vs_Normal - Significant genes (padj < 0.05):", sum(res_df_overall$padj < 0.05, na.rm=TRUE), "\n")

# Analyze each cell type separately
cell_types <- unique(meta$celltype)
for (ct in cell_types) {
  ct_samples <- rownames(meta)[meta$celltype == ct]
  ct_meta <- meta[ct_samples, ]

  # Skip if only one condition is available for this cell type
  if (length(unique(ct_meta$condition)) < 2) {
    cat("Skipping", ct, "- only one condition available\n")
    next
  }

  ct_counts <- counts[, ct_samples]

  # Create cell-type specific DESeq dataset
  ct_dds <- DESeqDataSetFromMatrix(
    countData = ct_counts,
    colData = ct_meta,
    design = ~ condition
  )

  # Use gene-wise dispersions for cell-type specific analysis too
  ct_dds <- estimateSizeFactors(ct_dds)
  ct_dds <- estimateDispersionsGeneEst(ct_dds)
  dispersions(ct_dds) <- mcols(ct_dds)$dispGeneEst
  ct_dds <- nbinomWaldTest(ct_dds)

  # Get results
  ct_res <- results(ct_dds, contrast=c("condition", "OpioidDependence", "Normal"))
  ct_res_df <- as.data.frame(ct_res)
  ct_res_df$gene <- rownames(ct_res_df)

  out_path <- file.path(output_dir, paste0("deseq2_results_", ct, "_OUD_vs_Normal.csv"))
  write.csv(ct_res_df, file = out_path, row.names = FALSE)

  cat(ct, "_OUD_vs_Normal - Significant genes (padj < 0.05):", sum(ct_res_df$padj < 0.05, na.rm=TRUE), "\n")
}

# Export normalized counts
norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, file = file.path(output_dir, "normalized_counts.csv"))

# Add gene symbols if mapping available
gene_mapping_path <- file.path(base_dir, "GSE207128/pseudobulk_results/GSE207128_pseudobulk_gene_mapping.csv")
if (file.exists(gene_mapping_path)) {
  gene_mapping <- read_csv(gene_mapping_path, col_names = FALSE)
  colnames(gene_mapping) <- c("ensembl_id", "gene_symbol")

  # Process all result files
  result_files <- list.files(output_dir, pattern="deseq2_results_.*\\.csv", full.names=TRUE)
  for (res_file in result_files) {
    res_df <- read_csv(res_file)
    gene_id_to_symbol <- setNames(gene_mapping$gene_symbol, gene_mapping$ensembl_id)
    res_df$gene_symbol <- gene_id_to_symbol[res_df$gene]
    write.csv(res_df, res_file, row.names = FALSE)
  }
}

cat("DESeq2 analysis complete. Results saved in:", output_dir, "\n")