# GSE207128 - DESeq2 Analysis using improved pseudobulk data
# Fixed version addressing model matrix rank issues and package dependencies

library(DESeq2)
library(dplyr)
library(readr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(stringr)
library(tibble)  # Added for rownames_to_column function

# Set base directory
base_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD"
pseudobulk_dir <- file.path(base_dir, "GSE207128/pseudobulk_results")
output_dir <- file.path(base_dir, "GSE207128/DESeq2_outputs")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load the pseudobulk data
counts <- readRDS(file.path(pseudobulk_dir, "GSE207128_pseudobulk_counts_for_deseq2.rds"))
metadata <- readRDS(file.path(pseudobulk_dir, "GSE207128_pseudobulk_metadata_for_deseq2.rds"))

# Initialize diagnostics
cat("Counts matrix dimensions:", dim(counts), "\n")
cat("Metadata rows:", nrow(metadata), "\n")

# Fix the metadata rownames if needed
if ("sample_name" %in% colnames(metadata)) {
  rownames(metadata) <- metadata$sample_name
}

# Fix the condition mapping
metadata$condition <- ifelse(
  grepl("^GSM627882[2-4]", metadata$sample_id) | grepl("^D[1-3]_", metadata$celltype),
  "OpioidDependence",
  ifelse(grepl("^GSM627881[9-9]|^GSM627882[0-1]", metadata$sample_id) | grepl("^N[1-3]_", metadata$celltype),
         "Normal", NA)
)

print(table(metadata$condition, useNA = "ifany"))

# Only keep samples with valid condition assignment
metadata <- metadata[!is.na(metadata$condition), ]

# Extract the actual cell type without sample prefix
metadata$celltype_clean <- str_replace(metadata$celltype, "^[ND][1-3]_", "")
cat("\nCleaned cell types:", paste(unique(metadata$celltype_clean)[1:5], "..."), "\n")

# Remove adipocytes as requested (not essential for analysis)
metadata <- metadata[metadata$celltype_clean != "Adipocytes", ]
cat("\nRemoved adipocytes from the analysis\n")

# Ensure count matrix and metadata match
common_samples <- intersect(colnames(counts), rownames(metadata))
counts <- counts[, common_samples]
metadata <- metadata[common_samples, ]

# Print dataset info after cleanup
cat("Total samples after filtering:", nrow(metadata), "\n")
cat("Samples per condition:", table(metadata$condition), "\n")
cat("Cell types (after cleaning):", length(unique(metadata$celltype_clean)), "\n")
cat("Samples per cell type & condition:\n")
print(table(metadata$celltype_clean, metadata$condition))

# Setup DESeq2 dataset with cleaned cell types
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ celltype_clean + condition
)

# Pre-filter genes with at least 10 counts total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat("Kept", sum(keep), "of", length(keep), "genes after filtering\n")

# Run DESeq2
dds <- DESeq(dds)

# Get overall results
res_overall <- results(dds, contrast=c("condition", "OpioidDependence", "Normal"),
                      alpha = 0.05)

# Create a more complete result dataframe
res_df_overall <- as.data.frame(res_overall) %>%
  rownames_to_column("gene_id") %>%
  arrange(padj)

# Save overall results
write_csv(res_df_overall, file.path(output_dir, "deseq2_results_overall_OUD_vs_Normal.csv"))

# Summarize results
cat("\n=== Overall Analysis Results ===\n")
cat("Total tested genes:", nrow(res_df_overall), "\n")
cat("Significant genes (padj < 0.05):", sum(res_df_overall$padj < 0.05, na.rm=TRUE), "\n")
cat("Upregulated in OUD:", sum(res_df_overall$padj < 0.05 & res_df_overall$log2FoldChange > 0, na.rm=TRUE), "\n")
cat("Downregulated in OUD:", sum(res_df_overall$padj < 0.05 & res_df_overall$log2FoldChange < 0, na.rm=TRUE), "\n")

# Create MA plot
pdf(file.path(output_dir, "MA_plot_overall.pdf"), width=8, height=6)
plotMA(res_overall, main="MA Plot - Overall Analysis")
dev.off()

# Analyze each cell type separately
celltypes <- unique(metadata$celltype_clean)
significant_by_celltype <- data.frame(
  celltype = character(),
  total_DEGs = integer(),
  up = integer(),
  down = integer()
)

for (ct in celltypes) {
  # Get samples for this cell type
  ct_samples <- rownames(metadata)[metadata$celltype_clean == ct]
  if (length(ct_samples) < 3) {
    cat("Skipping", ct, "- too few samples (", length(ct_samples), ")\n")
    next
  }

  # Check if we have both conditions
  ct_meta <- metadata[ct_samples, ]
  if (length(unique(ct_meta$condition)) < 2) {
    cat("Skipping", ct, "- only one condition available\n")
    next
  }

  cat("\nAnalyzing", ct, "...\n")

  # Create subset counts and metadata
  ct_counts <- counts[, ct_samples]

  # Create DESeq2 dataset
  ct_dds <- DESeqDataSetFromMatrix(
    countData = ct_counts,
    colData = ct_meta,
    design = ~ condition
  )

  # Pre-filter for this cell type
  keep_ct <- rowSums(counts(ct_dds)) >= 10
  ct_dds <- ct_dds[keep_ct, ]

  # Run DESeq2
  ct_dds <- DESeq(ct_dds)

  # Get results
  ct_res <- results(ct_dds, contrast=c("condition", "OpioidDependence", "Normal"),
                   alpha = 0.05)

  # Create results dataframe
  ct_res_df <- as.data.frame(ct_res) %>%
    rownames_to_column("gene_id") %>%
    arrange(padj)

  # Save results
  out_path <- file.path(output_dir, paste0("deseq2_results_", gsub(" ", "_", ct), "_OUD_vs_Normal.csv"))
  write_csv(ct_res_df, out_path)

  # Count significant genes
  ct_sig <- sum(ct_res_df$padj < 0.05, na.rm=TRUE)
  ct_up <- sum(ct_res_df$padj < 0.05 & ct_res_df$log2FoldChange > 0, na.rm=TRUE)
  ct_down <- sum(ct_res_df$padj < 0.05 & ct_res_df$log2FoldChange < 0, na.rm=TRUE)

  significant_by_celltype <- rbind(significant_by_celltype,
                                  data.frame(celltype = ct,
                                            total_DEGs = ct_sig,
                                            up = ct_up,
                                            down = ct_down))

  cat(ct, "- Significant genes:", ct_sig,
      "(Up:", ct_up, "Down:", ct_down, ")\n")

  # Create MA plot for this cell type
  pdf(file.path(output_dir, paste0("MA_plot_", gsub(" ", "_", ct), ".pdf")), width=8, height=6)
  plotMA(ct_res, main=paste0("MA Plot - ", ct))
  dev.off()
}

# Save summary of DEGs by cell type
write_csv(significant_by_celltype, file.path(output_dir, "DEGs_summary_by_celltype.csv"))

# Add gene symbols to all result files
gene_mapping_path <- file.path(pseudobulk_dir, "GSE207128_pseudobulk_gene_mapping.csv")
if (file.exists(gene_mapping_path)) {
  gene_mapping <- read_csv(gene_mapping_path)
  colnames(gene_mapping) <- c("gene_id", "gene_symbol")

  # Process all result files
  result_files <- list.files(output_dir, pattern="deseq2_results_.*\\.csv", full.names=TRUE)
  for (res_file in result_files) {
    res_df <- read_csv(res_file)
    res_df <- left_join(res_df, gene_mapping, by="gene_id")

    # Reorder columns to have gene_symbol near the beginning
    if ("gene_symbol" %in% colnames(res_df)) {
      first_cols <- c("gene_id", "gene_symbol")
      other_cols <- setdiff(colnames(res_df), first_cols)
      res_df <- res_df[, c(first_cols, other_cols)]
    }

    write_csv(res_df, res_file)
  }
}

cat("\nDESeq2 analysis complete. Results saved in:", output_dir, "\n")