#===============================================================================
# BULK RNA-SEQ PREPROCESSING & QUALITY CONTROL PIPELINE
#===============================================================================
# 
# This script performs streamlined preprocessing and quality control analysis
# for bulk RNA-seq data using DESeq2 with batch correction.
#
# WORKFLOW:
# 1. Data loading and validation
# 2. Gene ID conversion (Ensembl to symbols)
# 3. DESeq2 filtering and normalization
# 4. Essential QC visualizations
# 5. Batch effect detection and correction
# 6. Export for downstream analysis
#
#===============================================================================

# Load required libraries
library(DESeq2)
library(limma)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Set paths
setwd("/Users/aumchampaneri/Complement-OUD/Multi-Omics Study")
data_dir <- "data/raw/bulkrna"
output_dir <- "data/processed/bulkrna"
plots_dir <- "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study/results/bulkrna"

# Create directories
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

#===============================================================================
# 1. DATA LOADING & VALIDATION
#===============================================================================

# Read count matrix
count_files <- list.files(data_dir, pattern = "*.csv", full.names = TRUE)
count_files <- count_files[!grepl("metadata", count_files, ignore.case = TRUE)]
counts <- read.csv(count_files[1], header = TRUE, row.names = 1, check.names = FALSE)

# Clean data
counts <- counts[rowSums(counts, na.rm = TRUE) > 0, ]
cat("Count matrix dimensions:", dim(counts), "\n")

# Convert Ensembl IDs to Gene Symbols
ensembl_pattern <- "^ENS[A-Z]*[0-9]+$|^ENS[A-Z]*[0-9]+\\.[0-9]+$"
if (any(grepl(ensembl_pattern, rownames(counts)[1:10]))) {
  ensembl_ids <- sub("\\.[0-9]+$", "", rownames(counts))
  gene_symbols <- mapIds(org.Hs.eg.db, keys = ensembl_ids, column = "SYMBOL", 
                        keytype = "ENSEMBL", multiVals = "first")
  
  # Create final gene names
  final_names <- ifelse(is.na(gene_symbols), rownames(counts), gene_symbols)
  # Handle duplicates
  final_names <- make.unique(final_names)
  rownames(counts) <- final_names
  
  cat("Successfully converted", sum(!is.na(gene_symbols)), "gene IDs\n")
}

# Read metadata
sample_info <- read.csv(file.path(data_dir, "sample_metadata.csv"), 
                       header = TRUE, row.names = 1)

# Align samples
common_samples <- intersect(colnames(counts), rownames(sample_info))
counts <- counts[, common_samples]
sample_info <- sample_info[common_samples, ]

cat("Final dimensions - Counts:", dim(counts), "Metadata:", dim(sample_info), "\n")

#===============================================================================
# 2. DESEQ2 FILTERING & NORMALIZATION
#===============================================================================

# Create DESeq2 object
sample_info$condition <- factor(sample_info$condition)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = sample_info, 
                             design = ~ condition)

# Filter low-expressed genes
keep_genes <- rowSums(counts(dds) >= 10) >= 3
dds_filtered <- dds[keep_genes, ]
cat("Genes after filtering:", nrow(dds_filtered), "\n")

# Normalize
dds_norm <- DESeq(dds_filtered)
vsd <- vst(dds_norm, blind = FALSE)

#===============================================================================
# 3. ESSENTIAL QC VISUALIZATIONS
#===============================================================================

# PCA plot
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "PCA Plot", 
       x = paste0("PC1: ", round(attr(pca_data, "percentVar")[1] * 100), "% variance"),
       y = paste0("PC2: ", round(attr(pca_data, "percentVar")[2] * 100), "% variance")) +
  theme_minimal()

ggsave(file.path(plots_dir, "pca_plot.png"), p_pca, width = 8, height = 6)

# Sample clustering heatmap
sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

png(file.path(plots_dir, "sample_clustering.png"), width = 800, height = 600)
pheatmap(sample_dist_matrix, col = colors,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists)
dev.off()

#===============================================================================
# 4. BATCH EFFECT ANALYSIS & CORRECTION
#===============================================================================

# Detect batch effects
cat("Analyzing batch effects...\n")
potential_batch_vars <- c("batch", "region", "sex")
available_batch_vars <- potential_batch_vars[potential_batch_vars %in% colnames(sample_info)]

# Create PCA plots colored by batch variables
for (batch_var in available_batch_vars) {
  if (length(unique(sample_info[[batch_var]])) > 1) {
    p_batch <- ggplot(pca_data, aes(x = PC1, y = PC2, 
                                   color = factor(sample_info[[batch_var]]))) +
      geom_point(size = 3, alpha = 0.8) +
      labs(title = paste("PCA colored by", batch_var), 
           x = paste0("PC1: ", round(attr(pca_data, "percentVar")[1] * 100), "% variance"),
           y = paste0("PC2: ", round(attr(pca_data, "percentVar")[2] * 100), "% variance"),
           color = batch_var) +
      theme_minimal()
    
    ggsave(file.path(plots_dir, paste0("pca_", batch_var, ".png")), 
           p_batch, width = 8, height = 6)
  }
}

# Apply batch correction if batch variable exists
if ("batch" %in% available_batch_vars) {
  cat("Applying batch correction...\n")
  
  # DESeq2 batch correction
  sample_info$batch <- factor(sample_info$batch)
  dds_batch <- DESeqDataSetFromMatrix(countData = counts, colData = sample_info,
                                     design = ~ batch + condition)
  keep_genes_batch <- rowSums(counts(dds_batch) >= 10) >= 3
  dds_batch_filtered <- dds_batch[keep_genes_batch, ]
  dds_batch_norm <- DESeq(dds_batch_filtered)
  vsd_batch <- vst(dds_batch_norm, blind = FALSE)
  
  # Limma batch correction for visualization
  design_bio <- model.matrix(~ condition, data = sample_info)
  vst_corrected <- removeBatchEffect(assay(vsd_batch), 
                                   batch = sample_info$batch,
                                   design = design_bio)
  
  # PCA on corrected data
  pca_corrected <- prcomp(t(vst_corrected), scale. = TRUE)
  pca_corrected_data <- data.frame(
    PC1 = pca_corrected$x[, 1],
    PC2 = pca_corrected$x[, 2],
    condition = sample_info$condition,
    batch = sample_info$batch
  )
  
  pc1_var <- round(summary(pca_corrected)$importance[2, 1] * 100, 1)
  pc2_var <- round(summary(pca_corrected)$importance[2, 2] * 100, 1)
  
  # Before/After comparison plots
  p_before <- ggplot(pca_data, aes(x = PC1, y = PC2, color = factor(sample_info$batch))) +
    geom_point(size = 3, alpha = 0.8) +
    labs(title = "Before Batch Correction", color = "Batch") +
    theme_minimal()
  
  p_after <- ggplot(pca_corrected_data, aes(x = PC1, y = PC2, color = batch)) +
    geom_point(size = 3, alpha = 0.8) +
    labs(title = "After Batch Correction", color = "Batch") +
    theme_minimal()
  
  batch_comparison <- grid.arrange(p_before, p_after, ncol = 2,
                                 top = "Batch Effect Correction")
  ggsave(file.path(plots_dir, "batch_correction_comparison.png"), 
         batch_comparison, width = 12, height = 5)
  
  cat("Batch correction completed!\n")
  
} else {
  dds_batch_norm <- dds_norm
  cat("No batch variable found - using standard normalization\n")
}

#===============================================================================
# 5. SAVE DATA FOR DOWNSTREAM ANALYSIS
#===============================================================================

save(dds_batch_norm, sample_info, 
     file = file.path(output_dir, "deseq2_for_DE_analysis.RData"))

cat("\nPreprocessing completed successfully!\n")
cat("Outputs:\n")
cat("- DESeq2 object ready for differential expression analysis\n")
cat("- QC and batch correction plots\n")
cat("- Ready for downstream analysis!\n")
