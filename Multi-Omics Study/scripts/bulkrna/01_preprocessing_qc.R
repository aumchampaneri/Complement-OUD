#!/usr/bin/env Rscript

# =============================================================================
# Comprehensive Bulk RNA-seq Analysis using DESeq2
# =============================================================================
# 
# Description: Complete workflow for bulk RNA-seq processing, QC, normalization,
#              and differential expression analysis using DESeq2
# 
# Author: Multi-Omics Study Team
# Date: 2024
# 
# Features:
# - Advanced quality control and outlier detection
# - DESeq2 differential expression analysis
# - Comprehensive visualization and reporting
# - Statistical rigor with multiple testing correction
# - Integration-ready outputs for comparison with snRNA-seq
# 
# =============================================================================

# Load required libraries
required_packages <- c(
  "DESeq2", "dplyr", "tidyr", "readr", "ggplot2", "pheatmap", 
  "RColorBrewer", "viridis", "corrplot", "VennDiagram", "UpSetR",
  "ggrepel", "plotly", "DT", "htmlwidgets", "reshape2", "stringr",
  "limma", "edgeR", "biomaRt", "org.Hs.eg.db", "clusterProfiler",
  "EnhancedVolcano", "ComplexHeatmap", "circlize", "PoiClaClu",
  "RUVSeq", "sva", "FactoMineR", "factoextra", "WGCNA"
)

# Function to install and load packages
load_packages <- function(packages) {
  failed_packages <- c()
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat("Installing package:", pkg, "\n")
      tryCatch({
        if (pkg %in% c("DESeq2", "limma", "edgeR", "biomaRt", "org.Hs.eg.db", 
                       "clusterProfiler", "EnhancedVolcano", "ComplexHeatmap",
                       "RUVSeq")) {
          BiocManager::install(pkg, quiet = TRUE, update = FALSE)
        } else {
          install.packages(pkg, quiet = TRUE)
        }
        library(pkg, character.only = TRUE)
        cat("Successfully loaded:", pkg, "\n")
      }, error = function(e) {
        cat("Failed to install package:", pkg, "-", e$message, "\n")
        failed_packages <<- c(failed_packages, pkg)
      })
    } else {
      cat("Package already loaded:", pkg, "\n")
    }
  }
  
  if (length(failed_packages) > 0) {
    warning("Failed to install packages: ", paste(failed_packages, collapse = ", "))
  }
}

# Install BiocManager if needed
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  library(BiocManager)
}

# Load packages
cat("Loading required packages for bulk RNA-seq analysis...\n")
load_packages(required_packages)

# Suppress dplyr messages
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

# Set up logging
log_file <- "bulkrna_deseq2_analysis.log"
log_messages <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- paste(timestamp, "-", message)
  cat(log_entry, "\n")
  write(log_entry, file = log_file, append = TRUE)
}

log_messages("Starting Comprehensive Bulk RNA-seq Analysis with DESeq2")

# Set seed for reproducibility
set.seed(42)

# Create output directories
output_dirs <- list(
  main = "../../results/bulkrna/deseq2_analysis",
  qc = "../../results/bulkrna/deseq2_analysis/quality_control",
  de = "../../results/bulkrna/deseq2_analysis/differential_expression", 
  plots = "../../results/bulkrna/deseq2_analysis/plots",
  tables = "../../results/bulkrna/deseq2_analysis/tables",
  normalized = "../../results/bulkrna/deseq2_analysis/normalized_data"
)

for (dir in output_dirs) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
}

# Color schemes
color_schemes <- list(
  sex = c("Male" = "#3498db", "Female" = "#e74c3c"),
  condition = c("Control" = "#95a5a6", "OUD" = "#e67e22"),
  brain_region = c("Caudate" = "#9b59b6", "Putamen" = "#1abc9c"),
  batch_colors = RColorBrewer::brewer.pal(8, "Set2")
)

# =============================================================================
# DATA LOADING
# =============================================================================

log_messages("Loading count data and metadata...")

# Define data paths (adjust these to your actual data locations)
count_file <- "../../data/raw/bulkrna/GSE174409_raw_counts_02102020.csv"
metadata_file <- "../../data/raw/bulkrna/sample_metadata.csv"

# Alternative paths to search for data
possible_count_paths <- c(
  "../../data/raw/bulkrna/GSE174409_raw_counts_02102020.csv",
  "../../data/raw/bulkrna/counts_matrix.csv",
  "../../data/raw/bulkrna/gene_counts.csv",
  "../../data/bulkrna/counts_matrix.csv",
  "../../data/counts.csv"
)

possible_metadata_paths <- c(
  "../../data/raw/bulkrna/sample_metadata.csv",
  "../../data/raw/bulkrna/metadata.csv", 
  "../../data/bulkrna/sample_metadata.csv",
  "../../data/metadata.csv",
  "data/sample_metadata.csv"
)

# Find count matrix file
count_file <- NULL
for (path in possible_count_paths) {
  if (file.exists(path)) {
    count_file <- path
    break
  }
}

# Find metadata file
metadata_file <- NULL
for (path in possible_metadata_paths) {
  if (file.exists(path)) {
    metadata_file <- path
    break
  }
}

if (is.null(count_file) || is.null(metadata_file)) {
  cat("Data files not found. Please update the file paths in the script.\n")
  cat("Looking for:\n")
  cat("Count matrix: One of", paste(possible_count_paths, collapse = ", "), "\n")
  cat("Metadata: One of", paste(possible_metadata_paths, collapse = ", "), "\n")
  stop("Required data files not found")
}

# Load count matrix
log_messages(paste("Loading count matrix from:", count_file))
if (str_detect(count_file, "\\.csv$")) {
  counts_raw <- read_csv(count_file, show_col_types = FALSE)
} else {
  counts_raw <- read_delim(count_file, delim = "\t", show_col_types = FALSE)
}

# Load metadata
log_messages(paste("Loading metadata from:", metadata_file))
if (str_detect(metadata_file, "\\.csv$")) {
  metadata <- read_csv(metadata_file, show_col_types = FALSE)
} else {
  metadata <- read_delim(metadata_file, delim = "\t", show_col_types = FALSE)
}

log_messages(paste("Loaded count matrix with", nrow(counts_raw), "genes and", 
                  ncol(counts_raw)-1, "samples"))
log_messages(paste("Loaded metadata for", nrow(metadata), "samples"))

# =============================================================================
# DATA PREPROCESSING
# =============================================================================

log_messages("Preprocessing count data and metadata...")

# Process count matrix
# Check if first column contains gene names (non-numeric)
first_col_sample <- counts_raw[[1]][1:5]  # Check first few values
if (any(is.na(suppressWarnings(as.numeric(first_col_sample))))) {
  # First column is likely gene names
  gene_names <- counts_raw[[1]]
  count_matrix <- counts_raw[, -1] %>% as.matrix()
  rownames(count_matrix) <- gene_names
} else {
  # No gene names column, use all data
  count_matrix <- as.matrix(counts_raw)
  rownames(count_matrix) <- paste0("Gene_", 1:nrow(count_matrix))
}

# Ensure all counts are integers
count_matrix <- round(count_matrix)
count_matrix[count_matrix < 0] <- 0

# Process metadata - the metadata is already properly formatted from our extraction
metadata_clean <- metadata

# Ensure sample_id is character and other factors are properly formatted
metadata_clean$sample_id <- as.character(metadata_clean$sample_id)
metadata_clean$condition <- factor(metadata_clean$condition, levels = c("Control", "OUD"))
metadata_clean$sex <- factor(metadata_clean$sex)
metadata_clean$region <- factor(metadata_clean$region)
metadata_clean$batch <- factor(metadata_clean$batch)

log_messages(paste("Metadata columns:", paste(names(metadata_clean), collapse = ", ")))

# Verify factor levels
log_messages(paste("Condition levels:", paste(levels(metadata_clean$condition), collapse = ", ")))
log_messages(paste("Region levels:", paste(levels(metadata_clean$region), collapse = ", ")))
log_messages(paste("Sex levels:", paste(levels(metadata_clean$sex), collapse = ", ")))

# Match samples between count matrix and metadata
common_samples <- intersect(colnames(count_matrix), metadata_clean$sample_id)
if (length(common_samples) == 0) {
  stop("No matching samples found between count matrix and metadata")
}

# Subset and reorder
count_matrix <- count_matrix[, common_samples]
metadata_clean <- metadata_clean[match(common_samples, metadata_clean$sample_id), ]
rownames(metadata_clean) <- metadata_clean$sample_id

log_messages(paste("Final dataset:", nrow(count_matrix), "genes,", 
                  ncol(count_matrix), "samples"))

# =============================================================================
# QUALITY CONTROL ANALYSIS
# =============================================================================

log_messages("Performing comprehensive quality control analysis...")

# 1. Basic count statistics
count_stats <- data.frame(
  sample_id = colnames(count_matrix),
  total_counts = colSums(count_matrix),
  detected_genes = colSums(count_matrix > 0),
  median_counts = apply(count_matrix, 2, median),
  mean_counts = colMeans(count_matrix)
) %>%
  left_join(metadata_clean, by = "sample_id")

# 2. Gene-level statistics
gene_stats <- data.frame(
  gene_id = rownames(count_matrix),
  total_counts = rowSums(count_matrix),
  n_samples_expressed = rowSums(count_matrix > 0),
  mean_counts = rowMeans(count_matrix),
  variance = apply(count_matrix, 1, var)
)

# 3. Filter low-count genes (expressed in at least 10% of samples)
min_samples <- ceiling(0.1 * ncol(count_matrix))
keep_genes <- gene_stats$n_samples_expressed >= min_samples & gene_stats$mean_counts >= 1

count_matrix_filtered <- count_matrix[keep_genes, ]
gene_stats_filtered <- gene_stats[keep_genes, ]

log_messages(paste("Filtered to", sum(keep_genes), "genes (", 
                  round(100 * sum(keep_genes)/nrow(count_matrix), 1), 
                  "% of total)"))

# =============================================================================
# QUALITY CONTROL VISUALIZATIONS
# =============================================================================

log_messages("Creating quality control visualizations...")

# 1. Library size distribution
p1 <- ggplot(count_stats, aes(x = sample_id, y = total_counts/1e6)) +
  geom_col(aes(fill = condition)) +
  scale_fill_manual(values = color_schemes$condition) +
  labs(title = "Library Sizes", 
       x = "Sample", y = "Total Counts (Millions)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dirs$qc, "library_sizes.pdf"), p1, width = 12, height = 6)
ggsave(file.path(output_dirs$qc, "library_sizes.png"), p1, width = 12, height = 6, dpi = 300)

# 2. Number of detected genes
p2 <- ggplot(count_stats, aes(x = sample_id, y = detected_genes)) +
  geom_col(aes(fill = condition)) +
  scale_fill_manual(values = color_schemes$condition) +
  labs(title = "Number of Detected Genes", 
       x = "Sample", y = "Detected Genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dirs$qc, "detected_genes.pdf"), p2, width = 12, height = 6)
ggsave(file.path(output_dirs$qc, "detected_genes.png"), p2, width = 12, height = 6, dpi = 300)

# 3. Count distribution per sample (boxplot)
count_long <- count_matrix_filtered %>%
  as.data.frame() %>%
  mutate(gene_id = rownames(.)) %>%
  pivot_longer(-gene_id, names_to = "sample_id", values_to = "counts") %>%
  left_join(metadata_clean, by = "sample_id") %>%
  mutate(log_counts = log10(counts + 1))

p3 <- ggplot(count_long, aes(x = sample_id, y = log_counts)) +
  geom_boxplot(aes(fill = condition), alpha = 0.7) +
  scale_fill_manual(values = color_schemes$condition) +
  labs(title = "Count Distribution per Sample", 
       x = "Sample", y = "Log10(Counts + 1)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dirs$qc, "count_distributions.pdf"), p3, width = 12, height = 6)
ggsave(file.path(output_dirs$qc, "count_distributions.png"), p3, width = 12, height = 6, dpi = 300)

# =============================================================================
# DESEQ2 DATASET CREATION AND NORMALIZATION
# =============================================================================

log_messages("Creating DESeq2 dataset and performing normalization...")

# Create DESeqDataSet
# Start with simple design focused on main comparison
design_formula <- ~ condition

# Check if we have balanced design for more complex models
condition_counts <- table(metadata_clean$condition)
region_counts <- table(metadata_clean$region)
sex_counts <- table(metadata_clean$sex)

# Only add region if we have enough samples in each combination
if (all(condition_counts >= 10) && all(region_counts >= 10)) {
  # Check for confounding
  condition_region_table <- table(metadata_clean$condition, metadata_clean$region)
  if (all(condition_region_table >= 3)) {
    design_formula <- ~ condition + region
  }
}

log_messages(paste("Using design formula:", deparse(design_formula)))

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_filtered,
  colData = metadata_clean,
  design = design_formula
)

# Add gene information
mcols(dds)$baseMean <- rowMeans(counts(dds))
mcols(dds)$baseVar <- apply(counts(dds), 1, var)

# Perform DESeq2 analysis
log_messages("Running DESeq2 analysis...")
dds <- DESeq(dds)

# Extract normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
log_normalized_counts <- log2(normalized_counts + 1)

# Variance stabilizing transformation for downstream analysis
vsd <- vst(dds, blind = FALSE)
rlog_counts <- rlog(dds, blind = FALSE)

log_messages("DESeq2 analysis completed successfully")

# =============================================================================
# SAMPLE QUALITY ASSESSMENT
# =============================================================================

log_messages("Performing sample quality assessment...")

# 1. PCA analysis
pca_data <- plotPCA(vsd, intgroup = design_factors, returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p4 <- ggplot(pca_data, aes(PC1, PC2)) +
  geom_point(aes(color = condition, shape = sex), size = 3) +
  scale_color_manual(values = color_schemes$condition) +
  labs(
    title = "PCA - Sample Clustering",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(output_dirs$qc, "pca_samples.pdf"), p4, width = 10, height = 8)
ggsave(file.path(output_dirs$qc, "pca_samples.png"), p4, width = 10, height = 8, dpi = 300)

# 2. Sample correlation heatmap
sample_cor <- cor(log_normalized_counts, method = "pearson")

# Create annotation for heatmap
annotation_df <- metadata_clean[, design_factors, drop = FALSE]
rownames(annotation_df) <- metadata_clean$sample_id

# Define annotation colors
ann_colors <- list()
for (factor in design_factors) {
  if (factor %in% names(color_schemes)) {
    ann_colors[[factor]] <- color_schemes[[factor]]
  } else {
    levels_count <- length(unique(annotation_df[[factor]]))
    ann_colors[[factor]] <- RColorBrewer::brewer.pal(min(levels_count, 8), "Set2")[1:levels_count]
    names(ann_colors[[factor]]) <- unique(annotation_df[[factor]])
  }
}

# Create heatmap
pdf(file.path(output_dirs$qc, "sample_correlation_heatmap.pdf"), width = 12, height = 10)
pheatmap(
  sample_cor,
  annotation_col = annotation_df,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = "Sample-Sample Correlation"
)
dev.off()

png(file.path(output_dirs$qc, "sample_correlation_heatmap.png"), 
    width = 12, height = 10, units = "in", res = 300)
pheatmap(
  sample_cor,
  annotation_col = annotation_df,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = "Sample-Sample Correlation"
)
dev.off()

# 3. Outlier detection using PCA
pca_scores <- prcomp(t(log_normalized_counts), scale. = TRUE)
pc_scores <- pca_scores$x[, 1:3]

# Calculate Mahalanobis distance
mahal_dist <- mahalanobis(pc_scores, colMeans(pc_scores), cov(pc_scores))
outlier_threshold <- qchisq(0.95, df = 3)  # 95% confidence
outliers <- names(mahal_dist)[mahal_dist > outlier_threshold]

if (length(outliers) > 0) {
  log_messages(paste("Potential outliers detected:", paste(outliers, collapse = ", ")))
  
  # Create outlier plot
  outlier_data <- data.frame(
    sample_id = names(mahal_dist),
    mahal_distance = mahal_dist,
    is_outlier = mahal_dist > outlier_threshold
  ) %>%
    left_join(metadata_clean, by = "sample_id")
  
  p5 <- ggplot(outlier_data, aes(x = sample_id, y = mahal_distance)) +
    geom_point(aes(color = is_outlier, shape = condition), size = 3) +
    geom_hline(yintercept = outlier_threshold, linetype = "dashed", color = "red") +
    scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
    labs(title = "Outlier Detection (Mahalanobis Distance)",
         x = "Sample", y = "Mahalanobis Distance") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(output_dirs$qc, "outlier_detection.pdf"), p5, width = 12, height = 6)
  ggsave(file.path(output_dirs$qc, "outlier_detection.png"), p5, width = 12, height = 6, dpi = 300)
} else {
  log_messages("No outliers detected")
}

# =============================================================================
# BATCH EFFECT ASSESSMENT
# =============================================================================

if ("batch" %in% names(metadata_clean)) {
  log_messages("Assessing batch effects...")
  
  # RLE plot for batch effects
  rle_data <- log_normalized_counts - rowMedians(log_normalized_counts)
  rle_long <- rle_data %>%
    as.data.frame() %>%
    mutate(gene_id = rownames(.)) %>%
    pivot_longer(-gene_id, names_to = "sample_id", values_to = "rle") %>%
    left_join(metadata_clean, by = "sample_id")
  
  p6 <- ggplot(rle_long, aes(x = sample_id, y = rle)) +
    geom_boxplot(aes(fill = batch), alpha = 0.7) +
    scale_fill_manual(values = color_schemes$batch_colors) +
    labs(title = "Relative Log Expression (RLE) by Batch",
         x = "Sample", y = "RLE") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(output_dirs$qc, "rle_batch_effects.pdf"), p6, width = 12, height = 6)
  ggsave(file.path(output_dirs$qc, "rle_batch_effects.png"), p6, width = 12, height = 6, dpi = 300)
}

# =============================================================================
# DIFFERENTIAL EXPRESSION ANALYSIS
# =============================================================================

log_messages("Performing differential expression analysis...")

# Define contrasts based on your experimental design
# Only test contrasts that are in the design formula
contrasts_to_test <- list()

# Basic OUD vs Control (always included in design)
if ("OUD" %in% levels(metadata_clean$condition) && "Control" %in% levels(metadata_clean$condition)) {
  contrasts_to_test[["OUD_vs_Control"]] <- c("condition", "OUD", "Control")
}

# Region-specific contrasts (only if region is in the design)
if (grepl("region", deparse(design_formula)) && "region" %in% names(metadata_clean)) {
  region_levels <- levels(metadata_clean$region)
  if (length(region_levels) >= 2) {
    contrasts_to_test[[paste(region_levels[1], "vs", region_levels[2], sep = "_")]] <- 
      c("region", region_levels[1], region_levels[2])
  }
}

# Store all results
de_results <- list()
de_summary <- data.frame()

# Run differential expression for each contrast
for (contrast_name in names(contrasts_to_test)) {
  log_messages(paste("Testing contrast:", contrast_name))
  
  contrast_vec <- contrasts_to_test[[contrast_name]]
  
  tryCatch({
    # Extract results
    res <- results(dds, contrast = contrast_vec, alpha = 0.05)
    
    # Shrink log2 fold changes
    res_shrunk <- lfcShrink(dds, contrast = contrast_vec, res = res, type = "ashr")
    
    # Convert to data frame and add gene info
    res_df <- res_shrunk %>%
      as.data.frame() %>%
      mutate(
        gene_id = rownames(.),
        contrast = contrast_name,
        significant_fdr05 = padj < 0.05 & !is.na(padj),
        significant_fdr01 = padj < 0.01 & !is.na(padj),
        significant_nominal = pvalue < 0.05 & !is.na(pvalue),
        direction = case_when(
          log2FoldChange > 0 & significant_fdr05 ~ "Up",
          log2FoldChange < 0 & significant_fdr05 ~ "Down", 
          TRUE ~ "NS"
        )
      ) %>%
      arrange(padj)
    
    de_results[[contrast_name]] <- res_df
    
    # Summary statistics
    n_up <- sum(res_df$direction == "Up", na.rm = TRUE)
    n_down <- sum(res_df$direction == "Down", na.rm = TRUE)
    n_total_sig <- n_up + n_down
    
    de_summary <- rbind(de_summary, data.frame(
      contrast = contrast_name,
      n_genes_tested = nrow(res_df),
      n_significant_fdr05 = n_total_sig,
      n_upregulated = n_up,
      n_downregulated = n_down,
      median_pvalue = median(res_df$pvalue, na.rm = TRUE),
      median_padj = median(res_df$padj, na.rm = TRUE)
    ))
    
    log_messages(paste("Contrast", contrast_name, "completed:", n_total_sig, "significant genes"))
    
  }, error = function(e) {
    log_messages(paste("Error in contrast", contrast_name, ":", e$message))
  })
}

log_messages(paste("Differential expression analysis completed for", length(de_results), "contrasts"))

# =============================================================================
# VISUALIZATION OF DE RESULTS
# =============================================================================

log_messages("Creating differential expression visualizations...")

# Create plots for each contrast
for (contrast_name in names(de_results)) {
  res_df <- de_results[[contrast_name]]
  
  # 1. Volcano plot
  volcano_data <- res_df %>%
    mutate(
      log_padj = -log10(padj),
      label = ifelse(significant_fdr05 & abs(log2FoldChange) > 1, gene_id, "")
    )
  
  p_volcano <- ggplot(volcano_data, aes(x = log2FoldChange, y = log_padj)) +
    geom_point(aes(color = direction), alpha = 0.6, size = 1) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
    labs(
      title = paste("Volcano Plot:", contrast_name),
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Add labels for top genes
  if (sum(volcano_data$label != "") > 0) {
    p_volcano <- p_volcano + 
      geom_text_repel(aes(label = label), size = 3, max.overlaps = 20)
  }
  
  ggsave(file.path(output_dirs$plots, paste0("volcano_", contrast_name, ".pdf")), 
         p_volcano, width = 10, height = 8)
  ggsave(file.path(output_dirs$plots, paste0("volcano_", contrast_name, ".png")), 
         p_volcano, width = 10, height = 8, dpi = 300)
  
  # 2. MA plot
  p_ma <- ggplot(res_df, aes(x = log10(baseMean), y = log2FoldChange)) +
    geom_point(aes(color = direction), alpha = 0.6, size = 1) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    labs(
      title = paste("MA Plot:", contrast_name),
      x = "Log10 Mean Expression",
      y = "Log2 Fold Change"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave(file.path(output_dirs$plots, paste0("ma_plot_", contrast_name, ".pdf")), 
         p_ma, width = 10, height = 8)
  ggsave(file.path(output_dirs$plots, paste0("ma_plot_", contrast_name, ".png")), 
         p_ma, width = 10, height = 8, dpi = 300)
}

# 3. Summary barplot
p_summary <- ggplot(de_summary, aes(x = contrast)) +
  geom_col(aes(y = n_upregulated), fill = "red", alpha = 0.7) +
  geom_col(aes(y = -n_downregulated), fill = "blue", alpha = 0.7) +
  geom_hline(yintercept = 0, color = "black") +
  labs(
    title = "Differential Expression Summary",
    x = "Contrast",
    y = "Number of Genes (Up: positive, Down: negative)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dirs$plots, "de_summary.pdf"), p_summary, width = 10, height = 6)
ggsave(file.path(output_dirs$plots, "de_summary.png"), p_summary, width = 10, height = 6, dpi = 300)

# =============================================================================
# HEATMAP OF TOP GENES
# =============================================================================

log_messages("Creating heatmaps of top differentially expressed genes...")

# Get top DE genes across all contrasts
all_sig_genes <- unique(unlist(lapply(de_results, function(x) {
  x %>% filter(significant_fdr05, abs(log2FoldChange) > 1) %>% pull(gene_id)
})))

if (length(all_sig_genes) > 0) {
  # Limit to top genes to avoid overcrowding
  top_genes <- all_sig_genes[1:min(100, length(all_sig_genes))]
  
  # Extract expression data for top genes
  heatmap_data <- log_normalized_counts[top_genes, , drop = FALSE]
  
  # Z-score normalization across samples
  heatmap_data_scaled <- t(scale(t(heatmap_data)))
  
  # Create heatmap
  pdf(file.path(output_dirs$plots, "top_genes_heatmap.pdf"), width = 12, height = 10)
  pheatmap(
    heatmap_data_scaled,
    annotation_col = annotation_df,
    annotation_colors = ann_colors,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    show_rownames = nrow(heatmap_data_scaled) <= 50,
    main = "Top Differentially Expressed Genes"
  )
  dev.off()
  
  png(file.path(output_dirs$plots, "top_genes_heatmap.png"), 
      width = 12, height = 10, units = "in", res = 300)
  pheatmap(
    heatmap_data_scaled,
    annotation_col = annotation_df,
    annotation_colors = ann_colors,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    show_rownames = nrow(heatmap_data_scaled) <= 50,
    main = "Top Differentially Expressed Genes"
  )
  dev.off()
}

# =============================================================================
# SAVE RESULTS
# =============================================================================

log_messages("Saving analysis results...")

# Save DE results
for (contrast_name in names(de_results)) {
  write_csv(de_results[[contrast_name]], 
           file.path(output_dirs$tables, paste0("de_results_", contrast_name, ".csv")))
}

# Save summary
write_csv(de_summary, file.path(output_dirs$tables, "de_summary.csv"))

# Save normalized counts
write_csv(as.data.frame(normalized_counts) %>% rownames_to_column("gene_id"), 
         file.path(output_dirs$normalized, "normalized_counts.csv"))

# Save VST counts
write_csv(as.data.frame(assay(vsd)) %>% rownames_to_column("gene_id"), 
         file.path(output_dirs$normalized, "vst_counts.csv"))

# Save metadata
write_csv(metadata_clean, file.path(output_dirs$tables, "sample_metadata_processed.csv"))

# Save QC stats
write_csv(count_stats, file.path(output_dirs$tables, "sample_qc_stats.csv"))
write_csv(gene_stats_filtered, file.path(output_dirs$tables, "gene_qc_stats.csv"))

# Save DESeq2 object
saveRDS(dds, file.path(output_dirs$main, "deseq2_dataset.rds"))

# =============================================================================
# GENERATE HTML REPORT
# =============================================================================

log_messages("Generating HTML summary report...")

create_html_report <- function() {
  
  # Calculate summary statistics
  total_samples <- ncol(count_matrix)
  total_genes_initial <- nrow(count_matrix)
  total_genes_filtered <- nrow(count_matrix_filtered)
  total_contrasts <- length(de_results)
  
  html_content <- paste0(
    "<html><head><title>Bulk RNA-seq Analysis Report</title>",
    "<style>body{font-family: Arial, sans-serif; margin: 40px;}",
    "table{border-collapse: collapse; width: 100%;}",
    "th, td{border: 1px solid #ddd; padding: 8px; text-align: left;}",
    "th{background-color: #f2f2f2;}</style></head><body>",
    
    "<h1>Bulk RNA-seq Analysis Report - DESeq2</h1>",
    "<h2>Analysis Summary</h2>",
    "<p>Generated on: ", Sys.time(), "</p>",
    "<p>Analysis method: DESeq2</p>",
    "<p>Total samples: ", total_samples, "</p>",
    "<p>Total genes (initial): ", total_genes_initial, "</p>",
    "<p>Total genes (after filtering): ", total_genes_filtered, "</p>",
    "<p>Total contrasts tested: ", total_contrasts, "</p>",
    
    "<h2>Differential Expression Summary</h2>",
    "<table>",
    "<tr><th>Contrast</th><th>Genes Tested</th><th>Significant (FDR<0.05)</th>",
    "<th>Upregulated</th><th>Downregulated</th></tr>"
  )
  
  for (i in 1:nrow(de_summary)) {
    html_content <- paste0(html_content,
      "<tr><td>", de_summary$contrast[i], "</td>",
      "<td>", de_summary$n_genes_tested[i], "</td>", 
      "<td>", de_summary$n_significant_fdr05[i], "</td>",
      "<td>", de_summary$n_upregulated[i], "</td>",
      "<td>", de_summary$n_downregulated[i], "</td></tr>"
    )
  }
  
  html_content <- paste0(html_content, "</table>")
  
  # Add QC summary
  html_content <- paste0(html_content,
    "<h2>Quality Control Summary</h2>",
    "<p>Mean library size: ", round(mean(count_stats$total_counts)/1e6, 2), " million reads</p>",
    "<p>Mean detected genes: ", round(mean(count_stats$detected_genes)), " genes</p>"
  )
  
  if (length(outliers) > 0) {
    html_content <- paste0(html_content,
      "<p>Potential outliers detected: ", paste(outliers, collapse = ", "), "</p>"
    )
  } else {
    html_content <- paste0(html_content, "<p>No outliers detected</p>")
  }
  
  html_content <- paste0(html_content, "</body></html>")
  
  # Write HTML file
  writeLines(html_content, file.path(output_dirs$main, "bulk_rnaseq_report.html"))
  
  log_messages("Created HTML summary report")
}

create_html_report()

# =============================================================================
# FINAL SUMMARY
# =============================================================================

log_messages("Bulk RNA-seq Analysis Complete!")
log_messages("===============================")
log_messages(paste("Total samples analyzed:", ncol(count_matrix)))
log_messages(paste("Total genes (filtered):", nrow(count_matrix_filtered)))
log_messages(paste("Total contrasts tested:", length(de_results)))
log_messages(paste("Total significant genes (any contrast):", length(all_sig_genes)))
log_messages(paste("Output directory:", output_dirs$main))

# Print summary of significant genes per contrast
log_messages("Significant genes per contrast:")
for (i in 1:nrow(de_summary)) {
  log_messages(paste(de_summary$contrast[i], ":", de_summary$n_significant_fdr05[i], "genes"))
}

cat("\n🎉 Bulk RNA-seq Analysis completed successfully!\n")
cat("📁 Results saved in:", output_dirs$main, "\n")
cat("📊 Check the HTML report for a summary of findings\n")
cat("📈 Plots available in:", output_dirs$plots, "\n")
cat("📋 Tables available in:", output_dirs$tables, "\n")
cat("🔬 Ready for pathway analysis and comparison with snRNA-seq!\n")

# End of script
