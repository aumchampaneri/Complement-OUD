# Comprehensive Bulk RNA-seq QC and Preprocessing
# ----------------------------------------------

# Load required libraries
library(edgeR)
library(limma)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
# Optional packages
library(PCAtools)
library(sva)

# Create QC directory
qc_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE289002/QC"
dir.create(qc_dir, showWarnings = FALSE)

# File paths
count_file <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE289002/GSE289002_mouse_raw_counts.csv"
metadata_file <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE289002/mouse_metadata.csv"

# ----------------------
# 1. DATA LOADING
# ----------------------
counts <- read.csv(count_file, row.names = 1)
metadata <- read.csv(metadata_file)

# Check if samples in counts match the 'title' column in metadata
cat("Count matrix column names (first 10):", head(colnames(counts), 10), "\n")
cat("Metadata title values (first 10):", head(metadata$title, 10), "\n")

# Find matching samples between counts and metadata
common_samples <- intersect(colnames(counts), metadata$title)
cat("Found", length(common_samples), "matching samples between counts and metadata.\n")

if (length(common_samples) == 0) {
  stop("No matching samples between counts and metadata. Check column names and title values.")
} else {
  # Subset metadata to only include samples present in count data
  metadata <- metadata[metadata$title %in% colnames(counts), ]

  # Subset counts to only include samples present in metadata
  counts <- counts[, metadata$title]

  # Verify match
  all_match <- all(colnames(counts) == metadata$title)
  cat("All sample names match between counts and metadata:", all_match, "\n")
}

# Add treatment and sex as factors
metadata$treatment <- factor(metadata$treatment,
                           levels = c("Sal", "Mor + 24h", "Mor + 2W", "Chronic mor"))
metadata$sex <- factor(metadata$sex)
metadata$region <- factor(metadata$region)

# Basic dataset summary
cat("Dataset dimensions:", dim(counts)[1], "genes,", dim(counts)[2], "samples\n")
cat("Treatment groups:", paste(levels(metadata$treatment), collapse=", "), "\n")
cat("Brain regions:", paste(levels(metadata$region), collapse=", "), "\n")
cat("Sex distribution:", table(metadata$sex), "\n")

# ----------------------
# 2. GENE-LEVEL QC
# ----------------------

# Calculate gene stats
gene_stats <- data.frame(
  gene_id = rownames(counts),
  mean_count = rowMeans(counts),
  median_count = apply(counts, 1, median),
  min_count = apply(counts, 1, min),
  max_count = apply(counts, 1, max),
  zero_samples = rowSums(counts == 0),
  pct_zeros = rowSums(counts == 0) / ncol(counts) * 100
)

# Save gene stats
write.csv(gene_stats, file.path(qc_dir, "gene_stats.csv"), row.names = FALSE)

# Plot gene expression distribution
png(file.path(qc_dir, "gene_expression_distribution.png"), width = 1000, height = 800)
hist(log10(gene_stats$mean_count + 1), breaks = 100,
     xlab = "log10(mean count + 1)", main = "Gene Expression Distribution")
abline(v = log10(10), col = "red", lwd = 2, lty = 2)
dev.off()

# ----------------------
# 3. SAMPLE-LEVEL QC
# ----------------------

# Calculate library sizes
lib_sizes <- colSums(counts)
metadata$lib_size <- lib_sizes

# Create DGEList object
dge <- DGEList(counts = counts)

# Calculate CPM (useful for QC)
logcpm <- cpm(dge, log = TRUE)

# Plot library sizes
png(file.path(qc_dir, "library_sizes.png"), width = 1000, height = 800)
barplot(lib_sizes/1e6, names.arg = metadata$title, las = 2, cex.names = 0.7,
        main = "Library Sizes", ylab = "Millions of Reads")
abline(h = mean(lib_sizes/1e6), col = "red", lwd = 2, lty = 2)
dev.off()

# Boxplot of log CPM values
png(file.path(qc_dir, "log_cpm_boxplot.png"), width = 1200, height = 800)
boxplot(logcpm, las = 2, col = as.numeric(metadata$treatment),
        main = "Log2 CPM Distribution", cex.axis = 0.7)
legend("topright", legend = levels(metadata$treatment),
       fill = 1:length(levels(metadata$treatment)))
dev.off()

# Density plot of log CPM values
png(file.path(qc_dir, "log_cpm_density.png"), width = 1000, height = 800)
plotDensities(logcpm, col = as.numeric(metadata$treatment),
              legend = "topright", main = "Log2 CPM Density")
legend("topright", legend = levels(metadata$treatment),
       fill = 1:length(levels(metadata$treatment)))
dev.off()

# Sample correlation heatmap
sample_cor <- cor(logcpm)
annotation_col <- data.frame(
  Treatment = metadata$treatment,
  Sex = metadata$sex,
  Region = metadata$region,
  row.names = colnames(logcpm)
)

# Define colors for annotation
ann_colors <- list(
  Treatment = c("Sal" = "lightblue", "Mor + 24h" = "orange",
                "Mor + 2W" = "red", "Chronic mor" = "darkred"),
  Sex = c("male" = "blue", "female" = "pink"),
  Region = c("PFC" = "purple", "NAc" = "green")
)

png(file.path(qc_dir, "sample_correlation_heatmap.png"), width = 1200, height = 1000)
pheatmap(sample_cor, annotation_col = annotation_col, annotation_colors = ann_colors,
         show_colnames = FALSE, main = "Sample Correlation")
dev.off()

# ----------------------
# 4. EXPLORATORY DATA ANALYSIS
# ----------------------

# MDS plot
png(file.path(qc_dir, "mds_plot.png"), width = 1000, height = 800)
par(mfrow = c(1, 2))
# MDS by treatment
plotMDS(dge, col = as.numeric(metadata$treatment), main = "MDS by Treatment")
legend("topright", legend = levels(metadata$treatment),
       col = 1:length(levels(metadata$treatment)), pch = 16)
# MDS by region
plotMDS(dge, col = as.numeric(metadata$region), main = "MDS by Region")
legend("topright", legend = levels(metadata$region),
       col = 1:length(levels(metadata$region)), pch = 16)
dev.off()

# PCA analysis
pca_res <- prcomp(t(logcpm))
variance_explained <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100

# PCA plot by treatment and region
pca_data <- data.frame(
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2],
  PC3 = pca_res$x[,3],
  Treatment = metadata$treatment,
  Sex = metadata$sex,
  Region = metadata$region
)

png(file.path(qc_dir, "pca_treatment_region.png"), width = 1200, height = 600)
par(mfrow = c(1, 2))
# PCA by treatment
plot(pca_data$PC1, pca_data$PC2, col = as.numeric(pca_data$Treatment),
     pch = as.numeric(pca_data$Sex), cex = 1.5,
     xlab = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
     ylab = paste0("PC2 (", round(variance_explained[2], 2), "%)"),
     main = "PCA - Treatment & Sex")
legend("topright", legend = c(levels(pca_data$Treatment),
                             paste0("Sex: ", levels(pca_data$Sex))),
       col = c(1:length(levels(pca_data$Treatment)), rep(1, length(levels(pca_data$Sex)))),
       pch = c(rep(1, length(levels(pca_data$Treatment))), 1:length(levels(pca_data$Sex))))

# PCA by region
plot(pca_data$PC1, pca_data$PC2, col = as.numeric(pca_data$Region),
     pch = 16, cex = 1.5,
     xlab = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
     ylab = paste0("PC2 (", round(variance_explained[2], 2), "%)"),
     main = "PCA - Brain Region")
legend("topright", legend = levels(pca_data$Region),
       col = 1:length(levels(pca_data$Region)), pch = 16)
dev.off()

# ----------------------
# 5. FILTERING AND NORMALIZATION
# ----------------------

# Filter genes by expression level - keep genes expressed in at least 3 samples with > 10 counts
keep <- filterByExpr(dge, min.count = 10, min.total.count = 15)
dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]
cat("Genes before filtering:", nrow(dge), "\n")
cat("Genes after filtering:", nrow(dge_filtered), "\n")

# Normalize using TMM method
dge_filtered <- calcNormFactors(dge_filtered, method = "TMM")
logcpm_filtered_norm <- cpm(dge_filtered, log = TRUE)

# Before-after normalization comparison
png(file.path(qc_dir, "before_after_normalization.png"), width = 1200, height = 600)
par(mfrow = c(1, 2))
boxplot(logcpm, main = "Before normalization", las = 2, cex.axis = 0.7)
boxplot(logcpm_filtered_norm, main = "After filtering & normalization", las = 2, cex.axis = 0.7)
dev.off()

# ----------------------
# 6. BATCH EFFECT ASSESSMENT
# ----------------------

# Check if batch is a significant factor using PCA
if("batch" %in% colnames(metadata)) {
  metadata$batch <- factor(metadata$batch)

  # PCA plot by batch
  png(file.path(qc_dir, "pca_batch_effect.png"), width = 800, height = 800)
  plot(pca_data$PC1, pca_data$PC2, col = as.numeric(factor(metadata$batch)),
       pch = 16, cex = 1.5,
       xlab = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
       ylab = paste0("PC2 (", round(variance_explained[2], 2), "%)"),
       main = "PCA - Batch Effect")
  legend("topright", legend = levels(factor(metadata$batch)),
         col = 1:length(levels(factor(metadata$batch))), pch = 16)
  dev.off()
}

# ----------------------
# 7. FINAL QC VISUALIZATION
# ----------------------

# MDS plot after filtering and normalization
png(file.path(qc_dir, "mds_after_normalization.png"), width = 1000, height = 800)
plotMDS(dge_filtered, col = as.numeric(metadata$treatment) + 4 * as.numeric(metadata$region),
        main = "MDS after filtering and normalization")
legend("topright",
       legend = c(paste("Treatment:", levels(metadata$treatment)),
                  paste("Region:", levels(metadata$region))),
       col = c(1:4, 5:6), pch = 19)
dev.off()

# PCA after normalization
pca_norm <- prcomp(t(logcpm_filtered_norm))
variance_explained_norm <- (pca_norm$sdev^2) / sum(pca_norm$sdev^2) * 100

# PCA plot
pca_norm_data <- data.frame(
  PC1 = pca_norm$x[,1],
  PC2 = pca_norm$x[,2],
  Treatment = metadata$treatment,
  Sex = metadata$sex,
  Region = metadata$region
)

png(file.path(qc_dir, "pca_after_normalization.png"), width = 1000, height = 800)
plot(pca_norm_data$PC1, pca_norm_data$PC2,
     col = as.numeric(pca_norm_data$Treatment) + 4 * as.numeric(pca_norm_data$Region),
     pch = as.numeric(pca_norm_data$Sex) + 2,
     cex = 1.5,
     xlab = paste0("PC1 (", round(variance_explained_norm[1], 2), "%)"),
     ylab = paste0("PC2 (", round(variance_explained_norm[2], 2), "%)"),
     main = "PCA - After Normalization")
legend("topright",
       legend = c(levels(pca_norm_data$Treatment),
                  levels(pca_norm_data$Region),
                  levels(pca_norm_data$Sex)),
       col = c(1:4, 5:6, rep(1,2)),
       pch = c(rep(1,4), rep(1,2), 3:4))
dev.off()

# Save the processed data for downstream analysis
saveRDS(dge_filtered, file.path(qc_dir, "dge_filtered_normalized.rds"))
saveRDS(metadata, file.path(qc_dir, "metadata.rds"))

# ----------------------
# 8. SUMMARY REPORT
# ----------------------
sink(file.path(qc_dir, "preprocessing_summary.txt"))
cat("=== BULK RNA-SEQ PREPROCESSING SUMMARY ===\n\n")
cat("Total samples:", ncol(counts), "\n")
cat("Original genes:", nrow(counts), "\n")
cat("Genes after filtering:", nrow(dge_filtered), "\n")
cat("Normalization method: TMM\n\n")

cat("=== SAMPLE GROUPS ===\n")
print(table(metadata$treatment, metadata$region))
cat("\n")
print(table(metadata$treatment, metadata$sex))
cat("\n")

cat("=== LIBRARY SIZE STATS ===\n")
print(summary(lib_sizes))
cat("\n")

cat("=== TOP PC VARIANCE EXPLAINED ===\n")
print(data.frame(PC = 1:5,
                 Variance = variance_explained_norm[1:5],
                 Cumulative = cumsum(variance_explained_norm)[1:5]))
sink()

cat("Preprocessing complete. Results saved in:", qc_dir, "\n")