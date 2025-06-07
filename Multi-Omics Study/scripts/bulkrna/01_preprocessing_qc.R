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
library(biomaRt)

# Set paths
setwd("/Users/aumchampaneri/Complement-OUD/Multi-Omics Study")
data_dir <- "data/raw/bulkrna"
output_dir <- "data/processed/bulkrna/preprocessing"
plots_dir <- "results/bulkrna/preprocessing"

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
cat("Converting Ensembl IDs to Gene Symbols...\n")

# Check if row names look like Ensembl IDs
ensembl_pattern <- "^ENS[A-Z]*[0-9]+$|^ENS[A-Z]*[0-9]+\\.[0-9]+$"
if (any(grepl(ensembl_pattern, rownames(counts)[1:10]))) {
  cat("Detected Ensembl IDs, converting to gene symbols...\n")
  
  # Remove version numbers from Ensembl IDs if present
  ensembl_ids_clean <- sub("\\.[0-9]+$", "", rownames(counts))
  
  tryCatch({
    # Method 1: Using org.Hs.eg.db
    cat("Trying org.Hs.eg.db conversion...\n")
    gene_symbols_org <- mapIds(org.Hs.eg.db,
                              keys = ensembl_ids_clean,
                              column = "SYMBOL",
                              keytype = "ENSEMBL",
                              multiVals = "first")
    
    # Method 2: Using biomaRt with multiple attributes for genes not found
    missing_genes <- is.na(gene_symbols_org)
    if (sum(missing_genes) > 0) {
      cat("Using biomaRt for", sum(missing_genes), "remaining genes...\n")
      
      # Connect to biomart
      ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
      
      # Get multiple gene name attributes
      biomart_results <- getBM(
        attributes = c('ensembl_gene_id', 'external_gene_name', 'hgnc_symbol', 'description'),
        filters = 'ensembl_gene_id',
        values = ensembl_ids_clean[missing_genes],
        mart = ensembl
      )
      
      # Update gene symbols with biomart results
      for (i in which(missing_genes)) {
        ensembl_id <- ensembl_ids_clean[i]
        biomart_row <- biomart_results[biomart_results$ensembl_gene_id == ensembl_id, ]
        
        if (nrow(biomart_row) > 0) {
          # Try external_gene_name first, then hgnc_symbol
          if (biomart_row$external_gene_name[1] != "" && !is.na(biomart_row$external_gene_name[1])) {
            gene_symbols_org[i] <- biomart_row$external_gene_name[1]
          } else if (biomart_row$hgnc_symbol[1] != "" && !is.na(biomart_row$hgnc_symbol[1])) {
            gene_symbols_org[i] <- biomart_row$hgnc_symbol[1]
          }
        }
      }
    }
    
    # Method 3: Try gene synonyms for remaining missing genes using ENTREZID
    still_missing <- is.na(gene_symbols_org) | gene_symbols_org == ""
    if (sum(still_missing) > 0) {
      cat("Trying alternative gene mapping for", sum(still_missing), "remaining genes...\n")
      
      tryCatch({
        # First try to get ENTREZ IDs for the remaining Ensembl IDs
        entrez_ids <- mapIds(org.Hs.eg.db,
                           keys = ensembl_ids_clean[still_missing],
                           column = "ENTREZID",
                           keytype = "ENSEMBL",
                           multiVals = "first")
        
        # Then get symbols from ENTREZ IDs
        for (i in which(still_missing)) {
          idx_in_missing <- which(which(still_missing) == i)
          if (!is.na(entrez_ids[idx_in_missing])) {
            symbol_from_entrez <- mapIds(org.Hs.eg.db,
                                       keys = entrez_ids[idx_in_missing],
                                       column = "SYMBOL",
                                       keytype = "ENTREZID",
                                       multiVals = "first")
            if (!is.na(symbol_from_entrez) && symbol_from_entrez != "") {
              gene_symbols_org[i] <- symbol_from_entrez
            }
          }
        }
      }, error = function(e) {
        cat("Alternative mapping method failed:", e$message, "\n")
      })
    }
    
    # Method 4: Try to get gene types/biotypes for remaining genes
    final_missing <- is.na(gene_symbols_org) | gene_symbols_org == ""
    if (sum(final_missing) > 0) {
      cat("Adding gene type information for", sum(final_missing), "remaining genes...\n")
      
      tryCatch({
        # Get gene biotypes from biomaRt for categorization
        biotype_results <- getBM(
          attributes = c('ensembl_gene_id', 'gene_biotype'),
          filters = 'ensembl_gene_id',
          values = ensembl_ids_clean[final_missing],
          mart = ensembl
        )
        
        # Create informative names for genes without symbols
        for (i in which(final_missing)) {
          ensembl_id <- ensembl_ids_clean[i]
          biotype_row <- biotype_results[biotype_results$ensembl_gene_id == ensembl_id, ]
          
          if (nrow(biotype_row) > 0 && !is.na(biotype_row$gene_biotype[1])) {
            # Create a more informative name
            gene_symbols_org[i] <- paste0(ensembl_id, "_", biotype_row$gene_biotype[1])
          }
        }
      }, error = function(e) {
        cat("Gene type annotation failed:", e$message, "\n")
      })
    }

    # Create final mapping
    gene_mapping <- data.frame(
      ensembl_id = ensembl_ids_clean,
      original_id = rownames(counts),
      gene_symbol = gene_symbols_org,
      stringsAsFactors = FALSE
    )
    
    # For genes without symbols, keep original ID
    gene_mapping$final_name <- ifelse(is.na(gene_mapping$gene_symbol) | gene_mapping$gene_symbol == "", 
                                     gene_mapping$original_id, 
                                     gene_mapping$gene_symbol)
    
    # Handle duplicated gene symbols
    gene_mapping$final_name <- make.unique(gene_mapping$final_name)
    
    # Update row names
    rownames(counts) <- gene_mapping$final_name
    
    # Save mapping
    write.table(gene_mapping, file.path(output_dir, "gene_id_mapping.txt"), 
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Report conversion statistics
    converted_count <- sum(!is.na(gene_symbols_org) & gene_symbols_org != "")
    cat("Successfully converted", converted_count, "out of", length(gene_symbols_org), "Ensembl IDs\n")
    cat("Conversion rate:", round(converted_count/length(gene_symbols_org)*100, 1), "%\n")
    
    # Report final gene name statistics
    ensembl_remaining <- sum(grepl("^ENS", gene_mapping$final_name))
    cat("Genes with symbol names:", nrow(gene_mapping) - ensembl_remaining, "\n")
    cat("Genes remaining as Ensembl IDs:", ensembl_remaining, "\n")
    
  }, error = function(e) {
    cat("Error in gene symbol conversion:", e$message, "\n")
    cat("Proceeding with original Ensembl IDs\n")
  })
  
} else {
  cat("Row names do not appear to be Ensembl IDs, keeping original names\n")
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

# Library size and distribution plots
library_sizes <- data.frame(
  sample = colnames(counts),
  total_counts = colSums(counts),
  condition = sample_info$condition,
  region = sample_info$region,
  sex = sample_info$sex
)

# Library size boxplot
p_libsize <- ggplot(library_sizes, aes(x = condition, y = total_counts, fill = condition)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  scale_y_log10() +
  labs(title = "Library Sizes by Condition",
       x = "Condition", y = "Total Counts (log10)") +
  theme_minimal()

ggsave(file.path(plots_dir, "library_sizes.png"), p_libsize, width = 8, height = 6)

# Density plot of normalized counts
vsd_mat <- assay(vsd)
vsd_df <- data.frame(
  values = as.vector(vsd_mat),
  sample = rep(colnames(vsd_mat), each = nrow(vsd_mat)),
  condition = rep(sample_info$condition, each = nrow(vsd_mat))
)

p_density <- ggplot(vsd_df, aes(x = values, color = condition)) +
  geom_density(alpha = 0.7) +
  labs(title = "Density of VST Normalized Counts",
       x = "VST Normalized Counts", y = "Density") +
  theme_minimal()

ggsave(file.path(plots_dir, "count_distribution.png"), p_density, width = 10, height = 6)

# Mean-variance relationship plot
mean_var_data <- data.frame(
  mean = rowMeans(counts(dds_norm)),
  variance = apply(counts(dds_norm), 1, var)
)

p_meanvar <- ggplot(mean_var_data, aes(x = mean, y = variance)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = "loess", color = "red") +
  labs(title = "Mean-Variance Relationship",
       x = "Mean Count (log10)", y = "Variance (log10)") +
  theme_minimal()

ggsave(file.path(plots_dir, "mean_variance.png"), p_meanvar, width = 8, height = 6)

# Cook's distance plot to identify outlier samples
cooks_distances <- apply(assays(dds_norm)[["cooks"]], 2, max, na.rm = TRUE)
cooks_df <- data.frame(
  sample = names(cooks_distances),
  cooks_distance = cooks_distances,
  condition = sample_info$condition
)

p_cooks <- ggplot(cooks_df, aes(x = reorder(sample, cooks_distance), y = cooks_distance, fill = condition)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = qf(0.99, df1 = ncol(sample_info), df2 = ncol(counts) - ncol(sample_info)), 
             linetype = "dashed", color = "red") +
  labs(title = "Cook's Distance by Sample",
       x = "Sample", y = "Cook's Distance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plots_dir, "cooks_distance.png"), p_cooks, width = 12, height = 6)

# Gene detection rates
detection_rates <- data.frame(
  sample = colnames(counts),
  genes_detected = colSums(counts > 0),
  condition = sample_info$condition
)

p_detection <- ggplot(detection_rates, aes(x = condition, y = genes_detected, fill = condition)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  labs(title = "Gene Detection Rates by Condition",
       x = "Condition", y = "Number of Genes Detected") +
  theme_minimal()

ggsave(file.path(plots_dir, "gene_detection.png"), p_detection, width = 8, height = 6)

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
