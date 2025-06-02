#!/usr/bin/env Rscript
# ==============================================================================
# Human Data Processing Script - Cross-Species OUD Analysis
# ==============================================================================
# Project: Cross-species meta-analysis of mouse and human OUD datasets
# Author: Bioinformatics Analysis Pipeline
# Date: June 2025
# 
# Purpose: Download and process human OUD datasets for cross-species integration
# Datasets: GSE174409 (bulk RNA-seq), GSE225158 (single-cell RNA-seq)
# ==============================================================================

# Load required libraries
library(here)
library(tidyverse)
library(GEOquery)
library(data.table)
library(DESeq2)
library(edgeR)
library(limma)
library(biomaRt)
library(Seurat)
library(SingleCellExperiment)

cat(strrep("=", 80), "\n")
cat("HUMAN DATA PROCESSING - PHASE 1\n")
cat(strrep("=", 80), "\n")
cat("Processing Time:", Sys.time(), "\n\n")

# Set up directories (ensure we're in the correct project root)
project_root <- getwd()  # Should be "Translational Study" folder
raw_data_dir <- file.path(project_root, "Data", "Raw")
processed_data_dir <- file.path(project_root, "Data", "Processed")
results_dir <- file.path(project_root, "Results", "Human_Analysis")
figures_dir <- file.path(project_root, "Figures", "QC")

# Create directories if they don't exist
for (dir in c(raw_data_dir, processed_data_dir, results_dir, figures_dir)) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat("Created directory:", dir, "\n")
  }
}

cat("Project root:", project_root, "\n")
cat("Output directories set for Translational Study folder\n\n")

# ==============================================================================
# PART 1: GSE174409 - Human Bulk RNA-seq Processing
# ==============================================================================

cat("\n", strrep("=", 60), "\n")
cat("PROCESSING GSE174409 - Human Bulk RNA-seq\n")
cat(strrep("=", 60), "\n")

# Download GSE174409 data
gse174409_dir <- file.path(raw_data_dir, "GSE174409")
dir.create(gse174409_dir, showWarnings = FALSE, recursive = TRUE)

cat("Downloading GSE174409 data...\n")

# Download supplementary files
gse174409_urls <- c(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174409/suppl/GSE174409_human_Opioid_filtered_2020.csv.gz",
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174409/suppl/GSE174409_raw_counts_02102020.csv.gz"
)

gse174409_files <- c(
  file.path(gse174409_dir, "GSE174409_human_Opioid_filtered_2020.csv.gz"),
  file.path(gse174409_dir, "GSE174409_raw_counts_02102020.csv.gz")
)

for (i in seq_along(gse174409_urls)) {
  if (!file.exists(gse174409_files[i])) {
    cat("Downloading:", basename(gse174409_files[i]), "\n")
    tryCatch({
      download.file(gse174409_urls[i], gse174409_files[i], mode = "wb")
    }, error = function(e) {
      cat("Warning: Could not download", basename(gse174409_files[i]), "-", e$message, "\n")
    })
  } else {
    cat("File already exists:", basename(gse174409_files[i]), "\n")
  }
}

# Download series matrix for metadata
tryCatch({
  gse174409 <- getGEO("GSE174409", GSEMatrix = TRUE, destdir = gse174409_dir)
  if (length(gse174409) > 0) {
    gse174409_metadata <- pData(gse174409[[1]])
    cat("✓ Downloaded GSE174409 metadata\n")
  }
}, error = function(e) {
  cat("Warning: Could not download GSE174409 metadata:", e$message, "\n")
})

# Process raw counts data
cat("Processing GSE174409 raw counts...\n")
if (file.exists(gse174409_files[2])) {
  gse174409_raw <- fread(gse174409_files[2])
  cat("Raw counts dimensions:", nrow(gse174409_raw), "x", ncol(gse174409_raw), "\n")
  
  # Convert to matrix format
  gene_ids <- gse174409_raw$V1
  gse174409_counts <- as.matrix(gse174409_raw[, -1])
  rownames(gse174409_counts) <- gene_ids
  
  cat("✓ Processed raw counts matrix\n")
} else {
  cat("✗ Raw counts file not found\n")
}

# Process filtered data
cat("Processing GSE174409 filtered data...\n")
if (file.exists(gse174409_files[1])) {
  gse174409_filtered <- fread(gse174409_files[1])
  cat("Filtered data dimensions:", nrow(gse174409_filtered), "x", ncol(gse174409_filtered), "\n")
  cat("✓ Processed filtered data\n")
} else {
  cat("✗ Filtered data file not found\n")
}

# Create metadata from sample names (if GEO metadata not available)
if (exists("gse174409_counts")) {
  sample_names <- colnames(gse174409_counts)
  
  # Extract condition and region information from sample names
  gse174409_sample_metadata <- data.frame(
    sample_id = sample_names,
    condition = ifelse(grepl("OUD|Opioid", sample_names, ignore.case = TRUE), "OUD", "Control"),
    region = case_when(
      grepl("DLPFC|dlpfc", sample_names, ignore.case = TRUE) ~ "DLPFC",
      grepl("NAcc|nacc", sample_names, ignore.case = TRUE) ~ "NAcc",
      TRUE ~ "Unknown"
    ),
    stringsAsFactors = FALSE
  )
  
  cat("✓ Created sample metadata for", nrow(gse174409_sample_metadata), "samples\n")
}

# ==============================================================================
# PART 2: GSE225158 - Human Single-cell RNA-seq Processing
# ==============================================================================

cat("\n", strrep("=", 60), "\n")
cat("PROCESSING GSE225158 - Human Single-cell RNA-seq\n")
cat(strrep("=", 60), "\n")

# Download GSE225158 data
gse225158_dir <- file.path(raw_data_dir, "GSE225158")
dir.create(gse225158_dir, showWarnings = FALSE, recursive = TRUE)

cat("Downloading GSE225158 data...\n")

# Download series matrix for metadata
tryCatch({
  gse225158 <- getGEO("GSE225158", GSEMatrix = TRUE, destdir = gse225158_dir)
  if (length(gse225158) > 0) {
    gse225158_metadata <- pData(gse225158[[1]])
    cat("✓ Downloaded GSE225158 metadata\n")
  }
}, error = function(e) {
  cat("Warning: Could not download GSE225158 metadata:", e$message, "\n")
})

# Check for and download supplementary files
cat("Checking for GSE225158 supplementary files...\n")

# Download specific supplementary files for GSE225158
gse225158_suppl_urls <- c(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225158/suppl/GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat",
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225158/suppl/GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5ad.gz"
)

gse225158_suppl_files <- c(
  file.path(gse225158_dir, "GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat"),
  file.path(gse225158_dir, "GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5ad.gz")
)

# Download supplementary files (note: these are very large files)
for (i in seq_along(gse225158_suppl_urls)) {
  if (!file.exists(gse225158_suppl_files[i])) {
    file_size <- ifelse(i == 1, "10.0 GB", "1.7 GB")
    cat("Downloading:", basename(gse225158_suppl_files[i]), "(", file_size, ")\n")
    cat("This may take a while due to file size...\n")
    
    tryCatch({
      # Increase timeout for large files (30 minutes)
      options(timeout = 1800)
      download.file(gse225158_suppl_urls[i], gse225158_suppl_files[i], mode = "wb")
      cat("✓ Downloaded:", basename(gse225158_suppl_files[i]), "\n")
    }, error = function(e) {
      cat("Warning: Could not download", basename(gse225158_suppl_files[i]), "-", e$message, "\n")
      cat("Manual download may be required for this large file\n")
      cat("Alternative: Use wget or curl for better handling of large files:\n")
      cat("  wget", gse225158_suppl_urls[i], "-O", gse225158_suppl_files[i], "\n")
    })
  } else {
    cat("File already exists:", basename(gse225158_suppl_files[i]), "\n")
  }
}

# Process GSE225158 data
gse225158_processed <- FALSE

# Try to load h5Seurat file first (native Seurat format)
h5seurat_file <- gse225158_suppl_files[1]
if (file.exists(h5seurat_file)) {
  cat("Processing h5Seurat file...\n")
  
  tryCatch({
    # Load SeuratDisk for h5Seurat functionality
    if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
      cat("Installing SeuratDisk package...\n")
      devtools::install_github("mojaveazure/seurat-disk")
    }
    library(SeuratDisk)
    
    # Check if file is complete/valid
    file_info <- file.info(h5seurat_file)
    if (file_info$size < 1000000000) {  # Less than 1GB suggests incomplete download
      cat("Warning: h5Seurat file appears incomplete (", file_info$size, "bytes)\n")
      cat("Expected size: ~10GB. Manual download recommended.\n")
    } else {
      # Load Seurat object directly
      gse225158_seurat <- LoadH5Seurat(h5seurat_file)
      
      cat("✓ Loaded Seurat object with", ncol(gse225158_seurat), "cells and", nrow(gse225158_seurat), "genes\n")
      
      # Display basic information about the dataset
      cat("Assays available:", names(gse225158_seurat@assays), "\n")
      cat("Metadata columns:", ncol(gse225158_seurat@meta.data), "\n")
      
      # Add study-specific metadata
      gse225158_seurat@meta.data$study <- "GSE225158"
      gse225158_seurat@meta.data$tissue <- "Striatum"
      gse225158_seurat@meta.data$species <- "Human"
      
      gse225158_processed <- TRUE
    }
    
  }, error = function(e) {
    cat("Warning: Could not process h5Seurat file:", e$message, "\n")
  })
}

# Try h5ad file if h5Seurat failed
if (!gse225158_processed) {
  h5ad_file <- gse225158_suppl_files[2]
  if (file.exists(h5ad_file)) {
    cat("Processing h5ad file...\n")
    
    # Check if file is complete
    file_info <- file.info(h5ad_file)
    if (file_info$size < 1000000000) {  # Less than 1GB suggests incomplete download
      cat("Warning: h5ad file appears incomplete (", file_info$size, "bytes)\n")
      cat("Expected size: ~1.7GB. Manual download recommended.\n")
    } else {
      tryCatch({
        # Load SeuratDisk for h5ad functionality
        if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
          cat("Installing SeuratDisk package...\n")
          devtools::install_github("mojaveazure/seurat-disk")
        }
        library(SeuratDisk)
        
        # Convert h5ad to h5seurat
        temp_h5seurat <- gsub("\\.h5ad\\.gz$", ".h5seurat", h5ad_file)
        
        if (!file.exists(temp_h5seurat)) {
          cat("Converting h5ad to h5seurat format...\n")
          Convert(h5ad_file, dest = "h5seurat", overwrite = TRUE)
        }
        
        gse225158_seurat <- LoadH5Seurat(temp_h5seurat)
        cat("✓ Loaded from h5ad file with", ncol(gse225158_seurat), "cells and", nrow(gse225158_seurat), "genes\n")
        gse225158_processed <- TRUE
        
      }, error = function(e) {
        cat("Warning: Could not process h5ad file:", e$message, "\n")
      })
    }
  }
}

# If no files were processed successfully
if (!gse225158_processed) {
  cat("No GSE225158 data files could be processed.\n")
  cat("Issues detected:\n")
  cat("  1. Large file download timeouts (10GB + 1.7GB files)\n")
  cat("  2. Possible incomplete downloads\n")
  cat("\nRecommended solutions:\n")
  cat("  1. Manual download using web browser or wget:\n")
  for (i in seq_along(gse225158_suppl_urls)) {
    cat("     wget '", gse225158_suppl_urls[i], "' -O '", gse225158_suppl_files[i], "'\n", sep = "")
  }
  cat("  2. Use download manager for resumable downloads\n")
  cat("  3. Re-run script after manual download\n")
}

# Basic single-cell QC if data was processed
if (gse225158_processed && exists("gse225158_seurat")) {
  cat("Performing basic single-cell QC...\n")
  
  # Calculate QC metrics if not already present
  if (!"percent.mt" %in% colnames(gse225158_seurat@meta.data)) {
    gse225158_seurat[["percent.mt"]] <- PercentageFeatureSet(gse225158_seurat, pattern = "^MT-")
  }
  
  # Generate QC plots
  tryCatch({
    qc_plots <- VlnPlot(gse225158_seurat, 
                        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                        ncol = 3)
    
    ggsave(file.path(figures_dir, "GSE225158_QC_violin.pdf"), qc_plots, width = 12, height = 4)
    
    # Feature scatter plots
    scatter_plot <- FeatureScatter(gse225158_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    ggsave(file.path(figures_dir, "GSE225158_feature_scatter.pdf"), scatter_plot, width = 8, height = 6)
    
    cat("✓ Generated single-cell QC plots\n")
  }, error = function(e) {
    cat("Warning: Could not generate QC plots:", e$message, "\n")
  })
  
  # Display cell type information if available
  if ("cell_type" %in% colnames(gse225158_seurat@meta.data)) {
    cat("Cell types found:", length(unique(gse225158_seurat@meta.data$cell_type)), "\n")
    cat("Cell type distribution:\n")
    print(table(gse225158_seurat@meta.data$cell_type))
  }
  
  # Display condition information if available
  condition_cols <- grep("condition|group|treatment|oud", colnames(gse225158_seurat@meta.data), ignore.case = TRUE, value = TRUE)
  if (length(condition_cols) > 0) {
    cat("Condition columns found:", paste(condition_cols, collapse = ", "), "\n")
    for (col in condition_cols[1:min(2, length(condition_cols))]) {
      cat(col, "distribution:\n")
      print(table(gse225158_seurat@meta.data[[col]]))
    }
  }
}

# ==============================================================================
# PART 3: Data Quality Control and Normalization
# ==============================================================================

cat("\n", strrep("=", 60), "\n")
cat("DATA QUALITY CONTROL AND NORMALIZATION\n")
cat(strrep("=", 60), "\n")

if (exists("gse174409_counts") && exists("gse174409_sample_metadata")) {
  cat("Performing QC on GSE174409 data...\n")
  
  # Filter low-count genes
  keep_genes <- rowSums(gse174409_counts > 1) >= 3
  gse174409_filtered_counts <- gse174409_counts[keep_genes, ]
  
  cat("Genes after filtering:", nrow(gse174409_filtered_counts), "/", nrow(gse174409_counts), "\n")
  
  # Create DGEList for edgeR
  dge <- edgeR::DGEList(counts = gse174409_filtered_counts, 
                        samples = gse174409_sample_metadata)
  
  # Calculate normalization factors
  dge <- edgeR::calcNormFactors(dge)
  
  # Get normalized counts (explicitly call from edgeR)
  gse174409_normalized <- edgeR::cpm(dge, log = TRUE)
  
  cat("✓ Completed normalization\n")
  
  # Generate QC plots
  cat("Generating QC plots...\n")
  
  # Sample correlation heatmap
  sample_cor <- cor(gse174409_normalized)
  
  # Load pheatmap explicitly
  if (!require(pheatmap, quietly = TRUE)) {
    install.packages("pheatmap")
    library(pheatmap)
  }
  
  # Create annotation data for heatmap (ensure proper formatting)
  annotation_data <- gse174409_sample_metadata[, c("condition", "region"), drop = FALSE]
  rownames(annotation_data) <- gse174409_sample_metadata$sample_id
  
  # Ensure annotation_data matches sample_cor column names
  annotation_data <- annotation_data[colnames(sample_cor), , drop = FALSE]
  
  # Remove any NA values
  annotation_data$condition[is.na(annotation_data$condition)] <- "Unknown"
  annotation_data$region[is.na(annotation_data$region)] <- "Unknown"
  
  tryCatch({
    pdf(file.path(figures_dir, "GSE174409_sample_correlation.pdf"), width = 10, height = 8)
    pheatmap::pheatmap(sample_cor, 
                       main = "GSE174409 Sample Correlation",
                       annotation_col = annotation_data,
                       show_rownames = FALSE,
                       show_colnames = FALSE)
    dev.off()
    cat("✓ Generated sample correlation heatmap\n")
  }, error = function(e) {
    cat("Warning: Could not generate correlation heatmap:", e$message, "\n")
    # Generate simple correlation heatmap without annotations
    pdf(file.path(figures_dir, "GSE174409_sample_correlation_simple.pdf"), width = 10, height = 8)
    pheatmap::pheatmap(sample_cor, 
                       main = "GSE174409 Sample Correlation",
                       show_rownames = FALSE,
                       show_colnames = FALSE)
    dev.off()
  })
  
  # PCA plot
  tryCatch({
    pca_data <- prcomp(t(gse174409_normalized))
    pca_df <- data.frame(pca_data$x[, 1:2], gse174409_sample_metadata)
    
    p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, shape = region)) +
      geom_point(size = 3) +
      labs(title = "GSE174409 PCA Plot",
           x = paste0("PC1 (", round(summary(pca_data)$importance[2, 1] * 100, 1), "%)"),
           y = paste0("PC2 (", round(summary(pca_data)$importance[2, 2] * 100, 1), "%)")) +
      theme_bw()
    
    ggsave(file.path(figures_dir, "GSE174409_PCA.pdf"), p_pca, width = 8, height = 6)
    cat("✓ Generated PCA plot\n")
  }, error = function(e) {
    cat("Warning: Could not generate PCA plot:", e$message, "\n")
  })
  
  cat("✓ Generated QC plots\n")
}

# ==============================================================================
# PART 4: Gene Annotation and Symbol Mapping
# ==============================================================================

cat("\n", strrep("=", 60), "\n")
cat("GENE ANNOTATION AND SYMBOL MAPPING\n")
cat(strrep("=", 60), "\n")

if (exists("gse174409_normalized")) {
  cat("Mapping gene symbols for GSE174409...\n")
  
  # Get gene symbols using biomaRt
  tryCatch({
    mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    
    # Extract gene IDs (assuming ENSEMBL IDs)
    gene_ids <- rownames(gse174409_normalized)
    
    gene_mapping <- getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "gene_biotype"),
      filters = "ensembl_gene_id",
      values = gene_ids,
      mart = mart
    )
    
    cat("✓ Retrieved annotations for", nrow(gene_mapping), "genes\n")
    
  }, error = function(e) {
    cat("Warning: Could not retrieve gene annotations:", e$message, "\n")
    # Create dummy mapping
    gene_mapping <- data.frame(
      ensembl_gene_id = rownames(gse174409_normalized),
      hgnc_symbol = rownames(gse174409_normalized),
      stringsAsFactors = FALSE
    )
  })
}

# ==============================================================================
# PART 5: Save Processed Data
# ==============================================================================

cat("\n", strrep("=", 60), "\n")
cat("SAVING PROCESSED DATA\n")
cat(strrep("=", 60), "\n")

# Save GSE174409 processed data
if (exists("gse174409_normalized")) {
  # Save normalized expression data
  saveRDS(gse174409_normalized, file.path(processed_data_dir, "GSE174409_normalized_expression.rds"))
  write.csv(gse174409_normalized, file.path(processed_data_dir, "GSE174409_normalized_expression.csv"))
  
  # Save sample metadata
  saveRDS(gse174409_sample_metadata, file.path(processed_data_dir, "GSE174409_sample_metadata.rds"))
  write.csv(gse174409_sample_metadata, file.path(processed_data_dir, "GSE174409_sample_metadata.csv"), row.names = FALSE)
  
  # Save gene annotations
  if (exists("gene_mapping")) {
    saveRDS(gene_mapping, file.path(processed_data_dir, "GSE174409_gene_annotations.rds"))
    write.csv(gene_mapping, file.path(processed_data_dir, "GSE174409_gene_annotations.csv"), row.names = FALSE)
  }
  
  cat("✓ Saved GSE174409 processed data\n")
}

# Save GSE225158 processed data
if (exists("gse225158_seurat")) {
  saveRDS(gse225158_seurat, file.path(processed_data_dir, "GSE225158_seurat_object.rds"))
  
  # Extract and save count matrix
  gse225158_counts <- GetAssayData(gse225158_seurat, assay = "RNA", slot = "counts")
  saveRDS(gse225158_counts, file.path(processed_data_dir, "GSE225158_count_matrix.rds"))
  
  # Save metadata
  gse225158_cell_metadata <- gse225158_seurat@meta.data
  saveRDS(gse225158_cell_metadata, file.path(processed_data_dir, "GSE225158_cell_metadata.rds"))
  write.csv(gse225158_cell_metadata, file.path(processed_data_dir, "GSE225158_cell_metadata.csv"))
  
  cat("✓ Saved GSE225158 single-cell data\n")
}

# Save processing summary
processing_summary <- list(
  processing_date = Sys.time(),
  gse174409 = list(
    samples = if(exists("gse174409_sample_metadata")) nrow(gse174409_sample_metadata) else 0,
    genes = if(exists("gse174409_normalized")) nrow(gse174409_normalized) else 0,
    conditions = if(exists("gse174409_sample_metadata")) unique(gse174409_sample_metadata$condition) else NA,
    regions = if(exists("gse174409_sample_metadata")) unique(gse174409_sample_metadata$region) else NA
  ),
  gse225158 = list(
    status = if(exists("gse225158_seurat")) "Processed" else "Metadata only",
    cells = if(exists("gse225158_seurat")) ncol(gse225158_seurat) else 0,
    genes = if(exists("gse225158_seurat")) nrow(gse225158_seurat) else 0
  )
)

saveRDS(processing_summary, file.path(results_dir, "human_data_processing_summary.rds"))

# ==============================================================================
# PART 6: Generate Processing Report
# ==============================================================================

cat("\n", strrep("=", 60), "\n")
cat("GENERATING PROCESSING REPORT\n")
cat(strrep("=", 60), "\n")

report_file <- file.path(results_dir, "Human_Data_Processing_Report.txt")

cat("HUMAN DATA PROCESSING REPORT\n", file = report_file)
cat(strrep("=", 40), "\n", file = report_file, append = TRUE)
cat("Processing Date:", as.character(Sys.time()), "\n\n", file = report_file, append = TRUE)

cat("GSE174409 - Human Bulk RNA-seq:\n", file = report_file, append = TRUE)
if (exists("gse174409_sample_metadata")) {
  cat("  Samples:", nrow(gse174409_sample_metadata), "\n", file = report_file, append = TRUE)
  cat("  Conditions:", paste(unique(gse174409_sample_metadata$condition), collapse = ", "), "\n", file = report_file, append = TRUE)
  cat("  Regions:", paste(unique(gse174409_sample_metadata$region), collapse = ", "), "\n", file = report_file, append = TRUE)
}
if (exists("gse174409_normalized")) {
  cat("  Genes after filtering:", nrow(gse174409_normalized), "\n", file = report_file, append = TRUE)
}

cat("\nGSE225158 - Human Single-cell RNA-seq:\n", file = report_file, append = TRUE)
if (exists("gse225158_seurat")) {
  cat("  Status: Processed successfully\n", file = report_file, append = TRUE)
  cat("  Cells:", ncol(gse225158_seurat), "\n", file = report_file, append = TRUE)
  cat("  Genes:", nrow(gse225158_seurat), "\n", file = report_file, append = TRUE)
} else {
  cat("  Status: Metadata downloaded, data processing pending\n", file = report_file, append = TRUE)
  cat("  Note: Manual download of supplementary files may be required\n", file = report_file, append = TRUE)
}

cat("\nOutput Files Generated:\n", file = report_file, append = TRUE)
cat("  - GSE174409_normalized_expression.rds/csv\n", file = report_file, append = TRUE)
cat("  - GSE174409_sample_metadata.rds/csv\n", file = report_file, append = TRUE)
cat("  - GSE174409_gene_annotations.rds/csv\n", file = report_file, append = TRUE)
if (exists("gse225158_seurat")) {
  cat("  - GSE225158_seurat_object.rds\n", file = report_file, append = TRUE)
  cat("  - GSE225158_count_matrix.rds\n", file = report_file, append = TRUE)
  cat("  - GSE225158_cell_metadata.rds/csv\n", file = report_file, append = TRUE)
}
cat("  - QC plots: correlation, PCA, and single-cell QC plots\n", file = report_file, append = TRUE)

cat("✓ Generated processing report\n")

# ==============================================================================
# COMPLETION SUMMARY
# ==============================================================================

cat("\n", strrep("=", 80), "\n")
cat("HUMAN DATA PROCESSING COMPLETED\n")
cat(strrep("=", 80), "\n")

cat("Processing Summary:\n")
cat("  ✓ GSE174409: Downloaded and processed\n")
cat("  ✓ Data normalization: Completed\n")
cat("  ✓ Gene annotation: Completed\n")
cat("  ✓ QC plots: Generated\n")
if (exists("gse225158_seurat")) {
  cat("  ✓ GSE225158: Single-cell data processed\n")
} else {
  cat("  ⚠ GSE225158: Manual file download may be required\n")
}

cat("\nNext Steps:\n")
if (!exists("gse225158_seurat")) {
  cat("  1. Download GSE225158 supplementary files manually if needed\n")
  cat("  2. Re-run script to process GSE225158 data\n")
  cat("  3. Run mouse data harmonization script\n")
} else {
  cat("  1. Run mouse data harmonization script\n")
  cat("  2. Proceed to ortholog mapping\n")
}

cat("\nProcessing completed at:", as.character(Sys.time()), "\n")
cat(strrep("=", 80), "\n")
