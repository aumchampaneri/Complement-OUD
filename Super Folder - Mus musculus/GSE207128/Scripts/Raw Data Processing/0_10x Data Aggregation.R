#===============================================================================
# 10x Genomics Data Aggregation Script
#===============================================================================
# Purpose: Read raw 10x Genomics single-cell RNA-seq data from multiple samples
#          and conditions, then aggregate into a single Seurat object for analysis
#
# Author: Aum Champaneri
# Created: 2025MAY21
# Last Modified: 2025MAY24
#
# Input Data Structure:
#   Raw Data/
#   ├── Naive/          (*_N1_*, *_N2_*, *_N3_* files)
#   ├── Dependent/      (*_D1_*, *_D2_*, *_D3_* files)
#   └── Withdrawal/     (*_W1_*, *_W2_*, *_W3_* files)
#
# Each condition folder contains:
#   - *_matrix.mtx.gz    (gene expression count matrix)
#   - *_features.tsv.gz  (gene information)
#   - *_barcodes.tsv.gz  (cell barcodes)
#
# Output: combined_seurat.rds - Merged Seurat object with all samples
#===============================================================================

# Load required libraries
library(Seurat)      # Single-cell RNA-seq analysis
library(dplyr)       # Data manipulation
library(Matrix)      # Sparse matrix operations
library(ggplot2)     # Plotting (for potential future use)
library(fs)          # Modern file system operations
library(R.utils)     # File compression utilities (gunzip)

#===============================================================================
# SETUP AND CONFIGURATION
#===============================================================================

# Set working directories
setwd("/Users/aumchampaneri/PycharmProjects/Complement-OUD/Super Folder - Mus musculus/GSE207128")
input_dir <- "Raw Data"                    # Directory containing raw 10x data
output_dir <- "Outputs/Processed Data"     # Output directory for processed data

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define experimental conditions and their file prefixes
conditions <- c("Naive", "Dependent", "Withdrawal")  # Three experimental conditions
condition_prefixes <- c("N", "D", "W")               # Corresponding file prefixes

#===============================================================================
# HELPER FUNCTIONS
#===============================================================================

# Function: decompress_if_needed
# Purpose: Automatically decompress .gz files only if uncompressed version doesn't exist
# Input: file_path (string) - path to potentially compressed file
# Output: path to uncompressed file
# Note: Preserves original compressed file
decompress_if_needed <- function(file_path) {
  if (grepl("\\.gz$", file_path)) {
    # Generate path for decompressed file
    decompressed <- sub("\\.gz$", "", file_path)
    
    # Only decompress if uncompressed version doesn't exist
    if (!file.exists(decompressed)) {
      cat("    Decompressing:", basename(file_path), "\n")
      gunzip(file_path, remove = FALSE, overwrite = TRUE)
    }
    return(decompressed)
  } else {
    # File is already uncompressed
    return(file_path)
  }
}

#===============================================================================
# DATA PROCESSING INITIALIZATION
#===============================================================================

# Initialize containers for processing
seurat_list <- list()    # List to store individual Seurat objects
sample_counter <- 1      # Counter for indexing Seurat objects
temp_files <- c()        # Track temporary decompressed files for cleanup

cat("Starting 10x Genomics data aggregation...\n")
cat("Processing", length(conditions), "conditions:", paste(conditions, collapse=", "), "\n\n")

#===============================================================================
# MAIN PROCESSING LOOP
#===============================================================================

# Process each experimental condition
for (i in seq_along(conditions)) {
  condition <- conditions[i]
  prefix <- condition_prefixes[i]
  condition_path <- path(input_dir, condition)
  
  cat(sprintf("Processing condition: %s (prefix: %s)\n", condition, prefix))
  
  # Auto-detect sample numbers for this condition
  # Search for files with pattern: *_[PREFIX][NUMBER]_*
  all_files <- dir_ls(condition_path)
  sample_nums <- unique(gsub(paste0(".*_", prefix, "(\\d+)_.*"), "\\1", basename(all_files)))
  sample_nums <- sample_nums[grep("^\\d+$", sample_nums)]  # Keep only numeric matches
  
  cat("  Found", length(sample_nums), "samples:", paste(sort(sample_nums), collapse=", "), "\n")

  # Process each sample within the condition
  for (sample_num in sample_nums) {
    sample_id <- paste0(prefix, sample_num)  # Create sample ID (e.g., N1, D2, W3)
    cat(sprintf("  Processing sample: %s\n", sample_id))

    # Wrap sample processing in error handling to continue if one sample fails
    tryCatch({
      
      #-------------------------------------------------------------------------
      # STEP 1: Locate required 10x files for this sample
      #-------------------------------------------------------------------------
      
      # Define file patterns for 10x Genomics standard output
      matrix_pattern <- paste0("*_", sample_id, "_matrix.mtx.gz")
      features_pattern <- paste0("*_", sample_id, "_features.tsv.gz")
      barcodes_pattern <- paste0("*_", sample_id, "_barcodes.tsv.gz")

      # Find files matching patterns
      matrix_gz <- dir_ls(condition_path, glob = matrix_pattern)
      if (length(matrix_gz) == 0) stop(paste("No matrix file found for", sample_id))
      matrix_gz <- matrix_gz[1]  # Take first match if multiple files found

      features_gz <- dir_ls(condition_path, glob = features_pattern)
      if (length(features_gz) == 0) stop(paste("No features file found for", sample_id))
      features_gz <- features_gz[1]

      barcodes_gz <- dir_ls(condition_path, glob = barcodes_pattern)
      if (length(barcodes_gz) == 0) stop(paste("No barcodes file found for", sample_id))
      barcodes_gz <- barcodes_gz[1]

      #-------------------------------------------------------------------------
      # STEP 2: Decompress files if needed
      #-------------------------------------------------------------------------
      
      matrix_file <- decompress_if_needed(matrix_gz)
      features_file <- decompress_if_needed(features_gz)
      barcodes_file <- decompress_if_needed(barcodes_gz)

      # Track temporary files for later cleanup
      if (matrix_file != matrix_gz) temp_files <- c(temp_files, matrix_file)
      if (features_file != features_gz) temp_files <- c(temp_files, features_file)
      if (barcodes_file != barcodes_gz) temp_files <- c(temp_files, barcodes_file)

      #-------------------------------------------------------------------------
      # STEP 3: Read 10x data into R
      #-------------------------------------------------------------------------
      
      # Read the sparse matrix and associated metadata
      # ReadMtx expects: matrix file, features/genes file, cell barcodes file
      expression_matrix <- ReadMtx(
        mtx = matrix_file,      # Gene expression count matrix
        features = features_file, # Gene names and IDs
        cells = barcodes_file    # Cell barcodes
      )

      #-------------------------------------------------------------------------
      # STEP 4: Create Seurat object with quality control
      #-------------------------------------------------------------------------
      
      # Create Seurat object with basic filtering
      seurat_obj <- CreateSeuratObject(
        counts = expression_matrix,
        project = paste0(condition, "_", sample_id),  # Project name
        min.cells = 3,      # Include features detected in at least 3 cells
        min.features = 200  # Include cells with at least 200 detected features
      )

      #-------------------------------------------------------------------------
      # STEP 5: Add metadata annotations
      #-------------------------------------------------------------------------
      
      # Add experimental condition (Naive, Dependent, Withdrawal)
      seurat_obj$condition <- condition
      
      # Add sample identifier (N1, N2, D1, etc.)
      seurat_obj$sample_id <- sample_id

      #-------------------------------------------------------------------------
      # STEP 6: Store processed object
      #-------------------------------------------------------------------------
      
      # Add to list of successfully processed samples
      seurat_list[[sample_counter]] <- seurat_obj
      sample_counter <- sample_counter + 1

      cat(sprintf("  Successfully processed %s: %d cells, %d features\n", 
                  sample_id, ncol(seurat_obj), nrow(seurat_obj)))

    }, error = function(e) {
      # Log error but continue processing other samples
      cat(sprintf("  Error processing %s: %s\n", sample_id, e$message))
    })
  }
  cat("\n")  # Add spacing between conditions
}

#===============================================================================
# DATA MERGING AND OUTPUT
#===============================================================================

# Merge all successfully processed samples into a single Seurat object
if (length(seurat_list) > 0) {
  cat("Merging", length(seurat_list), "samples into a single Seurat object...\n")

  # Create unique cell identifiers for merging
  # Format: {condition}_{sample_id} (e.g., Naive_N1, Dependent_D2)
  # This prevents cell barcode conflicts between samples
  cell_ids <- sapply(seurat_list, function(x) paste0(x$condition[1], "_", x$sample_id[1]))

  # Merge Seurat objects
  # First object is the base, remaining objects are merged into it
  combined_seurat <- merge(seurat_list[[1]], 
                         y = seurat_list[-1],          # All objects except the first
                         add.cell.ids = cell_ids)      # Add prefixes to cell names

  #-----------------------------------------------------------------------------
  # Save processed data
  #-----------------------------------------------------------------------------
  
  output_file <- file.path(output_dir, "combined_seurat.rds")
  saveRDS(combined_seurat, file = output_file)
  cat("Seurat object saved to:", output_file, "\n\n")

  #-----------------------------------------------------------------------------
  # Generate summary report
  #-----------------------------------------------------------------------------
  
  cat("=== PROCESSING SUMMARY ===\n")
  cat("Structure of combined Seurat object:\n")
  print(combined_seurat)
  
  cat("\nDimensions (features x cells):", paste(dim(combined_seurat), collapse=" x "), "\n")
  cat("Available assays:", paste(Assays(combined_seurat), collapse=", "), "\n")
  
  cat("\nSample distribution across conditions:\n")
  print(table(combined_seurat$condition, combined_seurat$sample_id))
  
  cat("\nCell counts by condition:\n")
  print(table(combined_seurat$condition))
  
} else {
  cat("ERROR: No samples were successfully processed.\n")
  cat("Please check:\n")
  cat("  1. File paths and directory structure\n")
  cat("  2. File naming conventions\n")
  cat("  3. File permissions\n")
}

#===============================================================================
# CLEANUP
#===============================================================================

# Remove temporary decompressed files to save disk space
if (length(temp_files) > 0) {
  cat("\nCleaning up temporary files...\n")
  files_removed <- 0
  
  for (file in temp_files) {
    if (file.exists(file)) {
      file.remove(file)
      files_removed <- files_removed + 1
    }
  }
  
  cat("Cleanup complete. Removed", files_removed, "temporary files.\n")
}

cat("\n=== SCRIPT COMPLETED ===\n")

#===============================================================================
# USAGE NOTES
#===============================================================================
# 
# To use the resulting merged object:
#   combined_seurat <- readRDS("Outputs/Processed Data/combined_seurat.rds")
#
# The object contains:
#   - Raw count data in the "RNA" assay
#   - Metadata: condition (Naive/Dependent/Withdrawal) and sample_id (N1, D2, etc.)
#   - Cell names formatted as: {condition}_{sample_id}_{original_barcode}
#
# Next steps typically include:
#   1. Quality control analysis
#   2. Normalization
#   3. Dimensionality reduction
#   4. Clustering
#   5. Differential expression analysis
#
#===============================================================================