# GSE207128 Data Processing Script
# This script processes features, barcodes, and matrix files from the GSE207128 dataset

# Load required libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(fs)

# Set working directory to the GSE207128 folder
# Uncomment and modify if needed
# setwd("/path/to/GSE207128")

# Function to read 10X format data from a directory
read_10X_data <- function(data_dir) {
  cat("Processing directory:", data_dir, "\n")

  # Check if files exist
  feature_path <- file.path(data_dir, "features.tsv.gz")
  barcode_path <- file.path(data_dir, "barcodes.tsv.gz")
  matrix_path <- file.path(data_dir, "matrix.mtx.gz")

  if (!file.exists(feature_path) || !file.exists(barcode_path) || !file.exists(matrix_path)) {
    warning("Required files not found in directory: ", data_dir)
    return(NULL)
  }

  # Read the data
  features <- read.delim(feature_path, header = FALSE, stringsAsFactors = FALSE)
  barcodes <- read.delim(barcode_path, header = FALSE, stringsAsFactors = FALSE)
  matrix <- readMM(matrix_path)

  # Create row names for the matrix (genes)
  rownames(matrix) <- features$V1

  # Create column names for the matrix (cells)
  colnames(matrix) <- barcodes$V1

  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = matrix,
                                  project = basename(data_dir),
                                  min.cells = 3,
                                  min.features = 200)

  # Add metadata
  seurat_obj$run <- basename(data_dir)

  return(seurat_obj)
}

# Find all run directories in the GSE207128 folder
# Assumes directories are organized as GSE207128/[run_id]/
run_dirs <- list.dirs(".", recursive = FALSE)
cat("Found", length(run_dirs), "potential run directories\n")

# Process each run directory
seurat_objects <- list()
for (dir in run_dirs) {
  seurat_obj <- read_10X_data(dir)
  if (!is.null(seurat_obj)) {
    seurat_objects[[basename(dir)]] <- seurat_obj
    cat("Successfully processed", basename(dir), "\n")
  }
}

# Report summary
cat("Successfully processed", length(seurat_objects), "run directories\n")

# Save the processed data
saveRDS(seurat_objects, "GSE207128_processed_seurat_objects.rds")

# If you want to merge all objects into one (optional)
if (length(seurat_objects) > 0) {
  merged_seurat <- merge(seurat_objects[[1]],
                        y = seurat_objects[-1],
                        add.cell.ids = names(seurat_objects),
                        project = "GSE207128")

  # Save the merged object
  saveRDS(merged_seurat, "GSE207128_merged_seurat_object.rds")

  # Quick QC summary
  cat("Merged data contains", nrow(merged_seurat), "genes and", ncol(merged_seurat), "cells\n")

  # Plot distribution of cells per sample
  p <- ggplot(data.frame(Sample = merged_seurat$run), aes(x = Sample)) +
    geom_bar() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Number of cells per sample", y = "Cell count")

  ggsave("cells_per_sample.pdf", p, width = 10, height = 6)
}

cat("Processing complete!\n")