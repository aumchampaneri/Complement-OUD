# Analysis Outputs_0 Data Processing Script
# This script processes features, barcodes, and matrix files from the Analysis Outputs_0 dataset

# Load required libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(fs)

# Set working directory to the Analysis Outputs_0 folder
setwd("/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE207128/GSE207128_RAW")

# Create output directory if it doesn't exist
output_dir <- "../processed_data"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Modify the read_10X_data function to handle the specific file naming pattern
read_10X_data <- function(sample_id) {
  cat("Processing sample:", sample_id, "\n")

  # Construct file paths with the sample ID prefix
  feature_path <- file.path("", paste0(sample_id, "_features.tsv.gz"))
  barcode_path <- file.path("", paste0(sample_id, "_barcodes.tsv.gz"))
  matrix_path <- file.path("", paste0(sample_id, "_matrix.mtx.gz"))

  if (!file.exists(feature_path) || !file.exists(barcode_path) || !file.exists(matrix_path)) {
    warning("Required files not found for sample: ", sample_id)
    return(NULL)
  }

  # Read the data
  features <- read.delim(feature_path, header = FALSE, stringsAsFactors = FALSE)
  barcodes <- read.delim(barcode_path, header = FALSE, stringsAsFactors = FALSE)
  matrix <- Matrix::readMM(matrix_path)

  # Create row names for the matrix (genes)
  rownames(matrix) <- features$V1

  # Create column names for the matrix (cells)
  colnames(matrix) <- barcodes$V1

  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = matrix,
                                  project = sample_id,
                                  min.cells = 3,
                                  min.features = 200)

  # Add metadata
  seurat_obj$sample <- sample_id

  # Free up memory
  rm(matrix, features, barcodes)
  gc()

  return(seurat_obj)
}

# Get list of all files in the directory
all_files <- list.files("")

# Extract unique sample IDs by looking at barcode files
barcode_files <- grep("_barcodes.tsv.gz$", all_files, value = TRUE)
sample_ids <- gsub("_barcodes.tsv.gz$", "", barcode_files)

cat("Found", length(sample_ids), "potential samples\n")

# Process each sample
seurat_objects <- list()
for (sample_id in sample_ids) {
  tryCatch({
    seurat_obj <- read_10X_data(sample_id)
    if (!is.null(seurat_obj)) {
      seurat_objects[[sample_id]] <- seurat_obj
      cat("Successfully processed", sample_id, "\n")

      # Save each object individually to free up memory
      saveRDS(seurat_obj, file.path(output_dir, paste0(sample_id, "_seurat_object.rds")))

      # Remove objects to free memory
      rm(seurat_obj)
      gc()
    }
  }, error = function(e) {
    cat("Error processing sample", sample_id, ":", e$message, "\n")
  })
}

# Save the sample IDs that were successfully processed
saveRDS(names(seurat_objects), file.path(output_dir, "processed_sample_ids.rds"))

cat("Processing complete! Processed", length(seurat_objects), "samples\n")