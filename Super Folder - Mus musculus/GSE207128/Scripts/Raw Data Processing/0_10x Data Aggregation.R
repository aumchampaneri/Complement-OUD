# Read the raw data from the 10x Genomics output and aggregate it into a single Seurat object
# Load required libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(fs)
library(R.utils)  # for gunzip

# Set working directories
setwd("/Users/aumchampaneri/PycharmProjects/Complement-OUD/Super Folder - Mus musculus/GSE207128")
input_dir <- "Raw Data"
output_dir <- "Outputs/Processed Data"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define the three conditions
conditions <- c("Naive", "Dependent", "Withdrawal")
condition_prefixes <- c("N", "D", "W")  # Prefixes in filenames

# Helper to decompress if not already done
decompress_if_needed <- function(file_path) {
  if (grepl("\\.gz$", file_path)) {
    decompressed <- sub("\\.gz$", "", file_path)
    if (!file.exists(decompressed)) {
      gunzip(file_path, remove = FALSE, overwrite = TRUE)
    }
    return(decompressed)
  } else {
    return(file_path)
  }
}

# Initialize Seurat list
seurat_list <- list()
sample_counter <- 1
temp_files <- c() # Track temporary files for cleanup

for (i in seq_along(conditions)) {
  condition <- conditions[i]
  prefix <- condition_prefixes[i]
  condition_path <- path(input_dir, condition)
  cat(sprintf("Processing condition: %s (prefix: %s)\n", condition, prefix))

  # Find sample numbers (1, 2, 3) for this condition
  all_files <- dir_ls(condition_path)
  sample_nums <- unique(gsub(paste0(".*_", prefix, "(\\d+)_.*"), "\\1", basename(all_files)))
  sample_nums <- sample_nums[grep("^\\d+$", sample_nums)]  # Keep only numeric matches

  for (sample_num in sample_nums) {
    sample_id <- paste0(prefix, sample_num)
    cat(sprintf("  Processing sample: %s\n", sample_id))

    tryCatch({
      # Get specific files for this sample
      matrix_pattern <- paste0("*_", sample_id, "_matrix.mtx.gz")
      features_pattern <- paste0("*_", sample_id, "_features.tsv.gz")
      barcodes_pattern <- paste0("*_", sample_id, "_barcodes.tsv.gz")

      matrix_gz <- dir_ls(condition_path, glob = matrix_pattern)
      if (length(matrix_gz) == 0) stop(paste("No matrix file found for", sample_id))
      matrix_gz <- matrix_gz[1]  # Take first if multiple matches

      features_gz <- dir_ls(condition_path, glob = features_pattern)
      if (length(features_gz) == 0) stop(paste("No features file found for", sample_id))
      features_gz <- features_gz[1]

      barcodes_gz <- dir_ls(condition_path, glob = barcodes_pattern)
      if (length(barcodes_gz) == 0) stop(paste("No barcodes file found for", sample_id))
      barcodes_gz <- barcodes_gz[1]

      # Decompress files if needed
      matrix_file <- decompress_if_needed(matrix_gz)
      features_file <- decompress_if_needed(features_gz)
      barcodes_file <- decompress_if_needed(barcodes_gz)

      # Track temp files
      if (matrix_file != matrix_gz) temp_files <- c(temp_files, matrix_file)
      if (features_file != features_gz) temp_files <- c(temp_files, features_file)
      if (barcodes_file != barcodes_gz) temp_files <- c(temp_files, barcodes_file)

      # Read matrix
      expression_matrix <- ReadMtx(
        mtx = matrix_file,
        features = features_file,
        cells = barcodes_file
      )

      # Create Seurat object
      seurat_obj <- CreateSeuratObject(
        counts = expression_matrix,
        project = paste0(condition, "_", sample_id),
        min.cells = 3,
        min.features = 200
      )

      # Annotate metadata
      seurat_obj$condition <- condition
      seurat_obj$sample_id <- sample_id

      # Store in list
      seurat_list[[sample_counter]] <- seurat_obj
      sample_counter <- sample_counter + 1

      cat(sprintf("  Successfully processed %s\n", sample_id))

    }, error = function(e) {
      cat(sprintf("  Error processing %s: %s\n", sample_id, e$message))
    })
  }
}

# Merge all into a single Seurat object
if (length(seurat_list) > 0) {
  cat("Merging", length(seurat_list), "samples into a single Seurat object...\n")

  # Create proper cell IDs for merging
  cell_ids <- sapply(seurat_list, function(x) paste0(x$condition[1], "_", x$sample_id[1]))

  combined_seurat <- merge(seurat_list[[1]], y = seurat_list[-1],
                         add.cell.ids = cell_ids)

  # Save Seurat object
  saveRDS(combined_seurat, file = file.path(output_dir, "combined_seurat.rds"))
  cat("Seurat object saved to:", file.path(output_dir, "combined_seurat.rds"), "\n")

  # Print structure of the combined object
  cat("Structure of combined Seurat object:\n")
  print(combined_seurat)
  cat("Dimensions (features x cells):", paste(dim(combined_seurat), collapse=" x "), "\n")
  cat("Available assays:", paste(Assays(combined_seurat), collapse=", "), "\n")
  cat("Sample distribution:\n")
  print(table(combined_seurat$condition, combined_seurat$sample_id))
} else {
  cat("Error: No samples were successfully processed.\n")
}

# Clean up temporary files
if (length(temp_files) > 0) {
  cat("Cleaning up temporary files...\n")
  for (file in temp_files) {
    if (file.exists(file)) {
      file.remove(file)
    }
  }
  cat("Cleanup complete.\n")
}