# Script to explore data structure before pseudobulking

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# Load the filtered Seurat object
seurat_filtered <- readRDS("/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE207128/harmony_results/singleR_annotation/seurat_filtered.rds")

# Basic object dimensions
cat("==== Basic Information ====\n")
cat("Number of cells:", ncol(seurat_filtered), "\n")
cat("Number of genes:", nrow(seurat_filtered), "\n")

# Examine metadata structure
cat("\n==== Metadata Structure ====\n")
cat("Available metadata columns:\n")
print(colnames(seurat_filtered@meta.data))

# Print first few rows of metadata
cat("\nSample of metadata (first 5 rows):\n")
print(head(seurat_filtered@meta.data, 5))

# Cell type information
cat("\n==== Cell Type Information ====\n")
cat("Cell types from SingleR annotation:\n")
print(table(seurat_filtered$singler_label))

# Check for sample or condition information
cat