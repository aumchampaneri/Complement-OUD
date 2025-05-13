# Script to add gene symbols to Seurat object with Ensembl IDs

# Load required libraries
library(Seurat)
library(biomaRt)
library(dplyr)

# Load the filtered Seurat object
seurat_filtered <- readRDS("/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE207128/harmony_results/singleR_annotation/seurat_filtered.rds")

# Basic data inspection
cat("Data dimensions:", dim(seurat_filtered), "\n")
cat("First few gene IDs:", head(rownames(seurat_filtered)), "\n")

# Check if gene names are Ensembl IDs
ensembl_pattern <- "^ENSMUS"
gene_ids <- rownames(seurat_filtered)
is_ensembl <- any(grepl(ensembl_pattern, gene_ids))

if(!is_ensembl) {
  cat("Gene IDs don't appear to be Ensembl IDs. No conversion needed.\n")
} else {
  cat("Converting Ensembl IDs to gene symbols...\n")

  # Connect to BioMart
  cat("Connecting to BioMart...\n")
  mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

  # Get mapping between Ensembl IDs and gene symbols
  cat("Fetching gene symbol mapping...\n")
  mapping <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "description"),
    filters = "ensembl_gene_id",
    values = gene_ids,
    mart = mart
  )

  cat("Found mappings for", nrow(mapping), "out of", length(gene_ids), "Ensembl IDs\n")

  # Create mapping dataframe
  gene_mapping <- data.frame(
    ensembl_id = gene_ids,
    gene_symbol = NA,
    description = NA,
    stringsAsFactors = FALSE
  )
  rownames(gene_mapping) <- gene_ids

  # Fill in gene symbols and descriptions
  matched_idx <- match(mapping$ensembl_gene_id, gene_ids)
  gene_mapping$gene_symbol[matched_idx] <- mapping$external_gene_name
  gene_mapping$description[matched_idx] <- mapping$description

  # Use Ensembl ID for genes without a symbol
  missing_symbols <- is.na(gene_mapping$gene_symbol) | gene_mapping$gene_symbol == ""
  gene_mapping$gene_symbol[missing_symbols] <- gene_ids[missing_symbols]

  # Add the mapping to the Seurat object
  seurat_filtered@misc$gene_mapping <- gene_mapping

  # Save outputs
  output_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE207128/harmony_results/singleR_annotation"
  write.csv(gene_mapping, file.path(output_dir, "gene_id_symbol_mapping.csv"))
  saveRDS(seurat_filtered, file.path(output_dir, "seurat_filtered_with_mapping.rds"))

  # Print sample of mapping
  cat("\nSample of gene ID to symbol mapping:\n")
  print(head(gene_mapping, 10))
}

# Check other metadata
cat("\n==== Sample Metadata ====\n")
meta_cols <- colnames(seurat_filtered@meta.data)
cat("Available metadata columns:", paste(meta_cols, collapse=", "), "\n\n")

# Check for potential pseudobulking variables
if("orig.ident" %in% meta_cols) {
  cat("Sample identifiers (orig.ident):\n")
  print(table(seurat_filtered$orig.ident))
}

cat("\nData ready for pseudobulking with gene symbols now accessible!\n")