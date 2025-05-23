# Analysis Outputs_0 - Pseudobulk Analysis by Treatment Condition
# Optimized for DESeq2 compatibility

library(Seurat)
library(dplyr)
library(readr)
library(Matrix)

#----- 1. Set up paths and load data -----#
base_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE207128"
seurat_path <- file.path(base_dir, "harmony_results/singleR_annotation/seurat_filtered_with_mapping.rds")
output_dir <- file.path(base_dir, "pseudobulk_results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load Seurat object
cat("Loading Seurat object...\n")
seurat_obj <- readRDS(seurat_path)
gene_mapping <- seurat_obj@misc$gene_mapping

#----- 2. Create treatment condition variable -----#
cat("Creating treatment condition variable...\n")

# Map samples to conditions
if("sample_names" %in% colnames(seurat_obj@meta.data)) {
  treatment_map <- list(
    # Normal/control samples
    "GSM6278819_N1" = "Normal",
    "GSM6278820_N2" = "Normal",
    "GSM6278821_N3" = "Normal",
    # Opioid dependence samples
    "GSM6278822_D1" = "OpioidDependence",
    "GSM6278823_D2" = "OpioidDependence",
    "GSM6278824_D3" = "OpioidDependence"
  )

  # Add treatment and ensure sample ID is preserved
  seurat_obj$treatment_condition <- sapply(seurat_obj$sample_names, function(x) treatment_map[[x]])
  seurat_obj$sample_id <- seurat_obj$sample_names
} else if("orig.ident" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj$treatment_condition <- ifelse(
    grepl("D[1-3]$", seurat_obj$orig.ident),
    "OpioidDependence",
    ifelse(grepl("N[1-3]$", seurat_obj$orig.ident), "Normal", NA)
  )
  seurat_obj$sample_id <- seurat_obj$orig.ident
} else {
  stop("Could not find sample identifiers")
}

# Show distribution
cat("\nTreatment condition distribution:\n")
print(table(seurat_obj$treatment_condition))
cat("\nSample distribution:\n")
print(table(seurat_obj$sample_id))

#----- 3. Filter very low-expressed genes -----#
cat("\nFiltering low-expressed genes...\n")
counts <- GetAssayData(seurat_obj, slot = "counts", assay = "RNA")
genes_to_keep <- rownames(counts)[rowSums(counts > 0) >= 10]
cat("Keeping", length(genes_to_keep), "of", nrow(counts), "genes after filtering\n")

#----- 4. Define improved pseudobulk function -----#
create_pseudobulk_for_deseq2 <- function(seurat_obj, genes_to_keep = NULL) {
  # Get raw count data
  counts <- GetAssayData(seurat_obj, slot = "counts", assay = "RNA")

  # Filter genes if requested
  if (!is.null(genes_to_keep)) {
    counts <- counts[genes_to_keep, ]
  }

  # Prepare aggregation groups - by sample_id + cell type
  groups <- paste(seurat_obj$sample_id, seurat_obj$singler_label, sep = "_")
  unique_groups <- unique(groups)

  # Initialize pseudobulk matrix
  pseudobulk <- matrix(0, nrow = nrow(counts), ncol = length(unique_groups))
  rownames(pseudobulk) <- rownames(counts)
  colnames(pseudobulk) <- unique_groups

  # Sum counts for each group
  cat("Creating pseudobulk profiles by summing counts within groups:\n")
  for(group in unique_groups) {
    cells <- which(groups == group)
    pseudobulk[, group] <- rowSums(counts[, cells, drop = FALSE])
    cat(sprintf("Group %s: %d cells, total counts: %.1f million\n",
               group, length(cells), sum(pseudobulk[, group])/1e6))
  }

  # Create metadata for the pseudobulk samples
  meta <- data.frame(
    sample_name = unique_groups,
    sample_id = gsub("^(.+?)_(.+)$", "\\1", unique_groups),
    celltype = gsub("^(.+?)_(.+)$", "\\2", unique_groups)
  )

  # Add treatment condition
  treatment_lookup <- setNames(
    seurat_obj$treatment_condition,
    seurat_obj$sample_id
  )
  meta$condition <- treatment_lookup[meta$sample_id]

  # Ensure pseudobulk columns and metadata rows match
  meta <- meta[match(colnames(pseudobulk), meta$sample_name), ]

  return(list(
    counts = pseudobulk,
    metadata = meta
  ))
}

#----- 5. Generate pseudobulk data -----#
cat("\nCreating pseudobulk data optimized for DESeq2...\n")
pseudobulk_data <- create_pseudobulk_for_deseq2(seurat_obj, genes_to_keep)

# Check the structure
cat("\nPseudobulk matrix dimensions:", dim(pseudobulk_data$counts), "\n")
cat("Pseudobulk metadata rows:", nrow(pseudobulk_data$metadata), "\n")

# Create additional columns for DESeq2
pseudobulk_data$metadata$celltype_condition <- paste(
  pseudobulk_data$metadata$celltype,
  pseudobulk_data$metadata$condition,
  sep = "_"
)

#----- 6. Map gene IDs to symbols -----#
# Map Ensembl IDs to gene symbols for retained genes
gene_symbols <- gene_mapping$gene_symbol[match(rownames(pseudobulk_data$counts),
                                            gene_mapping$ensembl_id)]
id_symbol_mapping <- data.frame(
  ensembl_id = rownames(pseudobulk_data$counts),
  gene_symbol = gene_symbols,
  stringsAsFactors = FALSE
)

#----- 7. Save results -----#
# Save the pseudobulk data
saveRDS(pseudobulk_data$counts,
       file.path(output_dir, "GSE207128_pseudobulk_counts_for_deseq2.rds"))

# Save the metadata
saveRDS(pseudobulk_data$metadata,
       file.path(output_dir, "GSE207128_pseudobulk_metadata_for_deseq2.rds"))
write_csv(pseudobulk_data$metadata,
         file.path(output_dir, "GSE207128_pseudobulk_metadata_for_deseq2.csv"))

# Save the gene mapping
write_csv(id_symbol_mapping,
         file.path(output_dir, "GSE207128_pseudobulk_gene_mapping.csv"))

# Also create backward-compatible versions for existing scripts
saveRDS(pseudobulk_data$counts,
       file.path(output_dir, "GSE207128_pseudobulk_celltype_by_treatment.rds"))

cat("\nPseudobulk generation complete. Files saved to:", output_dir, "\n")
cat("Key improvements for DESeq2 compatibility:\n")
cat("- Used raw counts instead of normalized data\n")
cat("- Summed counts rather than averaging\n")
cat("- Preserved biological replicates (per sample)\n")
cat("- Pre-filtered very low-expressed genes\n")
cat("- Created proper metadata for DESeq2 design formula\n")
cat("- Maintained integer counts required by DESeq2\n")