# GSE207128 - Pseudobulk Analysis by Treatment Condition
# This script creates pseudobulk profiles by opioid dependence status

library(Seurat)
library(dplyr)
library(readr)

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
# Extract sample information from sample_names
cat("Creating treatment condition variable...\n")

# Check if sample_names column exists
if("sample_names" %in% colnames(seurat_obj@meta.data)) {
  # Map samples to conditions based on the information provided
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

  # Create treatment condition variable
  seurat_obj$treatment_condition <- sapply(seurat_obj$sample_names,
                                          function(x) treatment_map[[x]])

  # Show distribution
  cat("\nTreatment condition distribution:\n")
  print(table(seurat_obj$treatment_condition))
} else {
  # If sample_names doesn't exist, look for other potential identifiers
  # like orig.ident that might contain the sample information
  if("orig.ident" %in% colnames(seurat_obj@meta.data)) {
    cat("\nSample identifiers (orig.ident):\n")
    print(table(seurat_obj$orig.ident))

    # Try to infer treatment from orig.ident
    seurat_obj$treatment_condition <- ifelse(
      grepl("D[1-3]$", seurat_obj$orig.ident),
      "OpioidDependence",
      ifelse(grepl("N[1-3]$", seurat_obj$orig.ident), "Normal", NA)
    )

    cat("\nInferred treatment condition distribution:\n")
    print(table(seurat_obj$treatment_condition, useNA = "ifany"))
  } else {
    stop("Could not find sample_names or orig.ident columns to determine treatment conditions")
  }
}

#----- 3. Define pseudobulk function -----#
create_pseudobulk <- function(seurat_obj, group_by_col) {
  # Get normalized expression data
  expr_matrix <- GetAssayData(seurat_obj, slot = "data", assay = "RNA")

  # Get grouping values
  groups <- seurat_obj@meta.data[[group_by_col]]

  # Initialize pseudobulk matrix
  unique_groups <- unique(groups)
  pseudobulk <- matrix(0, nrow = nrow(expr_matrix), ncol = length(unique_groups))
  rownames(pseudobulk) <- rownames(expr_matrix)
  colnames(pseudobulk) <- unique_groups

  # Calculate mean expression for each group
  for(group in unique_groups) {
    cells <- which(groups == group)
    pseudobulk[, group] <- rowMeans(expr_matrix[, cells, drop = FALSE])
    cat(sprintf("Group %s: %d cells\n", group, length(cells)))
  }

  return(pseudobulk)
}

#----- 4. Generate pseudobulk matrices -----#
# Create pseudobulk by treatment condition
cat("\nCreating pseudobulk profiles by treatment condition...\n")
pseudobulk_condition <- create_pseudobulk(seurat_obj, "treatment_condition")
cat("Pseudobulk matrix dimensions (by treatment condition):", dim(pseudobulk_condition), "\n")

# Optional: Create pseudobulk by cell type within each treatment condition
cat("\nCreating pseudobulk profiles by cell type within each treatment condition...\n")
seurat_obj$celltype_treatment <- paste(seurat_obj$singler_label,
                                      seurat_obj$treatment_condition, sep="_")
pseudobulk_celltype_condition <- create_pseudobulk(seurat_obj, "celltype_treatment")
cat("Pseudobulk matrix dimensions (celltype within treatment):",
    dim(pseudobulk_celltype_condition), "\n")

#----- 5. Map gene IDs to symbols -----#
# Map Ensembl IDs to gene symbols
gene_symbols <- gene_mapping$gene_symbol[match(rownames(pseudobulk_condition),
                                             gene_mapping$ensembl_id)]
names(gene_symbols) <- rownames(pseudobulk_condition)

# Create a data frame with both IDs and symbols
id_symbol_mapping <- data.frame(
  ensembl_id = rownames(pseudobulk_condition),
  gene_symbol = gene_symbols,
  stringsAsFactors = FALSE
)

#----- 6. Save results -----#
# Save the condition pseudobulk matrix
saveRDS(pseudobulk_condition,
        file.path(output_dir, "GSE207128_pseudobulk_by_treatment.rds"))

# Save the cell type within condition pseudobulk matrix
saveRDS(pseudobulk_celltype_condition,
        file.path(output_dir, "GSE207128_pseudobulk_celltype_by_treatment.rds"))

# Save the gene mapping
write_csv(id_symbol_mapping,
          file.path(output_dir, "GSE207128_pseudobulk_gene_mapping.csv"))

cat("\nPseudobulk generation complete. Files saved to:", output_dir, "\n")
cat("Files saved:\n")
cat("- GSE207128_pseudobulk_by_treatment.rds\n")
cat("- GSE207128_pseudobulk_celltype_by_treatment.rds\n")
cat("- GSE207128_pseudobulk_gene_mapping.csv\n")