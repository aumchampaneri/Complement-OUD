library(Seurat)
library(SingleR)
library(celldex)
library(ggplot2)
library(pheatmap)
library(biomaRt)  # Add this for ID conversion

# Add at the top of your script after the library() calls
if (!requireNamespace("scran", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("scran")
}

# Also make sure other potential dependencies are installed
required_packages <- c("scuttle", "scater", "BiocNeighbors")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# Load your annotated Seurat object
seurat_annotated <- readRDS("/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE207128/harmony_results/cell_annotation/seurat_annotated.rds")

# Create output directory for validation
validation_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE207128/harmony_results/annotation_validation"
dir.create(validation_dir, showWarnings = FALSE, recursive = TRUE)

# Get expression matrix (log-normalized data)
expr_matrix <- GetAssayData(seurat_annotated, layer = "data")  # Fixed deprecated 'slot' parameter

# Check if we need to convert Ensembl IDs to gene symbols
if(grepl("^ENSMUS", rownames(expr_matrix)[1])) {
  cat("Converting Ensembl IDs to gene symbols for SingleR compatibility...\n")

  # Set up biomaRt connection
  mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

  # Get mappings for all genes in the dataset
  ensembl_ids <- rownames(expr_matrix)

  # Get mappings in batches to avoid timeout
  batch_size <- 1000
  all_mappings <- data.frame(ensembl_gene_id=character(), external_gene_name=character())

  for(i in seq(1, length(ensembl_ids), by=batch_size)) {
    end_idx <- min(i + batch_size - 1, length(ensembl_ids))
    cat("Mapping batch", i, "to", end_idx, "of", length(ensembl_ids), "genes\n")

    batch_ids <- ensembl_ids[i:end_idx]
    batch_mapping <- getBM(
      attributes = c("ensembl_gene_id", "external_gene_name"),
      filters = "ensembl_gene_id",
      values = batch_ids,
      mart = mart
    )

    all_mappings <- rbind(all_mappings, batch_mapping)
  }

  # Remove duplicate gene symbols (keep first occurrence)
  all_mappings <- all_mappings[!duplicated(all_mappings$external_gene_name),]

  # Filter out empty gene symbols
  all_mappings <- all_mappings[all_mappings$external_gene_name != "",]

  cat("Found mappings for", nrow(all_mappings), "of", length(ensembl_ids), "genes\n")

  # Create a new matrix with gene symbols as rownames
  matched_indices <- match(rownames(expr_matrix), all_mappings$ensembl_gene_id)
  genes_with_symbols <- !is.na(matched_indices)

  # Create new matrix with only genes that have symbol mappings
  converted_matrix <- expr_matrix[genes_with_symbols,]
  rownames(converted_matrix) <- all_mappings$external_gene_name[matched_indices[genes_with_symbols]]

  # Remove duplicated gene symbols (keep first occurrence)
  converted_matrix <- converted_matrix[!duplicated(rownames(converted_matrix)),]

  expr_matrix <- converted_matrix
  cat("Expression matrix converted to use gene symbols. Dimensions:", dim(expr_matrix), "\n")
}

# 1. Use Mouse Brain Atlas as reference
mouse_ref <- celldex::MouseRNAseqData()

# Check for common genes
common_genes <- intersect(rownames(expr_matrix), rownames(mouse_ref))
cat("Number of common genes between dataset and reference:", length(common_genes), "\n")

if(length(common_genes) < 100) {
  stop("Too few common genes for reliable annotation. Check gene ID formats.")
}

# Run SingleR with common genes
singler_results <- SingleR(
  test = expr_matrix,
  ref = mouse_ref,
  labels = mouse_ref$label.main,
  de.method = "wilcox"
)

# Add SingleR labels to Seurat object
seurat_annotated$singler_label <- singler_results$labels

# Compare annotations
comparison_table <- table(
  "Marker-based" = seurat_annotated$predicted_celltype,
  "SingleR" = seurat_annotated$singler_label
)

# Save comparison results
write.csv(comparison_table, file.path(validation_dir, "annotation_comparison.csv"))

# Create a heatmap of the comparison
pheatmap(log1p(comparison_table),
         filename = file.path(validation_dir, "annotation_comparison_heatmap.png"),
         main = "Marker-based vs Reference-based Annotation")

# Plot UMAP with SingleR labels
p1 <- DimPlot(seurat_annotated, reduction = "umap", group.by = "predicted_celltype",
              label = TRUE, repel = TRUE) + ggtitle("Marker-based")
p2 <- DimPlot(seurat_annotated, reduction = "umap", group.by = "singler_label",
              label = TRUE, repel = TRUE) + ggtitle("SingleR Reference-based")
combined <- p1 + p2
ggsave(file.path(validation_dir, "annotation_method_comparison.png"), combined, width = 16, height = 8)

# Calculate agreement percentage
agreement <- sum(diag(comparison_table)) / sum(comparison_table) * 100
cat("Overall agreement between methods:", round(agreement, 2), "%\n")

# Save validation results
saveRDS(singler_results, file.path(validation_dir, "singler_results.rds"))
saveRDS(seurat_annotated, file.path(validation_dir, "seurat_dual_annotated.rds"))