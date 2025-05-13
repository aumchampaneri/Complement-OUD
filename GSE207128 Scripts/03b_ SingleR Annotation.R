library(Seurat)
library(SingleR)
library(celldex)
library(ggplot2)
library(pheatmap)
library(BiocParallel)
library(biomaRt)

# Set CRAN mirror for package installation
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install required packages
required_packages <- c("scran", "scuttle", "scater", "BiocNeighbors", "BiocParallel")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(pkg)
  }
}

# Set up paths
base_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE207128"
validation_dir <- file.path(base_dir, "harmony_results/annotation_validation")
dir.create(validation_dir, showWarnings = FALSE, recursive = TRUE)

# Load your annotated Seurat object
seurat_annotated <- readRDS(file.path(base_dir, "harmony_results/cell_annotation/seurat_annotated.rds"))

# Get expression matrix (log-normalized data)
expr_matrix <- GetAssayData(seurat_annotated, slot = "data")

# Better Ensembl to symbol conversion function
convert_ensembl_to_symbol <- function(matrix) {
  if(!grepl("^ENSMUS", rownames(matrix)[1])) {
    cat("Data already using gene symbols. Skipping conversion.\n")
    return(matrix)
  }

  cat("Converting Ensembl IDs to gene symbols...\n")

  # Set up biomaRt connection
  mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

  # Process in batches to avoid timeouts
  batch_size <- 1000
  total_genes <- nrow(matrix)
  all_mappings <- data.frame(ensembl_gene_id=character(), external_gene_name=character())

  for(i in seq(1, total_genes, by=batch_size)) {
    end_idx <- min(i + batch_size - 1, total_genes)
    cat("Mapping batch", i, "to", end_idx, "of", total_genes, "genes\n")

    batch_ids <- rownames(matrix)[i:end_idx]
    batch_mapping <- getBM(
      attributes = c("ensembl_gene_id", "external_gene_name"),
      filters = "ensembl_gene_id",
      values = batch_ids,
      mart = mart
    )

    all_mappings <- rbind(all_mappings, batch_mapping)
  }

  # Keep only rows with valid gene symbols
  all_mappings <- all_mappings[all_mappings$external_gene_name != "", ]
  # Deal with multiple mappings (take first one)
  all_mappings <- all_mappings[!duplicated(all_mappings$ensembl_gene_id), ]

  # Create a new matrix with gene symbols as rownames
  idx <- match(rownames(matrix), all_mappings$ensembl_gene_id)
  valid <- !is.na(idx)

  new_matrix <- matrix[valid, ]
  rownames(new_matrix) <- all_mappings$external_gene_name[idx[valid]]
  # Remove duplicate gene symbols (keep first occurrence)
  new_matrix <- new_matrix[!duplicated(rownames(new_matrix)), ]

  cat("Converted", nrow(new_matrix), "genes of", nrow(matrix), "total\n")
  return(new_matrix)
}

# Apply the conversion
expr_matrix <- convert_ensembl_to_symbol(expr_matrix)

# Print some diagnostics about the expression matrix
cat("Expression matrix dimensions:", dim(expr_matrix), "\n")
cat("Sample gene symbols:", head(rownames(expr_matrix)), "\n")

# Try multiple reference datasets for better matching
cat("\n=== Testing different reference datasets ===\n")

# 1. MouseRNAseqData (Brain-focused)
brain_ref <- MouseRNAseqData()
brain_common <- intersect(rownames(expr_matrix), rownames(brain_ref))
cat("MouseRNAseqData: Found", length(brain_common), "common genes\n")

# 2. ImmGen (Immune-focused)
immgen_ref <- ImmGenData()
immgen_common <- intersect(rownames(expr_matrix), rownames(immgen_ref))
cat("ImmGenData: Found", length(immgen_common), "common genes\n")

# Create standardization function to normalize cell type names
standardize_cell_type <- function(cell_type) {
  cell_type <- tolower(cell_type)

  # Brain cells
  if(grepl("neuron|excit|inhib", cell_type)) return("Neuron")
  if(grepl("astro", cell_type)) return("Astrocyte")
  if(grepl("oligo|opc", cell_type)) return("Oligodendrocyte")
  if(grepl("micro|glia", cell_type)) return("Microglia")
  if(grepl("endo", cell_type)) return("Endothelial")

  # Immune cells
  if(grepl("^b.?cell|b220|cd19", cell_type)) return("B cell")
  if(grepl("^t.?cell|cd4|cd8|lympho", cell_type)) return("T cell")
  if(grepl("macro|phag", cell_type)) return("Macrophage")
  if(grepl("mono", cell_type)) return("Monocyte")
  if(grepl("dc|dendri", cell_type)) return("Dendritic cell")

  # Return original if no match
  return(cell_type)
}

# Run SingleR on the brain reference dataset
cat("\n=== Running SingleR with MouseRNAseqData ===\n")
singler_brain <- SingleR(
  test = expr_matrix[brain_common, ],
  ref = brain_ref[brain_common, ],
  labels = brain_ref$label.main,
  BPPARAM = MulticoreParam(4)  # Use 4 cores for faster processing
)

# Run SingleR on ImmGen dataset
cat("\n=== Running SingleR with ImmGenData ===\n")
singler_immgen <- SingleR(
  test = expr_matrix[immgen_common, ],
  ref = immgen_ref[immgen_common, ],
  labels = immgen_ref$label.main,
  BPPARAM = MulticoreParam(4)
)

# Add SingleR labels to Seurat object
seurat_annotated$singler_brain_label <- singler_brain$labels
seurat_annotated$singler_immgen_label <- singler_immgen$labels
seurat_annotated$singler_label <- singler_brain$labels  # Use brain as default

# Add standardized labels
seurat_annotated$marker_standard <- sapply(seurat_annotated$predicted_celltype, standardize_cell_type)
seurat_annotated$singler_standard <- sapply(seurat_annotated$singler_label, standardize_cell_type)

# Compare annotations (with original labels)
raw_comparison <- table(
  "Marker-based" = seurat_annotated$predicted_celltype,
  "SingleR" = seurat_annotated$singler_label
)

# Compare annotations (with standardized labels)
standard_comparison <- table(
  "Marker-based" = seurat_annotated$marker_standard,
  "SingleR" = seurat_annotated$singler_standard
)

# Calculate agreement for both raw and standardized labels
raw_agreement <- sum(diag(raw_comparison)) / sum(raw_comparison) * 100
standard_agreement <- sum(diag(standard_comparison)) / sum(standard_comparison) * 100

cat("\n=== Annotation Agreement ===\n")
cat("Raw label agreement:", round(raw_agreement, 2), "%\n")
cat("Standardized label agreement:", round(standard_agreement, 2), "%\n")

# Save comparison results
write.csv(raw_comparison, file.path(validation_dir, "raw_annotation_comparison.csv"))
write.csv(standard_comparison, file.path(validation_dir, "standardized_annotation_comparison.csv"))

# Create heatmaps of the comparisons
pheatmap(log1p(raw_comparison),
         filename = file.path(validation_dir, "raw_annotation_heatmap.png"),
         main = "Raw Label Comparison",
         fontsize_row = 8,
         fontsize_col = 8)

pheatmap(log1p(standard_comparison),
         filename = file.path(validation_dir, "standardized_annotation_heatmap.png"),
         main = "Standardized Label Comparison",
         fontsize_row = 10,
         fontsize_col = 10)

# Plot UMAP with different label sets
p1 <- DimPlot(seurat_annotated, reduction = "umap", group.by = "predicted_celltype",
              label = TRUE, repel = TRUE) + ggtitle("Marker-based")
p2 <- DimPlot(seurat_annotated, reduction = "umap", group.by = "singler_label",
              label = TRUE, repel = TRUE) + ggtitle("SingleR (Brain)")
p3 <- DimPlot(seurat_annotated, reduction = "umap", group.by = "marker_standard",
              label = TRUE, repel = TRUE) + ggtitle("Standardized Marker")
p4 <- DimPlot(seurat_annotated, reduction = "umap", group.by = "singler_standard",
              label = TRUE, repel = TRUE) + ggtitle("Standardized SingleR")

# Save plots individually to avoid viewport issues
ggsave(file.path(validation_dir, "marker_annotation.png"), p1, width = 10, height = 8)
ggsave(file.path(validation_dir, "singler_annotation.png"), p2, width = 10, height = 8)
ggsave(file.path(validation_dir, "standardized_marker.png"), p3, width = 10, height = 8)
ggsave(file.path(validation_dir, "standardized_singler.png"), p4, width = 10, height = 8)

# Save validation results
saveRDS(list(brain = singler_brain, immgen = singler_immgen),
        file.path(validation_dir, "singler_results.rds"))
saveRDS(seurat_annotated, file.path(validation_dir, "seurat_dual_annotated.rds"))

cat("\nSingleR annotation complete! Results saved to", validation_dir, "\n")
cat("Use the standardized labels for better agreement between methods.\n")