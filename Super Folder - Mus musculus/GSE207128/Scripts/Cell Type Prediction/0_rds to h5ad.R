# Load required libraries
library(Seurat)
library(SeuratDisk)

# Set paths
input_path <- "/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus/GSE207128/Outputs/Processed Data/processed_seurat.rds"
output_dir <- "/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus/GSE207128/Outputs/Processed Data/"

# Load the Seurat object
seurat_obj <- readRDS(input_path)

# Print basic information about the object
cat("Seurat object loaded successfully!\n")
cat("Number of cells:", ncol(seurat_obj), "\n")
cat("Number of genes:", nrow(seurat_obj), "\n")
cat("Assays available:", names(seurat_obj@assays), "\n")

# Check if this is a v5 object and join layers if needed
if (inherits(seurat_obj[["RNA"]], "Assay5")) {
  cat("Detected Seurat v5 object with multiple layers. Joining layers...\n")
  seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
}

# Clean the object by keeping only essential data
# Create a minimal Seurat object with just the RNA counts
cleaned_obj <- CreateSeuratObject(
  counts = LayerData(seurat_obj, assay = "RNA", layer = "counts"),
  meta.data = seurat_obj@meta.data
)

# Add normalized data if available
if ("data" %in% Layers(seurat_obj[["RNA"]])) {
  cleaned_obj <- SetAssayData(cleaned_obj, layer = "data", 
                             new.data = LayerData(seurat_obj, assay = "RNA", layer = "data"))
}

# Try to add embeddings if they exist
if ("pca" %in% names(seurat_obj@reductions)) {
  cleaned_obj[["pca"]] <- seurat_obj[["pca"]]
}
if ("umap" %in% names(seurat_obj@reductions)) {
  cleaned_obj[["umap"]] <- seurat_obj[["umap"]]
}

# Convert to h5Seurat format
h5seurat_path <- file.path(output_dir, "processed_seurat_clean.h5Seurat")
SaveH5Seurat(cleaned_obj, filename = h5seurat_path, overwrite = TRUE)

# Convert h5Seurat to h5ad format
h5ad_path <- file.path(output_dir, "processed_seurat_clean.h5ad")
Convert(h5seurat_path, dest = "h5ad", overwrite = TRUE)

cat("Conversion completed!\n")
cat("H5AD file saved at:", h5ad_path, "\n")

# Clean up intermediate file
file.remove(h5seurat_path)
cat("Intermediate h5Seurat file removed.\n")