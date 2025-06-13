# Export Seurat object to H5AD format for scVI/scANVI

# Load required packages
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("sceasy", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("cellgeni/sceasy")
}

library(Seurat)
library(sceasy)

# Setup Python environment
conda_envs <- reticulate::conda_list()
if ("scvi_env" %in% conda_envs$name) {
  reticulate::use_condaenv("scvi_env", required = TRUE)
} else {
  reticulate::py_install(c("anndata", "h5py", "numpy"))
}

# File paths
input_path <- "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study/data/raw/HuMicA/HuMicA.rds"
output_path <- "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study/data/raw/HuMicA/HuMicA_for_scVI.h5ad"

# Load and convert
cat("Loading Seurat object...\n")
seurat_obj <- readRDS(input_path)
print(seurat_obj)

cat("\nMetadata columns (for batch_key/labels_key):\n")
print(colnames(seurat_obj@meta.data))

cat("\nConverting to H5AD...\n")
if (file.exists(output_path)) file.remove(output_path)

sceasy::convertFormat(
  obj = seurat_obj,
  from = "seurat",
  to = "anndata",
  outFile = output_path,
  assay = "RNA",
  main_layer = "counts",
  drop_single_values = FALSE
)

if (file.exists(output_path)) {
  cat("\nSuccess! H5AD file created at:\n", output_path, "\n")
  cat("\nNext steps in Python:\n")
  cat("import scanpy as sc\n")
  cat("import scvi\n")
  cat("adata = sc.read_h5ad('HuMicA_for_scVI.h5ad')\n")
  cat("scvi.model.SCVI.setup_anndata(adata, batch_key='Sample_ID')\n")
  cat("model = scvi.model.SCVI(adata)\n")
  cat("model.train()\n")
} else {
  cat("\nError: Conversion failed.\n")
}
