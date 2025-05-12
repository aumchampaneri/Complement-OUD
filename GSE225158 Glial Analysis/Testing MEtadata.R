library(SeuratDisk)

# Path to h5seurat file
filepath_h5seurat <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/Super Series - GSE233279/GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5seurat"

# Extract metadata and save as RDS for quicker access later
cat("Extracting metadata from h5seurat file...\n")
seurat_obj_preview <- LoadH5Seurat(filepath_h5seurat, load.counts = FALSE, load.assays = FALSE)
metadata <- seurat_obj_preview@meta.data
saveRDS(metadata, file = "/Users/aumchampaneri/PycharmProjects/Complement-OUD/Super Series - GSE233279/metadata.rds")
cat("Metadata successfully extracted and saved!\n")