# Install necessary packages (uncomment if needed)
# install.packages("Seurat")
# install.packages("remotes")
# remotes::install_github("mojaveazure/seurat-disk")
# remotes::install_github("sqjin/CellChat")
# BiocManager::install("limma")

# Load libraries
library(Seurat)
library(SeuratDisk)
library(CellChat)
library(patchwork)
library(ggplot2)

# Create output directory
output_dir <- "Glial Analysis/CellChat output"
dir.create(output_dir, showWarnings = FALSE)

# Directly load the h5seurat file
h5seurat_file <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/Super Series - GSE233279/GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat"
seurat_obj <- SeuratDisk::LoadH5Seurat(h5seurat_file)

# Debugging: Print available metadata columns
print("Available metadata columns:")
print(colnames(seurat_obj@meta.data))

# Debugging: Check for disease and sex columns
print("First few rows of metadata:")
print(head(seurat_obj@meta.data))

# Debugging: Check celltype3 distribution
print("Celltype3 distribution:")
print(table(seurat_obj$celltype3))

# Write detailed metadata to CSV for reference
write.csv(seurat_obj@meta.data, file.path(output_dir, "seurat_metadata.csv"))

# Use celltype3 as identifier
Idents(seurat_obj) <- seurat_obj$celltype3

# Extract data for CellChat
data.input <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
labels <- Idents(seurat_obj)
meta <- data.frame(group = labels, row.names = names(labels))

# Create CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

# Set database (choose appropriate one)
cellchat@DB <- CellChatDB.human # or CellChatDB.mouse

# Preprocessing
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Compute communication probability
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Compute and analyze network
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Visualizations with updated paths
# Network overview - Creates chord diagrams
pdf(file.path(output_dir, "cellchat_networks.pdf"), width=12, height=8)
par(mfrow = c(1,2))
netVisual_circle(cellchat@net$count, vertex.weight = table(cellchat@idents),
                weight.scale = TRUE, label.edge= FALSE,
                title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = table(cellchat@idents),
                weight.scale = TRUE, label.edge= FALSE,
                title.name = "Interaction strength")
dev.off()

# Interaction heatmap
pdf(file.path(output_dir, "cellchat_heatmap.pdf"), width=10, height=8)
netVisual_heatmap(cellchat, measure = "count")
dev.off()

# Top pathways
top_pathways <- showSignalingPathway(cellchat, sort.by = "pval", return.data = TRUE)
write.csv(top_pathways, file.path(output_dir, "top_signaling_pathways.csv"))

# Plot specific pathways
pdf(file.path(output_dir, "cellchat_top_pathways.pdf"), width=10, height=8)
pathways.show <- top_pathways$pathway[1:min(5, length(top_pathways$pathway))]
for (i in 1:length(pathways.show)) {
  netVisual_aggregate(cellchat, signaling = pathways.show[i], layout = "circle")
}
dev.off()

# Additional chord diagrams for specific pathways
pdf(file.path(output_dir, "cellchat_pathway_circles.pdf"), width=12, height=10)
for (i in 1:length(pathways.show)) {
  netVisual_aggregate(cellchat, signaling = pathways.show[i], layout = "circle")
  netVisual_chord_cell(cellchat, signaling = pathways.show[i], title.name = paste0(pathways.show[i], " signaling network"))
}
dev.off()

# Advanced analysis
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
pdf(file.path(output_dir, "cellchat_roles.pdf"), width=10, height=8)
netAnalysis_signalingRole_scatter(cellchat)
dev.off()

# Save results
saveRDS(cellchat, file = file.path(output_dir, "cellchat_results.rds"))