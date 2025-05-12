# Load libraries
library(Seurat)
library(SeuratDisk)
library(CellChat)
library(patchwork)
library(ggplot2)
library(future)

# Memory and parallel processing optimization
options(future.globals.maxSize = 4000 * 1024^2)  # Allow 4GB for parallelization
plan("multicore", workers = 4)  # Use 4 cores on M1 Mac
# Removed the problematic line: CellChat::setParallel(parallel = TRUE)

# Create output directory
output_dir <- "Glial Analysis/CellChat output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

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

# Memory optimization: Focus on glial cells only
glial_celltypes <- c("Astro", "Micro", "Oligo", "OPC")
seurat_glial <- subset(seurat_obj, celltype3 %in% glial_celltypes)
rm(seurat_obj); gc()  # Free memory

# Extract data for CellChat with gene filtering
data.input <- GetAssayData(seurat_glial, assay = "RNA", slot = "data")
# Keep only genes expressed in at least 10% of cells
min_cells <- 0.1 * min(table(Idents(seurat_glial)))
expressed_genes <- which(rowSums(data.input > 0) >= min_cells)
data.input <- data.input[expressed_genes, ]
rm(expressed_genes); gc()

# Create metadata
labels <- Idents(seurat_glial)
meta <- data.frame(group = labels, row.names = names(labels))
rm(seurat_glial, labels); gc()  # Free memory

# Create CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
rm(data.input, meta); gc()  # Free memory

# Focus on inflammatory pathways
CellChatDB <- CellChatDB.human
inflammatory_keywords <- c("inflammatory", "cytokine", "interferon", "immune",
                          "complement", "IL", "TNF", "NFkB", "toll-like",
                          "leukocyte", "chemokine")
pattern <- paste(inflammatory_keywords, collapse = "|")
all_pathways <- unique(CellChatDB$interaction$pathway_name)
inflammatory_pathways <- all_pathways[grepl(pattern, all_pathways, ignore.case = TRUE)]
print(paste("Found", length(inflammatory_pathways), "inflammatory pathways"))

# Subset database to inflammatory pathways
CellChatDB.inflammatory <- subsetDB(CellChatDB, search = inflammatory_pathways, key = "pathway_name")
rm(CellChatDB, inflammatory_keywords, pattern, inflammatory_pathways, all_pathways); gc()

# Set database
cellchat@DB <- CellChatDB.inflammatory
rm(CellChatDB.inflammatory); gc()

# Preprocessing with intermediate saves
cellchat <- subsetData(cellchat)
saveRDS(cellchat, file.path(output_dir, "cellchat_step1.rds")); gc()

cellchat <- identifyOverExpressedGenes(cellchat)
saveRDS(cellchat, file.path(output_dir, "cellchat_step2.rds")); gc()

cellchat <- identifyOverExpressedInteractions(cellchat)
saveRDS(cellchat, file.path(output_dir, "cellchat_step3.rds")); gc()

# Compute communication probability with stricter filtering and parallelization
cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE, future.plan = "multicore", future.globals.maxSize = 4000 * 1024^2)
saveRDS(cellchat, file.path(output_dir, "cellchat_step4.rds")); gc()

cellchat <- filterCommunication(cellchat, min.cells = 10, min.prob = 0.1)
saveRDS(cellchat, file.path(output_dir, "cellchat_step5.rds")); gc()

# Compute and analyze network with parallelization
cellchat <- computeCommunProbPathway(cellchat, future.plan = "multicore", future.globals.maxSize = 4000 * 1024^2)
saveRDS(cellchat, file.path(output_dir, "cellchat_step6.rds")); gc()

cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, file.path(output_dir, "cellchat_complete.rds")); gc()

# Visualizations with memory cleanup between steps
pdf(file.path(output_dir, "cellchat_networks.pdf"), width=12, height=8)
par(mfrow = c(1,2))
netVisual_circle(cellchat@net$count, vertex.weight = table(cellchat@idents),
                weight.scale = TRUE, label.edge= FALSE,
                title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = table(cellchat@idents),
                weight.scale = TRUE, label.edge= FALSE,
                title.name = "Interaction strength")
dev.off(); gc()

# Interaction heatmap
pdf(file.path(output_dir, "cellchat_heatmap.pdf"), width=10, height=8)
netVisual_heatmap(cellchat, measure = "count")
dev.off(); gc()

# Top pathways
top_pathways <- showSignalingPathway(cellchat, sort.by = "pval", return.data = TRUE)
write.csv(top_pathways, file.path(output_dir, "top_signaling_pathways.csv"))

# Plot specific pathways - limit to top 3 to save memory
pdf(file.path(output_dir, "cellchat_top_pathways.pdf"), width=10, height=8)
pathways.show <- top_pathways$pathway[1:min(3, length(top_pathways$pathway))]
for (i in 1:length(pathways.show)) {
  netVisual_aggregate(cellchat, signaling = pathways.show[i], layout = "circle")
  gc()
}
dev.off(); gc()

# Chord diagrams for specific pathways
pdf(file.path(output_dir, "cellchat_pathway_circles.pdf"), width=12, height=10)
for (i in 1:length(pathways.show)) {
  netVisual_aggregate(cellchat, signaling = pathways.show[i], layout = "circle")
  netVisual_chord_cell(cellchat, signaling = pathways.show[i],
                      title.name = paste0(pathways.show[i], " signaling network"))
  gc()
}
dev.off(); gc()

# Advanced analysis
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
pdf(file.path(output_dir, "cellchat_roles.pdf"), width=10, height=8)
netAnalysis_signalingRole_scatter(cellchat)
dev.off(); gc()

# Save final results
saveRDS(cellchat, file = file.path(output_dir, "cellchat_results.rds"))