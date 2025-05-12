library(nichenetr)
library(tidyverse)
library(Seurat)
library(Matrix)
library(ggplot2)

# ─── Load the data with optimization for large file ───────────────────
library(SeuratDisk)

# Define the cell types you're interested in
glia <- c('Microglia', 'Oligos', 'Oligos_Pre', 'Astrocytes')

# Path to h5seurat file
filepath_h5seurat <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/Super Series - GSE233279/GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5seurat"

cat("Loading h5seurat metadata to identify cells of interest...\n")
# First just get the metadata
h5_meta <- H5Seurat(filepath_h5seurat, mode = "r")
cell_meta <- h5_meta$meta.data$get()
h5_meta$close_all()

# Filter cells based on cell type
cat("Identifying cells belonging to glial cell types...\n")
cells_to_load <- rownames(cell_meta)[cell_meta$celltype3 %in% glia &
                                    cell_meta$Dx_OUD %in% c('OUD', 'None')]

# Now load only the cells you need
cat("Loading only the glial cells from the h5seurat file...\n")
seurat_obj <- LoadH5Seurat(filepath_h5seurat, cells = cells_to_load)
cat("Seurat object loaded successfully with only glial cells!\n")

# ─── Define sender and receiver cell types ───────────────────────
sender_celltypes <- c("Microglia")
receiver_celltypes <- c("Oligodendrocytes", "OPCs", "Astrocytes")

# Map from celltype3 to simplified cell types
seurat_obj$celltype <- plyr::mapvalues(
  seurat_obj$celltype3,
  from = c('Microglia', 'Oligos', 'Oligos_Pre', 'Astrocytes'),
  to = c('Microglia', 'Oligodendrocytes', 'OPCs', 'Astrocytes')
)

# ─── Preprocess expression data for NicheNet ───────────────────────
# Split by diagnosis
Idents(seurat_obj) <- "Dx_OUD"
oud_obj <- subset(seurat_obj, idents = "OUD")
ctl_obj <- subset(seurat_obj, idents = "CTL")

# Load prior knowledge networks from NicheNet
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
lr_network <- readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
weighted_networks <- readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))

# Get expressed genes in sender and receiver
expressed_genes_sender <- get_expressed_genes(oud_obj, sender_celltypes, "celltype")
expressed_genes_receiver <- get_expressed_genes(oud_obj, receiver_celltypes, "celltype")

# Define potential ligands
ligands <- lr_network %>%
  pull(from) %>%
  unique() %>%
  intersect(expressed_genes_sender)

# ─── Define DE genes between OUD and CTL in receivers ───────────────────
# Function to get DE genes
get_de_genes <- function(seurat_obj, ctl_obj, receiver_type) {
  Idents(seurat_obj) <- "celltype"
  cells_receiver <- WhichCells(seurat_obj, idents = receiver_type)
  receiver_obj <- seurat_obj[, cells_receiver]

  Idents(ctl_obj) <- "celltype"
  cells_receiver_ctl <- WhichCells(ctl_obj, idents = receiver_type)
  receiver_ctl_obj <- ctl_obj[, cells_receiver_ctl]

  # Merge objects and find DE genes
  combined <- merge(receiver_obj, receiver_ctl_obj)
  Idents(combined) <- "Dx_OUD"
  de_genes <- FindMarkers(combined, ident.1 = "OUD", ident.2 = "CTL",
                         min.pct = 0.1, logfc.threshold = 0.25)

  # Filter for significant genes
  de_genes <- de_genes %>%
    filter(p_val_adj < 0.05) %>%
    rownames_to_column("gene")

  return(de_genes$gene)
}

# Get DE genes for each receiver cell type
de_genes_oligo <- get_de_genes(seurat_obj, ctl_obj, "Oligodendrocytes")
de_genes_opc <- get_de_genes(seurat_obj, ctl_obj, "OPCs")
de_genes_astro <- get_de_genes(seurat_obj, ctl_obj, "Astrocytes")

# ─── Define inflammatory gene sets ───────────────────────────────
# Load MSigDB gene sets (you can download from MSigDB website)
# Example with inflammatory gene sets
inflammatory_go <- c("GO:0006954", "GO:0002526", "GO:0050729") # Inflammatory response gene sets
inflammatory_kegg <- c("hsa04668", "hsa04750") # TNF signaling, Inflammatory mediator regulation

# Get genes from these sets (simplified example, you would need to map GO IDs to genes)
inflammatory_genes <- c(de_genes_oligo, de_genes_opc, de_genes_astro) %>%
  intersect(expressed_genes_receiver)

# ─── Run NicheNet analysis ───────────────────────────────────────
# Perform ligand activity analysis
ligand_activities <- predict_ligand_activities(
  geneset = inflammatory_genes,
  background_expressed_genes = expressed_genes_receiver,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = ligands
)

# Get the best ligands explaining inflammatory changes
best_ligands <- ligand_activities %>%
  top_n(20, pearson) %>%
  arrange(desc(pearson)) %>%
  pull(test_ligand)

# Calculate target genes of top ligands
active_ligand_target_links <- best_ligands %>%
  lapply(function(lig) {
    get_weighted_ligand_target_links(lig, ligand_target_matrix, 0.33)
  }) %>%
  bind_rows()

# ─── Create Visualizations ───────────────────────────────────────
# Ligand activity plot
ligand_activity_plot <- ligand_activities %>%
  top_n(30, pearson) %>%
  mutate(test_ligand = fct_reorder(test_ligand, pearson)) %>%
  ggplot(aes(x = test_ligand, y = pearson)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Ligands", y = "Pearson correlation",
       title = "Ligand activity in OUD inflammatory response")

# Save results
output_dir <- "NicheNet_output"
dir.create(output_dir, showWarnings = FALSE)

# Save plots
ggsave(file.path(output_dir, "ligand_activity.pdf"), ligand_activity_plot, width = 10, height = 6)

# Save ligand activity data
write.csv(ligand_activities, file.path(output_dir, "ligand_activities.csv"), row.names = FALSE)
write.csv(active_ligand_target_links, file.path(output_dir, "ligand_target_links.csv"), row.names = FALSE)

# ─── Create circos plot of top ligand-target interactions ───────────────
# (This requires the circlize package)
library(circlize)

# Create circos plot for top ligands and their targets
top_links <- active_ligand_target_links %>%
  filter(ligand %in% best_ligands[1:5]) %>%
  top_n(50, weight)

# Save data for circos
write.csv(top_links, file.path(output_dir, "top_ligand_target_links.csv"), row.names = FALSE)

# Output completion message
cat("NicheNet analysis complete. Results saved to", output_dir, "\n")