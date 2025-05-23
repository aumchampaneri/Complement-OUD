# Load required libraries
library(STRINGdb)
library(ggplot2)
library(VennDiagram)
library(plotly)
library(htmlwidgets)
library(dplyr)
library(biomaRt)

# Set timeout to 300 seconds
options(timeout = 300)

# Define input and output folder paths
input_folder <- '/Users/aumchampaneri/PycharmProjects/Complement-OUD/Super Folder - GSE225158/GSE225158/DESeq2 outputs'
output_folder <- '/Users/aumchampaneri/PycharmProjects/Complement-OUD/Super Folder - GSE225158/GSE225158/STRING outputs'

# Step 1: Load DESeq2 output CSVs
deg1 <- read.csv(file.path(input_folder, "deseq2_results_M_OUD_vs_M_None.csv"))
deg2 <- read.csv(file.path(input_folder, "deseq2_results_F_OUD_vs_F_None.csv"))

# Step 2: Filter for significant genes
sig1 <- subset(deg1, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
sig2 <- subset(deg2, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)

# Extract gene names
genes1 <- data.frame(gene = sig1$gene)
genes2 <- data.frame(gene = sig2$gene)

# Step 3: Initialize STRINGdb
string_db <- STRINGdb$new(
  version = "11.5",
  species = 9606,  # Human
  score_threshold = 700,
  input_directory = ""
)

# Step 4: Map gene names to STRING IDs
mapped1 <- string_db$map(genes1, "gene", removeUnmappedRows = FALSE)
mapped2 <- string_db$map(genes2, "gene", removeUnmappedRows = FALSE)

# Save unmapped genes for debugging
unmapped_genes1 <- mapped1[is.na(mapped1$STRING_id), ]
unmapped_genes2 <- mapped2[is.na(mapped2$STRING_id), ]
write.csv(unmapped_genes1, file.path(output_folder, "unmapped_genes_condition1.csv"), row.names = FALSE)
write.csv(unmapped_genes2, file.path(output_folder, "unmapped_genes_condition2.csv"), row.names = FALSE)

# Remove unmapped rows
mapped1 <- mapped1[!is.na(mapped1$STRING_id), ]
mapped2 <- mapped2[!is.na(mapped2$STRING_id), ]

# Save mapped genes for debugging
write.csv(mapped1, file.path(output_folder, "mapped_genes_condition1.csv"), row.names = FALSE)
write.csv(mapped2, file.path(output_folder, "mapped_genes_condition2.csv"), row.names = FALSE)

# Step 5: Get enrichment results
enrichment1 <- string_db$get_enrichment(mapped1$STRING_id)
enrichment2 <- string_db$get_enrichment(mapped2$STRING_id)

# Save enrichment results
write.csv(enrichment1, file.path(output_folder, "enrichment1.csv"), row.names = FALSE)
write.csv(enrichment2, file.path(output_folder, "enrichment2.csv"), row.names = FALSE)

# Step 6: Overlay log2FoldChange on STRING network
if (!"preferred_name" %in% colnames(mapped1)) {
  colnames(mapped1)[colnames(mapped1) == "gene"] <- "preferred_name"
}
if (!"preferred_name" %in% colnames(mapped2)) {
  colnames(mapped2)[colnames(mapped2) == "gene"] <- "preferred_name"
}

mapped1_fc <- merge(mapped1, sig1[, c("gene", "log2FoldChange")], by.x = "preferred_name", by.y = "gene")
mapped2_fc <- merge(mapped2, sig2[, c("gene", "log2FoldChange")], by.x = "preferred_name", by.y = "gene")

# Plot STRING network with fold change coloring
string_db$plot_network(mapped1_fc$STRING_id, values = mapped1_fc$log2FoldChange)
string_db$plot_network(mapped2_fc$STRING_id, values = mapped2_fc$log2FoldChange)

# Step 7: Save PPI networks
ppi1 <- string_db$get_interactions(mapped1$STRING_id)
ppi2 <- string_db$get_interactions(mapped2$STRING_id)

write.csv(ppi1, file.path(output_folder, "ppi_network_condition1.csv"), row.names = FALSE)
write.csv(ppi2, file.path(output_folder, "ppi_network_condition2.csv"), row.names = FALSE)

# Step 8: Export Interaction Tables for Cytoscape
write.table(ppi1, file.path(output_folder, "ppi_condition1_cytoscape.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ppi2, file.path(output_folder, "ppi_condition2_cytoscape.txt"), sep = "\t", quote = FALSE, row.names = FALSE)