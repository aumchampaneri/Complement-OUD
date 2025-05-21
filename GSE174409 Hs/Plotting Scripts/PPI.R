# ===== Protein-Protein Interaction (PPI) Network Analysis =====
# This script creates a STRING-based PPI network for complement genes,
# with options to highlight sex differences in expression

# ===== 1. Install and load required packages =====
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("STRINGdb", quietly = TRUE))
  BiocManager::install("STRINGdb")
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
if (!requireNamespace("RColorBrewer", quietly = TRUE))
  install.packages("RColorBrewer")
if (!requireNamespace("igraph", quietly = TRUE))
  install.packages("igraph")

library(STRINGdb)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(igraph)

# ===== 2. Load data =====
# Load DEG results by sex
cat("Loading DEG results by sex and region...\n")
deg_files <- list.files("GSE174409/NeuroinflammationResults/DEG_by_groups",
                        pattern = "DEG_.*\\.csv", full.names = TRUE)

# Create a list to store all results
all_deg_results <- list()
for (file in deg_files) {
  # Extract sex and region from filename
  filename <- basename(file)
  parts <- strsplit(gsub("DEG_|.csv", "", filename), "_")[[1]]
  sex <- parts[1]
  region <- parts[2]

  # Load the data
  results <- read.csv(file)
  all_deg_results[[paste(sex, region, sep="_")]] <- results
}

# Load complement gene list
cat("Loading complement gene list...\n")
comp_genes_file <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/Super Folder - GSE225158/GSE225158/KEGG outputs/kegg_complement_unique_genes.csv"
comp_genes <- read.csv(comp_genes_file)$gene
comp_genes <- toupper(comp_genes)

# Load gene mapping (ENSEMBL to Symbol)
gene_map <- read.csv('GSE174409/NeuroinflammationResults/ensembl_to_symbol_mapping.csv')
ensembl_to_symbol <- setNames(gene_map$Symbol, gene_map$EnsemblID)
symbol_to_ensembl <- setNames(gene_map$EnsemblID, toupper(gene_map$Symbol))

# ===== 3. Determine genes with sex differences =====
cat("Identifying complement genes with sex differences...\n")

# Method 1: Find complement genes that are significant in at least one sex-region combination
sig_comp_genes <- c()
for (key in names(all_deg_results)) {
  results <- all_deg_results[[key]]
  # Filter for significant complement genes
  sig_genes <- results %>%
    filter(adj.P.Val < 0.05 & is_complement == TRUE) %>%
    pull(gene_symbol)

  sig_comp_genes <- c(sig_comp_genes, sig_genes)
}
sig_comp_genes <- unique(sig_comp_genes)

# If we don't have enough significant genes, use all complement genes we can map
if (length(sig_comp_genes) < 5) {
  cat("Few significant complement genes found. Using full complement gene list.\n")
  # Map the complement genes we have symbols for
  mapped_comp_genes <- comp_genes[comp_genes %in% toupper(gene_map$Symbol)]
  if (length(mapped_comp_genes) > 0) {
    sig_comp_genes <- mapped_comp_genes
  } else {
    # Fallback to a standard list of complement genes
    sig_comp_genes <- c("C1QA", "C1QB", "C1QC", "C3", "C4A", "C4B", "C1R", "C1S",
                        "CFB", "CFH", "CFI", "C5AR1", "CR1", "VSIG4")
  }
}

cat("Selected", length(sig_comp_genes), "complement genes for PPI analysis.\n")
cat("Genes:", paste(head(sig_comp_genes, 10), collapse=", "),
    ifelse(length(sig_comp_genes) > 10, "...", ""), "\n")

# ===== 4. Initialize STRING database =====
cat("Initializing STRINGdb...\n")
string_db <- STRINGdb$new(
  version = "11.5",
  species = 9606,  # human
  score_threshold = 400,  # medium confidence
  input_directory = ""
)

# ===== 5. Map genes to STRING IDs =====
# Create data frame for mapping
df_genes <- data.frame(gene = sig_comp_genes)

# Map gene symbols to STRING IDs
cat("Mapping genes to STRING IDs...\n")
mapped_genes <- string_db$map(df_genes, "gene", removeUnmappedRows = TRUE)
cat("Successfully mapped", nrow(mapped_genes), "out of", length(sig_comp_genes), "genes\n")

# ===== 6. Get logFC values for each gene by sex and region =====
# Create a data frame to store logFC values
logfc_data <- data.frame(
  gene = mapped_genes$gene,
  Female_NAC = NA,
  Female_DLPFC = NA,
  Male_NAC = NA,
  Male_DLPFC = NA
)

# Fill in logFC values
for (key in names(all_deg_results)) {
  results <- all_deg_results[[key]]
  # Extract logFC for each gene
  for (i in 1:nrow(logfc_data)) {
    gene <- logfc_data$gene[i]
    idx <- which(results$gene_symbol == gene)
    if (length(idx) > 0) {
      logfc_data[[key]][i] <- results$logFC[idx[1]]
    }
  }
}

# Merge with mapped genes
annotated_genes <- merge(mapped_genes, logfc_data, by = "gene")

# ===== 7. Plot basic PPI network =====
cat("Generating basic PPI network...\n")
# Create output directory if it doesn't exist
dir.create("GSE174409/figures/ppi_networks", showWarnings = FALSE, recursive = TRUE)

pdf("GSE174409/figures/ppi_networks/complement_basic_ppi.pdf", width = 10, height = 8)
string_db$plot_network(annotated_genes$STRING_id)
dev.off()

# ===== 8. Generate IMAGES with different coloring schemes =====
cat("Generating colored PPI networks for each condition...\n")

# Get the interaction network data
interactions <- string_db$get_interactions(annotated_genes$STRING_id)

# Create an igraph network object
g <- graph_from_data_frame(interactions[, c("from", "to")], directed = FALSE)

# Add gene names as node labels
idx <- match(V(g)$name, annotated_genes$STRING_id)
V(g)$label <- annotated_genes$gene[idx]

# Function to create a colored network using igraph
create_colored_network <- function(condition, color_column) {
  # Map colors to logFC values
  idx <- match(V(g)$name, annotated_genes$STRING_id)
  logfc_values <- annotated_genes[[color_column]][idx]

  # Create a color gradient
  max_abs_logfc <- max(abs(logfc_values), na.rm = TRUE)
  color_gradient <- colorRampPalette(c("blue", "white", "red"))(100)

  # Scale logFC values to color indices (1-100)
  scaled_values <- ((logfc_values / max_abs_logfc) + 1) / 2 * 99 + 1
  scaled_values[is.na(scaled_values)] <- 50  # neutral color for NA

  # Assign colors to nodes
  V(g)$color <- color_gradient[round(scaled_values)]

  # Save to PDF
  pdf(paste0("GSE174409/figures/ppi_networks/complement_ppi_", condition, ".pdf"),
      width = 10, height = 10)

  # Plot the network
  plot(g,
       layout = layout_with_fr(g),
       vertex.size = 10,
       vertex.label.cex = 0.8,
       vertex.label.dist = 0.5,
       edge.width = 1,
       main = paste("PPI Network -", condition))

  # Add a legend for the color scale
  legend_colors <- color_gradient[c(1, 25, 50, 75, 100)]
  legend_labels <- c(
    paste0("Down: ", round(-max_abs_logfc, 2)),
    paste0("-", round(max_abs_logfc/2, 2)),
    "0",
    paste0("+", round(max_abs_logfc/2, 2)),
    paste0("Up: ", round(max_abs_logfc, 2))
  )

  legend("bottomright", legend = legend_labels, fill = legend_colors,
         title = "logFC", cex = 0.8, bty = "n")

  dev.off()
}

# Create networks for each condition
create_colored_network("Female_NAC", "Female_NAC")
create_colored_network("Female_DLPFC", "Female_DLPFC")
create_colored_network("Male_NAC", "Male_NAC")
create_colored_network("Male_DLPFC", "Male_DLPFC")

# ===== 9. Export data for Cytoscape =====
cat("Exporting data for Cytoscape visualization...\n")

# Write out node and edge files for Cytoscape
write.csv(annotated_genes, "GSE174409/figures/ppi_networks/complement_nodes.csv", row.names = FALSE)
write.csv(interactions, "GSE174409/figures/ppi_networks/complement_edges.csv", row.names = FALSE)

# ===== 10. Generate a README file with instructions for Cytoscape =====
cat("Creating README with Cytoscape instructions...\n")

readme_text <- "# Complement PPI Network Visualization Instructions

## Files
- complement_nodes.csv: Node information including gene symbols and logFC values
- complement_edges.csv: Edge information from STRING database
- complement_basic_ppi.pdf: Basic PPI network visualization
- colored network PDFs: Networks colored by logFC values for each condition

## Instructions for Cytoscape Visualization
1. Open Cytoscape
2. Install the STRING app if not already installed (Apps -> App Manager)
3. Import the node and edge files (File -> Import -> Network from Files)
4. Map the logFC values to node colors (Style tab -> Fill Color -> Map Column)
5. Set node size based on degree (Style tab -> Size -> Map Column)
6. Adjust layout as desired (Layout -> yFiles Layouts -> Organic)
7. Export as high-resolution PNG/PDF (File -> Export -> Network to Image)

## Additional Analyses
For hub gene analysis, use the NetworkAnalyzer tool (Tools -> NetworkAnalyzer -> Analyze Network)
"

writeLines(readme_text, "GSE174409/figures/ppi_networks/README.md")

cat("\nPPI analysis complete. Results saved in 'GSE174409/figures/ppi_networks/'\n")