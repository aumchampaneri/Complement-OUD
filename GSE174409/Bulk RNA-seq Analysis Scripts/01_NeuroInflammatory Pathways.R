# 01_NeuroInflammatory Pathways.R
# Analysis of neuroinflammation-related genes in the OUD RNA-seq dataset using KEGG pathways

# ----------------------
# 1. SETUP AND DATA LOADING
# ----------------------

# Load required libraries
suppressMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(pheatmap)
  library(dplyr)
  library(RColorBrewer)
  library(biomaRt)
  library(reshape2)
  library(tidyr)
  library(KEGGREST)
  library(plyr)
  library(clusterProfiler)
  library(org.Hs.eg.db) # For GO enrichment
})

# Set directories
base_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE174409"
qc_dir <- file.path(base_dir, "QC")
output_dir <- file.path(base_dir, "NeuroinflammationResults")
dir.create(output_dir, showWarnings = FALSE)

# Load preprocessed data
dge_filtered <- readRDS(file.path(qc_dir, "dge_filtered_normalized.rds"))
metadata_df <- readRDS(file.path(qc_dir, "metadata.rds"))
logcpm_filtered_norm <- readRDS(file.path(qc_dir, "logcpm_filtered_normalized.rds"))

# ----------------------
# 2. MAP ENSEMBL IDS TO GENE SYMBOLS
# ----------------------

# Get Ensembl IDs from data
ensembl_ids <- rownames(logcpm_filtered_norm)
cat("First few gene IDs in dataset:", head(ensembl_ids), "\n")

# Set up biomaRt to convert Ensembl IDs to gene symbols
cat("Converting Ensembl IDs to gene symbols using biomaRt...\n")
ensembl <- tryCatch({
  useMart("ensembl", dataset = "hsapiens_gene_ensembl")
}, error = function(e) {
  cat("Error connecting to Ensembl: ", e$message, "\n")
  cat("Using archived version from April 2023...\n")
  useMart("ENSEMBL_MART_ENSEMBL", host="https://apr2023.archive.ensembl.org",
          dataset="hsapiens_gene_ensembl")
})

mapping_result <- tryCatch({
  getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id", "description"),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = ensembl
  )
}, error = function(e) {
  cat("Error retrieving data from biomaRt: ", e$message, "\n")
  data.frame(
    ensembl_gene_id = character(0),
    external_gene_name = character(0),
    entrezgene_id = character(0),
    description = character(0),
    stringsAsFactors = FALSE
  )
})

# Create mapping dictionary
gene_symbols <- setNames(
  mapping_result$external_gene_name,
  mapping_result$ensembl_gene_id
)

# Count how many IDs were mapped
mapped_count <- sum(ensembl_ids %in% names(gene_symbols))
cat("Successfully mapped", mapped_count, "out of", length(ensembl_ids), "Ensembl IDs to gene symbols\n")

# Save the mapping for reference
mapping_df <- data.frame(
  EnsemblID = ensembl_ids,
  Symbol = gene_symbols[ensembl_ids],
  stringsAsFactors = FALSE
)
write.csv(mapping_df, file.path(output_dir, "ensembl_to_symbol_mapping.csv"), row.names = FALSE)

# ----------------------
# 3. IDENTIFY NEUROINFLAMMATION GENES USING KEGG PATHWAYS
# ----------------------

# Define relevant KEGG pathways for neuroinflammation
neuro_pathways <- c(
  "hsa04060", "hsa04061", "hsa04062", "hsa04620", "hsa04621", "hsa04622",
  "hsa04623", "hsa04625", "hsa04650", "hsa04657", "hsa04658", "hsa04659",
  "hsa04660", "hsa04662", "hsa04668", "hsa04672", "hsa04610", "hsa04640", "hsa04380"
)

# Retrieve genes from KEGG pathways
cat("Retrieving genes from KEGG pathways...\n")
neuro_genes_kegg <- list()

for (pathway in neuro_pathways) {
  tryCatch({
    pathway_info <- keggGet(pathway)
    if(length(pathway_info) > 0 && "GENE" %in% names(pathway_info[[1]])) {
      genes <- pathway_info[[1]]$GENE
      gene_ids <- genes[seq(1, length(genes), 2)]
      gene_symbols <- sapply(strsplit(genes[seq(2, length(genes), 2)], ";"), `[`, 1)
      gene_symbols <- trimws(gene_symbols)

      cat("Retrieved", length(gene_ids), "genes from pathway", pathway, "\n")

      neuro_genes_kegg[[pathway]] <- data.frame(
        EntrezID = gene_ids,
        Symbol = gene_symbols,
        Pathway = pathway,
        stringsAsFactors = FALSE
      )
    } else {
      cat("No genes found in pathway", pathway, "\n")
    }
  }, error = function(e) {
    cat("Error processing pathway", pathway, ":", e$message, "\n")
  })
}

# Combine all pathway genes
all_pathway_genes <- rbind.fill(neuro_genes_kegg)
all_pathway_genes <- all_pathway_genes[!duplicated(all_pathway_genes$Symbol), ]

# Get Entrez to Ensembl mapping
entrez_to_ensembl <- tryCatch({
  getBM(
    attributes = c("entrezgene_id", "ensembl_gene_id", "external_gene_name"),
    filters = "entrezgene_id",
    values = all_pathway_genes$EntrezID,
    mart = ensembl
  )
}, error = function(e) {
  cat("Error retrieving Entrez to Ensembl mapping: ", e$message, "\n")
  getBM(
    attributes = c("external_gene_name", "ensembl_gene_id"),
    filters = "external_gene_name",
    values = all_pathway_genes$Symbol,
    mart = ensembl
  )
})

# Find neuroinflammation genes in our dataset
neuro_ensembl <- entrez_to_ensembl$ensembl_gene_id[entrez_to_ensembl$ensembl_gene_id %in% ensembl_ids]
neuro_ensembl <- unique(neuro_ensembl)

# Create a mapping for analysis with pathway information
neuro_mapping <- data.frame(
  EnsemblID = neuro_ensembl,
  Symbol = mapping_df$Symbol[match(neuro_ensembl, mapping_df$EnsemblID)],
  Source = "KEGG pathways",
  stringsAsFactors = FALSE
)

# Save pathway details for reference
pathway_details <- data.frame(
  PathwayID = rep(names(neuro_genes_kegg), sapply(neuro_genes_kegg, nrow)),
  rbind.fill(neuro_genes_kegg)
)
write.csv(pathway_details, file.path(output_dir, "kegg_pathway_details.csv"), row.names = FALSE)

cat("Found", length(neuro_ensembl), "neuroinflammation-related genes in dataset\n")

# Save neuroinflammation gene list
write.csv(neuro_mapping, file.path(output_dir, "neuroinflammation_gene_list.csv"), row.names = FALSE)

# ----------------------
# 4. EXPLORATORY ANALYSIS
# ----------------------

# Extract expression data for neuroinflammation genes
neuro_expr <- logcpm_filtered_norm[neuro_ensembl, ]
neuro_expr_symbols <- neuro_expr
rownames(neuro_expr_symbols) <- neuro_mapping$Symbol

# Save expression data
write.csv(neuro_expr_symbols, file.path(output_dir, "neuroinflammation_expression.csv"))

# Create annotation for heatmap
anno <- data.frame(
  Diagnosis = metadata_df$diagnosis,
  Region = metadata_df$region,
  row.names = colnames(neuro_expr)
)

anno_colors <- list(
  Diagnosis = c("CONT" = "blue", "OUD" = "red"),
  Region = c("NAC" = "purple", "DLPFC" = "green")
)

# Generate heatmap
tryCatch({
  png(file.path(output_dir, "neuroinflammation_heatmap.png"), width = 1200, height = 1600, res = 120)
  pheatmap(
    neuro_expr_symbols,
    annotation_col = anno,
    annotation_colors = anno_colors,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_colnames = FALSE,
    fontsize_row = 10,
    main = "Neuroinflammation-Related Gene Expression"
  )
  dev.off()
}, error = function(e) {
  cat("Error generating heatmap:", e$message, "\n")
})

# PCA analysis
tryCatch({
  pca_neuro <- prcomp(t(neuro_expr))
  var_explained <- (pca_neuro$sdev^2) / sum(pca_neuro$sdev^2) * 100

  pca_data <- data.frame(
    PC1 = pca_neuro$x[,1],
    PC2 = pca_neuro$x[,2],
    Diagnosis = metadata_df$diagnosis,
    Region = metadata_df$region
  )

  png(file.path(output_dir, "neuroinflammation_pca.png"), width = 900, height = 700, res = 100)
  ggplot(pca_data, aes(x = PC1, y = PC2, color = Diagnosis, shape = Region)) +
    geom_point(size = 3) +
    labs(
      title = "PCA of Neuroinflammation Genes",
      x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
      y = paste0("PC2 (", round(var_explained[2], 1), "%)")
    ) +
    theme_bw() +
    scale_color_manual(values = c("CONT" = "blue", "OUD" = "red")) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.position = "right"
    )
  dev.off()
}, error = function(e) {
  cat("Error generating PCA plot:", e$message, "\n")
})

cat("✓ Neuroinflammation analysis complete. Output saved to:", output_dir, "\n")