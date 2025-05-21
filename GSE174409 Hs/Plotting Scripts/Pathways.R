# Load required libraries
suppressMessages({
  library(clusterProfiler)
  library(ggplot2)
  library(enrichplot)
  library(dplyr)
  library(biomaRt)
})

# Set file paths
base_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE174409/NeuroinflammationResults"
gene_list_file <- file.path(base_dir, "neuroinflammation_gene_list.csv")
metadata_file <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE174409/QC/metadata.rds"
output_dir <- base_dir

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Read the gene list and metadata
gene_list <- read.csv(gene_list_file, header = FALSE, stringsAsFactors = FALSE)
metadata <- readRDS(metadata_file)  # Load metadata from RDS file
ensembl_ids <- gene_list$V1  # Assuming the first column contains Ensembl IDs

# Convert Ensembl IDs to Entrez IDs
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
entrez_ids <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)
gene_ids <- na.omit(entrez_ids$entrezgene_id)  # Remove NA values

# Function to perform enrichment analysis with retry mechanism
perform_enrichment <- function(group_name, group_gene_ids) {
  retry_count <- 3
  enrich_res <- NULL

  for (i in 1:retry_count) {
    tryCatch({
      enrich_res <- enrichKEGG(
        gene = group_gene_ids,
        organism = "hsa",
        keyType = "kegg",
        pvalueCutoff = 0.05,
        use_internal_data = FALSE
      )
      break
    }, error = function(e) {
      cat("Attempt", i, "failed for group:", group_name, "-", e$message, "\n")
      if (i == retry_count) stop("Failed after", retry_count, "attempts.")
    })
  }

  if (is.null(enrich_res) || nrow(as.data.frame(enrich_res)) == 0) {
    cat("No pathways enriched for group:", group_name, "\n")
    return(NULL)
  }

  # Save enrichment results
  write.csv(as.data.frame(enrich_res), file.path(output_dir, paste0("pathway_enrichment_results_", group_name, ".csv")))

  # Generate dot plot
  tryCatch({
    png(file.path(output_dir, paste0("pathway_enrichment_dotplot_", group_name, ".png")), width = 1200, height = 800, res = 120)
    print(dotplot(enrich_res, showCategory = 20) +
      ggtitle(paste("Pathway Enrichment -", group_name)) +
      theme(plot.title = element_text(hjust = 0.5, size = 14)))
    dev.off()
  }, error = function(e) {
    cat("Error generating dotplot for group:", group_name, "-", e$message, "\n")
  })

  # Generate bar plot
  tryCatch({
    png(file.path(output_dir, paste0("pathway_enrichment_barplot_", group_name, ".png")), width = 1200, height = 800, res = 120)
    print(barplot(enrich_res, showCategory = 20, title = paste("Pathway Enrichment -", group_name)) +
      theme(plot.title = element_text(hjust = 0.5, size = 14)))
    dev.off()
  }, error = function(e) {
    cat("Error generating barplot for group:", group_name, "-", e$message, "\n")
  })

  cat("✓ Pathway enrichment plots generated for group:", group_name, "\n")
}

# Perform enrichment analysis by sex and brain region
for (sex in unique(metadata$sex)) {
  for (region in unique(metadata$region)) {
    group_name <- paste(sex, region, sep = "_")
    group_gene_ids <- gene_ids[metadata$sex == sex & metadata$region == region]
    if (length(group_gene_ids) > 0) {
      perform_enrichment(group_name, group_gene_ids)
    } else {
      cat("No genes found for group:", group_name, "\n")
    }
  }
}