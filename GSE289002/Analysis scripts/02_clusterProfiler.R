# Load required libraries
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)  # Mouse genome annotation
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(pheatmap)
library(RColorBrewer)

# Set up directories
project_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE289002"
de_results_dir <- file.path(project_dir, "DE_results")
enrichment_dir <- file.path(project_dir, "Enrichment_results")
dir.create(enrichment_dir, showWarnings = FALSE)

# Load gene mapping if available
gene_mapping_file <- file.path(de_results_dir, "ensembl_to_gene_symbol_mapping.csv")
if (file.exists(gene_mapping_file)) {
  gene_mapping <- read.csv(gene_mapping_file)
  # Create lookup functions
  ensembl_to_symbol <- setNames(gene_mapping$external_gene_name, gene_mapping$ensembl_gene_id)
  ensembl_to_entrez <- NULL

  # Map Ensembl IDs to Entrez IDs using org.Mm.eg.db
  ensembl_ids <- gene_mapping$ensembl_gene_id
  ensembl2entrez <- bitr(ensembl_ids, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  ensembl_to_entrez <- setNames(ensembl2entrez$ENTREZID, ensembl2entrez$ENSEMBL)

  cat("Loaded gene mappings for", length(ensembl_to_entrez), "genes\n")
} else {
  stop("Gene mapping file not found. Run the edgeR script first.")
}

# Define function to perform enrichment analysis
perform_enrichment <- function(de_file, output_dir, contrast_name, region) {
  cat("Processing", contrast_name, "in region", region, "\n")

  # Create output directory
  contrast_dir <- file.path(output_dir, paste0(region, "_", contrast_name))
  dir.create(contrast_dir, showWarnings = FALSE)

  # Load DE results
  de_results <- read.csv(de_file)

  # Extract significant genes (FDR < 0.05)
  sig_genes <- de_results[de_results$FDR < 0.05, ]

  # Skip if no significant genes
  if (nrow(sig_genes) == 0) {
    cat("No significant genes for", contrast_name, "in", region, "\n")
    return(NULL)
  }

  # Convert Ensembl IDs to Entrez IDs
  sig_entrez <- ensembl_to_entrez[sig_genes$ensembl_id]
  sig_entrez <- sig_entrez[!is.na(sig_entrez)]

  if (length(sig_entrez) == 0) {
    cat("No Entrez IDs found for significant genes in", contrast_name, "\n")
    return(NULL)
  }

  # Create ranked gene list for GSEA (optional)
  all_genes <- de_results$ensembl_id
  gene_list <- de_results$logFC
  names(gene_list) <- ensembl_to_entrez[all_genes]
  gene_list <- gene_list[!is.na(names(gene_list))]
  gene_list <- sort(gene_list, decreasing = TRUE)

  # GO enrichment analysis - Biological Process
  go_bp <- enrichGO(gene = sig_entrez,
                    OrgDb = org.Mm.eg.db,
                    keyType = "ENTREZID",
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2)

  if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
    # Save results
    write.csv(go_bp@result, file.path(contrast_dir, "GO_BP_results.csv"), row.names = FALSE)

    # Plot top results
    png(file.path(contrast_dir, "GO_BP_dotplot.png"), width = 1000, height = 800, res = 100)
    print(dotplot(go_bp, showCategory = 20, title = paste0("GO BP - ", region, " - ", contrast_name)))
    dev.off()

    png(file.path(contrast_dir, "GO_BP_enrichment_map.png"), width = 1200, height = 1000, res = 100)
    print(emapplot(pairwise_termsim(go_bp), showCategory = 30))
    dev.off()
  }

  # GO enrichment analysis - Molecular Function
  go_mf <- enrichGO(gene = sig_entrez,
                    OrgDb = org.Mm.eg.db,
                    keyType = "ENTREZID",
                    ont = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2)

  if (!is.null(go_mf) && nrow(go_mf@result) > 0) {
    write.csv(go_mf@result, file.path(contrast_dir, "GO_MF_results.csv"), row.names = FALSE)

    png(file.path(contrast_dir, "GO_MF_dotplot.png"), width = 1000, height = 800, res = 100)
    print(dotplot(go_mf, showCategory = 20, title = paste0("GO MF - ", region, " - ", contrast_name)))
    dev.off()
  }

  # GO enrichment analysis - Cellular Component
  go_cc <- enrichGO(gene = sig_entrez,
                    OrgDb = org.Mm.eg.db,
                    keyType = "ENTREZID",
                    ont = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2)

  if (!is.null(go_cc) && nrow(go_cc@result) > 0) {
    write.csv(go_cc@result, file.path(contrast_dir, "GO_CC_results.csv"), row.names = FALSE)

    png(file.path(contrast_dir, "GO_CC_dotplot.png"), width = 1000, height = 800, res = 100)
    print(dotplot(go_cc, showCategory = 20, title = paste0("GO CC - ", region, " - ", contrast_name)))
    dev.off()
  }

  # KEGG pathway analysis
  kegg <- enrichKEGG(gene = sig_entrez,
                     organism = "mmu",  # Mouse
                     keyType = "ncbi-geneid",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05)

  if (!is.null(kegg) && nrow(kegg@result) > 0) {
    write.csv(kegg@result, file.path(contrast_dir, "KEGG_results.csv"), row.names = FALSE)

    png(file.path(contrast_dir, "KEGG_dotplot.png"), width = 1000, height = 800, res = 100)
    print(dotplot(kegg, showCategory = 20, title = paste0("KEGG pathways - ", region, " - ", contrast_name)))
    dev.off()

    # Create pathway heatmap for top KEGG pathways
    if (nrow(kegg@result) >= 5) {
      top_pathways <- head(kegg@result$ID, 5)

      tryCatch({
        # Download KEGG pathway data for visualization
        pathway_data <- download_KEGG(top_pathways, "mmu")

        for (i in 1:length(pathway_data)) {
          pathway_id <- names(pathway_data)[i]
          pathway_name <- kegg@result$Description[kegg@result$ID == pathway_id]

          if (length(pathway_name) > 0) {
            png(file.path(contrast_dir, paste0("KEGG_pathway_", pathway_id, ".png")),
                width = 1200, height = 1000, res = 100)
            print(pathview(gene.data = gene_list,
                           pathway.id = pathway_id,
                           species = "mmu",
                           limit = list(gene = 2, cpd = 1)))
            dev.off()
          }
        }
      }, error = function(e) {
        cat("Error in pathway visualization:", conditionMessage(e), "\n")
      })
    }
  }

  # Return the enrichment results
  return(list(go_bp = go_bp, go_mf = go_mf, go_cc = go_cc, kegg = kegg))
}

# Process all regions and contrasts
regions <- list.dirs(de_results_dir, full.names = FALSE, recursive = FALSE)
regions <- regions[grep("^region_", regions)]

# Process each brain region
for (region_folder in regions) {
  region <- gsub("^region_", "", region_folder)
  region_path <- file.path(de_results_dir, region_folder)
  region_output <- file.path(enrichment_dir, region)
  dir.create(region_output, showWarnings = FALSE)

  # Get all CSV files (DE results)
  de_files <- list.files(region_path, pattern = "\\.csv$", full.names = TRUE)

  # Process each contrast
  for (file in de_files) {
    # Skip files that are not contrast results
    if (!grepl("vs|Interaction", basename(file))) {
      next
    }

    # Extract contrast name
    contrast_name <- gsub("\\.csv$", "", basename(file))

    # Perform enrichment analysis
    perform_enrichment(file, region_output, contrast_name, region)
  }
}

# Analyze complement-related genes specifically
complement_dir <- file.path(de_results_dir, "complement_genes")
if (dir.exists(complement_dir)) {
  # Load complement gene mapping
  complement_mapping_file <- file.path(complement_dir, "complement_genes_mapping.csv")
  if (file.exists(complement_mapping_file)) {
    complement_mapping <- read.csv(complement_mapping_file)
    complement_ensembl <- complement_mapping$ensembl_gene_id

    # Create a directory for complement gene analysis results
    complement_results <- file.path(enrichment_dir, "complement_analysis")
    dir.create(complement_results, showWarnings = FALSE)

    # Extract expression data for complement genes across all contrasts and regions
    complement_summary <- data.frame()

    for (region_folder in regions) {
      region <- gsub("^region_", "", region_folder)
      region_path <- file.path(de_results_dir, region_folder)

      # Get all CSV files (DE results)
      de_files <- list.files(region_path, pattern = "\\.csv$", full.names = TRUE)

      # Process each contrast
      for (file in de_files) {
        if (!grepl("vs|Interaction", basename(file))) {
          next
        }

        contrast_name <- gsub("\\.csv$", "", basename(file))
        de_results <- read.csv(file)

        # Filter for complement genes
        complement_de <- de_results[de_results$ensembl_id %in% complement_ensembl, ]

        if (nrow(complement_de) > 0) {
          # Add region and contrast info
          complement_de$region <- region
          complement_de$contrast <- contrast_name

          # Append to summary
          complement_summary <- rbind(complement_summary, complement_de)
        }
      }
    }

    # Save complement gene summary
    if (nrow(complement_summary) > 0) {
      write.csv(complement_summary, file.path(complement_results, "complement_DE_summary.csv"), row.names = FALSE)

      # Create a heatmap of logFC values
      complement_heatmap_data <- complement_summary %>%
        select(ensembl_id, gene_symbol, region, contrast, logFC) %>%
        pivot_wider(names_from = c(region, contrast), values_from = logFC)

      # Replace NA with 0 for visualization
      matrix_data <- as.matrix(complement_heatmap_data[, -c(1, 2)])
      rownames(matrix_data) <- complement_heatmap_data$gene_symbol

      # Generate heatmap
      png(file.path(complement_results, "complement_genes_heatmap.png"), width = 1200, height = 1000, res = 100)
      pheatmap(matrix_data,
               color = colorRampPalette(c("blue", "white", "red"))(100),
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               display_numbers = FALSE,
               fontsize_row = 10,
               fontsize_col = 8,
               angle_col = 45,
               main = "Complement Genes - Log2 Fold Change Across Conditions")
      dev.off()
    }
  }
}

# Generate a summary of all enrichment results
all_results <- data.frame()

for (region_folder in list.dirs(enrichment_dir, full.names = FALSE, recursive = FALSE)) {
  if (!file.exists(file.path(enrichment_dir, region_folder))) next

  for (contrast_folder in list.dirs(file.path(enrichment_dir, region_folder), full.names = FALSE, recursive = FALSE)) {
    # KEGG results
    kegg_file <- file.path(enrichment_dir, region_folder, contrast_folder, "KEGG_results.csv")
    if (file.exists(kegg_file)) {
      kegg_data <- read.csv(kegg_file)
      if (nrow(kegg_data) > 0) {
        kegg_data$region <- region_folder
        kegg_data$contrast <- contrast_folder
        kegg_data$database <- "KEGG"
        kegg_summary <- kegg_data[, c("region", "contrast", "database", "ID", "Description", "pvalue", "p.adjust", "Count")]
        all_results <- rbind(all_results, kegg_summary)
      }
    }

    # GO BP results
    go_bp_file <- file.path(enrichment_dir, region_folder, contrast_folder, "GO_BP_results.csv")
    if (file.exists(go_bp_file)) {
      go_data <- read.csv(go_bp_file)
      if (nrow(go_data) > 0) {
        go_data$region <- region_folder
        go_data$contrast <- contrast_folder
        go_data$database <- "GO_BP"
        go_summary <- go_data[, c("region", "contrast", "database", "ID", "Description", "pvalue", "p.adjust", "Count")]
        all_results <- rbind(all_results, go_summary)
      }
    }
  }
}

# Save overall summary
if (nrow(all_results) > 0) {
  write.csv(all_results, file.path(enrichment_dir, "all_enrichment_summary.csv"), row.names = FALSE)
}

cat("Enrichment analysis complete. Results saved in:", enrichment_dir, "\n")