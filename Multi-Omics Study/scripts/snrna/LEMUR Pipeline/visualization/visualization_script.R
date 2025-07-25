# Visualization Script for LEMUR DE/Enrichment Results
# ---------------------------------------------------
# Standardized plots for each stratification (sex, region, cell type)
# Outputs organized by stratification and plot type

suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
  library(ComplexHeatmap)
  library(clusterProfiler)
  library(viridis)
  library(cowplot)
  library(dplyr)
})

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path("outputs", "visualization", timestamp)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Helper: ensure stratification subdirectory exists
ensure_strat_dir <- function(strat_name) {
  strat_dir <- file.path(output_dir, strat_name)
  if (!dir.exists(strat_dir)) dir.create(strat_dir, recursive = TRUE)
  return(strat_dir)
}

# 1. Load DE and enrichment results
de_files <- list.files("../de_enrichment", pattern = "^DE_.*\\.rds$", full.names = TRUE)
enrich_files <- list.files("../de_enrichment", pattern = "^enrichment_.*\\.rds$", full.names = TRUE)
gsea_files <- list.files("../de_enrichment", pattern = "^gsea_.*\\.rds$", full.names = TRUE)
rrho_files <- list.files("../de_enrichment", pattern = "^rrho_.*\\.rds$", full.names = TRUE)

# 2. Define plotting functions

plot_volcano <- function(de_table, strat_name, top_n = 10) {
  de_table$signif <- de_table$adj_pval < 0.05
  de_table$label <- ""
  top_genes <- rownames(de_table)[order(de_table$adj_pval)][1:top_n]
  de_table$label[rownames(de_table) %in% top_genes] <- rownames(de_table)[rownames(de_table) %in% top_genes]
  ggplot(de_table, aes(x = logFC, y = -log10(adj_pval), color = signif)) +
    geom_point(alpha = 0.7, size = 1) +
    scale_color_manual(values = c("grey70", "firebrick")) +
    geom_text_repel(aes(label = label), max.overlaps = Inf, size = 3, box.padding = 0.5) +
    labs(
      title = paste("Volcano Plot:", strat_name),
      x = "log2 Fold Change", y = "-log10 Adjusted p-value"
    ) +
    theme_minimal(base_size = 12)
}

plot_heatmap <- function(de_table, strat_name, top_n = 20) {
  top_genes <- rownames(de_table)[order(de_table$adj_pval)][1:top_n]
  # Placeholder: random expression matrix for heatmap (replace with real data)
  set.seed(42)
  expr_mat <- matrix(rnorm(top_n * 8), nrow = top_n)
  rownames(expr_mat) <- top_genes
  colnames(expr_mat) <- paste0("Sample", 1:8)
  pheatmap(expr_mat,
    cluster_rows = TRUE, cluster_cols = TRUE,
    color = viridis(100), main = paste("Top DE Genes Heatmap:", strat_name),
    fontsize_row = 8, fontsize_col = 8
  )
}

plot_bar_enrichment <- function(enrich_obj, strat_name, enrich_type, top_n = 10) {
  if (is.null(enrich_obj) || nrow(enrich_obj@result) == 0) {
    return(NULL)
  }
  df <- enrich_obj@result %>%
    arrange(p.adjust) %>%
    head(top_n)
  ggplot(df, aes(x = reorder(Description, -log10(p.adjust)), y = -log10(p.adjust), fill = p.adjust)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_viridis_c(option = "C") +
    labs(
      title = paste("Top", enrich_type, "Enrichment:", strat_name),
      x = enrich_type, y = "-log10 Adjusted p-value"
    ) +
    theme_minimal(base_size = 12)
}

plot_gsea <- function(gsea_obj, strat_name, pathway = NULL) {
  if (is.null(gsea_obj) || nrow(gsea_obj@result) == 0) {
    return(NULL)
  }
  # If pathway not specified, use top NES
  if (is.null(pathway)) pathway <- gsea_obj@result$ID[which.max(abs(gsea_obj@result$NES))]
  gseaplot2(gsea_obj, geneSetID = pathway, title = paste("GSEA:", pathway, strat_name))
}

plot_rrho <- function(rrho_obj, strat_name) {
  if (is.null(rrho_obj)) {
    return(NULL)
  }
  # Placeholder: RRHO heatmap (replace with real RRHO plot)
  heatmap(rrho_obj$hypermat, main = paste("RRHO Heatmap:", strat_name), col = viridis(100))
}

# 3. Loop through stratifications and generate plots
for (de_file in de_files) {
  strat_name <- gsub("DE_|\\.rds", "", basename(de_file))
  strat_dir <- ensure_strat_dir(strat_name)
  de_table <- readRDS(de_file)
  # Volcano
  pdf(file.path(strat_dir, paste0("volcano_", strat_name, ".pdf")), width = 7, height = 5)
  print(plot_volcano(de_table, strat_name))
  dev.off()
  # Heatmap
  pdf(file.path(strat_dir, paste0("heatmap_", strat_name, ".pdf")), width = 7, height = 7)
  plot_heatmap(de_table, strat_name)
  dev.off()
  # Bar plots for enrichment
  for (enrich_type in c("GO", "Reactome", "MSigDB")) {
    enrich_file <- grep(paste0("enrichment_", strat_name, "_", enrich_type), enrich_files, value = TRUE)
    if (length(enrich_file) > 0) {
      enrich_obj <- readRDS(enrich_file[1])
      pdf(file.path(strat_dir, paste0("barplot_", enrich_type, "_", strat_name, ".pdf")), width = 7, height = 5)
      print(plot_bar_enrichment(enrich_obj, strat_name, enrich_type))
      dev.off()
    }
  }
  # GSEA
  gsea_file <- grep(paste0("gsea_", strat_name), gsea_files, value = TRUE)
  if (length(gsea_file) > 0) {
    gsea_obj <- readRDS(gsea_file[1])
    pdf(file.path(strat_dir, paste0("gsea_", strat_name, ".pdf")), width = 7, height = 5)
    print(plot_gsea(gsea_obj, strat_name))
    dev.off()
  }
  # RRHO (if available)
  rrho_file <- grep(paste0("rrho_", strat_name), rrho_files, value = TRUE)
  if (length(rrho_file) > 0) {
    rrho_obj <- readRDS(rrho_file[1])
    pdf(file.path(strat_dir, paste0("rrho_", strat_name, ".pdf")), width = 7, height = 7)
    plot_rrho(rrho_obj, strat_name)
    dev.off()
  }
  # Save summary tables (top DE genes, top pathways)
  top_de <- de_table[order(de_table$adj_pval), ][1:20, ]
  write.csv(top_de, file.path(strat_dir, paste0("top_DE_genes_", strat_name, ".csv")), row.names = TRUE)
  # For enrichment, save top terms if available
  for (enrich_type in c("GO", "Reactome", "MSigDB")) {
    enrich_file <- grep(paste0("enrichment_", strat_name, "_", enrich_type), enrich_files, value = TRUE)
    if (length(enrich_file) > 0) {
      enrich_obj <- readRDS(enrich_file[1])
      if (!is.null(enrich_obj) && nrow(enrich_obj@result) > 0) {
        top_terms <- enrich_obj@result %>%
          arrange(p.adjust) %>%
          head(20)
        write.csv(top_terms, file.path(strat_dir, paste0("top_", enrich_type, "_terms_", strat_name, ".csv")), row.names = FALSE)
      }
    }
  }
}

# 4. Log session info
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

# 5. Note: All outputs are organized by stratification and plot type in the visualization output directory.
#         Figures are saved as PDF, summary tables as CSV, and session info for reproducibility.
