#!/usr/bin/env Rscript
#' @title LEMUR Analysis Utility Functions
#' @description Helper functions for LEMUR analysis pipeline
#' @author Generated Analysis Pipeline
#' @date 2024

# =============================================================================
# UTILITY FUNCTIONS FOR LEMUR ANALYSIS
# =============================================================================

#' Check and install required packages
#' @param packages Character vector of package names
#' @return NULL (installs packages silently)
check_and_install_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      message("Installing package: ", pkg)
      if (pkg %in% rownames(available.packages())) {
        install.packages(pkg, dependencies = TRUE, quiet = TRUE)
      } else {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager", quiet = TRUE)
        }
        BiocManager::install(pkg, dependencies = TRUE, quiet = TRUE)
      }
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    }
  }
}

#' Comprehensive QC plotting function
#' @param sce SingleCellExperiment object
#' @param output_dir Output directory for plots
#' @param prefix File prefix for plots
#' @return NULL (saves plots)
plot_qc_metrics <- function(sce, output_dir, prefix = "qc") {
  library(ggplot2)
  library(patchwork)

  # Extract QC metrics
  qc_data <- data.frame(
    total_counts = sce$total_counts,
    n_genes = sce$n_genes,
    pct_mito = sce$subsets_Mito_percent,
    condition = sce$condition,
    donor_id = sce$donor_id
  )

  # 1. Distribution plots
  p1 <- ggplot(qc_data, aes(x = log10(total_counts), fill = condition)) +
    geom_density(alpha = 0.7) +
    scale_fill_viridis_d() +
    theme_minimal() +
    labs(title = "Total UMI Counts Distribution", x = "log10(Total Counts)")

  p2 <- ggplot(qc_data, aes(x = n_genes, fill = condition)) +
    geom_density(alpha = 0.7) +
    scale_fill_viridis_d() +
    theme_minimal() +
    labs(title = "Number of Genes Distribution", x = "Number of Genes")

  p3 <- ggplot(qc_data, aes(x = pct_mito, fill = condition)) +
    geom_density(alpha = 0.7) +
    scale_fill_viridis_d() +
    theme_minimal() +
    labs(title = "Mitochondrial Gene %", x = "Mitochondrial %")

  # 2. Scatter plots
  p4 <- ggplot(qc_data, aes(x = log10(total_counts), y = n_genes, color = condition)) +
    geom_point(alpha = 0.6, size = 0.5) +
    scale_color_viridis_d() +
    theme_minimal() +
    labs(title = "Genes vs UMI", x = "log10(Total Counts)", y = "Number of Genes")

  p5 <- ggplot(qc_data, aes(x = log10(total_counts), y = pct_mito, color = condition)) +
    geom_point(alpha = 0.6, size = 0.5) +
    scale_color_viridis_d() +
    theme_minimal() +
    labs(title = "Mitochondrial % vs UMI", x = "log10(Total Counts)", y = "Mitochondrial %")

  # 3. Box plots by condition
  p6 <- ggplot(qc_data, aes(x = condition, y = log10(total_counts), fill = condition)) +
    geom_boxplot(alpha = 0.7) +
    scale_fill_viridis_d() +
    theme_minimal() +
    labs(title = "UMI by Condition", x = "Condition", y = "log10(Total Counts)")

  # Combine plots
  combined <- (p1 + p2 + p3) / (p4 + p5 + p6)

  # Save
  ggsave(
    filename = file.path(output_dir, paste0(prefix, "_metrics.png")),
    plot = combined,
    width = 15, height = 10, dpi = 300
  )

  return(combined)
}

#' Create LEMUR diagnostic plots
#' @param fit LEMUR fit object
#' @param output_dir Output directory
#' @param prefix File prefix
#' @return NULL (saves plots)
plot_lemur_diagnostics <- function(fit, output_dir, prefix = "lemur") {
  library(ggplot2)
  library(patchwork)

  # Extract embedding
  embedding <- fit$embedding

  # 1. Embedding scatter plots (first few dimensions)
  plot_list <- list()
  for (i in 1:min(4, ncol(embedding))) {
    for (j in (i + 1):min(5, ncol(embedding))) {
      if (j <= ncol(embedding)) {
        df <- data.frame(
          x = embedding[, i],
          y = embedding[, j],
          condition = fit$colData$condition
        )

        p <- ggplot(df, aes(x = x, y = y, color = condition)) +
          geom_point(alpha = 0.6, size = 0.5) +
          scale_color_viridis_d() +
          theme_minimal() +
          labs(
            title = paste("LEMUR Embedding:", "Dim", i, "vs Dim", j),
            x = paste("Dimension", i),
            y = paste("Dimension", j)
          )

        plot_list[[paste0("dim_", i, "_", j)]] <- p
      }
    }
  }

  # 2. Explained variance plot
  if ("explained_variance" %in% names(fit)) {
    var_df <- data.frame(
      dimension = 1:length(fit$explained_variance),
      variance = fit$explained_variance,
      cumulative = cumsum(fit$explained_variance)
    )

    p_var <- ggplot(var_df, aes(x = dimension, y = variance)) +
      geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
      theme_minimal() +
      labs(
        title = "LEMUR Embedding - Explained Variance",
        x = "Dimension",
        y = "Explained Variance"
      )

    plot_list[["variance"]] <- p_var
  }

  # Combine and save
  if (length(plot_list) > 0) {
    combined <- wrap_plots(plot_list, ncol = 2)
    ggsave(
      filename = file.path(output_dir, paste0(prefix, "_diagnostics.png")),
      plot = combined,
      width = 12, height = 8, dpi = 300
    )
  }
}

#' Enhanced volcano plot with gene labels
#' @param results_table Data frame with DE results
#' @param output_dir Output directory
#' @param prefix File prefix
#' @param top_n Number of top genes to label
#' @param fdr_threshold FDR threshold for significance
#' @return ggplot object
enhanced_volcano_plot <- function(results_table, output_dir, prefix = "enhanced_volcano",
                                  top_n = 20, fdr_threshold = 0.05) {
  library(ggplot2)
  library(ggrepel)

  # Prepare data
  volcano_data <- results_table[!is.na(results_table$log_fc) & !is.na(results_table$adj_pvalue), ]
  volcano_data$neg_log10_padj <- -log10(volcano_data$adj_pvalue)
  volcano_data$significant <- volcano_data$adj_pvalue < fdr_threshold

  # Identify top genes to label
  sig_genes <- volcano_data[volcano_data$significant, ]
  if (nrow(sig_genes) > 0) {
    top_genes <- head(sig_genes[order(sig_genes$adj_pvalue), ], top_n)
    volcano_data$label <- ifelse(volcano_data$gene %in% top_genes$gene, volcano_data$gene, "")
  } else {
    volcano_data$label <- ""
  }

  # Create enhanced volcano plot
  p <- ggplot(volcano_data, aes(x = log_fc, y = neg_log10_padj)) +
    geom_point(aes(color = significant), alpha = 0.6, size = 0.8) +
    scale_color_manual(
      values = c("FALSE" = "grey70", "TRUE" = "red"),
      name = "Significant"
    ) +
    geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "blue", alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
    geom_text_repel(
      aes(label = label),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.3,
      point.padding = 0.3
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "Enhanced Volcano Plot - LEMUR DE Results",
      subtitle = paste("OUD vs Control | FDR <", fdr_threshold, "| Top", top_n, "genes labeled"),
      x = "Log2 Fold Change",
      y = "-log10(Adjusted P-value)"
    )

  # Save plot
  ggsave(
    filename = file.path(output_dir, paste0(prefix, ".png")),
    plot = p,
    width = 10, height = 8, dpi = 300
  )

  return(p)
}

#' Create gene expression heatmap for top DE genes
#' @param sce SingleCellExperiment object
#' @param results_table DE results table
#' @param output_dir Output directory
#' @param top_n Number of top genes to include
#' @param prefix File prefix
#' @return ComplexHeatmap object
plot_de_heatmap <- function(sce, results_table, output_dir, top_n = 50, prefix = "de_heatmap") {
  library(ComplexHeatmap)
  library(circlize)

  # Get top DE genes
  sig_results <- results_table[results_table$significant, ]
  if (nrow(sig_results) == 0) {
    message("No significant DE genes found for heatmap")
    return(NULL)
  }

  top_genes <- head(sig_results[order(sig_results$adj_pvalue), ], top_n)

  # Extract expression matrix
  gene_indices <- match(top_genes$gene, rownames(sce))
  gene_indices <- gene_indices[!is.na(gene_indices)]

  if (length(gene_indices) == 0) {
    message("No matching genes found in expression matrix")
    return(NULL)
  }

  expr_matrix <- logcounts(sce)[gene_indices, ]

  # Prepare annotations
  col_annotation <- data.frame(
    Condition = sce$condition,
    Donor = sce$donor_id
  )

  # Create color mappings
  condition_colors <- c("Control" = "blue", "OUD" = "red")
  donor_colors <- rainbow(length(unique(sce$donor_id)))
  names(donor_colors) <- unique(sce$donor_id)

  ha_column <- HeatmapAnnotation(
    Condition = col_annotation$Condition,
    Donor = col_annotation$Donor,
    col = list(
      Condition = condition_colors,
      Donor = donor_colors
    )
  )

  # Create heatmap
  ht <- Heatmap(
    expr_matrix,
    name = "log2(counts+1)",
    top_annotation = ha_column,
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_title = paste("Top", nrow(expr_matrix), "DE Genes"),
    row_title = "Genes",
    clustering_distance_rows = "pearson",
    clustering_distance_columns = "pearson"
  )

  # Save heatmap
  png(
    filename = file.path(output_dir, paste0(prefix, ".png")),
    width = 12, height = 10, units = "in", res = 300
  )
  draw(ht)
  dev.off()

  return(ht)
}

#' Generate pathway enrichment analysis
#' @param de_genes Character vector of DE gene names
#' @param background_genes Character vector of background genes
#' @param organism Organism for analysis ("human" or "mouse")
#' @return Data frame with enrichment results
run_pathway_enrichment <- function(de_genes, background_genes = NULL, organism = "human") {
  # This is a placeholder function - implement with your preferred enrichment method
  # e.g., clusterProfiler, fgsea, etc.

  message("Pathway enrichment analysis would be implemented here")
  message("DE genes provided: ", length(de_genes))

  # Return mock results for now
  data.frame(
    pathway = c("Mock Pathway 1", "Mock Pathway 2"),
    pvalue = c(0.001, 0.01),
    adj_pvalue = c(0.01, 0.05),
    genes = c("GENE1,GENE2", "GENE3,GENE4"),
    stringsAsFactors = FALSE
  )
}

#' Save LEMUR model and results
#' @param fit LEMUR fit object
#' @param results_table DE results table
#' @param output_dir Output directory
#' @return NULL (saves files)
save_lemur_model <- function(fit, results_table, output_dir) {
  # Save LEMUR fit object
  saveRDS(fit, file = file.path(output_dir, "lemur_fit.rds"))

  # Save results table
  write.csv(results_table, file = file.path(output_dir, "de_results_complete.csv"), row.names = FALSE)

  # Save embedding
  write.csv(fit$embedding, file = file.path(output_dir, "lemur_embedding.csv"))

  # Save harmony embedding if available
  if ("harmony_embedding" %in% names(fit)) {
    write.csv(fit$harmony_embedding, file = file.path(output_dir, "harmony_embedding.csv"))
  }

  message("LEMUR model and results saved to: ", output_dir)
}

#' Load and validate H5AD file
#' @param file_path Path to H5AD file
#' @param required_cols Required metadata columns
#' @return SingleCellExperiment object
load_and_validate_h5ad <- function(file_path, required_cols = c("donor_id", "condition")) {
  library(zellkonverter)

  # Check file exists
  if (!file.exists(file_path)) {
    stop("H5AD file not found: ", file_path)
  }

  # Load file
  message("Loading H5AD file: ", file_path)
  sce <- readH5AD(file_path, use_hdf5 = TRUE)

  # Validate required columns
  missing_cols <- setdiff(required_cols, colnames(colData(sce)))
  if (length(missing_cols) > 0) {
    stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
  }

  message("✅ H5AD file loaded successfully")
  message("📊 Dimensions: ", nrow(sce), " genes × ", ncol(sce), " cells")
  message("🏷️  Available metadata: ", paste(colnames(colData(sce)), collapse = ", "))

  return(sce)
}

#' Print analysis summary
#' @param sce SingleCellExperiment object
#' @param results_table DE results table
#' @param neighborhoods Neighborhood results
#' @param params Analysis parameters
#' @return NULL (prints summary)
print_analysis_summary <- function(sce, results_table, neighborhoods = NULL, params = NULL) {
  cat("\n" + rep("=", 60) + "\n")
  cat("📊 LEMUR ANALYSIS SUMMARY\n")
  cat(rep("=", 60) + "\n\n")

  cat("📈 Data Overview:\n")
  cat("  • Genes:", nrow(sce), "\n")
  cat("  • Cells:", ncol(sce), "\n")
  cat("  • Conditions:", paste(unique(sce$condition), collapse = ", "), "\n")
  cat("  • Donors:", length(unique(sce$donor_id)), "\n\n")

  cat("🔬 DE Analysis Results:\n")
  cat("  • Total tests:", nrow(results_table), "\n")
  cat("  • Significant genes:", sum(results_table$significant, na.rm = TRUE), "\n")
  cat("  • FDR threshold:", ifelse(is.null(params), "0.05", params$fdr_threshold), "\n\n")

  if (!is.null(neighborhoods)) {
    cat("🏘️  Neighborhood Analysis:\n")
    cat("  • DE neighborhoods found:", length(neighborhoods$de_neighborhoods), "\n\n")
  }

  cat("🎯 Top 5 DE Genes:\n")
  top_5 <- head(results_table[results_table$significant, ], 5)
  if (nrow(top_5) > 0) {
    for (i in 1:nrow(top_5)) {
      cat(sprintf(
        "  %d. %s (log2FC: %.3f, adj.p: %.2e)\n",
        i, top_5$gene[i], top_5$log_fc[i], top_5$adj_pvalue[i]
      ))
    }
  } else {
    cat("  No significant genes found\n")
  }

  cat("\n" + rep("=", 60) + "\n")
}
