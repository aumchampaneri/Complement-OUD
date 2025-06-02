# ==============================================================================
# Cross-Dataset Integration: GSE174409 vs GSE225158
# ==============================================================================
# Purpose: Integrate and compare differential expression results between datasets
# Dataset 1: GSE174409 - Bulk RNA-seq (ACC + NAc)
# Dataset 2: GSE225158 - Pseudobulk snRNA-seq (Caudate + Putamen)
# Methods: Gene overlap analysis, effect size correlation, meta-analysis
# ==============================================================================

# Configuration
BASE_DIR <- "/Users/aumchampaneri/Complement-OUD/Translational_Study/Homo_sapiens"
GSE174409_DIR <- file.path(BASE_DIR, "Results", "GSE174409")
GSE225158_DIR <- file.path(BASE_DIR, "Results", "GSE225158")
INTEGRATION_DIR <- file.path(BASE_DIR, "Results", "Integration")

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(VennDiagram)
  library(biomaRt)
})

# ==============================================================================
# GENE ID HARMONIZATION
# ==============================================================================

#' Convert Ensembl IDs to gene symbols using biomaRt
convert_ensembl_to_symbols <- function(ensembl_ids) {
  cat("Converting Ensembl IDs to gene symbols...\n")
  
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  gene_map <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = ensembl
  )
  
  # Remove empty symbols
  gene_map <- gene_map[gene_map$hgnc_symbol != "", ]
  
  cat(sprintf("Mapped %d/%d Ensembl IDs to gene symbols\n", 
              nrow(gene_map), length(ensembl_ids)))
  
  return(gene_map)
}

#' Harmonize gene identifiers between datasets
harmonize_gene_ids <- function(gse174409_results, gse225158_results) {
  cat("\n=== Harmonizing Gene Identifiers ===\n")
  
  # Get all Ensembl IDs from GSE174409
  all_ensembl_ids <- unique(unlist(lapply(gse174409_results, function(method) {
    rownames(method$results)
  })))
  
  # Convert to gene symbols
  gene_map <- convert_ensembl_to_symbols(all_ensembl_ids)
  
  # Convert GSE174409 results to use gene symbols
  gse174409_harmonized <- list()
  for (method_name in names(gse174409_results)) {
    method_results <- gse174409_results[[method_name]]$results
    
    # Merge with gene mapping
    method_results$ensembl_gene_id <- rownames(method_results)
    harmonized <- merge(method_results, gene_map, by = "ensembl_gene_id", all.x = FALSE)
    
    # Remove duplicates (keep best p-value)
    pval_col <- ifelse("adj.P.Val" %in% colnames(harmonized), "adj.P.Val", "padj")
    harmonized <- harmonized %>%
      group_by(hgnc_symbol) %>%
      slice_min(order_by = get(pval_col), n = 1, with_ties = FALSE) %>%
      ungroup()
    
    # Set gene symbols as row names
    harmonized <- as.data.frame(harmonized)
    rownames(harmonized) <- harmonized$hgnc_symbol
    harmonized$ensembl_gene_id <- NULL
    harmonized$hgnc_symbol <- NULL
    
    gse174409_harmonized[[method_name]] <- list(
      results = harmonized,
      fit = gse174409_results[[method_name]]$fit,
      correlation = gse174409_results[[method_name]]$correlation
    )
    
    cat(sprintf("Method %s: %d genes after harmonization\n", 
                method_name, nrow(harmonized)))
  }
  
  cat("✓ Gene ID harmonization complete\n")
  
  return(list(
    gse174409 = gse174409_harmonized,
    gse225158 = gse225158_results
  ))
}

# ==============================================================================
# DATA LOADING
# ==============================================================================

#' Load results from both datasets
load_integration_data <- function() {
  cat("=== Loading Cross-Dataset Results ===\n")
  
  # Check files exist
  gse174409_file <- file.path(GSE174409_DIR, "GSE174409_region_analysis_results.rds")
  gse225158_file <- file.path(GSE225158_DIR, "GSE225158_region_analysis_results.rds")
  
  if (!file.exists(gse174409_file)) {
    stop("GSE174409 results file not found: ", gse174409_file)
  }
  if (!file.exists(gse225158_file)) {
    stop("GSE225158 results file not found: ", gse225158_file)
  }
  
  # Load results
  gse174409_results <- readRDS(gse174409_file)
  gse225158_results <- readRDS(gse225158_file)
  
  cat("✓ GSE174409 loaded:", length(gse174409_results), "methods\n")
  cat("✓ GSE225158 loaded:", length(gse225158_results), "methods\n")
  
  # Harmonize gene identifiers
  harmonized_data <- harmonize_gene_ids(gse174409_results, gse225158_results)
  
  return(harmonized_data)
}

# ==============================================================================
# GENE OVERLAP ANALYSIS
# ==============================================================================

#' Analyze gene overlap between datasets
analyze_gene_overlap <- function(data) {
  cat("\n=== Gene Overlap Analysis ===\n")
  
  if (!dir.exists(INTEGRATION_DIR)) dir.create(INTEGRATION_DIR, recursive = TRUE)
  
  methods <- c("paired_limma", "mixed_effects", "deseq2")
  overlap_results <- list()
  
  for (method in methods) {
    cat("Analyzing method:", method, "\n")
    
    # Get results for this method
    gse174409_res <- data$gse174409[[method]]$results
    gse225158_res <- data$gse225158[[method]]$results
    
    # Get p-value column names
    pval_col_174409 <- ifelse("adj.P.Val" %in% colnames(gse174409_res), "adj.P.Val", "padj")
    pval_col_225158 <- ifelse("adj.P.Val" %in% colnames(gse225158_res), "adj.P.Val", "padj")
    
    # Get significant genes
    sig_174409 <- rownames(gse174409_res)[gse174409_res[[pval_col_174409]] < 0.05 & !is.na(gse174409_res[[pval_col_174409]])]
    sig_225158 <- rownames(gse225158_res)[gse225158_res[[pval_col_225158]] < 0.05 & !is.na(gse225158_res[[pval_col_225158]])]
    
    # Calculate overlap
    common_genes <- intersect(sig_174409, sig_225158)
    total_tested <- length(intersect(rownames(gse174409_res), rownames(gse225158_res)))
    
    overlap_results[[method]] <- list(
      gse174409_sig = sig_174409,
      gse225158_sig = sig_225158,
      common_sig = common_genes,
      total_tested = total_tested,
      overlap_count = length(common_genes),
      overlap_pct = round(length(common_genes) / min(length(sig_174409), length(sig_225158)) * 100, 2)
    )
    
    cat(sprintf("  GSE174409 significant: %d genes\n", length(sig_174409)))
    cat(sprintf("  GSE225158 significant: %d genes\n", length(sig_225158)))
    cat(sprintf("  Overlap: %d genes (%.1f%%)\n", length(common_genes), overlap_results[[method]]$overlap_pct))
  }
  
  return(overlap_results)
}

#' Create Venn diagrams for gene overlaps
create_venn_diagrams <- function(overlap_results) {
  cat("\n=== Creating Venn Diagrams ===\n")
  
  for (method in names(overlap_results)) {
    venn_file <- file.path(INTEGRATION_DIR, paste0("venn_", method, ".png"))
    
    png(venn_file, width = 800, height = 600)
    venn.plot <- venn.diagram(
      x = list(
        GSE174409 = overlap_results[[method]]$gse174409_sig,
        GSE225158 = overlap_results[[method]]$gse225158_sig
      ),
      category.names = c("GSE174409\n(Bulk RNA-seq)", "GSE225158\n(snRNA-seq)"),
      filename = NULL,
      fill = c("#FF6B6B", "#4ECDC4"),
      alpha = 0.7,
      cex = 1.5,
      cat.cex = 1.2,
      main = paste("Gene Overlap:", method),
      main.cex = 1.8
    )
    grid.draw(venn.plot)
    dev.off()
    
    cat("✓ Venn diagram saved:", venn_file, "\n")
  }
}

# ==============================================================================
# EFFECT SIZE CORRELATION
# ==============================================================================

#' Correlate effect sizes between datasets
correlate_effect_sizes <- function(data) {
  cat("\n=== Effect Size Correlation Analysis ===\n")
  
  methods <- c("paired_limma", "mixed_effects", "deseq2")
  correlation_results <- list()
  
  for (method in methods) {
    cat("Analyzing method:", method, "\n")
    
    # Get results
    gse174409_res <- data$gse174409[[method]]$results
    gse225158_res <- data$gse225158[[method]]$results
    
    # Get common genes
    common_genes <- intersect(rownames(gse174409_res), rownames(gse225158_res))
    
    if (length(common_genes) < 100) {
      cat("  Warning: Only", length(common_genes), "common genes found\n")
      next
    }
    
    # Extract effect sizes
    fc_col_174409 <- ifelse("logFC" %in% colnames(gse174409_res), "logFC", "log2FoldChange")
    fc_col_225158 <- ifelse("logFC" %in% colnames(gse225158_res), "logFC", "log2FoldChange")
    
    # Create correlation data
    cor_data <- data.frame(
      gene = common_genes,
      gse174409_fc = gse174409_res[common_genes, fc_col_174409],
      gse225158_fc = gse225158_res[common_genes, fc_col_225158]
    ) %>%
      filter(!is.na(gse174409_fc) & !is.na(gse225158_fc))
    
    # Calculate correlation
    cor_result <- cor.test(cor_data$gse174409_fc, cor_data$gse225158_fc)
    
    correlation_results[[method]] <- list(
      data = cor_data,
      correlation = cor_result$estimate,
      p_value = cor_result$p.value,
      n_genes = nrow(cor_data)
    )
    
    cat(sprintf("  Correlation: r = %.3f (p = %.2e, n = %d)\n", 
                cor_result$estimate, cor_result$p.value, nrow(cor_data)))
    
    # Create scatter plot
    p <- ggplot(cor_data, aes(x = gse174409_fc, y = gse225158_fc)) +
      geom_point(alpha = 0.6, color = "steelblue") +
      geom_smooth(method = "lm", color = "red", se = TRUE) +
      labs(
        title = paste("Effect Size Correlation:", method),
        subtitle = sprintf("r = %.3f, p = %.2e, n = %d genes", 
                          cor_result$estimate, cor_result$p.value, nrow(cor_data)),
        x = "GSE174409 Log2 Fold Change",
        y = "GSE225158 Log2 Fold Change"
      ) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    
    plot_file <- file.path(INTEGRATION_DIR, paste0("correlation_", method, ".png"))
    ggsave(plot_file, p, width = 8, height = 6, dpi = 300)
    cat("✓ Correlation plot saved:", plot_file, "\n")
  }
  
  return(correlation_results)
}

# ==============================================================================
# SUMMARY AND REPORTING
# ==============================================================================

#' Create comprehensive integration summary
create_integration_summary <- function(overlap_results, correlation_results) {
  cat("\n=== Creating Integration Summary ===\n")
  
  # Create summary table
  summary_df <- data.frame(
    Method = names(overlap_results),
    GSE174409_Significant = sapply(overlap_results, function(x) length(x$gse174409_sig)),
    GSE225158_Significant = sapply(overlap_results, function(x) length(x$gse225158_sig)),
    Overlap_Count = sapply(overlap_results, function(x) x$overlap_count),
    Overlap_Percent = sapply(overlap_results, function(x) x$overlap_pct),
    Effect_Correlation = sapply(names(overlap_results), function(m) {
      if (m %in% names(correlation_results)) {
        round(correlation_results[[m]]$correlation, 3)
      } else {
        NA
      }
    }),
    Correlation_PValue = sapply(names(overlap_results), function(m) {
      if (m %in% names(correlation_results)) {
        correlation_results[[m]]$p_value
      } else {
        NA
      }
    })
  )
  
  # Save summary
  summary_file <- file.path(INTEGRATION_DIR, "integration_summary.csv")
  write.csv(summary_df, summary_file, row.names = FALSE)
  
  cat("✓ Integration summary saved:", summary_file, "\n")
  cat("\n=== Integration Summary ===\n")
  print(summary_df)
  
  return(summary_df)
}

# ==============================================================================
# MAIN PIPELINE
# ==============================================================================

#' Run complete cross-dataset integration analysis
run_integration_analysis <- function() {
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("CROSS-DATASET INTEGRATION ANALYSIS\n")
  cat("GSE174409 (Bulk RNA-seq) vs GSE225158 (snRNA-seq)\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  tryCatch({
    # Load and harmonize data
    integration_data <- load_integration_data()
    
    # Analyze gene overlaps
    overlap_results <- analyze_gene_overlap(integration_data)
    create_venn_diagrams(overlap_results)
    
    # Correlate effect sizes
    correlation_results <- correlate_effect_sizes(integration_data)
    
    # Create summary
    summary_df <- create_integration_summary(overlap_results, correlation_results)
    
    cat("\n", paste(rep("=", 70), collapse = ""), "\n")
    cat("SUCCESS: Integration analysis complete!\n")
    cat("Results saved to:", INTEGRATION_DIR, "\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")
    
    return(list(
      overlap = overlap_results,
      correlation = correlation_results,
      summary = summary_df
    ))
    
  }, error = function(e) {
    cat("\nERROR:", e$message, "\n")
    stop(e)
  })
}

# ==============================================================================
# EXECUTION
# ==============================================================================

if (!exists("SOURCED")) {
  integration_results <- run_integration_analysis()
}
