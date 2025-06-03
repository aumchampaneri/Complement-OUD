# ==============================================================================
# Quick Cross-Dataset Integration Analysis
# ==============================================================================

# Load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(VennDiagram)
  library(grid)
})

# Configuration
BASE_DIR <- "/Users/aumchampaneri/Complement-OUD/Translational_Study/Homo_sapiens"
GSE174409_DIR <- file.path(BASE_DIR, "Results", "GSE174409")
GSE225158_DIR <- file.path(BASE_DIR, "Results", "GSE225158")
INTEGRATION_DIR <- file.path(BASE_DIR, "Results", "Integration")

if (!dir.exists(INTEGRATION_DIR)) dir.create(INTEGRATION_DIR, recursive = TRUE)

#' Quick integration analysis
run_quick_integration <- function() {
  cat("=== Quick Cross-Dataset Integration Analysis ===\n")
  
  # Load both datasets
  gse174409_file <- file.path(GSE174409_DIR, "GSE174409_region_analysis_results.rds")
  gse225158_file <- file.path(GSE225158_DIR, "GSE225158_region_analysis_results.rds")
  
  if (!file.exists(gse174409_file) || !file.exists(gse225158_file)) {
    stop("Required result files not found")
  }
  
  gse174409_data <- readRDS(gse174409_file)
  gse225158_data <- readRDS(gse225158_file)
  
  cat("✓ Data loaded successfully\n")
  
  # Quick overlap analysis
  methods <- c("paired_limma", "mixed_effects", "deseq2")
  overlap_summary <- data.frame()
  
  for (method in methods) {
    # Extract significant genes
    res174409 <- gse174409_data[[method]]$results
    res225158 <- gse225158_data[[method]]$results
    
    pval_col_174409 <- ifelse("adj.P.Val" %in% colnames(res174409), "adj.P.Val", "padj")
    pval_col_225158 <- ifelse("adj.P.Val" %in% colnames(res225158), "adj.P.Val", "padj")
    
    sig_174409 <- rownames(res174409)[res174409[[pval_col_174409]] < 0.05 & !is.na(res174409[[pval_col_174409]])]
    sig_225158 <- rownames(res225158)[res225158[[pval_col_225158]] < 0.05 & !is.na(res225158[[pval_col_225158]])]
    
    overlap_genes <- intersect(sig_174409, sig_225158)
    overlap_pct <- round(length(overlap_genes) / min(length(sig_174409), length(sig_225158)) * 100, 1)
    
    overlap_summary <- rbind(overlap_summary, data.frame(
      Method = method,
      GSE174409_Significant = length(sig_174409),
      GSE225158_Significant = length(sig_225158),
      Overlap_Count = length(overlap_genes),
      Overlap_Percent = overlap_pct
    ))
    
    cat(sprintf("%s: %d vs %d genes, %d overlap (%.1f%%)\n", 
                method, length(sig_174409), length(sig_225158), 
                length(overlap_genes), overlap_pct))
    
    # Create Venn diagram
    venn_file <- file.path(INTEGRATION_DIR, paste0("venn_", method, ".png"))
    
    png(venn_file, width = 800, height = 600)
    venn.plot <- venn.diagram(
      x = list(
        GSE174409 = sig_174409,
        GSE225158 = sig_225158
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
    
    cat("✓ Venn diagram saved:", basename(venn_file), "\n")
  }
  
  # Save summary
  write.csv(overlap_summary, file.path(INTEGRATION_DIR, "quick_integration_summary.csv"), row.names = FALSE)
  
  # Create summary plot
  p_summary <- ggplot(overlap_summary, aes(x = Method)) +
    geom_col(aes(y = GSE174409_Significant), alpha = 0.7, fill = "#FF6B6B", width = 0.4, position = position_nudge(x = -0.2)) +
    geom_col(aes(y = GSE225158_Significant), alpha = 0.7, fill = "#4ECDC4", width = 0.4, position = position_nudge(x = 0.2)) +
    geom_text(aes(y = GSE174409_Significant, label = GSE174409_Significant), 
              position = position_nudge(x = -0.2), vjust = -0.5, size = 3) +
    geom_text(aes(y = GSE225158_Significant, label = GSE225158_Significant), 
              position = position_nudge(x = 0.2), vjust = -0.5, size = 3) +
    labs(title = "Cross-Dataset Integration: Significant Genes by Method",
         subtitle = "GSE174409 (Bulk RNA-seq) vs GSE225158 (snRNA-seq)",
         x = "Analysis Method", y = "Number of Significant Genes") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(file.path(INTEGRATION_DIR, "method_comparison.png"), 
         p_summary, width = 10, height = 6, dpi = 300)
  
  cat("✓ Summary plot saved\n")
  cat("✓ Quick integration analysis complete!\n")
  cat("Results saved to:", INTEGRATION_DIR, "\n")
  
  return(overlap_summary)
}

# Run if sourced
if (!exists("SOURCED")) {
  quick_results <- run_quick_integration()
}
