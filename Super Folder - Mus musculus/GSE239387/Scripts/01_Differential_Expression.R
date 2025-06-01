# ========================================================================
# GSE239387 Differential Expression Analysis
# Control vs Morphine in Nucleus Accumbens - FPKM Data Analysis
# ========================================================================

# Load required libraries
library(limma)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano)

# Set working directory
setwd("/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus/GSE239387")

# Create output directories
if(!dir.exists("Outputs/02_Differential_Expression")) {
  dir.create("Outputs/02_Differential_Expression", recursive = TRUE)
}

# Load processed data
log_fpkm <- readRDS("Outputs/01_Processing_QC/Data/expression_log2fpkm.rds")
sample_metadata <- readRDS("Outputs/01_Processing_QC/Data/sample_metadata.rds")
annotation_data <- readRDS("Outputs/01_Processing_QC/Data/annotation_data.rds")

# Differential expression analysis using limma (appropriate for FPKM data)
design <- model.matrix(~ treatment, data = sample_metadata)
fit <- lmFit(log_fpkm, design)
fit <- eBayes(fit)

# Extract results
results <- topTable(fit, coef = "treatmentMorphine", number = Inf, sort.by = "P")
results$gene_id <- rownames(results)

# Add gene annotations
results <- results %>%
  left_join(annotation_data[, c("id", "Symbol", "Description")], 
            by = c("gene_id" = "id"))

# Save results
write.csv(results, "Outputs/02_Differential_Expression/GSE239387_DE_results.csv", 
          row.names = FALSE)

cat("Differential expression analysis complete for GSE239387\n")
cat("Significant genes (p < 0.05):", sum(results$P.Value < 0.05), "\n")
cat("Significant genes (FDR < 0.05):", sum(results$adj.P.Val < 0.05), "\n")

# ========================================================================
# ENHANCED VISUALIZATIONS AND ANALYSIS
# ========================================================================

# Create additional output directories
dir.create("Outputs/02_Differential_Expression/Plots", recursive = TRUE, showWarnings = FALSE)
dir.create("Outputs/02_Differential_Expression/Data", recursive = TRUE, showWarnings = FALSE)

# Function for safe plotting
safe_plot <- function(plot_func, filename, width = 10, height = 8) {
  tryCatch({
    png(paste0("Outputs/02_Differential_Expression/Plots/", filename, ".png"), 
        width = width, height = height, units = "in", res = 300)
    plot_func()
    dev.off()
    cat("✓ Plot saved:", filename, "\n")
  }, error = function(e) {
    dev.off()
    cat("✗ Error creating plot:", filename, ":", e$message, "\n")
  })
}

# 1. Volcano Plot
safe_plot(function() {
  # Basic volcano plot
  plot(results$logFC, -log10(results$P.Value),
       xlab = "Log2 Fold Change (Morphine vs Control)",
       ylab = "-Log10 P-value",
       main = "Volcano Plot: GSE239387 Differential Expression",
       pch = 16, cex = 0.6, col = "gray70")
  
  # Highlight significant genes
  sig_genes <- results$P.Value < 0.05 & abs(results$logFC) > 0.5
  points(results$logFC[sig_genes], -log10(results$P.Value[sig_genes]),
         col = ifelse(results$logFC[sig_genes] > 0, "red", "blue"), pch = 16, cex = 0.8)
  
  # Add reference lines
  abline(h = -log10(0.05), col = "gray", lty = 2)
  abline(v = c(-0.5, 0.5), col = "gray", lty = 2)
  
  # Add legend
  legend("topright", c("Upregulated", "Downregulated", "Non-significant"),
         col = c("red", "blue", "gray70"), pch = 16, cex = 0.8)
  
  # Label top genes
  top_genes <- head(results[order(results$P.Value), ], 10)
  text(top_genes$logFC, -log10(top_genes$P.Value), top_genes$Symbol,
       pos = 3, cex = 0.7, col = "black")
}, "01_Volcano_Plot")

# 2. MA Plot
safe_plot(function() {
  # Calculate average expression
  results$AveExpr <- rowMeans(log_fpkm[results$gene_id, ], na.rm = TRUE)
  
  plot(results$AveExpr, results$logFC,
       xlab = "Average Log2 Expression",
       ylab = "Log2 Fold Change (Morphine vs Control)",
       main = "MA Plot: GSE239387",
       pch = 16, cex = 0.6, col = "gray70")
  
  # Highlight significant genes
  sig_genes <- results$P.Value < 0.05 & abs(results$logFC) > 0.5
  points(results$AveExpr[sig_genes], results$logFC[sig_genes],
         col = ifelse(results$logFC[sig_genes] > 0, "red", "blue"), pch = 16, cex = 0.8)
  
  # Add reference line
  abline(h = 0, col = "gray", lty = 2)
  
  # Add legend
  legend("topright", c("Upregulated", "Downregulated", "Non-significant"),
         col = c("red", "blue", "gray70"), pch = 16, cex = 0.8)
}, "02_MA_Plot")

# 3. Expression Heatmap of Top DE Genes
safe_plot(function() {
  # Get top 50 DE genes
  top_de_genes <- head(results[order(results$P.Value), ], 50)
  
  # Extract expression data
  heatmap_data <- log_fpkm[top_de_genes$gene_id, ]
  
  # Scale by row (z-score)
  heatmap_data_scaled <- t(scale(t(heatmap_data)))
  
  # Set row names to gene symbols
  rownames(heatmap_data_scaled) <- top_de_genes$Symbol
  
  # Create annotation for samples
  col_annotation <- data.frame(
    Treatment = factor(ifelse(grepl("Ctrl", colnames(heatmap_data_scaled)), "Control", "Morphine"))
  )
  rownames(col_annotation) <- colnames(heatmap_data_scaled)
  
  # Create heatmap
  if(require(pheatmap, quietly = TRUE)) {
    pheatmap::pheatmap(heatmap_data_scaled,
             main = "Top 50 Differentially Expressed Genes",
             annotation_col = col_annotation,
             scale = "none",
             clustering_distance_rows = "correlation",
             clustering_distance_cols = "correlation",
             fontsize_row = 6,
             fontsize_col = 10,
             color = colorRampPalette(c("blue", "white", "red"))(100))
  } else {
    heatmap(heatmap_data_scaled, main = "Top 50 DE Genes", cexRow = 0.6, cexCol = 0.8)
  }
}, "03_Top_DE_Genes_Heatmap", width = 12, height = 16)

# ========================================================================
# COMPLEMENT-SPECIFIC ANALYSIS
# ========================================================================

cat("\n=== COMPLEMENT GENE DIFFERENTIAL EXPRESSION ANALYSIS ===\n")

# Load complement genes from processing script
if(file.exists("Outputs/01_Processing_QC/Data/complement_genes_with_expression.csv")) {
  complement_genes <- read.csv("Outputs/01_Processing_QC/Data/complement_genes_with_expression.csv")
  
  # Match with DE results
  complement_de <- results %>%
    filter(gene_id %in% complement_genes$id) %>%
    left_join(complement_genes[, c("id", "pathway_class")], by = c("gene_id" = "id")) %>%
    arrange(P.Value)
  
  cat("Found", nrow(complement_de), "complement genes in differential expression results\n")
  
  if(nrow(complement_de) > 0) {
    # Summary statistics
    complement_summary <- complement_de %>%
      summarise(
        total_genes = n(),
        significant_p05 = sum(P.Value < 0.05),
        significant_fdr05 = sum(adj.P.Val < 0.05),
        upregulated = sum(logFC > 0 & P.Value < 0.05),
        downregulated = sum(logFC < 0 & P.Value < 0.05),
        mean_logFC = mean(logFC),
        median_logFC = median(logFC)
      )
    
    print(complement_summary)
    
    # Pathway-specific analysis
    pathway_de_summary <- complement_de %>%
      group_by(pathway_class) %>%
      summarise(
        total_genes = n(),
        significant_genes = sum(P.Value < 0.05),
        mean_logFC = round(mean(logFC), 3),
        upregulated = sum(logFC > 0 & P.Value < 0.05),
        downregulated = sum(logFC < 0 & P.Value < 0.05),
        .groups = 'drop'
      ) %>%
      arrange(desc(significant_genes))
    
    cat("\nComplement pathway differential expression summary:\n")
    print(pathway_de_summary)
    
    # Save complement DE results
    write.csv(complement_de, "Outputs/02_Differential_Expression/Data/complement_DE_results.csv", row.names = FALSE)
    write.csv(pathway_de_summary, "Outputs/02_Differential_Expression/Data/complement_pathway_DE_summary.csv", row.names = FALSE)
    
    # Complement-specific volcano plot
    safe_plot(function() {
      plot(results$logFC, -log10(results$P.Value),
           xlab = "Log2 Fold Change (Morphine vs Control)",
           ylab = "-Log10 P-value",
           main = "Volcano Plot: Complement Genes Highlighted",
           pch = 16, cex = 0.6, col = "gray90")
      
      # Highlight complement genes
      points(complement_de$logFC, -log10(complement_de$P.Value),
             col = "red", pch = 16, cex = 1.2)
      
      # Add reference lines
      abline(h = -log10(0.05), col = "gray", lty = 2)
      abline(v = c(-0.5, 0.5), col = "gray", lty = 2)
      
      # Label significant complement genes
      sig_complement <- complement_de[complement_de$P.Value < 0.05, ]
      if(nrow(sig_complement) > 0) {
        text(sig_complement$logFC, -log10(sig_complement$P.Value), sig_complement$Symbol,
             pos = 3, cex = 0.8, col = "darkred", font = 2)
      }
      
      legend("topright", c("Complement genes", "Other genes"),
             col = c("red", "gray90"), pch = 16, cex = 0.8)
    }, "04_Complement_Volcano_Plot")
    
    # Complement heatmap
    if(nrow(complement_de) > 1) {
      safe_plot(function() {
        comp_expr_data <- log_fpkm[complement_de$gene_id, ]
        comp_expr_scaled <- t(scale(t(comp_expr_data)))
        rownames(comp_expr_scaled) <- complement_de$Symbol
        
        # Color by pathway class
        pathway_colors <- rainbow(length(unique(complement_de$pathway_class)))
        names(pathway_colors) <- unique(complement_de$pathway_class)
        
        row_annotation <- data.frame(
          Pathway = complement_de$pathway_class,
          Significant = ifelse(complement_de$P.Value < 0.05, "Yes", "No")
        )
        rownames(row_annotation) <- complement_de$Symbol
        
        col_annotation <- data.frame(
          Treatment = factor(ifelse(grepl "Ctrl", colnames(comp_expr_scaled)), "Control", "Morphine"))
        )
        rownames(col_annotation) <- colnames(comp_expr_scaled)
        
        if(require(pheatmap, quietly = TRUE)) {
          pheatmap::pheatmap(comp_expr_scaled,
                   main = "Complement Genes Differential Expression",
                   annotation_col = col_annotation,
                   annotation_row = row_annotation,
                   scale = "none",
                   clustering_distance_rows = "correlation",
                   clustering_distance_cols = "correlation",
                   fontsize_row = 8,
                   fontsize_col = 10,
                   color = colorRampPalette(c("blue", "white", "red"))(100))
        } else {
          heatmap(comp_expr_scaled, main = "Complement Genes DE", cexRow = 0.8, cexCol = 0.8)
        }
      }, "05_Complement_DE_Heatmap", width = 10, height = max(8, nrow(complement_de) * 0.3))
    }
  }
} else {
  cat("Complement gene data not found. Run processing script first.\n")
}

# ========================================================================
# SUMMARY STATISTICS
# ========================================================================

# Generate comprehensive summary
de_summary <- list(
  total_genes = nrow(results),
  significant_p05 = sum(results$P.Value < 0.05),
  significant_fdr05 = sum(results$adj.P.Val < 0.05),
  upregulated_p05 = sum(results$logFC > 0 & results$P.Value < 0.05),
  downregulated_p05 = sum(results$logFC < 0 & results$P.Value < 0.05),
  upregulated_fdr05 = sum(results$logFC > 0 & results$adj.P.Val < 0.05),
  downregulated_fdr05 = sum(results$logFC < 0 & results$adj.P.Val < 0.05),
  max_upregulation = max(results$logFC),
  max_downregulation = min(results$logFC),
  top_upregulated = results$Symbol[which.max(results$logFC)],
  top_downregulated = results$Symbol[which.min(results$logFC)]
)

# Save summary
write.csv(data.frame(Metric = names(de_summary), Value = unlist(de_summary)),
          "Outputs/02_Differential_Expression/Data/DE_summary_statistics.csv", row.names = FALSE)

# Print summary
cat("\n=== DIFFERENTIAL EXPRESSION SUMMARY ===\n")
cat("Total genes analyzed:", de_summary$total_genes, "\n")
cat("Significant genes (p < 0.05):", de_summary$significant_p05, "\n")
cat("Significant genes (FDR < 0.05):", de_summary$significant_fdr05, "\n")
cat("Upregulated (p < 0.05):", de_summary$upregulated_p05, "\n")
cat("Downregulated (p < 0.05):", de_summary$downregulated_p05, "\n")
cat("Most upregulated gene:", de_summary$top_upregulated, "(FC =", round(2^de_summary$max_upregulation, 2), ")\n")
cat("Most downregulated gene:", de_summary$top_downregulated, "(FC =", round(2^de_summary$max_downregulation, 2), ")\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Results saved to: Outputs/02_Differential_Expression/\n")
cat("- DE_results.csv: Complete differential expression results\n")
cat("- Plots/: Visualization files\n")
cat("- Data/: Summary statistics and complement-specific results\n")
