#===============================================================================
# WITHIN-CONDITION DIFFERENTIAL EXPRESSION ANALYSIS
#===============================================================================
# 
# This script performs differential expression analysis for within-condition
# comparisons (sex and region differences within OUD and Control groups).
# These analyses require subset-based approaches with simplified designs.
#
# WORKFLOW:
# 1. Load preprocessed DESeq2 object
# 2. Define subset-based contrasts for within-condition comparisons
# 3. Perform DE analysis on subsetted data
# 4. Create visualization plots
# 5. Save results to same directory as main DE analysis
#
#===============================================================================

# Load required libraries
library(DESeq2)
library(ggplot2)
library(dplyr)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)

# Set paths
setwd("/Users/aumchampaneri/Complement-OUD/Multi-Omics Study")
input_dir <- "data/processed/bulkrna/preprocessing"
output_dir <- "data/processed/bulkrna/differential_expression"
plots_dir <- "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study/results/bulkrna/differential_expression"

# Create directories (they should already exist from main script)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

#===============================================================================
# 1. LOAD PREPROCESSED DATA
#===============================================================================

# Load batch-corrected DESeq2 object
load(file.path(input_dir, "deseq2_for_DE_analysis.RData"))

cat("Loaded DESeq2 object with", ncol(dds_batch_norm), "samples and", 
    nrow(dds_batch_norm), "genes for within-condition analysis\n")

# Convert character columns to factors
colData(dds_batch_norm)$sex <- factor(colData(dds_batch_norm)$sex)
colData(dds_batch_norm)$region <- factor(colData(dds_batch_norm)$region)

cat("Sample distribution:\n")
print(table(colData(dds_batch_norm)$condition, colData(dds_batch_norm)$sex))
print(table(colData(dds_batch_norm)$condition, colData(dds_batch_norm)$region))

#===============================================================================
# 2. DEFINE WITHIN-CONDITION CONTRASTS
#===============================================================================

# Define the 4 within-condition contrasts
within_contrasts <- list(
  "OUD_Male_vs_Female" = list(
    contrast_var = "sex",
    contrast_levels = c("Male", "Female"),
    subset_condition = "condition == 'OUD'",
    description = "OUD Males vs OUD Females"
  ),
  
  "Control_Male_vs_Female" = list(
    contrast_var = "sex", 
    contrast_levels = c("Male", "Female"),
    subset_condition = "condition == 'Control'",
    description = "Control Males vs Control Females"
  ),
  
  "OUD_DLPFC_vs_NAC" = list(
    contrast_var = "region",
    contrast_levels = c("DLPFC", "NAC"),
    subset_condition = "condition == 'OUD'",
    description = "OUD DLPFC vs OUD NAC"
  ),
  
  "Control_DLPFC_vs_NAC" = list(
    contrast_var = "region",
    contrast_levels = c("DLPFC", "NAC"), 
    subset_condition = "condition == 'Control'",
    description = "Control DLPFC vs Control NAC"
  )
)

cat("Analyzing", length(within_contrasts), "within-condition contrasts\n")

#===============================================================================
# 3. PERFORM SUBSET-BASED DE ANALYSIS
#===============================================================================

de_results_within <- list()

for (contrast_name in names(within_contrasts)) {
  cat("\n=== Analyzing contrast:", contrast_name, "===\n")
  cat("Description:", within_contrasts[[contrast_name]]$description, "\n")
  
  tryCatch({
    # Get subset information
    subset_condition <- within_contrasts[[contrast_name]]$subset_condition
    contrast_var <- within_contrasts[[contrast_name]]$contrast_var
    contrast_levels <- within_contrasts[[contrast_name]]$contrast_levels
    
    cat("Subsetting data for:", subset_condition, "\n")
    cat("Comparing:", contrast_levels[1], "vs", contrast_levels[2], "in", contrast_var, "\n")
    
    # Create subset
    subset_samples <- with(as.data.frame(colData(dds_batch_norm)), eval(parse(text = subset_condition)))
    cat("Subset contains", sum(subset_samples), "samples\n")
    
    if (sum(subset_samples) < 6) {
      cat("Warning: Very few samples (", sum(subset_samples), ") - results may be unreliable\n")
    }
    
    # Create subset DESeq object
    dds_subset <- dds_batch_norm[, subset_samples]
    
    # Check distribution in subset
    cat("Distribution in subset:\n")
    if (contrast_var == "sex") {
      print(table(colData(dds_subset)$sex))
    } else if (contrast_var == "region") {
      print(table(colData(dds_subset)$region))
    }
    
    # Remove unused factor levels
    for (col in colnames(colData(dds_subset))) {
      if (is.factor(colData(dds_subset)[[col]])) {
        colData(dds_subset)[[col]] <- droplevels(colData(dds_subset)[[col]])
      }
    }
    
    # Try different designs for the subset
    designs_to_try <- list()
    designs_to_try[[paste0(contrast_var, "_only")]] <- as.formula(paste("~", contrast_var))
    designs_to_try[[paste0("batch_", contrast_var)]] <- as.formula(paste("~ batch +", contrast_var))
    
    dds_working <- NULL
    working_design <- NULL
    
    for (design_name in names(designs_to_try)) {
      cat("Trying design:", design_name, ":", as.character(designs_to_try[[design_name]]), "\n")
      
      tryCatch({
        design(dds_subset) <- designs_to_try[[design_name]]
        dds_test <- DESeq(dds_subset)
        cat("Success with design:", design_name, "\n")
        dds_working <- dds_test
        working_design <- design_name
        break
      }, error = function(e) {
        cat("Failed with design:", design_name, "Error:", e$message, "\n")
      })
    }
    
    if (is.null(dds_working)) {
      cat("No design worked for", contrast_name, "- skipping\n")
      next
    }
    
    cat("Using design:", working_design, "\n")
    cat("Available coefficients:", resultsNames(dds_working), "\n")
    
    # Extract results
    res <- results(dds_working, 
                   contrast = c(contrast_var, contrast_levels[1], contrast_levels[2]),
                   alpha = 0.05)
    
    # Convert to data frame and add gene names
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    res_df <- res_df[!is.na(res_df$padj), ]
    
    # Add significance classification
    res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 0.5
    res_df$direction <- ifelse(res_df$log2FoldChange > 0, "Up", "Down")
    
    # Store results
    de_results_within[[contrast_name]] <- res_df
    
    # Print summary
    cat("Total genes tested:", nrow(res_df), "\n")
    cat("Significant genes (padj < 0.05, |logFC| > 0.5):", sum(res_df$significant), "\n")
    cat("Upregulated:", sum(res_df$significant & res_df$direction == "Up"), "\n")
    cat("Downregulated:", sum(res_df$significant & res_df$direction == "Down"), "\n")
    
  }, error = function(e) {
    cat("Error analyzing contrast", contrast_name, ":", e$message, "\n")
    de_results_within[[contrast_name]] <- NULL
  })
}

#===============================================================================
# 4. CREATE VISUALIZATION PLOTS
#===============================================================================

# Function to create volcano plot
create_volcano_plot <- function(res_df, contrast_name, output_dir) {
  plot_df <- res_df
  plot_df$log_padj <- -log10(plot_df$padj)
  plot_df$color <- "Not Significant"
  plot_df$color[plot_df$significant & plot_df$direction == "Up"] <- "Upregulated"
  plot_df$color[plot_df$significant & plot_df$direction == "Down"] <- "Downregulated"
  
  p <- ggplot(plot_df, aes(x = log2FoldChange, y = log_padj, color = color)) +
    geom_point(alpha = 0.6, size = 0.8) +
    scale_color_manual(values = c("Not Significant" = "grey", 
                                 "Upregulated" = "red", 
                                 "Downregulated" = "blue")) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
    labs(title = paste("Volcano Plot:", gsub("_", " ", contrast_name)),
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-value",
         color = "Regulation") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          legend.position = "bottom")
  
  ggsave(file.path(output_dir, paste0("volcano_plot_", contrast_name, ".png")), 
         plot = p, width = 10, height = 8, dpi = 300)
  
  return(p)
}

# Function to create MA plot
create_ma_plot <- function(res_df, contrast_name, output_dir) {
  plot_df <- res_df
  plot_df$color <- "Not Significant"
  plot_df$color[plot_df$significant & plot_df$direction == "Up"] <- "Upregulated"
  plot_df$color[plot_df$significant & plot_df$direction == "Down"] <- "Downregulated"
  
  p <- ggplot(plot_df, aes(x = baseMean, y = log2FoldChange, color = color)) +
    geom_point(alpha = 0.6, size = 0.8) +
    scale_x_log10() +
    scale_color_manual(values = c("Not Significant" = "grey", 
                                 "Upregulated" = "red", 
                                 "Downregulated" = "blue")) +
    geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "solid", alpha = 0.8) +
    labs(title = paste("MA Plot:", gsub("_", " ", contrast_name)),
         x = "Mean of Normalized Counts",
         y = "Log2 Fold Change",
         color = "Regulation") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          legend.position = "bottom")
  
  ggsave(file.path(output_dir, paste0("ma_plot_", contrast_name, ".png")), 
         plot = p, width = 10, height = 8, dpi = 300)
  
  return(p)
}

# Generate plots for each contrast
for (contrast_name in names(de_results_within)) {
  if (!is.null(de_results_within[[contrast_name]])) {
    create_volcano_plot(de_results_within[[contrast_name]], contrast_name, plots_dir)
    create_ma_plot(de_results_within[[contrast_name]], contrast_name, plots_dir)
    cat("Created plots for", contrast_name, "\n")
  }
}

#===============================================================================
# 5. SAVE RESULTS
#===============================================================================

# Save individual contrast results
for (contrast_name in names(de_results_within)) {
  if (!is.null(de_results_within[[contrast_name]])) {
    
    # Save full results
    write.csv(de_results_within[[contrast_name]], 
              file.path(output_dir, paste0(contrast_name, "_full_results.csv")), 
              row.names = FALSE)
    
    # Save significant genes only
    sig_genes <- de_results_within[[contrast_name]][de_results_within[[contrast_name]]$significant, ]
    write.csv(sig_genes, 
              file.path(output_dir, paste0(contrast_name, "_significant_genes.csv")), 
              row.names = FALSE)
    
    cat("Saved results for", contrast_name, "\n")
  }
}

# Create summary table for within-condition contrasts
summary_table_within <- data.frame()
for (contrast_name in names(de_results_within)) {
  if (!is.null(de_results_within[[contrast_name]])) {
    res_df <- de_results_within[[contrast_name]]
    summary_row <- data.frame(
      Contrast = contrast_name,
      Description = within_contrasts[[contrast_name]]$description,
      Total_Genes = nrow(res_df),
      Significant_Genes = sum(res_df$significant),
      Upregulated = sum(res_df$significant & res_df$direction == "Up"),
      Downregulated = sum(res_df$significant & res_df$direction == "Down"),
      Percent_Significant = round(100 * sum(res_df$significant) / nrow(res_df), 2),
      Max_LogFC_Up = ifelse(any(res_df$direction == "Up"), 
                           round(max(res_df$log2FoldChange[res_df$direction == "Up"], na.rm = TRUE), 2), 0),
      Min_LogFC_Down = ifelse(any(res_df$direction == "Down"), 
                             round(min(res_df$log2FoldChange[res_df$direction == "Down"], na.rm = TRUE), 2), 0),
      Min_Padj = ifelse(nrow(res_df) > 0, 
                       format(min(res_df$padj, na.rm = TRUE), scientific = TRUE, digits = 3), "NA")
    )
    summary_table_within <- rbind(summary_table_within, summary_row)
  }
}

# Save within-condition summary
write.csv(summary_table_within, file.path(output_dir, "DE_analysis_within_condition_summary.csv"), 
          row.names = FALSE)

#===============================================================================
# 6. CREATE GENE LISTS FOR PATHWAY ANALYSIS
#===============================================================================

# Create directory for gene lists
genelist_dir <- file.path(output_dir, "gene_lists")
dir.create(genelist_dir, recursive = TRUE, showWarnings = FALSE)

# Function to save gene lists
save_gene_lists <- function(contrast_name, res_df, output_dir) {
  sig_genes <- res_df[res_df$significant, ]
  up_genes <- sig_genes[sig_genes$direction == "Up", ]
  down_genes <- sig_genes[sig_genes$direction == "Down", ]
  
  # Save all significant genes
  if (nrow(sig_genes) > 0) {
    write.table(sig_genes$gene, 
                file.path(output_dir, paste0(contrast_name, "_all_significant.txt")), 
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  # Save upregulated genes
  if (nrow(up_genes) > 0) {
    write.table(up_genes$gene, 
                file.path(output_dir, paste0(contrast_name, "_upregulated.txt")), 
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  # Save downregulated genes
  if (nrow(down_genes) > 0) {
    write.table(down_genes$gene, 
                file.path(output_dir, paste0(contrast_name, "_downregulated.txt")), 
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  # Save ranked list for GSEA
  ranked_genes <- res_df[order(-res_df$stat), c("gene", "stat")]
  write.table(ranked_genes, 
              file.path(output_dir, paste0(contrast_name, "_ranked_genes.rnk")), 
              row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  
  cat("Saved gene lists for", contrast_name, "\n")
}

# Generate gene lists for each contrast
for (contrast_name in names(de_results_within)) {
  if (!is.null(de_results_within[[contrast_name]])) {
    save_gene_lists(contrast_name, de_results_within[[contrast_name]], genelist_dir)
  }
}

#===============================================================================
# 7. SAVE WORKSPACE AND FINAL REPORT
#===============================================================================

# Save within-condition results
save(de_results_within, within_contrasts, summary_table_within,
     file = file.path(output_dir, "within_condition_differential_expression.RData"))

# Print final summary
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("WITHIN-CONDITION DIFFERENTIAL EXPRESSION ANALYSIS COMPLETED\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

cat("\nWITHIN-CONDITION CONTRASTS ANALYZED:\n")
for (i in 1:nrow(summary_table_within)) {
  cat(sprintf("%d. %s: %d significant genes (%d up, %d down)\n", 
              i, 
              summary_table_within$Contrast[i],
              summary_table_within$Significant_Genes[i],
              summary_table_within$Upregulated[i],
              summary_table_within$Downregulated[i]))
}

print(summary_table_within)

cat("\nFILES SAVED TO SAME DIRECTORIES AS MAIN ANALYSIS:\n")
cat("• Results:", output_dir, "\n")
cat("• Plots:", plots_dir, "\n")
cat("• Gene lists:", genelist_dir, "\n")

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("WITHIN-CONDITION ANALYSIS COMPLETE\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Print session info
cat("SESSION INFO:\n")
sessionInfo()
