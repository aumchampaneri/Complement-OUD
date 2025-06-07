#===============================================================================
# DIFFERENTIAL EXPRESSION ANALYSIS
#===============================================================================
# 
# This script performs differential expression analysis using DESeq2 with 
# batch-corrected data from the preprocessing pipeline.
#
# WORKFLOW:
# 1. Load batch-corrected DESeq2 object
# 2. Define contrasts and perform DE analysis
# 3. Extract and filter results
# 4. Create visualization plots
# 5. Save results for pathway enrichment
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

# Create directories
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

#===============================================================================
# 1. LOAD PREPROCESSED DATA
#===============================================================================

# Load batch-corrected DESeq2 object
load(file.path(input_dir, "deseq2_for_DE_analysis.RData"))

cat("Loaded DESeq2 object with", ncol(dds_batch_norm), "samples and", 
    nrow(dds_batch_norm), "genes\n")
cat("Sample conditions:\n")
print(table(sample_info$condition))

#===============================================================================
# 2. DEFINE CONTRASTS AND PERFORM DE ANALYSIS
#===============================================================================

# Check the design and sample information
cat("Design formula:", as.character(design(dds_batch_norm)), "\n")
cat("Available factor levels:\n")
print(colData(dds_batch_norm))

# Check what columns are available in sample metadata
cat("Available columns in sample_info:\n")
print(colnames(sample_info))
cat("Unique values in key columns:\n")
if("condition" %in% colnames(sample_info)) {
  cat("Condition:", paste(unique(sample_info$condition), collapse = ", "), "\n")
}
if("sex" %in% colnames(sample_info)) {
  cat("Sex:", paste(unique(sample_info$sex), collapse = ", "), "\n")
}
if("region" %in% colnames(sample_info)) {
  cat("Region:", paste(unique(sample_info$region), collapse = ", "), "\n")
}

# Convert character columns to factors for DESeq2 design
colData(dds_batch_norm)$sex <- factor(colData(dds_batch_norm)$sex)
colData(dds_batch_norm)$region <- factor(colData(dds_batch_norm)$region)

# Check factor levels
cat("Factor levels after conversion:\n")
cat("Sex levels:", levels(colData(dds_batch_norm)$sex), "\n")
cat("Region levels:", levels(colData(dds_batch_norm)$region), "\n")

# Check sample distribution to avoid rank deficiency
cat("Sample distribution:\n")
print(table(colData(dds_batch_norm)$condition, colData(dds_batch_norm)$sex))
print(table(colData(dds_batch_norm)$condition, colData(dds_batch_norm)$region))
print(table(colData(dds_batch_norm)$batch, colData(dds_batch_norm)$condition))

# Check for perfect confounding
cat("Checking for confounding:\n")
print("Batch vs Sex:")
print(table(colData(dds_batch_norm)$batch, colData(dds_batch_norm)$sex))
print("Batch vs Region:")
print(table(colData(dds_batch_norm)$batch, colData(dds_batch_norm)$region))

# Start with the simplest design that works
# Try different designs in order of complexity
designs_to_try <- list(
  "condition_only" = ~ condition,
  "batch_condition" = ~ batch + condition,
  "condition_sex_region" = ~ condition + sex + region,
  "batch_condition_sex" = ~ batch + condition + sex,
  "batch_condition_region" = ~ batch + condition + region
)

dds_working <- NULL
working_design <- NULL

for (design_name in names(designs_to_try)) {
  cat("Trying design:", design_name, ":", as.character(designs_to_try[[design_name]]), "\n")
  
  tryCatch({
    design(dds_batch_norm) <- designs_to_try[[design_name]]
    dds_test <- DESeq(dds_batch_norm)
    cat("Success with design:", design_name, "\n")
    dds_working <- dds_test
    working_design <- design_name
    break
  }, error = function(e) {
    cat("Failed with design:", design_name, "Error:", e$message, "\n")
  })
}

if (is.null(dds_working)) {
  stop("No design formula worked. Check for perfect confounding in your data.")
}

dds_batch_norm <- dds_working
cat("Using design:", working_design, "\n")

# Check what coefficients are available
cat("Available coefficients:\n")
print(resultsNames(dds_batch_norm))

# Define contrasts based on the working design
# Adjust contrasts based on what design worked
contrasts_list <- list()

# Always include the main condition comparison
contrasts_list[["Pooled_OUD_vs_Control"]] <- list(
  contrast = c("condition", "OUD", "Control"),
  description = "Pooled OUD vs Control"
)

# Add subset-based contrasts that will work regardless of design
contrasts_list[["Male_OUD_vs_Control"]] <- list(
  contrast = c("condition", "OUD", "Control"),
  description = "Male OUD vs Control",
  subset = "sex == 'Male'"
)

contrasts_list[["Female_OUD_vs_Control"]] <- list(
  contrast = c("condition", "OUD", "Control"),
  description = "Female OUD vs Control", 
  subset = "sex == 'Female'"
)

contrasts_list[["DLPFC_OUD_vs_Control"]] <- list(
  contrast = c("condition", "OUD", "Control"),
  description = "DLPFC OUD vs Control",
  subset = "region == 'DLPFC'"
)

contrasts_list[["NAC_OUD_vs_Control"]] <- list(
  contrast = c("condition", "OUD", "Control"),
  description = "NAC OUD vs Control",
  subset = "region == 'NAC'"
)

# Add within-condition comparisons only if sex is in the design
if ("sex" %in% all.vars(design(dds_batch_norm))) {
  contrasts_list[["OUD_Male_vs_Female"]] <- list(
    contrast = c("sex", "Male", "Female"),
    description = "OUD Males vs OUD Females",
    subset = "condition == 'OUD'"
  )
  
  contrasts_list[["Control_Male_vs_Female"]] <- list(
    contrast = c("sex", "Male", "Female"),
    description = "Control Males vs Control Females",
    subset = "condition == 'Control'"
  )
} else {
  cat("Sex not in design formula - skipping sex-based comparisons within conditions\n")
}

# Add region comparisons only if region is in the design
if ("region" %in% all.vars(design(dds_batch_norm))) {
  contrasts_list[["OUD_DLPFC_vs_NAC"]] <- list(
    contrast = c("region", "NAC", "DLPFC"), 
    description = "OUD DLPFC vs OUD NAC",
    subset = "condition == 'OUD'"
  )
  
  contrasts_list[["Control_DLPFC_vs_NAC"]] <- list(
    contrast = c("region", "NAC", "DLPFC"),
    description = "Control DLPFC vs Control NAC",
    subset = "condition == 'Control'"
  )
} else {
  cat("Region not in design formula - skipping region-based comparisons within conditions\n")
}

cat("Final contrast list contains", length(contrasts_list), "contrasts\n")

#===============================================================================
# 3. PERFORM DE ANALYSIS FOR EACH CONTRAST
#===============================================================================

de_results <- list()

for (contrast_name in names(contrasts_list)) {
  cat("\n=== Analyzing contrast:", contrast_name, "===\n")
  cat("Description:", contrasts_list[[contrast_name]]$description, "\n")
  
  tryCatch({
    # Check if this contrast requires subsetting
    if ("subset" %in% names(contrasts_list[[contrast_name]])) {
      # Create subset DESeq object
      subset_condition <- contrasts_list[[contrast_name]]$subset
      cat("Subsetting data for:", subset_condition, "\n")
      
      # Evaluate the subset condition
      subset_samples <- with(as.data.frame(colData(dds_batch_norm)), eval(parse(text = subset_condition)))
      
      if (sum(subset_samples) < 4) {
        cat("Not enough samples after subsetting (", sum(subset_samples), ") - skipping\n")
        next
      }
      
      dds_subset <- dds_batch_norm[, subset_samples]
      
      # Check if the contrast variables are still meaningful after subsetting
      contrast_var <- contrasts_list[[contrast_name]]$contrast[1]
      if (contrast_var %in% colnames(colData(dds_subset))) {
        unique_levels <- length(unique(colData(dds_subset)[[contrast_var]]))
        if (unique_levels < 2) {
          cat("Only", unique_levels, "level(s) of", contrast_var, "after subsetting - skipping\n")
          next
        }
      }
      
      # Remove unused factor levels and re-run DESeq
      for (col in colnames(colData(dds_subset))) {
        if (is.factor(colData(dds_subset)[[col]])) {
          colData(dds_subset)[[col]] <- droplevels(colData(dds_subset)[[col]])
        }
      }
      
      dds_subset <- DESeq(dds_subset)
      
      # Extract results
      res <- results(dds_subset, 
                     contrast = contrasts_list[[contrast_name]]$contrast,
                     alpha = 0.05)
    } else {
      # Use full dataset
      res <- results(dds_batch_norm, 
                     contrast = contrasts_list[[contrast_name]]$contrast,
                     alpha = 0.05)
    }
    
    # Convert to data frame and add gene names
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    res_df <- res_df[!is.na(res_df$padj), ]
    
    # Add significance classification
    res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 0.5
    res_df$direction <- ifelse(res_df$log2FoldChange > 0, "Up", "Down")
    
    # Store results
    de_results[[contrast_name]] <- res_df
    
    # Print summary
    cat("Total genes tested:", nrow(res_df), "\n")
    cat("Significant genes (padj < 0.05, |logFC| > 0.5):", sum(res_df$significant), "\n")
    cat("Upregulated:", sum(res_df$significant & res_df$direction == "Up"), "\n")
    cat("Downregulated:", sum(res_df$significant & res_df$direction == "Down"), "\n")
    
  }, error = function(e) {
    cat("Error analyzing contrast", contrast_name, ":", e$message, "\n")
    if (exists("dds_subset")) {
      cat("Available result names in subset:", resultsNames(dds_subset), "\n")
    } else {
      cat("Available result names:", resultsNames(dds_batch_norm), "\n")
    }
    
    # Skip this contrast and continue
    de_results[[contrast_name]] <- NULL
  })
}

#===============================================================================
# 4. CREATE VISUALIZATION PLOTS
#===============================================================================

# Function to create volcano plot
create_volcano_plot <- function(res_df, contrast_name, output_dir) {
  # Prepare data for plotting
  plot_df <- res_df
  plot_df$log_padj <- -log10(plot_df$padj)
  plot_df$color <- "Not Significant"
  plot_df$color[plot_df$significant & plot_df$direction == "Up"] <- "Upregulated"
  plot_df$color[plot_df$significant & plot_df$direction == "Down"] <- "Downregulated"
  
  # Create volcano plot
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
  
  # Save plot
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
for (contrast_name in names(de_results)) {
  if (!is.null(de_results[[contrast_name]])) {
    create_volcano_plot(de_results[[contrast_name]], contrast_name, plots_dir)
    create_ma_plot(de_results[[contrast_name]], contrast_name, plots_dir)
    cat("Created plots for", contrast_name, "\n")
  }
}

#===============================================================================
# 5. SAVE RESULTS
#===============================================================================

# Save individual contrast results
for (contrast_name in names(de_results)) {
  if (!is.null(de_results[[contrast_name]])) {
    
    # Save full results
    write.csv(de_results[[contrast_name]], 
              file.path(output_dir, paste0(contrast_name, "_full_results.csv")), 
              row.names = FALSE)
    
    # Save significant genes only
    sig_genes <- de_results[[contrast_name]][de_results[[contrast_name]]$significant, ]
    write.csv(sig_genes, 
              file.path(output_dir, paste0(contrast_name, "_significant_genes.csv")), 
              row.names = FALSE)
    
    cat("Saved results for", contrast_name, "\n")
  }
}

# Create comprehensive summary table
summary_table <- data.frame()
for (contrast_name in names(de_results)) {
  if (!is.null(de_results[[contrast_name]])) {
    res_df <- de_results[[contrast_name]]
    summary_row <- data.frame(
      Contrast = contrast_name,
      Description = contrasts_list[[contrast_name]]$description,
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
    summary_table <- rbind(summary_table, summary_row)
  }
}

# Save summary
write.csv(summary_table, file.path(output_dir, "DE_analysis_summary.csv"), 
          row.names = FALSE)

# Print final summary
cat("\n=== FINAL SUMMARY ===\n")
print(summary_table)

#===============================================================================
# 6. CREATE HEATMAPS FOR TOP DE GENES
#===============================================================================

# Function to create heatmap of top DE genes
create_heatmap <- function(contrast_name, res_df, dds_object, output_dir, top_n = 50) {
  if (sum(res_df$significant) < 5) {
    cat("Not enough significant genes for heatmap in", contrast_name, "\n")
    return(NULL)
  }
  
  # Get top significant genes
  sig_genes <- res_df[res_df$significant, ]
  sig_genes <- sig_genes[order(sig_genes$padj), ]
  top_genes <- head(sig_genes, top_n)
  
  if (nrow(top_genes) == 0) {
    cat("No significant genes found for", contrast_name, "\n")
    return(NULL)
  }
  
  # Get variance stabilized data
  vsd <- vst(dds_object, blind = FALSE)
  mat <- assay(vsd)[rownames(top_genes), ]
  
  # Z-score normalization
  mat <- t(scale(t(mat)))
  
  # Create annotation
  annotation_col <- data.frame(
    Condition = colData(dds_object)$condition,
    Sex = colData(dds_object)$sex,
    Region = colData(dds_object)$region,
    row.names = colnames(mat)
  )
  
  # Create heatmap
  png(file.path(output_dir, paste0("heatmap_top_", nrow(top_genes), "_", contrast_name, ".png")), 
      width = 12, height = 10, units = "in", res = 300)
  
  pheatmap(mat,
           annotation_col = annotation_col,
           scale = "none",  # Already scaled
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "complete",
           color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100),
           main = paste("Top", nrow(top_genes), "DE Genes:", gsub("_", " ", contrast_name)),
           fontsize = 8,
           fontsize_row = 6,
           show_rownames = TRUE,
           show_colnames = TRUE)
  
  dev.off()
  
  cat("Created heatmap for", contrast_name, "with", nrow(top_genes), "genes\n")
}

# Generate heatmaps for each contrast
for (contrast_name in names(de_results)) {
  if (!is.null(de_results[[contrast_name]])) {
    create_heatmap(contrast_name, de_results[[contrast_name]], dds_batch_norm, plots_dir)
  }
}

#===============================================================================
# 7. CREATE GENE LISTS FOR PATHWAY ANALYSIS
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
  
  # Save ranked list for GSEA (all genes ranked by stat)
  ranked_genes <- res_df[order(-res_df$stat), c("gene", "stat")]
  write.table(ranked_genes, 
              file.path(output_dir, paste0(contrast_name, "_ranked_genes.rnk")), 
              row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  
  cat("Saved gene lists for", contrast_name, "\n")
}

# Generate gene lists for each contrast
for (contrast_name in names(de_results)) {
  if (!is.null(de_results[[contrast_name]])) {
    save_gene_lists(contrast_name, de_results[[contrast_name]], genelist_dir)
  }
}

#===============================================================================
# 8. SAVE WORKSPACE AND OBJECTS
#===============================================================================

# Save all DE results and objects
save(dds_batch_norm, de_results, contrasts_list, summary_table,
     file = file.path(output_dir, "differential_expression_complete.RData"))

# Create a master gene list combining all contrasts
all_sig_genes <- unique(unlist(lapply(de_results, function(x) {
  if (!is.null(x)) {
    x$gene[x$significant]
  }
})))

write.table(all_sig_genes, 
            file.path(genelist_dir, "all_contrasts_significant_genes.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

#===============================================================================
# 9. FINAL REPORT
#===============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("DIFFERENTIAL EXPRESSION ANALYSIS COMPLETED\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

cat("\nCONTRASTS ANALYZED:\n")
for (i in 1:nrow(summary_table)) {
  cat(sprintf("%d. %s: %d significant genes (%d up, %d down)\n", 
              i, 
              summary_table$Contrast[i],
              summary_table$Significant_Genes[i],
              summary_table$Upregulated[i],
              summary_table$Downregulated[i]))
}

cat("\nOUTPUT FILES GENERATED:\n")
cat("• Full results CSV files for each contrast\n")
cat("• Significant genes CSV files for each contrast\n")
cat("• Volcano plots for each contrast\n")
cat("• MA plots for each contrast\n")
cat("• Heatmaps of top DE genes for each contrast\n")
cat("• Gene lists for pathway analysis (.txt and .rnk files)\n")
cat("• Summary table of all contrasts\n")
cat("• Complete R workspace with all objects\n")

cat("\nFILES SAVED TO:\n")
cat("• Results:", output_dir, "\n")
cat("• Plots:", plots_dir, "\n")
cat("• Gene lists:", genelist_dir, "\n")

cat("\nTOTAL UNIQUE SIGNIFICANT GENES ACROSS ALL CONTRASTS:", length(all_sig_genes), "\n")

cat("\nREADY FOR DOWNSTREAM ANALYSIS:\n")
cat("• Pathway enrichment analysis\n")
cat("• Gene set enrichment analysis (GSEA)\n")
cat("• Functional annotation\n")
cat("• Network analysis\n")

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("ANALYSIS COMPLETE - CHECK OUTPUT DIRECTORIES FOR RESULTS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Print session info for reproducibility
cat("SESSION INFO:\n")
sessionInfo()