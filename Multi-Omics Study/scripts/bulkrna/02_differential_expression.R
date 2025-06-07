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
# 2. DIFFERENTIAL EXPRESSION ANALYSIS
#===============================================================================

# Function to perform DE analysis and save results
perform_DE_analysis <- function(dds, contrast_name, contrast_vector, output_suffix) {
  cat("Performing DE analysis:", contrast_name, "\n")
  
  # Extract results
  res <- results(dds, contrast = contrast_vector, alpha = 0.05)
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  # Add statistics
  res_df <- res_df %>%
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    mutate(
      significant = padj < 0.05 & abs(log2FoldChange) > 1,
      direction = case_when(
        log2FoldChange > 1 & padj < 0.05 ~ "Upregulated",
        log2FoldChange < -1 & padj < 0.05 ~ "Downregulated",
        TRUE ~ "Not significant"
      )
    )
  
  # Summary statistics
  cat("Results for", contrast_name, ":\n")
  cat("Total genes tested:", nrow(res_df), "\n")
  cat("Significant genes (padj < 0.05, |log2FC| > 1):", sum(res_df$significant), "\n")
  cat("Upregulated:", sum(res_df$direction == "Upregulated"), "\n")
  cat("Downregulated:", sum(res_df$direction == "Downregulated"), "\n\n")
  
  # Save results
  write.table(res_df, file.path(output_dir, paste0("DE_results_", output_suffix, ".txt")), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  return(res_df)
}

# Define all contrasts matching snRNA-seq analysis
contrasts_list <- list(
  "Pooled_OUD_vs_Control" = c("condition", "OUD", "Control"),
  "OUD_vs_Control_DLPFC" = c("condition", "OUD", "Control"),  # Will subset to DLPFC
  "OUD_vs_Control_NAC" = c("condition", "OUD", "Control"),    # Will subset to NAC
  "OUD_vs_Control_Male" = c("condition", "OUD", "Control"),   # Will subset to Male
  "OUD_vs_Control_Female" = c("condition", "OUD", "Control")  # Will subset to Female
)

# Store all results
all_results <- list()

# 1. Pooled analysis (all samples)
all_results[["Pooled_OUD_vs_Control"]] <- perform_DE_analysis(
  dds_batch_norm, "Pooled OUD vs Control", 
  contrasts_list[["Pooled_OUD_vs_Control"]], "Pooled_OUD_vs_Control"
)

# 2. Region-specific analyses
# DLPFC samples
dlpfc_samples <- rownames(sample_info)[sample_info$region == "DLPFC"]
dds_dlpfc <- dds_batch_norm[, dlpfc_samples]
dds_dlpfc$condition <- droplevels(dds_dlpfc$condition)

all_results[["OUD_vs_Control_DLPFC"]] <- perform_DE_analysis(
  dds_dlpfc, "OUD vs Control in DLPFC", 
  contrasts_list[["OUD_vs_Control_DLPFC"]], "OUD_vs_Control_DLPFC"
)

# NAC samples
nac_samples <- rownames(sample_info)[sample_info$region == "NAC"]
dds_nac <- dds_batch_norm[, nac_samples]
dds_nac$condition <- droplevels(dds_nac$condition)

all_results[["OUD_vs_Control_NAC"]] <- perform_DE_analysis(
  dds_nac, "OUD vs Control in NAC", 
  contrasts_list[["OUD_vs_Control_NAC"]], "OUD_vs_Control_NAC"
)

# 3. Sex-specific analyses
# Male samples
male_samples <- rownames(sample_info)[sample_info$sex == "Male"]
dds_male <- dds_batch_norm[, male_samples]
dds_male$condition <- droplevels(dds_male$condition)

all_results[["OUD_vs_Control_Male"]] <- perform_DE_analysis(
  dds_male, "OUD vs Control in Males", 
  contrasts_list[["OUD_vs_Control_Male"]], "OUD_vs_Control_Male"
)

# Female samples
female_samples <- rownames(sample_info)[sample_info$sex == "Female"]
dds_female <- dds_batch_norm[, female_samples]
dds_female$condition <- droplevels(dds_female$condition)

all_results[["OUD_vs_Control_Female"]] <- perform_DE_analysis(
  dds_female, "OUD vs Control in Females", 
  contrasts_list[["OUD_vs_Control_Female"]], "OUD_vs_Control_Female"
)

# 4. Regional effect analyses (OUD samples only)
oud_samples <- rownames(sample_info)[sample_info$condition == "OUD"]
sample_info_oud <- sample_info[oud_samples, ]

# Create DESeq2 object for region comparison in OUD samples
if (length(unique(sample_info_oud$region)) > 1) {
  dds_oud_region <- DESeqDataSetFromMatrix(
    countData = counts(dds_batch_norm)[, oud_samples],
    colData = sample_info_oud,
    design = ~ region
  )
  dds_oud_region <- DESeq(dds_oud_region)
  
  all_results[["OUD_Effect_DLPFC_vs_NAC"]] <- perform_DE_analysis(
    dds_oud_region, "OUD Effect: DLPFC vs NAC", 
    c("region", "DLPFC", "NAC"), "OUD_Effect_DLPFC_vs_NAC"
  )
}

# 5. Sex effect analyses (OUD samples only)
if (length(unique(sample_info_oud$sex)) > 1) {
  dds_oud_sex <- DESeqDataSetFromMatrix(
    countData = counts(dds_batch_norm)[, oud_samples],
    colData = sample_info_oud,
    design = ~ sex
  )
  dds_oud_sex <- DESeq(dds_oud_sex)
  
  all_results[["OUD_Effect_Male_vs_Female"]] <- perform_DE_analysis(
    dds_oud_sex, "OUD Effect: Male vs Female", 
    c("sex", "Male", "Female"), "OUD_Effect_Male_vs_Female"
  )
}

#===============================================================================
# 3. VISUALIZATION (for all contrasts)
#===============================================================================

# Create plots for all contrasts
for (contrast_name in names(all_results)) {
  results_df <- all_results[[contrast_name]]
  
  # Volcano plot
  p_volcano <- EnhancedVolcano(results_df,
                             lab = results_df$gene,
                             x = 'log2FoldChange',
                             y = 'padj',
                             title = paste('Volcano Plot:', gsub("_", " ", contrast_name)),
                             subtitle = 'Bulk RNA-seq Differential Expression',
                             pCutoff = 0.05,
                             FCcutoff = 1,
                             pointSize = 2,
                             labSize = 3,
                             colAlpha = 0.7,
                             legendPosition = 'right',
                             drawConnectors = TRUE,
                             widthConnectors = 0.3,
                             max.overlaps = 20)
  
  ggsave(file.path(plots_dir, paste0("volcano_plot_", contrast_name, ".png")), 
         p_volcano, width = 12, height = 8)
  
  # MA plot
  p_ma <- ggplot(results_df, aes(x = baseMean, y = log2FoldChange)) +
    geom_point(aes(color = direction), alpha = 0.6) +
    scale_color_manual(values = c("Upregulated" = "red", 
                                 "Downregulated" = "blue", 
                                 "Not significant" = "grey")) +
    scale_x_log10() +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
    labs(title = paste("MA Plot:", gsub("_", " ", contrast_name)),
         x = "Mean Expression (log10)",
         y = "Log2 Fold Change",
         color = "Regulation") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12))
  
  ggsave(file.path(plots_dir, paste0("ma_plot_", contrast_name, ".png")), 
         p_ma, width = 10, height = 6)
  
  cat("Created plots for:", contrast_name, "\n")
}

# Heatmap of top differentially expressed genes (for main pooled contrast)
main_results <- all_results[["Pooled_OUD_vs_Control"]]
top_genes <- main_results %>% 
  filter(significant) %>% 
  arrange(padj) %>% 
  head(50) %>% 
  pull(gene)

if (length(top_genes) > 0) {
  # Get VSD-transformed counts for heatmap
  vsd_counts <- assay(vst(dds_batch_norm, blind = FALSE))
  heatmap_data <- vsd_counts[top_genes, ]
  
  # Create annotation for samples
  sample_annotation <- data.frame(
    Condition = sample_info$condition,
    row.names = colnames(heatmap_data)
  )
  
  # Create heatmap
  png(file.path(plots_dir, "heatmap_top_DE_genes_Pooled.png"), width = 1000, height = 800)
  pheatmap(heatmap_data,
           scale = "row",
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           annotation_col = sample_annotation,
           show_rownames = TRUE,
           show_colnames = FALSE,
           main = "Top 50 Differentially Expressed Genes (Pooled Analysis)")
  dev.off()
}

#===============================================================================
# 4. SAVE RESULTS FOR PATHWAY ENRICHMENT
#===============================================================================

# Save gene lists for each contrast
for (contrast_name in names(all_results)) {
  results_df <- all_results[[contrast_name]]
  
  # Upregulated genes
  upregulated_genes <- results_df %>% 
    filter(direction == "Upregulated") %>% 
    pull(gene)
  
  # Downregulated genes
  downregulated_genes <- results_df %>% 
    filter(direction == "Downregulated") %>% 
    pull(gene)
  
  # Significant genes
  significant_genes <- results_df %>% filter(significant)
  
  # Save files
  writeLines(upregulated_genes, 
             file.path(output_dir, paste0("upregulated_genes_", contrast_name, ".txt")))
  writeLines(downregulated_genes, 
             file.path(output_dir, paste0("downregulated_genes_", contrast_name, ".txt")))
  write.table(significant_genes, 
              file.path(output_dir, paste0("significant_genes_", contrast_name, ".txt")), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Ranked gene list for GSEA
  ranked_genes <- results_df %>%
    filter(!is.na(log2FoldChange)) %>%
    arrange(desc(log2FoldChange)) %>%
    dplyr::select(gene, log2FoldChange)
  
  write.table(ranked_genes, 
              file.path(output_dir, paste0("ranked_gene_list_", contrast_name, ".txt")), 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# Save all results objects
save(dds_batch_norm, all_results, sample_info, 
     file = file.path(output_dir, "DE_analysis_all_contrasts.RData"))

cat("\nDifferential expression analysis completed for all contrasts!\n")
cat("Contrasts analyzed:\n")
for (name in names(all_results)) {
  cat("-", name, "\n")
}
cat("\nOutputs saved for each contrast:\n")
cat("- Complete DE results table\n")
cat("- Significant genes table\n") 
cat("- Gene lists for pathway enrichment\n")
cat("- Ranked gene lists for GSEA\n")
cat("- Ready for pathway enrichment analysis (script 05)!\n")
