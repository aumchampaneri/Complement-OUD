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
output_dir <- "data/processed/bulkrna"
plots_dir <- "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study/results/bulkrna"

#===============================================================================
# 1. LOAD PREPROCESSED DATA
#===============================================================================

# Load batch-corrected DESeq2 object
load(file.path(output_dir, "deseq2_for_DE_analysis.RData"))

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
# 3. VISUALIZATION (for main contrast)
#===============================================================================

# Create plots for the main pooled analysis
main_results <- all_results[["Pooled_OUD_vs_Control"]]

# Volcano plot
p_volcano <- EnhancedVolcano(main_results,
                           lab = main_results$gene,
                           x = 'log2FoldChange',
                           y = 'padj',
                           title = 'Pooled: OUD vs Control',
                           subtitle = 'Bulk RNA-seq Differential Expression',
                           pCutoff = 0.05,
                           FCcutoff = 1,
                           pointSize = 2,
                           labSize = 3,
                           colAlpha = 0.7,
                           legendPosition = 'right',
                           drawConnectors = TRUE,
                           widthConnectors = 0.3)

ggsave(file.path(plots_dir, "volcano_plot_Pooled_OUD_vs_Control.png"), 
       p_volcano, width = 12, height = 8)

# MA plot
p_ma <- ggplot(main_results, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = direction), alpha = 0.6) +
  scale_color_manual(values = c("Upregulated" = "red", 
                               "Downregulated" = "blue", 
                               "Not significant" = "grey")) +
  scale_x_log10() +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
  labs(title = "MA Plot: Pooled OUD vs Control",
       x = "Mean Expression (log10)",
       y = "Log2 Fold Change",
       color = "Regulation") +
  theme_minimal()

ggsave(file.path(plots_dir, "ma_plot_Pooled_OUD_vs_Control.png"), 
       p_ma, width = 10, height = 6)

# Heatmap of top differentially expressed genes
top_genes <- main_results %>% 
  filter(significant) %>% 
  arrange(padj) %>% 
  head(50) %>% 
  pull(gene)