# Comprehensive Differential Expression Analysis for Morphine OUD Study
# ===================================================================

# Set reproducibility seed
set.seed(42)

# Load required libraries
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  library(RColorBrewer)
  library(VennDiagram)
  library(gridExtra)
  # Optional packages
  tryCatch(library(clusterProfiler), error = function(e) cat("clusterProfiler not available\n"))
  tryCatch(library(org.Mm.eg.db), error = function(e) cat("org.Mm.eg.db not available\n"))
  tryCatch(library(DOSE), error = function(e) cat("DOSE not available\n"))
  tryCatch(library(enrichplot), error = function(e) cat("enrichplot not available\n"))
})

# ----------------------
# SETUP AND DATA LOADING
# ----------------------

# Base directory and paths
base_dir <- "/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus/GSE289002"
input_dir <- file.path(base_dir, "Outputs/Enhanced_QC/Processed_Data")

# Create output directories
output_dirs <- list(
  base = file.path(base_dir, "Outputs/Differential_Expression"),
  results = file.path(base_dir, "Outputs/Differential_Expression/Results"),
  plots = file.path(base_dir, "Outputs/Differential_Expression/Plots"),
  tables = file.path(base_dir, "Outputs/Differential_Expression/Tables"),
  pathway = file.path(base_dir, "Outputs/Differential_Expression/Pathway_Analysis"),
  region_specific = file.path(base_dir, "Outputs/Differential_Expression/Region_Specific"),
  sex_specific = file.path(base_dir, "Outputs/Differential_Expression/Sex_Specific"),
  time_course = file.path(base_dir, "Outputs/Differential_Expression/Time_Course")
)

# Create all directories
lapply(output_dirs, function(x) dir.create(x, showWarnings = FALSE, recursive = TRUE))

cat("=== COMPREHENSIVE DIFFERENTIAL EXPRESSION ANALYSIS ===\n")
cat("Analysis Date:", as.character(Sys.time()), "\n")
cat("Output Directory:", output_dirs$base, "\n\n")

# Load processed data
cat("Loading preprocessed data...\n")
dge <- readRDS(file.path(input_dir, "dge_normalized_final.rds"))
metadata <- readRDS(file.path(input_dir, "metadata_enhanced.rds"))

# Verify data
cat("Dataset dimensions:", nrow(dge), "genes x", ncol(dge), "samples\n")
cat("Treatment groups:", paste(levels(metadata$treatment), collapse = ", "), "\n")
cat("Brain regions:", paste(levels(metadata$region), collapse = ", "), "\n")
cat("Sex distribution:", paste(names(table(metadata$sex)), ":", table(metadata$sex), collapse = "; "), "\n\n")

# Display experimental design
cat("Experimental Design Summary:\n")
design_summary <- table(metadata$treatment, metadata$region, metadata$sex)
print(design_summary)
cat("\n")

# ----------------------
# 1. MAIN DESIGN MATRIX AND MODEL SETUP
# ----------------------

cat("=== SETTING UP DESIGN MATRIX AND MODEL ===\n")

# Create comprehensive design matrix
design_full <- model.matrix(~ 0 + treatment + sex + region, data = metadata)

# FIX: Clean up column names to be R-valid
colnames(design_full) <- gsub("treatment", "", colnames(design_full))
colnames(design_full) <- gsub("sex", "Sex_", colnames(design_full))
colnames(design_full) <- gsub("region", "Region_", colnames(design_full))

# Make names syntactically valid
colnames(design_full) <- make.names(colnames(design_full))

cat("Design matrix columns:", paste(colnames(design_full), collapse = ", "), "\n")

# Estimate dispersions
cat("Estimating dispersions...\n")
dge <- estimateDisp(dge, design_full)

# Plot dispersion estimates
png(file.path(output_dirs$plots, "dispersion_estimates.png"), width = 1200, height = 800)
plotBCV(dge, main = "Biological Coefficient of Variation")
dev.off()

# Fit quasi-likelihood GLM
cat("Fitting GLM model...\n")
fit_full <- glmQLFit(dge, design_full)

# Plot quasi-likelihood dispersions
png(file.path(output_dirs$plots, "ql_dispersions.png"), width = 1200, height = 800)
plotQLDisp(fit_full, main = "Quasi-likelihood Dispersions")
dev.off()

# ----------------------
# 2. DEFINE ALL CONTRASTS
# ----------------------

cat("\n=== DEFINING CONTRASTS FOR ALL COMPARISONS ===\n")

# Print actual column names for reference
cat("Available design matrix columns:\n")
print(colnames(design_full))

# Main treatment contrasts - using the actual cleaned names
main_contrasts <- makeContrasts(
  # Primary comparisons vs Control (Saline)
  Acute_vs_Control = Mor...24h - Sal,
  Short_vs_Control = Mor...2W - Sal,
  Chronic_vs_Control = Chronic.mor - Sal,
  
  # Time course comparisons
  Short_vs_Acute = Mor...2W - Mor...24h,
  Chronic_vs_Acute = Chronic.mor - Mor...24h,
  Chronic_vs_Short = Chronic.mor - Mor...2W,
  
  levels = design_full
)

cat("Main contrasts defined:", ncol(main_contrasts), "comparisons\n")
print(colnames(main_contrasts))

# ----------------------
# 3. PERFORM DIFFERENTIAL EXPRESSION TESTS
# ----------------------

cat("\n=== PERFORMING DIFFERENTIAL EXPRESSION TESTS ===\n")

# Function to perform DE analysis and extract results
perform_de_analysis <- function(fit, contrast, contrast_name, fdr_cutoff = 0.05, lfc_cutoff = 1) {
  
  cat("Analyzing:", contrast_name, "\n")
  
  # Perform test
  qlf <- glmQLFTest(fit, contrast = contrast)
  
  # Extract results
  results <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  results$gene_id <- rownames(results)
  
  # Add significance classification
  results$significant <- (results$FDR < fdr_cutoff) & (abs(results$logFC) > log2(lfc_cutoff))
  results$direction <- ifelse(results$logFC > 0, "Up", "Down")
  results$regulation <- ifelse(results$significant, 
                              paste(results$direction, "regulated"), 
                              "Not significant")
  
  # Summary statistics
  summary_stats <- list(
    total_genes = nrow(results),
    significant_genes = sum(results$significant),
    upregulated = sum(results$significant & results$logFC > 0),
    downregulated = sum(results$significant & results$logFC < 0),
    max_logfc_up = ifelse(any(results$logFC > 0), max(results$logFC), 0),
    max_logfc_down = ifelse(any(results$logFC < 0), min(results$logFC), 0),
    min_fdr = min(results$FDR)
  )
  
  cat("  - Total genes tested:", summary_stats$total_genes, "\n")
  cat("  - Significant genes:", summary_stats$significant_genes, "\n")
  cat("  - Upregulated:", summary_stats$upregulated, "\n")
  cat("  - Downregulated:", summary_stats$downregulated, "\n")
  
  return(list(
    results = results,
    qlf = qlf,
    summary = summary_stats,
    contrast_name = contrast_name
  ))
}

# Perform all main DE analyses
de_results <- list()
for (i in 1:ncol(main_contrasts)) {
  contrast_name <- colnames(main_contrasts)[i]
  de_results[[contrast_name]] <- perform_de_analysis(
    fit_full, 
    main_contrasts[, i], 
    contrast_name
  )
}

# ----------------------
# 4. SAVE RESULTS TABLES
# ----------------------

cat("\n=== SAVING RESULTS TABLES ===\n")

# Save individual results
for (contrast_name in names(de_results)) {
  
  # Full results
  write.csv(de_results[[contrast_name]]$results, 
            file.path(output_dirs$tables, paste0(contrast_name, "_full_results.csv")), 
            row.names = FALSE)
  
  # Significant genes only
  sig_genes <- de_results[[contrast_name]]$results[de_results[[contrast_name]]$results$significant, ]
  write.csv(sig_genes, 
            file.path(output_dirs$tables, paste0(contrast_name, "_significant_genes.csv")), 
            row.names = FALSE)
  
  cat("Saved results for:", contrast_name, "\n")
}

# Create comprehensive summary table
summary_table <- data.frame(
  Comparison = names(de_results),
  Total_Genes = sapply(de_results, function(x) x$summary$total_genes),
  Significant_Genes = sapply(de_results, function(x) x$summary$significant_genes),
  Upregulated = sapply(de_results, function(x) x$summary$upregulated),
  Downregulated = sapply(de_results, function(x) x$summary$downregulated),
  Max_LogFC_Up = round(sapply(de_results, function(x) x$summary$max_logfc_up), 3),
  Max_LogFC_Down = round(sapply(de_results, function(x) x$summary$max_logfc_down), 3),
  Min_FDR = sapply(de_results, function(x) x$summary$min_fdr),
  stringsAsFactors = FALSE
)

write.csv(summary_table, file.path(output_dirs$tables, "DE_analysis_summary.csv"), row.names = FALSE)
cat("Comprehensive summary saved\n")

# ----------------------
# 5. VISUALIZATION: VOLCANO PLOTS
# ----------------------

cat("\n=== CREATING VOLCANO PLOTS ===\n")

create_volcano_plot <- function(results, title, output_file) {
  
  # Prepare data
  volcano_data <- results$results
  volcano_data$neg_log10_fdr <- -log10(volcano_data$FDR)
  
  # Create plot
  p <- ggplot(volcano_data, aes(x = logFC, y = neg_log10_fdr)) +
    geom_point(aes(color = regulation), alpha = 0.6, size = 1) +
    scale_color_manual(values = c("Up regulated" = "red", 
                                 "Down regulated" = "blue", 
                                 "Not significant" = "grey")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    labs(title = title,
         x = "Log2 Fold Change",
         y = "-Log10 FDR",
         color = "Regulation") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.position = "top") +
    annotate("text", x = Inf, y = Inf, 
             label = paste("Significant genes:", results$summary$significant_genes),
             hjust = 1.1, vjust = 2, size = 4)
  
  ggsave(output_file, p, width = 10, height = 8, dpi = 300)
  return(p)
}

# Create volcano plots for all comparisons
volcano_plots <- list()
for (contrast_name in names(de_results)) {
  output_file <- file.path(output_dirs$plots, paste0("volcano_", contrast_name, ".png"))
  volcano_plots[[contrast_name]] <- create_volcano_plot(
    de_results[[contrast_name]], 
    paste("Volcano Plot:", contrast_name),
    output_file
  )
  cat("Created volcano plot for:", contrast_name, "\n")
}

# ----------------------
# 6. HEATMAP OF TOP DE GENES
# ----------------------

cat("\n=== CREATING HEATMAPS ===\n")

create_de_heatmap <- function(dge, metadata, de_results_list, top_n = 50) {
  
  # Get top DE genes from each comparison
  all_top_genes <- c()
  for (contrast_name in names(de_results_list)) {
    sig_genes <- de_results_list[[contrast_name]]$results[de_results_list[[contrast_name]]$results$significant, ]
    if (nrow(sig_genes) > 0) {
      top_genes <- head(sig_genes[order(sig_genes$FDR), ], top_n)$gene_id
      all_top_genes <- c(all_top_genes, top_genes)
    }
  }
  
  # Remove duplicates and get final gene set
  unique_top_genes <- unique(all_top_genes)
  
  if (length(unique_top_genes) == 0) {
    cat("No significant genes found for heatmap\n")
    return(NULL)
  }
  
  cat("Creating heatmap with", length(unique_top_genes), "genes\n")
  
  # Get log-CPM values
  logcpm <- cpm(dge, log = TRUE)
  heatmap_data <- logcpm[unique_top_genes, ]
  
  # Scale by row (z-score)
  heatmap_data_scaled <- t(scale(t(heatmap_data)))
  
  # Create annotation
  annotation_col <- data.frame(
    Treatment = metadata$treatment,
    Region = metadata$region,
    Sex = metadata$sex,
    row.names = colnames(heatmap_data_scaled)
  )
  
  # Color schemes
  ann_colors <- list(
    Treatment = RColorBrewer::brewer.pal(length(levels(metadata$treatment)), "Set1"),
    Region = c("NAc" = "lightblue", "PFC" = "lightcoral"),
    Sex = c("female" = "pink", "male" = "lightblue")
  )
  names(ann_colors$Treatment) <- levels(metadata$treatment)
  
  # Create heatmap
  png(file.path(output_dirs$plots, "top_DE_genes_heatmap.png"), width = 1600, height = 1200, res = 300)
  pheatmap(heatmap_data_scaled,
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           show_rownames = FALSE,
           show_colnames = FALSE,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           scale = "none",
           color = colorRampPalette(c("blue", "white", "red"))(100),
           main = paste("Top DE Genes Across All Comparisons (n =", length(unique_top_genes), ")"),
           fontsize = 10)
  dev.off()
  
  return(unique_top_genes)
}

# Create comprehensive heatmap
top_de_genes <- create_de_heatmap(dge, metadata, de_results)

# ----------------------
# 7. QUESTION 1: ACUTE vs CHRONIC RESPONSE GENES
# ----------------------

cat("\n=== QUESTION 1: ACUTE vs CHRONIC RESPONSE ANALYSIS ===\n")

analyze_acute_vs_chronic <- function(de_results) {
  
  # Get significant genes from acute and chronic comparisons
  acute_genes <- de_results$Acute_vs_Control$results[de_results$Acute_vs_Control$results$significant, ]
  chronic_genes <- de_results$Chronic_vs_Control$results[de_results$Chronic_vs_Control$results$significant, ]
  
  # Gene sets
  acute_only <- setdiff(acute_genes$gene_id, chronic_genes$gene_id)
  chronic_only <- setdiff(chronic_genes$gene_id, acute_genes$gene_id)
  shared_genes <- intersect(acute_genes$gene_id, chronic_genes$gene_id)
  
  # Create comparison summary
  acute_chronic_summary <- data.frame(
    Category = c("Acute only", "Chronic only", "Shared", "Total acute", "Total chronic"),
    Count = c(length(acute_only), length(chronic_only), length(shared_genes),
              nrow(acute_genes), nrow(chronic_genes))
  )
  
  cat("Acute vs Chronic Response Analysis:\n")
  print(acute_chronic_summary)
  
  # Save gene lists
  write.csv(data.frame(gene_id = acute_only), 
            file.path(output_dirs$results, "acute_specific_genes.csv"), row.names = FALSE)
  write.csv(data.frame(gene_id = chronic_only), 
            file.path(output_dirs$results, "chronic_specific_genes.csv"), row.names = FALSE)
  write.csv(data.frame(gene_id = shared_genes), 
            file.path(output_dirs$results, "acute_chronic_shared_genes.csv"), row.names = FALSE)
  
  # Create Venn diagram
  png(file.path(output_dirs$plots, "acute_vs_chronic_venn.png"), width = 800, height = 800)
  venn_plot <- venn.diagram(
    x = list(
      Acute = acute_genes$gene_id,
      Chronic = chronic_genes$gene_id
    ),
    category.names = c("Acute (24h)", "Chronic"),
    filename = NULL,
    fill = c("lightblue", "lightcoral"),
    alpha = 0.7,
    cex = 2,
    fontfamily = "sans",
    cat.cex = 1.5,
    cat.fontfamily = "sans"
  )
  grid.draw(venn_plot)
  dev.off()
  
  return(list(
    acute_only = acute_only,
    chronic_only = chronic_only,
    shared = shared_genes,
    summary = acute_chronic_summary
  ))
}

acute_chronic_analysis <- analyze_acute_vs_chronic(de_results)

# ----------------------
# 8. QUESTION 2: REGION-SPECIFIC ANALYSIS (PFC vs NAc)
# ----------------------

cat("\n=== QUESTION 2: REGION-SPECIFIC ANALYSIS ===\n")

perform_region_specific_analysis <- function(dge, metadata) {
  
  cat("Performing region-specific differential expression analysis...\n")
  
  region_results <- list()
  
  for (region in c("PFC", "NAc")) {
    
    cat("Analyzing region:", region, "\n")
    
    # Subset data for this region
    region_samples <- metadata$title[metadata$region == region]
    region_metadata <- metadata[metadata$region == region, ]
    region_dge <- dge[, region_samples]
    
    # Create design matrix for this region
    region_design <- model.matrix(~ 0 + treatment + sex, data = region_metadata)
    colnames(region_design) <- gsub("treatment", "", colnames(region_design))
    colnames(region_design) <- gsub("sex", "Sex_", colnames(region_design))
    
    # FIX: Make names syntactically valid for region-specific analysis
    colnames(region_design) <- make.names(colnames(region_design))
    
    cat("Region design matrix columns for", region, ":", paste(colnames(region_design), collapse = ", "), "\n")
    
    # Fit model
    region_dge <- estimateDisp(region_dge, region_design)
    region_fit <- glmQLFit(region_dge, region_design)
    
    # Define contrasts - using cleaned names
    region_contrasts <- makeContrasts(
      Acute_vs_Control = Mor...24h - Sal,
      Short_vs_Control = Mor...2W - Sal,
      Chronic_vs_Control = Chronic.mor - Sal,
      levels = region_design
    )
    
    # Perform DE analysis for each contrast
    region_de_results <- list()
    for (i in 1:ncol(region_contrasts)) {
      contrast_name <- paste(region, colnames(region_contrasts)[i], sep = "_")
      region_de_results[[contrast_name]] <- perform_de_analysis(
        region_fit, 
        region_contrasts[, i], 
        contrast_name
      )
    }
    
    region_results[[region]] <- region_de_results
  }
  
  return(region_results)
}

region_specific_results <- perform_region_specific_analysis(dge, metadata)

# Save region-specific results
cat("Saving region-specific results...\n")
for (region in names(region_specific_results)) {
  for (contrast_name in names(region_specific_results[[region]])) {
    
    # Save full results
    write.csv(region_specific_results[[region]][[contrast_name]]$results,
              file.path(output_dirs$region_specific, paste0(contrast_name, "_full_results.csv")),
              row.names = FALSE)
    
    # Save significant genes
    sig_genes <- region_specific_results[[region]][[contrast_name]]$results[
      region_specific_results[[region]][[contrast_name]]$results$significant, ]
    write.csv(sig_genes,
              file.path(output_dirs$region_specific, paste0(contrast_name, "_significant_genes.csv")),
              row.names = FALSE)
  }
}

# Compare PFC vs NAc responses
compare_region_responses <- function(region_results) {
  
  cat("Comparing PFC vs NAc responses...\n")
  
  comparisons <- c("Acute_vs_Control", "Short_vs_Control", "Chronic_vs_Control")
  region_comparison_summary <- data.frame()
  
  for (comp in comparisons) {
    
    pfc_genes <- region_results$PFC[[paste("PFC", comp, sep = "_")]]$results[
      region_results$PFC[[paste("PFC", comp, sep = "_")]]$results$significant, ]$gene_id
    
    nac_genes <- region_results$NAc[[paste("NAc", comp, sep = "_")]]$results[
      region_results$NAc[[paste("NAc", comp, sep = "_")]]$results$significant, ]$gene_id
    
    pfc_only <- setdiff(pfc_genes, nac_genes)
    nac_only <- setdiff(nac_genes, pfc_genes)
    shared <- intersect(pfc_genes, nac_genes)
    
    summary_row <- data.frame(
      Comparison = comp,
      PFC_only = length(pfc_only),
      NAc_only = length(nac_only),
      Shared = length(shared),
      Total_PFC = length(pfc_genes),
      Total_NAc = length(nac_genes)
    )
    
    region_comparison_summary <- rbind(region_comparison_summary, summary_row)
    
    # Create Venn diagram for this comparison
    png(file.path(output_dirs$region_specific, paste0("region_comparison_", comp, "_venn.png")), 
        width = 800, height = 800)
    venn_plot <- venn.diagram(
      x = list(
        PFC = pfc_genes,
        NAc = nac_genes
      ),
      category.names = c("PFC", "NAc"),
      filename = NULL,
      fill = c("lightblue", "lightgreen"),
      alpha = 0.7,
      cex = 2,
      fontfamily = "sans",
      cat.cex = 1.5,
      cat.fontfamily = "sans",
      main = paste("Region Comparison:", comp)
    )
    grid.draw(venn_plot)
    dev.off()
  }
  
  write.csv(region_comparison_summary, 
            file.path(output_dirs$region_specific, "region_comparison_summary.csv"), 
            row.names = FALSE)
  
  cat("Region comparison analysis completed\n")
  print(region_comparison_summary)
  
  return(region_comparison_summary)
}

region_comparison_summary <- compare_region_responses(region_specific_results)

# ----------------------
# 9. QUESTION 3: SEX-SPECIFIC ANALYSIS
# ----------------------

cat("\n=== QUESTION 3: SEX-SPECIFIC ANALYSIS ===\n")

perform_sex_specific_analysis <- function(dge, metadata) {
  
  cat("Performing sex-specific differential expression analysis...\n")
  
  sex_results <- list()
  
  for (sex in c("male", "female")) {
    
    cat("Analyzing sex:", sex, "\n")
    
    # Subset data for this sex
    sex_samples <- metadata$title[metadata$sex == sex]
    sex_metadata <- metadata[metadata$sex == sex, ]
    sex_dge <- dge[, sex_samples]
    
    # Create design matrix for this sex
    sex_design <- model.matrix(~ 0 + treatment + region, data = sex_metadata)
    colnames(sex_design) <- gsub("treatment", "", colnames(sex_design))
    colnames(sex_design) <- gsub("region", "Region_", colnames(sex_design))
    
    # FIX: Make names syntactically valid for sex-specific analysis
    colnames(sex_design) <- make.names(colnames(sex_design))
    
    cat("Sex design matrix columns for", sex, ":", paste(colnames(sex_design), collapse = ", "), "\n")
    
    # Fit model
    sex_dge <- estimateDisp(sex_dge, sex_design)
    sex_fit <- glmQLFit(sex_dge, sex_design)
    
    # Define contrasts - using cleaned names
    sex_contrasts <- makeContrasts(
      Acute_vs_Control = Mor...24h - Sal,
      Short_vs_Control = Mor...2W - Sal,
      Chronic_vs_Control = Chronic.mor - Sal,
      levels = sex_design
    )
    
    # Perform DE analysis for each contrast
    sex_de_results <- list()
    for (i in 1:ncol(sex_contrasts)) {
      contrast_name <- paste(sex, colnames(sex_contrasts)[i], sep = "_")
      sex_de_results[[contrast_name]] <- perform_de_analysis(
        sex_fit, 
        sex_contrasts[, i], 
        contrast_name
      )
    }
    
    sex_results[[sex]] <- sex_de_results
  }
  
  return(sex_results)
}

sex_specific_results <- perform_sex_specific_analysis(dge, metadata)

# Save sex-specific results
cat("Saving sex-specific results...\n")
for (sex in names(sex_specific_results)) {
  for (contrast_name in names(sex_specific_results[[sex]])) {
    
    # Save full results
    write.csv(sex_specific_results[[sex]][[contrast_name]]$results,
              file.path(output_dirs$sex_specific, paste0(contrast_name, "_full_results.csv")),
              row.names = FALSE)
    
    # Save significant genes
    sig_genes <- sex_specific_results[[sex]][[contrast_name]]$results[
      sex_specific_results[[sex]][[contrast_name]]$results$significant, ]
    write.csv(sig_genes,
              file.path(output_dirs$sex_specific, paste0(contrast_name, "_significant_genes.csv")),
              row.names = FALSE)
  }
}

# Compare male vs female responses
compare_sex_responses <- function(sex_results) {
  
  cat("Comparing male vs female responses...\n")
  
  comparisons <- c("Acute_vs_Control", "Short_vs_Control", "Chronic_vs_Control")
  sex_comparison_summary <- data.frame()
  
  for (comp in comparisons) {
    
    male_genes <- sex_results$male[[paste("male", comp, sep = "_")]]$results[
      sex_results$male[[paste("male", comp, sep = "_")]]$results$significant, ]$gene_id
    
    female_genes <- sex_results$female[[paste("female", comp, sep = "_")]]$results[
      sex_results$female[[paste("female", comp, sep = "_")]]$results$significant, ]$gene_id
    
    male_only <- setdiff(male_genes, female_genes)
    female_only <- setdiff(female_genes, male_genes)
    shared <- intersect(male_genes, female_genes)
    
    summary_row <- data.frame(
      Comparison = comp,
      Male_only = length(male_only),
      Female_only = length(female_only),
      Shared = length(shared),
      Total_Male = length(male_genes),
      Total_Female = length(female_genes)
    )
    
    sex_comparison_summary <- rbind(sex_comparison_summary, summary_row)
    
    # Create Venn diagram for this comparison
    png(file.path(output_dirs$sex_specific, paste0("sex_comparison_", comp, "_venn.png")), 
        width = 800, height = 800)
    venn_plot <- venn.diagram(
      x = list(
        Male = male_genes,
        Female = female_genes
      ),
      category.names = c("Male", "Female"),
      filename = NULL,
      fill = c("lightblue", "pink"),
      alpha = 0.7,
      cex = 2,
      fontfamily = "sans",
      cat.cex = 1.5,
      cat.fontfamily = "sans",
      main = paste("Sex Comparison:", comp)
    )
    grid.draw(venn_plot)
    dev.off()
  }
  
  write.csv(sex_comparison_summary, 
            file.path(output_dirs$sex_specific, "sex_comparison_summary.csv"), 
            row.names = FALSE)
  
  cat("Sex comparison analysis completed\n")
  print(sex_comparison_summary)
  
  return(sex_comparison_summary)
}

sex_comparison_summary <- compare_sex_responses(sex_specific_results)

# ----------------------
# 10. QUESTION 4: TIME COURSE ANALYSIS
# ----------------------

cat("\n=== QUESTION 4: TIME COURSE ANALYSIS ===\n")

analyze_time_course <- function(de_results) {
  
  cat("Analyzing temporal patterns of gene expression...\n")
  
  # Get significant genes from each time point
  acute_genes <- de_results$Acute_vs_Control$results[de_results$Acute_vs_Control$results$significant, ]
  short_genes <- de_results$Short_vs_Control$results[de_results$Short_vs_Control$results$significant, ]
  chronic_genes <- de_results$Chronic_vs_Control$results[de_results$Chronic_vs_Control$results$significant, ]
  
  # Create time course patterns
  all_time_genes <- unique(c(acute_genes$gene_id, short_genes$gene_id, chronic_genes$gene_id))
  
  time_course_patterns <- data.frame(
    gene_id = all_time_genes,
    acute_sig = all_time_genes %in% acute_genes$gene_id,
    short_sig = all_time_genes %in% short_genes$gene_id,
    chronic_sig = all_time_genes %in% chronic_genes$gene_id,
    stringsAsFactors = FALSE
  )
  
  # Define temporal patterns
  time_course_patterns$pattern <- with(time_course_patterns, 
    ifelse(acute_sig & !short_sig & !chronic_sig, "Acute_only",
    ifelse(!acute_sig & short_sig & !chronic_sig, "Short_only",
    ifelse(!acute_sig & !short_sig & chronic_sig, "Chronic_only",
    ifelse(acute_sig & short_sig & !chronic_sig, "Early_response",
    ifelse(acute_sig & !short_sig & chronic_sig, "Acute_and_chronic",
    ifelse(!acute_sig & short_sig & chronic_sig, "Late_response",
    ifelse(acute_sig & short_sig & chronic_sig, "Sustained_response",
           "Other"))))))))
  
  # Count patterns
  pattern_counts <- table(time_course_patterns$pattern)
  
  cat("Temporal expression patterns:\n")
  print(pattern_counts)
  
  # Save time course analysis
  write.csv(time_course_patterns, 
            file.path(output_dirs$time_course, "temporal_expression_patterns.csv"), 
            row.names = FALSE)
  
  # Create pattern visualization
  pattern_summary <- data.frame(
    Pattern = names(pattern_counts),
    Count = as.numeric(pattern_counts)
  )
  
  p_patterns <- ggplot(pattern_summary, aes(x = reorder(Pattern, Count), y = Count)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    coord_flip() +
    labs(title = "Temporal Expression Patterns",
         x = "Expression Pattern",
         y = "Number of Genes") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  ggsave(file.path(output_dirs$time_course, "temporal_patterns_barplot.png"), 
         p_patterns, width = 10, height = 6, dpi = 300)
  
  return(list(
    patterns = time_course_patterns,
    pattern_counts = pattern_counts,
    pattern_summary = pattern_summary
  ))
}

time_course_analysis <- analyze_time_course(de_results)

# ----------------------
# 11. COMPREHENSIVE SUMMARY REPORT
# ----------------------

cat("\n=== GENERATING COMPREHENSIVE SUMMARY REPORT ===\n")

create_de_summary_report <- function(output_file) {
  
  sink(output_file)
  
  cat("================================================================\n")
  cat("COMPREHENSIVE DIFFERENTIAL EXPRESSION ANALYSIS REPORT\n")
  cat("================================================================\n\n")
  
  cat("Analysis Information:\n")
  cat("- Date:", as.character(Sys.time()), "\n")
  cat("- Dataset:", nrow(dge), "genes x", ncol(dge), "samples\n")
  cat("- Treatment groups:", paste(levels(metadata$treatment), collapse = ", "), "\n")
  cat("- Brain regions:", paste(levels(metadata$region), collapse = ", "), "\n")
  cat("- Sex groups:", paste(levels(metadata$sex), collapse = ", "), "\n\n")
  
  cat("=== MAIN DIFFERENTIAL EXPRESSION RESULTS ===\n\n")
  
  cat("Summary of all comparisons:\n")
  print(summary_table)
  cat("\n")
  
  cat("=== QUESTION 1: ACUTE vs CHRONIC RESPONSE GENES ===\n\n")
  cat("Acute vs Chronic Response Analysis:\n")
  print(acute_chronic_analysis$summary)
  cat("\nKey findings:\n")
  cat("- Acute-specific genes:", length(acute_chronic_analysis$acute_only), "\n")
  cat("- Chronic-specific genes:", length(acute_chronic_analysis$chronic_only), "\n")
  cat("- Shared acute/chronic genes:", length(acute_chronic_analysis$shared), "\n\n")
  
  cat("=== QUESTION 2: REGION-SPECIFIC RESPONSES ===\n\n")
  cat("PFC vs NAc Response Comparison:\n")
  print(region_comparison_summary)
  cat("\nKey findings:\n")
  cat("- Region-specific responses vary by treatment duration\n")
  cat("- See individual Venn diagrams for detailed overlaps\n\n")
  
  cat("=== QUESTION 3: SEX-SPECIFIC RESPONSES ===\n\n")
  cat("Male vs Female Response Comparison:\n")
  print(sex_comparison_summary)
  cat("\nKey findings:\n")
  cat("- Sex differences in morphine response patterns\n")
  cat("- See individual analyses for sex-specific gene lists\n\n")
  
  cat("=== QUESTION 4: TEMPORAL EXPRESSION PATTERNS ===\n\n")
  cat("Time Course Expression Patterns:\n")
  print(time_course_analysis$pattern_counts)
  cat("\nKey findings:\n")
  if("Sustained_response" %in% names(time_course_analysis$pattern_counts)) {
    cat("- Sustained response genes:", time_course_analysis$pattern_counts["Sustained_response"], "\n")
  }
  if("Acute_only" %in% names(time_course_analysis$pattern_counts)) {
    cat("- Acute-only response genes:", time_course_analysis$pattern_counts["Acute_only"], "\n")
  }
  if("Chronic_only" %in% names(time_course_analysis$pattern_counts)) {
    cat("- Chronic-only response genes:", time_course_analysis$pattern_counts["Chronic_only"], "\n")
  }
  cat("\n")
  
  cat("=== FILES GENERATED ===\n\n")
  
  cat("Main Results:\n")
  cat("- DE_analysis_summary.csv: Overview of all comparisons\n")
  cat("- Individual comparison results in Tables/ directory\n")
  cat("- Volcano plots for all comparisons in Plots/ directory\n\n")
  
  cat("Question-Specific Analyses:\n")
  cat("- acute_specific_genes.csv, chronic_specific_genes.csv\n")
  cat("- Region-specific results in Region_Specific/ directory\n")
  cat("- Sex-specific results in Sex_Specific/ directory\n")
  cat("- Temporal patterns in Time_Course/ directory\n\n")
  
  cat("Visualizations:\n")
  cat("- Volcano plots for all comparisons\n")
  cat("- Venn diagrams for overlaps\n")
  cat("- Heatmap of top DE genes\n")
  cat("- Temporal pattern plots\n\n")
  
  cat("=== RECOMMENDATIONS FOR FURTHER ANALYSIS ===\n\n")
  cat("1. Pathway enrichment analysis on gene sets\n")
  cat("2. Gene ontology analysis for biological processes\n")
  cat("3. Network analysis of co-expressed genes\n")
  cat("4. Validation of top candidates by qPCR\n")
  cat("5. Investigation of transcription factor networks\n\n")
  
  cat("================================================================\n")
  cat("ANALYSIS COMPLETE\n")
  cat("================================================================\n")
  
  sink()
}

# Generate comprehensive report
create_de_summary_report(file.path(output_dirs$results, "comprehensive_DE_analysis_report.txt"))

# Save session information
writeLines(capture.output(sessionInfo()), 
           file.path(output_dirs$results, "session_info.txt"))

# ----------------------
# 12. FINAL STATUS SUMMARY
# ----------------------

cat("\n")
cat("================================================================\n")
cat("DIFFERENTIAL EXPRESSION ANALYSIS COMPLETE\n")
cat("================================================================\n")
cat("Analysis completed at:", as.character(Sys.time()), "\n\n")

cat("KEY FINDINGS SUMMARY:\n")
cat("- Total comparisons performed:", length(de_results), "\n")
cat("- Most significant comparison:", summary_table$Comparison[which.max(summary_table$Significant_Genes)], 
    "with", max(summary_table$Significant_Genes), "DE genes\n")
cat("- Acute-specific genes:", length(acute_chronic_analysis$acute_only), "\n")
cat("- Chronic-specific genes:", length(acute_chronic_analysis$chronic_only), "\n")
cat("- Shared acute/chronic genes:", length(acute_chronic_analysis$shared), "\n\n")

cat("OUTPUT DIRECTORIES:\n")
cat("- Main results:", output_dirs$results, "\n")
cat("- All tables:", output_dirs$tables, "\n")
cat("- All plots:", output_dirs$plots, "\n")
cat("- Region-specific:", output_dirs$region_specific, "\n")
cat("- Sex-specific:", output_dirs$sex_specific, "\n")
cat("- Time course:", output_dirs$time_course, "\n\n")

cat("READY FOR:\n")
cat("- Pathway enrichment analysis\n")
cat("- Gene ontology analysis\n")
cat("- Network analysis\n")
cat("- Publication figure generation\n")
cat("- Experimental validation planning\n")

cat("================================================================\n")