# ==============================================================================
# Sex-Specific Neuroinflammatory Analysis
# ==============================================================================
# Purpose: Analyze sex-specific differential expression and interactions
# Datasets: GSE174409 (Bulk RNA-seq) + GSE225158 (snRNA-seq pseudobulk)
# Focus: Sex main effects, sex×region interactions, sex×OUD interactions
# ==============================================================================

# Configuration
BASE_DIR <- "/Users/aumchampaneri/Complement-OUD/Translational_Study/Homo_sapiens"
GSE174409_DIR <- file.path(BASE_DIR, "Results", "GSE174409")
GSE225158_DIR <- file.path(BASE_DIR, "Results", "GSE225158")
SEX_ANALYSIS_DIR <- file.path(BASE_DIR, "Results", "Sex_Specific_Analysis")

# Load libraries
suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(VennDiagram)
  library(grid)
})

# Create output directory
if (!dir.exists(SEX_ANALYSIS_DIR)) dir.create(SEX_ANALYSIS_DIR, recursive = TRUE)

# ==============================================================================
# DATA LOADING FOR SEX ANALYSIS
# ==============================================================================

#' Load datasets with focus on sex variables
load_sex_analysis_data <- function() {
  cat("=== Loading Data for Sex-Specific Analysis ===\n")
  
  # Load GSE174409 results
  gse174409_file <- file.path(GSE174409_DIR, "GSE174409_region_analysis_results.rds")
  gse225158_file <- file.path(GSE225158_DIR, "GSE225158_region_analysis_results.rds")
  
  if (!file.exists(gse174409_file) || !file.exists(gse225158_file)) {
    stop("Previous analysis results not found. Run individual analyses first.")
  }
  
  gse174409_data <- readRDS(gse174409_file)
  gse225158_data <- readRDS(gse225158_file)
  
  # Extract expression data and metadata
  gse174409_expr <- gse174409_data$expression_data
  gse225158_expr <- gse225158_data$expression_data
  
  cat("GSE174409 sex distribution:", table(gse174409_expr$metadata$Sex), "\n")
  cat("GSE225158 sex distribution:", table(gse225158_expr$metadata$Sex), "\n")
  
  return(list(
    gse174409 = gse174409_expr,
    gse225158 = gse225158_expr
  ))
}

# ==============================================================================
# SEX MAIN EFFECTS ANALYSIS
# ==============================================================================

#' Analyze sex main effects in each dataset
analyze_sex_main_effects <- function(data) {
  cat("\n=== Analyzing Sex Main Effects ===\n")
  
  sex_results <- list()
  
  for (dataset_name in names(data)) {
    cat("\nAnalyzing", dataset_name, "\n")
    
    expr_data <- data[[dataset_name]]
    counts <- expr_data$filtered_counts
    metadata <- expr_data$metadata
    
    # Create DGEList
    dge <- DGEList(counts = counts)
    keep <- filterByExpr(dge, group = metadata$Sex)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- calcNormFactors(dge)
    
    # Design matrix with sex as main effect
    design <- model.matrix(~ Sex + Region + OUD_Status, data = metadata)
    
    # Voom transformation
    v <- voom(dge, design, plot = FALSE)
    
    # Account for subject pairing if available
    if ("Subject_ID" %in% colnames(metadata)) {
      corfit <- duplicateCorrelation(v, design, block = metadata$Subject_ID)
      fit <- lmFit(v, design, block = metadata$Subject_ID, correlation = corfit$consensus)
    } else {
      fit <- lmFit(v, design)
    }
    
    fit <- eBayes(fit)
    
    # Extract sex effect
    sex_coef <- grep("Sex", colnames(design), value = TRUE)[1]
    if (length(sex_coef) > 0) {
      results <- topTable(fit, coef = sex_coef, number = Inf, sort.by = "P")
      results$Dataset <- dataset_name
      sex_results[[dataset_name]] <- results
      
      n_sig <- sum(results$adj.P.Val < 0.05, na.rm = TRUE)
      cat("  Sex main effect:", n_sig, "significant genes (FDR < 0.05)\n")
    }
  }
  
  return(sex_results)
}

# ==============================================================================
# SEX × REGION INTERACTION ANALYSIS
# ==============================================================================

#' Analyze sex × region interactions
analyze_sex_region_interactions <- function(data) {
  cat("\n=== Analyzing Sex × Region Interactions ===\n")
  
  interaction_results <- list()
  
  for (dataset_name in names(data)) {
    cat("\nAnalyzing", dataset_name, "\n")
    
    expr_data <- data[[dataset_name]]
    counts <- expr_data$filtered_counts
    metadata <- expr_data$metadata
    
    # Create DGEList
    dge <- DGEList(counts = counts)
    keep <- filterByExpr(dge, group = interaction(metadata$Sex, metadata$Region))
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- calcNormFactors(dge)
    
    # Design matrix with interaction
    design <- model.matrix(~ Sex * Region + OUD_Status, data = metadata)
    
    # Voom transformation
    v <- voom(dge, design, plot = FALSE)
    
    # Account for subject pairing
    if ("Subject_ID" %in% colnames(metadata)) {
      corfit <- duplicateCorrelation(v, design, block = metadata$Subject_ID)
      fit <- lmFit(v, design, block = metadata$Subject_ID, correlation = corfit$consensus)
    } else {
      fit <- lmFit(v, design)
    }
    
    fit <- eBayes(fit)
    
    # Extract interaction effect
    interaction_coef <- grep("Sex.*Region|Region.*Sex", colnames(design), value = TRUE)
    if (length(interaction_coef) > 0) {
      results <- topTable(fit, coef = interaction_coef[1], number = Inf, sort.by = "P")
      results$Dataset <- dataset_name
      interaction_results[[dataset_name]] <- results
      
      n_sig <- sum(results$adj.P.Val < 0.05, na.rm = TRUE)
      cat("  Sex × Region interaction:", n_sig, "significant genes (FDR < 0.05)\n")
    }
  }
  
  return(interaction_results)
}

# ==============================================================================
# SEX-STRATIFIED REGIONAL ANALYSIS
# ==============================================================================

#' Perform region comparisons stratified by sex
analyze_sex_stratified_regions <- function(data) {
  cat("\n=== Sex-Stratified Regional Analysis ===\n")
  
  stratified_results <- list()
  
  for (dataset_name in names(data)) {
    cat("\nAnalyzing", dataset_name, "\n")
    
    expr_data <- data[[dataset_name]]
    counts <- expr_data$filtered_counts
    metadata <- expr_data$metadata
    
    dataset_results <- list()
    
    # Analyze each sex separately
    for (sex in unique(metadata$Sex)) {
      cat("  Sex:", sex, "\n")
      
      # Filter to specific sex
      sex_mask <- metadata$Sex == sex
      sex_counts <- counts[, sex_mask]
      sex_metadata <- metadata[sex_mask, ]
      
      # Skip if too few samples
      if (ncol(sex_counts) < 6) {
        cat("    Skipping - too few samples:", ncol(sex_counts), "\n")
        next
      }
      
      # Create DGEList
      dge <- DGEList(counts = sex_counts)
      keep <- filterByExpr(dge, group = sex_metadata$Region)
      dge <- dge[keep, , keep.lib.sizes = FALSE]
      dge <- calcNormFactors(dge)
      
      # Design matrix
      design <- model.matrix(~ Region + OUD_Status, data = sex_metadata)
      
      # Voom transformation
      v <- voom(dge, design, plot = FALSE)
      
      # Account for subject pairing
      if ("Subject_ID" %in% colnames(sex_metadata)) {
        corfit <- duplicateCorrelation(v, design, block = sex_metadata$Subject_ID)
        fit <- lmFit(v, design, block = sex_metadata$Subject_ID, correlation = corfit$consensus)
      } else {
        fit <- lmFit(v, design)
      }
      
      fit <- eBayes(fit)
      
      # Extract region effect
      region_coef <- grep("Region", colnames(design), value = TRUE)[1]
      if (length(region_coef) > 0) {
        results <- topTable(fit, coef = region_coef, number = Inf, sort.by = "P")
        results$Sex <- sex
        results$Dataset <- dataset_name
        
        n_sig <- sum(results$adj.P.Val < 0.05, na.rm = TRUE)
        cat("    Regional differences:", n_sig, "significant genes (FDR < 0.05)\n")
        
        dataset_results[[paste0("Sex_", sex)]] <- results
      }
    }
    
    stratified_results[[dataset_name]] <- dataset_results
  }
  
  return(stratified_results)
}

# ==============================================================================
# SEX-SPECIFIC VISUALIZATION
# ==============================================================================

#' Create sex-specific visualizations
create_sex_visualizations <- function(sex_main_effects, sex_interactions, stratified_results) {
  cat("\n=== Creating Sex-Specific Visualizations ===\n")
  
  # 1. Sex main effects comparison
  if (length(sex_main_effects) >= 2) {
    gse174409_sex <- sex_main_effects$gse174409
    gse225158_sex <- sex_main_effects$gse225158
    
    # Compare significant genes between datasets
    gse174409_sig <- rownames(gse174409_sex)[gse174409_sex$adj.P.Val < 0.05 & !is.na(gse174409_sex$adj.P.Val)]
    gse225158_sig <- rownames(gse225158_sex)[gse225158_sex$adj.P.Val < 0.05 & !is.na(gse225158_sex$adj.P.Val)]
    
    # Venn diagram for sex effects
    png(file.path(SEX_ANALYSIS_DIR, "sex_effects_venn.png"), width = 800, height = 600)
    venn.plot <- venn.diagram(
      x = list(
        GSE174409_Sex = gse174409_sig,
        GSE225158_Sex = gse225158_sig
      ),
      category.names = c("GSE174409\n(Sex Effects)", "GSE225158\n(Sex Effects)"),
      filename = NULL,
      fill = c("#E74C3C", "#9B59B6"),
      alpha = 0.7,
      main = "Sex Main Effects: Cross-Dataset Overlap"
    )
    grid.draw(venn.plot)
    dev.off()
    
    cat("✓ Sex effects Venn diagram saved\n")
  }
  
  # 2. Sex-stratified results comparison
  if (length(stratified_results) >= 2) {
    # Create summary plot
    summary_data <- data.frame()
    
    for (dataset in names(stratified_results)) {
      for (sex_group in names(stratified_results[[dataset]])) {
        results <- stratified_results[[dataset]][[sex_group]]
        n_sig <- sum(results$adj.P.Val < 0.05, na.rm = TRUE)
        
        summary_data <- rbind(summary_data, data.frame(
          Dataset = dataset,
          Sex = gsub("Sex_", "", sex_group),
          N_Significant = n_sig
        ))
      }
    }
    
    if (nrow(summary_data) > 0) {
      p_sex_stratified <- ggplot(summary_data, aes(x = Sex, y = N_Significant, fill = Dataset)) +
        geom_col(position = "dodge", alpha = 0.8) +
        labs(title = "Sex-Stratified Regional Differences",
             subtitle = "Number of significant genes (region effect within each sex)",
             x = "Sex", y = "Significant Genes (FDR < 0.05)") +
        theme_minimal() +
        scale_fill_manual(values = c("gse174409" = "#3498DB", "gse225158" = "#E67E22"))
      
      ggsave(file.path(SEX_ANALYSIS_DIR, "sex_stratified_regional_differences.png"),
             p_sex_stratified, width = 10, height = 6, dpi = 300)
      
      cat("✓ Sex-stratified comparison plot saved\n")
    }
  }
  
  # 3. Sex interaction effects summary
  if (length(sex_interactions) >= 2) {
    interaction_summary <- data.frame(
      Dataset = names(sex_interactions),
      N_Significant_Interactions = sapply(sex_interactions, function(x) {
        sum(x$adj.P.Val < 0.05, na.rm = TRUE)
      })
    )
    
    p_interactions <- ggplot(interaction_summary, aes(x = Dataset, y = N_Significant_Interactions)) +
      geom_col(fill = "#F39C12", alpha = 0.8) +
      geom_text(aes(label = N_Significant_Interactions), vjust = -0.5) +
      labs(title = "Sex × Region Interaction Effects",
           subtitle = "Genes showing differential regional effects by sex",
           x = "Dataset", y = "Significant Interactions (FDR < 0.05)") +
      theme_minimal()
    
    ggsave(file.path(SEX_ANALYSIS_DIR, "sex_region_interactions.png"),
           p_interactions, width = 8, height = 6, dpi = 300)
    
    cat("✓ Sex interaction plot saved\n")
  }
}

# ==============================================================================
# SEX-SPECIFIC RESULTS EXPORT
# ==============================================================================

#' Export sex-specific analysis results
export_sex_results <- function(sex_main_effects, sex_interactions, stratified_results) {
  cat("\n=== Exporting Sex-Specific Results ===\n")
  
  # Export sex main effects
  for (dataset in names(sex_main_effects)) {
    filename <- file.path(SEX_ANALYSIS_DIR, paste0(dataset, "_sex_main_effects.csv"))
    write.csv(sex_main_effects[[dataset]], filename, row.names = TRUE)
    cat("✓ Sex main effects saved:", basename(filename), "\n")
  }
  
  # Export sex interactions
  for (dataset in names(sex_interactions)) {
    filename <- file.path(SEX_ANALYSIS_DIR, paste0(dataset, "_sex_region_interactions.csv"))
    write.csv(sex_interactions[[dataset]], filename, row.names = TRUE)
    cat("✓ Sex interactions saved:", basename(filename), "\n")
  }
  
  # Export stratified results
  for (dataset in names(stratified_results)) {
    for (sex_group in names(stratified_results[[dataset]])) {
      filename <- file.path(SEX_ANALYSIS_DIR, paste0(dataset, "_", sex_group, "_regional_effects.csv"))
      write.csv(stratified_results[[dataset]][[sex_group]], filename, row.names = TRUE)
      cat("✓ Stratified results saved:", basename(filename), "\n")
    }
  }
  
  # Create comprehensive summary
  sex_summary <- data.frame(
    Analysis_Type = c("Sex Main Effects", "Sex × Region Interactions", "Male Regional Effects", "Female Regional Effects"),
    GSE174409_Significant = c(
      if ("gse174409" %in% names(sex_main_effects)) sum(sex_main_effects$gse174409$adj.P.Val < 0.05, na.rm = TRUE) else 0,
      if ("gse174409" %in% names(sex_interactions)) sum(sex_interactions$gse174409$adj.P.Val < 0.05, na.rm = TRUE) else 0,
      if ("gse174409" %in% names(stratified_results) && "Sex_Male" %in% names(stratified_results$gse174409)) sum(stratified_results$gse174409$Sex_Male$adj.P.Val < 0.05, na.rm = TRUE) else 0,
      if ("gse174409" %in% names(stratified_results) && "Sex_Female" %in% names(stratified_results$gse174409)) sum(stratified_results$gse174409$Sex_Female$adj.P.Val < 0.05, na.rm = TRUE) else 0
    ),
    GSE225158_Significant = c(
      if ("gse225158" %in% names(sex_main_effects)) sum(sex_main_effects$gse225158$adj.P.Val < 0.05, na.rm = TRUE) else 0,
      if ("gse225158" %in% names(sex_interactions)) sum(sex_interactions$gse225158$adj.P.Val < 0.05, na.rm = TRUE) else 0,
      if ("gse225158" %in% names(stratified_results) && "Sex_Male" %in% names(stratified_results$gse225158)) sum(stratified_results$gse225158$Sex_Male$adj.P.Val < 0.05, na.rm = TRUE) else 0,
      if ("gse225158" %in% names(stratified_results) && "Sex_Female" %in% names(stratified_results$gse225158)) sum(stratified_results$gse225158$Sex_Female$adj.P.Val < 0.05, na.rm = TRUE) else 0
    )
  )
  
  write.csv(sex_summary, file.path(SEX_ANALYSIS_DIR, "sex_analysis_summary.csv"), row.names = FALSE)
  cat("✓ Sex analysis summary saved\n")
  
  return(sex_summary)
}

# ==============================================================================
# MAIN SEX ANALYSIS PIPELINE
# ==============================================================================

#' Run comprehensive sex-specific analysis
run_sex_specific_analysis <- function() {
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("COMPREHENSIVE SEX-SPECIFIC NEUROINFLAMMATORY ANALYSIS\n")
  cat("Cross-dataset sex effects, interactions, and stratified comparisons\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  tryCatch({
    # Load data
    data <- load_sex_analysis_data()
    
    # Run analyses
    sex_main_effects <- analyze_sex_main_effects(data)
    sex_interactions <- analyze_sex_region_interactions(data)
    stratified_results <- analyze_sex_stratified_regions(data)
    
    # Create visualizations
    create_sex_visualizations(sex_main_effects, sex_interactions, stratified_results)
    
    # Export results
    sex_summary <- export_sex_results(sex_main_effects, sex_interactions, stratified_results)
    
    cat("\n", paste(rep("=", 70), collapse = ""), "\n")
    cat("SUCCESS: Sex-specific analysis complete!\n")
    cat("✓ Sex main effects analyzed\n")
    cat("✓ Sex × Region interactions tested\n")
    cat("✓ Sex-stratified regional comparisons performed\n")
    cat("✓ Results exported with visualizations\n")
    cat("Results saved to:", SEX_ANALYSIS_DIR, "\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")
    
    cat("\n=== SEX ANALYSIS SUMMARY ===\n")
    print(sex_summary)
    
    return(list(
      sex_main_effects = sex_main_effects,
      sex_interactions = sex_interactions,
      stratified_results = stratified_results,
      summary = sex_summary
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
  sex_analysis_results <- run_sex_specific_analysis()
}
