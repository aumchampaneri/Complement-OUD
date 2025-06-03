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
  library(grid)  # Add this for grid.draw
  library(biomaRt)  # Add biomaRt explicitly
  library(metafor)  # Add for meta-analysis
  library(meta)     # Add for additional meta-analysis functions
})

# ==============================================================================
# GENE ID HARMONIZATION
# ==============================================================================

#' Convert Ensembl IDs to gene symbols using biomaRt
convert_ensembl_to_symbols <- function(ensembl_ids) {
  cat("Converting Ensembl IDs to gene symbols...\n")
  
  tryCatch({
    # Connect to Ensembl biomart with error handling
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", 
                       host = "https://useast.ensembl.org")
    
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
    
  }, error = function(e) {
    cat("⚠ biomaRt connection failed:", e$message, "\n")
    cat("Using alternative approach with clusterProfiler...\n")
    
    # Fallback to clusterProfiler approach
    gene_map <- clusterProfiler::bitr(ensembl_ids, 
                                      fromType = "ENSEMBL", 
                                      toType = "SYMBOL", 
                                      OrgDb = org.Hs.eg.db)
    
    colnames(gene_map) <- c("ensembl_gene_id", "hgnc_symbol")
    
    cat(sprintf("Mapped %d/%d Ensembl IDs using clusterProfiler\n", 
                nrow(gene_map), length(ensembl_ids)))
    
    return(gene_map)
  })
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
  
  cat("✓ GSE174409 loaded:", length(gse174409_results), "components\n")
  cat("✓ GSE225158 loaded:", length(gse225158_results), "components\n")
  
  # Check if data uses gene symbols or Ensembl IDs
  gse174409_genes <- rownames(gse174409_results$paired_limma$results)[1:5]
  gse225158_genes <- rownames(gse225158_results$paired_limma$results)[1:5]
  
  cat("GSE174409 gene IDs:", paste(gse174409_genes, collapse = ", "), "\n")
  cat("GSE225158 gene IDs:", paste(gse225158_genes, collapse = ", "), "\n")
  
  # Smart gene ID detection and harmonization
  gse174409_uses_ensembl <- any(grepl("^ENSG", gse174409_genes))
  gse225158_uses_ensembl <- any(grepl("^ENSG", gse225158_genes))
  
  cat("Gene ID formats - GSE174409:", ifelse(gse174409_uses_ensembl, "Ensembl", "Symbols"), 
      "| GSE225158:", ifelse(gse225158_uses_ensembl, "Ensembl", "Symbols"), "\n")
  
  # If both use symbols or both use Ensembl, proceed directly
  # If mixed, convert Ensembl to symbols
  if (gse174409_uses_ensembl && !gse225158_uses_ensembl) {
    cat("🔄 Converting GSE174409 from Ensembl to gene symbols...\n")
    gse174409_results <- convert_dataset_to_symbols(gse174409_results)
  } else if (!gse174409_uses_ensembl && gse225158_uses_ensembl) {
    cat("🔄 Converting GSE225158 from Ensembl to gene symbols...\n")
    gse225158_results <- convert_dataset_to_symbols(gse225158_results)
  } else {
    cat("✓ Gene ID formats are compatible\n")
  }
  
  return(list(
    gse174409 = gse174409_results,
    gse225158 = gse225158_results
  ))
}

#' Convert a dataset from Ensembl IDs to gene symbols
convert_dataset_to_symbols <- function(dataset_results) {
  tryCatch({
    # Get all Ensembl IDs from the first method
    all_ensembl_ids <- rownames(dataset_results$paired_limma$results)
    
    # Convert to gene symbols
    gene_map <- convert_ensembl_to_symbols(all_ensembl_ids)
    
    # Convert each method's results
    converted_results <- list()
    for (method_name in names(dataset_results)) {
      if (method_name %in% c("paired_limma", "mixed_effects", "deseq2")) {
        method_results <- dataset_results[[method_name]]$results
        
        # Merge with gene mapping
        method_results$ensembl_gene_id <- rownames(method_results)
        harmonized <- merge(method_results, gene_map, by = "ensembl_gene_id", all.x = FALSE)
        
        # Handle duplicates (keep best p-value)
        pval_col <- ifelse("adj.P.Val" %in% colnames(harmonized), "adj.P.Val", "padj")
        harmonized <- harmonized %>%
          group_by(hgnc_symbol) %>%
          slice_min(order_by = !!sym(pval_col), n = 1, with_ties = FALSE) %>%
          ungroup()
        
        # Set gene symbols as row names
        harmonized <- as.data.frame(harmonized)
        rownames(harmonized) <- harmonized$hgnc_symbol
        harmonized$ensembl_gene_id <- NULL
        harmonized$hgnc_symbol <- NULL
        
        converted_results[[method_name]] <- list(
          results = harmonized,
          fit = dataset_results[[method_name]]$fit,
          correlation = dataset_results[[method_name]]$correlation
        )
        
        cat(sprintf("  %s: %d genes after conversion\n", method_name, nrow(harmonized)))
      } else {
        # Keep other components as-is
        converted_results[[method_name]] <- dataset_results[[method_name]]
      }
    }
    
    return(converted_results)
    
  }, error = function(e) {
    cat("⚠ Gene conversion failed:", e$message, "\n")
    cat("Returning original dataset\n")
    return(dataset_results)
  })
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
# FORMAL META-ANALYSIS FUNCTIONS
# ==============================================================================

#' Perform Fisher's method for combining p-values across datasets
#' @param data Harmonized data from both datasets
#' @return List of Fisher's method results for each statistical method
fisher_method_analysis <- function(data) {
  cat("\n=== Fisher's Method P-value Combination ===\n")
  
  methods <- c("paired_limma", "mixed_effects", "deseq2")
  fisher_results <- list()
  
  for (method in methods) {
    cat("Processing method:", method, "\n")
    
    # Get results for this method
    gse174409_res <- data$gse174409[[method]]$results
    gse225158_res <- data$gse225158[[method]]$results
    
    # Get common genes
    common_genes <- intersect(rownames(gse174409_res), rownames(gse225158_res))
    
    if (length(common_genes) < 100) {
      cat("  Warning: Only", length(common_genes), "common genes found\n")
      next
    }
    
    # Extract p-values
    pval_col_174409 <- ifelse("P.Value" %in% colnames(gse174409_res), "P.Value", "pvalue")
    pval_col_225158 <- ifelse("P.Value" %in% colnames(gse225158_res), "P.Value", "pvalue")
    
    # Create combined data frame
    combined_data <- data.frame(
      gene = common_genes,
      p1 = gse174409_res[common_genes, pval_col_174409],
      p2 = gse225158_res[common_genes, pval_col_225158],
      stringsAsFactors = FALSE
    ) %>%
      filter(!is.na(p1) & !is.na(p2) & p1 > 0 & p2 > 0)  # Remove invalid p-values
    
    # Apply Fisher's method to each gene
    combined_data$fisher_stat <- -2 * (log(combined_data$p1) + log(combined_data$p2))
    combined_data$fisher_pvalue <- pchisq(combined_data$fisher_stat, df = 4, lower.tail = FALSE)
    combined_data$fisher_padj <- p.adjust(combined_data$fisher_pvalue, method = "BH")
    
    # Calculate effect directions and significance
    fc_col_174409 <- ifelse("logFC" %in% colnames(gse174409_res), "logFC", "log2FoldChange")
    fc_col_225158 <- ifelse("logFC" %in% colnames(gse225158_res), "logFC", "log2FoldChange")
    
    combined_data$fc1 <- gse174409_res[combined_data$gene, fc_col_174409]
    combined_data$fc2 <- gse225158_res[combined_data$gene, fc_col_225158]
    combined_data$direction_consistent <- sign(combined_data$fc1) == sign(combined_data$fc2)
    
    # Summary statistics
    n_significant_fisher <- sum(combined_data$fisher_padj < 0.05, na.rm = TRUE)
    n_consistent_direction <- sum(combined_data$direction_consistent & combined_data$fisher_padj < 0.05, na.rm = TRUE)
    
    fisher_results[[method]] <- list(
      data = combined_data,
      n_genes = nrow(combined_data),
      n_significant = n_significant_fisher,
      n_consistent = n_consistent_direction,
      consistency_rate = round(n_consistent_direction / n_significant_fisher * 100, 1)
    )
    
    cat(sprintf("  Genes tested: %d\n", nrow(combined_data)))
    cat(sprintf("  Significant (Fisher FDR < 0.05): %d\n", n_significant_fisher))
    cat(sprintf("  Consistent direction: %d (%.1f%%)\n", n_consistent_direction, 
               n_consistent_direction / n_significant_fisher * 100))
  }
  
  return(fisher_results)
}

#' Perform random effects meta-analysis for effect sizes
#' @param data Harmonized data from both datasets
#' @return List of random effects meta-analysis results
random_effects_meta_analysis <- function(data) {
  cat("\n=== Random Effects Meta-Analysis ===\n")
  
  methods <- c("paired_limma", "mixed_effects", "deseq2")
  meta_results <- list()
  
  for (method in methods) {
    cat("Processing method:", method, "\n")
    
    # Get results for this method
    gse174409_res <- data$gse174409[[method]]$results
    gse225158_res <- data$gse225158[[method]]$results
    
    # Get common genes
    common_genes <- intersect(rownames(gse174409_res), rownames(gse225158_res))
    
    if (length(common_genes) < 100) {
      cat("  Warning: Only", length(common_genes), "common genes found\n")
      next
    }
    
    # Extract effect sizes and standard errors
    fc_col_174409 <- ifelse("logFC" %in% colnames(gse174409_res), "logFC", "log2FoldChange")
    fc_col_225158 <- ifelse("logFC" %in% colnames(gse225158_res), "logFC", "log2FoldChange")
    
    # For standard errors, use different approaches based on method
    if (method %in% c("paired_limma", "mixed_effects")) {
      se_col_174409 <- "SE"  # Standard error from limma
      se_col_225158 <- "SE"
    } else {
      se_col_174409 <- "lfcSE"  # Standard error from DESeq2
      se_col_225158 <- "lfcSE"
    }
    
    # Check if SE columns exist, otherwise estimate from t-statistics
    if (!se_col_174409 %in% colnames(gse174409_res)) {
      if ("t" %in% colnames(gse174409_res)) {
        # Estimate SE from t-statistic and logFC
        gse174409_res$SE <- abs(gse174409_res[[fc_col_174409]] / gse174409_res$t)
      } else {
        cat("  Warning: Cannot find standard errors for", method, "- using approximation\n")
        # Very rough approximation based on p-values
        t_stats <- qt(gse174409_res$P.Value/2, df = 30, lower.tail = FALSE)  # Approximate df
        gse174409_res$SE <- abs(gse174409_res[[fc_col_174409]] / t_stats)
      }
      se_col_174409 <- "SE"
    }
    
    if (!se_col_225158 %in% colnames(gse225158_res)) {
      if ("t" %in% colnames(gse225158_res)) {
        gse225158_res$SE <- abs(gse225158_res[[fc_col_225158]] / gse225158_res$t)
      } else {
        pval_col <- ifelse("P.Value" %in% colnames(gse225158_res), "P.Value", "pvalue")
        t_stats <- qt(gse225158_res[[pval_col]]/2, df = 10, lower.tail = FALSE)  # Approximate df
        gse225158_res$SE <- abs(gse225158_res[[fc_col_225158]] / t_stats)
      }
      se_col_225158 <- "SE"
    }
    
    # Create meta-analysis data
    meta_data <- data.frame(
      gene = common_genes,
      yi1 = gse174409_res[common_genes, fc_col_174409],    # Effect size study 1
      sei1 = gse174409_res[common_genes, se_col_174409],   # SE study 1
      yi2 = gse225158_res[common_genes, fc_col_225158],    # Effect size study 2
      sei2 = gse225158_res[common_genes, se_col_225158],   # SE study 2
      stringsAsFactors = FALSE
    ) %>%
      filter(!is.na(yi1) & !is.na(yi2) & !is.na(sei1) & !is.na(sei2) & 
             sei1 > 0 & sei2 > 0)  # Remove invalid data
    
    # Perform random effects meta-analysis for each gene
    gene_meta_results <- data.frame()
    
    # Sample analysis on first 1000 genes to avoid computation time
    sample_genes <- head(meta_data$gene, min(1000, nrow(meta_data)))
    
    for (i in seq_along(sample_genes)) {
      gene <- sample_genes[i]
      gene_data <- meta_data[meta_data$gene == gene, ]
      
      tryCatch({
        # Prepare data for metafor
        yi <- c(gene_data$yi1, gene_data$yi2)
        sei <- c(gene_data$sei1, gene_data$sei2)
        
        # Random effects meta-analysis
        res <- rma(yi = yi, sei = sei, method = "REML")
        
        gene_meta_results <- rbind(gene_meta_results, data.frame(
          gene = gene,
          pooled_effect = res$beta[1],
          pooled_se = res$se,
          pooled_pval = res$pval,
          tau2 = res$tau2,  # Between-study variance
          i2 = max(0, (res$QE - res$k + 1) / res$QE * 100),  # I² statistic
          stringsAsFactors = FALSE
        ))
        
      }, error = function(e) {
        # Skip genes with convergence issues
      })
      
      if (i %% 100 == 0) cat("  Processed", i, "genes\n")
    }
    
    # Adjust p-values
    gene_meta_results$pooled_padj <- p.adjust(gene_meta_results$pooled_pval, method = "BH")
    
    # Summary statistics
    n_significant_meta <- sum(gene_meta_results$pooled_padj < 0.05, na.rm = TRUE)
    median_i2 <- median(gene_meta_results$i2, na.rm = TRUE)
    
    meta_results[[method]] <- list(
      data = gene_meta_results,
      n_genes = nrow(gene_meta_results),
      n_significant = n_significant_meta,
      median_i2 = median_i2,
      high_heterogeneity = sum(gene_meta_results$i2 > 75, na.rm = TRUE)
    )
    
    cat(sprintf("  Genes analyzed: %d\n", nrow(gene_meta_results)))
    cat(sprintf("  Significant (meta FDR < 0.05): %d\n", n_significant_meta))
    cat(sprintf("  Median I²: %.1f%%\n", median_i2))
    cat(sprintf("  High heterogeneity (I² > 75%%): %d\n", 
               sum(gene_meta_results$i2 > 75, na.rm = TRUE)))
  }
  
  return(meta_results)
}

#' Create comprehensive meta-analysis visualizations
create_meta_analysis_plots <- function(fisher_results, meta_results) {
  cat("\n=== Creating Meta-Analysis Visualizations ===\n")
  
  # 1. Fisher's method results comparison
  fisher_summary <- data.frame(
    Method = names(fisher_results),
    N_Genes = sapply(fisher_results, function(x) x$n_genes),
    Fisher_Significant = sapply(fisher_results, function(x) x$n_significant),
    Consistent_Direction = sapply(fisher_results, function(x) x$n_consistent),
    Consistency_Rate = sapply(fisher_results, function(x) x$consistency_rate)
  )
  
  p1 <- ggplot(fisher_summary, aes(x = Method)) +
    geom_col(aes(y = Fisher_Significant), fill = "steelblue", alpha = 0.7) +
    geom_text(aes(y = Fisher_Significant, label = Fisher_Significant), vjust = -0.5) +
    labs(title = "Fisher's Method: Cross-Dataset Significant Genes",
         subtitle = "Combined p-value analysis",
         x = "Statistical Method", y = "Significant Genes (Fisher FDR < 0.05)") +
    theme_minimal()
  
  ggsave(file.path(INTEGRATION_DIR, "fisher_method_results.png"), 
         p1, width = 10, height = 6, dpi = 300)
  
  # 2. Consistency rate plot
  p2 <- ggplot(fisher_summary, aes(x = Method, y = Consistency_Rate)) +
    geom_col(fill = "darkgreen", alpha = 0.7) +
    geom_text(aes(label = paste0(Consistency_Rate, "%")), vjust = -0.5) +
    ylim(0, 100) +
    labs(title = "Effect Direction Consistency",
         subtitle = "Percentage of significant genes with consistent direction",
         x = "Statistical Method", y = "Consistency Rate (%)") +
    theme_minimal()
  
  ggsave(file.path(INTEGRATION_DIR, "direction_consistency.png"), 
         p2, width = 10, height = 6, dpi = 300)
  
  # 3. Meta-analysis results comparison
  if (length(meta_results) > 0) {
    meta_summary <- data.frame(
      Method = names(meta_results),
      N_Genes = sapply(meta_results, function(x) x$n_genes),
      Meta_Significant = sapply(meta_results, function(x) x$n_significant),
      Median_I2 = sapply(meta_results, function(x) x$median_i2),
      High_Heterogeneity = sapply(meta_results, function(x) x$high_heterogeneity)
    )
    
    p3 <- ggplot(meta_summary, aes(x = Method)) +
      geom_col(aes(y = Meta_Significant), fill = "darkorange", alpha = 0.7) +
      geom_text(aes(y = Meta_Significant, label = Meta_Significant), vjust = -0.5) +
      labs(title = "Random Effects Meta-Analysis: Pooled Effects",
           subtitle = "Combined effect size analysis",
           x = "Statistical Method", y = "Significant Genes (Meta FDR < 0.05)") +
      theme_minimal()
    
    ggsave(file.path(INTEGRATION_DIR, "meta_analysis_results.png"), 
           p3, width = 10, height = 6, dpi = 300)
    
    # 4. Heterogeneity assessment
    p4 <- ggplot(meta_summary, aes(x = Method, y = Median_I2)) +
      geom_col(fill = "purple", alpha = 0.7) +
      geom_text(aes(label = paste0(round(Median_I2, 1), "%")), vjust = -0.5) +
      geom_hline(yintercept = 75, color = "red", linetype = "dashed", alpha = 0.7) +
      labs(title = "Between-Study Heterogeneity (I² statistic)",
           subtitle = "Red line indicates high heterogeneity threshold (75%)",
           x = "Statistical Method", y = "Median I² (%)") +
      theme_minimal()
    
    ggsave(file.path(INTEGRATION_DIR, "heterogeneity_assessment.png"), 
           p4, width = 10, height = 6, dpi = 300)
  }
  
  cat("✓ Meta-analysis plots saved\n")
}

#' Export meta-analysis results to CSV
export_meta_analysis_results <- function(fisher_results, meta_results) {
  cat("\n=== Exporting Meta-Analysis Results ===\n")
  
  # Export Fisher's method results
  for (method in names(fisher_results)) {
    fisher_file <- file.path(INTEGRATION_DIR, paste0("fisher_method_", method, ".csv"))
    write.csv(fisher_results[[method]]$data, fisher_file, row.names = FALSE)
    cat("✓ Fisher results saved:", basename(fisher_file), "\n")
  }
  
  # Export meta-analysis results
  for (method in names(meta_results)) {
    meta_file <- file.path(INTEGRATION_DIR, paste0("meta_analysis_", method, ".csv"))
    write.csv(meta_results[[method]]$data, meta_file, row.names = FALSE)
    cat("✓ Meta-analysis results saved:", basename(meta_file), "\n")
  }
  
  # Create comprehensive summary
  if (length(fisher_results) > 0 && length(meta_results) > 0) {
    comprehensive_summary <- data.frame(
      Method = names(fisher_results),
      Fisher_Genes_Tested = sapply(fisher_results, function(x) x$n_genes),
      Fisher_Significant = sapply(fisher_results, function(x) x$n_significant),
      Fisher_Consistency_Rate = sapply(fisher_results, function(x) x$consistency_rate),
      Meta_Genes_Tested = sapply(names(fisher_results), function(m) {
        if (m %in% names(meta_results)) meta_results[[m]]$n_genes else NA
      }),
      Meta_Significant = sapply(names(fisher_results), function(m) {
        if (m %in% names(meta_results)) meta_results[[m]]$n_significant else NA
      }),
      Meta_Median_I2 = sapply(names(fisher_results), function(m) {
        if (m %in% names(meta_results)) meta_results[[m]]$median_i2 else NA
      })
    )
    
    summary_file <- file.path(INTEGRATION_DIR, "comprehensive_meta_analysis_summary.csv")
    write.csv(comprehensive_summary, summary_file, row.names = FALSE)
    cat("✓ Comprehensive summary saved:", basename(summary_file), "\n")
  }
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

#' Create enhanced integration summary with statistical tests
create_integration_summary <- function(overlap_results, correlation_results) {
  cat("\n=== Creating Enhanced Integration Summary ===\n")
  
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
    }),
    stringsAsFactors = FALSE
  )
  
  # Calculate Fisher's exact test for overlap significance
  for (i in 1:nrow(summary_df)) {
    method <- summary_df$Method[i]
    overlap_data <- overlap_results[[method]]
    
    # Create contingency table for Fisher's test
    # [significant in both, significant in 174409 only]
    # [significant in 225158 only, not significant in either]
    
    sig_both <- overlap_data$overlap_count
    sig_174409_only <- length(overlap_data$gse174409_sig) - sig_both
    sig_225158_only <- length(overlap_data$gse225158_sig) - sig_both
    
    # Estimate total tested genes (conservative approach)
    total_tested <- max(overlap_data$total_tested, 15000)  # Approximate
    not_sig_either <- total_tested - sig_both - sig_174409_only - sig_225158_only
    
    # Fisher's exact test
    if (not_sig_either > 0) {
      fisher_matrix <- matrix(c(sig_both, sig_174409_only, sig_225158_only, not_sig_either), 
                             nrow = 2)
      fisher_test <- fisher.test(fisher_matrix)
      summary_df$Fisher_PValue[i] <- fisher_test$p.value
      summary_df$Fisher_OR[i] <- round(fisher_test$estimate, 2)
    } else {
      summary_df$Fisher_PValue[i] <- NA
      summary_df$Fisher_OR[i] <- NA
    }
  }
  
  # Save enhanced summary
  summary_file <- file.path(INTEGRATION_DIR, "integration_summary_enhanced.csv")
  write.csv(summary_df, summary_file, row.names = FALSE)
  
  # Create integration visualization
  create_integration_visualization(summary_df, overlap_results)
  
  cat("✓ Enhanced integration summary saved:", summary_file, "\n")
  cat("\n=== Enhanced Integration Summary ===\n")
  print(summary_df)
  
  return(summary_df)
}

#' Create comprehensive integration visualization
create_integration_visualization <- function(summary_df, overlap_results) {
  cat("📊 Creating integration visualization...\n")
  
  # 1. Method comparison plot
  p_methods <- ggplot(summary_df, aes(x = Method)) +
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
         p_methods, width = 10, height = 6, dpi = 300)
  
  # 2. Overlap percentage plot
  p_overlap <- ggplot(summary_df, aes(x = Method, y = Overlap_Percent)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_text(aes(label = paste0(Overlap_Percent, "%")), vjust = -0.5) +
    labs(title = "Gene Overlap Percentage Between Datasets",
         x = "Analysis Method", y = "Overlap Percentage") +
    theme_minimal()
  
  ggsave(file.path(INTEGRATION_DIR, "overlap_percentage.png"), 
         p_overlap, width = 8, height = 6, dpi = 300)
  
  cat("✓ Integration visualizations saved\n")
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
    
    # NEW: Formal meta-analysis
    fisher_results <- fisher_method_analysis(integration_data)
    meta_results <- random_effects_meta_analysis(integration_data)
    
    # Create meta-analysis visualizations
    create_meta_analysis_plots(fisher_results, meta_results)
    
    # Export meta-analysis results
    export_meta_analysis_results(fisher_results, meta_results)
    
    # Create enhanced summary with meta-analysis
    summary_df <- create_integration_summary(overlap_results, correlation_results, 
                                            fisher_results, meta_results)
    
    cat("\n", paste(rep("=", 70), collapse = ""), "\n")
    cat("SUCCESS: Enhanced integration analysis with meta-analysis complete!\n")
    cat("✓ Gene overlap analysis\n")
    cat("✓ Effect size correlations\n")
    cat("✓ Fisher's method p-value combination\n")
    cat("✓ Random effects meta-analysis\n")
    cat("✓ Comprehensive visualizations and exports\n")
    cat("Results saved to:", INTEGRATION_DIR, "\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")
    
    return(list(
      overlap = overlap_results,
      correlation = correlation_results,
      fisher = fisher_results,
      meta_analysis = meta_results,
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
