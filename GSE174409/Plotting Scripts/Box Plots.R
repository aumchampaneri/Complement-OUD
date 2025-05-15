# Complement_Gene_Plots.R
# Script with improved gene symbol matching for complement genes

# ----------------------
# 1. SETUP AND DATA LOADING
# ----------------------

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggpubr)
  library(RColorBrewer)
  library(gridExtra)
})

# Set directories
base_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE174409"
qc_dir <- file.path(base_dir, "QC")
results_dir <- file.path(base_dir, "NeuroinflammationResults")
complement_plot_dir <- file.path(base_dir, "ComplementGenePlots")
dir.create(complement_plot_dir, showWarnings = FALSE)

# Load preprocessed data
logcpm_filtered_norm <- readRDS(file.path(qc_dir, "logcpm_filtered_normalized.rds"))
metadata_df <- readRDS(file.path(qc_dir, "metadata.rds"))

# Load gene mapping
gene_mapping <- read.csv(file.path(results_dir, "ensembl_to_symbol_mapping.csv"))

# ----------------------
# 2. IMPROVED GENE SYMBOL MATCHING
# ----------------------

#' Match complement genes correctly
#' @param symbol Character: gene symbol to find
#' @param mapping Data frame: gene mapping data
#' @return Character: matching Ensembl ID or NULL if not found
find_gene_id <- function(symbol, mapping) {
  # Exact match with word boundary (to prevent partial matches)
  exact_matches <- mapping$EnsemblID[grep(paste0("^", symbol, "$"), mapping$Symbol)]

  if(length(exact_matches) > 0 && !is.na(exact_matches[1])) {
    if(length(exact_matches) > 1) {
      cat("Multiple exact matches found for", symbol, "- using first match\n")
    }
    cat("Found exact match for", symbol, "\n")
    return(exact_matches[1])
  }

  # Manual mapping for complement genes (to ensure correct matching)
  complement_map <- list(
    "C1QA" = "C1QA",
    "C3" = "C3",          # Exact match for complement component 3
    "C4A" = "C4A",         # Exact match for complement component 4A
    "C5AR1" = "C5AR1",
    "CFH" = "CFH",
    "C3AR1" = "C3AR1"
  )

  if(symbol %in% names(complement_map)) {
    exact_symbol <- complement_map[[symbol]]
    exact_matches <- mapping$EnsemblID[mapping$Symbol == exact_symbol]

    if(length(exact_matches) > 0 && !is.na(exact_matches[1])) {
      cat("Found complement gene", symbol, "using manual mapping\n")
      return(exact_matches[1])
    }
  }

  # For short gene names, be very strict with matching
  if(nchar(symbol) <= 3) {
    # Try to match with word boundaries to prevent substring matches
    bounded_matches <- mapping$EnsemblID[grep(paste0("\\b", symbol, "\\b"), mapping$Symbol)]
    if(length(bounded_matches) > 0 && !is.na(bounded_matches[1])) {
      matched_symbol <- mapping$Symbol[grep(paste0("\\b", symbol, "\\b"), mapping$Symbol)][1]
      cat("Found short gene name", symbol, "with strict matching:", matched_symbol, "\n")
      return(bounded_matches[1])
    }
  }

  cat("Could not find match for", symbol, "- checking all symbols in mapping data\n")

  # Print all symbols that contain this gene name to help diagnose
  contains_matches <- grep(symbol, mapping$Symbol, fixed=TRUE)
  if(length(contains_matches) > 0) {
    cat("  Symbols containing '", symbol, "': ",
        paste(head(mapping$Symbol[contains_matches], 10), collapse=", "), "\n")
  }

  return(NULL)
}

#' Plot gene expression by sex and diagnosis
#' @param gene_symbol Character: gene symbol to plot
#' @param plot_type Character: "boxplot" or "violin"
#' @param facet_by_region Logical: whether to facet by region
#' @param region Character: specific region to plot if not faceting
#' @param add_points Logical: whether to add individual data points
#' @return A ggplot object or NULL if gene not found
plot_gene_by_sex_diagnosis <- function(gene_symbol,
                                      plot_type = "boxplot",
                                      facet_by_region = TRUE,
                                      region = NULL,
                                      add_points = TRUE) {

  # Find Ensembl ID with improved matching
  ensembl_id <- find_gene_id(gene_symbol, gene_mapping)

  if(is.null(ensembl_id)) {
    cat("Could not find Ensembl ID for", gene_symbol, "\n")
    return(NULL)
  }

  # Check if gene is in expression data
  if(!(ensembl_id %in% rownames(logcpm_filtered_norm))) {
    cat("Gene ID", ensembl_id, "not found in expression data\n")
    return(NULL)
  }

  # Get the actual symbol used in the mapping for the title
  actual_symbol <- gene_mapping$Symbol[gene_mapping$EnsemblID == ensembl_id][1]

  # Create data frame for plotting
  df <- data.frame(
    Expression = logcpm_filtered_norm[ensembl_id, ],
    Sex = metadata_df$sex,
    Diagnosis = metadata_df$diagnosis,
    Region = metadata_df$region
  )

  # Filter by region if specified and not faceting
  if(!facet_by_region && !is.null(region)) {
    df <- df[df$Region == region, ]
    if(nrow(df) == 0) {
      cat("No data for gene", gene_symbol, "in region", region, "\n")
      return(NULL)
    }
  }

  # Create plot
  p <- ggplot(df, aes(x = Sex, y = Expression, fill = Diagnosis))

  # Add appropriate geom based on plot type
  if(plot_type == "violin") {
    p <- p + geom_violin(alpha = 0.7, trim = FALSE)
    if(add_points) {
      p <- p + geom_jitter(width = 0.1, height = 0, alpha = 0.6, size = 1)
    }
  } else {
    p <- p + geom_boxplot(alpha = 0.7, outlier.shape = NA)
    if(add_points) {
      p <- p + geom_jitter(width = 0.2, height = 0, alpha = 0.6, size = 1)
    }
  }

  # Add faceting if requested
  if(facet_by_region) {
    p <- p + facet_wrap(~Region)
    title <- paste("Expression of", gene_symbol, "by Sex and Diagnosis")
  } else if(!is.null(region)) {
    title <- paste("Expression of", gene_symbol, "by Sex and Diagnosis in", region)
  } else {
    title <- paste("Expression of", gene_symbol, "by Sex and Diagnosis")
  }

  # Add styling
  p <- p +
    labs(title = title, y = "Log2 CPM Expression") +
    scale_fill_manual(values = c("CONT" = "blue", "OUD" = "red")) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title.x = element_blank(),
      legend.position = "bottom"
    )

  # Add statistics comparing OUD vs CONT within each sex
  p <- p + stat_compare_means(
    aes(group = Diagnosis),
    method = "t.test",
    label = "p.signif",
    label.y.npc = 0.95
  )

  return(p)
}

# ----------------------
# 3. ANALYZE GENE MAPPING DATA
# ----------------------

# Print information about the gene mapping data
cat("Total gene mappings available:", nrow(gene_mapping), "\n")
cat("Number of unique gene symbols:", length(unique(gene_mapping$Symbol)), "\n")
cat("Number of genes with missing symbols:", sum(is.na(gene_mapping$Symbol)), "\n")

# Look for exact matches of complement genes
complement_genes <- c("CFH", "C3AR1", "C1QA", "C3", "C4A", "C5AR1")
cat("\nChecking for exact matches of complement genes:\n")
for(gene in complement_genes) {
  exact_matches <- gene_mapping$Symbol[gene_mapping$Symbol == gene]
  if(length(exact_matches) > 0) {
    cat("✓", gene, "found with exact match\n")

    # Get the Ensembl ID and check if it's in the expression data
    ensembl_id <- gene_mapping$EnsemblID[gene_mapping$Symbol == gene][1]
    if(ensembl_id %in% rownames(logcpm_filtered_norm)) {
      cat("  ✓ Gene is in expression data\n")
    } else {
      cat("  ✗ Gene not found in expression data\n")
    }
  } else {
    cat("✗", gene, "not found with exact match\n")
  }
}

# ----------------------
# 4. GENERATE PLOTS FOR COMPLEMENT GENES
# ----------------------

# Plot each gene with boxplot
for(gene in complement_genes) {
  box_facet <- plot_gene_by_sex_diagnosis(gene, plot_type = "boxplot", facet_by_region = TRUE)
  if(!is.null(box_facet)) {
    ggsave(
      file.path(complement_plot_dir, paste0(gene, "_boxplot_faceted.png")),
      box_facet,
      width = 8,
      height = 6,
      dpi = 300
    )
    cat("Created boxplot for", gene, "\n")
  }
}

# Create a list of valid plots
valid_plots <- list()
for(gene in complement_genes) {
  p <- plot_gene_by_sex_diagnosis(gene, plot_type = "boxplot", facet_by_region = TRUE)
  if(!is.null(p)) {
    valid_plots[[gene]] <- p
  }
}

# Create a combined plot of all valid genes
if(length(valid_plots) > 0) {
  ncol_val <- min(2, length(valid_plots))
  combined_plot <- gridExtra::grid.arrange(grobs = valid_plots, ncol = ncol_val)
  ggsave(
    file.path(complement_plot_dir, "complement_genes_combined.png"),
    combined_plot,
    width = 8 * ncol_val,
    height = 6 * ceiling(length(valid_plots)/ncol_val),
    dpi = 300
  )
  cat("Created combined plot with", length(valid_plots), "complement genes\n")
}

cat("✓ Complement gene plots generated successfully. Output saved to:", complement_plot_dir, "\n")