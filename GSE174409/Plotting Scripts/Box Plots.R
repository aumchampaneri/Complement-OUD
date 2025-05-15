# Complement_Gene_Plots.R
# Script with improved gene symbol matching and enhanced plotting aesthetics

# ----------------------
# 1. SETUP AND DATA LOADING
# ----------------------

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggpubr)
  library(RColorBrewer)
  library(ggbeeswarm)
  library(patchwork)
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
  exact_matches <- mapping$EnsemblID[grep(paste0("^", symbol, "$"), mapping$Symbol)]
  if (length(exact_matches) > 0 && !is.na(exact_matches[1])) {
    return(exact_matches[1])
  }

  complement_map <- list(
    "C1QA" = "C1QA",
    "C3" = "C3",
    "C4A" = "C4A",
    "C5AR1" = "C5AR1",
    "CFH" = "CFH",
    "C3AR1" = "C3AR1"
  )

  if (symbol %in% names(complement_map)) {
    exact_symbol <- complement_map[[symbol]]
    exact_matches <- mapping$EnsemblID[mapping$Symbol == exact_symbol]
    if (length(exact_matches) > 0 && !is.na(exact_matches[1])) {
      return(exact_matches[1])
    }
  }

  return(NULL)
}

# ----------------------
# 3. PLOTTING FUNCTION
# ----------------------

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

  ensembl_id <- find_gene_id(gene_symbol, gene_mapping)
  if (is.null(ensembl_id)) {
    cat("Could not find Ensembl ID for", gene_symbol, "\n")
    return(NULL)
  }

  if (!(ensembl_id %in% rownames(logcpm_filtered_norm))) {
    cat("Gene ID", ensembl_id, "not found in expression data\n")
    return(NULL)
  }

  df <- data.frame(
    Expression = logcpm_filtered_norm[ensembl_id, ],
    Sex = metadata_df$sex,
    Diagnosis = metadata_df$diagnosis,
    Region = metadata_df$region
  )

  if (!facet_by_region && !is.null(region)) {
    df <- df[df$Region == region, ]
    if (nrow(df) == 0) {
      cat("No data for gene", gene_symbol, "in region", region, "\n")
      return(NULL)
    }
  }

  y_max <- max(df$Expression, na.rm = TRUE)

  counts_df <- df %>%
    group_by(Sex, Diagnosis) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(label = paste0("n=", n))

  p <- ggplot(df, aes(x = Sex, y = Expression, fill = Diagnosis))

  if (plot_type == "violin") {
    p <- p + geom_violin(alpha = 0.7, trim = FALSE)
    if (add_points) {
      p <- p + geom_quasirandom(dodge.width = 0.75, alpha = 0.6, size = 1)
    }
  } else {
    p <- p + geom_boxplot(alpha = 0.7, outlier.shape = NA)
    if (add_points) {
      p <- p + geom_quasirandom(dodge.width = 0.75, alpha = 0.6, size = 1)
    }
  }

  if (facet_by_region) {
    p <- p + facet_wrap(~Region)
    title <- paste("Expression of", gene_symbol, "by Sex and Diagnosis")
  } else if (!is.null(region)) {
    title <- paste("Expression of", gene_symbol, "by Sex and Diagnosis in", region)
  } else {
    title <- paste("Expression of", gene_symbol, "by Sex and Diagnosis")
  }

  p <- p +
    labs(title = title, y = "Log2 CPM Expression") +
    scale_fill_brewer(palette = "Set1") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.text = element_text(size = 12),
      axis.title.y = element_text(size = 14),
      legend.title = element_blank(),
      legend.position = "bottom",
      strip.text = element_text(size = 12, face = "bold")
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    geom_text(data = counts_df, aes(x = Sex, y = -1, label = label), inherit.aes = FALSE) +
    stat_compare_means(
      aes(group = Diagnosis),
      method = "wilcox.test",
      label = "p.signif",
      label.y = y_max + 0.5
    )

  return(p)
}

# ----------------------
# 4. GENERATE PLOTS FOR COMPLEMENT GENES
# ----------------------

complement_genes <- c("CFH", "C3AR1", "C1QA", "C1QB", "C3", "C4A", "C5", "C5AR1")

valid_plots <- list()
for (gene in complement_genes) {
  p <- plot_gene_by_sex_diagnosis(gene, plot_type = "boxplot", facet_by_region = TRUE)
  if (!is.null(p)) {
    valid_plots[[gene]] <- p
    ggsave(
      file.path(complement_plot_dir, paste0(gene, "_boxplot_faceted.png")),
      p,
      width = 8,
      height = 6,
      dpi = 300,
      device = "png"
    )
    cat("Created boxplot for", gene, "\n")
  }
}

if (length(valid_plots) > 0) {
  ncol_val <- min(2, length(valid_plots))
  combined_plot <- wrap_plots(valid_plots, ncol = ncol_val) +
    plot_annotation(title = "Complement Gene Expression by Sex and Diagnosis")
  ggsave(
    file.path(complement_plot_dir, "complement_genes_combined.pdf"),
    combined_plot,
    width = 8 * ncol_val,
    height = 6 * ceiling(length(valid_plots) / ncol_val),
    dpi = 300,
    device = "pdf"
  )
  cat("Created combined plot with", length(valid_plots), "complement genes\n")
}

cat("✓ Complement gene plots generated successfully. Output saved to:", complement_plot_dir, "\n")