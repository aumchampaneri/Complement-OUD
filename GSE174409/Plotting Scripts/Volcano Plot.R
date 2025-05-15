library(limma)
library(edgeR)
library(readr)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(stringr)
library(purrr)
library(gridExtra)
library(cowplot) # Add this to make sure we have the get_legend function

# ===== PART 1: Load Data =====
cat("==== 1. Loading data ====\n")
expr <- readRDS('/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE174409/QC/logcpm_filtered_normalized.rds')
meta <- readRDS('/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE174409/QC/metadata.rds')

# Find sample ID column
potential_id_cols <- c("title", "sample_name", "sample_id", "geo_accession", "run")
found_cols <- intersect(potential_id_cols, colnames(meta))
id_col <- found_cols[1]
meta$sample_id <- meta[[id_col]]

# Standardize column names
req_cols <- c("diagnosis", "sex", "region")
for (col in req_cols) {
  idx <- which(tolower(colnames(meta)) == col)
  if (length(idx) > 0) colnames(meta)[idx] <- col
}

# Match samples
common_samples <- intersect(colnames(expr), meta$sample_id)
cat("Matched samples:", length(common_samples), "\n")
expr_matched <- expr[, common_samples]
meta_matched <- meta[meta$sample_id %in% common_samples, ]
meta_matched <- meta_matched[match(colnames(expr_matched), meta_matched$sample_id), ]

# Load gene mappings and complement genes
gene_map <- read_csv('GSE174409/NeuroinflammationResults/ensembl_to_symbol_mapping.csv', show_col_types = FALSE)
ensembl_to_symbol <- setNames(gene_map$Symbol, gene_map$EnsemblID)
symbol_to_ensembl <- setNames(gene_map$EnsemblID, toupper(gene_map$Symbol))
comp_genes <- read_csv('/Users/aumchampaneri/PycharmProjects/Complement-OUD/Super Folder - GSE225158/GSE225158/KEGG outputs/kegg_complement_unique_genes.csv', show_col_types = FALSE)$gene
comp_genes <- toupper(comp_genes)
ensembl_ids <- symbol_to_ensembl[comp_genes]
ensembl_ids <- ensembl_ids[!is.na(ensembl_ids)]

# ===== PART 2: Perform DEG Analysis by Groups =====
cat("==== 2. Running DEG analysis for each sex and region ====\n")

# Create output directories
dir.create("GSE174409/NeuroinflammationResults/DEG_by_groups", showWarnings=FALSE, recursive=TRUE)
dir.create("GSE174409/figures/volcano_plots", showWarnings=FALSE, recursive=TRUE)

# Get unique sex and region values
sexes <- unique(meta_matched$sex)
regions <- unique(meta_matched$region)

# Create a function to run DEG analysis for a subset
run_deg_analysis <- function(sex_val, region_val) {
  # Filter data for this combination
  subset_idx <- meta_matched$sex == sex_val & meta_matched$region == region_val
  if (sum(subset_idx) < 6) {
    cat("Skipping", sex_val, region_val, "- too few samples (", sum(subset_idx), ")\n")
    return(NULL)
  }

  subset_meta <- meta_matched[subset_idx, ]
  subset_expr <- expr_matched[, subset_idx]

  # Check if we have enough samples per group
  subset_meta$diagnosis <- factor(subset_meta$diagnosis)
  group_counts <- table(subset_meta$diagnosis)
  if (min(group_counts) < 2) {
    cat("Skipping", sex_val, region_val, "- insufficient replicates in at least one group\n")
    return(NULL)
  }

  cat("Running DEG analysis for", sex_val, "in", region_val,
      "(n =", sum(subset_idx), ", groups:", paste(names(group_counts), collapse="/"), ")\n")

  # Create design matrix and run limma
  design <- model.matrix(~ diagnosis, data=subset_meta)
  fit <- lmFit(subset_expr, design)
  fit <- eBayes(fit)

  # Get results
  results <- topTable(fit, coef=2, number=Inf, sort.by="P") %>%
    rownames_to_column("gene_id")

  # Add gene symbols and flag complement genes
  results$gene_symbol <- ensembl_to_symbol[results$gene_id]
  results$is_complement <- results$gene_id %in% ensembl_ids

  # Replace NA symbols with gene_id
  results$gene_symbol[is.na(results$gene_symbol)] <- results$gene_id[is.na(results$gene_symbol)]

  # Save results
  filename <- paste0("GSE174409/NeuroinflammationResults/DEG_by_groups/DEG_",
                    sex_val, "_", region_val, ".csv")
  write_csv(results, filename)

  # Also create and save volcano plot
  create_volcano_plot(results, sex_val, region_val)

  return(results)
}

# Function to create volcano plot
create_volcano_plot <- function(deg_results, sex_val, region_val) {
  # Set significance thresholds
  sig_threshold <- 0.05
  fc_threshold <- 1  # log2 scale (2-fold change)

  # Standardize column names
  deg_results <- deg_results %>%
    rename(log2FoldChange = logFC, pvalue = P.Value, padj = adj.P.Val)

  # Label significant complement genes
  genes_to_label <- deg_results %>%
    filter(is_complement == TRUE & padj < sig_threshold & abs(log2FoldChange) > fc_threshold)

  # If no significant complement genes, label top complement genes by p-value
  if (nrow(genes_to_label) == 0) {
    genes_to_label <- deg_results %>%
      filter(is_complement == TRUE) %>%
      arrange(pvalue) %>%
      head(10)
  }

  # Count significant genes
  sig_count <- sum(deg_results$padj < sig_threshold, na.rm=TRUE)
  comp_sig_count <- sum(deg_results$is_complement & deg_results$padj < sig_threshold, na.rm=TRUE)

  # Create volcano plot
  volcano_plot <- ggplot(deg_results, aes(x = log2FoldChange, y = -log10(pvalue))) +
    # Add points with different colors for complement genes
    geom_point(aes(color = is_complement, alpha = padj < sig_threshold), size = 1.5) +
    # Add threshold lines
    geom_hline(yintercept = -log10(sig_threshold), linetype = "dashed", color = "darkgrey") +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "darkgrey") +
    # Add labels for complement genes
    geom_label_repel(
      data = genes_to_label,
      aes(label = gene_symbol),
      box.padding = 0.5,
      max.overlaps = 15,
      size = 3
    ) +
    # Set colors and transparency
    scale_color_manual(values = c("grey30", "red"), labels = c("Other genes", "Complement genes")) +
    scale_alpha_manual(values = c(0.4, 0.8)) +
    # Add labels
    labs(
      title = paste0("Volcano Plot: ", sex_val, " - ", region_val),
      subtitle = paste0("Significant DEGs: ", sig_count,
                      " (", comp_sig_count, " complement genes)"),
      x = "Log2 Fold Change",
      y = "-Log10 P-value",
      color = "Gene Type"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )

  # Save plot
  filename <- paste0("GSE174409/figures/volcano_plots/volcano_", sex_val, "_", region_val)
  ggsave(paste0(filename, ".pdf"), volcano_plot, width = 8, height = 7)
  ggsave(paste0(filename, ".png"), volcano_plot, width = 8, height = 7, dpi = 300)

  return(volcano_plot)
}

# Run analysis for each combination
all_results <- list()
all_plots <- list()

for (sex_val in sexes) {
  for (region_val in regions) {
    results <- run_deg_analysis(sex_val, region_val)
    if (!is.null(results)) {
      key <- paste(sex_val, region_val, sep="_")
      all_results[[key]] <- results
    }
  }
}

# ===== PART 3: Create summary figure of all plots =====
cat("\n==== 3. Creating summary figure ====\n")

# Recreate all plots with consistent limits for comparison
all_plots <- list()
max_y <- 0
max_x <- 0

# First pass to determine common axis limits
for (sex_val in sexes) {
  for (region_val in regions) {
    key <- paste(sex_val, region_val, sep="_")
    if (key %in% names(all_results)) {
      results <- all_results[[key]]
      max_y <- max(max_y, -log10(min(results$P.Value[results$P.Value > 0], na.rm=TRUE)), na.rm=TRUE)
      max_x <- max(max_x, max(abs(results$logFC), na.rm=TRUE))
    }
  }
}

# Add a small buffer to the limits
max_x <- max_x * 1.05
max_y <- max_y * 1.05

# Second pass to create plots with common limits
for (sex_val in sexes) {
  for (region_val in regions) {
    key <- paste(sex_val, region_val, sep="_")
    if (key %in% names(all_results)) {
      # Standardize column names
      deg_results <- all_results[[key]] %>%
        rename(log2FoldChange = logFC, pvalue = P.Value, padj = adj.P.Val)

      # Set significance thresholds
      sig_threshold <- 0.05
      fc_threshold <- 1

      # Label significant complement genes
      genes_to_label <- deg_results %>%
        filter(is_complement == TRUE & padj < sig_threshold & abs(log2FoldChange) > fc_threshold)

      # If no significant complement genes, label top complement genes by p-value
      if (nrow(genes_to_label) == 0) {
        genes_to_label <- deg_results %>%
          filter(is_complement == TRUE) %>%
          arrange(pvalue) %>%
          head(10)
      }

      # Count significant genes
      sig_count <- sum(deg_results$padj < sig_threshold, na.rm=TRUE)
      comp_sig_count <- sum(deg_results$is_complement & deg_results$padj < sig_threshold, na.rm=TRUE)

      # Create volcano plot with consistent limits
      p <- ggplot(deg_results, aes(x = log2FoldChange, y = -log10(pvalue))) +
        geom_point(aes(color = is_complement, alpha = padj < sig_threshold), size = 1.2) +
        geom_hline(yintercept = -log10(sig_threshold), linetype = "dashed", color = "darkgrey") +
        geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "darkgrey") +
        geom_label_repel(
          data = genes_to_label,
          aes(label = gene_symbol),
          box.padding = 0.5,
          max.overlaps = 10,
          size = 2.5
        ) +
        scale_color_manual(values = c("grey30", "red"), labels = c("Other", "Complement")) +
        scale_alpha_manual(values = c(0.4, 0.8)) +
        xlim(-max_x, max_x) +
        ylim(0, max_y) +
        labs(
          title = paste0(sex_val, " - ", region_val),
          subtitle = paste0("DEGs: ", sig_count, " (", comp_sig_count, " complement)"),
          x = "Log2 Fold Change",
          y = "-Log10 P-value"
        ) +
        theme_bw() +
        theme(
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 9),
          axis.title = element_text(size = 10)
        )

      all_plots[[key]] <- p
    }
  }
}

# Arrange plots in a grid - FIX FOR THE ERROR
if (length(all_plots) > 0) {
  # Print the number of plots for debugging
  cat("Number of plots to arrange:", length(all_plots), "\n")

  # Create a common legend
  dummy_plot <- ggplot(data.frame(x = 1, y = 1, is_complement = c(FALSE, TRUE), sig = c(FALSE, TRUE))) +
    geom_point(aes(x = x, y = y, color = is_complement, alpha = sig)) +
    scale_color_manual(values = c("grey30", "red"), labels = c("Other genes", "Complement genes")) +
    scale_alpha_manual(values = c(0.4, 0.8), labels = c("Not significant", "Significant (FDR < 0.05)")) +
    theme(legend.position = "bottom")

  legend <- cowplot::get_legend(dummy_plot)

  # Calculate grid dimensions - FIXED calculation
  n_plots <- length(all_plots)
  n_cols <- min(2, n_plots)  # Use max 2 columns instead of 3
  n_rows <- ceiling(n_plots / n_cols)

  # Ensure the product of rows and columns is at least equal to the number of plots
  if (n_rows * n_cols < n_plots) {
    n_cols <- ceiling(sqrt(n_plots))
    n_rows <- ceiling(n_plots / n_cols)
  }

  cat("Grid dimensions:", n_rows, "rows by", n_cols, "columns\n")

  # Convert all_plots list to a list with names
  plot_list <- all_plots

  # Add error handling for grid arrangement
  tryCatch({
    # Create the grid
    grid_plot <- gridExtra::grid.arrange(
      grobs = plot_list,
      ncol = n_cols,
      nrow = n_rows,
      bottom = legend
    )

    # Save the summary figure
    ggsave("GSE174409/figures/volcano_plots/all_volcano_plots_summary.pdf",
           grid_plot, width = min(12, 6*n_cols), height = min(15, 5*n_rows + 1))
    ggsave("GSE174409/figures/volcano_plots/all_volcano_plots_summary.png",
           grid_plot, width = min(12, 6*n_cols), height = min(15, 5*n_rows + 1), dpi = 300)

    cat("Created summary figure with", n_plots, "volcano plots\n")
  }, error = function(e) {
    cat("Error creating grid:", as.character(e), "\n")
    cat("Saving individual plots instead\n")

    # Alternative: save a separate comparison figure
    pdf("GSE174409/figures/volcano_plots/all_volcano_plots_separate.pdf", width = 12, height = 10)
    for (p in plot_list) {
      print(p)
    }
    dev.off()
  })
}

cat("\nAnalysis complete. Results saved in 'GSE174409/NeuroinflammationResults/DEG_by_groups/'\n")
cat("Volcano plots saved in 'GSE174409/figures/volcano_plots/'\n")