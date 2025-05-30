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
library(cowplot)

# ===== PART 1: Load Data =====
message("==== 1. Loading data ====")
expr <- readRDS('/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE174409/QC/logcpm_filtered_normalized.rds')
meta <- readRDS('/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE174409/QC/metadata.rds')

# Strip Ensembl version suffixes if present
rownames(expr) <- gsub("\\.\\d+$", "", rownames(expr))

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
message("Matched samples: ", length(common_samples))
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
message("==== 2. Running DEG analysis for each sex and region ====")

# Create output directories
dir.create("GSE174409/NeuroinflammationResults/DEG_by_groups", showWarnings=FALSE, recursive=TRUE)
dir.create("GSE174409/figures/volcano_plots", showWarnings=FALSE, recursive=TRUE)
dir.create("GSE174409/NeuroinflammationResults/complement_genes", showWarnings=FALSE, recursive=TRUE)

# Get unique sex and region values
sexes <- unique(meta_matched$sex)
regions <- unique(meta_matched$region)

# Define list of key genes you want to always highlight
genes_of_interest <- c("C1QA", "C1QB", "C3", "C3AR1", "C4A", "C5", "C5AR1", "CFH")

# Function to create volcano plot with flexible parameters and additional gene labeling
create_volcano_plot <- function(deg_results, sex_val, region_val, xlim = NULL, ylim = NULL,
                              save_plot = TRUE, detailed_colors = FALSE, additional_genes = NULL,
                              label_top_deg = 15) {
  # Set significance thresholds - more lenient for visualization
  sig_threshold <- 0.05  # FDR threshold
  fc_threshold <- 0.5    # Lower for visualization (log2FC of 0.5 = ~1.4 fold)

  # Standardize column names
  deg_results <- deg_results %>%
    rename(log2FoldChange = logFC, pvalue = P.Value, padj = adj.P.Val)

  # Filter out NA values which would break the plot
  deg_results <- deg_results %>%
    filter(!is.na(pvalue) & !is.na(log2FoldChange) & !is.na(padj))

  # Count significant genes
  sig_count <- sum(deg_results$padj < sig_threshold, na.rm=TRUE)
  comp_sig_count <- sum(deg_results$is_complement & deg_results$padj < sig_threshold, na.rm=TRUE)

  # If no data to plot, create an empty plot with a message
  if (nrow(deg_results) == 0) {
    p <- ggplot() +
      annotate("text", x = 0, y = 0, label = "No data to plot") +
      theme_minimal() +
      labs(title = paste0("Volcano Plot: ", sex_val, " - ", region_val),
           subtitle = "No genes passed filtering criteria")
    return(p)
  }

  # Create color grouping for more detailed visualization
  if (detailed_colors) {
    deg_results$color_group <- case_when(
      deg_results$is_complement & deg_results$log2FoldChange > 0 & deg_results$padj < sig_threshold ~ "Complement Up",
      deg_results$is_complement & deg_results$log2FoldChange < 0 & deg_results$padj < sig_threshold ~ "Complement Down",
      deg_results$is_complement ~ "Complement (NS)",
      deg_results$padj < sig_threshold & deg_results$log2FoldChange > 0 ~ "Other DEG Up",
      deg_results$padj < sig_threshold & deg_results$log2FoldChange < 0 ~ "Other DEG Down",
      TRUE ~ "Non-significant"
    )

    # Initialize genes to label
    genes_to_label <- deg_results[0,]

    # Add complement genes that are significant
    complement_sig <- deg_results %>%
      filter(is_complement == TRUE & padj < sig_threshold)
    if (nrow(complement_sig) > 0) {
      genes_to_label <- bind_rows(genes_to_label, complement_sig)
    }

    # Add top complement genes even if not significant
    top_complement <- deg_results %>%
      filter(is_complement == TRUE) %>%
      arrange(pvalue) %>%
      head(10)
    if (nrow(top_complement) > 0) {
      genes_to_label <- bind_rows(genes_to_label, top_complement)
    }

    # Add top DEGs by p-value
    top_degs <- deg_results %>%
      filter(padj < sig_threshold | (is.na(padj) & pvalue < 0.01)) %>%
      arrange(pvalue) %>%
      head(label_top_deg)
    if (nrow(top_degs) > 0) {
      genes_to_label <- bind_rows(genes_to_label, top_degs)
    }

    # If still no genes to label, get top genes by fold change
    if (nrow(genes_to_label) == 0) {
      top_fc <- deg_results %>%
        arrange(desc(abs(log2FoldChange))) %>%
        head(label_top_deg)
      genes_to_label <- bind_rows(genes_to_label, top_fc)
    }

    # Add user-specified genes
    if (!is.null(additional_genes)) {
      additional_to_label <- deg_results %>%
        filter(gene_symbol %in% additional_genes | gene_id %in% additional_genes)
      if (nrow(additional_to_label) > 0) {
        genes_to_label <- bind_rows(genes_to_label, additional_to_label)
      }
    }

    # Remove duplicates
    genes_to_label <- distinct(genes_to_label)

    # Create volcano plot with detailed colors
    p <- ggplot(deg_results, aes(x = log2FoldChange, y = -log10(pvalue))) +
      geom_point(aes(color = color_group), size = 1.5, alpha = 0.8) +
      geom_hline(yintercept = -log10(sig_threshold), linetype = "dashed", color = "darkgrey") +
      geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "darkgrey") +
      geom_label_repel(
        data = genes_to_label,
        aes(label = gene_symbol),
        box.padding = 0.5,
        max.overlaps = 30,
        size = 3
      ) +
      scale_color_manual(values = c(
        "Complement Up" = "#D7263D",
        "Complement Down" = "#1B998B",
        "Complement (NS)" = "#FFA07A",
        "Other DEG Up" = "#2E86C1",
        "Other DEG Down" = "#7D3C98",
        "Non-significant" = "grey80"
      ), name = "Gene Type") +
      labs(
        title = paste0("Volcano Plot: ", sex_val, " - ", region_val),
        subtitle = paste0("Significant DEGs: ", sig_count,
                          " (", comp_sig_count, " complement genes)"),
        x = "Log2 Fold Change",
        y = "-Log10 P-value"
      )
  } else {
    # Initialize genes to label
    genes_to_label <- deg_results[0,]

    # Add complement genes and top DEGs
    complement_genes <- deg_results %>%
      filter(is_complement == TRUE) %>%
      arrange(pvalue) %>%
      head(10)

    top_degs <- deg_results %>%
      filter(padj < sig_threshold) %>%
      arrange(pvalue) %>%
      head(label_top_deg)

    genes_to_label <- bind_rows(genes_to_label, complement_genes, top_degs)

    # Add user-specified genes
    if (!is.null(additional_genes)) {
      additional_to_label <- deg_results %>%
        filter(gene_symbol %in% additional_genes | gene_id %in% additional_genes)
      genes_to_label <- bind_rows(genes_to_label, additional_to_label)
    }

    # Remove duplicates
    genes_to_label <- distinct(genes_to_label)

    # Original color scheme
    p <- ggplot(deg_results, aes(x = log2FoldChange, y = -log10(pvalue))) +
      geom_point(aes(color = is_complement, alpha = padj < sig_threshold), size = 1.5) +
      geom_hline(yintercept = -log10(sig_threshold), linetype = "dashed", color = "darkgrey") +
      geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "darkgrey") +
      geom_label_repel(
        data = genes_to_label,
        aes(label = gene_symbol),
        box.padding = 0.5,
        max.overlaps = 30,
        size = 3
      ) +
      scale_color_manual(values = c("grey30", "red"), labels = c("Other genes", "Complement genes")) +
      scale_alpha_manual(values = c(0.4, 0.8)) +
      labs(
        title = paste0("Volcano Plot: ", sex_val, " - ", region_val),
        subtitle = paste0("Significant DEGs: ", sig_count,
                          " (", comp_sig_count, " complement genes)"),
        x = "Log2 Fold Change",
        y = "-Log10 P-value",
        color = "Gene Type"
      )
  }

  # Apply axis limits if provided
  if (!is.null(xlim)) p <- p + xlim(xlim[1], xlim[2])
  if (!is.null(ylim)) p <- p + ylim(ylim[1], ylim[2])

  # Apply theme
  p <- p + theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )

  # Save plot if requested
  if (save_plot) {
    filename <- paste0("GSE174409/figures/volcano_plots/volcano_", sex_val, "_", region_val)

    # Check if files already exist
    if (file.exists(paste0(filename, ".pdf")) || file.exists(paste0(filename, ".png"))) {
      message("Overwriting existing volcano plot files for ", sex_val, " - ", region_val)
    }

    ggsave(paste0(filename, ".pdf"), p, width = 8, height = 7)
    ggsave(paste0(filename, ".png"), p, width = 8, height = 7, dpi = 300)
  }

  return(p)
}

# Create a function to run DEG analysis for a subset
run_deg_analysis <- function(sex_val, region_val, additional_genes = NULL,
                             use_treat = TRUE, lfc_threshold = 0.5) {
  # Filter data for this combination
  subset_idx <- meta_matched$sex == sex_val & meta_matched$region == region_val
  if (sum(subset_idx) < 6) {
    message("Skipping ", sex_val, " ", region_val, " - too few samples (", sum(subset_idx), ")")
    return(NULL)
  }

  subset_meta <- meta_matched[subset_idx, ]
  subset_expr <- expr_matched[, subset_idx]

  # Check if we have enough samples per group
  subset_meta$diagnosis <- factor(subset_meta$diagnosis)
  group_counts <- table(subset_meta$diagnosis)
  if (min(group_counts) < 2) {
    message("Skipping ", sex_val, " ", region_val, " - insufficient replicates in at least one group")
    return(NULL)
  }

  message("Running DEG analysis for ", sex_val, " in ", region_val,
          " (n = ", sum(subset_idx), ", groups: ", paste(names(group_counts), collapse="/"), ")")

  # Create design matrix and run limma
  # Check if we have enough samples to include covariates
  design <- model.matrix(~ diagnosis, data=subset_meta)

  # Fit the linear model
  fit <- lmFit(subset_expr, design)

  # Apply empirical Bayes and fold change threshold (treat)
  if (use_treat) {
    message("  Using 'treat' method with log2FC threshold ", lfc_threshold)
    fit <- treat(fit, lfc=lfc_threshold)
    # Get results using topTreat to enforce fold-change threshold
    results <- topTreat(fit, coef=2, number=Inf, sort.by="p") %>%
      rownames_to_column("gene_id")
  } else {
    # Standard eBayes
    message("  Using standard 'eBayes' method")
    fit <- eBayes(fit)
    # Get results
    results <- topTable(fit, coef=2, number=Inf, sort.by="P") %>%
      rownames_to_column("gene_id")
  }

  # Add gene symbols and flag complement genes
  results$gene_symbol <- ensembl_to_symbol[results$gene_id]
  results$is_complement <- results$gene_id %in% ensembl_ids

  # Replace NA symbols with gene_id
  results$gene_symbol[is.na(results$gene_symbol)] <- results$gene_id[is.na(results$gene_symbol)]

  # Report summary statistics
  sig_genes <- sum(results$adj.P.Val < 0.05, na.rm = TRUE)
  sig_up <- sum(results$adj.P.Val < 0.05 & results$logFC > 0, na.rm = TRUE)
  sig_down <- sum(results$adj.P.Val < 0.05 & results$logFC < 0, na.rm = TRUE)
  message("  Found ", sig_genes, " significant DEGs (", sig_up, " up, ", sig_down, " down)")

  # If using treat but finding no DEGs, run again with standard approach
  if (use_treat && sig_genes == 0) {
    message("  No significant DEGs found with treat method, running standard analysis")
    return(run_deg_analysis(sex_val, region_val, additional_genes, use_treat = FALSE))
  }

  # Save results
  method_prefix <- ifelse(use_treat, "treatDEG_", "DEG_")
  filename <- paste0("GSE174409/NeuroinflammationResults/DEG_by_groups/", method_prefix,
                     sex_val, "_", region_val, ".csv")
  write_csv(results, filename)

  # Export complement genes specifically to separate CSV
  complement_results <- results %>%
    filter(is_complement == TRUE) %>%
    arrange(P.Value)

  # Save complement genes results
  comp_filename <- paste0("GSE174409/NeuroinflammationResults/complement_genes/complement_",
                          method_prefix, sex_val, "_", region_val, ".csv")
  write_csv(complement_results, comp_filename)

  # Create and save volcano plot with detailed colors and additional genes
  create_volcano_plot(results, sex_val, region_val, detailed_colors = TRUE,
                      additional_genes = additional_genes)

  return(results)
}

# Run analysis for each combination
all_results <- list()
all_plots <- list()

for (sex_val in sexes) {
  for (region_val in regions) {
    # Try with a lower fold change threshold (0.5 = ~1.4 fold)
    results <- run_deg_analysis(sex_val, region_val,
                               additional_genes = genes_of_interest,
                               use_treat = TRUE,
                               lfc_threshold = 0.5)

    if (!is.null(results)) {
      key <- paste(sex_val, region_val, sep="_")
      all_results[[key]] <- results
    }
  }
}

# ===== PART 3: Create summary figure of all plots =====
message("\n==== 3. Creating summary figure ====")

# Only create summary if we have results
if (length(all_results) > 0) {
  # First pass to determine common axis limits
  max_y <- 0
  max_x <- 0

  for (sex_val in sexes) {
    for (region_val in regions) {
      key <- paste(sex_val, region_val, sep="_")
      if (key %in% names(all_results)) {
        results <- all_results[[key]]
        # Filter out NA values
        valid_results <- results %>% filter(!is.na(P.Value) & P.Value > 0 & !is.na(logFC))
        if (nrow(valid_results) > 0) {
          max_y <- max(max_y, -log10(min(valid_results$P.Value, na.rm=TRUE)), na.rm=TRUE)
          max_x <- max(max_x, max(abs(valid_results$logFC), na.rm=TRUE))
        }
      }
    }
  }

  # Add a small buffer to the limits
  max_x <- max_x * 1.05
  max_y <- max_y * 1.05

  # Second pass to create plots with common limits
  plot_list <- list()
  for (sex_val in sexes) {
    for (region_val in regions) {
      key <- paste(sex_val, region_val, sep="_")
      if (key %in% names(all_results)) {
        # Create simplified volcano plot for summary with consistent limits
        p <- create_volcano_plot(
          all_results[[key]],
          sex_val,
          region_val,
          xlim = c(-max_x, max_x),
          ylim = c(0, max_y),
          save_plot = FALSE,
          additional_genes = genes_of_interest
        ) +
          theme(
            legend.position = "none",
            plot.title = element_text(hjust = 0.5, size = 11),
            plot.subtitle = element_text(hjust = 0.5, size = 9),
            axis.title = element_text(size = 10)
          )

        plot_list[[key]] <- p
      }
    }
  }

  # Arrange plots in a grid
  if (length(plot_list) > 0) {
    # Print the number of plots for debugging
    message("Number of plots to arrange: ", length(plot_list))

    # Create a common legend
    dummy_plot <- ggplot(data.frame(x = 1, y = 1, is_complement = c(FALSE, TRUE), sig = c(FALSE, TRUE))) +
      geom_point(aes(x = x, y = y, color = is_complement, alpha = sig)) +
      scale_color_manual(values = c("grey30", "red"), labels = c("Other genes", "Complement genes")) +
      scale_alpha_manual(values = c(0.4, 0.8), labels = c("Not significant", "Significant (FDR < 0.05)")) +
      theme(legend.position = "bottom")

    legend <- cowplot::get_legend(dummy_plot)

    # Calculate grid dimensions
    n_plots <- length(plot_list)
    n_cols <- min(2, n_plots)
    n_rows <- ceiling(n_plots / n_cols)

    # Ensure the product of rows and columns is at least equal to the number of plots
    if (n_rows * n_cols < n_plots) {
      n_cols <- ceiling(sqrt(n_plots))
      n_rows <- ceiling(n_plots / n_cols)
    }

    message("Grid dimensions: ", n_rows, " rows by ", n_cols, " columns")

    # Add error handling for grid arrangement
    tryCatch({
      # Use cowplot instead of gridExtra for better layout control
      main_grid <- cowplot::plot_grid(
        plotlist = plot_list,
        ncol = n_cols,
        nrow = n_rows,
        align = "hv"
      )

      # Add the legend at the bottom
      grid_plot <- cowplot::plot_grid(
        main_grid,
        legend,
        ncol = 1,
        rel_heights = c(1, 0.1)
      )

      # Save the summary figure
      ggsave("GSE174409/figures/volcano_plots/all_volcano_plots_summary.pdf",
             grid_plot, width = min(12, 6*n_cols), height = min(15, 5*n_rows + 1))
      ggsave("GSE174409/figures/volcano_plots/all_volcano_plots_summary.png",
             grid_plot, width = min(12, 6*n_cols), height = min(15, 5*n_rows + 1), dpi = 300)

      message("Created summary figure with ", n_plots, " volcano plots")
    }, error = function(e) {
      message("Error creating grid: ", as.character(e))
      message("Saving individual plots instead")

      # Alternative: save a separate comparison figure
      pdf("GSE174409/figures/volcano_plots/all_volcano_plots_separate.pdf", width = 12, height = 10)
      for (p in plot_list) {
        print(p)
      }
      dev.off()
    })
  }
}

# Create a combined CSV of all complement genes from all analyses
message("\n==== 4. Creating summary complement genes file ====")
if (length(all_results) > 0) {
  all_complement_genes <- bind_rows(lapply(names(all_results), function(key) {
    parts <- strsplit(key, "_")[[1]]
    sex_val <- parts[1]
    region_val <- parts[2]

    results <- all_results[[key]]
    complement_results <- results %>%
      filter(is_complement == TRUE) %>%
      mutate(
        sex = sex_val,
        region = region_val,
        significant = adj.P.Val < 0.05,
        direction = ifelse(logFC > 0, "Up", "Down")
      ) %>%
      select(gene_id, gene_symbol, logFC, P.Value, adj.P.Val, sex, region, significant, direction)

    return(complement_results)
  }))

  # Save the combined complement genes table
  write_csv(all_complement_genes, "GSE174409/NeuroinflammationResults/complement_genes/all_complement_genes_summary.csv")
  message("Created summary file of all complement genes across conditions")
} else {
  message("No results to summarize")
}

message("\nAnalysis complete. Results saved in 'GSE174409/NeuroinflammationResults/DEG_by_groups/'")
message("Volcano plots saved in 'GSE174409/figures/volcano_plots/'")
message("Complement gene results saved in 'GSE174409/NeuroinflammationResults/complement_genes/'")