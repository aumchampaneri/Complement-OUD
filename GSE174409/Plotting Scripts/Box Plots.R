# ===== SETUP =====
# Install required packages if not already installed
required_packages <- c("ggplot2", "ggbeeswarm", "ggpubr", "dplyr", "rstatix", "patchwork")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# ===== DATA LOADING =====
expression_data_path <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE174409/QC/logcpm_filtered_normalized.rds"
metadata_path <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE174409/QC/metadata.rds"
gene_mapping_path <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE174409/NeuroinflammationResults/ensembl_to_symbol_mapping.csv"

# Load data files
if (file.exists(expression_data_path)) {
  logcpm_filtered_norm <- readRDS(expression_data_path)
} else {
  stop("Expression data file not found at: ", expression_data_path)
}

if (file.exists(metadata_path)) {
  metadata_df <- readRDS(metadata_path)
} else {
  stop("Metadata file not found at: ", metadata_path)
}

if (file.exists(gene_mapping_path)) {
  gene_mapping <- read.csv(gene_mapping_path)
} else {
  stop("Gene mapping file not found at: ", gene_mapping_path)
}

# Define gene list to plot
complement_genes <- c("C1QA", "C1QB", "C3", "C3AR1", "C4A", "CFH", "C5", "C5AR1")

# Set output directory
complement_plot_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE174409/plots/complement_genes"

# ===== HELPER FUNCTIONS =====
find_gene_id <- function(gene_symbol, gene_mapping) {
  row_idx <- which(gene_mapping$Symbol == gene_symbol)
  if (length(row_idx) > 0) {
    return(gene_mapping$EnsemblID[row_idx[1]])
  } else {
    row_idx <- which(toupper(gene_mapping$Symbol) == toupper(gene_symbol))
    if (length(row_idx) > 0) {
      return(gene_mapping$EnsemblID[row_idx[1]])
    } else {
      return(NULL)
    }
  }
}

# ===== PLOTTING FUNCTIONS =====
# For plots comparing Sex and Diagnosis by Region
plot_gene_with_significance <- function(gene_symbol,
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

  # Create combined group variable
  df$Group <- interaction(df$Sex, df$Diagnosis)

  # Filter by region if needed
  if (!facet_by_region && !is.null(region)) {
    df <- df[df$Region == region, ]
    if (nrow(df) == 0) {
      cat("No data for gene", gene_symbol, "in region", region, "\n")
      return(NULL)
    }
  }

  if (facet_by_region) {
    # For faceted plots, create a separate plot for each region with significance testing
    regions <- unique(df$Region)
    plot_list <- list()

    for (reg in regions) {
      reg_df <- df[df$Region == reg, ]

      p <- ggplot(reg_df, aes(x = Group, y = Expression, fill = Group)) +
        labs(subtitle = reg) +
        theme(plot.subtitle = element_text(hjust = 0.5, face = "bold"))

      if (plot_type == "violin") {
        p <- p + geom_violin(alpha = 0.7, trim = FALSE)
        if (add_points) {
          p <- p + geom_quasirandom(alpha = 0.6, size = 1)
        }
      } else {
        p <- p + geom_boxplot(alpha = 0.7, outlier.shape = NA)
        if (add_points) {
          p <- p + geom_quasirandom(alpha = 0.6, size = 1)
        }
      }

      # Add significance for this region
      if (length(unique(reg_df$Group)) >= 2) {
        # Create all pairwise comparisons
        comparisons <- list()
        groups <- as.character(unique(reg_df$Group))
        for (i in 1:(length(groups)-1)) {
          for (j in (i+1):length(groups)) {
            comparisons[[length(comparisons) + 1]] <- c(groups[i], groups[j])
          }
        }

        p <- p + stat_compare_means(
          comparisons = comparisons,
          method = "wilcox.test",
          label = "p.signif",
          p.adjust.method = "BH"
        )
      }

      # Format region plot
      p <- p +
        labs(y = NULL, x = NULL) +
        scale_fill_brewer(palette = "Set1") +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none"
        )

      plot_list[[reg]] <- p
    }

    # Combine all region plots
    combined_plot <- patchwork::wrap_plots(plot_list, ncol = 2) +
      plot_annotation(
        title = paste("Expression of", gene_symbol, "by Sex and Diagnosis"),
        theme = theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
        )
      ) +
      patchwork::plot_layout(guides = "collect") &
      theme(legend.position = "none") &
      ylab("Log2 CPM Expression")

    return(combined_plot)

  } else {
    # For non-faceted plots (single region)
    p <- ggplot(df, aes(x = Group, y = Expression, fill = Group))

    if (plot_type == "violin") {
      p <- p + geom_violin(alpha = 0.7, trim = FALSE)
      if (add_points) {
        p <- p + geom_quasirandom(alpha = 0.6, size = 1)
      }
    } else {
      p <- p + geom_boxplot(alpha = 0.7, outlier.shape = NA)
      if (add_points) {
        p <- p + geom_quasirandom(alpha = 0.6, size = 1)
      }
    }

    # Add significance testing for all pairs
    if (length(unique(df$Group)) >= 2) {
      # Create all pairwise comparisons
      comparisons <- list()
      groups <- as.character(unique(df$Group))
      for (i in 1:(length(groups)-1)) {
        for (j in (i+1):length(groups)) {
          comparisons[[length(comparisons) + 1]] <- c(groups[i], groups[j])
        }
      }

      p <- p + stat_compare_means(
        comparisons = comparisons,
        method = "wilcox.test",
        label = "p.signif",
        p.adjust.method = "BH"
      )
    }

    # Format single region plot
    reg_text <- if(!is.null(region)) paste(" in", region) else ""
    p <- p +
      labs(
        title = paste("Expression of", gene_symbol, "by Sex and Diagnosis", reg_text),
        y = "Log2 CPM Expression",
        x = ""
      ) +
      scale_fill_brewer(palette = "Set1") +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(size = 14),
        legend.position = "none"
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.30)))

    return(p)
  }
}

# For plots comparing regions with faceting by Sex and Diagnosis
plot_gene_by_region <- function(gene_symbol,
                               plot_type = "boxplot",
                               add_points = TRUE) {

  ensembl_id <- find_gene_id(gene_symbol, gene_mapping)
  if (is.null(ensembl_id)) {
    cat("Could not find Ensembl ID for", gene_symbol, "\n")
    return(NULL)
  }

  df <- data.frame(
    Expression = logcpm_filtered_norm[ensembl_id, ],
    Sex = metadata_df$sex,
    Diagnosis = metadata_df$diagnosis,
    Region = metadata_df$region
  )

  # Create combined group for faceting
  df$SexDiagnosis <- paste(df$Sex, df$Diagnosis, sep = "_")

  # Create separate plots for each combination of Sex and Diagnosis
  plot_list <- list()
  facet_groups <- unique(df$SexDiagnosis)

  for (group in facet_groups) {
    group_df <- df[df$SexDiagnosis == group, ]
    sex_diag <- strsplit(group, "_")[[1]]

    p <- ggplot(group_df, aes(x = Region, y = Expression, fill = Region)) +
      labs(subtitle = paste(sex_diag[1], sex_diag[2])) +
      theme(plot.subtitle = element_text(hjust = 0.5, face = "bold"))

    if (plot_type == "violin") {
      p <- p + geom_violin(alpha = 0.7, trim = FALSE)
      if (add_points) {
        p <- p + geom_quasirandom(alpha = 0.6, size = 1)
      }
    } else {
      p <- p + geom_boxplot(alpha = 0.7, outlier.shape = NA)
      if (add_points) {
        p <- p + geom_quasirandom(alpha = 0.6, size = 1)
      }
    }

    # Add significance testing for regions
    if (length(unique(group_df$Region)) >= 2) {
      # Create all pairwise comparisons between regions
      comparisons <- list()
      regions <- unique(group_df$Region)
      for (i in 1:(length(regions)-1)) {
        for (j in (i+1):length(regions)) {
          comparisons[[length(comparisons) + 1]] <- c(regions[i], regions[j])
        }
      }

      # Add significance testing
      p <- p + stat_compare_means(
        comparisons = comparisons,
        method = "wilcox.test",
        label = "p.signif",
        p.adjust.method = "BH"
      )
    }

    # Format region plot
    p <- p +
      labs(y = NULL, x = NULL) +
      scale_fill_brewer(palette = "Set2") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )

    plot_list[[group]] <- p
  }

  # Combine all Sex/Diagnosis plots
  combined_plot <- patchwork::wrap_plots(plot_list, ncol = 2) +
    plot_annotation(
      title = paste("Expression of", gene_symbol, "by Region"),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
      )
    ) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "none") &
    ylab("Log2 CPM Expression")

  return(combined_plot)
}

# ===== MAIN FUNCTION =====
generate_all_plots <- function(gene_list, output_dir) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Loop through all genes and create plots
  for (gene in gene_list) {
    cat("Processing gene:", gene, "\n")

    # 1. Create plot comparing Sex and Diagnosis with faceting by Region
    p1 <- plot_gene_with_significance(
      gene_symbol = gene,
      plot_type = "boxplot",
      facet_by_region = TRUE,
      add_points = TRUE
    )

    if (!is.null(p1)) {
      filename <- file.path(output_dir, paste0(gene, "_by_sex_diagnosis.png"))
      ggsave(filename, p1, width = 12, height = 10, dpi = 300)
      cat("  Saved:", filename, "\n")
    }

    # 2. Create plot comparing regions with faceting by Sex and Diagnosis
    p2 <- plot_gene_by_region(
      gene_symbol = gene,
      plot_type = "boxplot",
      add_points = TRUE
    )

    if (!is.null(p2)) {
      filename <- file.path(output_dir, paste0(gene, "_by_region.png"))
      ggsave(filename, p2, width = 12, height = 10, dpi = 300)
      cat("  Saved:", filename, "\n")
    }

    # 3. For each region, create a non-faceted plot that shows all pairwise comparisons
    regions <- unique(metadata_df$region)
    for (reg in regions) {
      p3 <- plot_gene_with_significance(
        gene_symbol = gene,
        plot_type = "boxplot",
        facet_by_region = FALSE,
        region = reg,
        add_points = TRUE
      )

      if (!is.null(p3)) {
        filename <- file.path(output_dir, paste0(gene, "_", reg, ".png"))
        ggsave(filename, p3, width = 8, height = 6, dpi = 300)
        cat("  Saved:", filename, "\n")
      }
    }
  }

  cat("All plots generated successfully!\n")
}

# ===== RUN THE SCRIPT =====
generate_all_plots(complement_genes, complement_plot_dir)