#!/usr/bin/env Rscript

# =============================================================================
# Transcription Factor Activity Inference using DecoupleR + DoRothEA
# =============================================================================
# 
# Description: Comprehensive TF activity analysis of differential expression 
#              results from edgeR analysis using DoRothEA regulons and decoupleR
# 
# Author: Multi-Omics Study Team
# Date: 2024
# 
# Input: edgeR differential expression results
# Output: TF activity scores, statistical tests, and visualization
# 
# Methods: DoRothEA regulons + decoupleR (VIPER, ULM, MLM)
# =============================================================================

# Load required libraries with error handling
required_packages <- c(
  "decoupleR", "dplyr", "tidyr", "readr", "ggplot2", "ComplexHeatmap", 
  "circlize", "pheatmap", "viridis", "RColorBrewer", "reshape2",
  "broom", "purrr", "tibble", "stringr", "corrplot", "igraph",
  "networkD3", "plotly", "DT", "htmlwidgets", "org.Hs.eg.db",
  "ggrepel", "rlang"
)

# Suppress dplyr messages and conflicts
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(purrr)
  library(stringr)
  library(readr)
  library(rlang)
})

# Function to install and load packages
load_packages <- function(packages) {
  failed_packages <- c()
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat("Installing package:", pkg, "\n")
      tryCatch({
        if (pkg %in% c("decoupleR", "org.Hs.eg.db", "ComplexHeatmap")) {
          BiocManager::install(pkg, quiet = TRUE, update = FALSE)
        } else {
          install.packages(pkg, quiet = TRUE)
        }
        library(pkg, character.only = TRUE)
        cat("Successfully installed and loaded:", pkg, "\n")
      }, error = function(e) {
        cat("Failed to install package:", pkg, "-", e$message, "\n")
        failed_packages <<- c(failed_packages, pkg)
      })
    } else {
      cat("Package already loaded:", pkg, "\n")
    }
  }
  
  if (length(failed_packages) > 0) {
    warning("Failed to install packages: ", paste(failed_packages, collapse = ", "))
    cat("Continuing with available packages...\n")
  }
}

# Install BiocManager if not available
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  library(BiocManager)
}

# Load packages
cat("Loading required packages for TF activity analysis...\n")
load_packages(required_packages)

# Set up logging
log_file <- "tf_activity_analysis.log"
log_messages <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- paste(timestamp, "-", message)
  cat(log_entry, "\n")
  write(log_entry, file = log_file, append = TRUE)
}

log_messages("Starting TF Activity Analysis with DecoupleR + DoRothEA")

# =============================================================================
# CONFIGURATION
# =============================================================================

# Set seed for reproducibility
set.seed(42)

# Create output directories
output_dirs <- list(
  main = "../../results/snrna_scvi/tf_activity_analysis",
  plots = "../../results/snrna_scvi/tf_activity_analysis/plots",
  tables = "../../results/snrna_scvi/tf_activity_analysis/tables",
  networks = "../../results/snrna_scvi/tf_activity_analysis/networks"
)

for (dir in output_dirs) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
}

# Color schemes for OUD analysis
color_schemes <- list(
  sex = c("Male" = "#3498db", "Female" = "#e74c3c"),
  condition = c("Control" = "#95a5a6", "OUD" = "#e67e22"),
  brain_region = c("Caudate" = "#9b59b6", "Putamen" = "#1abc9c"),
  tf_activity = viridis::viridis(100)
)

# =============================================================================
# DATA LOADING AND PREPARATION
# =============================================================================

log_messages("Loading edgeR differential expression results...")

# Define expected edgeR results files (adjust paths as needed)
edger_results_dir <- "../../results/snrna_scvi/differential_expression_edgeR"

# Alternative paths to search for edgeR results
possible_paths <- c(
  "../../results/snrna_scvi/differential_expression_edgeR",
  "results/snrna_scvi/differential_expression_edgeR",
  "results/differential_expression",
  "results/edgeR",
  "outputs/differential_expression",
  "outputs/edgeR",
  "data/edgeR_results",
  "."
)

# Find edgeR results files
# Look for edgeR results files
edger_files <- c()
for (path in possible_paths) {
  if (dir.exists(path)) {
    files <- list.files(path, pattern = "(edgeR.*results|DE|differential).*\\.(csv|tsv|txt)$", 
                       full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
    # Filter out summary files and focus on individual contrast results
    files <- files[!str_detect(basename(files), "(summary|all_contrasts|significant)")]
    edger_files <- c(edger_files, files)
  }
}

if (length(edger_files) == 0) {
  cat("No edgeR results files found. Please specify the correct path.\n")
  cat("Searched in:\n")
  for (path in possible_paths) {
    cat("  -", path, "\n")
  }
  cat("\nPlease update the 'edger_results_dir' variable in the script.\n")
  stop("No edgeR results files found")
}

log_messages(paste("Found", length(edger_files), "potential edgeR results files"))

# Load all edgeR results
edger_results <- list()
contrast_names <- c()

for (file in edger_files[seq_len(min(10, length(edger_files)))]) {  # Limit to first 10 files
  # Extract contrast name from filename
  contrast_name <- basename(file) %>% 
    str_remove("(?i)(edger?_?|de_?|differential_?)") %>% 
    str_remove("\\.(csv|tsv|txt)$") %>%
    str_replace_all("[^A-Za-z0-9_]", "_")
  
  if (contrast_name == "" || contrast_name %in% contrast_names) {
    contrast_name <- paste0("contrast_", length(contrast_names) + 1)
  }
  
  contrast_names <- c(contrast_names, contrast_name)
  
  # Load results
  tryCatch({
    if (str_detect(file, "\\.(csv)$")) {
      result <- read_csv(file, show_col_types = FALSE)
    } else {
      result <- read_delim(file, delim = "\t", show_col_types = FALSE)
    }
    
    # Check for required columns (flexible naming)
    logfc_cols <- names(result)[str_detect(names(result), "(?i)(logfc|log2fc|log_fc|fold_change)")]
    pval_cols <- names(result)[str_detect(names(result), "(?i)(pvalue|p_value|pval|p\\.value)")]
    fdr_cols <- names(result)[str_detect(names(result), "(?i)(fdr|adj|padj|p\\.adj|qvalue)")]
    
    if (length(logfc_cols) == 0 || length(pval_cols) == 0) {
      log_messages(paste("Warning: Missing required columns in", file))
      next
    }
    
    # Standardize column names - check if they already exist
    if (!"logFC" %in% names(result) && length(logfc_cols) > 0) {
      result <- result %>%
        dplyr::rename(logFC = !!rlang::sym(logfc_cols[1]))
    }
    
    if (!"PValue" %in% names(result) && length(pval_cols) > 0) {
      result <- result %>%
        dplyr::rename(PValue = !!rlang::sym(pval_cols[1]))
    }
    
    if (!"FDR" %in% names(result)) {
      if (length(fdr_cols) > 0) {
        result <- result %>% 
          dplyr::rename(FDR = !!rlang::sym(fdr_cols[1]))
      } else {
        result$FDR <- p.adjust(result$PValue, method = "fdr")
      }
    }
    
    # Add gene symbols if not present
    if (!"gene_symbol" %in% colnames(result)) {
      # Try different possible gene identifier columns
      gene_cols <- names(result)[str_detect(names(result), "(?i)(gene|symbol|name)")]
      ensembl_cols <- names(result)[str_detect(names(result), "(?i)(ensembl|ensg)")]
      
      if (length(gene_cols) > 0) {
        result$gene_symbol <- result[[gene_cols[1]]]
      } else if (length(ensembl_cols) > 0) {
        # Convert ENSEMBL to gene symbols
        library(org.Hs.eg.db)
        gene_symbols <- mapIds(org.Hs.eg.db, 
                             keys = result[[ensembl_cols[1]]],
                             column = "SYMBOL", 
                             keytype = "ENSEMBL",
                             multiVals = "first")
        result$gene_symbol <- gene_symbols
      } else if ("gene" %in% names(result)) {
        result$gene_symbol <- result$gene
      } else if (!"X1" %in% names(result) && !is.null(rownames(result))) {
        result$gene_symbol <- rownames(result)
      } else {
        result$gene_symbol <- result[[1]]  # Use first column as gene names
      }
    }
    
    edger_results[[contrast_name]] <- result
    log_messages(paste("Loaded", contrast_name, "with", nrow(result), "genes"))
    
  }, error = function(e) {
    log_messages(paste("Error loading", file, ":", e$message))
  })
}

if (length(edger_results) == 0) {
  stop("No valid edgeR results files could be loaded")
}

log_messages(paste("Successfully loaded", length(edger_results), "contrasts"))

# =============================================================================
# DOROTHEA REGULON PREPARATION
# =============================================================================

log_messages("Preparing DoRothEA transcription factor regulons...")

# Get DoRothEA regulons
# Using confidence levels A, B, C (high to medium confidence)
dorothea_regulons <- get_dorothea(organism = "human", levels = c("A", "B", "C"))

log_messages(paste("Raw DoRothEA data loaded with", nrow(dorothea_regulons), "interactions"))
log_messages(paste("DoRothEA columns:", paste(colnames(dorothea_regulons), collapse = ", ")))

# Check and fix DoRothEA structure for decoupleR compatibility
# The DoRothEA data structure might use 'source' instead of 'tf'
if ("source" %in% colnames(dorothea_regulons) && !"tf" %in% colnames(dorothea_regulons)) {
  dorothea_regulons <- dorothea_regulons %>%
    dplyr::rename(tf = source)
}

if ("weight" %in% colnames(dorothea_regulons) && !"mor" %in% colnames(dorothea_regulons)) {
  dorothea_regulons <- dorothea_regulons %>%
    dplyr::rename(mor = weight)
}

if (!"likelihood" %in% colnames(dorothea_regulons)) {
  dorothea_regulons <- dorothea_regulons %>%
    dplyr::mutate(likelihood = 1)
}

# Verify we have the correct columns now
required_cols <- c("tf", "target", "mor")
missing_cols <- required_cols[!required_cols %in% colnames(dorothea_regulons)]

if (length(missing_cols) > 0) {
  log_messages(paste("ERROR: Missing required columns:", paste(missing_cols, collapse = ", ")))
  log_messages(paste("Available columns:", paste(colnames(dorothea_regulons), collapse = ", ")))
  stop("DoRothEA data structure incompatible")
}

log_messages(paste("Loaded DoRothEA with", 
                  length(unique(dorothea_regulons$tf)), "TFs and",
                  nrow(dorothea_regulons), "interactions"))

# Show top TFs by number of targets
tf_summary <- dorothea_regulons %>%
  dplyr::group_by(.data$tf) %>%
  dplyr::summarise(
    n_targets = dplyr::n(),
    confidence_levels = paste(unique(.data$confidence), collapse = ","),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(.data$n_targets))

log_messages("Top 10 TFs by target count:")
print(head(tf_summary, 10))

# =============================================================================
# TF ACTIVITY INFERENCE
# =============================================================================

log_messages("Running TF activity inference with multiple methods...")

# Storage for results
tf_activity_results <- list()
tf_statistics <- list()

# Methods to use (note: decoupleR uses different method names)
inference_methods <- c("viper", "ulm", "mlm")

for (contrast_name in names(edger_results)) {
  log_messages(paste("Processing contrast:", contrast_name))
  
  result <- edger_results[[contrast_name]]
  
  # Prepare gene statistics matrix
  # Create a matrix with genes as rows and statistics as columns
  gene_stats <- result %>%
    dplyr::filter(!is.na(gene_symbol) & gene_symbol != "" & gene_symbol != "NA") %>%
    dplyr::filter(!is.na(logFC) & !is.na(PValue)) %>%
    dplyr::select(gene_symbol, logFC, PValue, FDR) %>%
    dplyr::distinct(gene_symbol, .keep_all = TRUE) %>%
    tibble::column_to_rownames("gene_symbol")
  
  log_messages(paste("Prepared gene statistics for", nrow(gene_stats), "genes"))
  
  # Run TF activity inference with multiple methods
  contrast_activities <- list()
  
  # Convert to matrix format required by decoupleR
  gene_matrix <- as.matrix(gene_stats[, "logFC", drop = FALSE])
  colnames(gene_matrix) <- contrast_name
  
  for (method in inference_methods) {
    log_messages(paste("Running", method, "method..."))
    
    tryCatch({
      # Run decoupleR inference based on method
      if (method == "viper") {
        tf_acts <- decoupleR::run_viper(
          mat = gene_matrix,
          network = dorothea_regulons,
          .source = "tf",
          .target = "target",
          .mor = "mor"
        )
      } else if (method == "ulm") {
        tf_acts <- decoupleR::run_ulm(
          mat = gene_matrix,
          network = dorothea_regulons,
          .source = "tf",
          .target = "target",
          .mor = "mor"
        )
      } else if (method == "mlm") {
        tf_acts <- decoupleR::run_mlm(
          mat = gene_matrix,
          network = dorothea_regulons,
          .source = "tf",
          .target = "target",
          .mor = "mor"
        )
      }
      
      # Store results
      contrast_activities[[method]] <- tf_acts
      
      log_messages(paste("Completed", method, "- found activities for", 
                        nrow(tf_acts), "TFs"))
      
    }, error = function(e) {
      log_messages(paste("Error in", method, ":", e$message))
    })
  }
  
  tf_activity_results[[contrast_name]] <- contrast_activities
}

log_messages("Completed TF activity inference for all contrasts")

# =============================================================================
# STATISTICAL ANALYSIS AND FILTERING
# =============================================================================

log_messages("Performing statistical analysis of TF activities...")

# Function to extract significant TFs
extract_significant_tfs <- function(tf_results, p_cutoff = 0.05, activity_cutoff = 1) {
  significant_tfs <- list()
  
  for (contrast in names(tf_results)) {
    for (method in names(tf_results[[contrast]])) {
      tf_data <- tf_results[[contrast]][[method]]
      
      # Filter for significant TFs
      if (nrow(tf_data) > 0 && all(c("statistic", "score", "p_value") %in% names(tf_data))) {
        sig_tfs <- tf_data %>%
          dplyr::filter(statistic != 0, 
                 abs(score) > activity_cutoff,
                 p_value < p_cutoff) %>%
          dplyr::arrange(p_value)
      } else {
        sig_tfs <- data.frame()
      }
      
      significant_tfs[[paste(contrast, method, sep = "_")]] <- sig_tfs
    }
  }
  
  return(significant_tfs)
}

# Extract significant TFs
significant_tfs <- extract_significant_tfs(tf_activity_results)

# Create summary statistics
tf_summary_stats <- purrr::map_dfr(significant_tfs, function(x) {
  data.frame(
    n_significant = nrow(x),
    n_activated = sum(x$score > 0),
    n_repressed = sum(x$score < 0),
    mean_activity = mean(abs(x$score)),
    median_pvalue = median(x$p_value)
  )
}, .id = "analysis")

log_messages("TF Activity Summary:")
print(tf_summary_stats)

# =============================================================================
# DATA INTEGRATION AND COMPARISON
# =============================================================================

log_messages("Integrating TF activities across contrasts and methods...")

# Create a comprehensive TF activity matrix
create_tf_activity_matrix <- function(tf_results, method = "viper") {
  tf_matrices <- list()
  
  for (contrast in names(tf_results)) {
    if (method %in% names(tf_results[[contrast]])) {
      tf_data <- tf_results[[contrast]][[method]] %>%
        dplyr::select(.data$source, .data$score) %>%
        dplyr::rename(!!rlang::sym(contrast) := .data$score)
      
      tf_matrices[[contrast]] <- tf_data
    }
  }
  
  # Merge all contrasts
  if (length(tf_matrices) > 1) {
    tf_matrix <- purrr::reduce(tf_matrices, dplyr::full_join, by = "source")
  } else if (length(tf_matrices) == 1) {
    tf_matrix <- tf_matrices[[1]]
  } else {
    return(NULL)
  }
  
  # Convert to matrix format
  tf_matrix <- tf_matrix %>%
    tibble::column_to_rownames("source") %>%
    as.matrix()
  
  return(tf_matrix)
}

# Create activity matrices for each method
tf_activity_matrices <- purrr::map(inference_methods, function(method) {
  create_tf_activity_matrix(tf_activity_results, method)
})
names(tf_activity_matrices) <- inference_methods

# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================

# Function to create TF activity heatmaps
create_tf_heatmap <- function(tf_matrix, title, filename) {
  if (is.null(tf_matrix) || nrow(tf_matrix) == 0) {
    log_messages(paste("No data for heatmap:", title))
    return(NULL)
  }
  
  # Filter for TFs with substantial activity
  tf_filtered <- tf_matrix[apply(abs(tf_matrix), 1, max, na.rm = TRUE) > 1, , drop = FALSE]
  
  if (nrow(tf_filtered) == 0) {
    log_messages(paste("No TFs pass activity threshold for:", title))
    return(NULL)
  }
  
  # Create heatmap - PDF version
  pdf(file.path(output_dirs$plots, paste0(filename, ".pdf")), 
      width = 12, height = max(8, nrow(tf_filtered) * 0.3))
  
  heatmap_obj <- pheatmap::pheatmap(
    tf_filtered,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    scale = "none",
    color = colorRampPalette(c("blue", "white", "red"))(100),
    breaks = seq(-max(abs(tf_filtered), na.rm = TRUE), 
                 max(abs(tf_filtered), na.rm = TRUE), 
                 length.out = 101),
    main = title,
    fontsize = 10,
    fontsize_row = 8,
    fontsize_col = 10,
    angle_col = 45
  )
  
  dev.off()
  
  # Create heatmap - PNG version
  png(file.path(output_dirs$plots, paste0(filename, ".png")), 
      width = 12, height = max(8, nrow(tf_filtered) * 0.3), units = "in", res = 300)
  
  pheatmap::pheatmap(
    tf_filtered,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    scale = "none",
    color = colorRampPalette(c("blue", "white", "red"))(100),
    breaks = seq(-max(abs(tf_filtered), na.rm = TRUE), 
                 max(abs(tf_filtered), na.rm = TRUE), 
                 length.out = 101),
    main = title,
    fontsize = 10,
    fontsize_row = 8,
    fontsize_col = 10,
    angle_col = 45
  )
  
  dev.off()
  
  log_messages(paste("Created heatmap:", filename))
  return(heatmap_obj)
}

# Function to create TF activity volcano plots
create_tf_volcano <- function(tf_data, contrast_name, method, filename) {
  if (nrow(tf_data) == 0) return(NULL)
  
  # Add significance categories
  tf_data <- tf_data %>%
    dplyr::mutate(
      significance = dplyr::case_when(
        .data$p_value < 0.001 & abs(.data$score) > 2 ~ "Highly Significant",
        .data$p_value < 0.05 & abs(.data$score) > 1 ~ "Significant", 
        .data$p_value < 0.1 ~ "Trend",
        TRUE ~ "Not Significant"
      ),
      log_p = -log10(.data$p_value),
      label = ifelse(.data$p_value < 0.05 & abs(.data$score) > 1, .data$source, "")
    )
  
  # Create volcano plot
  p <- ggplot(tf_data, aes(x = .data$score, y = .data$log_p)) +
    geom_point(aes(color = .data$significance), alpha = 0.7, size = 2) +
    ggrepel::geom_text_repel(aes(label = .data$label), size = 3, max.overlaps = 20) +
    scale_color_manual(values = c(
      "Highly Significant" = "#d73027",
      "Significant" = "#fc8d59", 
      "Trend" = "#fee08b",
      "Not Significant" = "#999999"
    )) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
    labs(
      title = paste("TF Activity Volcano Plot:", contrast_name, "-", method),
      x = "TF Activity Score",
      y = "-log10(p-value)",
      color = "Significance"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "bottom"
    )
  
  ggsave(file.path(output_dirs$plots, paste0(filename, ".pdf")), 
         plot = p, width = 10, height = 8, device = "pdf")
  
  ggsave(file.path(output_dirs$plots, paste0(filename, ".png")), 
         plot = p, width = 10, height = 8, device = "png", dpi = 300)
  
  log_messages(paste("Created volcano plot:", filename))
  return(p)
}

# Function to create TF activity barplots
create_tf_barplot <- function(tf_data, contrast_name, method, top_n = 20, filename) {
  if (nrow(tf_data) == 0) return(NULL)
  
  # Get top activated and repressed TFs
  top_activated <- tf_data %>%
    dplyr::filter(.data$score > 0, .data$p_value < 0.05) %>%
    dplyr::top_n(top_n/2, .data$score) %>%
    dplyr::mutate(direction = "Activated")
  
  top_repressed <- tf_data %>%
    dplyr::filter(.data$score < 0, .data$p_value < 0.05) %>%
    dplyr::top_n(top_n/2, -.data$score) %>%
    dplyr::mutate(direction = "Repressed")
  
  plot_data <- dplyr::bind_rows(top_activated, top_repressed) %>%
    dplyr::mutate(source = factor(.data$source, levels = .data$source[order(.data$score)]))
  
  if (nrow(plot_data) == 0) {
    log_messages(paste("No significant TFs for barplot:", contrast_name, method))
    return(NULL)
  }
  
  # Create barplot
  p <- ggplot(plot_data, aes(x = .data$source, y = .data$score, fill = .data$direction)) +
    geom_col() +
    scale_fill_manual(values = c("Activated" = "#e74c3c", "Repressed" = "#3498db")) +
    coord_flip() +
    labs(
      title = paste("Top TF Activities:", contrast_name, "-", method),
      x = "Transcription Factor",
      y = "Activity Score",
      fill = "Direction"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.y = element_text(size = 8)
    )
  
  ggplot2::ggsave(file.path(output_dirs$plots, paste0(filename, ".pdf")), 
         plot = p, width = 10, height = max(6, nrow(plot_data) * 0.3), device = "pdf")
  
  ggplot2::ggsave(file.path(output_dirs$plots, paste0(filename, ".png")), 
         plot = p, width = 10, height = max(6, nrow(plot_data) * 0.3), device = "png", dpi = 300)
  
  log_messages(paste("Created barplot:", filename))
  return(p)
}

# =============================================================================
# GENERATE VISUALIZATIONS
# =============================================================================

log_messages("Creating TF activity visualizations...")

# 1. Create heatmaps for each method
for (method in names(tf_activity_matrices)) {
  if (!is.null(tf_activity_matrices[[method]])) {
    create_tf_heatmap(
      tf_activity_matrices[[method]], 
      paste("TF Activity Heatmap -", toupper(method)),
      paste0("tf_activity_heatmap_", method)
    )
  }
}

# 2. Create volcano plots and barplots for each contrast and method
for (contrast in names(tf_activity_results)) {
  for (method in names(tf_activity_results[[contrast]])) {
    tf_data <- tf_activity_results[[contrast]][[method]]
    
    # Volcano plot
    create_tf_volcano(
      tf_data, contrast, method,
      paste0("tf_volcano_", contrast, "_", method)
    )
    
    # Barplot
    create_tf_barplot(
      tf_data, contrast, method, 20,
      paste0("tf_barplot_", contrast, "_", method)
    )
  }
}

# =============================================================================
# NETWORK ANALYSIS
# =============================================================================

log_messages("Performing TF network analysis...")

# Function to create TF-target network
create_tf_network <- function(significant_tfs_list, dorothea_regulons, 
                             min_targets = 5, top_tfs = 30) {
  
  # Get top TFs across all analyses
  if (length(significant_tfs_list) > 0) {
    all_sig_tfs <- purrr::map_dfr(significant_tfs_list, identity, .id = "analysis")
    
    if (nrow(all_sig_tfs) > 0 && "source" %in% names(all_sig_tfs)) {
      all_sig_tfs <- all_sig_tfs %>%
        dplyr::group_by(source) %>%
        dplyr::summarise(
          mean_activity = mean(abs(score)),
          max_activity = max(abs(score)),
          n_analyses = dplyr::n(),
          mean_pvalue = mean(p_value),
          .groups = "drop"
        ) %>%
        dplyr::filter(n_analyses >= 1) %>%  # Must be significant in at least 1 analysis
        dplyr::top_n(top_tfs, mean_activity)
    } else {
      all_sig_tfs <- data.frame()
    }
  } else {
    all_sig_tfs <- data.frame()
  }
  
  if (nrow(all_sig_tfs) == 0) {
    log_messages("No TFs found for network analysis")
    return(NULL)
  }
  
  # Filter DoRothEA for these TFs
  if (nrow(all_sig_tfs) > 0) {
    network_regulons <- dorothea_regulons %>%
      dplyr::filter(tf %in% all_sig_tfs$source) %>%
      dplyr::group_by(tf) %>%
      dplyr::filter(dplyr::n() >= min_targets) %>%  # Minimum target requirement
      dplyr::ungroup()
  } else {
    network_regulons <- data.frame()
  }
  
  # Create igraph network
  if (nrow(network_regulons) == 0) {
    log_messages("No regulons pass network filters")
    return(NULL)
  }
  
  network <- igraph::graph_from_data_frame(
    network_regulons[, c("tf", "target")], 
    directed = TRUE
  )
  
  # Add TF activity as vertex attributes
  igraph::V(network)$type <- ifelse(igraph::V(network)$name %in% all_sig_tfs$source, "TF", "Target")
  
  if (nrow(all_sig_tfs) > 0) {
    tf_activities <- all_sig_tfs %>%
      dplyr::select(source, mean_activity) %>%
      tibble::deframe()
  } else {
    tf_activities <- c()
  }
  
  igraph::V(network)$activity <- ifelse(
    igraph::V(network)$name %in% names(tf_activities),
    tf_activities[igraph::V(network)$name],
    0
  )
  
  log_messages(paste("Created network with", igraph::vcount(network), "nodes and", 
                    igraph::ecount(network), "edges"))
  
  return(list(network = network, tf_data = all_sig_tfs))
}

# Create TF network
tf_network_result <- create_tf_network(significant_tfs, dorothea_regulons)

# Function to plot TF network
plot_tf_network <- function(network_result, filename) {
  if (is.null(network_result)) return(NULL)
  
  network <- network_result$network
  
  # Set layout
  layout <- igraph::layout_with_fr(network)
  
  # Set node properties
  igraph::V(network)$size <- ifelse(igraph::V(network)$type == "TF", 8, 4)
  igraph::V(network)$color <- ifelse(igraph::V(network)$type == "TF", "#e74c3c", "#95a5a6")
  igraph::V(network)$label.cex <- ifelse(igraph::V(network)$type == "TF", 0.8, 0.6)
  
  # Plot network - PDF version
  pdf(file.path(output_dirs$networks, paste0(filename, ".pdf")), 
      width = 16, height = 12)
  
  plot(network,
       layout = layout,
       vertex.label = ifelse(igraph::V(network)$type == "TF", igraph::V(network)$name, ""),
       vertex.label.color = "black",
       vertex.label.dist = 1,
       edge.arrow.size = 0.3,
       edge.color = "#7f8c8d",
       main = "TF-Target Regulatory Network",
       sub = paste("Top TFs from OUD Analysis")
  )
  
  # Add legend
  legend("bottomleft", 
         legend = c("Transcription Factor", "Target Gene"),
         col = c("#e74c3c", "#95a5a6"),
         pch = 19,
         cex = 0.8)
  
  dev.off()
  
  # Plot network - PNG version
  png(file.path(output_dirs$networks, paste0(filename, ".png")), 
      width = 16, height = 12, units = "in", res = 300)
  
  plot(network,
       layout = layout,
       vertex.label = ifelse(igraph::V(network)$type == "TF", igraph::V(network)$name, ""),
       vertex.label.color = "black",
       vertex.label.dist = 1,
       edge.arrow.size = 0.3,
       edge.color = "#7f8c8d",
       main = "TF-Target Regulatory Network",
       sub = paste("Top TFs from OUD Analysis")
  )
  
  # Add legend
  legend("bottomleft", 
         legend = c("Transcription Factor", "Target Gene"),
         col = c("#e74c3c", "#95a5a6"),
         pch = 19,
         cex = 0.8)
  
  dev.off()
  
  log_messages(paste("Created network plot:", filename))
}

# Plot TF network
if (!is.null(tf_network_result)) {
  plot_tf_network(tf_network_result, "tf_regulatory_network")
}

# =============================================================================
# COMPARATIVE ANALYSIS
# =============================================================================

log_messages("Performing comparative analysis across contrasts...")

# Function to compare TF activities between contrasts
compare_tf_activities <- function(tf_results, method = "viper") {
  # Extract TF activities for the specified method
  tf_activities <- purrr::map(tf_results, function(x) {
    if (method %in% names(x) && nrow(x[[method]]) > 0) {
      x[[method]] %>% dplyr::select(source, score, p_value)
    } else {
      NULL
    }
  })
  
  # Remove NULL results
  tf_activities <- tf_activities[!purrr::map_lgl(tf_activities, is.null)]
  
  if (length(tf_activities) < 2) {
    log_messages("Not enough contrasts for comparison")
    return(NULL)
  }
  
  # Find common TFs
  common_tfs <- Reduce(intersect, purrr::map(tf_activities, ~ .x$source))
  
  if (length(common_tfs) == 0) {
    log_messages("No common TFs found across contrasts")
    return(NULL)
  }
  
  # Create comparison matrix
  comparison_data <- purrr::map_dfr(tf_activities, function(x) {
    x %>% dplyr::filter(source %in% common_tfs)
  }, .id = "contrast")
  
  return(comparison_data)
}

# Compare TF activities
tf_comparison <- compare_tf_activities(tf_activity_results)

# Create correlation plots if comparison data exists
if (!is.null(tf_comparison)) {
  
  # Reshape for correlation analysis
  tf_wide <- tf_comparison %>%
    dplyr::select(contrast, source, score) %>%
    tidyr::pivot_wider(names_from = contrast, values_from = score) %>%
    tibble::column_to_rownames("source")
  
  # Calculate correlations
  tf_correlations <- cor(tf_wide, use = "complete.obs")
  
  # Plot correlation heatmap - PDF version
  pdf(file.path(output_dirs$plots, "tf_activity_correlations.pdf"), 
      width = 10, height = 8)
  
  corrplot::corrplot(tf_correlations, 
           method = "color",
           order = "hclust",
           tl.cex = 0.8,
           title = "TF Activity Correlations Across Contrasts",
           mar = c(0,0,2,0))
  
  dev.off()
  
  # Plot correlation heatmap - PNG version
  png(file.path(output_dirs$plots, "tf_activity_correlations.png"), 
      width = 10, height = 8, units = "in", res = 300)
  
  corrplot::corrplot(tf_correlations, 
           method = "color",
           order = "hclust",
           tl.cex = 0.8,
           title = "TF Activity Correlations Across Contrasts",
           mar = c(0,0,2,0))
  
  dev.off()
  
  log_messages("Created TF activity correlation plot")
}

# =============================================================================
# OUTPUT GENERATION
# =============================================================================

log_messages("Generating output files and reports...")

# Save TF activity matrices
for (method in names(tf_activity_matrices)) {
  if (!is.null(tf_activity_matrices[[method]])) {
    write.csv(tf_activity_matrices[[method]], 
              file.path(output_dirs$tables, paste0("tf_activity_matrix_", method, ".csv")))
  }
}

# Save significant TFs
significant_tfs_df <- purrr::map_dfr(significant_tfs, identity, .id = "analysis")
write.csv(significant_tfs_df, 
          file.path(output_dirs$tables, "significant_tfs_all_analyses.csv"),
          row.names = FALSE)

# Save summary statistics
write.csv(tf_summary_stats, 
          file.path(output_dirs$tables, "tf_activity_summary_stats.csv"),
          row.names = TRUE)

# Save network data if available
if (!is.null(tf_network_result)) {
  write.csv(tf_network_result$tf_data, 
            file.path(output_dirs$networks, "network_tf_data.csv"),
            row.names = FALSE)
}

# =============================================================================
# GENERATE HTML REPORT
# =============================================================================

log_messages("Generating HTML summary report...")

# Create summary report
create_html_report <- function() {
  
  # Count significant TFs per analysis
  if (nrow(significant_tfs_df) > 0) {
    sig_counts <- significant_tfs_df %>%
      dplyr::group_by(analysis) %>%
      dplyr::summarise(
        total_tfs = dplyr::n(),
        activated = sum(score > 0),
        repressed = sum(score < 0),
        .groups = "drop"
      )
    
    # Top TFs across all analyses
    top_tfs_global <- significant_tfs_df %>%
      dplyr::group_by(source) %>%
      dplyr::summarise(
        n_analyses = dplyr::n(),
        mean_activity = mean(abs(score)),
        max_activity = max(abs(score)),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(mean_activity)) %>%
      head(20)
  } else {
    sig_counts <- data.frame(analysis = character(0), total_tfs = integer(0), 
                           activated = integer(0), repressed = integer(0))
    top_tfs_global <- data.frame(source = character(0), n_analyses = integer(0),
                               mean_activity = numeric(0), max_activity = numeric(0))
  }
  
  # Create HTML content
  html_content <- paste0(
    "<html><head><title>TF Activity Analysis Report</title>",
    "<style>body{font-family: Arial, sans-serif; margin: 40px;}",
    "table{border-collapse: collapse; width: 100%;}",
    "th, td{border: 1px solid #ddd; padding: 8px; text-align: left;}",
    "th{background-color: #f2f2f2;}</style></head><body>",
    
    "<h1>Transcription Factor Activity Analysis Report</h1>",
    "<h2>Analysis Summary</h2>",
    "<p>Generated on: ", Sys.time(), "</p>",
    "<p>Method: DoRothEA + DecoupleR</p>",
    "<p>Total contrasts analyzed: ", length(edger_results), "</p>",
    "<p>TF confidence levels: A, B, C</p>",
    
    "<h2>Significant TFs per Analysis</h2>",
    "<table>",
    "<tr><th>Analysis</th><th>Total TFs</th><th>Activated</th><th>Repressed</th></tr>"
  )
  
  if (nrow(sig_counts) > 0) {
    for (i in seq_len(nrow(sig_counts))) {
      html_content <- paste0(html_content,
        "<tr><td>", sig_counts$analysis[i], "</td>",
        "<td>", sig_counts$total_tfs[i], "</td>",
        "<td>", sig_counts$activated[i], "</td>",
        "<td>", sig_counts$repressed[i], "</td></tr>"
      )
    }
  } else {
    html_content <- paste0(html_content,
      "<tr><td colspan='4'>No significant TFs found</td></tr>"
    )
  }
  
  html_content <- paste0(html_content, "</table>")
  
  html_content <- paste0(html_content,
    "<h2>Top 20 TFs by Mean Activity</h2>",
    "<table>",
    "<tr><th>TF</th><th>N Analyses</th><th>Mean Activity</th><th>Max Activity</th></tr>"
  )
  
  if (nrow(top_tfs_global) > 0) {
    for (i in seq_len(nrow(top_tfs_global))) {
      html_content <- paste0(html_content,
        "<tr><td>", top_tfs_global$source[i], "</td>",
        "<td>", top_tfs_global$n_analyses[i], "</td>",
        "<td>", round(top_tfs_global$mean_activity[i], 3), "</td>",
        "<td>", round(top_tfs_global$max_activity[i], 3), "</td></tr>"
      )
    }
  } else {
    html_content <- paste0(html_content,
      "<tr><td colspan='4'>No significant TFs found</td></tr>"
    )
  }
  
  html_content <- paste0(html_content, "</table></body></html>")
  
  # Write HTML file
  writeLines(html_content, file.path(output_dirs$main, "tf_activity_report.html"))
  
  log_messages("Created HTML summary report")
}

create_html_report()

# =============================================================================
# FINAL SUMMARY
# =============================================================================

log_messages("TF Activity Analysis Complete!")
log_messages("=================================")
log_messages(paste("Total contrasts analyzed:", length(edger_results)))
log_messages(paste("Total significant TF-contrast combinations:", nrow(significant_tfs_df)))
log_messages(paste("Unique TFs identified:", length(unique(significant_tfs_df$source))))
log_messages(paste("Output directory:", output_dirs$main))

# Print top 10 most frequently significant TFs
if (nrow(significant_tfs_df) > 0) {
  top_frequent_tfs <- significant_tfs_df %>%
    dplyr::count(source, sort = TRUE) %>%
    head(10)

  log_messages("Top 10 most frequently significant TFs:")
  for (i in seq_len(nrow(top_frequent_tfs))) {
    log_messages(paste(i, ".", top_frequent_tfs$source[i], "(", top_frequent_tfs$n[i], "analyses)"))
  }
} else {
  log_messages("No significant TFs found across all analyses")
}

# Session info
log_messages("Session Info:")
log_messages(paste("R version:", R.version.string))
log_messages(paste("decoupleR version:", packageVersion("decoupleR")))

cat("\n🎉 TF Activity Analysis completed successfully!\n")
cat("📁 Results saved in:", output_dirs$main, "\n")
cat("📊 Check the HTML report for a summary of findings\n")
cat("📈 Plots available in:", output_dirs$plots, "\n")
cat("📋 Tables available in:", output_dirs$tables, "\n")

# End of script