# =============================================================================
# GSEA Analysis Utility Functions
# =============================================================================
#
# Purpose: Utility functions for GSEA pathway analysis
# Author: Multi-Omics OUD Study
# Date: 2024
# =============================================================================

# =============================================================================
# PACKAGE MANAGEMENT FUNCTIONS
# =============================================================================

#' Install and load required packages
#' @param packages Character vector of package names
#' @param source Package source ("bioc", "cran", or "auto")
load_packages <- function(packages, source = "auto") {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat("📦 Installing", pkg, "...\n")

      if (source == "bioc" || pkg %in% c(
        "clusterProfiler", "org.Hs.eg.db", "msigdbr", "enrichplot",
        "DOSE", "ReactomePA", "pathview"
      )) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
        }
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
      } else {
        install.packages(pkg, dependencies = TRUE)
      }
    }

    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    cat("✅", pkg, "loaded\n")
  }
}

#' Check if all required packages are available
check_packages <- function(required_packages) {
  missing <- c()
  available <- c()

  all_packages <- c(required_packages$bioc, required_packages$cran)

  for (pkg in all_packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      available <- c(available, pkg)
    } else {
      missing <- c(missing, pkg)
    }
  }

  list(
    available = available,
    missing = missing,
    all_available = length(missing) == 0
  )
}

# =============================================================================
# FILE AND DIRECTORY MANAGEMENT
# =============================================================================

#' Create output directories safely
#' @param paths List of directory paths to create
create_output_dirs <- function(paths) {
  for (path in paths) {
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE, showWarnings = FALSE)
      cat("📁 Created directory:", basename(path), "\n")
    }
  }
}

#' Safe data loading with error handling
#' @param file_path Path to the data file
#' @param show_summary Logical, whether to show data summary
safe_load_data <- function(file_path, show_summary = TRUE) {
  if (!file.exists(file_path)) {
    cat("⚠️  File not found:", basename(file_path), "\n")
    return(NULL)
  }

  tryCatch(
    {
      data <- readr::read_csv(file_path, show_col_types = FALSE)

      if (nrow(data) == 0) {
        cat("⚠️  Empty data file:", basename(file_path), "\n")
        return(NULL)
      }

      if (show_summary) {
        cat("✅ Loaded:", basename(file_path), "(", nrow(data), "rows,", ncol(data), "cols)\n")
      }

      return(data)
    },
    error = function(e) {
      cat("❌ Error loading", basename(file_path), ":", e$message, "\n")
      return(NULL)
    }
  )
}

#' Generate timestamped filename
#' @param prefix File prefix
#' @param suffix File suffix/extension
#' @param use_timestamp Logical, whether to include timestamp
timestamped_filename <- function(prefix, suffix, use_timestamp = TRUE) {
  if (use_timestamp) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    paste0(prefix, "_", timestamp, ".", suffix)
  } else {
    paste0(prefix, ".", suffix)
  }
}

# =============================================================================
# LOGGING FUNCTIONS
# =============================================================================

#' Setup logging system
#' @param log_dir Directory for log files
#' @param log_level Logging level
setup_logging <- function(log_dir, log_level = "INFO") {
  log_file <- file.path(log_dir, timestamped_filename("gsea_pipeline", "log"))
  log_con <- file(log_file, "w")

  list(
    log_file = log_file,
    log_connection = log_con
  )
}

#' Log message with timestamp
#' @param msg Message to log
#' @param level Log level (DEBUG, INFO, WARNING, ERROR)
#' @param log_con Log file connection (optional)
log_message <- function(msg, level = "INFO", log_con = NULL) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- paste0("[", timestamp, "] [", level, "] ", msg)

  cat(log_entry, "\n")

  if (!is.null(log_con)) {
    writeLines(log_entry, log_con)
    flush(log_con)
  }
}

# =============================================================================
# DATA PROCESSING FUNCTIONS
# =============================================================================

#' Validate differential expression data
#' @param data Differential expression data frame
#' @param required_cols Required column names
validate_de_data <- function(data, required_cols = c("gene", "effect_size", "pval", "adj_pval")) {
  if (is.null(data) || nrow(data) == 0) {
    return(list(valid = FALSE, message = "Empty or null data"))
  }

  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    return(list(
      valid = FALSE,
      message = paste("Missing columns:", paste(missing_cols, collapse = ", "))
    ))
  }

  # Check for missing values in critical columns
  critical_na <- sapply(required_cols, function(col) sum(is.na(data[[col]])))
  if (any(critical_na > 0)) {
    return(list(
      valid = FALSE,
      message = paste("Missing values found in:", names(critical_na)[critical_na > 0])
    ))
  }

  list(valid = TRUE, message = "Data validation passed")
}

#' Clean and prepare gene list for analysis
#' @param effect_sizes Named vector of effect sizes
#' @param remove_na Logical, whether to remove NA values
prepare_gene_list <- function(effect_sizes, remove_na = TRUE) {
  if (remove_na) {
    effect_sizes <- effect_sizes[!is.na(effect_sizes)]
  }

  # Ensure numeric and sort
  effect_sizes <- as.numeric(effect_sizes)
  effect_sizes <- sort(effect_sizes, decreasing = TRUE)

  return(effect_sizes)
}

#' Convert gene symbols to Entrez IDs
#' @param gene_symbols Character vector of gene symbols
#' @param org_db Organism database (default: org.Hs.eg.db)
convert_gene_ids <- function(gene_symbols, org_db = org.Hs.eg.db::org.Hs.eg.db) {
  tryCatch(
    {
      conversion <- clusterProfiler::bitr(
        gene_symbols,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = org_db
      )

      conversion_rate <- nrow(conversion) / length(gene_symbols) * 100

      list(
        conversion_table = conversion,
        conversion_rate = conversion_rate,
        success = TRUE
      )
    },
    error = function(e) {
      list(
        conversion_table = NULL,
        conversion_rate = 0,
        success = FALSE,
        error = e$message
      )
    }
  )
}

# =============================================================================
# PATHWAY FILTERING FUNCTIONS
# =============================================================================

#' Filter pathways for neuroscience relevance
#' @param pathway_data MSigDB pathway data
#' @param keywords Vector of filtering keywords
filter_neuro_pathways <- function(pathway_data, keywords) {
  if (is.null(pathway_data) || nrow(pathway_data) == 0) {
    return(NULL)
  }

  pattern <- paste(keywords, collapse = "|")
  filtered <- pathway_data[grepl(pattern, pathway_data$gs_name, ignore.case = TRUE), ]

  return(filtered)
}

#' Load and filter MSigDB collections
#' @param collections List of collection specifications
#' @param neuro_keywords Keywords for filtering
#' @param species Species name (default: "Homo sapiens")
load_msigdb_collections <- function(collections, neuro_keywords, species = "Homo sapiens") {
  result <- list()

  for (name in names(collections)) {
    collection_info <- collections[[name]]

    tryCatch(
      {
        if (is.null(collection_info$subcollection)) {
          msigdb_data <- msigdbr::msigdbr(
            species = species,
            collection = collection_info$collection
          )
        } else {
          msigdb_data <- msigdbr::msigdbr(
            species = species,
            collection = collection_info$collection,
            subcollection = collection_info$subcollection
          )
        }

        # Filter for neuroscience relevance
        filtered_data <- filter_neuro_pathways(msigdb_data, neuro_keywords)

        if (!is.null(filtered_data) && nrow(filtered_data) > 0) {
          result[[name]] <- filtered_data
          log_message(paste(
            collection_info$name, "loaded:",
            nrow(filtered_data), "neuro-relevant pathways"
          ))
        } else {
          log_message(paste(collection_info$name, "- no relevant pathways found"), "WARNING")
        }
      },
      error = function(e) {
        log_message(paste("Error loading", collection_info$name, ":", e$message), "ERROR")
      }
    )
  }

  return(result)
}

# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

#' Run GSEA analysis safely
#' @param gene_list Named numeric vector (gene symbols -> effect sizes)
#' @param gene_sets Gene set data frame
#' @param analysis_name Name of the analysis
#' @param params Analysis parameters
run_gsea_safe <- function(gene_list, gene_sets, analysis_name, params) {
  if (is.null(gene_sets) || nrow(gene_sets) == 0) {
    log_message(paste("No gene sets available for", analysis_name), "WARNING")
    return(NULL)
  }

  tryCatch(
    {
      result <- clusterProfiler::GSEA(
        geneList = gene_list,
        TERM2GENE = gene_sets[, c("gs_name", "entrez_gene")],
        pvalueCutoff = params$pvalue_cutoff,
        pAdjustMethod = "BH",
        minGSSize = params$min_gs_size,
        maxGSSize = params$max_gs_size,
        eps = params$eps,
        verbose = FALSE
      )

      if (!is.null(result) && nrow(result@result) > 0) {
        log_message(paste(analysis_name, "GSEA: found", nrow(result@result), "enriched pathways"))
        return(result)
      } else {
        log_message(paste(analysis_name, "GSEA: no significant pathways"), "WARNING")
        return(NULL)
      }
    },
    error = function(e) {
      log_message(paste(analysis_name, "GSEA error:", e$message), "ERROR")
      return(NULL)
    }
  )
}

#' Run Over-Representation Analysis safely
#' @param genes Character vector of gene symbols
#' @param universe_genes Character vector of universe genes
#' @param gene_sets Gene set data frame
#' @param analysis_name Name of the analysis
#' @param params Analysis parameters
run_ora_safe <- function(genes, universe_genes, gene_sets, analysis_name, params) {
  if (is.null(gene_sets) || nrow(gene_sets) == 0) {
    log_message(paste("No gene sets available for", analysis_name), "WARNING")
    return(NULL)
  }

  tryCatch(
    {
      result <- clusterProfiler::enricher(
        gene = genes,
        universe = universe_genes,
        TERM2GENE = gene_sets[, c("gs_name", "gene_symbol")],
        pvalueCutoff = params$pvalue_cutoff,
        pAdjustMethod = "BH",
        minGSSize = params$min_gs_size,
        maxGSSize = params$max_gs_size
      )

      if (!is.null(result) && nrow(result@result) > 0) {
        log_message(paste(analysis_name, "ORA: found", nrow(result@result), "enriched pathways"))
        return(result)
      } else {
        log_message(paste(analysis_name, "ORA: no significant pathways"), "WARNING")
        return(NULL)
      }
    },
    error = function(e) {
      log_message(paste(analysis_name, "ORA error:", e$message), "ERROR")
      return(NULL)
    }
  )
}

# =============================================================================
# PLOTTING UTILITY FUNCTIONS
# =============================================================================

#' Clean pathway names for display
#' @param names Character vector of pathway names
#' @param max_length Maximum character length
clean_pathway_names <- function(names, max_length = 50) {
  # Remove common prefixes
  names <- stringr::str_replace(names, "^HALLMARK_", "")
  names <- stringr::str_replace(names, "^GOBP_", "")
  names <- stringr::str_replace(names, "^GOCC_", "")
  names <- stringr::str_replace(names, "^GOMF_", "")
  names <- stringr::str_replace(names, "^KEGG_", "")
  names <- stringr::str_replace(names, "^REACTOME_", "")

  # Replace underscores with spaces
  names <- stringr::str_replace_all(names, "_", " ")

  # Truncate if too long
  names <- ifelse(
    nchar(names) > max_length,
    paste0(substr(names, 1, max_length - 3), "..."),
    names
  )

  return(names)
}

#' Get appropriate color palette
#' @param n Number of colors needed
#' @param type Type of palette ("continuous", "categorical", "diverging")
get_color_palette <- function(n, type = "categorical") {
  if (type == "continuous") {
    return(viridis::viridis(n))
  } else if (type == "diverging") {
    return(RColorBrewer::brewer.pal(min(n, 11), "RdBu"))
  } else {
    # Categorical
    if (n <= 8) {
      return(RColorBrewer::brewer.pal(max(3, n), "Set2"))
    } else if (n <= 12) {
      return(RColorBrewer::brewer.pal(n, "Set3"))
    } else {
      return(viridis::viridis(n))
    }
  }
}

#' Save plot safely with error handling
#' @param plot_obj ggplot object
#' @param filename Output filename
#' @param width Plot width
#' @param height Plot height
#' @param dpi Resolution
save_plot_safe <- function(plot_obj, filename, width = 12, height = 8, dpi = 300) {
  if (is.null(plot_obj)) {
    return(FALSE)
  }

  tryCatch(
    {
      ggplot2::ggsave(filename, plot = plot_obj, width = width, height = height, dpi = dpi)
      log_message(paste("Saved plot:", basename(filename)))
      return(TRUE)
    },
    error = function(e) {
      log_message(paste("Plot save error for", basename(filename), ":", e$message), "ERROR")
      return(FALSE)
    }
  )
}

# =============================================================================
# REPORT GENERATION FUNCTIONS
# =============================================================================

#' Generate pathway analysis summary
#' @param results_list List of analysis results
#' @param gene_info Gene information list
generate_analysis_summary <- function(results_list, gene_info) {
  summary <- list(
    total_genes = gene_info$total_genes,
    converted_genes = gene_info$converted_genes,
    conversion_rate = gene_info$conversion_rate,
    databases_analyzed = length(results_list),
    total_pathways = 0,
    significant_pathways = 0
  )

  for (result in results_list) {
    if (!is.null(result)) {
      summary$total_pathways <- summary$total_pathways + nrow(result@result)
      summary$significant_pathways <- summary$significant_pathways +
        sum(result@result$p.adjust < 0.05)
    }
  }

  return(summary)
}

#' Create markdown report
#' @param summary Analysis summary
#' @param output_file Output file path
#' @param include_details Logical, whether to include detailed results
create_markdown_report <- function(summary, output_file, include_details = TRUE) {
  report_content <- paste0(
    "# GSEA Analysis Report\n\n",
    "**Analysis Date:** ", format(Sys.Date(), "%B %d, %Y"), "\n",
    "**Generated:** ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n",
    "## Summary Statistics\n\n",
    "- **Total Input Genes:** ", summary$total_genes, "\n",
    "- **Successfully Converted:** ", summary$converted_genes,
    " (", round(summary$conversion_rate, 1), "%)\n",
    "- **Databases Analyzed:** ", summary$databases_analyzed, "\n",
    "- **Total Pathways Found:** ", summary$total_pathways, "\n",
    "- **Significant Pathways:** ", summary$significant_pathways, "\n\n"
  )

  writeLines(report_content, output_file)
  log_message(paste("Report saved to:", basename(output_file)))
}

# =============================================================================
# VALIDATION FUNCTIONS
# =============================================================================

#' Validate analysis configuration
#' @param config Configuration list
validate_analysis_config <- function(config) {
  errors <- c()

  # Check statistical parameters
  if (config$STATS$pvalue_cutoff <= 0 || config$STATS$pvalue_cutoff >= 1) {
    errors <- c(errors, "Invalid p-value cutoff")
  }

  if (config$STATS$min_gs_size <= 0) {
    errors <- c(errors, "Minimum gene set size must be positive")
  }

  # Check paths
  if (!dir.exists(dirname(config$PATHS$input_file))) {
    errors <- c(errors, "Input directory does not exist")
  }

  # Check plot parameters
  if (config$PLOT_PARAMS$width <= 0 || config$PLOT_PARAMS$height <= 0) {
    errors <- c(errors, "Invalid plot dimensions")
  }

  if (length(errors) > 0) {
    stop("Configuration validation failed:\n", paste(errors, collapse = "\n"))
  }

  return(TRUE)
}

#' System requirements check
check_system_requirements <- function() {
  requirements <- list(
    r_version = numeric_version(paste(R.version$major, R.version$minor, sep = ".")),
    min_r_version = numeric_version("4.0.0"),
    memory_gb = as.numeric(object.size(1:1e6)) * 1000 / 1024^3,
    min_memory_gb = 4
  )

  issues <- c()

  if (requirements$r_version < requirements$min_r_version) {
    issues <- c(issues, paste("R version", requirements$r_version, "< required", requirements$min_r_version))
  }

  if (requirements$memory_gb < requirements$min_memory_gb) {
    issues <- c(issues, paste("Available memory", round(requirements$memory_gb, 1), "GB < required", requirements$min_memory_gb, "GB"))
  }

  list(
    passed = length(issues) == 0,
    issues = issues,
    requirements = requirements
  )
}

# =============================================================================
# UTILITY HELPER FUNCTIONS
# =============================================================================

#' Print progress bar
#' @param current Current progress
#' @param total Total steps
#' @param width Width of progress bar
print_progress <- function(current, total, width = 50) {
  percent <- round(current / total * 100)
  filled <- round(width * current / total)
  bar <- paste0(
    "[", paste(rep("=", filled), collapse = ""),
    paste(rep(" ", width - filled), collapse = ""), "]"
  )
  cat("\r", bar, " ", percent, "% (", current, "/", total, ")")
  if (current == total) cat("\n")
}

#' Calculate enrichment statistics
#' @param result_obj Enrichment result object
calculate_enrichment_stats <- function(result_obj) {
  if (is.null(result_obj) || nrow(result_obj@result) == 0) {
    return(list(
      total_pathways = 0,
      significant_pathways = 0,
      top_pathway = NULL,
      mean_fold_enrichment = NA
    ))
  }

  result_df <- result_obj@result

  list(
    total_pathways = nrow(result_df),
    significant_pathways = sum(result_df$p.adjust < 0.05),
    top_pathway = result_df$Description[1],
    mean_fold_enrichment = mean(as.numeric(result_df$FoldEnrichment), na.rm = TRUE)
  )
}

#' Memory cleanup function
cleanup_memory <- function() {
  invisible(gc())
  cat("🧹 Memory cleanup completed\n")
}

# Message when functions are loaded
cat("🔧 GSEA utility functions loaded successfully\n")
