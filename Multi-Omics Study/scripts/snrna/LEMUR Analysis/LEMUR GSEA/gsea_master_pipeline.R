# =============================================================================
# GSEA Master Pipeline - Consolidated Analysis Script
# =============================================================================
#
# Purpose: Comprehensive Gene Set Enrichment Analysis for LEMUR results
# Features:
#   - Automated setup and validation
#   - Multiple pathway databases (MSigDB, KEGG, Reactome)
#   - Enhanced visualization with publication-ready plots
#   - Comprehensive reporting and documentation
#
# Author: Multi-Omics OUD Study
# Date: 2024
# Version: 2.0 (Consolidated)
# =============================================================================

# =============================================================================
# INITIALIZATION
# =============================================================================

cat("🧬 GSEA Master Pipeline - Consolidated Analysis\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("🕐 Analysis started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# Start timing
pipeline_start_time <- Sys.time()

# Load configuration and functions
source("gsea_config.R")
source("gsea_functions.R")

# Global variables for logging
LOG_CONNECTION <- NULL
ANALYSIS_SUMMARY <- list()

# =============================================================================
# PHASE 1: SETUP AND VALIDATION
# =============================================================================

cat("🔧 PHASE 1: Setup and Validation\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

# Setup logging
setup_analysis_logging <- function() {
  log_setup <- setup_logging(PATHS$log_dir)
  LOG_CONNECTION <<- log_setup$log_connection
  log_message("GSEA Master Pipeline started", "INFO", LOG_CONNECTION)
  log_message(paste("Configuration loaded from:", "gsea_config.R"), "INFO", LOG_CONNECTION)
}

# Create output directories
create_output_dirs(unlist(PATHS))
setup_analysis_logging()

# Validate configuration
tryCatch(
  {
    validate_analysis_config(list(STATS = STATS, PATHS = PATHS, PLOT_PARAMS = PLOT_PARAMS))
    log_message("Configuration validation passed", "INFO", LOG_CONNECTION)
  },
  error = function(e) {
    log_message(paste("Configuration validation failed:", e$message), "ERROR", LOG_CONNECTION)
    stop("❌ Configuration validation failed")
  }
)

# Check system requirements
system_check <- check_system_requirements()
if (!system_check$passed) {
  for (issue in system_check$issues) {
    log_message(issue, "WARNING", LOG_CONNECTION)
  }
  cat("⚠️  System requirement warnings detected (check logs)\n")
}

# Check and load packages
cat("📦 Checking and loading required packages...\n")
package_status <- check_packages(REQUIRED_PACKAGES)

if (!package_status$all_available) {
  cat("🔄 Installing missing packages...\n")
  if (length(REQUIRED_PACKAGES$bioc) > 0) {
    load_packages(REQUIRED_PACKAGES$bioc, "bioc")
  }
  if (length(REQUIRED_PACKAGES$cran) > 0) {
    load_packages(REQUIRED_PACKAGES$cran, "cran")
  }
} else {
  # Load packages silently
  suppressPackageStartupMessages({
    load_packages(c(REQUIRED_PACKAGES$bioc, REQUIRED_PACKAGES$cran))
  })
}

cat("✅ Setup and validation complete\n\n")

# =============================================================================
# PHASE 2: DATA LOADING AND PREPARATION
# =============================================================================

cat("📊 PHASE 2: Data Loading and Preparation\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

# Load differential expression data
de_data <- safe_load_data(PATHS$input_file, show_summary = TRUE)
if (is.null(de_data)) {
  log_message("Failed to load input data", "ERROR", LOG_CONNECTION)
  stop("❌ Input data not found. Please run LEMUR analysis first.")
}

# Validate data structure
validation_result <- validate_de_data(de_data)
if (!validation_result$valid) {
  log_message(paste("Data validation failed:", validation_result$message), "ERROR", LOG_CONNECTION)
  stop("❌ Invalid input data structure")
}

# Data summary
n_genes <- nrow(de_data)
effect_range <- range(de_data$effect_size, na.rm = TRUE)
n_significant <- sum(de_data$adj_pval < STATS$fdr_threshold, na.rm = TRUE)

log_message(paste("Loaded", n_genes, "genes"), "INFO", LOG_CONNECTION)
log_message(paste("Effect size range:", round(effect_range[1], 3), "to", round(effect_range[2], 3)), "INFO", LOG_CONNECTION)
log_message(paste("Significant genes (FDR <", STATS$fdr_threshold, "):", n_significant), "INFO", LOG_CONNECTION)

cat("📊 Dataset summary:\n")
cat("  - Total genes:", n_genes, "\n")
cat("  - Effect size range:", round(effect_range[1], 3), "to", round(effect_range[2], 3), "\n")
cat("  - Significant genes:", n_significant, "\n")

# Prepare gene lists
cat("🧬 Preparing gene lists for analysis...\n")

# Create ranked gene list for GSEA
gene_list <- de_data$effect_size
names(gene_list) <- de_data$gene
gene_list <- prepare_gene_list(gene_list, remove_na = TRUE)

# Extract gene symbols for ORA
gene_symbols <- de_data$gene[!is.na(de_data$effect_size)]

# Convert gene symbols to Entrez IDs
conversion_result <- convert_gene_ids(gene_symbols)

if (!conversion_result$success) {
  log_message(paste("Gene ID conversion failed:", conversion_result$error), "ERROR", LOG_CONNECTION)
  stop("❌ Failed to convert gene IDs")
}

gene_entrez_table <- conversion_result$conversion_table
conversion_rate <- conversion_result$conversion_rate

# Create Entrez ID ranked list for GSEA
gene_list_entrez <- gene_list[gene_entrez_table$SYMBOL]
names(gene_list_entrez) <- gene_entrez_table$ENTREZID
gene_list_entrez <- prepare_gene_list(gene_list_entrez, remove_na = TRUE)

log_message(paste(
  "Gene ID conversion:", nrow(gene_entrez_table), "/", length(gene_symbols),
  "(", round(conversion_rate, 1), "%)"
), "INFO", LOG_CONNECTION)

cat(
  "🔄 Gene ID conversion:", nrow(gene_entrez_table), "/", length(gene_symbols),
  "(", round(conversion_rate, 1), "%)\n"
)

if (conversion_rate < GENE_FILTER$conversion_rate_threshold * 100) {
  log_message("Low gene ID conversion rate", "WARNING", LOG_CONNECTION)
  cat("⚠️  Low gene ID conversion rate\n")
}

cat("✅ Data preparation complete\n\n")

# Store gene information for reporting
ANALYSIS_SUMMARY$genes <- list(
  total_genes = n_genes,
  significant_genes = n_significant,
  converted_genes = nrow(gene_entrez_table),
  conversion_rate = conversion_rate
)

# =============================================================================
# PHASE 3: PATHWAY DATABASE PREPARATION
# =============================================================================

cat("🗃️ PHASE 3: Pathway Database Preparation\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

# Load MSigDB collections
if (WORKFLOW$load_msigdb) {
  cat("📚 Loading MSigDB collections...\n")

  msigdb_collections <- load_msigdb_collections(
    MSIGDB_COLLECTIONS,
    NEURO_KEYWORDS,
    species = "Homo sapiens"
  )

  ANALYSIS_SUMMARY$msigdb <- sapply(msigdb_collections, function(x) {
    if (!is.null(x)) nrow(x) else 0
  })
} else {
  msigdb_collections <- list()
  ANALYSIS_SUMMARY$msigdb <- list()
}

cat("✅ Database preparation complete\n\n")

# =============================================================================
# PHASE 4: GENE SET ENRICHMENT ANALYSIS (GSEA)
# =============================================================================

if (WORKFLOW$run_gsea) {
  cat("🔬 PHASE 4: Gene Set Enrichment Analysis (GSEA)\n")
  cat(paste(rep("-", 40), collapse = ""), "\n")

  gsea_results <- list()

  for (collection_name in names(msigdb_collections)) {
    collection_data <- msigdb_collections[[collection_name]]
    analysis_name <- MSIGDB_COLLECTIONS[[collection_name]]$name

    cat("🔄 Running GSEA for", analysis_name, "...\n")

    result <- run_gsea_safe(
      gene_list_entrez,
      collection_data,
      analysis_name,
      STATS
    )

    if (!is.null(result)) {
      gsea_results[[collection_name]] <- result
    }
  }

  cat("✅ GSEA analysis complete (", length(gsea_results), "successful analyses)\n\n")
  ANALYSIS_SUMMARY$gsea_results <- length(gsea_results)
} else {
  gsea_results <- list()
  ANALYSIS_SUMMARY$gsea_results <- 0
  cat("⏭️  GSEA analysis skipped\n\n")
}

# =============================================================================
# PHASE 5: OVER-REPRESENTATION ANALYSIS (ORA)
# =============================================================================

if (WORKFLOW$run_ora) {
  cat("🎯 PHASE 5: Over-Representation Analysis (ORA)\n")
  cat(paste(rep("-", 40), collapse = ""), "\n")

  # Get universe genes
  all_human_genes <- keys(org.Hs.eg.db, keytype = "SYMBOL")

  ora_results <- list()

  for (collection_name in names(msigdb_collections)) {
    collection_data <- msigdb_collections[[collection_name]]
    analysis_name <- MSIGDB_COLLECTIONS[[collection_name]]$name

    cat("🔄 Running ORA for", analysis_name, "...\n")

    result <- run_ora_safe(
      gene_symbols,
      all_human_genes,
      collection_data,
      analysis_name,
      STATS
    )

    if (!is.null(result)) {
      ora_results[[collection_name]] <- result
    }
  }

  cat("✅ ORA analysis complete (", length(ora_results), "successful analyses)\n\n")
  ANALYSIS_SUMMARY$ora_results <- length(ora_results)
} else {
  ora_results <- list()
  ANALYSIS_SUMMARY$ora_results <- 0
  cat("⏭️  ORA analysis skipped\n\n")
}

# =============================================================================
# PHASE 6: ADDITIONAL PATHWAY DATABASES
# =============================================================================

cat("🔍 PHASE 6: Additional Pathway Databases\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

additional_results <- list()

# KEGG pathway analysis
if (WORKFLOW$run_kegg && PATHWAY_DATABASES$kegg$enabled) {
  cat("🔄 Running KEGG pathway analysis...\n")

  kegg_result <- tryCatch(
    {
      enrichKEGG(
        gene = gene_entrez_table$ENTREZID,
        organism = PATHWAY_DATABASES$kegg$organism,
        pvalueCutoff = STATS$pvalue_cutoff,
        pAdjustMethod = "BH"
      )
    },
    error = function(e) {
      log_message(paste("KEGG analysis error:", e$message), "ERROR", LOG_CONNECTION)
      NULL
    }
  )

  if (!is.null(kegg_result) && nrow(kegg_result@result) > 0) {
    additional_results$kegg <- kegg_result
    cat("✅ KEGG:", nrow(kegg_result@result), "pathways\n")
  } else {
    cat("⚠️  KEGG: no significant pathways\n")
  }
}

# Reactome pathway analysis
if (WORKFLOW$run_reactome && PATHWAY_DATABASES$reactome$enabled) {
  cat("🔄 Running Reactome pathway analysis...\n")

  reactome_result <- tryCatch(
    {
      enrichPathway(
        gene = gene_entrez_table$ENTREZID,
        pvalueCutoff = STATS$pvalue_cutoff,
        pAdjustMethod = "BH",
        readable = TRUE
      )
    },
    error = function(e) {
      log_message(paste("Reactome analysis error:", e$message), "ERROR", LOG_CONNECTION)
      NULL
    }
  )

  if (!is.null(reactome_result) && nrow(reactome_result@result) > 0) {
    additional_results$reactome <- reactome_result
    cat("✅ Reactome:", nrow(reactome_result@result), "pathways\n")
  } else {
    cat("⚠️  Reactome: no significant pathways\n")
  }
}

cat("✅ Additional database analysis complete\n\n")

# =============================================================================
# PHASE 7: VISUALIZATION
# =============================================================================

if (WORKFLOW$create_enhanced_plots) {
  cat("🎨 PHASE 7: Enhanced Visualization\n")
  cat(paste(rep("-", 40), collapse = ""), "\n")

  plot_count <- 0

  # Function to create enhanced plots for each result set
  create_result_plots <- function(result_obj, result_name, result_type = "ORA") {
    if (is.null(result_obj) || nrow(result_obj@result) == 0) {
      return(0)
    }

    plots_created <- 0
    result_data <- result_obj@result

    # Clean pathway names
    result_data$Description_clean <- clean_pathway_names(
      result_data$Description,
      PLOT_PARAMS$max_pathway_name_length
    )

    # Enhanced dot plot
    if (PLOT_TYPES$enhanced_dotplot) {
      tryCatch(
        {
          plot_data <- result_data %>%
            arrange(p.adjust) %>%
            head(PLOT_PARAMS$max_pathways_plot) %>%
            mutate(
              neg_log_p = -log10(p.adjust),
              Count = as.numeric(Count)
            )

          p_dot <- ggplot(plot_data, aes(x = Count, y = reorder(Description_clean, Count))) +
            geom_point(aes(size = Count, color = neg_log_p), alpha = 0.8) +
            scale_color_viridis_c(name = "-log10(FDR)", direction = -1) +
            scale_size_continuous(name = "Gene Count", range = PLOT_PARAMS$point_size_range) +
            labs(
              title = paste(result_name, "- Enhanced Dot Plot"),
              subtitle = paste("Top", nrow(plot_data), "pathways by gene count"),
              x = "Gene Count",
              y = "Pathway"
            ) +
            theme_minimal() +
            theme(
              plot.title = element_text(size = PLOT_PARAMS$text_size_title, face = "bold"),
              plot.subtitle = element_text(size = PLOT_PARAMS$text_size_subtitle),
              axis.text.y = element_text(size = PLOT_PARAMS$text_size_axis),
              legend.position = "right"
            )

          filename <- file.path(PATHS$plot_dir, paste0(
            tolower(gsub(" ", "_", result_name)),
            "_enhanced_dotplot.png"
          ))
          if (save_plot_safe(p_dot, filename, PLOT_PARAMS$width, PLOT_PARAMS$height, PLOT_PARAMS$dpi)) {
            plots_created <- plots_created + 1
          }
        },
        error = function(e) {
          log_message(paste("Dot plot error for", result_name, ":", e$message), "ERROR", LOG_CONNECTION)
        }
      )
    }

    # Enhanced bar plot
    if (PLOT_TYPES$enhanced_barplot) {
      tryCatch(
        {
          plot_data <- result_data %>%
            arrange(p.adjust) %>%
            head(PLOT_PARAMS$max_pathways_plot) %>%
            mutate(
              FoldEnrichment = as.numeric(FoldEnrichment),
              Count = as.numeric(Count)
            )

          p_bar <- ggplot(plot_data, aes(x = Count, y = reorder(Description_clean, Count))) +
            geom_col(aes(fill = FoldEnrichment), alpha = 0.8) +
            scale_fill_viridis_c(name = "Fold\nEnrichment", option = "plasma") +
            labs(
              title = paste(result_name, "- Enhanced Bar Plot"),
              subtitle = paste("Top", nrow(plot_data), "pathways by significance"),
              x = "Gene Count",
              y = "Pathway"
            ) +
            theme_minimal() +
            theme(
              plot.title = element_text(size = PLOT_PARAMS$text_size_title, face = "bold"),
              plot.subtitle = element_text(size = PLOT_PARAMS$text_size_subtitle),
              axis.text.y = element_text(size = PLOT_PARAMS$text_size_axis),
              legend.position = "right"
            )

          filename <- file.path(PATHS$plot_dir, paste0(
            tolower(gsub(" ", "_", result_name)),
            "_enhanced_barplot.png"
          ))
          if (save_plot_safe(p_bar, filename, PLOT_PARAMS$width, PLOT_PARAMS$height, PLOT_PARAMS$dpi)) {
            plots_created <- plots_created + 1
          }
        },
        error = function(e) {
          log_message(paste("Bar plot error for", result_name, ":", e$message), "ERROR", LOG_CONNECTION)
        }
      )
    }

    return(plots_created)
  }

  # Create plots for ORA results
  for (name in names(ora_results)) {
    analysis_name <- MSIGDB_COLLECTIONS[[name]]$name
    plots_created <- create_result_plots(ora_results[[name]], analysis_name, "ORA")
    plot_count <- plot_count + plots_created
    cat("🎨", analysis_name, ":", plots_created, "plots created\n")
  }

  # Create plots for additional database results
  for (name in names(additional_results)) {
    db_name <- PATHWAY_DATABASES[[name]]$name
    plots_created <- create_result_plots(additional_results[[name]], db_name, "Pathway")
    plot_count <- plot_count + plots_created
    cat("🎨", db_name, ":", plots_created, "plots created\n")
  }

  # Create summary visualizations
  if (PLOT_TYPES$database_summary) {
    tryCatch(
      {
        summary_data <- data.frame(
          Database = character(0),
          Pathways = numeric(0),
          Top_FDR = numeric(0)
        )

        # Collect data from all results
        all_results <- c(ora_results, additional_results)
        for (name in names(all_results)) {
          result <- all_results[[name]]
          if (!is.null(result) && nrow(result@result) > 0) {
            db_name <- if (name %in% names(MSIGDB_COLLECTIONS)) {
              MSIGDB_COLLECTIONS[[name]]$name
            } else {
              PATHWAY_DATABASES[[name]]$name
            }

            summary_data <- rbind(summary_data, data.frame(
              Database = db_name,
              Pathways = nrow(result@result),
              Top_FDR = min(result@result$p.adjust)
            ))
          }
        }

        if (nrow(summary_data) > 0) {
          p_summary <- ggplot(summary_data, aes(x = reorder(Database, Pathways), y = Pathways)) +
            geom_col(aes(fill = -log10(Top_FDR)), alpha = 0.8) +
            geom_text(aes(label = Pathways), hjust = -0.1, size = 4) +
            scale_fill_viridis_c(name = "Top Pathway\n-log10(FDR)") +
            coord_flip() +
            labs(
              title = "Pathway Analysis Database Summary",
              subtitle = "Number of enriched pathways per database",
              x = "Database",
              y = "Number of Enriched Pathways"
            ) +
            theme_minimal() +
            theme(
              plot.title = element_text(size = PLOT_PARAMS$text_size_title, face = "bold"),
              plot.subtitle = element_text(size = PLOT_PARAMS$text_size_subtitle)
            )

          filename <- file.path(PATHS$plot_dir, "database_summary.png")
          if (save_plot_safe(p_summary, filename, PLOT_PARAMS$width, 6, PLOT_PARAMS$dpi)) {
            plot_count <- plot_count + 1
          }
        }
      },
      error = function(e) {
        log_message(paste("Database summary plot error:", e$message), "ERROR", LOG_CONNECTION)
      }
    )
  }

  cat("✅ Visualization complete (", plot_count, "plots created)\n\n")
  ANALYSIS_SUMMARY$plots_created <- plot_count
} else {
  cat("⏭️  Visualization skipped\n\n")
  ANALYSIS_SUMMARY$plots_created <- 0
}

# =============================================================================
# PHASE 8: SAVE RESULTS
# =============================================================================

cat("💾 PHASE 8: Saving Results\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

files_saved <- 0

# Save ORA results
for (name in names(ora_results)) {
  filename <- file.path(PATHS$table_dir, paste0("ora_", name, "_results.csv"))
  readr::write_csv(ora_results[[name]]@result, filename)
  files_saved <- files_saved + 1
}

# Save additional database results
for (name in names(additional_results)) {
  filename <- file.path(PATHS$table_dir, paste0(name, "_results.csv"))
  readr::write_csv(additional_results[[name]]@result, filename)
  files_saved <- files_saved + 1
}

# Save input gene list
input_gene_data <- data.frame(
  Gene = names(gene_list),
  Effect_Size = gene_list,
  Entrez_ID = gene_entrez_table$ENTREZID[match(names(gene_list), gene_entrez_table$SYMBOL)]
)
readr::write_csv(input_gene_data, file.path(PATHS$table_dir, "input_gene_list.csv"))
files_saved <- files_saved + 1

cat("✅ Results saved (", files_saved, "files)\n\n")
ANALYSIS_SUMMARY$files_saved <- files_saved

# =============================================================================
# PHASE 9: GENERATE REPORTS
# =============================================================================

if (WORKFLOW$generate_reports) {
  cat("📋 PHASE 9: Generate Reports\n")
  cat(paste(rep("-", 40), collapse = ""), "\n")

  # Generate analysis summary
  analysis_summary <- generate_analysis_summary(
    c(ora_results, additional_results),
    ANALYSIS_SUMMARY$genes
  )

  # Create comprehensive report
  report_file <- file.path(PATHS$report_dir, "gsea_analysis_report.md")
  create_markdown_report(analysis_summary, report_file, include_details = TRUE)

  # Create summary statistics file
  summary_stats <- data.frame(
    Metric = c(
      "Total Genes", "Converted Genes", "Conversion Rate (%)",
      "Databases Analyzed", "Total Pathways", "Plots Created"
    ),
    Value = c(
      ANALYSIS_SUMMARY$genes$total_genes,
      ANALYSIS_SUMMARY$genes$converted_genes,
      round(ANALYSIS_SUMMARY$genes$conversion_rate, 1),
      length(c(ora_results, additional_results)),
      sum(sapply(c(ora_results, additional_results), function(x) nrow(x@result))),
      ANALYSIS_SUMMARY$plots_created
    )
  )

  readr::write_csv(summary_stats, file.path(PATHS$table_dir, "analysis_summary.csv"))

  cat("✅ Reports generated\n\n")
}

# =============================================================================
# PHASE 10: CLEANUP AND FINALIZATION
# =============================================================================

cat("🔧 PHASE 10: Cleanup and Finalization\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

# Save workspace
if (WORKFLOW$save_workspace) {
  workspace_file <- file.path(PATHS$output_dir, "gsea_master_workspace.RData")
  save.image(workspace_file)
  log_message(paste("Workspace saved to:", workspace_file), "INFO", LOG_CONNECTION)
}

# Calculate runtime
pipeline_end_time <- Sys.time()
total_runtime <- as.numeric(difftime(pipeline_end_time, pipeline_start_time, units = "mins"))

# Final logging
log_message(paste("Pipeline completed in", round(total_runtime, 2), "minutes"), "INFO", LOG_CONNECTION)
log_message(
  paste(
    "Total pathways found:",
    sum(sapply(c(ora_results, additional_results), function(x) nrow(x@result)))
  ),
  "INFO", LOG_CONNECTION
)

# Close log file
if (!is.null(LOG_CONNECTION)) {
  close(LOG_CONNECTION)
}

# Memory cleanup
cleanup_memory()

cat("✅ Cleanup complete\n\n")

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("🎉 GSEA MASTER PIPELINE COMPLETE!\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Print comprehensive summary
cat("📊 ANALYSIS SUMMARY:\n")
cat("  - Runtime:", round(total_runtime, 2), "minutes\n")
cat("  - Input genes:", ANALYSIS_SUMMARY$genes$total_genes, "\n")
cat("  - Gene conversion rate:", round(ANALYSIS_SUMMARY$genes$conversion_rate, 1), "%\n")
cat("  - ORA analyses:", ANALYSIS_SUMMARY$ora_results, "\n")
cat("  - Additional databases:", length(additional_results), "\n")

# Calculate total pathways
total_pathways <- sum(sapply(c(ora_results, additional_results), function(x) {
  if (!is.null(x)) nrow(x@result) else 0
}))
cat("  - Total enriched pathways:", total_pathways, "\n")
cat("  - Visualizations created:", ANALYSIS_SUMMARY$plots_created, "\n")
cat("  - Files saved:", ANALYSIS_SUMMARY$files_saved, "\n")

cat("\n📁 OUTPUT LOCATIONS:\n")
cat("  📊 Plots:", PATHS$plot_dir, "\n")
cat("  📋 Tables:", PATHS$table_dir, "\n")
cat("  📄 Reports:", PATHS$report_dir, "\n")
cat("  📜 Logs:", PATHS$log_dir, "\n")

# Highlight key findings
cat("\n🔑 KEY FINDINGS:\n")
all_results <- c(ora_results, additional_results)
for (name in names(all_results)) {
  result <- all_results[[name]]
  if (!is.null(result) && nrow(result@result) > 0) {
    top_pathway <- result@result[1, ]
    analysis_name <- if (name %in% names(MSIGDB_COLLECTIONS)) {
      MSIGDB_COLLECTIONS[[name]]$name
    } else {
      PATHWAY_DATABASES[[name]]$name
    }

    cat("🧠", analysis_name, "top pathway:\n")
    pathway_name <- if (nchar(top_pathway$Description) > 60) {
      paste0(substr(top_pathway$Description, 1, 57), "...")
    } else {
      top_pathway$Description
    }
    cat("   ", pathway_name, "\n")
    cat(
      "   FDR:", formatC(top_pathway$p.adjust, format = "e", digits = 2),
      ", Genes:", top_pathway$Count, "\n"
    )
  }
}

cat("\n🎯 NEXT STEPS:\n")
cat("1. Review the analysis report in outputs/reports/\n")
cat("2. Examine pathway results in outputs/tables/\n")
cat("3. Explore visualizations in outputs/plots/\n")
cat("4. Consider functional validation of top pathways\n")
cat("5. Integrate findings with literature review\n")

cat("\n🔬 Analysis completed successfully! 🧬\n")
cat("📧 For questions, refer to the project documentation.\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("🕐 Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
