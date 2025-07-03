# =============================================================================
# Complete GSEA Pipeline Execution Script
# =============================================================================
#
# Purpose: Master script to run complete pathway analysis pipeline
# Features:
#   - Automated validation and setup
#   - Comprehensive GSEA and ORA analysis
#   - Multiple database integration
#   - Robust error handling and logging
#   - Automated report generation
#
# Author: Multi-Omics OUD Study
# Date: 2024
# =============================================================================

# Start timing
start_time <- Sys.time()

cat("🧬 Complete GSEA Pipeline for LEMUR Results\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("🕐 Analysis started:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# =============================================================================
# SETUP AND VALIDATION
# =============================================================================

cat("🔧 PHASE 1: Setup and Validation\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

# Function to safely load packages
load_package <- function(pkg, source = "auto") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("📦 Installing", pkg, "...\n")
    if (source == "bioc" || pkg %in% c(
      "clusterProfiler", "org.Hs.eg.db", "msigdbr",
      "enrichplot", "DOSE", "ReactomePA", "pathview"
    )) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
  }

  library(pkg, character.only = TRUE)
  cat("✅", pkg, "loaded\n")
}

# Load all required packages
required_packages <- list(
  bioc = c(
    "clusterProfiler", "org.Hs.eg.db", "msigdbr", "enrichplot",
    "DOSE", "ReactomePA", "pathview"
  ),
  cran = c(
    "ggplot2", "dplyr", "readr", "httr", "jsonlite", "stringr",
    "VennDiagram", "pheatmap", "RColorBrewer", "here", "knitr",
    "rmarkdown", "DT", "plotly"
  )
)

cat("📦 Loading required packages...\n")
for (pkg in c(required_packages$bioc, required_packages$cran)) {
  tryCatch(
    {
      load_package(pkg)
    },
    error = function(e) {
      cat("❌ Failed to load", pkg, ":", e$message, "\n")
    }
  )
}

# Set up directories
BASE_DIR <- "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study/scripts/snrna/LEMUR Analysis"
INPUT_FILE <- file.path(BASE_DIR, "outputs/tables/top_de_genes.csv")
OUTPUT_DIR <- file.path(BASE_DIR, "LEMUR GSEA/outputs")

# Create output directories
dirs_to_create <- c(
  OUTPUT_DIR,
  file.path(OUTPUT_DIR, "plots"),
  file.path(OUTPUT_DIR, "tables"),
  file.path(OUTPUT_DIR, "reports"),
  file.path(OUTPUT_DIR, "logs")
)

for (dir_path in dirs_to_create) {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

# Setup logging
log_file <- file.path(OUTPUT_DIR, "logs", paste0("gsea_pipeline_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con <- file(log_file, "w")

# Function to log messages
log_message <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- paste0("[", timestamp, "] [", level, "] ", msg)
  cat(log_entry, "\n")
  writeLines(log_entry, log_con)
  flush(log_con)
}

log_message("GSEA Pipeline Started")
log_message(paste("Base directory:", BASE_DIR))
log_message(paste("Output directory:", OUTPUT_DIR))

# =============================================================================
# DATA VALIDATION AND LOADING
# =============================================================================

cat("\n📊 PHASE 2: Data Validation and Loading\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

# Validate input file
if (!file.exists(INPUT_FILE)) {
  log_message(paste("Input file not found:", INPUT_FILE), "ERROR")
  stop("❌ Input file not found. Please run LEMUR analysis first.")
}

log_message("Loading differential expression results")
de_results <- tryCatch(
  {
    read_csv(INPUT_FILE, show_col_types = FALSE)
  },
  error = function(e) {
    log_message(paste("Error loading input file:", e$message), "ERROR")
    stop("❌ Failed to load input data")
  }
)

# Validate data structure
required_columns <- c("gene", "effect_size", "pval", "adj_pval")
missing_cols <- setdiff(required_columns, colnames(de_results))

if (length(missing_cols) > 0) {
  log_message(paste("Missing required columns:", paste(missing_cols, collapse = ", ")), "ERROR")
  stop("❌ Invalid data structure")
}

# Data summary
n_genes <- nrow(de_results)
effect_range <- range(de_results$effect_size, na.rm = TRUE)
n_significant <- sum(de_results$adj_pval < 0.05, na.rm = TRUE)

log_message(paste("Loaded", n_genes, "genes"))
log_message(paste("Effect size range:", round(effect_range[1], 3), "to", round(effect_range[2], 3)))
log_message(paste("Significant genes (FDR < 0.05):", n_significant))

cat("✅ Data validation complete\n")
cat("📊 Dataset summary:\n")
cat("  - Total genes:", n_genes, "\n")
cat("  - Effect size range:", round(effect_range[1], 3), "to", round(effect_range[2], 3), "\n")
cat("  - Significant genes:", n_significant, "\n")

# =============================================================================
# GENE LIST PREPARATION
# =============================================================================

cat("\n🧬 PHASE 3: Gene List Preparation\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

# Create ranked gene list
gene_list <- de_results$effect_size
names(gene_list) <- de_results$gene
gene_list <- sort(gene_list, decreasing = TRUE)

# Remove any NA values and ensure numeric
gene_list <- gene_list[!is.na(gene_list)]
gene_list <- as.numeric(gene_list)
names(gene_list) <- de_results$gene[!is.na(de_results$effect_size)]

# Extract gene symbols
gene_symbols <- de_results$gene

# Convert to Entrez IDs
log_message("Converting gene symbols to Entrez IDs")
gene_entrez <- tryCatch(
  {
    bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  },
  error = function(e) {
    log_message(paste("Gene ID conversion error:", e$message), "ERROR")
    stop("❌ Failed to convert gene IDs")
  }
)

# Create Entrez ID ranked list
gene_list_entrez <- gene_list[gene_entrez$SYMBOL]
names(gene_list_entrez) <- gene_entrez$ENTREZID

# Remove any NA values and ensure proper sorting
gene_list_entrez <- gene_list_entrez[!is.na(gene_list_entrez)]
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)

conversion_rate <- nrow(gene_entrez) / length(gene_symbols) * 100
log_message(paste(
  "Gene ID conversion:", nrow(gene_entrez), "/", length(gene_symbols),
  "(", round(conversion_rate, 1), "%)"
))

cat("✅ Gene list preparation complete\n")
cat(
  "🔄 Gene ID conversion:", nrow(gene_entrez), "/", length(gene_symbols),
  "(", round(conversion_rate, 1), "%)\n"
)

# =============================================================================
# GENE SET DATABASE PREPARATION
# =============================================================================

cat("\n🗃️ PHASE 4: Gene Set Database Preparation\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

# Load MSigDB collections
log_message("Loading MSigDB gene sets")

msigdb_collections <- list()
collection_info <- list(
  hallmark = list(category = "H", subcategory = NULL, name = "Hallmark"),
  canonical = list(category = "C2", subcategory = "CP", name = "Canonical Pathways"),
  go_bp = list(category = "C5", subcategory = "GO:BP", name = "GO Biological Process"),
  go_cc = list(category = "C5", subcategory = "GO:CC", name = "GO Cellular Component"),
  go_mf = list(category = "C5", subcategory = "GO:MF", name = "GO Molecular Function")
)

for (collection_name in names(collection_info)) {
  info <- collection_info[[collection_name]]
  log_message(paste("Loading", info$name, "gene sets"))

  msigdb_collections[[collection_name]] <- tryCatch(
    {
      if (is.null(info$subcategory)) {
        msigdbr(species = "Homo sapiens", category = info$category)
      } else {
        msigdbr(species = "Homo sapiens", category = info$category, subcategory = info$subcategory)
      }
    },
    error = function(e) {
      log_message(paste("Error loading", info$name, ":", e$message), "WARNING")
      NULL
    }
  )
}

# Filter for neuroscience/addiction relevant pathways
neuro_keywords <- c(
  "synaptic", "synapse", "neuron", "neural", "dopamine", "serotonin",
  "GABA", "glutamate", "addiction", "reward", "plasticity", "learning",
  "memory", "behavior", "channel", "receptor", "neurotransmitter",
  "signal", "calcium", "potassium", "sodium", "vesicle", "axon",
  "dendrite", "spine", "excitatory", "inhibitory"
)

filter_neuro_pathways <- function(msigdb_data) {
  if (is.null(msigdb_data) || nrow(msigdb_data) == 0) {
    return(NULL)
  }
  pattern <- paste(neuro_keywords, collapse = "|")
  msigdb_data[grepl(pattern, msigdb_data$gs_name, ignore.case = TRUE), ]
}

# Apply neuroscience filter
msigdb_neuro <- lapply(msigdb_collections, filter_neuro_pathways)
msigdb_neuro <- msigdb_neuro[!sapply(msigdb_neuro, is.null)]

# Log gene set statistics
for (collection_name in names(msigdb_neuro)) {
  n_pathways <- length(unique(msigdb_neuro[[collection_name]]$gs_name))
  log_message(paste(collection_info[[collection_name]]$name, "neuro pathways:", n_pathways))
  cat("📊", collection_info[[collection_name]]$name, ":", n_pathways, "neuro-relevant pathways\n")
}

# =============================================================================
# GSEA ANALYSIS
# =============================================================================

cat("\n🔬 PHASE 5: Gene Set Enrichment Analysis (GSEA)\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

# GSEA parameters
PARAMS <- list(
  pvalue_cutoff = 0.05,
  qvalue_cutoff = 0.2,
  min_gs_size = 10,
  max_gs_size = 500,
  eps = 1e-10
)

log_message("Starting GSEA analyses")

# Function to run GSEA safely
run_gsea_safe <- function(gene_list, gene_sets, analysis_name) {
  log_message(paste("Running GSEA for", analysis_name))

  if (is.null(gene_sets) || nrow(gene_sets) == 0) {
    log_message(paste("No gene sets available for", analysis_name), "WARNING")
    return(NULL)
  }

  tryCatch(
    {
      result <- GSEA(
        geneList = gene_list,
        TERM2GENE = gene_sets[, c("gs_name", "entrez_gene")],
        pvalueCutoff = PARAMS$pvalue_cutoff,
        pAdjustMethod = "BH",
        minGSSize = PARAMS$min_gs_size,
        maxGSSize = PARAMS$max_gs_size,
        eps = PARAMS$eps,
        verbose = FALSE
      )

      if (!is.null(result) && nrow(result@result) > 0) {
        n_enriched <- nrow(result@result)
        log_message(paste(analysis_name, "GSEA: found", n_enriched, "enriched pathways"))
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

# Run GSEA analyses
gsea_results <- list()
for (collection_name in names(msigdb_neuro)) {
  gsea_results[[collection_name]] <- run_gsea_safe(
    gene_list_entrez,
    msigdb_neuro[[collection_name]],
    collection_info[[collection_name]]$name
  )
}

# Remove NULL results
gsea_results <- gsea_results[!sapply(gsea_results, is.null)]

cat("✅ GSEA analysis complete\n")
cat("📊 Results:", length(gsea_results), "successful analyses\n")

# =============================================================================
# OVER-REPRESENTATION ANALYSIS (ORA)
# =============================================================================

cat("\n🎯 PHASE 6: Over-Representation Analysis (ORA)\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

log_message("Starting ORA analyses")

# Function to run ORA safely
run_ora_safe <- function(genes, universe_genes, gene_sets, analysis_name) {
  log_message(paste("Running ORA for", analysis_name))

  if (is.null(gene_sets) || nrow(gene_sets) == 0) {
    log_message(paste("No gene sets available for", analysis_name), "WARNING")
    return(NULL)
  }

  tryCatch(
    {
      result <- enricher(
        gene = genes,
        universe = universe_genes,
        TERM2GENE = gene_sets[, c("gs_name", "gene_symbol")],
        pvalueCutoff = PARAMS$pvalue_cutoff,
        pAdjustMethod = "BH",
        minGSSize = PARAMS$min_gs_size,
        maxGSSize = PARAMS$max_gs_size
      )

      if (!is.null(result) && nrow(result@result) > 0) {
        n_enriched <- nrow(result@result)
        log_message(paste(analysis_name, "ORA: found", n_enriched, "enriched pathways"))
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

# Get universe genes
all_human_genes <- keys(org.Hs.eg.db, keytype = "SYMBOL")

# Run ORA analyses
ora_results <- list()
for (collection_name in names(msigdb_neuro)) {
  ora_results[[collection_name]] <- run_ora_safe(
    gene_symbols,
    all_human_genes,
    msigdb_neuro[[collection_name]],
    collection_info[[collection_name]]$name
  )
}

# Remove NULL results
ora_results <- ora_results[!sapply(ora_results, is.null)]

cat("✅ ORA analysis complete\n")
cat("📊 Results:", length(ora_results), "successful analyses\n")

# =============================================================================
# ADDITIONAL PATHWAY DATABASES
# =============================================================================

cat("\n🔍 PHASE 7: Additional Pathway Databases\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

# KEGG pathway analysis
log_message("Running KEGG pathway analysis")
kegg_result <- tryCatch(
  {
    enrichKEGG(
      gene = gene_entrez$ENTREZID,
      organism = "hsa",
      pvalueCutoff = PARAMS$pvalue_cutoff,
      pAdjustMethod = "BH"
    )
  },
  error = function(e) {
    log_message(paste("KEGG analysis error:", e$message), "ERROR")
    NULL
  }
)

# Reactome pathway analysis
log_message("Running Reactome pathway analysis")
reactome_result <- tryCatch(
  {
    enrichPathway(
      gene = gene_entrez$ENTREZID,
      pvalueCutoff = PARAMS$pvalue_cutoff,
      pAdjustMethod = "BH",
      readable = TRUE
    )
  },
  error = function(e) {
    log_message(paste("Reactome analysis error:", e$message), "ERROR")
    NULL
  }
)

# Disease Ontology analysis
log_message("Running Disease Ontology analysis")
do_result <- tryCatch(
  {
    enrichDO(
      gene = gene_entrez$ENTREZID,
      ont = "DO",
      pvalueCutoff = PARAMS$pvalue_cutoff,
      pAdjustMethod = "BH",
      readable = TRUE
    )
  },
  error = function(e) {
    log_message(paste("Disease Ontology analysis error:", e$message), "ERROR")
    NULL
  }
)

# Report additional analyses
additional_results <- list(kegg = kegg_result, reactome = reactome_result, disease_ontology = do_result)
for (name in names(additional_results)) {
  result <- additional_results[[name]]
  if (!is.null(result) && nrow(result@result) > 0) {
    cat("✅", toupper(name), ":", nrow(result@result), "pathways\n")
  } else {
    cat("⚠️ ", toupper(name), ": no significant pathways\n")
  }
}

# =============================================================================
# VISUALIZATION
# =============================================================================

cat("\n📊 PHASE 8: Visualization Generation\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

log_message("Generating visualizations")

# Safe plotting function
safe_plot <- function(plot_func, filename, width = 12, height = 8) {
  tryCatch(
    {
      p <- plot_func()
      if (!is.null(p)) {
        ggsave(filename, plot = p, width = width, height = height, dpi = 300)
        log_message(paste("Saved plot:", basename(filename)))
        return(TRUE)
      }
      return(FALSE)
    },
    error = function(e) {
      log_message(paste("Plot error for", basename(filename), ":", e$message), "ERROR")
      return(FALSE)
    }
  )
}

plot_count <- 0

# GSEA plots
for (analysis_name in names(gsea_results)) {
  result <- gsea_results[[analysis_name]]

  # Dot plot
  if (safe_plot(
    function() {
      dotplot(result, showCategory = 20) +
        ggtitle(paste("GSEA", collection_info[[analysis_name]]$name, "- Dot Plot"))
    },
    file.path(OUTPUT_DIR, "plots", paste0("gsea_", analysis_name, "_dotplot.png"))
  )) {
    plot_count <- plot_count + 1
  }

  # Ridge plot
  if (safe_plot(
    function() {
      ridgeplot(result, showCategory = 15) +
        ggtitle(paste("GSEA", collection_info[[analysis_name]]$name, "- Ridge Plot"))
    },
    file.path(OUTPUT_DIR, "plots", paste0("gsea_", analysis_name, "_ridgeplot.png"))
  )) {
    plot_count <- plot_count + 1
  }
}

# ORA plots
for (analysis_name in names(ora_results)) {
  result <- ora_results[[analysis_name]]

  # Bar plot
  if (safe_plot(
    function() {
      barplot(result, showCategory = 20) +
        ggtitle(paste("ORA", collection_info[[analysis_name]]$name, "- Bar Plot"))
    },
    file.path(OUTPUT_DIR, "plots", paste0("ora_", analysis_name, "_barplot.png"))
  )) {
    plot_count <- plot_count + 1
  }

  # Dot plot
  if (safe_plot(
    function() {
      dotplot(result, showCategory = 20) +
        ggtitle(paste("ORA", collection_info[[analysis_name]]$name, "- Dot Plot"))
    },
    file.path(OUTPUT_DIR, "plots", paste0("ora_", analysis_name, "_dotplot.png"))
  )) {
    plot_count <- plot_count + 1
  }
}

# Additional database plots
if (!is.null(kegg_result) && nrow(kegg_result@result) > 0) {
  if (safe_plot(
    function() dotplot(kegg_result, showCategory = 20) + ggtitle("KEGG Pathways"),
    file.path(OUTPUT_DIR, "plots", "kegg_dotplot.png")
  )) {
    plot_count <- plot_count + 1
  }
}

if (!is.null(reactome_result) && nrow(reactome_result@result) > 0) {
  if (safe_plot(
    function() dotplot(reactome_result, showCategory = 20) + ggtitle("Reactome Pathways"),
    file.path(OUTPUT_DIR, "plots", "reactome_dotplot.png")
  )) {
    plot_count <- plot_count + 1
  }
}

cat("✅ Visualization complete\n")
cat("📊 Generated:", plot_count, "plots\n")

# =============================================================================
# SAVE RESULTS
# =============================================================================

cat("\n💾 PHASE 9: Saving Results\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

log_message("Saving analysis results")

file_count <- 0

# Save GSEA results
for (name in names(gsea_results)) {
  filename <- file.path(OUTPUT_DIR, "tables", paste0("gsea_", name, "_results.csv"))
  write_csv(gsea_results[[name]]@result, filename)
  file_count <- file_count + 1
}

# Save ORA results
for (name in names(ora_results)) {
  filename <- file.path(OUTPUT_DIR, "tables", paste0("ora_", name, "_results.csv"))
  write_csv(ora_results[[name]]@result, filename)
  file_count <- file_count + 1
}

# Save additional results
if (!is.null(kegg_result) && nrow(kegg_result@result) > 0) {
  write_csv(kegg_result@result, file.path(OUTPUT_DIR, "tables", "kegg_results.csv"))
  file_count <- file_count + 1
}

if (!is.null(reactome_result) && nrow(reactome_result@result) > 0) {
  write_csv(reactome_result@result, file.path(OUTPUT_DIR, "tables", "reactome_results.csv"))
  file_count <- file_count + 1
}

if (!is.null(do_result) && nrow(do_result@result) > 0) {
  write_csv(do_result@result, file.path(OUTPUT_DIR, "tables", "disease_ontology_results.csv"))
  file_count <- file_count + 1
}

# Save input gene list
write_csv(
  data.frame(Gene = names(gene_list), Effect_Size = gene_list),
  file.path(OUTPUT_DIR, "tables", "input_gene_list.csv")
)
file_count <- file_count + 1

cat("✅ Results saved\n")
cat("📁 Created:", file_count, "result files\n")

# =============================================================================
# GENERATE COMPREHENSIVE REPORT
# =============================================================================

cat("\n📋 PHASE 10: Generating Comprehensive Report\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

log_message("Generating comprehensive summary report")

# Calculate end time and duration
end_time <- Sys.time()
duration <- as.numeric(difftime(end_time, start_time, units = "mins"))

# Create detailed summary report
create_detailed_report <- function() {
  report_file <- file.path(OUTPUT_DIR, "reports", "comprehensive_GSEA_report.md")

  # Collect all significant pathways
  all_significant_pathways <- data.frame()

  # GSEA pathways
  for (name in names(gsea_results)) {
    if (!is.null(gsea_results[[name]])) {
      df <- gsea_results[[name]]@result
      df$Analysis_Type <- "GSEA"
      df$Database <- collection_info[[name]]$name
      df$Analysis_ID <- name
      all_significant_pathways <- rbind(
        all_significant_pathways,
        df[, c(
          "ID", "Description", "Analysis_Type", "Database",
          "Analysis_ID", "NES", "pvalue", "p.adjust"
        )]
      )
    }
  }

  # ORA pathways
  for (name in names(ora_results)) {
    if (!is.null(ora_results[[name]])) {
      df <- ora_results[[name]]@result
      df$Analysis_Type <- "ORA"
      df$Database <- collection_info[[name]]$name
      df$Analysis_ID <- name
      df$NES <- NA # ORA doesn't have NES
      all_significant_pathways <- rbind(
        all_significant_pathways,
        df[, c(
          "ID", "Description", "Analysis_Type", "Database",
          "Analysis_ID", "NES", "pvalue", "p.adjust"
        )]
      )
    }
  }

  # Generate report content
  report_content <- paste0(
    "# Comprehensive Gene Set Enrichment Analysis Report\n\n",
    "**Analysis Date:** ", format(Sys.Date(), "%B %d, %Y"), "\n",
    "**Analysis Duration:** ", round(duration, 2), " minutes\n",
    "**Software:** R ", R.version.string, "\n\n",
    "## Executive Summary\n\n",
    "This comprehensive pathway analysis identified functional enrichments among ",
    n_genes, " genes from LEMUR differential expression analysis, focusing on addiction-related ",
    "pathways, neurotransmitter signaling, and synaptic function.\n\n",
    "### Key Statistics\n",
    "- **Input Genes:** ", n_genes, " (", n_significant, " significant at FDR < 0.05)\n",
    "- **Effect Size Range:** ", round(effect_range[1], 3), " to ", round(effect_range[2], 3), "\n",
    "- **Gene ID Conversion:** ", nrow(gene_entrez), "/", length(gene_symbols),
    " (", round(conversion_rate, 1), "%)\n",
    "- **Total Enriched Pathways:** ", nrow(all_significant_pathways), "\n",
    "- **Analysis Methods:** GSEA, Over-Representation Analysis\n",
    "- **Databases:** MSigDB, KEGG, Reactome, Disease Ontology\n\n",
    "## Analysis Methods\n\n",
    "### Gene Set Enrichment Analysis (GSEA)\n",
    "- **Approach:** Rank-based enrichment using effect sizes\n",
    "- **Advantage:** Considers all genes, not just significant ones\n",
    "- **Metric:** Normalized Enrichment Score (NES)\n\n",
    "### Over-Representation Analysis (ORA)\n",
    "- **Approach:** Statistical over-representation of significant genes\n",
    "- **Advantage:** Direct interpretation of significant gene overlap\n",
    "- **Metric:** Fold enrichment and p-values\n\n",
    "## Results Summary\n\n"
  )

  # Add analysis-specific results
  if (length(gsea_results) > 0) {
    report_content <- paste0(report_content, "### GSEA Results\n\n")
    for (name in names(gsea_results)) {
      n_pathways <- nrow(gsea_results[[name]]@result)
      report_content <- paste0(
        report_content, "- **", collection_info[[name]]$name,
        ":** ", n_pathways, " enriched pathways\n"
      )
    }
    report_content <- paste0(report_content, "\n")
  }

  if (length(ora_results) > 0) {
    report_content <- paste0(report_content, "### Over-Representation Analysis Results\n\n")
    for (name in names(ora_results)) {
      n_pathways <- nrow(ora_results[[name]]@result)
      report_content <- paste0(
        report_content, "- **", collection_info[[name]]$name,
        ":** ", n_pathways, " enriched pathways\n"
      )
    }
    report_content <- paste0(report_content, "\n")
  }

  # Top pathways section
  if (nrow(all_significant_pathways) > 0) {
    # Get top 10 most significant pathways overall
    top_pathways <- all_significant_pathways[order(all_significant_pathways$p.adjust), ][1:min(10, nrow(all_significant_pathways)), ]

    report_content <- paste0(report_content, "### Top 10 Most Significant Pathways\n\n")
    report_content <- paste0(report_content, "| Rank | Pathway | Database | Analysis | P-value | FDR |\n")
    report_content <- paste0(report_content, "|------|---------|----------|-----------|---------|-----|\n")

    for (i in 1:nrow(top_pathways)) {
      report_content <- paste0(
        report_content, "| ", i, " | ",
        substr(top_pathways$Description[i], 1, 50),
        ifelse(nchar(top_pathways$Description[i]) > 50, "...", ""),
        " | ", top_pathways$Database[i],
        " | ", top_pathways$Analysis_Type[i],
        " | ", formatC(top_pathways$pvalue[i], format = "e", digits = 2),
        " | ", formatC(top_pathways$p.adjust[i], format = "e", digits = 2), " |\n"
      )
    }
    report_content <- paste0(report_content, "\n")
  }

  # Add pathway categories analysis
  if (nrow(all_significant_pathways) > 0) {
    # Identify addiction-related pathways
    addiction_keywords <- c("addiction", "reward", "dopamine", "substance", "alcohol", "drug")
    addiction_pattern <- paste(addiction_keywords, collapse = "|")
    addiction_pathways <- all_significant_pathways[grepl(addiction_pattern,
      all_significant_pathways$Description,
      ignore.case = TRUE
    ), ]

    # Identify synaptic pathways
    synaptic_keywords <- c("synaptic", "synapse", "neurotransmitter", "vesicle", "axon", "dendrite")
    synaptic_pattern <- paste(synaptic_keywords, collapse = "|")
    synaptic_pathways <- all_significant_pathways[grepl(synaptic_pattern,
      all_significant_pathways$Description,
      ignore.case = TRUE
    ), ]

    report_content <- paste0(report_content, "### Pathway Categories\n\n")
    report_content <- paste0(report_content, "- **Addiction-Related Pathways:** ", nrow(addiction_pathways), "\n")
    report_content <- paste0(report_content, "- **Synaptic Function Pathways:** ", nrow(synaptic_pathways), "\n")
    report_content <- paste0(
      report_content, "- **Other Neuro Pathways:** ",
      nrow(all_significant_pathways) - nrow(addiction_pathways) - nrow(synaptic_pathways), "\n\n"
    )
  }

  # Add file inventory
  report_content <- paste0(report_content, "## Output Files\n\n")
  report_content <- paste0(report_content, "### Analysis Results\n")
  for (name in names(gsea_results)) {
    report_content <- paste0(report_content, "- `tables/gsea_", name, "_results.csv`\n")
  }
  for (name in names(ora_results)) {
    report_content <- paste0(report_content, "- `tables/ora_", name, "_results.csv`\n")
  }
  report_content <- paste0(report_content, "\n### Visualizations\n")
  report_content <- paste0(report_content, "- Generated ", plot_count, " visualization files in `plots/` directory\n")
  report_content <- paste0(report_content, "- Includes dot plots, ridge plots, bar plots, and network plots\n\n")

  # Add methodology details
  report_content <- paste0(report_content, "## Technical Details\n\n")
  report_content <- paste0(report_content, "### Analysis Parameters\n")
  report_content <- paste0(report_content, "- **P-value cutoff:** ", PARAMS$pvalue_cutoff, "\n")
  report_content <- paste0(report_content, "- **FDR cutoff:** ", PARAMS$qvalue_cutoff, "\n")
  report_content <- paste0(report_content, "- **Min gene set size:** ", PARAMS$min_gs_size, "\n")
  report_content <- paste0(report_content, "- **Max gene set size:** ", PARAMS$max_gs_size, "\n\n")

  report_content <- paste0(report_content, "### Database Versions\n")
  report_content <- paste0(report_content, "- **MSigDB:** Filtered for neuroscience relevance\n")
  report_content <- paste0(report_content, "- **KEGG:** Human pathways\n")
  report_content <- paste0(report_content, "- **Reactome:** Human pathways\n")
  report_content <- paste0(report_content, "- **Gene Ontology:** Biological Process, Cellular Component, Molecular Function\n\n")

  # Add interpretation guidelines
  report_content <- paste0(report_content, "## Interpretation Guidelines\n\n")
  report_content <- paste0(report_content, "### GSEA Metrics\n")
  report_content <- paste0(report_content, "- **NES (Normalized Enrichment Score):** Magnitude and direction of enrichment\n")
  report_content <- paste0(report_content, "  - Positive NES: Pathway upregulated in condition\n")
  report_content <- paste0(report_content, "  - Negative NES: Pathway downregulated in condition\n")
  report_content <- paste0(report_content, "  - |NES| > 1.0: Strong enrichment\n\n")

  report_content <- paste0(report_content, "### ORA Metrics\n")
  report_content <- paste0(report_content, "- **Fold Enrichment:** Ratio of observed vs expected genes\n")
  report_content <- paste0(report_content, "- **Gene Ratio:** Proportion of input genes in pathway\n")
  report_content <- paste0(report_content, "- **Background Ratio:** Proportion of all genes in pathway\n\n")

  # Footer
  report_content <- paste0(report_content, "---\n")
  report_content <- paste0(report_content, "*Report generated by Complete GSEA Pipeline*\n")
  report_content <- paste0(report_content, "*Multi-Omics OUD Study*\n")
  report_content <- paste0(report_content, "*Generated on: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "*\n")

  # Write report
  writeLines(report_content, report_file)
  log_message(paste("Comprehensive report saved to:", report_file))

  return(report_file)
}

report_file <- create_detailed_report()

# Save workspace
workspace_file <- file.path(OUTPUT_DIR, "complete_gsea_workspace.RData")
save.image(workspace_file)
log_message(paste("Workspace saved to:", workspace_file))

cat("✅ Report generation complete\n")
cat("📋 Comprehensive report:", basename(report_file), "\n")

# =============================================================================
# FINAL SUMMARY AND CLEANUP
# =============================================================================

cat("\n🎉 PIPELINE COMPLETE!\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Close log file
close(log_con)

# Final summary
cat("📊 ANALYSIS SUMMARY:\n")
cat("  - Duration:", round(duration, 2), "minutes\n")
cat("  - Input genes:", n_genes, "\n")
cat("  - Gene ID conversion rate:", round(conversion_rate, 1), "%\n")
cat("  - GSEA analyses:", length(gsea_results), "\n")
cat("  - ORA analyses:", length(ora_results), "\n")

# Calculate total pathways safely
total_gsea <- if (length(gsea_results) > 0) {
  sum(sapply(gsea_results, function(x) if (!is.null(x)) nrow(x@result) else 0))
} else {
  0
}

total_ora <- if (length(ora_results) > 0) {
  sum(sapply(ora_results, function(x) if (!is.null(x)) nrow(x@result) else 0))
} else {
  0
}

cat("  - Total enriched pathways:", total_gsea + total_ora, "\n")
cat("  - Visualizations created:", plot_count, "\n")
cat("  - Result files saved:", file_count, "\n")

cat("\n📁 OUTPUT STRUCTURE:\n")
cat("  📂 outputs/\n")
cat("    📊 plots/          - ", plot_count, " visualization files\n")
cat("    📋 tables/         - ", file_count, " result CSV files\n")
cat("    📄 reports/        - Comprehensive analysis report\n")
cat("    📜 logs/           - Execution log file\n")

cat("\n🔑 KEY FINDINGS:\n")
if (length(gsea_results) > 0) {
  for (name in names(gsea_results)) {
    if (nrow(gsea_results[[name]]@result) > 0) {
      top_pathway <- gsea_results[[name]]@result[1, ]
      cat("🧠", collection_info[[name]]$name, "top pathway:\n")
      cat(
        "   ", substr(top_pathway$Description, 1, 60),
        ifelse(nchar(top_pathway$Description) > 60, "...", ""), "\n"
      )
      cat(
        "   NES:", round(top_pathway$NES, 3), ", FDR:",
        formatC(top_pathway$p.adjust, format = "e", digits = 2), "\n"
      )
    }
  }
}

cat("\n🎯 NEXT STEPS:\n")
cat("1. Review the comprehensive report in outputs/reports/\n")
cat("2. Examine individual pathway results in outputs/tables/\n")
cat("3. Explore visualizations in outputs/plots/\n")
cat("4. Consider functional validation of top pathways\n")
cat("5. Integrate findings with literature and clinical relevance\n")

cat("\n🔬 Happy pathway analysis! 🧬\n")
cat("📧 For questions, refer to the project documentation.\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Display completion time
cat("🕐 Analysis completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
