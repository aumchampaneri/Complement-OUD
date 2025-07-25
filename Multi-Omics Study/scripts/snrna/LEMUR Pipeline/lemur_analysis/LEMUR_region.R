#!/usr/bin/env Rscript

# =============================================================================
# Single-nucleus RNA-seq Analysis Pipeline using LEMUR
# =============================================================================
# Title: Single nuclei transcriptomics in human and non-human primate striatum
#        in opioid use disorder
# GEO Accession: GSE225158
# Author: Generated Analysis Pipeline
# Date: July 2025
# Description: Comprehensive snRNA-seq analysis pipeline using LEMUR package
#              for differential expression analysis with interaction effects
# =============================================================================

# Set up environment and reproducibility
set.seed(42)
options(stringsAsFactors = FALSE)

# Create timestamp for output files
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
cat("\n========== LEMUR snRNA-seq Pipeline ==========\n")
cat("Start time:", timestamp, "\n")

# =============================================================================
# 1. SETUP AND DEPENDENCIES
# =============================================================================

# Required packages
required_packages <- c(
  # Core single-cell packages
  "SingleCellExperiment", "scater", "scran", "zellkonverter",
  # LEMUR and dependencies
  "lemur", "glmGamPoi", "harmony",
  # Visualization
  "ggplot2", "cowplot", "ComplexHeatmap", "pheatmap",
  # Data manipulation
  "dplyr", "tidyr", "tibble", "readr",
  # Gene set enrichment
  "clusterProfiler", "org.Hs.eg.db", "org.Mmu.eg.db", "msigdbr",
  "ReactomePA", "DOSE", "enrichplot",
  # Utilities
  "BiocParallel", "BiocSingular", "Matrix"
)

# Install missing packages
missing_packages <- required_packages[!sapply(required_packages, require, character.only = TRUE)]
if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  if (!require("BiocManager", quietly = TRUE)) suppressMessages(install.packages("BiocManager"))
  suppressMessages(suppressPackageStartupMessages(BiocManager::install(missing_packages, update = FALSE, ask = FALSE)))
}

# Maximum suppression of package loading output
sink("/dev/null")
invisible(lapply(required_packages, function(pkg) suppressMessages(suppressPackageStartupMessages(library(pkg, character.only = TRUE)))))
sink()

# Set parallel processing using BiocParallel::MulticoreParam
library(BiocParallel)
num_workers <- max(1, parallel::detectCores() - 1)
BPPARAM <- BiocParallel::MulticoreParam(workers = num_workers)
BiocParallel::register(BPPARAM)

# =============================================================================
# 2. CONFIGURATION AND PATHS
# =============================================================================

# Input configuration
input_file <- "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study/data/raw/snrna/GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5ad"
# Set output directory to outputs/lemur_analysis_TIMESTAMP/
output_base <- "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study/scripts/snrna/LEMUR Pipeline/outputs"
output_dir <- file.path(output_base, paste0("lemur_analysis_", timestamp))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cat("Output directory:", output_dir, "\n")

# Create subdirectories for organization
dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "objects"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "gsea"), recursive = TRUE, showWarnings = FALSE)

# Initialize log file
log_file <- file.path(output_dir, paste0("analysis_log_", timestamp, ".txt"))
cat("Log file:", log_file, "\n")
log_conn <- file(log_file, "w")
writeLines(paste("Analysis started at:", Sys.time()), log_conn)
writeLines(paste("Input file:", input_file), log_conn)
writeLines(paste("Output directory:", output_dir), log_conn)
close(log_conn)

# Function to log messages
log_message <- function(message) {
  cat("[", format(Sys.time(), "%H:%M:%S"), "] ", message, "\n")
  cat(paste(Sys.time(), ":", message), "\n", file = log_file, append = TRUE)
}

# Modularized: Factor level check and logging
check_and_log_factors <- function(sce, stage = "post-filtering") {
  for (col in c("condition", "sex", "brain_region", "donor_id")) {
    if (col %in% colnames(colData(sce))) {
      colData(sce)[[col]] <- droplevels(as.factor(colData(sce)[[col]]))
      tbl <- table(colData(sce)[[col]])
      log_message(paste(stage, "counts for", col, ":", paste(names(tbl), tbl, sep = "=", collapse = ", ")))
      if (length(tbl) < 2) {
        stop(paste("Grouping variable", col, "has only", length(tbl), "level(s) after", stage, ". LEMUR requires at least 2."))
      }
      if (any(tbl == 0)) {
        stop(paste("Grouping variable", col, "has a group with zero cells after", stage, ". Please adjust filtering."))
      }
    }
  }
  return(sce)
}

# Checkpointing utility
save_checkpoint <- function(obj, name, timestamp, output_dir) {
  # Always use saveHDF5SummarizedExperiment for SingleCellExperiment objects
  if (inherits(obj, "SingleCellExperiment")) {
    saveHDF5SummarizedExperiment(obj, file.path(output_dir, "checkpoints", paste0(name, "_", timestamp)))
    log_message(paste("Checkpoint saved (HDF5):", name))
  } else {
    saveRDS(obj, file.path(output_dir, "checkpoints", paste0(name, "_", timestamp, ".rds")))
    log_message(paste("Checkpoint saved:", name))
  }
}

log_message("Starting LEMUR snRNA-seq analysis pipeline")

# =============================================================================
# 3. DATA LOADING AND INITIAL INSPECTION
# =============================================================================

log_message("STEP 1: Data loading started.")
if (!file.exists(input_file)) stop("Input file not found: ", input_file)
sce <- zellkonverter::readH5AD(input_file, use_hdf5 = TRUE)
log_message(paste("Loaded SCE object with", ncol(sce), "cells and", nrow(sce), "genes"))

# ---- Assay name correction ----
if (!"counts" %in% assayNames(sce)) {
  if ("matrix" %in% assayNames(sce)) {
    assay(sce, "counts") <- assay(sce, "matrix")
  } else if ("X" %in% assayNames(sce)) {
    assay(sce, "counts") <- assay(sce, "X")
  }
}

# ---- Metadata cleaning ----
meta_data <- colData(sce)
colnames(meta_data)[colnames(meta_data) == "ID"] <- "donor_id"
colnames(meta_data)[colnames(meta_data) == "Dx_OUD"] <- "condition"
colnames(meta_data)[colnames(meta_data) == "Region"] <- "brain_region"
colnames(meta_data)[colnames(meta_data) == "Sex"] <- "sex"
colData(sce) <- meta_data

# Ensure grouping variables are factors with at least two levels
for (col in c("condition", "sex", "brain_region", "donor_id")) {
  if (col %in% colnames(colData(sce))) {
    colData(sce)[[col]] <- as.factor(colData(sce)[[col]])
    log_message(paste("Levels for", col, ":", paste(levels(colData(sce)[[col]]), collapse = ", ")))
    if (length(levels(colData(sce)[[col]])) < 2) {
      stop(paste("Variable", col, "must have at least two levels for contrasts."))
    }
  }
}

# ---- Add species column for all samples ----
if (!"species" %in% colnames(colData(sce))) {
  colData(sce)$species <- "human"
}

# ---- Check for required columns ----
required_cols <- c("donor_id", "condition")
missing_cols <- setdiff(required_cols, colnames(colData(sce)))
if (length(missing_cols) > 0) {
  stop("❌ Missing required columns: ", paste(missing_cols, collapse = ", "))
}

library(HDF5Array)
saveHDF5SummarizedExperiment(sce, file.path(output_dir, "objects", paste0("sce_initial_", timestamp)))
log_message("STEP 1: Data loading complete.")

log_message("Initial data dimensions:")
log_message(paste("  Genes:", nrow(sce)))
log_message(paste("  Cells:", ncol(sce)))
log_message(paste("  Assays:", paste(assayNames(sce), collapse = ", ")))
log_message("Available metadata columns:")
log_message(paste("  ", paste(colnames(colData(sce)), collapse = ", ")))

################################################################################
# 2. QC & Filtering
################################################################################
log_message("STEP 2: QC/filtering started.")
is_mito <- grepl("^MT-|^mt-", rownames(sce))
is_ribo <- grepl("^RPS|^RPL|^rps|^rpl", rownames(sce))
sce <- scater::addPerCellQC(sce, subsets=list(Mito=is_mito, Ribo=is_ribo))

# Log pre-filtering cell and gene counts and group sizes
log_message(paste("Pre-filtering: Genes =", nrow(sce), "Cells =", ncol(sce)))
sce <- check_and_log_factors(sce, stage = "pre-filtering")

# Loosen filtering thresholds (user-configurable/adaptive)
gene_filter <- Matrix::rowSums(counts(sce) > 0) >= 3
cell_filter <- (
  sce$detected >= 100 & sce$detected <= 8000 &
  sce$sum >= 500 & sce$sum <= 150000 &
  sce$subsets_Mito_percent <= 10
)
sce <- sce[gene_filter, cell_filter]

# Log post-filtering cell and gene counts and group sizes
log_message(paste("Post-filtering: Genes =", nrow(sce), "Cells =", ncol(sce)))
sce <- check_and_log_factors(sce, stage = "post-filtering")

# Adaptive filtering warning
cell_loss <- 1 - (ncol(sce) / (ncol(sce) + sum(!cell_filter)))
gene_loss <- 1 - (nrow(sce) / (nrow(sce) + sum(!gene_filter)))
if (cell_loss > 0.8) log_message("WARNING: More than 80% of cells were filtered out. Consider loosening thresholds.")
if (gene_loss > 0.8) log_message("WARNING: More than 80% of genes were filtered out. Consider loosening thresholds.")

# Checkpoint: Save filtered SCE
dir.create(file.path(output_dir, "checkpoints"), showWarnings = FALSE, recursive = TRUE)
save_checkpoint(sce, "sce_filtered", timestamp, output_dir)

saveHDF5SummarizedExperiment(sce, file.path(output_dir, "objects", paste0("sce_filtered_", timestamp)))
log_message("STEP 2: QC/filtering complete.")

################################################################################
# 3. Normalization & HVG Selection
################################################################################
log_message("STEP 3: Normalization and HVG selection started.")
set.seed(42)
clusters <- scran::quickCluster(sce)
sce <- scran::computeSumFactors(sce, clusters=clusters)
sce <- scater::logNormCounts(sce)
gene_var <- scran::modelGeneVar(sce)
hvgs <- scran::getTopHVGs(gene_var, n=2000)
sce_hvg <- sce[hvgs, ]
saveHDF5SummarizedExperiment(sce_hvg, file.path(output_dir, "objects", paste0("sce_hvg_", timestamp)))
log_message(paste("Top HVGs selected:", length(hvgs)))
log_message("STEP 3: Normalization and HVG selection complete.")

################################################################################
# 4. Experimental Design
################################################################################
log_message("STEP 4: Experimental design setup started.")
n_embedding <- 15
test_fraction <- 0.5

cond_levels <- levels(colData(sce)$condition)
sex_levels <- levels(colData(sce)$sex)
region_levels <- levels(colData(sce)$brain_region)

# Modular LEMUR analysis setup
lemur_results <- list()

# 1. Sex-Stratified Analysis (Primary)
for (sex in sex_levels) {
  log_message(paste("Running sex-stratified LEMUR for:", sex))
  sce_sex <- sce[, colData(sce)$sex == sex]
  # Drop unused levels for safety
  sce_sex$condition <- droplevels(sce_sex$condition)
  sce_sex$brain_region <- droplevels(sce_sex$brain_region)
  sce_sex$donor_id <- droplevels(sce_sex$donor_id)
  # Remove cells with NA in any covariate
  covariates_sex <- c("condition", "brain_region", "Age", "PMI", "RIN")
  keep_cells <- complete.cases(as.data.frame(colData(sce_sex)[, covariates_sex, drop=FALSE]))
  sce_sex <- sce_sex[, keep_cells]
  # Design formula: condition + brain_region + Age + PMI + RIN (donor_id removed to avoid collinearity)
  design_formula_sex <- ~ condition + brain_region + Age + PMI + RIN
  fit_sex <- tryCatch({
    lemur(sce_sex, design = design_formula_sex, n_embedding = n_embedding, test_fraction = test_fraction, verbose = TRUE)
  }, error = function(e) {
    log_message(paste("LEMUR fitting error (sex-stratified,", sex, "):", e$message))
    NULL
  })
  if (!is.null(fit_sex)) {
    de_sex <- tryCatch({
      test_de(fit_sex, contrast = cond(condition = cond_levels[2]) - cond(condition = cond_levels[1]), ridge_penalty = 0)
    }, error = function(e) {
      log_message(paste("LEMUR DE error (sex-stratified,", sex, "):", e$message))
      log_message("LEMUR DE extraction failed for sex-stratified model. Check for invalid slot access or internal errors.")
      log_message("Printing structure of fit_sex for debugging:")
      log_message(capture.output(str(fit_sex)))
      log_message("Printing traceback for debugging:")
      log_message(capture.output(traceback()))
      NULL
    })
    lemur_results[[paste0("sex_", sex, "_OUD_vs_Control")]] <- list(fit = fit_sex, de = de_sex)
    # Save only in-memory components of lemur_fit object
    if (!is.null(fit_sex)) {
      saveRDS(fit_sex$coefficients, file.path(output_dir, "objects", paste0("lemur_fit_sex_", sex, "_coefficients_", timestamp, ".rds")))
      saveRDS(fit_sex$embedding, file.path(output_dir, "objects", paste0("lemur_fit_sex_", sex, "_embedding_", timestamp, ".rds")))
      if (!is.null(fit_sex$design)) {
        saveRDS(fit_sex$design, file.path(output_dir, "objects", paste0("lemur_fit_sex_", sex, "_design_", timestamp, ".rds")))
      }
    }
    # Save DE results as RDS
    saveRDS(de_sex, file.path(output_dir, "objects", paste0("lemur_de_sex_", sex, "_", timestamp, ".rds")))
  }
}

# 2. Region-Stratified Analysis (Secondary)
for (region in region_levels) {
  log_message(paste("Running region-stratified LEMUR for:", region))
  sce_region <- sce[, colData(sce)$brain_region == region]
  # Drop unused levels for safety
  sce_region$condition <- droplevels(sce_region$condition)
  sce_region$sex <- droplevels(sce_region$sex)
  sce_region$donor_id <- droplevels(sce_region$donor_id)
  # Remove cells with NA in any covariate
  covariates_region <- c("condition", "sex", "Age", "PMI", "RIN")
  keep_cells <- complete.cases(as.data.frame(colData(sce_region)[, covariates_region, drop=FALSE]))
  sce_region <- sce_region[, keep_cells]
  # Design formula: condition + sex + Age + PMI + RIN (donor_id removed to avoid collinearity)
  design_formula_region <- ~ condition + sex + Age + PMI + RIN
  fit_region <- tryCatch({
    lemur(sce_region, design = design_formula_region, n_embedding = n_embedding, test_fraction = test_fraction, verbose = TRUE)
  }, error = function(e) {
    log_message(paste("LEMUR fitting error (region-stratified,", region, "):", e$message))
    NULL
  })
  if (!is.null(fit_region)) {
    de_region <- tryCatch({
      test_de(fit_region, contrast = cond(condition = cond_levels[2]) - cond(condition = cond_levels[1]), ridge_penalty = 0)
    }, error = function(e) {
      log_message(paste("LEMUR DE error (region-stratified,", region, "):", e$message))
      log_message("LEMUR DE extraction failed for region-stratified model. Check for invalid slot access or internal errors.")
      log_message("Printing structure of fit_region for debugging:")
      log_message(capture.output(str(fit_region)))
      log_message("Printing traceback for debugging:")
      log_message(capture.output(traceback()))
      NULL
    })
    lemur_results[[paste0("region_", region, "_OUD_vs_Control")]] <- list(fit = fit_region, de = de_region)
    # Save only in-memory components of lemur_fit object
    if (!is.null(fit_region)) {
      saveRDS(fit_region$coefficients, file.path(output_dir, "objects", paste0("lemur_fit_region_", region, "_coefficients_", timestamp, ".rds")))
      saveRDS(fit_region$embedding, file.path(output_dir, "objects", paste0("lemur_fit_region_", region, "_embedding_", timestamp, ".rds")))
      if (!is.null(fit_region$design)) {
        saveRDS(fit_region$design, file.path(output_dir, "objects", paste0("lemur_fit_region_", region, "_design_", timestamp, ".rds")))
      }
    }
    # Save DE results as RDS
    saveRDS(de_region, file.path(output_dir, "objects", paste0("lemur_de_region_", region, "_", timestamp, ".rds")))
  }
}

log_message("STEP 4: Modular LEMUR analysis complete.")

################################################################################
# 5. Primary LEMUR Fitting
################################################################################
# (Removed: now handled by modular stratified fitting above)

################################################################################
# 6. Primary DE Analysis
################################################################################
# DE analysis is now handled within the modular LEMUR loops above.
# Remove previous DE analysis block.
log_message("STEP 6: Primary DE analysis complete.")

################################################################################
# 7. Enrichment Analysis
################################################################################
log_message("STEP 7: Enrichment analysis started.")
de_table_male <- as.data.frame(de_results1)
de_table_female <- as.data.frame(de_results2)
de_table_interaction <- as.data.frame(de_interaction)
de_table_male$gene <- rownames(de_table_male)
de_table_female$gene <- rownames(de_table_female)
de_table_interaction$gene <- rownames(de_table_interaction)
enrichment_male_up <- tryCatch({
  perform_enrichment(de_table_male$gene[de_table_male$lfc > 0 & de_table_male$pval < 0.05], rownames(sce_hvg), "male_up")
}, error = function(e) {
  log_message(paste("Enrichment error (male_up):", e$message))
  NULL
})
enrichment_male_down <- tryCatch({
  perform_enrichment(de_table_male$gene[de_table_male$lfc < 0 & de_table_male$pval < 0.05], rownames(sce_hvg), "male_down")
}, error = function(e) {
  log_message(paste("Enrichment error (male_down):", e$message))
  NULL
})
enrichment_female_up <- tryCatch({
  perform_enrichment(de_table_female$gene[de_table_female$lfc > 0 & de_table_female$pval < 0.05], rownames(sce_hvg), "female_up")
}, error = function(e) {
  log_message(paste("Enrichment error (female_up):", e$message))
  NULL
})
enrichment_female_down <- tryCatch({
  perform_enrichment(de_table_female$gene[de_table_female$lfc < 0 & de_table_female$pval < 0.05], rownames(sce_hvg), "female_down")
}, error = function(e) {
  log_message(paste("Enrichment error (female_down):", e$message))
  NULL
})
for (obj_name in c("enrichment_male_up", "enrichment_male_down", "enrichment_female_up", "enrichment_female_down")) {
  if (exists(obj_name) && !is.null(get(obj_name))) {
    save_enrichment_results(get(obj_name), gsub("enrichment_", "", obj_name))
  } else {
    log_message(paste("No enrichment results for", obj_name, "; skipping save/plot."))
  }
}
log_message("STEP 7: Enrichment analysis complete.")

################################################################################
# 8. Secondary Analysis (Region Interaction)
################################################################################
log_message("STEP 8: Secondary analysis (region interaction) started.")
fit2 <- tryCatch({
  lemur(sce_hvg, design = design_formula2, n_embedding = n_embedding, test_fraction = test_fraction, verbose = TRUE)
}, error = function(e) {
  log_message(paste("LEMUR secondary fitting error:", e$message))
  NULL
})
if (!is.null(fit2)) {
  fit2 <- align_harmony(fit2, grouping_vars = vars(donor_id, species))
  de_regional <- tryCatch({
    test_de(fit2,
      contrast =
        (cond(brain_region = region_levels[1], condition = cond_levels[2]) - cond(brain_region = region_levels[1], condition = cond_levels[1])) -
        (cond(brain_region = region_levels[2], condition = cond_levels[2]) - cond(brain_region = region_levels[2], condition = cond_levels[1]))
    )
  }, error = function(e) {
    log_message(paste("LEMUR secondary DE error:", e$message))
    NULL
  })
  saveRDS(fit2, file.path(output_dir, "objects", paste0("lemur_fit_secondary_", timestamp, ".rds")))
  saveRDS(de_regional, file.path(output_dir, "objects", paste0("de_results_regional_", timestamp, ".rds")))
  if (!is.null(de_regional)) {
    de_table_regional <- as.data.frame(de_regional)
    de_table_regional$gene <- rownames(de_table_regional)
    write.csv(de_table_regional, file.path(output_dir, "tables", paste0("de_table_regional_", timestamp, ".csv")), row.names = FALSE)
    volcano_regional <- create_volcano_plot(
      de_table_regional, "Brain Region × Condition Interaction",
      file.path(output_dir, "figures", paste0("volcano_regional_", timestamp, ".pdf"))
    )
  }
  log_message("STEP 8: Secondary analysis complete.")
} else {
  log_message("LEMUR secondary fitting failed; skipping downstream analysis.")
}

log_message("Pipeline completed successfully!")

# Save enrichment results
save_enrichment_results <- function(enrichment_list, prefix) {
  if (is.null(enrichment_list)) {
    return()
  }

  for (analysis_type in names(enrichment_list)) {
    if (!is.null(enrichment_list[[analysis_type]])) {
      result_df <- as.data.frame(enrichment_list[[analysis_type]])
      if (nrow(result_df) > 0) {
        write.csv(result_df,
          file.path(output_dir, "gsea", paste0(prefix, "_", analysis_type, "_", timestamp, ".csv")),
          row.names = FALSE
        )

        # Create enrichment plot
        if (nrow(result_df) >= 5) {
          tryCatch(
            {
              p <- dotplot(enrichment_list[[analysis_type]], showCategory = 15) +
                labs(title = paste(prefix, analysis_type, "Enrichment")) +
                theme_minimal()
              ggsave(file.path(output_dir, "gsea", paste0(prefix, "_", analysis_type, "_plot_", timestamp, ".pdf")),
                p,
                width = 10, height = 8
              )
            },
            error = function(e) {
              log_message(paste("Plot creation failed for", prefix, analysis_type, ":", e$message))
            }
          )
        }
      }
    }
  }
}

# Robustly save all enrichment results (check existence and non-null)
for (obj_name in c("enrichment_male_up", "enrichment_male_down", "enrichment_female_up", "enrichment_female_down")) {
  if (exists(obj_name) && !is.null(get(obj_name))) {
    save_enrichment_results(get(obj_name), gsub("enrichment_", "", obj_name))
  } else {
    log_message(paste("No enrichment results for", obj_name, "; skipping save/plot."))
  }
}

log_message("Gene set enrichment analysis completed")

# =============================================================================
# 12. SECONDARY ANALYSIS: BRAIN REGION INTERACTION
# =============================================================================

log_message("Starting secondary analysis: brain region interaction...")

# Fit secondary LEMUR model with brain_region * condition interaction
set.seed(42)
log_message("STEP 10: Secondary analysis (region interaction) started.")
fit2 <- tryCatch({
  lemur(sce_hvg,
    design = design_formula2, n_embedding = n_embedding,
    test_fraction = test_fraction, verbose = TRUE
  )
}, error = function(e) {
  log_message(paste("LEMUR secondary fitting error:", e$message))
  NULL
})
if (!is.null(fit2)) {
  log_message("LEMUR secondary model fitting complete.")
}

# Align with Harmony
if (!is.null(fit2)) {
  log_message("STEP 10: Secondary Harmony alignment started.")
  fit2 <- align_harmony(fit2, grouping_vars = vars(donor_id, species))
  log_message("LEMUR secondary Harmony alignment complete.")
}

log_message("Secondary LEMUR model fitted and aligned")

# Test regional interaction effects
if (!is.null(fit2)) {
  log_message("STEP 10: Secondary DE analysis started.")
  de_regional <- tryCatch({
    test_de(fit2,
      contrast =
        (cond(brain_region = region_levels[1], condition = cond_levels[2]) - cond(brain_region = region_levels[1], condition = cond_levels[1])) -
          (cond(brain_region = region_levels[2], condition = cond_levels[2]) - cond(brain_region = region_levels[2], condition = cond_levels[1]))
    )
  }, error = function(e) {
    log_message(paste("LEMUR secondary DE error:", e$message))
    NULL
  })
  log_message("LEMUR secondary DE analysis complete.")
}

# Save secondary analysis results
saveRDS(fit2, file.path(output_dir, "objects", paste0("lemur_fit_secondary_", timestamp, ".rds")))
saveRDS(de_regional, file.path(output_dir, "objects", paste0("de_results_regional_", timestamp, ".rds")))

# Create regional interaction DE table
de_table_regional <- as.data.frame(de_regional)
de_table_regional$gene <- rownames(de_table_regional)
write.csv(de_table_regional, file.path(output_dir, "tables", paste0("de_table_regional_", timestamp, ".csv")), row.names = FALSE)

# Regional volcano plot
volcano_regional <- create_volcano_plot(
  de_table_regional, "Brain Region × Condition Interaction",
  file.path(output_dir, "figures", paste0("volcano_regional_", timestamp, ".pdf"))
)

log_message("Secondary analysis completed")

# =============================================================================
# 13. CELL TYPE-SPECIFIC ANALYSIS (OPTIONAL)
# =============================================================================

log_message("Checking for cell type annotations...")

# Check if cell type information is available
cell_type_cols <- colnames(colData(sce_hvg))[grepl("cell_type|cluster|annotation", colnames(colData(sce_hvg)), ignore.case = TRUE)]

if (length(cell_type_cols) > 0) {
  log_message(paste("Found potential cell type columns:", paste(cell_type_cols, collapse = ", ")))

  # Use the first cell type column found
  cell_type_col <- cell_type_cols[1]
  cell_types <- unique(colData(sce_hvg)[[cell_type_col]])

  log_message(paste("Cell types found:", paste(cell_types, collapse = ", ")))

  # Perform cell type-specific analysis for major populations
  major_cell_types <- cell_types[table(colData(sce_hvg)[[cell_type_col]]) >= 100] # At least 100 cells

  if (length(major_cell_types) > 0) {
    log_message(paste("Performing cell type-specific analysis for:", paste(major_cell_types, collapse = ", ")))

    for (ct in major_cell_types) {
      log_message(paste("Analyzing cell type:", ct))

      # Subset to cell type
      cells_ct <- colData(sce_hvg)[[cell_type_col]] == ct
      sce_ct <- sce_hvg[, cells_ct]

      if (ncol(sce_ct) >= 50) { # Minimum 50 cells for analysis
        tryCatch(
          {
            # Fit cell type-specific LEMUR model
            set.seed(42)
            fit_ct <- lemur(sce_ct,
              design = design_formula1, n_embedding = min(10, ncol(sce_ct) / 10),
              test_fraction = test_fraction, verbose = FALSE
            )

            # Test DE
            de_ct <- test_de(fit_ct,
              contrast =
                (cond(sex = sex_levels[1], condition = cond_levels[2]) - cond(sex = sex_levels[1], condition = cond_levels[1])) -
                  (cond(sex = sex_levels[2], condition = cond_levels[2]) - cond(sex = sex_levels[2], condition = cond_levels[1]))
            )

            # Save results
            ct_name <- gsub("[^A-Za-z0-9]", "_", ct)
            saveRDS(fit_ct, file.path(output_dir, "objects", paste0("lemur_fit_", ct_name, "_", timestamp, ".rds")))
            saveRDS(de_ct, file.path(output_dir, "objects", paste0("de_results_", ct_name, "_", timestamp, ".rds")))

            # Create DE table
            de_table_ct <- as.data.frame(de_ct)
            de_table_ct$gene <- rownames(de_table_ct)
            write.csv(de_table_ct, file.path(output_dir, "tables", paste0("de_table_", ct_name, "_", timestamp, ".csv")), row.names = FALSE)

            log_message(paste("Cell type-specific analysis completed for", ct))
          },
          error = function(e) {
            log_message(paste("Cell type-specific analysis failed for", ct, ":", e$message))
          }
        )
      } else {
        log_message(paste("Skipping", ct, "- too few cells (", ncol(sce_ct), ")"))
      }
    }
  }
} else {
  log_message("No cell type annotations found - skipping cell type-specific analysis")
}

# =============================================================================
# 14. ROBUSTNESS CHECKS
# =============================================================================

log_message("Performing robustness checks...")

# Test different embedding dimensions
robustness_results <- list()

for (n_emb in c(10, 20)) {
  log_message(paste("Testing robustness with n_embedding =", n_emb))

  tryCatch(
    {
      set.seed(42)
      fit_robust <- lemur(sce_hvg,
        design = design_formula1, n_embedding = n_emb,
        test_fraction = test_fraction, verbose = FALSE
      )
      fit_robust <- align_harmony(fit_robust, grouping_vars = vars(donor_id, species))

      de_robust <- test_de(fit_robust,
        contrast =
          (cond(sex = sex_levels[1], condition = cond_levels[2]) - cond(sex = sex_levels[1], condition = cond_levels[1])) -
            (cond(sex = sex_levels[2], condition = cond_levels[2]) - cond(sex = sex_levels[2], condition = cond_levels[1]))
      )

      robustness_results[[paste0("n_emb_", n_emb)]] <- de_robust

      # Save robustness results
      saveRDS(de_robust, file.path(output_dir, "objects", paste0("de_results_robust_", n_emb, "_", timestamp, ".rds")))

      log_message(paste("Robustness check completed for n_embedding =", n_emb))
    },
    error = function(e) {
      log_message(paste("Robustness check failed for n_embedding =", n_emb, ":", e$message))
    }
  )
}

# Compare robustness results
if (length(robustness_results) > 0) {
  log_message("Comparing robustness results...")

  # Extract significant genes from each analysis
  robust_sig_genes <- list()
  for (analysis in names(robustness_results)) {
    de_table <- as.data.frame(robustness_results[[analysis]])
    sig_genes <- rownames(de_table)[de_table$adj_pval < 0.05]
    robust_sig_genes[[analysis]] <- sig_genes
  }

  # Add primary analysis results
  primary_sig_genes <- rownames(de_table_interaction)[de_table_interaction$adj_pval < 0.05]
  robust_sig_genes[["primary"]] <- primary_sig_genes

  # Create comparison table
  all_genes <- unique(unlist(robust_sig_genes))
  comparison_matrix <- matrix(FALSE, nrow = length(all_genes), ncol = length(robust_sig_genes))
  rownames(comparison_matrix) <- all_genes
  colnames(comparison_matrix) <- names(robust_sig_genes)

  for (i in 1:length(robust_sig_genes)) {
    comparison_matrix[robust_sig_genes[[i]], i] <- TRUE
  }

  # Save comparison
  write.csv(comparison_matrix, file.path(output_dir, "tables", paste0("robustness_comparison_", timestamp, ".csv")))

  # Calculate overlap statistics
  overlap_stats <- data.frame(
    analysis = names(robust_sig_genes),
    n_sig_genes = sapply(robust_sig_genes, length),
    overlap_with_primary = sapply(robust_sig_genes, function(x) length(intersect(x, primary_sig_genes)))
  )

  write.csv(overlap_stats, file.path(output_dir, "tables", paste0("robustness_overlap_stats_", timestamp, ".csv")), row.names = FALSE)
}

# =============================================================================
# 15. SUMMARY AND FINAL OUTPUTS
# =============================================================================

log_message("Creating final summary...")

# Create comprehensive summary
summary_stats <- data.frame(
  metric = c(
    "Initial genes", "Initial cells", "Filtered genes", "Filtered cells",
    "HVG genes", "LEMUR embedding dimensions", "Test fraction",
    "Male DE genes (adj.p < 0.05)", "Female DE genes (adj.p < 0.05)",
    "Interaction DE genes (adj.p < 0.05)", "Regional DE genes (adj.p < 0.05)"
  ),
  value = c(
    nrow(sce), ncol(sce), nrow(sce_filtered), ncol(sce_filtered),
    nrow(sce_hvg), n_embedding, test_fraction,
    sum(de_table_male$adj_pval < 0.05, na.rm = TRUE),
    sum(de_table_female$adj_pval < 0.05, na.rm = TRUE),
    sum(de_table_interaction$adj_pval < 0.05, na.rm = TRUE),
    sum(de_table_regional$adj_pval < 0.05, na.rm = TRUE)
  )
)

write.csv(summary_stats, file.path(output_dir, "tables", paste0("analysis_summary_", timestamp, ".csv")), row.names = FALSE)

# Create top DE genes summary
create_top_genes_summary <- function(de_table, analysis_name, n_top = 20) {
  # Top upregulated genes
  top_up <- de_table[order(de_table$lfc, decreasing = TRUE), ][1:min(n_top, nrow(de_table)), ]
  top_up$direction <- "Up"
  top_up$analysis <- analysis_name

  # Top downregulated genes
  top_down <- de_table[order(de_table$lfc, decreasing = FALSE), ][1:min(n_top, nrow(de_table)), ]
  top_down$direction <- "Down"
  top_down$analysis <- analysis_name

  return(rbind(top_up, top_down))
}

# Combine top genes from all analyses
top_genes_all <- rbind(
  create_top_genes_summary(de_table_male, "Male_OUD_vs_Control"),
  create_top_genes_summary(de_table_female, "Female_OUD_vs_Control"),
  create_top_genes_summary(de_table_interaction, "Sex_x_Condition_Interaction"),
  create_top_genes_summary(de_table_regional, "Region_x_Condition_Interaction")
)

write.csv(top_genes_all, file.path(output_dir, "tables", paste0("top_genes_summary_", timestamp, ".csv")), row.names = FALSE)

# Session info
session_info <- sessionInfo()
saveRDS(session_info, file.path(output_dir, "objects", paste0("session_info_", timestamp, ".rds")))

# Create final log summary
final_log <- c(
  paste("Analysis completed at:", Sys.time()),
  paste("Total runtime:", difftime(Sys.time(), as.POSIXct(timestamp, format = "%Y%m%d_%H%M%S"), units = "mins"), "minutes"),
  "",
  "=== ANALYSIS SUMMARY ===",
  paste("Input file:", input_file),
  paste("Output directory:", output_dir),
  paste("Initial dimensions:", nrow(sce), "genes x", ncol(sce), "cells"),
  paste("Final dimensions:", nrow(sce_hvg), "genes x", ncol(sce_hvg), "cells"),
  paste("LEMUR embedding dimensions:", n_embedding),
  paste("Test fraction:", test_fraction),
  "",
  "=== DIFFERENTIAL EXPRESSION RESULTS ===",
  paste("Male OUD vs Control DE genes (adj.p < 0.05):", sum(de_table_male$adj_pval < 0.05, na.rm = TRUE)),
  paste("Female OUD vs Control DE genes (adj.p < 0.05):", sum(de_table_female$adj_pval < 0.05, na.rm = TRUE)),
  paste("Sex x Condition interaction DE genes (adj.p < 0.05):", sum(de_table_interaction$adj_pval < 0.05, na.rm = TRUE)),
  paste("Region x Condition interaction DE genes (adj.p < 0.05):", sum(de_table_regional$adj_pval < 0.05, na.rm = TRUE)),
  "",
  "=== OUTPUT FILES ===",
  "Objects saved:",
  paste("  - Initial SCE object:", file.path("objects", paste0("sce_initial_", timestamp, ".rds"))),
  paste("  - Filtered SCE object:", file.path("objects", paste0("sce_filtered_", timestamp, ".rds"))),
  paste("  - Normalized SCE object:", file.path("objects", paste0("sce_normalized_", timestamp, ".rds"))),
  paste("  - HVG SCE object:", file.path("objects", paste0("sce_hvg_", timestamp, ".rds"))),
  paste("  - Primary LEMUR fit:", file.path("objects", paste0("lemur_fit_primary_", timestamp, ".rds"))),
  paste("  - Secondary LEMUR fit:", file.path("objects", paste0("lemur_fit_secondary_", timestamp, ".rds"))),
  "",
  "Tables saved:",
  paste("  - Analysis summary:", file.path("tables", paste0("analysis_summary_", timestamp, ".csv"))),
  paste("  - Top genes summary:", file.path("tables", paste0("top_genes_summary_", timestamp, ".csv"))),
  paste("  - Male DE results:", file.path("tables", paste0("de_table_male_", timestamp, ".csv"))),
  paste("  - Female DE results:", file.path("tables", paste0("de_table_female_", timestamp, ".csv"))),
  paste("  - Interaction DE results:", file.path("tables", paste0("de_table_interaction_", timestamp, ".csv"))),
  paste("  - Regional DE results:", file.path("tables", paste0("de_table_regional_", timestamp, ".csv"))),
  "",
  "Figures saved:",
  paste("  - QC plots:", file.path("figures", paste0("qc_plots_pre_filter_", timestamp, ".pdf"))),
  paste("  - HVG plot:", file.path("figures", paste0("hvg_plot_", timestamp, ".pdf"))),
  paste("  - UMAP plots:", file.path("figures", paste0("umap_plots_", timestamp, ".pdf"))),
  paste("  - Volcano plots:", file.path("figures", paste0("volcano_*_", timestamp, ".pdf"))),
  "",
  "=== R SESSION INFO ===",
  capture.output(print(session_info))
)

# Write final log
writeLines(final_log, file.path(output_dir, paste0("final_analysis_log_", timestamp, ".txt")))

# Print completion message
cat("\n")
cat("=============================================================================\n")
cat("LEMUR snRNA-seq Analysis Pipeline Completed Successfully!\n")
cat("=============================================================================\n")
cat("Output directory:", output_dir, "\n")
cat("Analysis timestamp:", timestamp, "\n")
cat("Total DE genes identified:\n")
cat("  - Male OUD vs Control:", sum(de_table_male$adj_pval < 0.05, na.rm = TRUE), "\n")
cat("  - Female OUD vs Control:", sum(de_table_female$adj_pval < 0.05, na.rm = TRUE), "\n")
cat("  - Sex x Condition interaction:", sum(de_table_interaction$adj_pval < 0.05, na.rm = TRUE), "\n")
cat("  - Region x Condition interaction:", sum(de_table_regional$adj_pval < 0.05, na.rm = TRUE), "\n")
cat("\nAll results saved to:", output_dir, "\n")
cat("=============================================================================\n")

log_message("STEP 2: Quality control and filtering started.")

saveRDS(sce, file.path(output_dir, "objects", paste0("sce_initial_", timestamp, ".rds")))

# Calculate QC metrics
# Mitochondrial genes (human: MT-, NHP: mt-)
is_mito_human <- grepl("^MT-", rownames(sce))
is_mito_nhp <- grepl("^mt-", rownames(sce))
is_mito <- is_mito_human | is_mito_nhp

# Ribosomal genes (RPS, RPL)
is_ribo <- grepl("^RPS|^RPL", rownames(sce))

# Add QC metrics to colData
colData(sce)$n_genes <- colSums(counts(sce) > 0)
colData(sce)$total_counts <- colSums(counts(sce))
colData(sce)$pct_mito <- colSums(counts(sce)[is_mito, ]) / colData(sce)$total_counts * 100
colData(sce)$pct_ribo <- colSums(counts(sce)[is_ribo, ]) / colData(sce)$total_counts * 100

# Add gene-level QC
rowData(sce)$n_cells <- rowSums(counts(sce) > 0)
rowData(sce)$mean_counts <- rowMeans(counts(sce))

# Pre-filtering summary
pre_filter_summary <- data.frame(
  metric = c("n_genes", "n_cells", "median_genes_per_cell", "median_UMI_per_cell",
             "median_pct_mito", "median_pct_ribo"),
  value = c(nrow(sce), ncol(sce), median(colData(sce)$n_genes),
            median(colData(sce)$total_counts), median(colData(sce)$pct_mito),
            median(colData(sce)$pct_ribo))
)

write.csv(pre_filter_summary, file.path(output_dir, "tables", paste0("pre_filter_summary_", timestamp, ".csv")), row.names = FALSE)
log_message("STEP 2: Quality control and filtering complete.")

# Create QC plots
log_message("Generating QC plots...")

# QC violin plots
qc_plots <- list()

qc_plots$genes <- ggplot(data.frame(colData(sce)), aes(x = "All cells", y = n_genes)) +
  geom_violin(fill = "lightblue", alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white") +
  labs(title = "Genes per cell", x = "", y = "Number of genes") +
  theme_minimal()

qc_plots$counts <- ggplot(data.frame(colData(sce)), aes(x = "All cells", y = total_counts)) +
  geom_violin(fill = "lightgreen", alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white") +
  scale_y_log10() +
  labs(title = "UMI counts per cell", x = "", y = "Total UMI counts (log10)") +
  theme_minimal()

qc_plots$mito <- ggplot(data.frame(colData(sce)), aes(x = "All cells", y = pct_mito)) +
  geom_violin(fill = "lightcoral", alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white") +
  labs(title = "Mitochondrial gene %", x = "", y = "% Mitochondrial reads") +
  theme_minimal()

qc_plots$ribo <- ggplot(data.frame(colData(sce)), aes(x = "All cells", y = pct_ribo)) +
  geom_violin(fill = "lightyellow", alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white") +
  labs(title = "Ribosomal gene %", x = "", y = "% Ribosomal reads") +
  theme_minimal()

# Scatter plots
qc_plots$scatter1 <- ggplot(data.frame(colData(sce)), aes(x = total_counts, y = n_genes, color = pct_mito)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(
    title = "Genes vs UMI counts", x = "Total UMI counts (log10)", y = "Number of genes",
    color = "% Mito"
  ) +
  theme_minimal()

qc_plots$scatter2 <- ggplot(data.frame(colData(sce)), aes(x = total_counts, y = pct_mito, color = n_genes)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(
    title = "Mito % vs UMI counts", x = "Total UMI counts (log10)", y = "% Mitochondrial reads",
    color = "N genes"
  ) +
  theme_minimal()

# Save QC plots
pdf(file.path(output_dir, "figures", paste0("qc_plots_pre_filter_", timestamp, ".pdf")), width = 12, height = 8)
plot_grid(plotlist = qc_plots[1:4], ncol = 2, nrow = 2)
plot_grid(plotlist = qc_plots[5:6], ncol = 2, nrow = 1)
dev.off()

# Apply filtering criteria
log_message("Applying filtering criteria...")

# Gene filtering: expressed in at least 3 cells
genes_keep <- rowData(sce)$n_cells >= 3

# Cell filtering
cells_keep <- colData(sce)$n_genes >= 200 &
  colData(sce)$n_genes <= 6000 &
  colData(sce)$total_counts >= 1000 &
  colData(sce)$total_counts <= 100000 &
  colData(sce)$pct_mito <= 5

# Apply filters
sce_filtered <- sce[genes_keep, cells_keep]
sce <- sce_filtered  # Update main object
gc() # Free memory after filtering

# Checkpoint: Save post-filter SCE
save_checkpoint(sce, "sce_post_filter", timestamp, output_dir)

# ---- Drop unused levels, log group counts, and check grouping variables after filtering ----
for (col in c("condition", "sex", "brain_region", "donor_id")) {
  if (col %in% colnames(colData(sce))) {
    colData(sce)[[col]] <- droplevels(as.factor(colData(sce)[[col]]))
    tbl <- table(colData(sce)[[col]])
    log_message(paste("Post-filtering counts for", col, ":", paste(names(tbl), tbl, sep = "=", collapse = ", ")))
    if (col == "donor_id") {
      # Parse donor_id for brain region assignment
      donor_ids <- names(tbl)
      region_from_id <- substr(donor_ids, 1, 1)
      log_message(paste("Donor brain region assignment from donor_id pattern:", paste(donor_ids, region_from_id, sep = "->", collapse = ", ")))
    }
    if (length(tbl) < 2) {
      stop(paste("Grouping variable", col, "has only", length(tbl), "level(s) after filtering. LEMUR requires at least 2."))
    }
    if (any(tbl == 0)) {
      stop(paste("Grouping variable", col, "has a group with zero cells after filtering. Please adjust filtering."))
    }
  }
}
gc() # Free memory after filtering

# ---- Re-factor and check grouping variables after filtering (now modularized above) ----
sce <- check_and_log_factors(sce, stage = "LEMUR pre-fit")

################################################################################
# 4. Normalization & HVG Selection
################################################################################
log_message("STEP 3: Normalization and HVG selection started.")
set.seed(42)
clusters <- scran::quickCluster(sce)
sce <- scran::computeSumFactors(sce, clusters=clusters)
sce <- scater::logNormCounts(sce)
gene_var <- scran::modelGeneVar(sce)
hvgs <- scran::getTopHVGs(gene_var, n=2000)
sce_hvg <- sce[hvgs, ]
saveRDS(sce, file.path(output_dir, "objects", paste0("sce_normalized_", timestamp, ".rds")))
saveHDF5SummarizedExperiment(sce_hvg, file.path(output_dir, "objects", paste0("sce_hvg_", timestamp)))
log_message(paste("Top HVGs selected:", length(hvgs)))
log_message("STEP 3: Normalization and HVG selection complete.")
saveRDS(sce, file.path(output_dir, "objects", paste0("sce_post_qc_pre_enrichment_", timestamp, ".rds")))

log_message("STEP 3: Normalization and HVG selection complete.")

# Save intermediate object for debugging
saveHDF5SummarizedExperiment(sce, file.path(output_dir, "objects", paste0("sce_post_qc_pre_enrichment_", timestamp)))

log_message("STEP 3: Normalization and HVG selection complete.")

# Check that HVG object exists before any downstream analysis
if (!exists("sce_hvg") || is.null(sce_hvg) || !inherits(sce_hvg, "SingleCellExperiment")) {
  log_message("Error: HVG-filtered SCE object 'sce_hvg' not found or invalid. Halting downstream analysis.")
  stop("Error: HVG-filtered SCE object 'sce_hvg' not found or invalid.")
}

# ---- Enrichment analysis (moved here, after HVG selection) ----
background_genes <- rownames(sce_hvg)
enrichment_male_up <- perform_enrichment(gene_lists_male$up, background_genes, "Male Up")
enrichment_male_down <- perform_enrichment(gene_lists_male$down, background_genes, "Male Down")
enrichment_female_up <- perform_enrichment(gene_lists_female$up, background_genes, "Female Up")
enrichment_female_down <- perform_enrichment(gene_lists_female$down, background_genes, "Female Down")
for (obj_name in c("enrichment_male_up", "enrichment_male_down", "enrichment_female_up", "enrichment_female_down")) {
  if (exists(obj_name) && !is.null(get(obj_name))) {
    save_enrichment_results(get(obj_name), gsub("enrichment_", "", obj_name))
  } else {
    log_message(paste("No enrichment results for", obj_name, "; skipping save/plot."))
  }
}
log_message("Gene set enrichment analysis completed")

top100_male_up <- head(gene_lists_male$up, 100)
top100_male_down <- head(gene_lists_male$down, 100)
top100_female_up <- head(gene_lists_female$up, 100)
top100_female_down <- head(gene_lists_female$down, 100)
for (obj_name in c("male_up_top100", "male_down_top100", "female_up_top100", "female_down_top100")) {
  enrich_obj <- perform_enrichment(get(paste0("top100_", obj_name)), background_genes, obj_name)
  if (!is.null(enrich_obj)) {
    save_enrichment_results(enrich_obj, obj_name)
  } else {
    log_message(paste("No enrichment results for", obj_name, "; skipping save/plot."))
  }
}
log_message("Top 100 up/down gene set enrichment results saved")

# Check that HVG object exists before any downstream analysis
if (!exists("sce_hvg") || is.null(sce_hvg) || !inherits(sce_hvg, "SingleCellExperiment")) {
  log_message("Error: HVG-filtered SCE object 'sce_hvg' not found or invalid. Halting downstream analysis.")
  stop("Error: HVG-filtered SCE object 'sce_hvg' not found or invalid.")
}

# =============================================================================
# 6. EXPERIMENTAL DESIGN SETUP
# =============================================================================

log_message("Setting up experimental design...")

# Check required columns in colData
required_cols <- c("donor_id", "species", "brain_region", "sex", "condition")
available_cols <- colnames(colData(sce_hvg))

# Try to match column names (case-insensitive, partial matching)
col_mapping <- list()
for (col in required_cols) {
  matches <- available_cols[grepl(col, available_cols, ignore.case = TRUE)]
  if (length(matches) > 0) {
    col_mapping[[col]] <- matches[1]
    log_message(paste("Mapped", col, "to", matches[1]))
  } else {
    log_message(paste("Warning: Could not find column for", col))
  }
}

# Create standardized column names if mapping exists
if (length(col_mapping) == length(required_cols)) {
  for (i in seq_along(col_mapping)) {
    old_name <- col_mapping[[i]]
    new_name <- names(col_mapping)[i]
    if (old_name != new_name) {
      colData(sce_hvg)[[new_name]] <- colData(sce_hvg)[[old_name]]
    }
  }
} else {
  log_message("Warning: Not all required columns found. Using available columns.")
  log_message(paste("Available columns:", paste(available_cols, collapse = ", ")))
}

# Create design summary
design_summary <- colData(sce_hvg) %>%
  as.data.frame() %>%
  select(any_of(required_cols)) %>%
  summary()

capture.output(design_summary, file = file.path(output_dir, "tables", paste0("design_summary_", timestamp, ".txt")))

# =============================================================================
# 7. LEMUR ANALYSIS PIPELINE
# =============================================================================

log_message("Starting LEMUR analysis...")

# Set LEMUR parameters
log_message("STEP 4: Experimental design setup started.")
n_embedding <- 15
test_fraction <- 0.5
design_formula1 <- ~ donor_id + species + brain_region + sex * condition
design_formula2 <- ~ donor_id + species + brain_region * condition
log_message("STEP 4: Experimental design setup complete.")

# Fit primary LEMUR model
set.seed(42)
fit1 <- lemur(sce_hvg,
  design = design_formula1, n_embedding = n_embedding,
  test_fraction = test_fraction, verbose = TRUE
)

# Align embedding with Harmony
if (!is.null(fit1)) {
  log_message("STEP 5: LEMUR Harmony alignment started.")
  fit1 <- align_harmony(fit1, grouping_vars = vars(donor_id, species))
  log_message("LEMUR Harmony alignment complete.")
}

log_message("Primary LEMUR model fitted and aligned")

# Perform differential expression analysis
de_male <- test_de(fit1, contrast = cond(sex = "male", condition = "OUD") - cond(sex = "male", condition = "Control"))
de_female <- test_de(fit1, contrast = cond(sex = "female", condition = "OUD") - cond(sex = "female", condition = "Control"))
de_interaction <- test_de(fit1,
  contrast =
    (cond(sex = "male", condition = "OUD") - cond(sex = "male", condition = "Control")) -
      (cond(sex = "female", condition = "OUD") - cond(sex = "female", condition = "Control"))
)

# Save model and DE results
saveRDS(fit1, file.path(output_dir, "objects", paste0("lemur_fit_primary_", timestamp, ".rds")))
# Convert to data frames
de_table_male <- as.data.frame(de_male)
de_table_female <- as.data.frame(de_female)
de_table_interaction <- as.data.frame(de_interaction)
de_table_male$gene <- rownames(de_table_male)
de_table_female$gene <- rownames(de_table_female)
de_table_interaction$gene <- rownames(de_table_interaction)

# Save tables
write.csv(de_table_male, file.path(output_dir, "tables", paste0("de_table_male_", timestamp, ".csv")), row.names = FALSE)
write.csv(de_table_female, file.path(output_dir, "tables", paste0("de_table_female_", timestamp, ".csv")), row.names = FALSE)
write.csv(de_table_interaction, file.path(output_dir, "tables", paste0("de_table_interaction_", timestamp, ".csv")), row.names = FALSE)

# Create volcano plots
volcano_male <- create_volcano_plot(
  de_table_male, "Male OUD vs Control",
  file.path(output_dir, "figures", paste0("volcano_male_", timestamp, ".pdf"))
)
volcano_female <- create_volcano_plot(
  de_table_female, "Female OUD vs Control",
  file.path(output_dir, "figures", paste0("volcano_female_", timestamp, ".pdf"))
)
volcano_interaction <- create_volcano_plot(
  de_table_interaction, "Sex × Condition Interaction",
  file.path(output_dir, "figures", paste0("volcano_interaction_", timestamp, ".pdf"))
)

log_message("Primary differential expression analysis completed")

# Perform enrichment analysis for primary DE
enrichment_male_up <- perform_enrichment(
  de_table_male[de_table_male$lfc > 0 & de_table_male$adj_pval < 0.1, "gene"],
  rownames(sce_hvg),
  "male_up"
)
enrichment_male_down <- perform_enrichment(
  de_table_male[de_table_male$lfc < 0 & de_table_male$adj_pval < 0.1, "gene"],
  rownames(sce_hvg),
  "male_down"
)
enrichment_female_up <- perform_enrichment(
  de_table_female[de_table_female$lfc > 0 & de_table_female$adj_pval < 0.1, "gene"],
  rownames(sce_hvg),
  "female_up"
)
enrichment_female_down <- perform_enrichment(
  de_table_female[de_table_female$lfc < 0 & de_table_female$adj_pval < 0.1, "gene"],
  rownames(sce_hvg),
  "female_down"
)

log_message("Primary enrichment analysis completed")
log_message("Fitting primary LEMUR model (sex * condition interaction)...")
set.seed(42)

# Check if we have logcounts
if (!"logcounts" %in% assayNames(sce_hvg)) {
  log_message("Warning: logcounts not found, using counts")
  logcounts(sce_hvg) <- log1p(counts(sce_hvg))
}

# Fit LEMUR model
fit1 <- lemur(sce_hvg,
  design = design_formula1, n_embedding = n_embedding,
  test_fraction = test_fraction, verbose = TRUE
)

log_message("LEMUR model fitted successfully")

# Align with Harmony
log_message("Aligning with Harmony...")
fit1 <- align_harmony(fit1, grouping_vars = vars(donor_id, species))

log_message("Harmony alignment completed")

# Save LEMUR fit
saveRDS(fit1, file.path(output_dir, "objects", paste0("lemur_fit_primary_", timestamp, ".rds")))

# =============================================================================
# 8. DIFFERENTIAL EXPRESSION TESTING
# =============================================================================

log_message("Performing differential expression testing...")

# Test primary contrast: sex-specific OUD effects
# Male OUD vs Male Control
# Female OUD vs Female Control
# Then compare the difference (interaction effect)

log_message("Testing sex-specific OUD effects...")

# Define contrasts for sex-specific effects
set.seed(42)
if (!is.null(fit1)) {
  log_message("STEP 6: LEMUR differential expression started.")
  de_results1 <- tryCatch({
    test_de(fit1, contrast = cond(sex = "male", condition = "OUD") - cond(sex = "male", condition = "Control"))
  }, error = function(e) {
    log_message(paste("LEMUR DE error (male):", e$message))
    NULL
  })
  de_results2 <- tryCatch({
    test_de(fit1, contrast = cond(sex = "female", condition = "OUD") - cond(sex = "female", condition = "Control"))
  }, error = function(e) {
    log_message(paste("LEMUR DE error (female):", e$message))
    NULL
  })
  log_message("LEMUR differential expression complete.")
}

# Test interaction effect (difference in differences)
de_interaction <- test_de(fit1,
  contrast =
    (cond(sex = "male", condition = "OUD") - cond(sex = "male", condition = "Control")) -
      (cond(sex = "female", condition = "OUD") - cond(sex = "female", condition = "Control"))
)

log_message("Differential expression testing completed")

# Save DE results


# =============================================================================
# 9. DE NEIGHBORHOODS ANALYSIS
# =============================================================================

log_message("Finding DE neighborhoods...")

# Find DE neighborhoods
set.seed(42)
if (!is.null(fit1)) {
  log_message("STEP 7: LEMUR DE neighborhoods started.")
  neighborhoods <- tryCatch({
    find_de_neighborhoods(fit1,
                         group_by = vars(donor_id, sex, condition),
                         selection_procedure = "contrast",
                         size_factor_method = "ratio",
                         verbose = TRUE)
  }, error = function(e) {
    log_message(paste("LEMUR neighborhoods error:", e$message))
    NULL
  })
  log_message("LEMUR DE neighborhoods complete.")
}
)

log_message("DE neighborhoods analysis completed")

# Save neighborhoods
saveRDS(neighborhoods, file.path(output_dir, "objects", paste0("de_neighborhoods_", timestamp, ".rds")))

# Extract and save DE tables
de_table_male <- as.data.frame(de_results1)
de_table_female <- as.data.frame(de_results2)
de_table_interaction <- as.data.frame(de_interaction)

# Add gene names and additional information
de_table_male$gene <- rownames(de_table_male)
de_table_female$gene <- rownames(de_table_female)
de_table_interaction$gene <- rownames(de_table_interaction)

# Save DE tables
write.csv(de_table_male, file.path(output_dir, "tables", paste0("de_table_male_", timestamp, ".csv")), row.names = FALSE)
write.csv(de_table_female, file.path(output_dir, "tables", paste0("de_table_female_", timestamp, ".csv")), row.names = FALSE)
write.csv(de_table_interaction, file.path(output_dir, "tables", paste0("de_table_interaction_", timestamp, ".csv")), row.names = FALSE)

# =============================================================================
# 10. VISUALIZATION
# =============================================================================

log_message("Creating visualization plots...")

# Volcano plots
create_volcano_plot <- function(de_table, title, filename) {
  log_message(paste("STEP 8: Creating volcano plot for", title))
  p <- ggplot(de_table, aes(x = lfc, y = -log10(adj_pval))) +
    geom_point(alpha = 0.6, size = 0.8) +
    geom_point(data = subset(de_table, adj_pval < 0.05 & abs(lfc) > 0.5),
               color = "red", size = 1) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "blue") +
    labs(title = title, x = "Log2 Fold Change", y = "-log10(Adjusted P-value)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(filename, p, width = 8, height = 6)
  log_message(paste("Volcano plot saved:", filename))
  return(p)
}

# ============================= PRIMARY LEMUR FITTING AND DE ANALYSIS =============================

log_message("STEP 4: Primary LEMUR fitting started.")
fit1 <- tryCatch({
  lemur(sce_hvg, design = design_formula1, n_embedding = n_embedding, test_fraction = test_fraction, verbose = TRUE)
}, error = function(e) {
  log_message(paste("LEMUR fitting error:", e$message))
  NULL
})
if (!is.null(fit1)) {
  log_message("STEP 4: Primary LEMUR fitting complete.")
  fit1 <- align_harmony(fit1, grouping_vars = vars(donor_id, species))
  log_message("STEP 4: Harmony alignment complete.")

  log_message("STEP 5: Primary DE analysis started.")
  de_results1 <- tryCatch({
    test_de(fit1, contrast = cond(sex = "male", condition = "OUD") - cond(sex = "male", condition = "Control"))
  }, error = function(e) {
    log_message(paste("LEMUR DE error (male):", e$message))
    NULL
  })
  de_results2 <- tryCatch({
    test_de(fit1, contrast = cond(sex = "female", condition = "OUD") - cond(sex = "female", condition = "Control"))
  }, error = function(e) {
    log_message(paste("LEMUR DE error (female):", e$message))
    NULL
  })
  de_interaction <- tryCatch({
    test_de(fit1,
      contrast =
        (cond(sex = "male", condition = "OUD") - cond(sex = "male", condition = "Control")) -
        (cond(sex = "female", condition = "OUD") - cond(sex = "female", condition = "Control"))
    )
  }, error = function(e) {
    log_message(paste("LEMUR DE error (interaction):", e$message))
    NULL
  })
  log_message("STEP 5: Primary DE analysis complete.")

  # Save model and DE results
  saveRDS(fit1, file.path(output_dir, "objects", paste0("lemur_fit_primary_", timestamp, ".rds")))

  # Convert to data frames
  de_table_male <- as.data.frame(de_results1)
  de_table_female <- as.data.frame(de_results2)
  de_table_interaction <- as.data.frame(de_interaction)
  de_table_male$gene <- rownames(de_table_male)
  de_table_female$gene <- rownames(de_table_female)
  de_table_interaction$gene <- rownames(de_table_interaction)

  # Save tables
  write.csv(de_table_male, file.path(output_dir, "tables", paste0("de_table_male_", timestamp, ".csv")), row.names = FALSE)
  write.csv(de_table_female, file.path(output_dir, "tables", paste0("de_table_female_", timestamp, ".csv")), row.names = FALSE)
  write.csv(de_table_interaction, file.path(output_dir, "tables", paste0("de_table_interaction_", timestamp, ".csv")), row.names = FALSE)

  # Volcano plots
  volcano_male <- create_volcano_plot(
    de_table_male, "Male OUD vs Control",
    file.path(output_dir, "figures", paste0("volcano_male_", timestamp, ".pdf"))
  )
  volcano_female <- create_volcano_plot(
    de_table_female, "Female OUD vs Control",
    file.path(output_dir, "figures", paste0("volcano_female_", timestamp, ".pdf"))
  )
  volcano_interaction <- create_volcano_plot(
    de_table_interaction, "Sex × Condition Interaction",
    file.path(output_dir, "figures", paste0("volcano_interaction_", timestamp, ".pdf"))
  )

  log_message("STEP 6: Primary enrichment started.")
  # Example enrichment (replace with your enrichment function as needed)
  enrichment_male_up <- tryCatch({
    perform_enrichment(de_table_male$gene[de_table_male$lfc > 0 & de_table_male$pval < 0.05], rownames(sce_hvg), "male_up")
  }, error = function(e) {
    log_message(paste("Enrichment error (male_up):", e$message))
    NULL
  })
  enrichment_male_down <- tryCatch({
    perform_enrichment(de_table_male$gene[de_table_male$lfc < 0 & de_table_male$pval < 0.05], rownames(sce_hvg), "male_down")
  }, error = function(e) {
    log_message(paste("Enrichment error (male_down):", e$message))
    NULL
  })
  enrichment_female_up <- tryCatch({
    perform_enrichment(de_table_female$gene[de_table_female$lfc > 0 & de_table_female$pval < 0.05], rownames(sce_hvg), "female_up")
  }, error = function(e) {
    log_message(paste("Enrichment error (female_up):", e$message))
    NULL
  })
  enrichment_female_down <- tryCatch({
    perform_enrichment(de_table_female$gene[de_table_female$lfc < 0 & de_table_female$pval < 0.05], rownames(sce_hvg), "female_down")
  }, error = function(e) {
    log_message(paste("Enrichment error (female_down):", e$message))
    NULL
  })
  for (obj_name in c("enrichment_male_up", "enrichment_male_down", "enrichment_female_up", "enrichment_female_down")) {
    if (exists(obj_name) && !is.null(get(obj_name))) {
      save_enrichment_results(get(obj_name), gsub("enrichment_", "", obj_name))
    } else {
      log_message(paste("No enrichment results for", obj_name, "; skipping save/plot."))
    }
  }
  log_message("STEP 6: Primary enrichment complete.")

} else {
  log_message("LEMUR fitting failed; skipping downstream analysis.")
}

# ============================= SECONDARY ANALYSIS (REGION INTERACTION) =============================

log_message("STEP 7: Secondary analysis (region interaction) started.")
fit2 <- tryCatch({
  lemur(sce_hvg,
    design = design_formula2, n_embedding = n_embedding,
    test_fraction = test_fraction, verbose = TRUE
  )
}, error = function(e) {
  log_message(paste("LEMUR secondary fitting error:", e$message))
  NULL
})
if (!is.null(fit2)) {
  log_message("LEMUR secondary model fitting complete.")
  fit2 <- align_harmony(fit2, grouping_vars = vars(donor_id, species))
  log_message("LEMUR secondary Harmony alignment complete.")

  de_regional <- tryCatch({
    test_de(fit2,
      contrast =
        (cond(brain_region = "caudate", condition = "OUD") - cond(brain_region = "caudate", condition = "Control")) -
        (cond(brain_region = "putamen", condition = "OUD") - cond(brain_region = "putamen", condition = "Control"))
    )
  }, error = function(e) {
    log_message(paste("LEMUR secondary DE error:", e$message))
    NULL
  })
  log_message("LEMUR secondary DE analysis complete.")

  saveRDS(fit2, file.path(output_dir, "objects", paste0("lemur_fit_secondary_", timestamp, ".rds")))
  saveRDS(de_regional, file.path(output_dir, "objects", paste0("de_results_regional_", timestamp, ".rds")))

  if (!is.null(de_regional)) {
    de_table_regional <- as.data.frame(de_regional)
    de_table_regional$gene <- rownames(de_table_regional)
    write.csv(de_table_regional, file.path(output_dir, "tables", paste0("de_table_regional_", timestamp, ".csv")), row.names = FALSE)
    volcano_regional <- create_volcano_plot(
      de_table_regional, "Brain Region × Condition Interaction",
      file.path(output_dir, "figures", paste0("volcano_regional_", timestamp, ".pdf"))
    )
  }
} else {
  log_message("LEMUR secondary fitting failed; skipping downstream analysis.")
}
log_message("STEP 7: Secondary analysis (region interaction) complete.")

# ============================= PRIMARY LEMUR FITTING AND DE ANALYSIS =============================

log_message("STEP 5: Primary LEMUR fitting started.")
fit1 <- tryCatch({
  lemur(sce_hvg, design = design_formula1, n_embedding = n_embedding, test_fraction = test_fraction, verbose = TRUE)
}, error = function(e) {
  log_message(paste("LEMUR fitting error:", e$message))
  NULL
})
if (!is.null(fit1)) {
  log_message("STEP 5: Primary LEMUR fitting complete.")
  fit1 <- align_harmony(fit1, grouping_vars = vars(donor_id, species))
  log_message("STEP 5: Harmony alignment complete.")
  saveRDS(fit1, file.path(output_dir, "objects", paste0("lemur_fit_primary_", timestamp, ".rds")))
} else {
  stop("LEMUR fitting failed; pipeline halted.")
}

# ============================= PRIMARY DE ANALYSIS =============================

log_message("STEP 6: Primary DE analysis started.")
de_results1 <- tryCatch({
  test_de(fit1, contrast = cond(sex = "male", condition = "OUD") - cond(sex = "male", condition = "Control"))
}, error = function(e) {
  log_message(paste("LEMUR DE error (male):", e$message))
  NULL
})
de_results2 <- tryCatch({
  test_de(fit1, contrast = cond(sex = "female", condition = "OUD") - cond(sex = "female", condition = "Control"))
}, error = function(e) {
  log_message(paste("LEMUR DE error (female):", e$message))
  NULL
})
de_interaction <- tryCatch({
  test_de(fit1,
    contrast =
      (cond(sex = "male", condition = "OUD") - cond(sex = "male", condition = "Control")) -
      (cond(sex = "female", condition = "OUD") - cond(sex = "female", condition = "Control"))
  )
}, error = function(e) {
  log_message(paste("LEMUR DE error (interaction):", e$message))
  NULL
})
log_message("STEP 6: Primary DE analysis complete.")

# ============================= ENRICHMENT ANALYSIS =============================

log_message("STEP 7: Enrichment analysis started.")
# Example enrichment (replace with your enrichment function as needed)
de_table_male <- as.data.frame(de_results1)
de_table_female <- as.data.frame(de_results2)
de_table_interaction <- as.data.frame(de_interaction)
de_table_male$gene <- rownames(de_table_male)
de_table_female$gene <- rownames(de_table_female)
de_table_interaction$gene <- rownames(de_table_interaction)

enrichment_male_up <- tryCatch({
  perform_enrichment(de_table_male$gene[de_table_male$lfc > 0 & de_table_male$pval < 0.05], rownames(sce_hvg), "male_up")
}, error = function(e) {
  log_message(paste("Enrichment error (male_up):", e$message))
  NULL
})
enrichment_male_down <- tryCatch({
  perform_enrichment(de_table_male$gene[de_table_male$lfc < 0 & de_table_male$pval < 0.05], rownames(sce_hvg), "male_down")
}, error = function(e) {
  log_message(paste("Enrichment error (male_down):", e$message))
  NULL
})
enrichment_female_up <- tryCatch({
  perform_enrichment(de_table_female$gene[de_table_female$lfc > 0 & de_table_female$pval < 0.05], rownames(sce_hvg), "female_up")
}, error = function(e) {
  log_message(paste("Enrichment error (female_up):", e$message))
  NULL
})
enrichment_female_down <- tryCatch({
  perform_enrichment(de_table_female$gene[de_table_female$lfc < 0 & de_table_female$pval < 0.05], rownames(sce_hvg), "female_down")
}, error = function(e) {
  log_message(paste("Enrichment error (female_down):", e$message))
  NULL
})
for (obj_name in c("enrichment_male_up", "enrichment_male_down", "enrichment_female_up", "enrichment_female_down")) {
  if (exists(obj_name) && !is.null(get(obj_name))) {
    save_enrichment_results(get(obj_name), gsub("enrichment_", "", obj_name))
  } else {
    log_message(paste("No enrichment results for", obj_name, "; skipping save/plot."))
  }
}
log_message("STEP 7: Enrichment analysis complete.")

# ============================= SECONDARY ANALYSIS (REGION INTERACTION) =============================

log_message("STEP 8: Secondary analysis (region interaction) started.")
fit2 <- tryCatch({
  lemur(sce_hvg, design = design_formula2, n_embedding = n_embedding, test_fraction = test_fraction, verbose = TRUE)
}, error = function(e) {
  log_message(paste("LEMUR secondary fitting error:", e$message))
  NULL
})
if (!is.null(fit2)) {
  fit2 <- align_harmony(fit2, grouping_vars = vars(donor_id, species))
  de_regional <- tryCatch({
    test_de(fit2,
      contrast =
        (cond(brain_region = "caudate", condition = "OUD") - cond(brain_region = "caudate", condition = "Control")) -
        (cond(brain_region = "putamen", condition = "OUD") - cond(brain_region = "putamen", condition = "Control"))
    )
  }, error = function(e) {
    log_message(paste("LEMUR secondary DE error:", e$message))
    NULL
  })
  saveRDS(fit2, file.path(output_dir, "objects", paste0("lemur_fit_secondary_", timestamp, ".rds")))
  saveRDS(de_regional, file.path(output_dir, "objects", paste0("de_results_regional_", timestamp, ".rds")))
  if (!is.null(de_regional)) {
    de_table_regional <- as.data.frame(de_regional)
    de_table_regional$gene <- rownames(de_table_regional)
    write.csv(de_table_regional, file.path(output_dir, "tables", paste0("de_table_regional_", timestamp, ".csv")), row.names = FALSE)
    volcano_regional <- create_volcano_plot(
      de_table_regional, "Brain Region × Condition Interaction",
      file.path(output_dir, "figures", paste0("volcano_regional_", timestamp, ".pdf"))
    )
  }
  log_message("STEP 8: Secondary analysis complete.")
} else {
  log_message("LEMUR secondary fitting failed; skipping downstream analysis.")
}

log_message("Pipeline completed successfully!")

# UMAP plots
log_message("Creating UMAP plots...")

# Generate UMAP embedding
set.seed(42)
sce_hvg <- runUMAP(sce_hvg, dimred = "PCA", n_dimred = 15)

# Create UMAP plots by metadata
if (all(c("sex", "condition", "brain_region") %in% colnames(colData(sce_hvg)))) {
  umap_plots <- list()

  umap_plots$sex <- plotReducedDim(sce_hvg, dimred = "UMAP", colour_by = "sex") +
    labs(title = "UMAP by Sex") +
    theme_minimal()

  umap_plots$condition <- plotReducedDim(sce_hvg, dimred = "UMAP", colour_by = "condition") +
    labs(title = "UMAP by Condition") +
    theme_minimal()

  umap_plots$region <- plotReducedDim(sce_hvg, dimred = "UMAP", colour_by = "brain_region") +
    labs(title = "UMAP by Brain Region") +
    theme_minimal()

  # Save UMAP plots
  pdf(file.path(output_dir, "figures", paste0("umap_plots_", timestamp, ".pdf")), width = 12, height = 4)
  plot_grid(plotlist = umap_plots, ncol = 3, nrow = 1)
  dev.off()
}

# Gene-specific UMAP plots for key genes
key_genes <- c("GAP43", "CXCL8", "RPS11")
available_key_genes <- intersect(key_genes, rownames(sce_hvg))

if (length(available_key_genes) > 0) {
  log_message(paste("Creating gene-specific UMAP plots for:", paste(available_key_genes, collapse = ", ")))

  pdf(file.path(output_dir, "figures", paste0("umap_genes_", timestamp, ".pdf")), width = 12, height = 4)
  for (gene in available_key_genes) {
    p <- plotReducedDim(sce_hvg, dimred = "UMAP", colour_by = gene) +
      labs(title = paste("UMAP -", gene, "expression")) +
      theme_minimal()
    print(p)
  }
  dev.off()
}

# =============================================================================
# 11. GENE SET ENRICHMENT ANALYSIS
# =============================================================================

log_message("Performing gene set enrichment analysis...")

# Prepare gene lists for enrichment
prepare_gene_list <- function(de_table, adj_pval_cutoff = 0.1, lfc_cutoff = 0) {
  sig_genes <- de_table[de_table$adj_pval < adj_pval_cutoff & abs(de_table$lfc) > lfc_cutoff, ]
  up_genes <- sig_genes[sig_genes$lfc > 0, ]$gene
  down_genes <- sig_genes[sig_genes$lfc < 0, ]$gene
  all_genes <- sig_genes$gene

  return(list(up = up_genes, down = down_genes, all = all_genes))
}

# Prepare gene lists
gene_lists_male <- prepare_gene_list(de_table_male)
gene_lists_female <- prepare_gene_list(de_table_female)
gene_lists_interaction <- prepare_gene_list(de_table_interaction)

log_message(paste(
  "Male DE genes - Up:", length(gene_lists_male$up),
  "Down:", length(gene_lists_male$down),
  "Total:", length(gene_lists_male$all)
))
log_message(paste(
  "Female DE genes - Up:", length(gene_lists_female$up),
  "Down:", length(gene_lists_female$down),
  "Total:", length(gene_lists_female$all)
))

# Function to perform enrichment analysis
perform_enrichment <- function(gene_list, background_genes, analysis_name) {
  if (length(gene_list) == 0) {
    log_message(paste("No genes for", analysis_name, "- skipping"))
    return(NULL)
  }

  results <- list()

  # GO Biological Process
  tryCatch(
    {
      results$GO_BP <- enrichGO(
        gene = gene_list,
        universe = background_genes,
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        keyType = "SYMBOL"
      )
    },
    error = function(e) {
      log_message(paste("GO BP enrichment failed for", analysis_name, ":", e$message))
    }
  )

  # KEGG
  tryCatch(
    {
      # Convert to Entrez IDs
      entrez_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
      if (nrow(entrez_ids) > 0) {
        results$KEGG <- enrichKEGG(
          gene = entrez_ids$ENTREZID,
          organism = "hsa",
          pvalueCutoff = 0.05,
          pAdjustMethod = "BH"
        )
      }
    },
    error = function(e) {
      log_message(paste("KEGG enrichment failed for", analysis_name, ":", e$message))
    }
  )

  # Hallmark gene sets
  tryCatch(
    {
      hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
      results$Hallmark <- enricher(
        gene = gene_list,
        TERM2GENE = hallmark_sets[, c("gs_name", "gene_symbol")],
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH"
      )
    },
    error = function(e) {
      log_message(paste("Hallmark enrichment failed for", analysis_name, ":", e$message))
    }
  )
# Robust DE filtering, logging, and safe enrichment with error handling

background_genes <- rownames(sce_hvg)

# Example: de_results is a named list of DE tables (male, female, interaction, etc.)
message("Available DE results: ", paste(names(de_results), collapse = ", "))
if (length(de_results) == 0) {
  stop("❌ No DE results found — check earlier DE step.")
}

for (test_name in names(de_results)) {
  tryCatch({
    de_table <- as.data.frame(de_results[[test_name]])
    message(test_name, " DE genes before filtering: ", nrow(de_table))

    # Less strict filtering
    de_filtered <- de_table[
      abs(de_table$effect_size) > 0.15 &
      de_table$pval < 0.05 &
      !is.na(de_table$gene),
    ]

    message(test_name, " DE genes after filtering: ", nrow(de_filtered))

    if (is.null(de_filtered) || nrow(de_filtered) == 0) {
      message("⚠️ No DE genes passed the filter for: ", test_name)
      next
    }

    # Save filtered DE table
    write.csv(de_filtered, file.path(output_dir, "tables", paste0("de_table_filtered_", test_name, "_", timestamp, ".csv")), row.names = FALSE)

    # Safe enrichment
    enrich_obj <- tryCatch({
      perform_enrichment(de_filtered$gene, background_genes, test_name)
    }, error = function(e) {
      message("❌ Error in enrichment for ", test_name, ": ", e$message)
      NULL
    })

    if (!is.null(enrich_obj)) {
      save_enrichment_results(enrich_obj, paste0(test_name, "_filtered"))
    } else {
      message("No enrichment results for ", test_name, " ; skipping save/plot.")
    }

  }, error = function(e) {
    message("❌ Error in DE/enrichment for ", test_name, ": ", e$message)
  })
}

log_message("DE filtering and enrichment analysis completed")

# Save all major objects and tables with timestamped filenames
saveRDS(sce, file.path(output_dir, "objects", paste0("sce_initial_", timestamp, ".rds")))
saveHDF5SummarizedExperiment(sce_filtered, file.path(output_dir, "objects", paste0("sce_filtered_", timestamp)))
saveRDS(sce_hvg, file.path(output_dir, "objects", paste0("sce_hvg_", timestamp, ".rds")))
saveRDS(fit1, file.path(output_dir, "objects", paste0("lemur_fit_primary_", timestamp, ".rds")))
saveRDS(fit2, file.path(output_dir, "objects", paste0("lemur_fit_secondary_", timestamp, ".rds")))
saveRDS(de_regional, file.path(output_dir, "objects", paste0("de_results_regional_", timestamp, ".rds")))
saveRDS(neighborhoods, file.path(output_dir, "objects", paste0("de_neighborhoods_", timestamp, ".rds")))
saveRDS(session_info, file.path(output_dir, "objects", paste0("session_info_", timestamp, ".rds")))

write.csv(de_table_male, file.path(output_dir, "tables", paste0("de_table_male_", timestamp, ".csv")), row.names = FALSE)
write.csv(de_table_female, file.path(output_dir, "tables", paste0("de_table_female_", timestamp, ".csv")), row.names = FALSE)
write.csv(de_table_interaction, file.path(output_dir, "tables", paste0("de_table_interaction_", timestamp, ".csv")), row.names = FALSE)
write.csv(de_table_regional, file.path(output_dir, "tables", paste0("de_table_regional_", timestamp, ".csv")), row.names = FALSE)
write.csv(top_genes_all, file.path(output_dir, "tables", paste0("top_genes_summary_", timestamp, ".csv")), row.names = FALSE)
write.csv(summary_stats, file.path(output_dir, "tables", paste0("analysis_summary_", timestamp, ".csv")), row.names = FALSE)
write.csv(filter_comparison, file.path(output_dir, "tables", paste0("filter_comparison_", timestamp, ".csv")), row.names = FALSE)
write.csv(pre_filter_summary, file.path(output_dir, "tables", paste0("pre_filter_summary_", timestamp, ".csv")), row.names = FALSE)
write.csv(post_filter_summary, file.path(output_dir, "tables", paste0("post_filter_summary_", timestamp, ".csv")), row.names = FALSE)

# Save robustness results
if (exists("robustness_results")) {
  for (n_emb in names(robustness_results)) {
    saveRDS(robustness_results[[n_emb]], file.path(output_dir, "objects", paste0("de_results_robust_", n_emb, "_", timestamp, ".rds")))
  }
  if (exists("comparison_matrix")) {
    write.csv(comparison_matrix, file.path(output_dir, "tables", paste0("robustness_comparison_", timestamp, ".csv")))
  }
  if (exists("overlap_stats")) {
    write.csv(overlap_stats, file.path(output_dir, "tables", paste0("robustness_overlap_stats_", timestamp, ".csv")), row.names = FALSE)
  }
}

# Save cell type-specific results if present
if (exists("de_table_ct")) {
  write.csv(de_table_ct, file.path(output_dir, "tables", paste0("de_table_", ct_name, "_", timestamp, ".csv")), row.names = FALSE)
  saveRDS(fit_ct, file.path(output_dir, "objects", paste0("lemur_fit_", ct_name, "_", timestamp, ".rds")))
  saveRDS(de_ct, file.path(output_dir, "objects", paste0("de_results_", ct_name, "_", timestamp, ".rds")))
}

# Save final log summary
writeLines(final_log, file.path(output_dir, paste0("final_analysis_log_", timestamp, ".txt")))

log_message("All outputs saved and pipeline completed successfully.")
