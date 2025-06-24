#!/usr/bin/env Rscript
#' @title H5AD Data Validation Script
#' @description Quick validation of H5AD file structure and metadata
#' @author Generated Validation Script
#' @date 2024

# Validate H5AD file structure and metadata for LEMUR analysis
cat("🔍 H5AD Data Validation for LEMUR Analysis\n")
cat("==========================================\n\n")

# Suppress warnings for cleaner output
suppressPackageStartupMessages({
  library(zellkonverter)
  library(SingleCellExperiment)
})

# File paths
BASE_DIR <- "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study"
INPUT_H5AD <- file.path(BASE_DIR, "data/raw/snrna/GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5ad")

cat("📁 Input file:", INPUT_H5AD, "\n\n")

# Check if file exists
if (!file.exists(INPUT_H5AD)) {
  cat("❌ ERROR: H5AD file not found!\n")
  cat("Expected location:", INPUT_H5AD, "\n")
  quit(status = 1)
}

# Get file info
file_info <- file.info(INPUT_H5AD)
cat("📊 File Information:\n")
cat("  Size:", round(file_info$size / 1024^3, 2), "GB\n")
cat("  Modified:", as.character(file_info$mtime), "\n\n")

# Load the H5AD file
cat("🔄 Loading H5AD file...\n")
tryCatch(
  {
    sce <- readH5AD(INPUT_H5AD, use_hdf5 = TRUE)
    cat("✅ File loaded successfully!\n\n")
  },
  error = function(e) {
    cat("❌ ERROR loading H5AD file:", e$message, "\n")
    quit(status = 1)
  }
)

# Basic data structure
cat("📈 Data Structure:\n")
cat("  Genes (features):", nrow(sce), "\n")
cat("  Cells (observations):", ncol(sce), "\n")
cat("  Assays available:", paste(assayNames(sce), collapse = ", "), "\n\n")

# Check metadata columns
cat("🏷️  Metadata Columns:\n")
metadata_cols <- colnames(colData(sce))
cat("  Total columns:", length(metadata_cols), "\n")
cat("  Available columns:\n")
for (col in metadata_cols) {
  cat("    -", col, "\n")
}
cat("\n")

# Check required columns
required_cols <- c("donor_id", "condition")
cat("✅ Required Column Check:\n")
for (col in required_cols) {
  if (col %in% metadata_cols) {
    cat("  ✅", col, "- FOUND\n")
  } else {
    cat("  ❌", col, "- MISSING\n")
  }
}
cat("\n")

# Analyze condition column
if ("condition" %in% metadata_cols) {
  cat("🎯 Condition Analysis:\n")
  condition_table <- table(sce$condition)
  cat("  Unique conditions:", length(unique(sce$condition)), "\n")
  cat("  Condition counts:\n")
  for (cond in names(condition_table)) {
    cat("    -", cond, ":", condition_table[cond], "cells\n")
  }

  # Check for expected conditions
  expected_conditions <- c("Control", "OUD")
  missing_conditions <- setdiff(expected_conditions, names(condition_table))
  if (length(missing_conditions) > 0) {
    cat("  ⚠️  Missing expected conditions:", paste(missing_conditions, collapse = ", "), "\n")
  } else {
    cat("  ✅ All expected conditions found\n")
  }
  cat("\n")
}

# Analyze donor_id column
if ("donor_id" %in% metadata_cols) {
  cat("👥 Donor Analysis:\n")
  donor_table <- table(sce$donor_id)
  cat("  Unique donors:", length(unique(sce$donor_id)), "\n")
  cat("  Cells per donor:\n")
  donor_summary <- summary(as.numeric(donor_table))
  cat("    Min:", donor_summary["Min."], "\n")
  cat("    Median:", donor_summary["Median"], "\n")
  cat("    Mean:", round(donor_summary["Mean"], 1), "\n")
  cat("    Max:", donor_summary["Max."], "\n")

  # Show donor distribution
  if (length(unique(sce$donor_id)) <= 20) {
    cat("  Donor distribution:\n")
    for (donor in names(sort(donor_table, decreasing = TRUE))) {
      cat("    -", donor, ":", donor_table[donor], "cells\n")
    }
  }
  cat("\n")
}

# Cross-tabulation of condition and donor
if ("condition" %in% metadata_cols && "donor_id" %in% metadata_cols) {
  cat("🔄 Condition × Donor Cross-tabulation:\n")
  cross_tab <- table(sce$condition, sce$donor_id)
  print(cross_tab)
  cat("\n")
}

# Check for additional useful columns
cat("🔍 Additional Metadata Exploration:\n")
useful_cols <- c("Sex", "Age", "cell_type", "cluster", "seurat_clusters", "level1", "level2")
found_useful <- intersect(useful_cols, metadata_cols)
if (length(found_useful) > 0) {
  cat("  Found useful columns:\n")
  for (col in found_useful) {
    n_unique <- length(unique(colData(sce)[[col]]))
    cat("    -", col, ":", n_unique, "unique values\n")

    # Show distribution for categorical variables with reasonable number of levels
    if (n_unique <= 10) {
      col_table <- table(colData(sce)[[col]])
      cat("      Distribution:", paste(names(col_table), "=", col_table, collapse = ", "), "\n")
    }
  }
} else {
  cat("  No additional useful columns found\n")
}
cat("\n")

# Gene information
cat("🧬 Gene Information:\n")
if ("Symbol" %in% colnames(rowData(sce))) {
  cat("  ✅ Gene symbols available\n")
} else {
  cat("  ⚠️  Gene symbols not found in rowData\n")
}

# Check for mitochondrial genes
mito_genes <- grep("^MT-|^mt-", rownames(sce), value = TRUE)
cat("  Mitochondrial genes found:", length(mito_genes), "\n")
if (length(mito_genes) <= 20) {
  cat("    Examples:", paste(head(mito_genes, 10), collapse = ", "), "\n")
}

# Check for ribosomal genes
ribo_genes <- grep("^RP[SL]|^Rp[sl]", rownames(sce), value = TRUE)
cat("  Ribosomal genes found:", length(ribo_genes), "\n")
cat("\n")

# Expression data check
cat("📊 Expression Data Check:\n")
# Check if counts are integers (raw counts)
first_100_cells <- counts(sce)[1:min(100, nrow(sce)), 1:min(100, ncol(sce))]
has_non_integers <- any(first_100_cells != round(first_100_cells))
cat("  Contains non-integer values:", has_non_integers, "\n")

# Check sparsity
zero_fraction <- sum(first_100_cells == 0) / length(first_100_cells)
cat("  Sparsity (% zeros in sample):", round(zero_fraction * 100, 1), "%\n")

# Check value ranges
cat("  Count range in sample: [", min(first_100_cells), ",", max(first_100_cells), "]\n")
cat("\n")

# Memory usage
cat("💾 Memory Information:\n")
object_size <- object.size(sce)
cat("  Object size:", round(as.numeric(object_size) / 1024^3, 2), "GB\n")
cat("\n")

# Summary and recommendations
cat("📋 VALIDATION SUMMARY\n")
cat("====================\n")

# Check readiness for LEMUR analysis
ready_for_lemur <- TRUE
issues <- c()

# Check required columns
if (!"condition" %in% metadata_cols) {
  ready_for_lemur <- FALSE
  issues <- c(issues, "Missing 'condition' column")
}

if (!"donor_id" %in% metadata_cols) {
  ready_for_lemur <- FALSE
  issues <- c(issues, "Missing 'donor_id' column")
}

# Check condition values
if ("condition" %in% metadata_cols) {
  conditions <- unique(sce$condition)
  if (!all(c("Control", "OUD") %in% conditions)) {
    issues <- c(issues, "Expected conditions 'Control' and 'OUD' not found")
  }
}

# Check data size
if (ncol(sce) < 1000) {
  issues <- c(issues, "Very few cells (< 1000) - results may be unreliable")
}

if (nrow(sce) < 10000) {
  issues <- c(issues, "Few genes (< 10,000) - consider gene filtering")
}

# Final assessment
if (ready_for_lemur && length(issues) == 0) {
  cat("🎉 READY FOR LEMUR ANALYSIS!\n")
  cat("✅ All required metadata columns present\n")
  cat("✅ Expected conditions found\n")
  cat("✅ Sufficient data for analysis\n\n")

  cat("🚀 Next steps:\n")
  cat("  1. Run: ./run_lemur.sh --install-packages\n")
  cat("  2. Or: Rscript run_LEMUR_GSE225158.R\n")
} else {
  cat("⚠️  ISSUES DETECTED:\n")
  for (issue in issues) {
    cat("  ❌", issue, "\n")
  }
  cat("\n")
  cat("🔧 Please address these issues before running LEMUR analysis\n")
}

cat("\n")
cat("📞 For help, see README_LEMUR.md\n")
cat("🔍 Validation completed at:", as.character(Sys.time()), "\n")
