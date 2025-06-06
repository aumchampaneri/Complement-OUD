#!/usr/bin/env Rscript

# =============================================================================
# Validation Script for edgeR Files
# =============================================================================
# 
# Description: Quick validation to test loading your specific edgeR files
#              and verify they have the correct structure for TF analysis
# 
# Author: Multi-Omics Study Team
# Date: 2024
# 
# Usage: Rscript validate_edger_files.R
# =============================================================================

cat("🔍 Validating edgeR Files for TF Activity Analysis...\n")
cat("===================================================\n\n")

# Load required packages
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

# Define the correct path to your edgeR results
edger_path <- "../../results/snrna_scvi/differential_expression_edgeR"

cat("📁 Checking directory:", edger_path, "\n")

if (!dir.exists(edger_path)) {
  cat("❌ Directory not found! Please check the path.\n")
  cat("   Current working directory:", getwd(), "\n")
  stop("Cannot proceed without edgeR results directory")
}

# Find edgeR results files
edger_files <- list.files(edger_path, 
                         pattern = "edgeR.*results.*\\.csv$", 
                         full.names = TRUE)

cat("📊 Found", length(edger_files), "edgeR results files:\n")
for (file in edger_files) {
  cat("  -", basename(file), "\n")
}

if (length(edger_files) == 0) {
  cat("❌ No edgeR results files found!\n")
  stop("No files to validate")
}

cat("\n🧪 Testing file loading and structure...\n")

# Test loading each file
validation_results <- list()

for (i in seq_along(edger_files)) {
  file <- edger_files[i]
  filename <- basename(file)
  
  cat("\n", i, ". Testing:", filename, "\n")
  
  # Try to load the file
  tryCatch({
    data <- read_csv(file, show_col_types = FALSE)
    
    # Check basic structure
    cat("   ✅ File loaded successfully\n")
    cat("   📏 Dimensions:", nrow(data), "rows x", ncol(data), "columns\n")
    
    # Check required columns
    required_cols <- c("gene", "logFC", "PValue", "FDR")
    present_cols <- required_cols[required_cols %in% names(data)]
    missing_cols <- required_cols[!required_cols %in% names(data)]
    
    cat("   📋 Required columns present:", paste(present_cols, collapse = ", "), "\n")
    if (length(missing_cols) > 0) {
      cat("   ⚠️  Missing columns:", paste(missing_cols, collapse = ", "), "\n")
    }
    
    # Check data types and ranges
    if ("logFC" %in% names(data)) {
      logfc_range <- range(data$logFC, na.rm = TRUE)
      cat("   📈 logFC range:", round(logfc_range[1], 3), "to", round(logfc_range[2], 3), "\n")
    }
    
    if ("PValue" %in% names(data)) {
      pval_range <- range(data$PValue, na.rm = TRUE)
      n_sig <- sum(data$PValue < 0.05, na.rm = TRUE)
      cat("   📊 P-value range:", signif(pval_range[1], 3), "to", signif(pval_range[2], 3), "\n")
      cat("   🎯 Significant genes (p < 0.05):", n_sig, "\n")
    }
    
    if ("gene" %in% names(data)) {
      n_genes <- length(unique(data$gene))
      cat("   🧬 Unique genes:", n_genes, "\n")
      
      # Show first few gene names
      sample_genes <- head(data$gene, 5)
      cat("   📝 Sample genes:", paste(sample_genes, collapse = ", "), "\n")
    }
    
    # Check for contrast information
    if ("contrast" %in% names(data)) {
      contrasts <- unique(data$contrast)
      cat("   🔬 Contrasts:", paste(contrasts, collapse = ", "), "\n")
    }
    
    # Store validation results
    validation_results[[filename]] <- list(
      status = "success",
      n_rows = nrow(data),
      n_cols = ncol(data),
      columns = names(data),
      required_present = present_cols,
      required_missing = missing_cols
    )
    
  }, error = function(e) {
    cat("   ❌ Error loading file:", e$message, "\n")
    validation_results[[filename]] <<- list(
      status = "error",
      error_message = e$message
    )
  })
}

cat("\n📋 Validation Summary:\n")
cat("=====================\n")

successful_files <- 0
failed_files <- 0

for (filename in names(validation_results)) {
  result <- validation_results[[filename]]
  if (result$status == "success") {
    successful_files <- successful_files + 1
    cat("✅", filename, "- VALID\n")
  } else {
    failed_files <- failed_files + 1
    cat("❌", filename, "- FAILED:", result$error_message, "\n")
  }
}

cat("\n📊 Overall Results:\n")
cat("Successful files:", successful_files, "/", length(validation_results), "\n")

# Test compatibility with TF analysis requirements
cat("\n🧬 TF Analysis Compatibility Check:\n")

if (successful_files > 0) {
  # Test with first successful file
  test_file <- NULL
  for (filename in names(validation_results)) {
    if (validation_results[[filename]]$status == "success") {
      test_file <- file.path(edger_path, filename)
      break
    }
  }
  
  if (!is.null(test_file)) {
    cat("Testing with:", basename(test_file), "\n")
    
    test_data <- read_csv(test_file, show_col_types = FALSE)
    
    # Check gene format compatibility
    gene_col <- "gene"
    if (gene_col %in% names(test_data)) {
      sample_genes <- head(test_data[[gene_col]], 10)
      
      # Check for ENSEMBL IDs
      has_ensembl <- any(str_detect(sample_genes, "^ENSG"))
      # Check for gene symbols
      has_symbols <- any(str_detect(sample_genes, "^[A-Z][A-Z0-9-]+$"))
      
      cat("   🔍 Gene ID format analysis:\n")
      cat("     ENSEMBL IDs detected:", has_ensembl, "\n")
      cat("     Gene symbols detected:", has_symbols, "\n")
      
      if (has_symbols || has_ensembl) {
        cat("   ✅ Gene identifiers compatible with DoRothEA\n")
      } else {
        cat("   ⚠️  Gene identifier format may need conversion\n")
      }
    }
    
    # Check statistical significance distribution
    if ("PValue" %in% names(test_data) && "logFC" %in% names(test_data)) {
      n_up <- sum(test_data$logFC > 0 & test_data$PValue < 0.05, na.rm = TRUE)
      n_down <- sum(test_data$logFC < 0 & test_data$PValue < 0.05, na.rm = TRUE)
      
      cat("   📈 Significant upregulated genes:", n_up, "\n")
      cat("   📉 Significant downregulated genes:", n_down, "\n")
      
      if (n_up > 10 && n_down > 10) {
        cat("   ✅ Good balance of up/down regulation for TF analysis\n")
      } else {
        cat("   ⚠️  Limited number of significant genes - TF analysis may be limited\n")
      }
    }
  }
}

# Final recommendations
cat("\n💡 Recommendations:\n")

if (successful_files == length(validation_results)) {
  cat("🎉 All files are valid! You're ready to run TF activity analysis.\n")
  cat("\n🚀 Next steps:\n")
  cat("1. Run: Rscript test_tf_setup.R\n")
  cat("2. Then: Rscript 05_tf_activity_decoupler.R\n")
} else {
  cat("⚠️  Some files have issues. Please check:\n")
  for (filename in names(validation_results)) {
    result <- validation_results[[filename]]
    if (result$status == "error") {
      cat("   -", filename, ":", result$error_message, "\n")
    }
  }
}

cat("\n📝 File Structure Expected by TF Analysis:\n")
cat("   - gene: Gene identifiers (symbols or ENSEMBL)\n")
cat("   - logFC: Log fold change values\n")
cat("   - PValue: Statistical p-values\n")
cat("   - FDR: False discovery rate (optional, will be calculated if missing)\n")

cat("\n✨ Validation complete!\n")