#!/usr/bin/env Rscript

# =============================================================================
# Test Script for TF Activity Analysis Setup
# =============================================================================
# 
# Description: Quick test to verify all required packages and dependencies
#              are properly installed for TF activity analysis
# 
# Author: Multi-Omics Study Team
# Date: 2024
# 
# Usage: Rscript test_tf_setup.R
# =============================================================================

cat("🧬 Testing TF Activity Analysis Setup...\n")
cat("=====================================\n\n")

# Test package installations
required_packages <- c(
  "decoupleR", "dplyr", "tidyr", "readr", "ggplot2", "ComplexHeatmap", 
  "circlize", "pheatmap", "viridis", "RColorBrewer", "reshape2",
  "broom", "purrr", "tibble", "stringr", "corrplot", "igraph",
  "networkD3", "plotly", "DT", "htmlwidgets", "org.Hs.eg.db",
  "ggrepel", "rlang"
)

cat("📦 Checking package installations...\n")
missing_packages <- c()
installed_packages <- c()

for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    installed_packages <- c(installed_packages, pkg)
    cat("✅", pkg, "\n")
  } else {
    missing_packages <- c(missing_packages, pkg)
    cat("❌", pkg, "- MISSING\n")
  }
}

cat("\n📊 Installation Summary:\n")
cat("Installed:", length(installed_packages), "/", length(required_packages), "packages\n")

if (length(missing_packages) > 0) {
  cat("\n⚠️  Missing packages:\n")
  for (pkg in missing_packages) {
    cat("  -", pkg, "\n")
  }
  cat("\n💡 To install missing packages, run:\n")
  
  # Separate Bioconductor and CRAN packages
  bioc_packages <- c("decoupleR", "org.Hs.eg.db", "ComplexHeatmap")
  missing_bioc <- intersect(missing_packages, bioc_packages)
  missing_cran <- setdiff(missing_packages, bioc_packages)
  
  if (length(missing_bioc) > 0) {
    cat("BiocManager::install(c(", paste0('"', missing_bioc, '"', collapse = ", "), "))\n")
  }
  if (length(missing_cran) > 0) {
    cat("install.packages(c(", paste0('"', missing_cran, '"', collapse = ", "), "))\n")
  }
}

# Test core functionality
cat("\n🧪 Testing core functionality...\n")

# Test 1: Load core packages
test_packages <- c("decoupleR", "dplyr", "ggplot2")
all_core_loaded <- TRUE

for (pkg in test_packages) {
  tryCatch({
    library(pkg, character.only = TRUE, quietly = TRUE)
    cat("✅ Loaded", pkg, "\n")
  }, error = function(e) {
    cat("❌ Failed to load", pkg, ":", e$message, "\n")
    all_core_loaded <<- FALSE
  })
}

# Test 2: Check decoupleR functionality
if ("decoupleR" %in% installed_packages) {
  cat("\n🔬 Testing decoupleR functionality...\n")
  
  tryCatch({
    # Test DoRothEA access
    dorothea_test <- decoupleR::get_dorothea(organism = "human", levels = c("A"))
    cat("✅ DoRothEA database accessible (", nrow(dorothea_test), "interactions)\n")
    
    # Test basic decoupleR functions
    test_matrix <- matrix(rnorm(100), nrow = 10, ncol = 10)
    rownames(test_matrix) <- paste0("Gene", 1:10)
    colnames(test_matrix) <- paste0("Sample", 1:10)
    
    cat("✅ Test matrix creation successful\n")
    
  }, error = function(e) {
    cat("❌ decoupleR test failed:", e$message, "\n")
    all_core_loaded <- FALSE
  })
}

# Test 3: Check file system access
cat("\n📁 Testing file system access...\n")

# Test directory creation
test_dir <- "test_output"
tryCatch({
  dir.create(test_dir, showWarnings = FALSE)
  if (dir.exists(test_dir)) {
    cat("✅ Directory creation successful\n")
    unlink(test_dir, recursive = TRUE)  # Clean up
  } else {
    cat("❌ Directory creation failed\n")
  }
}, error = function(e) {
  cat("❌ File system test failed:", e$message, "\n")
})

# Test 4: Check for sample data files
cat("\n📄 Checking for input data...\n")

# Look for potential edgeR results files
search_paths <- c(
  "../../results/snrna_scvi/differential_expression_edgeR",
  "results/snrna_scvi/differential_expression_edgeR",
  "results/differential_expression",
  "results/edgeR", 
  "outputs/differential_expression",
  "outputs/edgeR",
  "data/edgeR_results",
  "."
)

found_files <- c()
for (path in search_paths) {
  if (dir.exists(path)) {
    files <- list.files(path, pattern = "(edgeR.*results|DE|differential).*\\.(csv|tsv|txt)$", 
                       full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
    # Filter out summary files and focus on individual contrast results
    files <- files[!grepl("(summary|all_contrasts|significant)", basename(files))]
    if (length(files) > 0) {
      found_files <- c(found_files, files)
      cat("✅ Found", length(files), "potential input files in", path, "\n")
    }
  }
}

if (length(found_files) == 0) {
  cat("⚠️  No edgeR results files found in standard locations\n")
  cat("   You may need to update file paths in the main script\n")
} else {
  cat("✅ Total potential input files found:", length(found_files), "\n")
}

# Final summary
cat("\n🎯 Setup Summary:\n")
cat("================\n")

if (length(missing_packages) == 0 && all_core_loaded) {
  cat("🎉 All systems ready! You can run the TF activity analysis.\n")
  cat("\n🚀 To start analysis, run:\n")
  cat("   Rscript 05_tf_activity_decoupler.R\n")
} else {
  cat("⚠️  Setup incomplete. Please install missing packages and resolve issues.\n")
}

cat("\n📋 Next steps:\n")
cat("1. Ensure all required packages are installed\n")
cat("2. Verify your edgeR results files are in the correct location\n")
cat("3. Update file paths in 05_tf_activity_decoupler.R if needed\n")
cat("4. Run the main TF activity analysis script\n")

cat("\n💡 For help:\n")
cat("- Check package documentation: help(package = 'decoupleR')\n")
cat("- Verify file paths match your data location\n")
cat("- Ensure edgeR results have required columns: logFC, PValue, gene_symbol\n")

cat("\n✨ Setup test complete!\n")