#!/usr/bin/env Rscript
#' @title Test LEMUR Setup and Dependencies
#' @description Quick test script to validate R environment for LEMUR analysis
#' @author Generated Test Script
#' @date 2024

# Test script to validate LEMUR R implementation
cat("🧪 Testing LEMUR R Setup\n")
cat("========================\n\n")

# Test 1: Check R version
cat("📊 R Environment:\n")
cat("R version:", R.version.string, "\n")
cat("Platform:", R.version$platform, "\n\n")

# Test 2: Package availability check
cat("📦 Checking package availability...\n")

required_packages <- c(
  "SingleCellExperiment", "scran", "scater", "lemur",
  "tidyverse", "ggplot2", "viridis", "patchwork",
  "zellkonverter", "HDF5Array", "uwot", "harmony",
  "BiocNeighbors", "Matrix", "limma"
)

package_status <- list()
missing_packages <- c()

for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    version <- try(as.character(packageVersion(pkg)), silent = TRUE)
    if (inherits(version, "try-error")) {
      version <- "unknown"
    }
    package_status[[pkg]] <- paste("✅", version)
    cat("  ✅", pkg, "- version:", version, "\n")
  } else {
    package_status[[pkg]] <- "❌ Not installed"
    missing_packages <- c(missing_packages, pkg)
    cat("  ❌", pkg, "- NOT INSTALLED\n")
  }
}

cat("\n")

# Test 3: Check input file
cat("📁 Checking input data...\n")
input_file <- "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study/data/raw/snrna/GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5ad"

if (file.exists(input_file)) {
  file_info <- file.info(input_file)
  cat("  ✅ Input H5AD file found\n")
  cat("     Path:", input_file, "\n")
  cat("     Size:", round(file_info$size / 1024^3, 2), "GB\n")
  cat("     Modified:", as.character(file_info$mtime), "\n")
} else {
  cat("  ❌ Input H5AD file NOT FOUND\n")
  cat("     Expected:", input_file, "\n")
}

cat("\n")

# Test 4: Try loading key packages
cat("🔧 Testing key package functionality...\n")

# Test zellkonverter
if ("zellkonverter" %in% names(package_status) && !grepl("❌", package_status[["zellkonverter"]])) {
  tryCatch(
    {
      library(zellkonverter)
      cat("  ✅ zellkonverter loads successfully\n")
    },
    error = function(e) {
      cat("  ❌ zellkonverter loading error:", e$message, "\n")
    }
  )
} else {
  cat("  ⚠️  zellkonverter not available for testing\n")
}

# Test LEMUR
if ("lemur" %in% names(package_status) && !grepl("❌", package_status[["lemur"]])) {
  tryCatch(
    {
      library(lemur)
      cat("  ✅ LEMUR loads successfully\n")
    },
    error = function(e) {
      cat("  ❌ LEMUR loading error:", e$message, "\n")
    }
  )
} else {
  cat("  ⚠️  LEMUR not available for testing\n")
}

# Test SingleCellExperiment
if ("SingleCellExperiment" %in% names(package_status) && !grepl("❌", package_status[["SingleCellExperiment"]])) {
  tryCatch(
    {
      library(SingleCellExperiment)
      cat("  ✅ SingleCellExperiment loads successfully\n")
    },
    error = function(e) {
      cat("  ❌ SingleCellExperiment loading error:", e$message, "\n")
    }
  )
} else {
  cat("  ⚠️  SingleCellExperiment not available for testing\n")
}

cat("\n")

# Test 5: Memory and system info
cat("💾 System Information:\n")
if (Sys.info()["sysname"] == "Darwin") {
  cat("  OS: macOS\n")
} else if (Sys.info()["sysname"] == "Linux") {
  cat("  OS: Linux\n")
} else {
  cat("  OS:", Sys.info()["sysname"], "\n")
}

# Get memory info (if available)
tryCatch(
  {
    if (requireNamespace("pryr", quietly = TRUE)) {
      library(pryr)
      cat("  Available memory:", pryr::mem_used(), "\n")
    } else {
      cat("  Memory info: pryr package not available\n")
    }
  },
  error = function(e) {
    cat("  Memory info: unable to determine\n")
  }
)

cat("\n")

# Summary and recommendations
cat("📋 SUMMARY\n")
cat("==========\n")

if (length(missing_packages) == 0) {
  cat("🎉 All required packages are installed!\n")
  cat("✅ Ready to run LEMUR analysis\n\n")

  cat("🚀 To run the analysis:\n")
  cat("   ./run_lemur.sh\n")
  cat("   or\n")
  cat("   Rscript run_LEMUR_GSE225158.R\n")
} else {
  cat("⚠️  Missing packages detected\n")
  cat("❌ Missing:", length(missing_packages), "packages\n")
  cat("   -", paste(missing_packages, collapse = "\n   - "), "\n\n")

  cat("🔧 To install missing packages:\n")
  cat("   ./run_lemur.sh --install-packages --skip-analysis\n")
  cat("   or run the following R commands:\n\n")

  # Generate installation commands
  bioc_packages <- c(
    "SingleCellExperiment", "scran", "scater", "lemur",
    "zellkonverter", "HDF5Array", "BiocNeighbors", "limma"
  )
  cran_packages <- setdiff(missing_packages, bioc_packages)
  bioc_missing <- intersect(missing_packages, bioc_packages)

  if (length(cran_packages) > 0) {
    cat("   # CRAN packages\n")
    cat("   install.packages(c(", paste0('"', cran_packages, '"', collapse = ", "), "))\n\n")
  }

  if (length(bioc_missing) > 0) {
    cat("   # Bioconductor packages\n")
    cat("   if (!requireNamespace('BiocManager', quietly = TRUE))\n")
    cat("       install.packages('BiocManager')\n")
    cat("   BiocManager::install(c(", paste0('"', bioc_missing, '"', collapse = ", "), "))\n")
  }
}

if (!file.exists(input_file)) {
  cat("\n⚠️  Input file issue:\n")
  cat("   The expected H5AD file was not found.\n")
  cat("   Please ensure the data file is available at:\n")
  cat("   ", input_file, "\n")
}

cat("\n")
cat("📞 For troubleshooting, see README_LEMUR.md\n")
cat("🧪 Test completed at:", as.character(Sys.time()), "\n")
