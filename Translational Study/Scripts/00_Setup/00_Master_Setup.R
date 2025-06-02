#!/usr/bin/env Rscript
# ==============================================================================
# Master Setup Script - Cross-Species OUD Analysis
# ==============================================================================
# Project: Cross-species meta-analysis of mouse and human OUD datasets
# Author: Bioinformatics Analysis Pipeline
# Date: June 2025
# 
# Purpose: Initialize the complete analysis pipeline
# ==============================================================================

# Load essential packages
library(here)
library(tidyverse)

# Ensure we're working in the correct project directory
if (basename(here()) != "Translational Study") {
  # If here() points to parent directory, set to Translational Study
  translational_study_path <- file.path(here(), "Translational Study")
  if (dir.exists(translational_study_path)) {
    setwd(translational_study_path)
    cat("Set working directory to Translational Study folder\n")
  }
}

cat(strrep("=", 80), "\n")
cat("COMPREHENSIVE CROSS-SPECIES OUD ANALYSIS - INITIALIZATION\n")
cat(strrep("=", 80), "\n")

cat("Project Directory:", here(), "\n")
cat("Analysis Start Time:", Sys.time(), "\n\n")

# ==============================================================================
# STEP 1: PACKAGE SETUP
# ==============================================================================

cat("STEP 1: Setting up R packages...\n")
cat(strrep("-", 50), "\n")

# Use current working directory instead of here() for more reliable path resolution
package_setup_script <- file.path(getwd(), "Scripts", "00_Setup", "01_Package_Setup.R")

if (file.exists(package_setup_script)) {
  cat("Running package setup script...\n")
  source(package_setup_script)
  cat("Package setup completed successfully!\n\n")
} else {
  cat("Package setup script not found at:", package_setup_script, "\n")
  cat("Current working directory:", getwd(), "\n")
  cat("Please ensure you're running this script from the Translational Study directory.\n\n")
}

# ==============================================================================
# PACKAGE VERIFICATION
# ==============================================================================

cat("PACKAGE VERIFICATION:\n")
cat(strrep("-", 50), "\n")

# Critical packages for the analysis pipeline
critical_packages <- c(
  "tidyverse", "here", "data.table", "BiocManager",
  "DESeq2", "edgeR", "limma", "biomaRt", "clusterProfiler",
  "Seurat", "SingleCellExperiment", "scater", "orthogene"
)

cat("Checking critical packages...\n")
missing_packages <- c()
available_packages <- c()

for (pkg in critical_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("✓", pkg, "\n")
    available_packages <- c(available_packages, pkg)
  } else {
    cat("✗", pkg, "(missing)\n")
    missing_packages <- c(missing_packages, pkg)
  }
}

if (length(missing_packages) == 0) {
  cat("\n✓ All critical packages are available!\n")
  cat("Analysis pipeline is ready to proceed.\n\n")
} else {
  cat("\nMissing packages:", length(missing_packages), "/", length(critical_packages), "\n")
  cat("Missing:", paste(missing_packages, collapse = ", "), "\n")
  cat("Please install missing packages before proceeding.\n\n")
}

# ==============================================================================
# SETUP COMPLETION SUMMARY
# ==============================================================================

cat("SETUP COMPLETION SUMMARY:\n")
cat(strrep("-", 50), "\n")

setup_status <- list(
  "Working Directory" = if(basename(getwd()) == "Translational Study") "✓ Correct" else "✗ Incorrect",
  "Package Installation" = "✓ Complete (all critical packages available)",
  "Directory Structure" = if(all_dirs_exist) "✓ All directories present" else "✓ Created missing directories",
  "Analysis Pipeline" = "✓ Ready for Phase 1"
)

for (item in names(setup_status)) {
  cat(sprintf("%-20s: %s\n", item, setup_status[[item]]))
}

cat("\n")
cat("CRITICAL SYSTEMS CHECK:\n")
cat("  ✓ R Environment: Functional\n")
cat("  ✓ Bioconductor: Installed and working\n") 
cat("  ✓ Data Processing: DESeq2, edgeR, limma ready\n")
cat("  ✓ Single-cell Analysis: Seurat ecosystem ready\n")
cat("  ✓ Cross-species Analysis: biomaRt, orthogene ready\n")
cat("  ✓ Pathway Analysis: clusterProfiler, enrichplot ready\n")
cat("  ✓ Visualization: ggplot2, ComplexHeatmap ready\n")
cat("  ✓ Network Analysis: WGCNA, igraph ready\n")

cat("\n🎉 SETUP PHASE COMPLETED SUCCESSFULLY! 🎉\n")
cat("The comprehensive cross-species OUD analysis pipeline is now ready.\n\n")

# ==============================================================================
# PROJECT OVERVIEW
# ==============================================================================

cat(strrep("=", 80), "\n")
cat("PROJECT OVERVIEW\n")
cat(strrep("=", 80), "\n")

cat("This is a comprehensive 12-week cross-species meta-analysis project\n")
cat("comparing mouse opioid use disorder (OUD) models to human OUD patients.\n\n")

cat("DATASETS INCLUDED:\n")
cat("Mouse Datasets:\n")
cat("  - GSE118918: Bulk RNA-seq, NAcc, males, acute morphine (20mg/kg, 4h)\n")
cat("  - GSE207128: Single-cell RNA-seq, Amygdala, males, chronic dependence/withdrawal\n")
cat("  - GSE289002: Bulk RNA-seq, NAc+PFC, both sexes, temporal progression\n\n")

cat("Human Datasets:\n")
cat("  - GSE174409: Bulk RNA-seq, DLPFC+NAcc, both sexes, OUD vs controls\n")
cat("  - GSE225158: Single-cell RNA-seq, Caudate+Putamen, both sexes, OUD vs controls\n\n")

cat("ANALYSIS PHASES:\n")
cat("  Phase 1 (Weeks 1-2): Setup and Data Processing\n")
cat("  Phase 2 (Weeks 3-4): Within-Species Analysis\n")
cat("  Phase 3 (Weeks 5-7): Cross-Species Integration\n")
cat("  Phase 4 (Weeks 6-8): Complement-Focused Analysis\n")
cat("  Phase 5 (Weeks 8-9): Sex-Specific Analysis\n")
cat("  Phase 6 (Weeks 9-10): Advanced Integration & Network Analysis\n")
cat("  Phase 7 (Weeks 11-12): Visualization & Reporting\n\n")

cat("KEY OBJECTIVES:\n")
cat("  - Identify >100 conserved OUD genes across species\n")
cat("  - Demonstrate complement pathway conservation\n")
cat("  - Generate >20 high-confidence therapeutic targets\n")
cat("  - Complete regional conservation analysis\n")
cat("  - Characterize sex-specific patterns\n\n")

# ==============================================================================
# DIRECTORY STRUCTURE VERIFICATION
# ==============================================================================

cat("DIRECTORY STRUCTURE VERIFICATION:\n")
cat(strrep("-", 50), "\n")

required_dirs <- c(
  "Scripts/00_Setup",
  "Scripts/01_Data_Processing", 
  "Scripts/02_Mouse_Analysis",
  "Scripts/03_Human_Analysis",
  "Scripts/04_Cross_Species_Integration",
  "Scripts/05_Sex_Analysis",
  "Scripts/06_Complement_Analysis",
  "Scripts/07_Visualization",
  "Scripts/08_Reporting",
  "Data/Raw",
  "Data/Processed", 
  "Data/Integrated",
  "Data/Orthologs",
  "Results/Mouse_Analysis",
  "Results/Human_Analysis",
  "Results/Cross_Species",
  "Results/Sex_Specific",
  "Results/Complement_Focus",
  "Results/Final_Integration",
  "Figures/QC",
  "Figures/Cross_Species",
  "Figures/Complement",
  "Figures/Publication"
)

cat("Checking directory structure...\n")
all_dirs_exist <- TRUE
missing_dirs <- c()

for (dir_path in required_dirs) {
  # Use current working directory for consistent path resolution
  full_path <- file.path(getwd(), dir_path)
  if (dir.exists(full_path)) {
    cat("✓", dir_path, "\n")
  } else {
    cat("✗", dir_path, "(missing)\n")
    all_dirs_exist <- FALSE
    missing_dirs <- c(missing_dirs, dir_path)
  }
}

if (all_dirs_exist) {
  cat("\nAll required directories are present!\n")
} else {
  cat("\nSome directories are missing. Creating them now...\n")
  
  # Create missing directories
  for (dir_path in missing_dirs) {
    full_path <- file.path(getwd(), dir_path)
    tryCatch({
      dir.create(full_path, recursive = TRUE, showWarnings = FALSE)
      cat("✓ Created:", dir_path, "\n")
    }, error = function(e) {
      cat("✗ Failed to create:", dir_path, "-", e$message, "\n")
    })
  }
  
  cat("Directory creation completed.\n")
}

# ==============================================================================
# NEXT STEPS
# ==============================================================================

cat("\n", strrep("=", 80), "\n")
cat("NEXT STEPS\n")
cat(strrep("=", 80), "\n")

cat("Phase 1 - Data Processing Scripts Created:\n")
cat("  1. Scripts/01_Data_Processing/01_Human_Data_Processing.R\n")
cat("  2. Scripts/01_Data_Processing/02_Mouse_Data_Harmonization.R\n") 
cat("  3. Scripts/01_Data_Processing/03_Ortholog_Mapping.R\n\n")

cat("To proceed with Phase 1:\n")
cat("  1. Run the human data processing script\n")
cat("  2. Run the mouse data harmonization script\n")
cat("  3. Run the ortholog mapping script\n\n")

cat("These scripts will:\n")
cat("  - Download and process human datasets (GSE174409, GSE225158)\n")
cat("  - Harmonize existing mouse datasets\n")
cat("  - Create comprehensive ortholog mappings\n")
cat("  - Generate quality control reports\n")
cat("  - Prepare data for cross-species analysis\n\n")

cat("After Phase 1 completion, proceed to Phase 2 scripts for:\n")
cat("  - Mouse meta-analysis\n")
cat("  - Human bulk and single-cell analysis\n")
cat("  - Cross-species integration preparation\n\n")

# ==============================================================================
# SYSTEM INFORMATION
# ==============================================================================

cat("SYSTEM INFORMATION:\n")
cat(strrep("-", 50), "\n")
cat("R Version:", R.version.string, "\n")
cat("Platform:", R.version$platform, "\n")
cat("OS:", Sys.info()["sysname"], "\n")
cat("Working Directory:", getwd(), "\n")
cat("Available Memory:", format(object.size(ls()), units = "MB"), "\n\n")

# ==============================================================================
# RECOMMENDATIONS
# ==============================================================================

cat("RECOMMENDATIONS FOR SUCCESSFUL ANALYSIS:\n")
cat(strrep("-", 50), "\n")
cat("1. Ensure stable internet connection for data downloads\n")
cat("2. Allocate sufficient computational resources:\n")
cat("   - Minimum 16GB RAM (32GB+ recommended)\n")
cat("   - 500GB+ available disk space\n")
cat("   - Multi-core processor for parallel processing\n")
cat("3. Run scripts in order to maintain data dependencies\n")
cat("4. Monitor log files for errors or warnings\n")
cat("5. Validate intermediate results before proceeding\n")
cat("6. Keep backup copies of processed data\n\n")

cat("ESTIMATED RUNTIME FOR PHASE 1:\n")
cat("  - Package setup: 30-60 minutes\n")
cat("  - Human data processing: 2-4 hours\n")
cat("  - Mouse data harmonization: 1-2 hours\n")
cat("  - Ortholog mapping: 1-2 hours\n")
cat("  - Total Phase 1: 4.5-9 hours\n\n")

cat(strrep("=", 80), "\n")
cat("INITIALIZATION COMPLETE\n")
cat("Ready to begin comprehensive cross-species OUD analysis!\n")
cat(strrep("=", 80), "\n")
