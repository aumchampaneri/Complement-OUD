#!/usr/bin/env Rscript
# ==============================================================================
# Comprehensive Cross-Species OUD Analysis - Package Setup
# ==============================================================================
# Project: Cross-species meta-analysis of mouse and human OUD datasets
# Author: Bioinformatics Analysis Pipeline
# Date: June 2025
# 
# Purpose: Install and load all required packages for the comprehensive analysis
# ==============================================================================

# Set options for better package installation
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 300)  # Increase timeout for large downloads

# Function to install packages if not already installed
install_if_missing <- function(packages, repo = "CRAN") {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat("Installing package:", pkg, "\n")
      if (repo == "CRAN") {
        install.packages(pkg, dependencies = TRUE)
      } else if (repo == "Bioconductor") {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
        }
        BiocManager::install(pkg, dependencies = TRUE, update = FALSE)
      }
      library(pkg, character.only = TRUE)
    } else {
      cat("Package", pkg, "already installed and loaded\n")
    }
  }
}

# ==============================================================================
# CORE DATA MANIPULATION PACKAGES
# ==============================================================================
cat("Installing core data manipulation packages...\n")
core_packages <- c(
  "tidyverse",      # Data manipulation and visualization
  "data.table",     # Fast data manipulation
  "readxl",         # Excel file reading
  "openxlsx",       # Excel file writing
  "here",           # Reproducible file paths
  "glue",           # String interpolation
  "janitor",        # Data cleaning
  "scales",         # Scaling functions for ggplot2
  "RColorBrewer",   # Color palettes
  "viridis"         # Perceptually uniform color palettes
)

install_if_missing(core_packages, "CRAN")

# ==============================================================================
# BIOINFORMATICS CORE PACKAGES
# ==============================================================================
cat("Installing core bioinformatics packages...\n")
bioc_core <- c(
  "BiocManager",           # Bioconductor package manager
  "DESeq2",               # Differential expression analysis
  "edgeR",                # Differential expression analysis
  "limma",                # Linear models for microarray/RNA-seq
  "biomaRt",              # Access to Ensembl databases
  "AnnotationDbi",        # Annotation database interface
  "org.Hs.eg.db",         # Human gene annotations
  "org.Mm.eg.db",         # Mouse gene annotations
  "GenomicFeatures",      # Genomic features manipulation
  "rtracklayer",          # Interface to genome browsers
  "Biostrings",           # String objects for biological sequences
  "IRanges",              # Infrastructure for ranges
  "GenomicRanges",        # Genomic ranges manipulation
  "SummarizedExperiment"  # Container for experiment data
)

install_if_missing(bioc_core, "Bioconductor")

# ==============================================================================
# SINGLE-CELL ANALYSIS PACKAGES
# ==============================================================================
cat("Installing single-cell analysis packages...\n")
singlecell_packages <- c(
  "Seurat",                    # Single-cell analysis
  "SingleCellExperiment",      # Single-cell data container
  "scater",                    # Single-cell analysis toolkit
  "scran",                     # Single-cell RNA-seq analysis
  "celldex",                   # Cell type reference datasets
  "SingleR",                   # Cell type annotation
  "scuttle",                   # Single-cell utilities
  "bluster",                   # Clustering algorithms
  "batchelor"                  # Batch correction
)

install_if_missing(singlecell_packages, "Bioconductor")

# ==============================================================================
# CROSS-SPECIES ANALYSIS PACKAGES
# ==============================================================================
cat("Installing cross-species analysis packages...\n")
cross_species_packages <- c(
  "orthogene",        # Ortholog mapping
  "homologene",       # Homolog gene mapping
  "gprofiler2"        # Functional enrichment analysis
)

install_if_missing(cross_species_packages, "CRAN")

# ==============================================================================
# PATHWAY ANALYSIS PACKAGES
# ==============================================================================
cat("Installing pathway analysis packages...\n")
pathway_packages <- c(
  "clusterProfiler",    # Statistical analysis of functional profiles
  "enrichplot",         # Visualization of enrichment results
  "pathview",           # Pathway visualization
  "ReactomePA",         # Reactome pathway analysis
  "DOSE",               # Disease ontology semantic analysis
  "msigdbr",            # MSigDB gene sets
  "fgsea",              # Fast gene set enrichment analysis
  "GSEABase"            # Gene set enrichment analysis base
)

install_if_missing(pathway_packages, "Bioconductor")

# ==============================================================================
# NETWORK ANALYSIS PACKAGES
# ==============================================================================
cat("Installing network analysis packages...\n")
network_packages <- c(
  "WGCNA",          # Weighted gene co-expression network analysis
  "igraph",         # Network analysis and visualization
  "corrplot",       # Correlation plot visualization
  "ggraph",         # Grammar of graphics for networks
  "tidygraph",      # Tidy API for graph manipulation
  "networkD3"       # Interactive network visualizations
)

install_if_missing(network_packages, "CRAN")

# ==============================================================================
# VISUALIZATION PACKAGES
# ==============================================================================
cat("Installing visualization packages...\n")
viz_packages_bioc <- c(
  "ComplexHeatmap",     # Complex heatmap visualization
  "EnhancedVolcano",    # Enhanced volcano plots
  "regioneR"            # Association analysis of genomic regions
)

viz_packages_cran <- c(
  "pheatmap",           # Pretty heatmaps
  "VennDiagram",        # Venn diagram creation
  "UpSetR",             # UpSet plots for set intersections
  "ggsci",              # Scientific journal color palettes
  "ggpubr",             # Publication ready plots
  "cowplot",            # Plot arrangements
  "patchwork",          # Combine ggplot2 plots
  "plotly",             # Interactive plots
  "DT",                 # Interactive data tables
  "kableExtra",         # Enhanced tables
  "gridExtra",          # Grid arrangements
  "ggnewscale",         # Multiple scales in ggplot2
  "ggrepel",            # Non-overlapping text labels
  "circlize"            # Circular visualization
)

install_if_missing(viz_packages_bioc, "Bioconductor")
install_if_missing(viz_packages_cran, "CRAN")

# ==============================================================================
# QUALITY CONTROL AND BATCH EFFECTS PACKAGES
# ==============================================================================
cat("Installing quality control and batch effect packages...\n")
qc_packages <- c(
  "sva",              # Surrogate variable analysis
  "RUVSeq",           # Remove unwanted variation
  "limma",            # Already installed above
  "ComBatSeq"         # Batch effect correction for RNA-seq
)

# ComBatSeq is part of sva, so we don't need separate installation
qc_packages_to_install <- c("sva", "RUVSeq")
install_if_missing(qc_packages_to_install, "Bioconductor")

# Harmony for single-cell batch correction
install_if_missing("harmony", "CRAN")

# ==============================================================================
# META-ANALYSIS PACKAGES
# ==============================================================================
cat("Installing meta-analysis packages...\n")
meta_packages <- c(
  "metafor",          # Meta-analysis framework
  "RankProd",         # Rank product method
  "MetaVolcanoR",     # Meta-analysis visualization
  "meta"              # General meta-analysis
)

# RankProd is from Bioconductor
install_if_missing("RankProd", "Bioconductor")
install_if_missing(c("metafor", "meta"), "CRAN")

# MetaVolcanoR might need to be installed from GitHub
if (!require("MetaVolcanoR", character.only = TRUE, quietly = TRUE)) {
  if (!require("devtools", character.only = TRUE, quietly = TRUE)) {
    install.packages("devtools")
  }
  devtools::install_github("hyunyulhenry/MetaVolcanoR")
  library(MetaVolcanoR)
}

# ==============================================================================
# STATISTICAL ANALYSIS PACKAGES
# ==============================================================================
cat("Installing additional statistical packages...\n")
stats_packages <- c(
  "broom",            # Tidy statistical output
  "car",              # Companion to applied regression
  "emmeans",          # Estimated marginal means
  "multcomp",         # Multiple comparisons
  "coin",             # Conditional inference
  "exactRankTests",   # Exact rank tests
  "Hmisc",            # Harrell miscellaneous
  "psych",            # Psychological/psychometric functions
  "effsize",          # Effect size calculation
  "pwr"               # Power analysis
)

install_if_missing(stats_packages, "CRAN")

# ==============================================================================
# REPORT GENERATION PACKAGES
# ==============================================================================
cat("Installing report generation packages...\n")
report_packages <- c(
  "rmarkdown",        # Dynamic documents
  "knitr",            # Elegant, flexible and fast dynamic report generation
  "bookdown",         # Authoring books with R Markdown
  "flexdashboard",    # Interactive dashboards
  "crosstalk",        # Inter-widget interactivity
  "htmlwidgets"       # HTML widgets for R
)

install_if_missing(report_packages, "CRAN")

# ==============================================================================
# PARALLEL PROCESSING PACKAGES
# ==============================================================================
cat("Installing parallel processing packages...\n")
parallel_packages <- c(
  "parallel",         # Base R parallel processing
  "doParallel",       # Parallel backend for foreach
  "foreach",          # Parallel for loops
  "future",           # Unified parallel processing
  "future.apply",     # Apply functions in parallel
  "BiocParallel"      # Bioconductor parallel processing
)

install_if_missing(c("doParallel", "foreach", "future", "future.apply"), "CRAN")
install_if_missing("BiocParallel", "Bioconductor")

# ==============================================================================
# DATABASE ACCESS PACKAGES
# ==============================================================================
cat("Installing database access packages...\n")
db_packages <- c(
  "GEOquery",         # Access to GEO datasets
  "ArrayExpress",     # Access to ArrayExpress
  "rentrez",          # Interface to NCBI databases
  "annotate"          # Annotation for microarrays
)

install_if_missing(db_packages, "Bioconductor")

# ==============================================================================
# SPECIALIZED ANALYSIS PACKAGES
# ==============================================================================
cat("Installing specialized analysis packages...\n")

# For gene set variation analysis
install_if_missing("GSVA", "Bioconductor")

# For time series analysis (for temporal progression in GSE289002)
install_if_missing(c("forecast", "tseries", "changepoint"), "CRAN")

# For machine learning approaches
install_if_missing(c("randomForest", "glmnet", "caret", "e1071"), "CRAN")

# For survival analysis (if needed for clinical data)
install_if_missing("survival", "CRAN")

# ==============================================================================
# COMPLEMENT-SPECIFIC GENE SETS
# ==============================================================================
cat("Loading complement-specific gene sets...\n")

# We'll define these in a separate script, but load msigdbr for access to gene sets
if (require("msigdbr", character.only = TRUE, quietly = TRUE)) {
  cat("MSigDB gene sets available for complement pathway analysis\n")
}

# ==============================================================================
# SESSION INFO AND PACKAGE VERSIONS
# ==============================================================================
cat("Generating session information...\n")

# Create a session info report
session_info <- sessionInfo()
cat("Session Info saved. Key package versions:\n")
cat("R version:", as.character(session_info$R.version$version.string), "\n")

# Save session info
saveRDS(session_info, file = here::here("Scripts", "00_Setup", "session_info.rds"))

# Create a package version summary
package_versions <- data.frame(
  Package = names(sessionInfo()$otherPkgs),
  Version = sapply(sessionInfo()$otherPkgs, function(x) x$Version),
  stringsAsFactors = FALSE
)

write.csv(package_versions, 
          file = here::here("Scripts", "00_Setup", "package_versions.csv"),
          row.names = FALSE)

# ==============================================================================
# LOAD COMMONLY USED PACKAGES
# ==============================================================================
cat("Loading commonly used packages for analysis...\n")

# Load most frequently used packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(here)
  library(DESeq2)
  library(Seurat)
  library(biomaRt)
  library(clusterProfiler)
  library(ComplexHeatmap)
  library(ggplot2)
  library(viridis)
})

cat("Package setup complete!\n")
cat("Total packages installed and ready for cross-species OUD analysis.\n")
cat("Session info and package versions saved to Scripts/00_Setup/\n")

# ==============================================================================
# CLEANUP
# ==============================================================================
# Clean up installation messages
invisible(gc())

cat("\n===============================================\n")
cat("PACKAGE SETUP COMPLETED SUCCESSFULLY\n")
cat("===============================================\n")
cat("You can now proceed with the analysis pipeline.\n")
cat("All required packages are installed and loaded.\n")
