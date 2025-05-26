# GSE207128 Project Dependencies
# Install all required packages for the analysis pipeline

# Core packages
install.packages(c(
  "Seurat",
  "dplyr",
  "ggplot2",
  "patchwork",
  "Matrix"
))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "SingleR",
  "celldex"
))

# Integration packages (optional - will fallback if not available)
install.packages("harmony")  # Primary integration method

# Print session info
sessionInfo()
