#===============================================================================
# DIFFERENTIAL EXPRESSION ANALYSIS
#===============================================================================
# 
# This script performs differential expression analysis using DESeq2 with 
# batch-corrected data from the preprocessing pipeline.
#
# WORKFLOW:
# 1. Load batch-corrected DESeq2 object
# 2. Define contrasts and perform DE analysis
# 3. Extract and filter results
# 4. Create visualization plots
# 5. Save results for pathway enrichment
#
#===============================================================================

# Load required libraries
library(DESeq2)
library(ggplot2)
library(dplyr)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)

# Set paths
setwd("/Users/aumchampaneri/Complement-OUD/Multi-Omics Study")
input_dir <- "data/processed/bulkrna/preprocessing"
output_dir <- "data/processed/bulkrna/differential_expression"
plots_dir <- "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study/results/bulkrna/differential_expression"

# Create directories
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

#===============================================================================
# 1. LOAD PREPROCESSED DATA
#===============================================================================

# Load batch-corrected DESeq2 object
load(file.path(input_dir, "deseq2_for_DE_analysis.RData"))

cat("Loaded DESeq2 object with", ncol(dds_batch_norm), "samples and", 
    nrow(dds_batch_norm), "genes\n")
cat("Sample conditions:\n")
print(table(sample_info$condition))
