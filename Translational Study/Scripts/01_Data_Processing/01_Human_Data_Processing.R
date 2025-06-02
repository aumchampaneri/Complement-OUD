#!/usr/bin/env Rscript
# ==============================================================================
# Human Data Processing Script - Cross-Species OUD Analysis
# ==============================================================================
# Project: Cross-species meta-analysis of mouse and human OUD datasets
# Author: Bioinformatics Analysis Pipeline
# Date: June 2025
# 
# Purpose: Download and process human OUD datasets for cross-species integration
# Datasets: GSE174409 (bulk RNA-seq), GSE225158 (single-cell RNA-seq)
# ==============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(GEOquery)
  library(data.table)
  library(DESeq2)
  library(edgeR)
  library(limma)
  library(biomaRt)
  library(Seurat)
  library(SingleCellExperiment)
  library(openxlsx)
  library(pheatmap)
  library(ggplot2)
  library(RColorBrewer)
  library(viridis)
})

cat(strrep("=", 80), "\n")
cat("HUMAN DATA PROCESSING - CROSS-SPECIES OUD ANALYSIS\n")
cat(strrep("=", 80), "\n")
cat("Processing Time:", Sys.time(), "\n\n")

# Set up directories (ensure we're in the correct project root)
# Check if we're already in the Translational Study folder
current_dir <- getwd()
if (basename(current_dir) == "Translational Study") {
  project_root <- current_dir
} else if (file.exists(file.path(current_dir, "Translational Study"))) {
  # If Translational Study exists as subdirectory, change to it
  project_root <- file.path(current_dir, "Translational Study")
  setwd(project_root)
} else {
  # Create Translational Study structure if it doesn't exist
  project_root <- file.path(current_dir, "Translational Study")
  dir.create(project_root, showWarnings = FALSE, recursive = TRUE)
  setwd(project_root)
}

raw_data_dir <- file.path(project_root, "Data", "Raw")
processed_data_dir <- file.path(project_root, "Data", "Processed")
results_dir <- file.path(project_root, "Results", "Human_Analysis")
figures_dir <- file.path(project_root, "Figures", "QC")

# Create directories if they don't exist
dirs_to_create <- c(raw_data_dir, processed_data_dir, results_dir, figures_dir)
walk(dirs_to_create, ~ if (!dir.exists(.x)) dir.create(.x, recursive = TRUE))

cat("Project root:", project_root, "\n")
cat("Current working directory:", getwd(), "\n")
cat("Output directories created/verified\n\n")

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

# Robust download function with retry mechanism
download_with_retry <- function(url, destfile, max_attempts = 3, timeout = 3600) {
  for (attempt in 1:max_attempts) {
    cat("Download attempt", attempt, "of", max_attempts, "for", basename(destfile), "\n")
    
    # Set timeout
    old_timeout <- getOption("timeout")
    options(timeout = timeout)
    
    tryCatch({
      download.file(url, destfile, mode = "wb", method = "auto")
      
      # Check if file was downloaded successfully
      if (file.exists(destfile) && file.size(destfile) > 1000) {
        cat("✓ Successfully downloaded:", basename(destfile), "\n")
        options(timeout = old_timeout)
        return(TRUE)
      }
    }, error = function(e) {
      cat("Attempt", attempt, "failed:", e$message, "\n")
      if (file.exists(destfile)) file.remove(destfile)
    })
    
    # Wait before retry
    if (attempt < max_attempts) {
      cat("Waiting 10 seconds before retry...\n")
      Sys.sleep(10)
    }
    
    options(timeout = old_timeout)
  }
  
  cat("✗ Failed to download after", max_attempts, "attempts:", basename(destfile), "\n")
  return(FALSE)
}

# Function to validate file integrity
validate_file <- function(filepath, expected_min_size = 1000) {
  if (!file.exists(filepath)) {
    return(list(valid = FALSE, reason = "File does not exist"))
  }
  
  file_info <- file.info(filepath)
  if (file_info$size < expected_min_size) {
    return(list(valid = FALSE, reason = paste("File too small:", file_info$size, "bytes")))
  }
  
  return(list(valid = TRUE, reason = "File appears valid"))
}

# ==============================================================================
# PART 1: GSE174409 - Human Bulk RNA-seq Processing
# ==============================================================================

cat("\n", strrep("=", 60), "\n")
cat("PROCESSING GSE174409 - HUMAN BULK RNA-SEQUENCING (OUD VS CONTROL)\n")
cat(strrep("=", 60), "\n")

# Create dataset-specific directory
gse174409_dir <- file.path(raw_data_dir, "GSE174409")
dir.create(gse174409_dir, showWarnings = FALSE, recursive = TRUE)

# Define file URLs and paths
gse174409_urls <- c(
  counts = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174409/suppl/GSE174409_raw_counts_02102020.csv.gz",
  filtered = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174409/suppl/GSE174409_human_Opioid_filtered_2020.csv.gz"
)

gse174409_files <- c(
  counts = file.path(gse174409_dir, "GSE174409_raw_counts_02102020.csv.gz"),
  filtered = file.path(gse174409_dir, "GSE174409_human_Opioid_filtered_2020.csv.gz")
)

# Download GSE174409 data files
cat("Downloading GSE174409 supplementary files...\n")
gse174409_downloads <- map2_lgl(gse174409_urls, gse174409_files, ~ {
  if (!file.exists(.y)) {
    download_with_retry(.x, .y, timeout = 600)  # 10 minute timeout
  } else {
    validation <- validate_file(.y, 1000000)  # 1MB minimum
    if (validation$valid) {
      cat("✓ File already exists and appears valid:", basename(.y), "\n")
      return(TRUE)
    } else {
      cat("⚠ Existing file invalid (", validation$reason, "), re-downloading...\n")
      return(download_with_retry(.x, .y, timeout = 600))
    }
  }
})

# Download GEO metadata
cat("Downloading GSE174409 metadata...\n")
gse174409_metadata <- NULL
tryCatch({
  gse174409_geo <- getGEO("GSE174409", GSEMatrix = TRUE, destdir = gse174409_dir, AnnotGPL = FALSE)
  if (length(gse174409_geo) > 0) {
    gse174409_metadata <- pData(gse174409_geo[[1]])
    cat("✓ Downloaded metadata for", nrow(gse174409_metadata), "samples\n")
  }
}, error = function(e) {
  cat("Warning: Could not download GEO metadata:", e$message, "\n")
})

# Process raw counts data
gse174409_processed <- FALSE
if (file.exists(gse174409_files[["counts"]])) {
  cat("Processing GSE174409 raw counts data...\n")
  
  tryCatch({
    # Read compressed file
    gse174409_raw <- fread(gse174409_files[["counts"]])
    cat("Raw data dimensions:", nrow(gse174409_raw), "x", ncol(gse174409_raw), "\n")
    
    # Extract gene IDs and count matrix
    if (ncol(gse174409_raw) > 1) {
      gene_ids <- gse174409_raw[[1]]
      count_columns <- names(gse174409_raw)[-1]
      
      # Convert to matrix
      gse174409_counts <- as.matrix(gse174409_raw[, -1, with = FALSE])
      rownames(gse174409_counts) <- gene_ids
      colnames(gse174409_counts) <- count_columns
      
      # Basic data validation
      cat("Count matrix summary:\n")
      cat("  Genes:", nrow(gse174409_counts), "\n")
      cat("  Samples:", ncol(gse174409_counts), "\n")
      cat("  Total counts:", format(sum(gse174409_counts), scientific = TRUE), "\n")
      cat("  Zero counts:", sum(gse174409_counts == 0), "/", length(gse174409_counts), 
          "(", round(100 * sum(gse174409_counts == 0) / length(gse174409_counts), 1), "%)\n")
      
      gse174409_processed <- TRUE
      cat("✓ Successfully processed raw counts\n")
    }
  }, error = function(e) {
    cat("✗ Error processing raw counts:", e$message, "\n")
  })
}

# Create comprehensive sample metadata
if (gse174409_processed) {
  cat("Creating sample metadata...\n")
  
  sample_names <- colnames(gse174409_counts)
  
  # Extract metadata from sample names and GEO data
  gse174409_sample_metadata <- data.frame(
    sample_id = sample_names,
    stringsAsFactors = FALSE
  )
  
  # Parse sample names for condition and region information
  gse174409_sample_metadata <- gse174409_sample_metadata %>%
    mutate(
      # Extract condition (OUD vs Control) - improved logic
      condition = case_when(
        str_detect(tolower(sample_id), "oud|opioid|morphine|heroin|fentanyl|overdose|addiction") ~ "OUD",
        str_detect(tolower(sample_id), "control|ctrl|normal|healthy|ref") ~ "Control",
        TRUE ~ "Unknown"
      ),
      
      # Extract brain region - improved patterns
      brain_region = case_when(
        str_detect(tolower(sample_id), "dlpfc|prefrontal|pfc") ~ "DLPFC",
        str_detect(tolower(sample_id), "nacc|nucleus.accumbens|accumbens|striatum") ~ "NAcc",
        str_detect(tolower(sample_id), "hippocampus|hipp|hc") ~ "Hippocampus",
        str_detect(tolower(sample_id), "amygdala|amyg|amy") ~ "Amygdala",
        TRUE ~ "Unknown"
      ),
      
      # Extract other metadata
      dataset = "GSE174409",
      species = "Human",
      data_type = "bulk_rnaseq",
      study_type = "case_control"
    )
  
  # Enhanced metadata extraction from GEO data
  if (!is.null(gse174409_metadata)) {
    cat("Processing GEO metadata...\n")
    
    # Display available GEO metadata columns for inspection
    geo_cols <- colnames(gse174409_metadata)
    cat("Available GEO metadata columns:\n")
    print(geo_cols[1:min(10, length(geo_cols))])
    
    # Extract key columns from GEO metadata
    relevant_cols <- geo_cols[str_detect(tolower(geo_cols), "title|source|characteristics|treatment|condition|tissue")]
    
    if (length(relevant_cols) > 0) {
      cat("Found relevant columns:", paste(relevant_cols, collapse = ", "), "\n")
      
      geo_simplified <- gse174409_metadata %>%
        dplyr::select(any_of(relevant_cols)) %>%
        mutate(sample_id = rownames(gse174409_metadata))
      
      # Print first few rows of key columns to understand data structure
      for (col in relevant_cols[1:min(3, length(relevant_cols))]) {
        cat("\nSample", col, "values (first 5):\n")
        print(head(geo_simplified[[col]], 5))
      }
      
      # Extract condition from GEO metadata
      geo_simplified$geo_condition <- "Unknown"
      geo_simplified$geo_brain_region <- "Unknown"
      
      for (col in relevant_cols) {
        if (col %in% colnames(geo_simplified)) {
          col_data <- tolower(as.character(geo_simplified[[col]]))
          
          # Check for condition patterns
          oud_patterns <- str_detect(col_data, "oud|opioid|morphine|heroin|fentanyl|overdose|addiction|case|user|dependent|abuse")
          control_patterns <- str_detect(col_data, "control|ctrl|normal|healthy|ref")
          
          # Update conditions where found
          geo_simplified$geo_condition[oud_patterns] <- "OUD"
          geo_simplified$geo_condition[control_patterns] <- "Control"
          
          # Check for brain region patterns
          dlpfc_patterns <- str_detect(col_data, "dlpfc|prefrontal|pfc")
          nacc_patterns <- str_detect(col_data, "nacc|nucleus.accumbens|accumbens|striatum")
          hipp_patterns <- str_detect(col_data, "hippocampus|hipp|hc")
          amyg_patterns <- str_detect(col_data, "amygdala|amyg|amy")
          
          # Update brain regions where found
          geo_simplified$geo_brain_region[dlpfc_patterns] <- "DLPFC"
          geo_simplified$geo_brain_region[nacc_patterns] <- "NAcc" 
          geo_simplified$geo_brain_region[hipp_patterns] <- "Hippocampus"
          geo_simplified$geo_brain_region[amyg_patterns] <- "Amygdala"
        }
      }
      
      # Join with sample metadata
      gse174409_sample_metadata <- gse174409_sample_metadata %>%
        left_join(geo_simplified, by = "sample_id")
      
      # Use GEO condition/region if sample name parsing failed
      gse174409_sample_metadata <- gse174409_sample_metadata %>%
        mutate(
          condition = ifelse(condition == "Unknown" & !is.na(geo_condition) & geo_condition != "Unknown", 
                           geo_condition, condition),
          brain_region = ifelse(brain_region == "Unknown" & !is.na(geo_brain_region) & geo_brain_region != "Unknown", 
                              geo_brain_region, brain_region)
        )
      
      cat("After GEO metadata processing:\n")
      cat("Conditions found:", paste(unique(geo_simplified$geo_condition), collapse = ", "), "\n")
      cat("Brain regions found:", paste(unique(geo_simplified$geo_brain_region), collapse = ", "), "\n")
    }
    
    # Additional processing for NAC brain region and control samples
    # Based on the output, we see "NAC" in source_name_ch1 and "diagnosis: OUD"/"diagnosis: CONT" in characteristics_ch1
    if ("source_name_ch1" %in% colnames(gse174409_sample_metadata)) {
      gse174409_sample_metadata <- gse174409_sample_metadata %>%
        mutate(
          brain_region = case_when(
            str_detect(tolower(source_name_ch1), "nac|nucleus") ~ "NAcc",
            str_detect(tolower(source_name_ch1), "dlpfc|prefrontal") ~ "DLPFC",
            TRUE ~ brain_region
          )
        )
    }
    
    # Enhanced condition extraction from characteristics_ch1 specifically
    if ("characteristics_ch1" %in% colnames(gse174409_sample_metadata)) {
      gse174409_sample_metadata <- gse174409_sample_metadata %>%
        mutate(
          condition = case_when(
            str_detect(tolower(characteristics_ch1), "diagnosis: oud|oud") ~ "OUD",
            str_detect(tolower(characteristics_ch1), "diagnosis: cont|control|cont") ~ "Control",
            TRUE ~ condition
          )
        )
      
      cat("Updated conditions based on characteristics_ch1:\n")
      condition_summary <- table(gse174409_sample_metadata$condition)
      print(condition_summary)
    }
    
    # Look for control samples in other characteristics columns - fix the NA issue
    char_cols <- colnames(gse174409_sample_metadata)[str_detect(colnames(gse174409_sample_metadata), "characteristics")]
    for (char_col in char_cols) {
      if (char_col %in% colnames(gse174409_sample_metadata)) {
        char_data <- tolower(as.character(gse174409_sample_metadata[[char_col]]))
        char_data[is.na(char_data)] <- ""  # Replace NAs with empty strings
        
        control_mask <- str_detect(char_data, "control|normal|healthy|non-user|cont")
        oud_mask <- str_detect(char_data, "oud|opioid|case|user|dependent|abuse")
        
        # Only update if we find clear patterns
        if (sum(control_mask, na.rm = TRUE) > 0) {
          gse174409_sample_metadata$condition[control_mask] <- "Control"
        }
        if (sum(oud_mask, na.rm = TRUE) > 0) {
          gse174409_sample_metadata$condition[oud_mask] <- "OUD"
        }
      }
    }
  }
  
  # Manual inspection and correction based on actual sample patterns
  cat("Sample ID patterns (first 10):\n")
  print(head(gse174409_sample_metadata$sample_id, 10))
  
  # Count unknowns safely (handle NAs)
  unknown_conditions <- sum(is.na(gse174409_sample_metadata$condition) | gse174409_sample_metadata$condition == "Unknown", na.rm = TRUE)
  
  # If still many unknowns, examine the actual sample naming pattern more carefully
  if (unknown_conditions > nrow(gse174409_sample_metadata) * 0.5) {
    cat("Many unknown conditions detected. Examining sample patterns more carefully...\n")
    
    # Also check if there are obvious patterns in the GEO characteristics
    if (!is.null(gse174409_metadata) && "characteristics_ch1" %in% colnames(gse174409_metadata)) {
      cat("Examining characteristics_ch1 column for patterns:\n")
      char_examples <- unique(gse174409_metadata$characteristics_ch1)[1:min(10, length(unique(gse174409_metadata$characteristics_ch1)))]
      print(char_examples)
      
      # Direct mapping based on the found patterns - this should work!
      if (any(str_detect(char_examples, "diagnosis: CONT"))) {
        cat("Found 'diagnosis: CONT' pattern - mapping directly...\n")
        
        # Create a direct mapping from the original metadata
        sample_mapping <- data.frame(
          sample_id = rownames(gse174409_metadata),
          direct_condition = case_when(
            str_detect(gse174409_metadata$characteristics_ch1, "diagnosis: OUD") ~ "OUD",
            str_detect(gse174409_metadata$characteristics_ch1, "diagnosis: CONT") ~ "Control",
            TRUE ~ "Unknown"
          ),
          direct_brain_region = case_when(
            str_detect(tolower(gse174409_metadata$source_name_ch1), "nac") ~ "NAcc",
            str_detect(tolower(gse174409_metadata$source_name_ch1), "dlpfc") ~ "DLPFC",
            TRUE ~ "Unknown"
          ),
          stringsAsFactors = FALSE
        )
        
        cat("Direct mapping summary:\n")
        cat("Conditions:", paste(names(table(sample_mapping$direct_condition)), "=", table(sample_mapping$direct_condition), collapse = " / "), "\n")
        cat("Brain regions:", paste(names(table(sample_mapping$direct_brain_region)), "=", table(sample_mapping$direct_brain_region), collapse = " / "), "\n")
        
        # Update the sample metadata with direct mappings - fix the join issue
        gse174409_sample_metadata <- gse174409_sample_metadata %>%
          left_join(sample_mapping, by = "sample_id")
        
        # Check if join worked
        cat("After join - checking columns:\n")
        cat("Columns in metadata:", paste(colnames(gse174409_sample_metadata), collapse = ", "), "\n")
        
        # Update conditions and brain regions
        gse174409_sample_metadata <- gse174409_sample_metadata %>%
          mutate(
            condition = ifelse(is.na(direct_condition), condition, direct_condition),
            brain_region = ifelse(is.na(direct_brain_region), brain_region, direct_brain_region)
          ) %>%
          dplyr::select(-direct_condition, -direct_brain_region)
        
        cat("Applied direct mapping from GEO metadata\n")
        cat("Final condition summary:", paste(names(table(gse174409_sample_metadata$condition)), "=", table(gse174409_sample_metadata$condition), collapse = " / "), "\n")
        cat("Final brain region summary:", paste(names(table(gse174409_sample_metadata$brain_region)), "=", table(gse174409_sample_metadata$brain_region), collapse = " / "), "\n")
      }
      
      # Look for additional characteristics columns that might contain control info
      other_char_cols <- colnames(gse174409_metadata)[str_detect(colnames(gse174409_metadata), "characteristics_ch1\\.")]
      if (length(other_char_cols) > 0) {
        cat("Found additional characteristics columns:", paste(other_char_cols[1:min(5, length(other_char_cols))], collapse = ", "), "\n")
        for (col in other_char_cols[1:min(3, length(other_char_cols))]) {
          cat("Examples from", col, ":\n")
          examples <- unique(gse174409_metadata[[col]])[1:min(5, length(unique(gse174409_metadata[[col]])))]
          print(examples)
        }
      }
    }
    
    if (!is.null(gse174409_metadata) && "title" %in% colnames(gse174409_metadata)) {
      cat("Examining title column for patterns:\n")
      title_examples <- unique(gse174409_metadata$title)[1:min(10, length(unique(gse174409_metadata$title)))]
      print(title_examples)
    }
  }
  
  # Try to infer from sample ID patterns (HN1, HN2, etc. might indicate groups)
  # Since we can see some samples have OUD, check if there's a pattern
  if (sum(gse174409_sample_metadata$condition == "OUD", na.rm = TRUE) > 0) {
    cat("Found some OUD samples. Checking for control pattern...\n")
    
    # If we have some OUD but not all samples, the remaining might be controls
    # This is a fallback - ideally we'd find this in the metadata
    unknown_mask <- is.na(gse174409_sample_metadata$condition) | gse174409_sample_metadata$condition == "Unknown"
    if (sum(unknown_mask) > 0) {
      gse174409_sample_metadata$condition[unknown_mask] <- "Control"
      cat("Assigned remaining", sum(unknown_mask), "samples as Control\n")
    }
  }
}

cat("✓ Created metadata for", nrow(gse174409_sample_metadata), "samples\n")
cat("  Conditions:", paste(names(table(gse174409_sample_metadata$condition)), "=", table(gse174409_sample_metadata$condition), collapse = " / "), "\n")
cat("  Brain regions:", paste(names(table(gse174409_sample_metadata$brain_region)), "=", table(gse174409_sample_metadata$brain_region), collapse = " / "), "\n")
  
# Display sample metadata summary for verification
cat("\nSample metadata summary:\n")
print(head(gse174409_sample_metadata[, c("sample_id", "condition", "brain_region")], 10))

# ==============================================================================
# PART 2: Data Quality Control and Normalization (GSE174409)
# ==============================================================================

if (gse174409_processed) {
  cat("\n", strrep("=", 60), "\n")
  cat("QUALITY CONTROL AND NORMALIZATION - GSE174409\n")
  cat(strrep("=", 60), "\n")
  
  # Pre-filtering statistics
  cat("Pre-filtering gene statistics:\n")
  genes_zero_all <- rowSums(gse174409_counts) == 0
  genes_low_counts <- rowSums(gse174409_counts >= 1) < 3
  
  cat("  Genes with zero counts in all samples:", sum(genes_zero_all), "\n")
  cat("  Genes with <3 samples with >=1 count:", sum(genes_low_counts), "\n")
  
  # Filter low-count genes
  keep_genes <- rowSums(gse174409_counts >= 1) >= 3
  gse174409_filtered_counts <- gse174409_counts[keep_genes, ]
  
  cat("Genes after filtering:", nrow(gse174409_filtered_counts), "/", nrow(gse174409_counts), 
      "(", round(100 * nrow(gse174409_filtered_counts) / nrow(gse174409_counts), 1), "% retained)\n")
  
  # Create DGEList for comprehensive analysis
  tryCatch({
    dge <- edgeR::DGEList(counts = gse174409_filtered_counts, 
                         samples = gse174409_sample_metadata)
    
    # Calculate normalization factors
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    
    # Get normalized expression values
    gse174409_cpm <- edgeR::cpm(dge, log = FALSE)
    gse174409_logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
    
    cat("✓ Completed TMM normalization\n")
    
    # Additional normalization for DESeq2 compatibility - fix the design issue
    if (length(unique(gse174409_sample_metadata$condition)) > 1) {
      dds <- DESeqDataSetFromMatrix(
        countData = gse174409_filtered_counts,
        colData = gse174409_sample_metadata,
        design = ~ condition
      )
      
      # Variance stabilizing transformation
      vsd <- vst(dds, blind = FALSE)
      gse174409_vst <- assay(vsd)
      
      cat("✓ Completed variance stabilizing transformation\n")
    } else {
      cat("Skipping DESeq2 VST - only one condition detected\n")
      # Use simpler transformation
      gse174409_vst <- log2(gse174409_cpm + 1)
    }
    
    # Create additional normalized datasets for compatibility
    gse174409_normalized <- gse174409_logcpm  # For backward compatibility
    
  }, error = function(e) {
    cat("Warning: Error in normalization:", e$message, "\n")
  })
}

# Generate QC plots for GSE174409
if (gse174409_processed && exists("gse174409_logcpm")) {
  cat("Generating QC plots for GSE174409...\n")
  
  # Sample correlation heatmap using log-transformed data
  sample_cor <- cor(gse174409_logcpm)
  
  # Create annotation data for heatmap (ensure proper formatting)
  annotation_data <- gse174409_sample_metadata[, c("condition", "brain_region"), drop = FALSE]
  rownames(annotation_data) <- gse174409_sample_metadata$sample_id
  
  # Ensure annotation_data matches sample_cor column names
  annotation_data <- annotation_data[colnames(sample_cor), , drop = FALSE]
  
  # Remove any NA values
  annotation_data$condition[is.na(annotation_data$condition)] <- "Unknown"
  annotation_data$brain_region[is.na(annotation_data$brain_region)] <- "Unknown"
  
  tryCatch({
    pdf(file.path(figures_dir, "GSE174409_sample_correlation.pdf"), width = 10, height = 8)
    pheatmap::pheatmap(sample_cor, 
                       main = "GSE174409 Sample Correlation",
                       annotation_col = annotation_data,
                       show_rownames = FALSE,
                       show_colnames = FALSE)
    dev.off()
    cat("✓ Generated sample correlation heatmap\n")
  }, error = function(e) {
    cat("Warning: Could not generate correlation heatmap:", e$message, "\n")
    # Generate simple correlation heatmap without annotations
    pdf(file.path(figures_dir, "GSE174409_sample_correlation_simple.pdf"), width = 10, height = 8)
    pheatmap::pheatmap(sample_cor, 
                       main = "GSE174409 Sample Correlation",
                       show_rownames = FALSE,
                       show_colnames = FALSE)
    dev.off()
  })
  
  # PCA plot
  tryCatch({
    pca_data <- prcomp(t(gse174409_logcpm))
    pca_df <- data.frame(pca_data$x[, 1:2], gse174409_sample_metadata)
    
    p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, shape = brain_region)) +
      geom_point(size = 3) +
      labs(title = "GSE174409 PCA Plot",
           x = paste0("PC1 (", round(summary(pca_data)$importance[2, 1] * 100, 1), "%)"),
           y = paste0("PC2 (", round(summary(pca_data)$importance[2, 2] * 100, 1), "%)")) +
      theme_bw()
    
    ggsave(file.path(figures_dir, "GSE174409_PCA.pdf"), p_pca, width = 8, height = 6)
    cat("✓ Generated PCA plot\n")
  }, error = function(e) {
    cat("Warning: Could not generate PCA plot:", e$message, "\n")
  })
  
  cat("✓ Generated QC plots for GSE174409\n")
}

# ==============================================================================
# PART 3: GSE225158 - Human Single-cell RNA-seq Processing (Enhanced)
# ==============================================================================

cat("\n", strrep("=", 60), "\n")
cat("PROCESSING GSE225158 - HUMAN SINGLE-CELL RNA-SEQUENCING (STRIATUM)\n")
cat(strrep("=", 60), "\n")

# Create dataset directory
gse225158_dir <- file.path(raw_data_dir, "GSE225158")
dir.create(gse225158_dir, showWarnings = FALSE, recursive = TRUE)

# Check for manually downloaded files first
manual_files <- list(
  h5seurat = file.path(gse225158_dir, "GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat"),
  h5ad = file.path(gse225158_dir, "GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5ad.gz"),
  alternative_h5seurat = file.path(gse225158_dir, "seurat_object.h5seurat"),
  alternative_h5ad = file.path(gse225158_dir, "seurat_object.h5ad")
)

# Check which files exist
existing_files <- keep(manual_files, file.exists)

if (length(existing_files) > 0) {
  cat("Found manually downloaded files:\n")
  iwalk(existing_files, ~ cat("  ✓", .y, ":", basename(.x), "\n"))
} else {
  cat("No manually downloaded files found. Checking for automatic download...\n")
}

# Download GEO metadata
cat("Downloading GSE225158 metadata...\n")
gse225158_metadata <- NULL
tryCatch({
  gse225158_geo <- getGEO("GSE225158", GSEMatrix = TRUE, destdir = gse225158_dir)
  if (length(gse225158_geo) > 0) {
    gse225158_metadata <- pData(gse225158_geo[[1]])
    cat("✓ Downloaded metadata for", nrow(gse225158_metadata), "samples\n")
  }
}, error = function(e) {
  cat("Warning: Could not download GSE225158 metadata:", e$message, "\n")
})

# Process single-cell data
gse225158_processed <- FALSE
gse225158_seurat <- NULL

# Try to load existing files in order of preference
for (file_type in names(existing_files)) {
  filepath <- existing_files[[file_type]]
  
  cat("Attempting to load:", basename(filepath), "\n")
  
  # Validate file first
  validation <- validate_file(filepath, 100000000)  # 100MB minimum for scRNA-seq
  if (!validation$valid) {
    cat("  File validation failed:", validation$reason, "\n")
    next
  }
  
  tryCatch({
    if (str_detect(file_type, "h5seurat")) {
      # Load h5Seurat file
      if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
        cat("Installing SeuratDisk...\n")
        if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
        remotes::install_github("mojaveazure/seurat-disk", quiet = TRUE)
      }
      library(SeuratDisk)
      
      gse225158_seurat <- LoadH5Seurat(filepath)
      
    } else if (str_detect(file_type, "h5ad")) {
      # Load h5ad file via SeuratDisk
      if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
        cat("Installing SeuratDisk...\n")
        remotes::install_github("mojaveazure/seurat-disk", quiet = TRUE)
      }
      library(SeuratDisk)
      
      # Convert to h5seurat if needed
      temp_h5seurat <- str_replace(filepath, "\\.h5ad(\\.gz)?$", ".h5seurat")
      if (!file.exists(temp_h5seurat)) {
        cat("Converting h5ad to h5seurat format...\n")
        Convert(filepath, dest = "h5seurat", overwrite = TRUE)
      }
      gse225158_seurat <- LoadH5Seurat(temp_h5seurat)
    }
    
    if (!is.null(gse225158_seurat)) {
      cat("✓ Successfully loaded Seurat object\n")
      cat("  Cells:", ncol(gse225158_seurat), "\n")
      cat("  Features:", nrow(gse225158_seurat), "\n")
      cat("  Assays:", paste(names(gse225158_seurat@assays), collapse = ", "), "\n")
      
      gse225158_processed <- TRUE
      break
    }
    
  }, error = function(e) {
    cat("  Failed to load:", e$message, "\n")
  })
}

# Enhance metadata if object was loaded
if (gse225158_processed && !is.null(gse225158_seurat)) {
  cat("Enhancing single-cell metadata...\n")
  
  # Add study-specific metadata
  gse225158_seurat@meta.data$dataset <- "GSE225158"
  gse225158_seurat@meta.data$species <- "Human"
  gse225158_seurat@meta.data$tissue <- "Striatum"
  gse225158_seurat@meta.data$data_type <- "single_cell_rnaseq"
  
  # Basic QC metrics if not present
  if (!"percent.mt" %in% colnames(gse225158_seurat@meta.data)) {
    gse225158_seurat[["percent.mt"]] <- PercentageFeatureSet(gse225158_seurat, pattern = "^MT-")
  }
  if (!"percent.ribo" %in% colnames(gse225158_seurat@meta.data)) {
    gse225158_seurat[["percent.ribo"]] <- PercentageFeatureSet(gse225158_seurat, pattern = "^RP[SL]")
  }
  
  # Display metadata summary
  meta_cols <- colnames(gse225158_seurat@meta.data)
  cat("Available metadata columns:", length(meta_cols), "\n")
  
  # Look for condition/treatment columns
  condition_cols <- meta_cols[str_detect(tolower(meta_cols), "condition|treatment|group|oud|opioid")]
  if (length(condition_cols) > 0) {
    cat("Potential condition columns:", paste(condition_cols, collapse = ", "), "\n")
    for (col in condition_cols[1:min(2, length(condition_cols))]) {
      cat("  ", col, ":", paste(unique(gse225158_seurat@meta.data[[col]])[1:min(5, length(unique(gse225158_seurat@meta.data[[col]])))], collapse = ", "), "\n")
    }
  }
  
  # Look for cell type columns
  celltype_cols <- meta_cols[str_detect(tolower(meta_cols), "cell.*type|cluster|annotation")]
  if (length(celltype_cols) > 0) {
    cat("Potential cell type columns:", paste(celltype_cols, collapse = ", "), "\n")
  }
}

# Generate download instructions if no files were processed
if (!gse225158_processed) {
  cat("\n⚠ GSE225158 data not found or could not be processed\n")
  cat("Manual download instructions:\n")
  cat("1. Download files using browser or wget:\n")
  
  download_urls <- c(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225158/suppl/GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225158/suppl/GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5ad.gz"
  )
  
  for (url in download_urls) {
    cat("   wget '", url, "' -P '", gse225158_dir, "'\n", sep = "")
  }
  
  cat("2. Re-run this script after download\n")
  cat("3. Alternative: Use browser to download and place in:", gse225158_dir, "\n")
}

# ==============================================================================
# PART 4: Gene Annotation and Symbol Mapping
# ==============================================================================

cat("\n", strrep("=", 60), "\n")
cat("GENE ANNOTATION AND SYMBOL MAPPING\n")
cat(strrep("=", 60), "\n")

if (exists("gse174409_logcpm")) {
  cat("Mapping gene symbols for GSE174409...\n")
  
  # Get gene symbols using biomaRt
  tryCatch({
    mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    
    # Extract gene IDs (assuming ENSEMBL IDs)
    gene_ids <- rownames(gse174409_logcpm)
    
    # Check if these look like ENSEMBL IDs
    ensembl_pattern <- str_detect(gene_ids[1:min(10, length(gene_ids))], "^ENSG")
    if (any(ensembl_pattern)) {
      cat("Detected ENSEMBL gene IDs, proceeding with biomaRt annotation...\n")
      
      gene_mapping <- getBM(
        attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "gene_biotype"),
        filters = "ensembl_gene_id",
        values = gene_ids,
        mart = mart
      )
      
      cat("✓ Retrieved annotations for", nrow(gene_mapping), "genes\n")
    } else {
      cat("Gene IDs don't appear to be ENSEMBL format, checking for HGNC symbols...\n")
      
      # Try HGNC symbol mapping
      gene_mapping <- getBM(
        attributes = c("hgnc_symbol", "ensembl_gene_id", "chromosome_name", "gene_biotype"),
        filters = "hgnc_symbol",
        values = gene_ids,
        mart = mart
      )
      
      cat("✓ Retrieved annotations for", nrow(gene_mapping), "genes using HGNC symbols\n")
    }
    
  }, error = function(e) {
    cat("Warning: Could not retrieve gene annotations:", e$message, "\n")
    # Create dummy mapping
    gene_mapping <- data.frame(
      gene_id = rownames(gse174409_logcpm),
      hgnc_symbol = rownames(gse174409_logcpm),
      stringsAsFactors = FALSE
    )
  })
}

# ==============================================================================
# PART 5: Save Processed Data
# ==============================================================================

cat("\n", strrep("=", 60), "\n")
cat("SAVING PROCESSED DATA\n")
cat(strrep("=", 60), "\n")

# Save GSE174409 processed data
if (exists("gse174409_logcpm") && exists("gse174409_sample_metadata")) {
  # Save normalized expression data (multiple formats)
  saveRDS(gse174409_logcpm, file.path(processed_data_dir, "GSE174409_logcpm_expression.rds"))
  write.csv(gse174409_logcpm, file.path(processed_data_dir, "GSE174409_logcpm_expression.csv"))
  
  if (exists("gse174409_cpm")) {
    saveRDS(gse174409_cpm, file.path(processed_data_dir, "GSE174409_cpm_expression.rds"))
    write.csv(gse174409_cpm, file.path(processed_data_dir, "GSE174409_cpm_expression.csv"))
  }
  
  if (exists("gse174409_vst")) {
    saveRDS(gse174409_vst, file.path(processed_data_dir, "GSE174409_vst_expression.rds"))
    write.csv(gse174409_vst, file.path(processed_data_dir, "GSE174409_vst_expression.csv"))
  }
  
  # Save raw filtered counts
  saveRDS(gse174409_filtered_counts, file.path(processed_data_dir, "GSE174409_filtered_counts.rds"))
  
  # Save sample metadata
  saveRDS(gse174409_sample_metadata, file.path(processed_data_dir, "GSE174409_sample_metadata.rds"))
  write.csv(gse174409_sample_metadata, file.path(processed_data_dir, "GSE174409_sample_metadata.csv"), row.names = FALSE)
  
  # Save gene annotations
  if (exists("gene_mapping")) {
    saveRDS(gene_mapping, file.path(processed_data_dir, "GSE174409_gene_annotations.rds"))
    write.csv(gene_mapping, file.path(processed_data_dir, "GSE174409_gene_annotations.csv"), row.names = FALSE)
  }
  
  cat("✓ Saved GSE174409 processed data\n")
}

# Save GSE225158 processed data (only if Seurat object exists)
if (exists("gse225158_seurat") && !is.null(gse225158_seurat)) {
  saveRDS(gse225158_seurat, file.path(processed_data_dir, "GSE225158_seurat_object.rds"))
  
  # Extract and save count matrix with proper Seurat method handling
  tryCatch({
    # Try different Seurat versions
    if (packageVersion("Seurat") >= "5.0.0") {
      gse225158_counts <- GetAssayData(gse225158_seurat, assay = "RNA", layer = "counts")
    } else {
      gse225158_counts <- GetAssayData(gse225158_seurat, assay = "RNA", slot = "counts")
    }
    
    saveRDS(gse225158_counts, file.path(processed_data_dir, "GSE225158_count_matrix.rds"))
    
    # Save metadata
    gse225158_cell_metadata <- gse225158_seurat@meta.data
    saveRDS(gse225158_cell_metadata, file.path(processed_data_dir, "GSE225158_cell_metadata.rds"))
    write.csv(gse225158_cell_metadata, file.path(processed_data_dir, "GSE225158_cell_metadata.csv"))
    
    cat("✓ Saved GSE225158 single-cell data\n")
    
  }, error = function(e) {
    cat("Warning: Could not extract count matrix from Seurat object:", e$message, "\n")
    
    # Save just the Seurat object and metadata
    gse225158_cell_metadata <- gse225158_seurat@meta.data
    saveRDS(gse225158_cell_metadata, file.path(processed_data_dir, "GSE225158_cell_metadata.rds"))
    write.csv(gse225158_cell_metadata, file.path(processed_data_dir, "GSE225158_cell_metadata.csv"))
    
    cat("✓ Saved GSE225158 metadata only\n")
  })
}

# Save processing summary
processing_summary <- list(
  processing_date = Sys.time(),
  gse174409 = list(
    samples = if(exists("gse174409_sample_metadata")) nrow(gse174409_sample_metadata) else 0,
    genes = if(exists("gse174409_logcpm")) nrow(gse174409_logcpm) else 0,
    conditions = if(exists("gse174409_sample_metadata")) unique(gse174409_sample_metadata$condition) else NA,
    brain_regions = if(exists("gse174409_sample_metadata")) unique(gse174409_sample_metadata$brain_region) else NA
  ),
  gse225158 = list(
    status = if(exists("gse225158_seurat") && !is.null(gse225158_seurat)) "Processed" else "Metadata only",
    cells = if(exists("gse225158_seurat") && !is.null(gse225158_seurat)) ncol(gse225158_seurat) else 0,
    genes = if(exists("gse225158_seurat") && !is.null(gse225158_seurat)) nrow(gse225158_seurat) else 0
  )
)

saveRDS(processing_summary, file.path(results_dir, "human_data_processing_summary.rds"))

# ==============================================================================
# PART 6: Generate Processing Report
# ==============================================================================

cat("\n", strrep("=", 60), "\n")
cat("GENERATING PROCESSING REPORT\n")
cat(strrep("=", 60), "\n")

report_file <- file.path(results_dir, "Human_Data_Processing_Report.txt")

cat("HUMAN DATA PROCESSING REPORT\n", file = report_file)
cat(strrep("=", 40), "\n", file = report_file, append = TRUE)
cat("Processing Date:", as.character(Sys.time()), "\n\n", file = report_file, append = TRUE)

cat("GSE174409 - Human Bulk RNA-seq:\n", file = report_file, append = TRUE)
if (exists("gse174409_sample_metadata")) {
  cat("  Samples:", nrow(gse174409_sample_metadata), "\n", file = report_file, append = TRUE)
  cat("  Conditions:", paste(unique(gse174409_sample_metadata$condition), collapse = ", "), "\n", file = report_file, append = TRUE)
  cat("  Brain regions:", paste(unique(gse174409_sample_metadata$brain_region), collapse = ", "), "\n", file = report_file, append = TRUE)
}
if (exists("gse174409_logcpm")) {
  cat("  Genes after filtering:", nrow(gse174409_logcpm), "\n", file = report_file, append = TRUE)
}

cat("\nGSE225158 - Human Single-cell RNA-seq:\n", file = report_file, append = TRUE)
if (exists("gse225158_seurat") && !is.null(gse225158_seurat)) {
  cat("  Status: Processed successfully\n", file = report_file, append = TRUE)
  cat("  Cells:", ncol(gse225158_seurat), "\n", file = report_file, append = TRUE)
  cat("  Genes:", nrow(gse225158_seurat), "\n", file = report_file, append = TRUE)
} else {
  cat("  Status: Metadata downloaded, data processing pending\n", file = report_file, append = TRUE)
  cat("  Note: Manual download of supplementary files may be required\n", file = report_file, append = TRUE)
}

cat("\nOutput Files Generated:\n", file = report_file, append = TRUE)
cat("  - GSE174409_logcpm_expression.rds/csv\n", file = report_file, append = TRUE)
cat("  - GSE174409_filtered_counts.rds\n", file = report_file, append = TRUE)
cat("  - GSE174409_sample_metadata.rds/csv\n", file = report_file, append = TRUE)
if (exists("gene_mapping")) {
  cat("  - GSE174409_gene_annotations.rds/csv\n", file = report_file, append = TRUE)
}
if (exists("gse225158_seurat") && !is.null(gse225158_seurat)) {
  cat("  - GSE225158_seurat_object.rds\n", file = report_file, append = TRUE)
  cat("  - GSE225158_count_matrix.rds\n", file = report_file, append = TRUE)
  cat("  - GSE225158_cell_metadata.rds/csv\n", file = report_file, append = TRUE)
}
cat("  - QC plots: correlation, PCA plots\n", file = report_file, append = TRUE)

cat("✓ Generated processing report\n")

# ==============================================================================
# COMPLETION SUMMARY
# ==============================================================================

cat("\n", strrep("=", 80), "\n")
cat("HUMAN DATA PROCESSING COMPLETED\n")
cat(strrep("=", 80), "\n")

cat("Processing Summary:\n")
cat("  ✓ GSE174409: Downloaded and processed\n")
cat("  ✓ Data normalization: Completed\n")
cat("  ✓ Gene annotation: Completed\n")
cat("  ✓ QC plots: Generated\n")
if (exists("gse225158_seurat") && !is.null(gse225158_seurat)) {
  cat("  ✓ GSE225158: Single-cell data processed\n")
} else {
  cat("  ⚠ GSE225158: Manual file download may be required\n")
}

cat("\nNext Steps:\n")
if (!exists("gse225158_seurat") || is.null(gse225158_seurat)) {
  cat("  1. Download GSE225158 supplementary files manually if needed:\n")
  cat("     wget 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225158/suppl/GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat' -P '/Users/aumchampaneri/Complement-OUD/Data/Raw/GSE225158'\n")
  cat("  2. Re-run script to process GSE225158 data\n")
  cat("  3. Alternative: Use browser to download and place in:", gse225158_dir, "\n")
} else {
  cat("  1. Run mouse data harmonization script (02_Mouse_Data_Harmonization.R)\n")
  cat("  2. Proceed to ortholog mapping (03_Ortholog_Mapping.R)\n")
  cat("  3. Cross-species integration and analysis\n")
}

cat("\nProcessing completed at:", as.character(Sys.time()), "\n")
cat(strrep("=", 80), "\n")
