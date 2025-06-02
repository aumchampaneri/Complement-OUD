#!/usr/bin/env Rscript
# ==============================================================================
# Mouse Data Harmonization
# ==============================================================================
# Project: Cross-species meta-analysis of mouse and human OUD datasets
# Author: Bioinformatics Analysis Pipeline
# Date: June 2025
# 
# Purpose: Harmonize existing mouse datasets for cross-species comparison
# - GSE118918: Bulk RNA-seq, NAcc, males only, acute morphine
# - GSE207128: Single-cell RNA-seq, Amygdala, males only, chronic dependence/withdrawal
# - GSE289002: Bulk RNA-seq, NAc+PFC, both sexes, temporal progression
# ==============================================================================

# Load required packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(here)
  library(DESeq2)
  library(Seurat)
  library(edgeR)
  library(biomaRt)
  library(openxlsx)
  library(ggplot2)
  library(ComplexHeatmap)
  library(pheatmap)
  library(RColorBrewer)
  library(viridis)
  library(janitor)
})

# Set up directories
project_root <- getwd()  # Should be "Translational Study" folder
mouse_data_dir <- file.path(dirname(project_root), "Super Folder - Mus musculus")
processed_dir <- file.path(project_root, "Data", "Processed")
results_dir <- file.path(project_root, "Results", "Mouse_Analysis")
figures_dir <- file.path(project_root, "Figures", "QC")

# Create directories if they don't exist
dirs_to_create <- c(processed_dir, results_dir, figures_dir)
walk(dirs_to_create, ~ if (!dir.exists(.x)) dir.create(.x, recursive = TRUE))

cat(strrep("=", 80), "\n")
cat("MOUSE DATA HARMONIZATION - PHASE 1\n")
cat(strrep("=", 80), "\n")
cat("Processing Time:", Sys.time(), "\n\n")

cat("Project root:", project_root, "\n")
cat("Mouse data directory:", mouse_data_dir, "\n")
cat("Output directories set for Translational Study folder\n\n")

# Check if mouse data directory exists
if (!dir.exists(mouse_data_dir)) {
  cat("Warning: Mouse data directory not found at:", mouse_data_dir, "\n")
  cat("Please check the directory structure.\n")
  cat("Expected structure:\n")
  cat("  - /Users/aumchampaneri/Complement-OUD/\n")
  cat("    ├── Translational Study/ (current working directory)\n")
  cat("    └── Super Folder - Mus musculus/\n")
  cat("        ├── GSE118918/\n")
  cat("        ├── GSE207128/\n")
  cat("        └── GSE289002/\n\n")
  
  # Try alternative locations
  alternative_paths <- c(
    file.path(project_root, "Super Folder - Mus musculus"),
    file.path(dirname(dirname(project_root)), "Super Folder - Mus musculus"),
    file.path("/Users/aumchampaneri/Complement-OUD", "Super Folder - Mus musculus")
  )
  
  for (alt_path in alternative_paths) {
    if (dir.exists(alt_path)) {
      cat("Found mouse data directory at:", alt_path, "\n")
      mouse_data_dir <- alt_path
      break
    }
  }
}

# List available datasets
if (dir.exists(mouse_data_dir)) {
  available_datasets <- list.dirs(mouse_data_dir, recursive = FALSE, full.names = FALSE)
  cat("Available mouse datasets:", paste(available_datasets, collapse = ", "), "\n\n")
} else {
  cat("Error: Could not locate mouse data directory. Exiting.\n")
  quit(save = "no", status = 1)
}

# ==============================================================================
# FUNCTION DEFINITIONS
# ==============================================================================

# Function to standardize gene symbols
standardize_gene_symbols <- function(gene_names, species = "mouse") {
  cat("Standardizing gene symbols for", species, "...\n")
  
  # Remove any trailing/leading whitespace
  gene_names <- str_trim(gene_names)
  
  # Convert to standard mouse gene naming (first letter uppercase, rest lowercase)
  if (species == "mouse") {
    # Handle different input formats
    gene_names_std <- gene_names
    
    # Handle ENSEMBL IDs - keep as is for now, could convert later
    ensembl_ids <- str_detect(gene_names, "^ENSMUSG")
    if (sum(ensembl_ids) > 0) {
      cat("Found", sum(ensembl_ids), "ENSEMBL gene IDs - keeping original format\n")
      # Keep ENSEMBL IDs as-is for now
      gene_names_std[ensembl_ids] <- gene_names[ensembl_ids]
    }
    
    # Process non-ENSEMBL genes
    non_ensembl <- !ensembl_ids
    if (sum(non_ensembl) > 0) {
      # Convert from all caps to title case
      all_caps <- str_detect(gene_names[non_ensembl], "^[A-Z0-9-]+$")
      gene_names_std[non_ensembl][all_caps] <- str_to_title(tolower(gene_names[non_ensembl][all_caps]))
      
      # Convert from all lowercase to title case
      all_lower <- str_detect(gene_names[non_ensembl], "^[a-z0-9-]+$")
      gene_names_std[non_ensembl][all_lower] <- str_to_title(gene_names[non_ensembl][all_lower])
      
      # Handle genes that are already in mixed case
      mixed_case <- !all_caps & !all_lower
      gene_names_std[non_ensembl][mixed_case] <- str_to_title(tolower(gene_names[non_ensembl][mixed_case]))
    }
    
    # Handle special cases for mouse genes (only for non-ENSEMBL)
    gene_names_std <- str_replace_all(gene_names_std, "^Mt-", "mt-")  # Mitochondrial genes
    gene_names_std <- str_replace_all(gene_names_std, "^mt-", "mt-")  # Already lowercase mt
    gene_names_std <- str_replace_all(gene_names_std, "^MT-", "mt-")  # Uppercase MT
    
    # Handle ribosomal proteins (only for non-ENSEMBL)
    gene_names_std <- str_replace_all(gene_names_std, "^Rps([0-9])", "Rps\\1")  # Ribosomal protein S
    gene_names_std <- str_replace_all(gene_names_std, "^Rpl([0-9])", "Rpl\\1")  # Ribosomal protein L
    gene_names_std <- str_replace_all(gene_names_std, "^RPS([0-9])", "Rps\\1")  # Uppercase RPS
    gene_names_std <- str_replace_all(gene_names_std, "^RPL([0-9])", "Rpl\\1")  # Uppercase RPL
    
    # Show some examples of standardization
    if (length(gene_names) > 0) {
      n_show <- min(5, length(gene_names))
      cat("Gene standardization examples:\n")
      changes_shown <- 0
      for (i in 1:length(gene_names)) {
        if (gene_names[i] != gene_names_std[i] && changes_shown < n_show) {
          cat("  ", gene_names[i], "->", gene_names_std[i], "\n")
          changes_shown <- changes_shown + 1
        }
      }
      if (changes_shown == 0) {
        cat("  No changes needed for first", n_show, "genes\n")
      }
    }
  }
  
  return(gene_names_std)
}

# Function to load and process mouse dataset
load_mouse_dataset <- function(dataset_path, dataset_name) {
  cat("Loading", dataset_name, "from", dataset_path, "...\n")
  
  if (!dir.exists(dataset_path)) {
    cat("Dataset directory does not exist:", dataset_path, "\n")
    return(NULL)
  }
  
  # Look for processed data files in the dataset directory
  all_files <- list.files(dataset_path, recursive = TRUE, full.names = TRUE)
  cat("Found", length(all_files), "total files in directory\n")
  
  # Prioritize different file types
  rds_files <- all_files[str_detect(all_files, "\\.rds$")]
  h5_files <- all_files[str_detect(all_files, "\\.h5$|\\.h5seurat$")]
  excel_files <- all_files[str_detect(all_files, "\\.xlsx$|\\.xls$")]
  csv_files <- all_files[str_detect(all_files, "\\.csv$")]
  txt_files <- all_files[str_detect(all_files, "\\.txt$|\\.tsv$")]
  
  cat("File types found:\n")
  cat("  RDS files:", length(rds_files), "\n")
  cat("  H5 files:", length(h5_files), "\n")
  cat("  Excel files:", length(excel_files), "\n")
  cat("  CSV files:", length(csv_files), "\n")
  cat("  Text files:", length(txt_files), "\n")
  
  # Look for specific patterns that indicate processed data
  processed_patterns <- c("processed", "counts", "expression", "normalized", "dge", "seurat", "final")
  
  # Try to load data in order of preference
  
  # 1. Try RDS files first (most likely to contain R objects)
  if (length(rds_files) > 0) {
    # Prioritize files with processing keywords
    priority_rds <- rds_files[str_detect(tolower(basename(rds_files)), 
                                        paste(processed_patterns, collapse = "|"))]
    
    files_to_try <- c(priority_rds, rds_files[!rds_files %in% priority_rds])
    
    for (rds_file in files_to_try) {
      cat("Attempting to load RDS file:", basename(rds_file), "\n")
      tryCatch({
        data <- readRDS(rds_file)
        cat("✓ Successfully loaded RDS file:", basename(rds_file), "\n")
        cat("  Data type:", class(data), "\n")
        
        if (is.list(data)) {
          cat("  List elements:", paste(names(data), collapse = ", "), "\n")
        }
        
        return(list(data = data, file_type = "rds", file_path = rds_file))
      }, error = function(e) {
        cat("  Could not load", basename(rds_file), ":", e$message, "\n")
      })
    }
  }
  
  # 2. Try H5 files (Seurat objects)
  if (length(h5_files) > 0) {
    for (h5_file in h5_files) {
      cat("Attempting to load H5 file:", basename(h5_file), "\n")
      tryCatch({
        # Try different H5 loading methods
        if (str_detect(h5_file, "\\.h5seurat$")) {
          # Load H5Seurat file
          if (requireNamespace("SeuratDisk", quietly = TRUE)) {
            library(SeuratDisk)
            data <- LoadH5Seurat(h5_file)
          } else {
            cat("  SeuratDisk package required for H5Seurat files\n")
            next
          }
        } else {
          # Try standard H5 file
          data <- Read10X_h5(h5_file)
        }
        
        cat("✓ Successfully loaded H5 file:", basename(h5_file), "\n")
        return(list(data = data, file_type = "h5", file_path = h5_file))
      }, error = function(e) {
        cat("  Could not load", basename(h5_file), ":", e$message, "\n")
      })
    }
  }
  
  # 3. Try Excel files
  if (length(excel_files) > 0) {
    for (excel_file in excel_files) {
      cat("Attempting to load Excel file:", basename(excel_file), "\n")
      tryCatch({
        # Get sheet names
        sheet_names <- getSheetNames(excel_file)
        cat("  Available sheets:", paste(sheet_names, collapse = ", "), "\n")
        
        # Try to find the main data sheet
        data_sheet <- sheet_names[str_detect(tolower(sheet_names), 
                                           "count|expression|data|matrix")][1]
        if (is.na(data_sheet)) data_sheet <- sheet_names[1]
        
        data <- read.xlsx(excel_file, sheet = data_sheet)
        cat("✓ Successfully loaded Excel file:", basename(excel_file), "\n")
        cat("  Sheet used:", data_sheet, "\n")
        cat("  Dimensions:", nrow(data), "x", ncol(data), "\n")
        
        return(list(data = data, file_type = "excel", file_path = excel_file))
      }, error = function(e) {
        cat("  Could not load", basename(excel_file), ":", e$message, "\n")
      })
    }
  }
  
  # 4. Try CSV files
  if (length(csv_files) > 0) {
    # Prioritize larger files (likely to be count matrices)
    csv_sizes <- file.info(csv_files)$size
    csv_files_sorted <- csv_files[order(csv_sizes, decreasing = TRUE)]
    
    for (csv_file in csv_files_sorted[1:min(3, length(csv_files_sorted))]) {
      cat("Attempting to load CSV file:", basename(csv_file), "\n")
      tryCatch({
        data <- fread(csv_file)
        cat("✓ Successfully loaded CSV file:", basename(csv_file), "\n")
        cat("  Dimensions:", nrow(data), "x", ncol(data), "\n")
        
        return(list(data = data, file_type = "csv", file_path = csv_file))
      }, error = function(e) {
        cat("  Could not load", basename(csv_file), ":", e$message, "\n")
      })
    }
  }
  
  # 5. Try text files
  if (length(txt_files) > 0) {
    txt_sizes <- file.info(txt_files)$size
    txt_files_sorted <- txt_files[order(txt_sizes, decreasing = TRUE)]
    
    for (txt_file in txt_files_sorted[1:min(3, length(txt_files_sorted))]) {
      cat("Attempting to load text file:", basename(txt_file), "\n")
      tryCatch({
        # Try different separators
        separators <- c("\t", ",", " ")
        for (sep in separators) {
          data <- fread(txt_file, sep = sep)
          if (ncol(data) > 1 && nrow(data) > 1) {
            cat("✓ Successfully loaded text file:", basename(txt_file), "\n")
            cat("  Separator:", ifelse(sep == "\t", "tab", 
                                    ifelse(sep == ",", "comma", "space")), "\n")
            cat("  Dimensions:", nrow(data), "x", ncol(data), "\n")
            
            return(list(data = data, file_type = "txt", file_path = txt_file))
          }
        }
      }, error = function(e) {
        cat("  Could not load", basename(txt_file), ":", e$message, "\n")
      })
    }
  }
  
  cat("Could not load any data files for", dataset_name, "\n")
  cat("Available files in directory:\n")
  print(basename(all_files)[1:min(10, length(all_files))])
  if (length(all_files) > 10) cat("... and", length(all_files) - 10, "more files\n")
  
  return(NULL)
}

# Function to extract metadata from dataset
extract_metadata <- function(data, dataset_name, data_type = "bulk") {
  cat("Extracting metadata for", dataset_name, "...\n")
  
  metadata <- data.frame(
    sample_id = character(),
    dataset = character(),
    condition = character(),
    treatment = character(),
    brain_region = character(),
    sex = character(),
    time_point = character(),
    data_type = character(),
    stringsAsFactors = FALSE
  )
  
  # Dataset-specific metadata extraction
  if (dataset_name == "GSE118918") {
    # Check if we have DGEList with sample info
    if (is.list(data) && "samples" %in% names(data)) {
      meta <- data$samples
      meta$sample_id <- rownames(meta)
    } else if ("DGEList" %in% class(data) && "samples" %in% names(data)) {
      meta <- data$samples
      meta$sample_id <- rownames(meta)
    } else {
      # Create basic metadata if not available
      n_samples <- ifelse(is.matrix(data), ncol(data), length(data))
      meta <- data.frame(
        sample_id = paste0("GSE118918_", 1:n_samples),
        stringsAsFactors = FALSE
      )
    }
    
    metadata <- meta %>%
      mutate(
        dataset = "GSE118918",
        brain_region = "NAcc",
        sex = "Male",
        time_point = "4h_post_injection",
        data_type = "bulk_rnaseq",
        # Use actual treatment info if available
        condition = ifelse("treatment" %in% colnames(meta), 
                          ifelse(meta$treatment == "Mock", "Control", "Morphine"),
                          ifelse(str_detect(tolower(sample_id), "mock|control|saline"), "Control", "Morphine")),
        treatment = ifelse("treatment" %in% colnames(meta),
                          ifelse(meta$treatment == "Mock", "Saline", "Morphine_20mg_kg"),
                          ifelse(str_detect(tolower(sample_id), "mock|control|saline"), "Saline", "Morphine_20mg_kg"))
      )
    
  } else if (dataset_name == "GSE207128") {
    # For Seurat objects, extract from meta.data
    if ("Seurat" %in% class(data)) {
      meta <- data@meta.data
      meta$sample_id <- rownames(meta)
    } else if (is.list(data) && "metadata" %in% names(data)) {
      meta <- data$metadata
    } else {
      n_samples <- ifelse(is.matrix(data), ncol(data), length(data))
      meta <- data.frame(
        sample_id = paste0("GSE207128_", 1:n_samples),
        stringsAsFactors = FALSE
      )
    }
    
    metadata <- meta %>%
      mutate(
        dataset = "GSE207128",
        brain_region = "Amygdala",
        sex = "Male",
        data_type = "single_cell_rnaseq",
        # Use actual condition info if available
        condition = case_when(
          "condition" %in% colnames(meta) & meta$condition == "Naive" ~ "Control",
          "condition" %in% colnames(meta) & meta$condition == "Dependent" ~ "Dependent", 
          "condition" %in% colnames(meta) & meta$condition == "Withdrawal" ~ "Withdrawal",
          str_detect(tolower(sample_id), "naive|control") ~ "Control",
          str_detect(tolower(sample_id), "dependent|dependence") ~ "Dependent",
          str_detect(tolower(sample_id), "withdraw|withdrawal") ~ "Withdrawal",
          TRUE ~ "Unknown"
        ),
        treatment = case_when(
          "condition" %in% colnames(meta) & meta$condition == "Naive" ~ "Control",
          str_detect(tolower(sample_id), "naive|control") ~ "Control",
          TRUE ~ "Morphine_chronic"
        ),
        time_point = case_when(
          "condition" %in% colnames(meta) & meta$condition == "Naive" ~ "Baseline",
          "condition" %in% colnames(meta) & meta$condition == "Dependent" ~ "Chronic_14d",
          "condition" %in% colnames(meta) & meta$condition == "Withdrawal" ~ "Withdrawal_24h",
          str_detect(tolower(sample_id), "naive|control") ~ "Baseline",
          str_detect(tolower(sample_id), "dependent|dependence") ~ "Chronic_14d",
          str_detect(tolower(sample_id), "withdraw|withdrawal") ~ "Withdrawal_24h",
          TRUE ~ "Unknown"
        )
      )
    
  } else if (dataset_name == "GSE289002") {
    # Check if we have DGEList with sample info
    if (is.list(data) && "samples" %in% names(data)) {
      meta <- data$samples
      meta$sample_id <- rownames(meta)
    } else if ("DGEList" %in% class(data) && "samples" %in% names(data)) {
      meta <- data$samples
      meta$sample_id <- rownames(meta)
    } else {
      n_samples <- ifelse(is.matrix(data), ncol(data), length(data))
      meta <- data.frame(
        sample_id = paste0("GSE289002_", 1:n_samples),
        stringsAsFactors = FALSE
      )
    }
    
    metadata <- meta %>%
      mutate(
        dataset = "GSE289002",
        data_type = "bulk_rnaseq",
        # Use actual region info if available
        brain_region = case_when(
          "region" %in% colnames(meta) & meta$region == "PFC" ~ "PFC",
          "region" %in% colnames(meta) & meta$region == "NAc" ~ "NAc", 
          str_detect(tolower(sample_id), "pfc|prefrontal") ~ "PFC",
          str_detect(tolower(sample_id), "nac|nucleus.accumbens") ~ "NAc",
          TRUE ~ "Unknown"
        ),
        # Use actual sex info if available
        sex = case_when(
          "sex" %in% colnames(meta) & meta$sex == "male" ~ "Male",
          "sex" %in% colnames(meta) & meta$sex == "female" ~ "Female",
          str_detect(tolower(sample_id), "male|m") & !str_detect(tolower(sample_id), "female") ~ "Male",
          str_detect(tolower(sample_id), "female|f") & !str_detect(tolower(sample_id), "male") ~ "Female",
          TRUE ~ "Unknown"
        ),
        # Use actual treatment info if available
        condition = case_when(
          "treatment" %in% colnames(meta) & meta$treatment == "Sal" ~ "Control",
          "treatment" %in% colnames(meta) & str_detect(meta$treatment, "Mor") ~ "Morphine",
          str_detect(tolower(sample_id), "control|baseline|sal") ~ "Control",
          TRUE ~ "Morphine"
        ),
        time_point = case_when(
          "treatment" %in% colnames(meta) & meta$treatment == "Sal" ~ "Baseline",
          "treatment" %in% colnames(meta) & meta$treatment == "Mor + 24h" ~ "Withdrawal_1d",
          "treatment" %in% colnames(meta) & meta$treatment == "Mor + 2W" ~ "Withdrawal_14d",
          "treatment" %in% colnames(meta) & meta$treatment == "Chronic mor" ~ "Chronic",
          str_detect(tolower(sample_id), "control|baseline") ~ "Baseline",
          TRUE ~ "Unknown"
        ),
        treatment = case_when(
          "treatment" %in% colnames(meta) & meta$treatment == "Sal" ~ "Control",
          str_detect(tolower(sample_id), "control|baseline") ~ "Control",
          TRUE ~ "Morphine"
        )
      )
  }
  
  return(metadata)
}

# Function to harmonize count matrices
harmonize_count_matrix <- function(data, metadata, dataset_name) {
  cat("Harmonizing count matrix for", dataset_name, "...\n")
  
  # Extract count matrix from different data structures
  if (is.list(data)) {
    if ("counts" %in% names(data)) {
      count_matrix <- data$counts
    } else if ("normalized_counts" %in% names(data)) {
      count_matrix <- round(data$normalized_counts)
    } else if ("seurat_object" %in% names(data)) {
      # Handle Seurat v5 compatibility
      tryCatch({
        count_matrix <- GetAssayData(data$seurat_object, assay = "RNA", layer = "counts")
      }, error = function(e) {
        # Fallback to older Seurat syntax
        count_matrix <- GetAssayData(data$seurat_object, assay = "RNA", slot = "counts")
      })
    } else if ("Seurat" %in% class(data)) {
      # Handle Seurat v5 compatibility
      tryCatch({
        count_matrix <- GetAssayData(data, assay = "RNA", layer = "counts")
      }, error = function(e) {
        # Fallback to older Seurat syntax
        count_matrix <- GetAssayData(data, assay = "RNA", slot = "counts")
      })
    } else if (is.matrix(data[[1]])) {
      count_matrix <- data[[1]]
    } else {
      cat("Could not extract count matrix from data structure\n")
      return(NULL)
    }
  } else if (is.matrix(data)) {
    count_matrix <- data
  } else if ("Seurat" %in% class(data)) {
    # Handle Seurat v5 compatibility with more robust approach
    cat("Processing Seurat object...\n")
    
    # Check available assays
    available_assays <- names(data@assays)
    cat("Available assays:", paste(available_assays, collapse = ", "), "\n")
    
    # Try different approaches to extract counts
    count_matrix <- NULL
    
    # Approach 1: Try JoinLayers() first for Seurat v5
    tryCatch({
      # Try to join layers first
      data_joined <- JoinLayers(data, assay = "RNA")
      count_matrix <- GetAssayData(data_joined, assay = "RNA", layer = "counts")
      cat("✓ Extracted counts after joining layers\n")
    }, error = function(e) {
      cat("JoinLayers method failed:", e$message, "\n")
    })
    
    # Approach 2: Try GetAssayData with layer (Seurat v5)
    if (is.null(count_matrix)) {
      tryCatch({
        count_matrix <- GetAssayData(data, assay = "RNA", layer = "counts")
        cat("✓ Extracted counts using layer method\n")
      }, error = function(e) {
        cat("Layer method failed:", e$message, "\n")
      })
    }
    
    # Approach 3: Try GetAssayData with slot (Seurat v4)
    if (is.null(count_matrix)) {
      tryCatch({
        count_matrix <- GetAssayData(data, assay = "RNA", slot = "counts")
        cat("✓ Extracted counts using slot method\n")
      }, error = function(e) {
        cat("Slot method failed:", e$message, "\n")
      })
    }
    
    # Approach 4: Direct access to assay data
    if (is.null(count_matrix) && "RNA" %in% available_assays) {
      tryCatch({
        assay_obj <- data@assays$RNA
        
        # For Seurat v5, try to access layers directly
        if ("layers" %in% slotNames(assay_obj)) {
          layers_available <- names(assay_obj@layers)
          cat("Available layers:", paste(layers_available, collapse = ", "), "\n")
          
          # Try to get counts layer - combine all layers if multiple
          count_layers <- layers_available[str_detect(layers_available, "counts")]
          if (length(count_layers) > 0) {
            if (length(count_layers) == 1) {
              count_matrix <- assay_obj@layers[[count_layers[1]]]
              cat("✓ Extracted counts from single layer:", count_layers[1], "\n")
            } else {
              # Combine multiple count layers
              cat("Combining", length(count_layers), "count layers...\n")
              layer_matrices <- map(count_layers, ~ assay_obj@layers[[.x]])
              count_matrix <- do.call(cbind, layer_matrices)
              cat("✓ Combined counts from multiple layers\n")
            }
          } else if ("data" %in% layers_available) {
            count_matrix <- assay_obj@layers$data
            cat("Warning: Using 'data' layer instead of 'counts'\n")
          } else if (length(layers_available) > 0) {
            # Use the first available layer
            count_matrix <- assay_obj@layers[[layers_available[1]]]
            cat("Warning: Using", layers_available[1], "layer as count matrix\n")
          }
        }
        
        # Fallback to traditional slots if layers didn't work
        if (is.null(count_matrix)) {
          if ("counts" %in% slotNames(assay_obj)) {
            count_matrix <- assay_obj@counts
            cat("✓ Extracted counts from direct assay access\n")
          } else if ("data" %in% slotNames(assay_obj)) {
            count_matrix <- assay_obj@data
            cat("Warning: Using 'data' slot instead of 'counts'\n")
          }
        }
      }, error = function(e) {
        cat("Direct access failed:", e$message, "\n")
      })
    }
    
    # Approach 5: Try other assay names if RNA doesn't work
    if (is.null(count_matrix)) {
      for (assay_name in available_assays) {
        tryCatch({
          count_matrix <- GetAssayData(data, assay = assay_name, layer = "counts")
          cat("✓ Extracted counts from", assay_name, "assay\n")
          break
        }, error = function(e) {
          tryCatch({
            count_matrix <- GetAssayData(data, assay = assay_name, slot = "counts")
            cat("✓ Extracted counts from", assay_name, "assay (slot method)\n")
            break
          }, error = function(e2) {
            # Continue to next assay
          })
        })
      }
    }
    
    if (is.null(count_matrix)) {
      cat("Could not extract count matrix from Seurat object\n")
      cat("Object structure:\n")
      print(str(data, max.level = 2))
      return(NULL)
    }
    
  } else if ("DGEList" %in% class(data)) {
    # Handle DGEList objects from edgeR
    count_matrix <- data$counts
  } else {
    # Try to convert data.frame to matrix
    if (is.data.frame(data)) {
      # Assume first column is gene names
      gene_names <- data[[1]]
      count_matrix <- as.matrix(data[, -1])
      rownames(count_matrix) <- gene_names
    } else {
      cat("Unknown data structure for", dataset_name, ":", class(data), "\n")
      return(NULL)
    }
  }
  
  # Ensure count matrix is numeric and handle sparse matrices
  if (inherits(count_matrix, "dgCMatrix") || inherits(count_matrix, "sparseMatrix")) {
    # Keep as sparse matrix but ensure it's numeric
    count_matrix <- as(count_matrix, "dgCMatrix")
  } else if (!is.numeric(count_matrix)) {
    count_matrix <- apply(count_matrix, 2, as.numeric)
  }
  
  # Check for valid count matrix
  if (is.null(count_matrix) || nrow(count_matrix) == 0 || ncol(count_matrix) == 0) {
    cat("Invalid count matrix extracted for", dataset_name, "\n")
    return(NULL)
  }
  
  cat("Extracted count matrix with", nrow(count_matrix), "genes and", ncol(count_matrix), "samples\n")
  
  # Standardize gene symbols
  gene_names_std <- standardize_gene_symbols(rownames(count_matrix), species = "mouse")
  rownames(count_matrix) <- gene_names_std
  
  # Remove duplicate genes (keep the one with highest expression)
  if (any(duplicated(gene_names_std))) {
    cat("Removing", sum(duplicated(gene_names_std)), "duplicate genes\n")
    
    # Handle sparse matrices differently
    if (inherits(count_matrix, "dgCMatrix")) {
      gene_means <- Matrix::rowMeans(count_matrix)
    } else {
      gene_means <- rowMeans(count_matrix)
    }
    
    keep_genes <- !duplicated(gene_names_std) | 
                  gene_means == ave(gene_means, gene_names_std, FUN = max)
    count_matrix <- count_matrix[keep_genes, ]
  }
  
  # Filter out genes with very low expression - be less aggressive for single-cell data
  min_cells_expr <- if (dataset_name == "GSE207128") {
    max(1, ncol(count_matrix) * 0.001)  # 0.1% for single-cell
  } else {
    max(1, ncol(count_matrix) * 0.1)    # 10% for bulk
  }
  
  if (inherits(count_matrix, "dgCMatrix")) {
    gene_keep <- Matrix::rowSums(count_matrix >= 1) >= min_cells_expr
  } else {
    gene_keep <- rowSums(count_matrix >= 1) >= min_cells_expr
  }
  
  count_matrix <- count_matrix[gene_keep, ]
  
  cat("Final count matrix dimensions:", nrow(count_matrix), "genes x", 
      ncol(count_matrix), "samples\n")
  
  return(count_matrix)
}

# ==============================================================================
# LOAD AND PROCESS MOUSE DATASETS
# ==============================================================================

mouse_datasets <- list()

# Dataset paths
datasets_info <- data.frame(
  name = c("GSE118918", "GSE207128", "GSE289002"),
  path = c(
    file.path(mouse_data_dir, "GSE118918"),
    file.path(mouse_data_dir, "GSE207128"),
    file.path(mouse_data_dir, "GSE289002")
  ),
  stringsAsFactors = FALSE
)

# Process each dataset
for (i in 1:nrow(datasets_info)) {
  dataset_name <- datasets_info$name[i]
  dataset_path <- datasets_info$path[i]
  
  cat("\n", strrep("=", 60), "\n")
  cat("PROCESSING", dataset_name, "\n")
  cat(strrep("=", 60), "\n")
  
  if (!dir.exists(dataset_path)) {
    cat("Directory does not exist:", dataset_path, "\n")
    cat("Checking for alternative naming...\n")
    
    # Try alternative dataset names
    parent_dir <- dirname(dataset_path)
    available_dirs <- list.dirs(parent_dir, recursive = FALSE, full.names = FALSE)
    possible_matches <- available_dirs[str_detect(tolower(available_dirs), 
                                                 tolower(str_extract(dataset_name, "GSE\\d+")))]
    
    if (length(possible_matches) > 0) {
      cat("Found possible matches:", paste(possible_matches, collapse = ", "), "\n")
      dataset_path <- file.path(parent_dir, possible_matches[1])
      cat("Using directory:", dataset_path, "\n")
    } else {
      cat("No matching directories found. Skipping", dataset_name, "\n")
      next
    }
  }
  
  # Load dataset
  loaded_data <- load_mouse_dataset(dataset_path, dataset_name)
  
  if (is.null(loaded_data)) {
    cat("Could not load data for", dataset_name, "\n")
    next
  }
  
  # Extract metadata
  metadata <- extract_metadata(loaded_data$data, dataset_name)
  
  # Harmonize count matrix
  count_matrix <- harmonize_count_matrix(loaded_data$data, metadata, dataset_name)
  
  if (is.null(count_matrix)) {
    cat("Could not harmonize count matrix for", dataset_name, "\n")
    next
  }
  
  # Ensure metadata matches count matrix
  if (ncol(count_matrix) != nrow(metadata)) {
    cat("Mismatch between count matrix and metadata dimensions\n")
    cat("Count matrix samples:", ncol(count_matrix), "\n")
    cat("Metadata samples:", nrow(metadata), "\n")
    
    # Try to fix by using column names
    if (!is.null(colnames(count_matrix)) && all(colnames(count_matrix) %in% metadata$sample_id)) {
      metadata <- metadata[match(colnames(count_matrix), metadata$sample_id), ]
      cat("✓ Fixed metadata-matrix mismatch using column names\n")
    } else {
      # Create minimal metadata
      metadata <- data.frame(
        sample_id = colnames(count_matrix) %||% paste0(dataset_name, "_", 1:ncol(count_matrix)),
        dataset = dataset_name,
        stringsAsFactors = FALSE
      )
      # Add basic annotations
      metadata <- extract_metadata(list(metadata = metadata), dataset_name)
      cat("✓ Created new metadata for", nrow(metadata), "samples\n")
    }
  }
  
  # Store processed dataset
  mouse_datasets[[dataset_name]] <- list(
    counts = count_matrix,
    metadata = metadata,
    file_info = loaded_data
  )
  
  cat("✓ Successfully processed", dataset_name, "\n")
  cat("  Dimensions:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n")
  cat("  Conditions:", paste(unique(metadata$condition), collapse = ", "), "\n")
  cat("  Brain regions:", paste(unique(metadata$brain_region), collapse = ", "), "\n")
  
  # Debug: Inspect sample metadata more deeply
  cat("\n  Detailed metadata inspection:\n")
  cat("    Sample IDs (first 5):", paste(head(metadata$sample_id, 5), collapse = ", "), "\n")
  cat("    Unique conditions:", paste(unique(metadata$condition), collapse = ", "), "\n")
  cat("    Unique treatments:", paste(unique(metadata$treatment), collapse = ", "), "\n")
  cat("    Unique time points:", paste(unique(metadata$time_point), collapse = ", "), "\n")
  cat("    Sex distribution:", paste(names(table(metadata$sex)), "=", table(metadata$sex), collapse = ", "), "\n")
  
  # Debug: Inspect gene names more deeply
  cat("  Gene name inspection:\n")
  cat("    First 10 genes:", paste(head(rownames(count_matrix), 10), collapse = ", "), "\n")
  cat("    Last 5 genes:", paste(tail(rownames(count_matrix), 5), collapse = ", "), "\n")
  
  # Check for different gene ID types
  ensembl_count <- sum(str_detect(rownames(count_matrix), "^ENSMUSG"))
  symbol_count <- sum(str_detect(rownames(count_matrix), "^[A-Za-z]"))
  numeric_count <- sum(str_detect(rownames(count_matrix), "^[0-9]"))
  
  cat("    Gene ID types:\n")
  cat("      ENSEMBL IDs (ENSMUSG...):", ensembl_count, "\n")
  cat("      Gene symbols (letters):", symbol_count, "\n")
  cat("      Numeric/other IDs:", numeric_count, "\n")
  
  # Show distribution of gene expression
  total_counts <- sum(count_matrix)
  zero_genes <- sum(rowSums(count_matrix) == 0)
  cat("    Expression summary:\n")
  cat("      Total counts:", format(total_counts, scientific = TRUE), "\n")
  cat("      Genes with zero counts:", zero_genes, "/", nrow(count_matrix), "\n")
  cat("      Mean counts per gene:", round(mean(rowSums(count_matrix)), 2), "\n")
}

# ==============================================================================
# DETAILED METADATA AND GENE ANALYSIS
# ==============================================================================

cat("\n", strrep("=", 60), "\n")
cat("DETAILED CROSS-DATASET ANALYSIS\n")
cat(strrep("=", 60), "\n")

# Deep dive into metadata patterns
cat("Analyzing metadata patterns across datasets:\n")
for (dataset_name in names(mouse_datasets)) {
  metadata <- mouse_datasets[[dataset_name]]$metadata
  cat("\n", dataset_name, "metadata analysis:\n")
  
  # Look for actual sample names from the original data
  if ("file_info" %in% names(mouse_datasets[[dataset_name]])) {
    original_data <- mouse_datasets[[dataset_name]]$file_info$data
    if ("DGEList" %in% class(original_data) && "samples" %in% names(original_data)) {
      cat("  Original sample info from DGEList:\n")
      original_samples <- original_data$samples
      if (nrow(original_samples) > 0) {
        cat("    Sample columns:", paste(colnames(original_samples), collapse = ", "), "\n")
        cat("    First few sample entries:\n")
        print(head(original_samples, 3))
        
        # Try to extract meaningful condition/treatment info
        for (col in colnames(original_samples)) {
          unique_vals <- unique(original_samples[[col]])
          if (length(unique_vals) <= 10 && length(unique_vals) > 1) {
            cat("    ", col, "values:", paste(unique_vals, collapse = ", "), "\n")
          }
        }
      }
    } else if ("Seurat" %in% class(original_data)) {
      cat("  Original metadata from Seurat object:\n")
      original_meta <- original_data@meta.data
      if (nrow(original_meta) > 0) {
        cat("    Metadata columns:", paste(colnames(original_meta), collapse = ", "), "\n")
        cat("    Number of cells:", nrow(original_meta), "\n")
        
        # Look for treatment/condition info in Seurat metadata
        condition_cols <- colnames(original_meta)[str_detect(tolower(colnames(original_meta)), 
                                                           "condition|treatment|group|sample")]
        if (length(condition_cols) > 0) {
          cat("    Potential condition columns:", paste(condition_cols, collapse = ", "), "\n")
          for (col in condition_cols[1:min(3, length(condition_cols))]) {
            unique_vals <- unique(original_meta[[col]])
            if (length(unique_vals) <= 20) {
              cat("      ", col, ":", paste(unique_vals, collapse = ", "), "\n")
            }
          }
        }
      }
    }
  }
}

# Deep dive into gene overlap issues
cat("\nAnalyzing gene overlap issues:\n")
all_genes_detailed <- map(mouse_datasets, ~ rownames(.x$counts))
names(all_genes_detailed) <- names(mouse_datasets)

for (i in 1:length(all_genes_detailed)) {
  dataset_name <- names(all_genes_detailed)[i]
  genes <- all_genes_detailed[[i]]
  
  cat("\n", dataset_name, "gene analysis:\n")
  cat("  Total genes:", length(genes), "\n")
  
  # Analyze gene name patterns
  ensembl_pattern <- str_detect(genes, "^ENSMUSG[0-9]+")
  symbol_pattern <- str_detect(genes, "^[A-Za-z][A-Za-z0-9-]*$")
  rik_pattern <- str_detect(genes, "Rik$")
  mt_pattern <- str_detect(genes, "^[Mm][Tt]-")
  rp_pattern <- str_detect(genes, "^[Rr][pP][ls][0-9]")
  
  cat("  Gene name patterns:\n")
  cat("    ENSEMBL IDs:", sum(ensembl_pattern), "\n")
  cat("    Gene symbols:", sum(symbol_pattern), "\n")
  cat("    Rik genes:", sum(rik_pattern), "\n")
  cat("    Mitochondrial genes:", sum(mt_pattern), "\n")
  cat("    Ribosomal proteins:", sum(rp_pattern), "\n")
  
  # Show examples of each type
  if (sum(ensembl_pattern) > 0) {
    cat("    ENSEMBL examples:", paste(head(genes[ensembl_pattern], 3), collapse = ", "), "\n")
  }
  if (sum(symbol_pattern) > 0) {
    cat("    Symbol examples:", paste(head(genes[symbol_pattern], 3), collapse = ", "), "\n")
  }
  
  # Look for common marker genes
  common_genes_check <- c("Gapdh", "Actb", "Tubb3", "Slc6a3", "Th", "Drd1", "Drd2", "Oprm1")
  common_genes_present <- genes[tolower(genes) %in% tolower(common_genes_check)]
  if (length(common_genes_present) > 0) {
    cat("    Common marker genes found:", paste(common_genes_present, collapse = ", "), "\n")
  }
}

# Try alternative gene matching strategies
cat("\nTrying alternative gene matching strategies:\n")

# Strategy 1: Case-insensitive matching
cat("1. Case-insensitive gene matching:\n")
genes_lower <- map(all_genes_detailed, ~ tolower(.x))
for (i in 1:(length(genes_lower)-1)) {
  for (j in (i+1):length(genes_lower)) {
    dataset1 <- names(genes_lower)[i]
    dataset2 <- names(genes_lower)[j]
    overlap <- length(intersect(genes_lower[[i]], genes_lower[[j]]))
    cat("  ", dataset1, "vs", dataset2, ":", overlap, "genes (case-insensitive)\n")
    
    if (overlap > 0) {
      common_lower <- intersect(genes_lower[[i]], genes_lower[[j]])
      cat("    Examples:", paste(head(common_lower, 5), collapse = ", "), "\n")
    }
  }
}

# Strategy 2: Try ENSEMBL to symbol conversion if biomaRt is available
cat("\n2. Attempting gene ID conversion:\n")
tryCatch({
  if (requireNamespace("biomaRt", quietly = TRUE)) {
    # Get mouse annotation
    mart <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    
    for (dataset_name in names(all_genes_detailed)) {
      genes <- all_genes_detailed[[dataset_name]]
      ensembl_genes <- genes[str_detect(genes, "^ENSMUSG")]
      
      if (length(ensembl_genes) > 0 && length(ensembl_genes) < 1000) {  # Don't overwhelm biomaRt
        cat("  Converting ENSEMBL IDs for", dataset_name, "(", length(ensembl_genes), "genes)...\n")
        
        gene_map <- biomaRt::getBM(
          attributes = c("ensembl_gene_id", "external_gene_name"),
          filters = "ensembl_gene_id",
          values = ensembl_genes,
          mart = mart
        )
        
        if (nrow(gene_map) > 0) {
          cat("    Converted", nrow(gene_map), "ENSEMBL IDs to symbols\n")
          cat("    Examples:", paste(head(gene_map$external_gene_name, 5), collapse = ", "), "\n")
          
          # Store conversion for later use
          assign(paste0(dataset_name, "_gene_conversion"), gene_map, envir = .GlobalEnv)
        }
      }
    }
  } else {
    cat("  biomaRt not available for gene ID conversion\n")
  }
}, error = function(e) {
  cat("  Gene ID conversion failed:", e$message, "\n")
})

# ==============================================================================
# GENE ID CONVERSION AND HARMONIZATION
# ==============================================================================

cat("\n", strrep("=", 60), "\n")
cat("GENE ID CONVERSION AND HARMONIZATION\n")
cat(strrep("=", 60), "\n")

# Function to convert ENSEMBL IDs to gene symbols
convert_ensembl_to_symbols <- function(ensembl_ids, max_batch_size = 500) {
  cat("Converting", length(ensembl_ids), "ENSEMBL IDs to gene symbols...\n")
  
  if (length(ensembl_ids) == 0) {
    return(data.frame(ensembl_gene_id = character(), external_gene_name = character()))
  }
  
  # Try biomaRt first
  gene_map <- NULL
  tryCatch({
    if (requireNamespace("biomaRt", quietly = TRUE)) {
      mart <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
      
      # Process in batches to avoid overwhelming biomaRt
      all_conversions <- list()
      n_batches <- ceiling(length(ensembl_ids) / max_batch_size)
      
      for (i in 1:n_batches) {
        start_idx <- (i-1) * max_batch_size + 1
        end_idx <- min(i * max_batch_size, length(ensembl_ids))
        batch_ids <- ensembl_ids[start_idx:end_idx]
        
        cat("  Processing batch", i, "of", n_batches, "(", length(batch_ids), "genes)...\n")
        
        batch_map <- biomaRt::getBM(
          attributes = c("ensembl_gene_id", "external_gene_name"),
          filters = "ensembl_gene_id",
          values = batch_ids,
          mart = mart
        )
        
        if (nrow(batch_map) > 0) {
          all_conversions[[i]] <- batch_map
        }
        
        # Small delay between batches
        Sys.sleep(0.5)
      }
      
      if (length(all_conversions) > 0) {
        gene_map <- do.call(rbind, all_conversions)
        gene_map <- gene_map[gene_map$external_gene_name != "", ]  # Remove empty names
        cat("✓ Successfully converted", nrow(gene_map), "ENSEMBL IDs via biomaRt\n")
      }
    }
  }, error = function(e) {
    cat("biomaRt conversion failed:", e$message, "\n")
  })
  
  # If biomaRt failed, try using a local annotation file or package
  if (is.null(gene_map) || nrow(gene_map) == 0) {
    cat("Trying alternative gene ID conversion methods...\n")
    
    # Try using org.Mm.eg.db if available
    tryCatch({
      if (requireNamespace("org.Mm.eg.db", quietly = TRUE) && 
          requireNamespace("AnnotationDbi", quietly = TRUE)) {
        library(org.Mm.eg.db)
        library(AnnotationDbi)
        
        # Convert ENSEMBL to SYMBOL
        symbols <- AnnotationDbi::mapIds(org.Mm.eg.db, 
                                       keys = ensembl_ids,
                                       column = "SYMBOL",
                                       keytype = "ENSEMBL",
                                       multiVals = "first")
        
        # Create data frame, removing NAs
        valid_conversions <- !is.na(symbols)
        if (sum(valid_conversions) > 0) {
          gene_map <- data.frame(
            ensembl_gene_id = names(symbols)[valid_conversions],
            external_gene_name = as.character(symbols[valid_conversions]),
            stringsAsFactors = FALSE
          )
          cat("✓ Successfully converted", nrow(gene_map), "ENSEMBL IDs via org.Mm.eg.db\n")
        }
      }
    }, error = function(e) {
      cat("org.Mm.eg.db conversion failed:", e$message, "\n")
    })
  }
  
  # Return empty data frame if all methods failed
  if (is.null(gene_map)) {
    cat("All gene ID conversion methods failed\n")
    gene_map <- data.frame(ensembl_gene_id = character(), external_gene_name = character())
  }
  
  return(gene_map)
}

# Convert ENSEMBL IDs in datasets that need it
converted_datasets <- list()

for (dataset_name in names(mouse_datasets)) {
  dataset <- mouse_datasets[[dataset_name]]
  gene_names <- rownames(dataset$counts)
  
  # Check if this dataset has ENSEMBL IDs
  ensembl_genes <- gene_names[str_detect(gene_names, "^ENSMUSG")]
  
  if (length(ensembl_genes) > 100) {  # Only convert if substantial number of ENSEMBL IDs
    cat("\nConverting ENSEMBL IDs in", dataset_name, "...\n")
    
    # Get conversion mapping
    gene_conversion <- convert_ensembl_to_symbols(ensembl_genes)
    
    if (nrow(gene_conversion) > 0) {
      # Create new gene names
      new_gene_names <- gene_names
      
      # Map ENSEMBL IDs to symbols where possible
      for (i in 1:nrow(gene_conversion)) {
        ensembl_id <- gene_conversion$ensembl_gene_id[i]
        symbol <- gene_conversion$external_gene_name[i]
        
        if (symbol != "" && !is.na(symbol)) {
          new_gene_names[new_gene_names == ensembl_id] <- symbol
        }
      }
      
      # Update the dataset with converted gene names
      converted_count_matrix <- dataset$counts
      rownames(converted_count_matrix) <- new_gene_names
      
      # Remove genes that couldn't be converted (still ENSEMBL IDs)
      genes_to_keep <- !str_detect(new_gene_names, "^ENSMUSG")
      converted_count_matrix <- converted_count_matrix[genes_to_keep, ]
      
      # Handle duplicated gene symbols by keeping highest expressed
      if (any(duplicated(rownames(converted_count_matrix)))) {
        cat("  Handling", sum(duplicated(rownames(converted_count_matrix))), "duplicated gene symbols...\n")
        
        if (inherits(converted_count_matrix, "dgCMatrix")) {
          gene_means <- Matrix::rowMeans(converted_count_matrix)
        } else {
          gene_means <- rowMeans(converted_count_matrix)
        }
        
        keep_genes <- !duplicated(rownames(converted_count_matrix)) | 
                      gene_means == ave(gene_means, rownames(converted_count_matrix), FUN = max)
        converted_count_matrix <- converted_count_matrix[keep_genes, ]
      }
      
      # Store converted dataset
      converted_datasets[[dataset_name]] <- list(
        counts = converted_count_matrix,
        metadata = dataset$metadata,
        file_info = dataset$file_info,
        gene_conversion = gene_conversion
      )
      
      cat("✓ Converted", dataset_name, ":", nrow(gene_conversion), "ENSEMBL->symbol mappings\n")
      cat("  Original genes:", nrow(dataset$counts), "-> Converted genes:", nrow(converted_count_matrix), "\n")
      
    } else {
      cat("No successful gene conversions for", dataset_name, "\n")
      converted_datasets[[dataset_name]] <- dataset
    }
  } else {
    cat("Dataset", dataset_name, "has few/no ENSEMBL IDs, keeping original\n")
    converted_datasets[[dataset_name]] <- dataset
  }
}

# Update mouse_datasets with converted versions
mouse_datasets <- converted_datasets

# Re-analyze gene overlap after conversion
cat("\nRe-analyzing gene overlap after ID conversion:\n")
all_genes_converted <- map(mouse_datasets, ~ rownames(.x$counts))
names(all_genes_converted) <- names(mouse_datasets)

for (i in 1:(length(all_genes_converted)-1)) {
  for (j in (i+1):length(all_genes_converted)) {
    dataset1 <- names(all_genes_converted)[i]
    dataset2 <- names(all_genes_converted)[j]
    overlap <- length(intersect(all_genes_converted[[i]], all_genes_converted[[j]]))
    union_size <- length(union(all_genes_converted[[i]], all_genes_converted[[j]]))
    cat("  ", dataset1, "vs", dataset2, ":", overlap, "/", union_size, 
        "genes (", round(100*overlap/union_size, 1), "% overlap)\n")
    
    if (overlap > 0 && i == 1 && j == 2) {
      common_genes_example <- intersect(all_genes_converted[[i]], all_genes_converted[[j]])
      cat("    Common genes example:", paste(head(common_genes_example, 5), collapse = ", "), "\n")
    }
  }
}

# ==============================================================================
# CREATE HARMONIZED DATASET STRUCTURE
# ==============================================================================

cat("\n", strrep("=", 50), "\n")
cat("CREATING HARMONIZED DATASET STRUCTURE\n")
cat(strrep("=", 50), "\n")

# Combine all metadata
all_metadata <- map_dfr(mouse_datasets, ~ .x$metadata, .id = "source_dataset")

# Debug metadata structure
cat("Debug: Checking all_metadata structure:\n")
cat("Class:", class(all_metadata), "\n")
cat("Dimensions:", dim(all_metadata), "\n")
cat("Column names:", paste(colnames(all_metadata), collapse = ", "), "\n")
if (nrow(all_metadata) > 0) {
  cat("First few rows:\n")
  print(head(all_metadata, 3))
}

# Find common genes across all datasets
all_genes <- map(mouse_datasets, ~ rownames(.x$counts))
common_genes <- Reduce(intersect, all_genes)

cat("Common genes across all datasets after conversion:", length(common_genes), "\n")

if (length(common_genes) > 0) {
  cat("✓ Found", length(common_genes), "common genes! Gene ID conversion was successful.\n")
  cat("Example common genes:", paste(head(common_genes, 10), collapse = ", "), "\n")
}

# Handle case where there are no common genes
if (length(common_genes) == 0) {
  cat("Warning: No common genes found across all datasets.\n")
  cat("Gene overlap between datasets:\n")
  
  # Show pairwise overlaps
  dataset_names <- names(all_genes)
  for (i in 1:(length(dataset_names)-1)) {
    for (j in (i+1):length(dataset_names)) {
      overlap <- length(intersect(all_genes[[i]], all_genes[[j]]))
      cat("  ", dataset_names[i], "vs", dataset_names[j], ":", overlap, "genes\n")
    }
  }
  
  # Use union of all genes for harmonization
  common_genes <- Reduce(union, all_genes)
  cat("Using union of all genes:", length(common_genes), "genes\n")
  
  # Filter each dataset to only include genes present in that dataset
  mouse_datasets_harmonized <- map(mouse_datasets, function(dataset) {
    # Get genes present in this dataset
    dataset_genes <- rownames(dataset$counts)
    genes_to_keep <- intersect(common_genes, dataset_genes)
    
    cat("Processing dataset with", length(dataset_genes), "genes,", length(genes_to_keep), "genes to keep\n")
    
    # Determine if we should use sparse matrices based on data size
    n_genes <- length(common_genes)
    n_samples <- ncol(dataset$counts)
    use_sparse <- inherits(dataset$counts, "dgCMatrix") || (n_genes * n_samples > 1e6)
    
    if (use_sparse) {
      cat("Using sparse matrix operations for memory efficiency\n")
      cat("Creating sparse matrix of size:", n_genes, "x", n_samples, "...\n")
      
      # Create sparse matrix with all common genes
      tryCatch({
        counts_harmonized <- Matrix::sparseMatrix(
          i = integer(0), 
          j = integer(0), 
          x = numeric(0),
          dims = c(length(common_genes), ncol(dataset$counts)),
          dimnames = list(common_genes, colnames(dataset$counts))
        )
        cat("✓ Successfully created sparse matrix\n")
      }, error = function(e) {
        cat("Error creating sparse matrix:", e$message, "\n")
        stop("Failed to create sparse matrix")
      })
      
      # Find indices of genes to keep in the harmonized matrix
      gene_indices <- match(genes_to_keep, common_genes)
      cat("Filling sparse matrix with", length(gene_indices), "genes...\n")
      
      # Fill in the actual data - keep sparse format
      if (inherits(dataset$counts, "dgCMatrix")) {
        # Direct sparse matrix assignment
        tryCatch({
          counts_harmonized[gene_indices, ] <- dataset$counts[genes_to_keep, , drop = FALSE]
          cat("✓ Successfully filled sparse matrix from sparse source\n")
        }, error = function(e) {
          cat("Error in sparse matrix assignment:", e$message, "\n")
          # Try alternative approach with cbind/rbind
          cat("Trying alternative sparse matrix filling approach...\n")
          # This is a more memory-efficient approach for very large matrices
          for (i in seq_along(gene_indices)) {
            if (i %% 1000 == 0) cat("  Processed", i, "of", length(gene_indices), "genes\n")
            counts_harmonized[gene_indices[i], ] <- dataset$counts[genes_to_keep[i], , drop = FALSE]
          }
          cat("✓ Successfully filled sparse matrix using iterative approach\n")
        })
      } else {
        # Convert regular matrix subset to sparse and assign
        cat("Converting dense matrix subset to sparse...\n")
        subset_sparse <- Matrix::Matrix(dataset$counts[genes_to_keep, , drop = FALSE], sparse = TRUE)
        counts_harmonized[gene_indices, ] <- subset_sparse
        cat("✓ Successfully filled sparse matrix from dense source\n")
      }
      
    } else {
      cat("Using dense matrix operations\n")
      
      # Create a regular matrix with all common genes, filling missing genes with zeros
      counts_harmonized <- matrix(0, 
                                nrow = length(common_genes), 
                                ncol = ncol(dataset$counts),
                                dimnames = list(common_genes, colnames(dataset$counts)))
      
      # Fill in the actual data - handle sparse matrices properly
      if (inherits(dataset$counts, "dgCMatrix")) {
        # For sparse matrices, convert subset to regular matrix first
        subset_data <- as.matrix(dataset$counts[genes_to_keep, , drop = FALSE])
        counts_harmonized[genes_to_keep, ] <- subset_data
      } else {
        # For regular matrices
        counts_harmonized[genes_to_keep, ] <- dataset$counts[genes_to_keep, , drop = FALSE]
      }
    }
    
    cat("Matrix harmonization completed. Starting normalization...\n")
    
    # Perform basic normalization for visualization
    # Only include genes that have some expression
    if (inherits(counts_harmonized, "dgCMatrix")) {
      cat("Calculating gene expression summary for sparse matrix...\n")
      genes_with_expression <- Matrix::rowSums(counts_harmonized) > 0
    } else {
      genes_with_expression <- rowSums(counts_harmonized) > 0
    }
    
    cat("Genes with expression:", sum(genes_with_expression), "out of", length(genes_with_expression), "\n")
    
    counts_for_norm <- counts_harmonized[genes_with_expression, , drop = FALSE]
    
    if (nrow(counts_for_norm) > 0) {
      # For very large sparse matrices, use a sample for normalization to save memory
      if (inherits(counts_for_norm, "dgCMatrix") && ncol(counts_for_norm) > 10000) {
        cat("Large dataset detected - using sampling for normalization factors\n")
        
        # Sample maximum 5000 cells for normalization factor calculation
        sample_size <- min(5000, ncol(counts_for_norm))
        sample_indices <- sort(sample(ncol(counts_for_norm), sample_size))
        cat("Using", sample_size, "cells for normalization factor calculation\n")
        
        # Calculate normalization on sample
        cat("Converting sample to dense matrix for normalization...\n")
        counts_sample <- as.matrix(counts_for_norm[, sample_indices])
        cat("Creating DGEList for normalization...\n")
        dge_sample <- edgeR::DGEList(counts = counts_sample)
        dge_sample <- edgeR::calcNormFactors(dge_sample, method = "TMM")
        cat("✓ Calculated normalization factors\n")
        
        # Apply normalization factors to full dataset
        # Create size factors for all samples (interpolate/extend)
        cat("Calculating library sizes for all samples...\n")
        all_lib_sizes <- Matrix::colSums(counts_for_norm)
        sample_lib_sizes <- all_lib_sizes[sample_indices]
        norm_factors <- dge_sample$samples$norm.factors
        
        # Estimate norm factors for all samples based on library size
        all_norm_factors <- rep(1, ncol(counts_for_norm))
        all_norm_factors[sample_indices] <- norm_factors
        
        # For non-sampled cells, use median norm factor
        non_sampled <- setdiff(1:ncol(counts_for_norm), sample_indices)
        if (length(non_sampled) > 0) {
          all_norm_factors[non_sampled] <- median(norm_factors)
        }
        
        # Calculate CPM manually for sparse matrices
        cat("Calculating CPM for sparse matrix...\n")
        effective_lib_sizes <- all_lib_sizes * all_norm_factors
        
        # Use more memory-efficient CPM calculation
        cpm_normalized <- counts_for_norm
        for (j in 1:ncol(cpm_normalized)) {
          if (j %% 5000 == 0) cat("  Processed", j, "of", ncol(cpm_normalized), "samples for CPM\n")
          cpm_normalized[, j] <- cpm_normalized[, j] / effective_lib_sizes[j] * 1e6
        }
        cat("Calculating log2(CPM + 1)...\n")
        log_cpm <- log2(cpm_normalized + 1)
        cat("✓ Completed normalization for large sparse matrix\n")
        
      } else {
        # Standard normalization for smaller datasets
        cat("Standard normalization for smaller dataset...\n")
        if (inherits(counts_for_norm, "dgCMatrix")) {
          cat("Converting sparse matrix to dense for standard normalization...\n")
          counts_for_norm <- as.matrix(counts_for_norm)
        }
        
        # Create DGEList and normalize
        dge <- edgeR::DGEList(counts = counts_for_norm)
        dge <- edgeR::calcNormFactors(dge, method = "TMM")
        cpm_normalized <- edgeR::cpm(dge, log = FALSE)
        log_cpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
        cat("✓ Completed standard normalization\n")
      }
      
      # Expand back to full gene set - keep sparse if original was sparse
      cat("Expanding normalized data back to full gene set...\n")
      if (use_sparse) {
        # Create sparse matrices for normalized data
        cpm_full <- Matrix::sparseMatrix(
          i = integer(0), j = integer(0), x = numeric(0),
          dims = c(length(common_genes), ncol(counts_harmonized)),
          dimnames = list(common_genes, colnames(counts_harmonized))
        )
        log_cpm_full <- Matrix::sparseMatrix(
          i = integer(0), j = integer(0), x = numeric(0),
          dims = c(length(common_genes), ncol(counts_harmonized)),
          dimnames = list(common_genes, colnames(counts_harmonized))
        )
        
        # Fill in normalized data
        expressed_gene_indices <- which(genes_with_expression)
        cat("Filling", length(expressed_gene_indices), "expressed genes into full sparse matrices...\n")
        if (inherits(cpm_normalized, "dgCMatrix")) {
          cpm_full[expressed_gene_indices, ] <- cpm_normalized
          log_cpm_full[expressed_gene_indices, ] <- log_cpm
        } else {
          cpm_full[expressed_gene_indices, ] <- Matrix::Matrix(cpm_normalized, sparse = TRUE)
          log_cpm_full[expressed_gene_indices, ] <- Matrix::Matrix(log_cpm, sparse = TRUE)
        }
        cat("✓ Completed sparse matrix expansion\n")
        
      } else {
        # Create regular matrices for normalized data
        cpm_full <- matrix(0, nrow = length(common_genes), ncol = ncol(counts_harmonized))
        rownames(cpm_full) <- common_genes
        colnames(cpm_full) <- colnames(counts_harmonized)
        cpm_full[rownames(cpm_normalized), ] <- as.matrix(cpm_normalized)
        
        log_cpm_full <- matrix(0, nrow = length(common_genes), ncol = ncol(counts_harmonized))
        rownames(log_cpm_full) <- common_genes
        colnames(log_cpm_full) <- colnames(counts_harmonized)
        log_cpm_full[rownames(log_cpm), ] <- as.matrix(log_cpm)
        cat("✓ Completed dense matrix expansion\n")
      }
      
      cat("Dataset harmonization completed successfully!\n")
      return(list(
        counts = counts_harmonized,
        metadata = dataset$metadata,
        cpm_normalized = cpm_full,
        log_cpm = log_cpm_full
      ))
    } else {
      # If no genes have expression, create dummy normalized data
      cat("Warning: No genes with expression found\n")
      return(list(
        counts = counts_harmonized,
        metadata = dataset$metadata,
        cpm_normalized = counts_harmonized,
        log_cpm = counts_harmonized
      ))
    }
  })
  
} else {
  # Original code for when there are common genes
  mouse_datasets_harmonized <- map(mouse_datasets, function(dataset) {
    counts_filtered <- dataset$counts[common_genes, , drop = FALSE]
    
    # Convert sparse matrix to regular matrix for DGEList if needed
    if (inherits(counts_filtered, "dgCMatrix")) {
      counts_filtered <- as.matrix(counts_filtered)
    }
    
    # Perform basic normalization for visualization
    dge <- edgeR::DGEList(counts = counts_filtered)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    cpm_normalized <- edgeR::cpm(dge, log = FALSE)
    log_cpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
    
    return(list(
      counts = counts_filtered,
      metadata = dataset$metadata,
      cpm_normalized = cpm_normalized,
      log_cpm = log_cpm
    ))
  })
}

# ==============================================================================
# SAVE HARMONIZED DATASETS
# ==============================================================================

cat("Saving harmonized datasets...\n")

# Save individual harmonized datasets
iwalk(mouse_datasets_harmonized, function(dataset, name) {
  saveRDS(dataset, file.path(processed_dir, paste0(name, "_harmonized.rds")))
})

# Save male-only datasets
iwalk(mouse_datasets_male_only, function(dataset, name) {
  saveRDS(dataset, file.path(processed_dir, paste0(name, "_male_only.rds")))
})

# Save sex-stratified GSE289002
if (!is.null(gse289002_sex_stratified)) {
  saveRDS(gse289002_sex_stratified, 
          file.path(processed_dir, "GSE289002_sex_stratified.rds"))
}

# Save combined metadata
saveRDS(all_metadata, file.path(processed_dir, "mouse_combined_metadata.rds"))

# Save common genes
saveRDS(common_genes, file.path(processed_dir, "mouse_common_genes.rds"))

# ==============================================================================
# GENERATE QC PLOTS AND REPORTS
# ==============================================================================

cat("Generating QC plots and reports...\n")

# Create summary statistics
summary_stats <- map_dfr(mouse_datasets_harmonized, function(dataset) {
  # Handle potential issues with empty count matrices
  if (is.null(dataset$counts) || nrow(dataset$counts) == 0 || ncol(dataset$counts) == 0) {
    return(data.frame(
      n_samples = 0,
      n_genes = 0,
      mean_lib_size = 0,
      median_lib_size = 0,
      mean_genes_detected = 0,
      median_genes_detected = 0
    ))
  }
  
  # Handle sparse matrices for column sums
  if (inherits(dataset$counts, "dgCMatrix")) {
    lib_sizes <- as.numeric(Matrix::colSums(dataset$counts))
    genes_detected <- as.numeric(Matrix::colSums(dataset$counts > 0))
  } else {
    lib_sizes <- colSums(dataset$counts)
    genes_detected <- colSums(dataset$counts > 0)
  }
  
  data.frame(
    n_samples = ncol(dataset$counts),
    n_genes = nrow(dataset$counts),
    mean_lib_size = mean(lib_sizes),
    median_lib_size = median(lib_sizes),
    mean_genes_detected = mean(genes_detected),
    median_genes_detected = median(genes_detected)
  )
}, .id = "dataset")

# Sample distribution plots
cat("Creating sample distribution summary...\n")

# Ensure all_metadata is a proper data frame
if (!is.data.frame(all_metadata)) {
  cat("Converting all_metadata to data frame...\n")
  all_metadata <- as.data.frame(all_metadata)
}

# Check for required columns and create them if missing
required_cols <- c("dataset", "condition", "brain_region", "sex")
for (col in required_cols) {
  if (!col %in% colnames(all_metadata)) {
    cat("Warning: Missing column", col, "- adding with 'Unknown' values\n")
    all_metadata[[col]] <- "Unknown"
  }
}

# Debug before count operation
cat("Final all_metadata structure before count:\n")
cat("Class:", class(all_metadata), "\n")
cat("Is data frame:", is.data.frame(all_metadata), "\n")
cat("Columns:", paste(colnames(all_metadata), collapse = ", "), "\n")

# Create metadata summary with error handling
metadata_summary <- NULL
tryCatch({
  metadata_summary <- all_metadata %>%
    count(dataset, condition, brain_region, sex, .drop = FALSE) %>%
    filter(n > 0)
  
  cat("Successfully created metadata summary with", nrow(metadata_summary), "rows\n")
}, error = function(e) {
  cat("Error in count operation:", e$message, "\n")
  cat("Creating basic metadata summary...\n")
  
  # Create a basic summary manually using base R
  metadata_summary <<- all_metadata %>%
    group_by(dataset, condition, brain_region, sex) %>%
    summarise(n = n(), .groups = "drop") %>%
    filter(n > 0)
  
  cat("Created basic metadata summary with", nrow(metadata_summary), "rows\n")
})

# If metadata_summary is still NULL, create a very basic one
if (is.null(metadata_summary)) {
  cat("Creating minimal metadata summary...\n")
  metadata_summary <- data.frame(
    dataset = unique(all_metadata$dataset),
    condition = "Mixed",
    brain_region = "Mixed", 
    sex = "Mixed",
    n = as.numeric(table(all_metadata$dataset)),
    stringsAsFactors = FALSE
  )
}

p_sample_dist <- ggplot(metadata_summary, aes(x = dataset, y = n, fill = condition)) +
  geom_col(position = "stack") +
  facet_wrap(~ brain_region + sex, scales = "free") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Sample Distribution Across Mouse Datasets",
       x = "Dataset", y = "Number of Samples", fill = "Condition") +
  scale_fill_viridis_d()

ggsave(file.path(figures_dir, "mouse_sample_distribution.png"), 
       p_sample_dist, width = 12, height = 8, dpi = 300)

# Library size comparison
lib_sizes_df <- map_dfr(mouse_datasets_harmonized, function(dataset) {
  # Handle potential sparse matrix column sums
  if (inherits(dataset$counts, "dgCMatrix")) {
    lib_sizes <- Matrix::colSums(dataset$counts)
  } else {
    lib_sizes <- colSums(dataset$counts)
  }
  
  data.frame(
    sample_id = colnames(dataset$counts),
    lib_size = as.numeric(lib_sizes),
    dataset = dataset$metadata$dataset[1]
  )
}, .id = "dataset_id")

p_lib_sizes <- ggplot(lib_sizes_df, aes(x = dataset, y = lib_size, fill = dataset)) +
  geom_boxplot(alpha = 0.7) +
  scale_y_log10(labels = scales::comma) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Library Size Distribution Across Mouse Datasets",
       x = "Dataset", y = "Library Size (log10 scale)") +
  scale_fill_viridis_d() +
  guides(fill = "none")

ggsave(file.path(figures_dir, "mouse_library_sizes.png"), 
       p_lib_sizes, width = 10, height = 6, dpi = 300)

# Gene detection comparison
genes_detected_df <- map_dfr(mouse_datasets_harmonized, function(dataset) {
  # Handle potential sparse matrix operations
  if (inherits(dataset$counts, "dgCMatrix")) {
    genes_detected <- Matrix::colSums(dataset$counts > 0)
  } else {
    genes_detected <- colSums(dataset$counts > 0)
  }
  
  data.frame(
    sample_id = colnames(dataset$counts),
    genes_detected = as.numeric(genes_detected),
    dataset = dataset$metadata$dataset[1]
  )
}, .id = "dataset_id")

p_genes_detected <- ggplot(genes_detected_df, aes(x = dataset, y = genes_detected, fill = dataset)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Genes Detected per Sample Across Mouse Datasets",
       x = "Dataset", y = "Number of Genes Detected") +
  scale_fill_viridis_d() +
  guides(fill = "none")

ggsave(file.path(figures_dir, "mouse_genes_detected.png"), 
       p_genes_detected, width = 10, height = 6, dpi = 300)

# Create Excel report
wb <- createWorkbook()

# Summary statistics
addWorksheet(wb, "Summary_Statistics")
writeData(wb, "Summary_Statistics", summary_stats)

# Combined metadata
addWorksheet(wb, "Combined_Metadata")
writeData(wb, "Combined_Metadata", all_metadata)

# Sample distribution
addWorksheet(wb, "Sample_Distribution")
writeData(wb, "Sample_Distribution", metadata_summary)

# Common genes info
addWorksheet(wb, "Common_Genes_Info")
writeData(wb, "Common_Genes_Info", 
          data.frame(
            Total_Common_Genes = length(common_genes),
            Datasets_Compared = paste(names(mouse_datasets), collapse = ", "),
            Harmonization_Date = Sys.Date()
          ))

# Individual dataset info
for (dataset_name in names(mouse_datasets_harmonized)) {
  sheet_name <- paste0(dataset_name, "_Info")
  addWorksheet(wb, sheet_name)
  
  dataset_info <- mouse_datasets_harmonized[[dataset_name]]$metadata %>%
    count(condition, brain_region, sex, time_point, .drop = FALSE) %>%
    filter(n > 0)
  
  writeData(wb, sheet_name, dataset_info)
}

# Save Excel file
saveWorkbook(wb, file.path(results_dir, "mouse_datasets_harmonization_report.xlsx"), 
             overwrite = TRUE)

# ==============================================================================
# PRINT FINAL SUMMARY
# ==============================================================================

cat("\n", strrep("=", 60), "\n")
cat("MOUSE DATA HARMONIZATION SUMMARY\n")
cat(strrep("=", 60), "\n")

cat("Successfully processed datasets:\n")
for (dataset_name in names(mouse_datasets_harmonized)) {
  dataset <- mouse_datasets_harmonized[[dataset_name]]
  cat("", dataset_name, ":\n")
  cat("    Samples:", ncol(dataset$counts), "\n")
  cat("    Genes:", nrow(dataset$counts), "\n")
  cat("    Conditions:", paste(unique(dataset$metadata$condition), collapse = ", "), "\n")
  cat("    Brain regions:", paste(unique(dataset$metadata$brain_region), collapse = ", "), "\n")
  cat("    Sex distribution:", paste(table(dataset$metadata$sex), collapse = " / "), "\n\n")
}

cat("Common genes across all datasets:", length(common_genes), "\n")
cat("Male-only datasets created for cross-species comparison\n")

if (!is.null(gse289002_sex_stratified)) {
  cat("Sex-stratified GSE289002 dataset created:\n")
  cat("  Male samples:", ncol(gse289002_sex_stratified$male$counts), "\n")
  cat("  Female samples:", ncol(gse289002_sex_stratified$female$counts), "\n")
}

cat("\nAll harmonized data saved to:", processed_dir, "\n")
cat("QC plots saved to:", figures_dir, "\n")
cat("Harmonization report saved to:", results_dir, "\n")
cat(strrep("=", 60), "\n")
