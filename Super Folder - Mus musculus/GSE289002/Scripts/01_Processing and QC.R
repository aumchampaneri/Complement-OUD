# Enhanced Bulk RNA-seq QC and Preprocessing with Comprehensive Analysis
# =====================================================================

# Set reproducibility seed
set.seed(42)

# Load required libraries with proper error handling
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(RColorBrewer)
  library(pheatmap)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
  # Optional packages with error handling
  tryCatch(library(PCAtools), error = function(e) cat("PCAtools not available\n"))
  tryCatch(library(sva), error = function(e) cat("sva not available\n"))
  tryCatch(library(pwr), error = function(e) cat("pwr not available\n"))
  tryCatch(library(preprocessCore), error = function(e) cat("preprocessCore not available\n"))
  tryCatch(library(ggforce), error = function(e) cat("ggforce not available\n"))
})

# ----------------------
# DIRECTORY SETUP & ORGANIZATION
# ----------------------

# Base directory
base_dir <- "/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus/GSE289002"

# Create organized output directories
output_dirs <- list(
  base = file.path(base_dir, "Outputs/Enhanced_QC"),
  plots = file.path(base_dir, "Outputs/Enhanced_QC/Plots"),
  data = file.path(base_dir, "Outputs/Enhanced_QC/Processed_Data"),
  reports = file.path(base_dir, "Outputs/Enhanced_QC/Reports"),
  stats = file.path(base_dir, "Outputs/Enhanced_QC/Statistics"),
  outliers = file.path(base_dir, "Outputs/Enhanced_QC/Outlier_Analysis"),
  normalization = file.path(base_dir, "Outputs/Enhanced_QC/Normalization_Comparison"),
  batch = file.path(base_dir, "Outputs/Enhanced_QC/Batch_Analysis")
)

# Create all directories
lapply(output_dirs, function(x) dir.create(x, showWarnings = FALSE, recursive = TRUE))

# File paths
count_file <- file.path(base_dir, "Raw_Data/GSE289002_mouse_raw_counts.csv")
soft_file <- file.path(base_dir, "Raw_Data/GSE289002_family.soft")

cat("=== ENHANCED BULK RNA-SEQ PREPROCESSING PIPELINE ===\n")
cat("Analysis Date:", as.character(Sys.time()), "\n")
cat("Output Directory:", output_dirs$base, "\n\n")

# ----------------------
# 1. METADATA EXTRACTION FROM SOFT FILE
# ----------------------

parse_soft_file <- function(soft_file_path) {
  cat("Parsing SOFT file for metadata extraction...\n")
  
  samples <- list()
  current_sample <- NULL
  sample_count <- 0
  
  con <- file(soft_file_path, "r")
  
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    line <- trimws(line)
    
    # Start of a new sample
    if (grepl("^\\^SAMPLE =", line)) {
      if (!is.null(current_sample)) {
        samples[[length(samples) + 1]] <- current_sample
      }
      sample_id <- trimws(sub("^\\^SAMPLE =", "", line))
      current_sample <- list(geo_accession = sample_id)
      sample_count <- sample_count + 1
    }
    
    # Sample attributes
    else if (grepl("^!Sample_", line) && !is.null(current_sample)) {
      parts <- strsplit(line, "=", fixed = TRUE)[[1]]
      if (length(parts) >= 2) {
        key <- trimws(sub("^!Sample_", "", parts[1]))
        value <- trimws(paste(parts[-1], collapse = "="))
        
        # For characteristics, parse the field: value format
        if (grepl("^characteristics_ch1", key)) {
          char_match <- regexpr("(.+?):\\s*(.+)", value, perl = TRUE)
          if (char_match > 0) {
            char_parts <- regmatches(value, char_match, invert = FALSE)
            if (length(char_parts) > 0) {
              char_split <- strsplit(char_parts, ":\\s*", perl = TRUE)[[1]]
              if (length(char_split) >= 2) {
                char_name <- trimws(char_split[1])
                char_value <- trimws(char_split[2])
                current_sample[[char_name]] <- char_value
              }
            }
          }
        } else {
          current_sample[[key]] <- value
        }
      }
    }
  }
  
  # Add the last sample
  if (!is.null(current_sample)) {
    samples[[length(samples) + 1]] <- current_sample
  }
  
  close(con)
  
  cat("Total samples found:", length(samples), "\n")
  
  # Filter for mouse samples only
  mouse_samples <- samples[sapply(samples, function(s) {
    organism <- s[["organism_ch1"]]
    !is.null(organism) && organism == "Mus musculus"
  })]
  
  cat("Mouse samples found:", length(mouse_samples), "\n")
  return(mouse_samples)
}

create_metadata_df <- function(mouse_samples) {
  cat("Creating metadata dataframe...\n")
  
  metadata_list <- lapply(mouse_samples, function(sample) {
    data.frame(
      geo_accession = sample[["geo_accession"]] %||% NA,
      title = sample[["title"]] %||% NA,
      treatment = sample[["treatment"]] %||% NA,
      sex = sample[["Sex"]] %||% sample[["sex"]] %||% NA,
      region = sample[["region"]] %||% NA,
      strain = sample[["strain"]] %||% NA,
      batch = sample[["batch"]] %||% NA,
      stringsAsFactors = FALSE
    )
  })
  
  metadata_df <- do.call(rbind, metadata_list)
  return(metadata_df)
}

# Helper function for null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

# Extract metadata from SOFT file
if (file.exists(soft_file)) {
  mouse_samples <- parse_soft_file(soft_file)
  metadata <- create_metadata_df(mouse_samples)
  
  # Save metadata for future use
  write.csv(metadata, file.path(output_dirs$data, "original_metadata.csv"), row.names = FALSE)
  cat("Metadata saved to:", file.path(output_dirs$data, "original_metadata.csv"), "\n")
} else {
  # Fallback to existing metadata file if SOFT file not found
  metadata_file <- file.path(dirname(count_file), "mouse_metadata.csv")
  if (file.exists(metadata_file)) {
    cat("SOFT file not found, loading existing metadata...\n")
    metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
  } else {
    stop("Neither SOFT file nor metadata file found!")
  }
}

# ----------------------
# 2. DATA LOADING AND INITIAL ASSESSMENT
# ----------------------

cat("Loading count data...\n")
counts_original <- read.csv(count_file, row.names = 1, check.names = FALSE)

# Print metadata summary
cat("\n=== INITIAL DATA ASSESSMENT ===\n")
cat("Original count matrix dimensions:", dim(counts_original)[1], "genes x", dim(counts_original)[2], "samples\n")
cat("Metadata dimensions:", dim(metadata)[1], "samples x", dim(metadata)[2], "variables\n")
cat("\nTreatment values found:", paste(unique(metadata$treatment), collapse = ", "), "\n")
cat("Sex values found:", paste(unique(metadata$sex), collapse = ", "), "\n")
cat("Region values found:", paste(unique(metadata$region), collapse = ", "), "\n")

# ----------------------
# 3. COMPREHENSIVE OUTLIER DETECTION
# ----------------------

detect_outliers_comprehensive <- function(counts, metadata, output_dir) {
  cat("\n=== COMPREHENSIVE OUTLIER DETECTION ===\n")
  
  outliers <- list()
  
  # 1. Library size outliers (3 MAD rule)
  lib_sizes <- colSums(counts)
  median_lib <- median(lib_sizes)
  mad_lib <- mad(lib_sizes)
  lib_outliers <- names(lib_sizes)[abs(lib_sizes - median_lib) > 3 * mad_lib]
  
  # 2. PCA-based outliers (samples far from cluster centroid)
  pca_temp <- prcomp(t(log2(counts + 1)))
  pc_dist <- sqrt(rowSums(pca_temp$x[,1:5]^2))
  pca_outliers <- names(pc_dist)[pc_dist > quantile(pc_dist, 0.95)]
  
  # 3. Manual outliers (known problematic samples)
  manual_outliers <- c("PFC_10")
  
  # 4. Low gene detection outliers
  genes_detected <- colSums(counts > 0)
  low_gene_outliers <- names(genes_detected)[genes_detected < quantile(genes_detected, 0.05)]
  
  outliers <- list(
    library_size = lib_outliers,
    pca_based = pca_outliers,
    manual = manual_outliers,
    low_gene_detection = low_gene_outliers
  )
  
  # Comprehensive outlier report
  all_outliers <- unique(c(lib_outliers, pca_outliers, manual_outliers, low_gene_outliers))
  
  # Create outlier summary
  if (length(all_outliers) > 0) {
    outlier_summary <- data.frame(
      sample = all_outliers,
      library_size_outlier = all_outliers %in% lib_outliers,
      pca_outlier = all_outliers %in% pca_outliers,
      manual_outlier = all_outliers %in% manual_outliers,
      low_gene_outlier = all_outliers %in% low_gene_outliers,
      library_size = lib_sizes[all_outliers],
      genes_detected = genes_detected[all_outliers]
    )
  } else {
    outlier_summary <- data.frame()
  }
  
  # Save outlier analysis
  write.csv(outlier_summary, file.path(output_dir, "outlier_detection_summary.csv"), row.names = FALSE)
  
  # Plot outlier detection
  png(file.path(output_dirs$plots, "outlier_detection_plots.png"), width = 1600, height = 1200)
  par(mfrow = c(2, 3))
  
  # Library size plot
  plot(lib_sizes, main = "Library Sizes", ylab = "Total Counts", pch = 16)
  abline(h = median_lib + 3*mad_lib, col = "red", lty = 2)
  abline(h = median_lib - 3*mad_lib, col = "red", lty = 2)
  if (length(lib_outliers) > 0) {
    points(which(names(lib_sizes) %in% lib_outliers), lib_sizes[lib_outliers], col = "red", pch = 16, cex = 1.5)
  }
  
  # PCA distance plot
  plot(pc_dist, main = "PCA Distance from Center", ylab = "Euclidean Distance", pch = 16)
  abline(h = quantile(pc_dist, 0.95), col = "red", lty = 2)
  if (length(pca_outliers) > 0) {
    points(which(names(pc_dist) %in% pca_outliers), pc_dist[pca_outliers], col = "red", pch = 16, cex = 1.5)
  }
  
  # Gene detection plot
  plot(genes_detected, main = "Genes Detected per Sample", ylab = "Number of Genes", pch = 16)
  abline(h = quantile(genes_detected, 0.05), col = "red", lty = 2)
  if (length(low_gene_outliers) > 0) {
    points(which(names(genes_detected) %in% low_gene_outliers), genes_detected[low_gene_outliers], col = "red", pch = 16, cex = 1.5)
  }
  
  # PCA plot highlighting outliers
  plot(pca_temp$x[,1], pca_temp$x[,2], main = "PCA - Outliers Highlighted", 
       xlab = "PC1", ylab = "PC2", pch = 16)
  if (length(all_outliers) > 0) {
    points(pca_temp$x[all_outliers,1], pca_temp$x[all_outliers,2], col = "red", pch = 16, cex = 1.5)
    text(pca_temp$x[all_outliers,1], pca_temp$x[all_outliers,2], labels = all_outliers, pos = 3, col = "red")
  }
  
  dev.off()
  
  cat("Detected outliers:\n")
  cat("- Library size outliers:", length(lib_outliers), "\n")
  cat("- PCA-based outliers:", length(pca_outliers), "\n")
  cat("- Manual outliers:", length(manual_outliers), "\n")
  cat("- Low gene detection outliers:", length(low_gene_outliers), "\n")
  cat("- Total unique outliers:", length(all_outliers), "\n")
  
  return(list(outliers = outliers, all_outliers = all_outliers, summary = outlier_summary))
}

# Detect outliers
outlier_results <- detect_outliers_comprehensive(counts_original, metadata, output_dirs$outliers)
outliers_to_remove <- outlier_results$all_outliers

# ----------------------
# 4. SAMPLE REMOVAL AND DATA MATCHING
# ----------------------

cat("\n=== SAMPLE REMOVAL AND DATA MATCHING ===\n")

# Remove outlier samples
if (length(outliers_to_remove) > 0) {
  cat("Removing", length(outliers_to_remove), "outlier samples:", paste(outliers_to_remove, collapse = ", "), "\n")
  
  # Remove from counts
  counts <- counts_original[, !colnames(counts_original) %in% outliers_to_remove]
  
  # Remove from metadata
  metadata <- metadata[!metadata$title %in% outliers_to_remove & 
                      !metadata$geo_accession %in% outliers_to_remove, ]
} else {
  counts <- counts_original
  cat("No outliers to remove.\n")
}

# Match samples between counts and metadata
cat("Matching samples between counts and metadata...\n")
common_samples <- intersect(colnames(counts), metadata$title)

if (length(common_samples) == 0) {
  common_samples <- intersect(colnames(counts), metadata$geo_accession)
  if (length(common_samples) > 0) {
    cat("Using geo_accession for sample matching. Found", length(common_samples), "matches.\n")
    metadata <- metadata[metadata$geo_accession %in% colnames(counts), ]
    counts <- counts[, metadata$geo_accession]
    metadata$title <- metadata$geo_accession
  } else {
    stop("No matching samples between counts and metadata!")
  }
} else {
  metadata <- metadata[metadata$title %in% colnames(counts), ]
  counts <- counts[, metadata$title]
}

# Verify sample matching
all_match <- all(colnames(counts) == metadata$title)
cat("Sample matching successful:", all_match, "\n")
cat("Final dataset:", nrow(counts), "genes x", ncol(counts), "samples\n")

# ----------------------
# 5. ENHANCED METADATA PROCESSING
# ----------------------

cat("\n=== ENHANCED METADATA PROCESSING ===\n")

# Use geo_accession as batch identifier (each sample is its own batch)
metadata$batch <- metadata$geo_accession

# Clean and format metadata
metadata$treatment <- factor(metadata$treatment)
metadata$sex <- factor(metadata$sex)
metadata$region <- factor(metadata$region)

# Handle missing values
missing_treatment <- sum(is.na(metadata$treatment))
if (missing_treatment > 0) {
  cat("Warning: Found", missing_treatment, "samples with missing treatment information\n")
}

# Advanced QC metrics
calculate_advanced_qc <- function(counts, metadata) {
  cat("Calculating advanced QC metrics...\n")
  
  # Mitochondrial gene percentage
  mito_genes <- grep("^mt-|^Mt-", rownames(counts), value = TRUE)
  if(length(mito_genes) > 0) {
    mito_percent <- colSums(counts[mito_genes, ]) / colSums(counts) * 100
  } else {
    mito_percent <- rep(0, ncol(counts))
  }
  
  # Ribosomal gene percentage
  ribo_genes <- grep("^Rpl|^Rps", rownames(counts), value = TRUE)
  if(length(ribo_genes) > 0) {
    ribo_percent <- colSums(counts[ribo_genes, ]) / colSums(counts) * 100
  } else {
    ribo_percent <- rep(0, ncol(counts))
  }
  
  # Gene detection rate
  genes_detected <- colSums(counts > 0)
  
  # Library size
  lib_sizes <- colSums(counts)
  
  # Add to metadata
  metadata$lib_size <- lib_sizes
  metadata$mito_percent <- mito_percent
  metadata$ribo_percent <- ribo_percent
  metadata$genes_detected <- genes_detected
  
  return(metadata)
}

metadata <- calculate_advanced_qc(counts, metadata)

# Display experimental design
cat("Experimental design summary:\n")
cat("Treatments:", paste(levels(metadata$treatment), collapse = ", "), "\n")
cat("Brain regions:", paste(levels(metadata$region), collapse = ", "), "\n")
cat("Sex distribution:", paste(table(metadata$sex), collapse = " / "), "\n")

# ----------------------
# 6. EXPERIMENTAL DESIGN ASSESSMENT
# ----------------------

assess_experimental_design <- function(metadata, output_dir) {
  cat("\n=== EXPERIMENTAL DESIGN ASSESSMENT ===\n")
  
  design_table <- table(metadata$treatment, metadata$sex)
  region_table <- table(metadata$treatment, metadata$region)
  min_group_size <- min(design_table[design_table > 0])
  
  # Power analysis (only if pwr package is available)
  power_analysis <- NULL
  if ("pwr" %in% loadedNamespaces()) {
    power_small <- pwr.t.test(n = min_group_size, d = 0.2, sig.level = 0.05, type = "two.sample")
    power_medium <- pwr.t.test(n = min_group_size, d = 0.5, sig.level = 0.05, type = "two.sample")
    power_large <- pwr.t.test(n = min_group_size, d = 0.8, sig.level = 0.05, type = "two.sample")
    
    power_analysis <- data.frame(
      effect_size = c("small (0.2)", "medium (0.5)", "large (0.8)"),
      power = c(power_small$power, power_medium$power, power_large$power)
    )
  }
  
  # Test for confounding - updated for meaningful comparisons
  confounding_tests <- list()
  
  # Test treatment vs sex
  sex_levels <- length(unique(metadata$sex[!is.na(metadata$sex)]))
  treatment_levels <- length(unique(metadata$treatment[!is.na(metadata$treatment)]))
  
  if(sex_levels > 1 && treatment_levels > 1) {
    tryCatch({
      sex_treatment_test <- chisq.test(metadata$treatment, metadata$sex)
      confounding_tests$sex_treatment <- sex_treatment_test$p.value
    }, error = function(e) {
      cat("Warning: Could not perform sex-treatment chi-square test:", e$message, "\n")
      confounding_tests$sex_treatment <- NA
    })
  } else {
    confounding_tests$sex_treatment <- NA
  }
  
  # Test treatment vs region
  region_levels <- length(unique(metadata$region[!is.na(metadata$region)]))
  
  if(region_levels > 1 && treatment_levels > 1) {
    tryCatch({
      region_treatment_test <- chisq.test(metadata$treatment, metadata$region)
      confounding_tests$region_treatment <- region_treatment_test$p.value
    }, error = function(e) {
      cat("Warning: Could not perform region-treatment chi-square test:", e$message, "\n")
      confounding_tests$region_treatment <- NA
    })
  } else {
    confounding_tests$region_treatment <- NA
  }
  
  # Create design assessment report
  design_report <- list(
    design_table = design_table,
    region_table = region_table,
    min_group_size = min_group_size,
    power_analysis = power_analysis,
    confounding_tests = confounding_tests
  )
  
  # Save design assessment
  sink(file.path(output_dir, "experimental_design_assessment.txt"))
  cat("=== EXPERIMENTAL DESIGN ASSESSMENT ===\n\n")
  cat("Design Table (Treatment x Sex):\n")
  print(design_table)
  cat("\nDesign Table (Treatment x Region):\n")
  print(region_table)
  cat("\nMinimum group size:", min_group_size, "\n\n")
  
  if (!is.null(power_analysis)) {
    cat("Power Analysis:\n")
    print(power_analysis)
    cat("\n")
  } else {
    cat("Power Analysis: pwr package not available\n\n")
  }
  
  cat("Confounding Tests (p-values):\n")
  for(test in names(confounding_tests)) {
    if(!is.na(confounding_tests[[test]])) {
      cat(test, ":", round(confounding_tests[[test]], 4), "\n")
    } else {
      cat(test, ": Not applicable\n")
    }
  }
  cat("\nBatch Information:\n")
  cat("Using geo_accession as batch identifier\n")
  cat("Total samples (unique batches):", length(unique(metadata$batch)), "\n")
  cat("Note: Each sample represents its own batch\n")
  sink()
  
  return(design_report)
}

design_assessment <- assess_experimental_design(metadata, output_dirs$reports)

# ----------------------
# 7. COMPREHENSIVE GENE-LEVEL QC
# ----------------------

cat("\n=== COMPREHENSIVE GENE-LEVEL QC ===\n")

# Calculate comprehensive gene statistics
gene_stats <- data.frame(
  gene_id = rownames(counts),
  mean_count = rowMeans(counts),
  median_count = apply(counts, 1, median),
  min_count = apply(counts, 1, min),
  max_count = apply(counts, 1, max),
  var_count = apply(counts, 1, var),
  cv_count = apply(counts, 1, function(x) sd(x)/mean(x)),
  zero_samples = rowSums(counts == 0),
  pct_zeros = rowSums(counts == 0) / ncol(counts) * 100,
  stringsAsFactors = FALSE
)

# Add gene biotype information if available
gene_stats$log_mean <- log10(gene_stats$mean_count + 1)

# Save comprehensive gene stats
write.csv(gene_stats, file.path(output_dirs$stats, "comprehensive_gene_statistics.csv"), row.names = FALSE)

# Enhanced gene expression visualization
png(file.path(output_dirs$plots, "gene_expression_comprehensive.png"), width = 1600, height = 1200)
par(mfrow = c(2, 3))

# Expression distribution
hist(gene_stats$log_mean, breaks = 100, main = "Gene Expression Distribution",
     xlab = "log10(mean count + 1)", col = "lightblue", border = "white")
abline(v = log10(10), col = "red", lwd = 2, lty = 2)

# Zero inflation
hist(gene_stats$pct_zeros, breaks = 50, main = "Zero Inflation Distribution",
     xlab = "% Samples with Zero Counts", col = "lightcoral", border = "white")

# Mean-variance relationship
plot(gene_stats$log_mean, log10(gene_stats$var_count + 1), 
     main = "Mean-Variance Relationship", pch = 16, col = rgb(0,0,0,0.3),
     xlab = "log10(mean)", ylab = "log10(variance)")

# Coefficient of variation
cv_finite <- gene_stats$cv_count[is.finite(gene_stats$cv_count)]
if(length(cv_finite) > 0) {
  plot(gene_stats$log_mean[is.finite(gene_stats$cv_count)], cv_finite, 
       main = "Coefficient of Variation", pch = 16, col = rgb(0,0,0,0.3),
       xlab = "log10(mean)", ylab = "CV", ylim = c(0, quantile(cv_finite, 0.95, na.rm = TRUE)))
} else {
  plot.new()
  text(0.5, 0.5, "CV calculation failed", cex = 1.5)
}

# Detection frequency
plot(gene_stats$mean_count, ncol(counts) - gene_stats$zero_samples,
     main = "Gene Detection Frequency", pch = 16, col = rgb(0,0,0,0.3),
     xlab = "Mean Count", ylab = "Samples Detected")

dev.off()

# ----------------------
# 8. ENHANCED FILTERING
# ----------------------

cat("\n=== ENHANCED GENE FILTERING ===\n")

# Create DGEList object
dge <- DGEList(counts = counts, samples = metadata)

# Enhanced filtering function
filter_genes_enhanced <- function(dge, metadata, min_cpm = 1, min_samples = 3) {
  cat("Applying enhanced gene filtering...\n")
  
  cpm_vals <- cpm(dge)
  group_sizes <- table(metadata$treatment)
  min_group_size <- min(group_sizes)
  
  # Require expression in at least min_samples per group
  keep_expression <- rowSums(cpm_vals >= min_cpm) >= min(min_samples, min_group_size)
  
  # Remove genes with extremely low variance (likely noise)
  gene_vars <- apply(cpm_vals, 1, var)
  keep_variance <- gene_vars > quantile(gene_vars, 0.1, na.rm = TRUE)
  
  # Standard edgeR filtering
  keep_edger <- filterByExpr(dge, group = metadata$treatment, min.count = 10, min.total.count = 15)
  
  # Combine filters
  keep_final <- keep_expression & keep_variance & keep_edger
  
  cat("Genes passing expression filter:", sum(keep_expression), "\n")
  cat("Genes passing variance filter:", sum(keep_variance), "\n")
  cat("Genes passing edgeR filter:", sum(keep_edger), "\n")
  cat("Genes passing all filters:", sum(keep_final), "\n")
  
  return(keep_final)
}

# Apply enhanced filtering
keep_genes <- filter_genes_enhanced(dge, metadata)
dge_filtered <- dge[keep_genes, , keep.lib.sizes = FALSE]

cat("Original genes:", nrow(dge), "\n")
cat("Filtered genes:", nrow(dge_filtered), "\n")
cat("Filtering removed:", nrow(dge) - nrow(dge_filtered), "genes\n")

# ----------------------
# 9. NORMALIZATION COMPARISON
# ----------------------

cat("\n=== NORMALIZATION METHODS COMPARISON ===\n")

compare_normalization_methods <- function(dge, output_dir) {
  cat("Comparing normalization methods...\n")
  
  # TMM normalization
  dge_tmm <- calcNormFactors(dge, method = "TMM")
  
  # RLE normalization
  dge_rle <- calcNormFactors(dge, method = "RLE")
  
  # Upper quartile normalization
  dge_uq <- calcNormFactors(dge, method = "upperquartile")
  
  # Get log-CPM values
  logcpm_tmm <- cpm(dge_tmm, log = TRUE)
  logcpm_rle <- cpm(dge_rle, log = TRUE)
  logcpm_uq <- cpm(dge_uq, log = TRUE)
  
  # Plot comparisons
  png(file.path(output_dir, "normalization_comparison.png"), width = 1600, height = 800)
  par(mfrow = c(2, 3))
  
  # Boxplots
  boxplot(logcpm_tmm, main = "TMM Normalization", las = 2, cex.axis = 0.6)
  boxplot(logcpm_rle, main = "RLE Normalization", las = 2, cex.axis = 0.6)
  boxplot(logcpm_uq, main = "Upper Quartile Normalization", las = 2, cex.axis = 0.6)
  
  # Density plots
  plotDensities(logcpm_tmm, main = "TMM - Density", legend = FALSE)
  plotDensities(logcpm_rle, main = "RLE - Density", legend = FALSE)
  plotDensities(logcpm_uq, main = "Upper Quartile - Density", legend = FALSE)
  
  dev.off()
  
  # Return TMM as default (best for most cases)
  return(list(
    tmm = dge_tmm,
    rle = dge_rle,
    uq = dge_uq,
    logcpm_tmm = logcpm_tmm
  ))
}

normalization_results <- compare_normalization_methods(dge_filtered, output_dirs$normalization)
dge_normalized <- normalization_results$tmm
logcpm_normalized <- normalization_results$logcpm_tmm

# ----------------------
# 10. SAMPLE-LEVEL VARIABILITY ANALYSIS
# ----------------------

cat("\n=== SAMPLE-LEVEL VARIABILITY ANALYSIS ===\n")

analyze_sample_effects <- function(dge, metadata, logcpm, output_dir) {
  cat("Analyzing sample-level variability...\n")
  cat("Note: Since each sample is its own batch (geo_accession), analyzing individual sample effects\n")
  
  # PCA for sample-level assessment
  pca_samples <- prcomp(t(logcpm))
  variance_explained_samples <- (pca_samples$sdev^2) / sum(pca_samples$sdev^2) * 100
  
  # Plot sample-level effects
  png(file.path(output_dir, "sample_level_analysis.png"), width = 1600, height = 1200)
  par(mfrow = c(2, 2))
  
  # PCA by treatment and region
  plot(pca_samples$x[,1], pca_samples$x[,2], 
       col = rainbow(length(levels(metadata$treatment)))[metadata$treatment],
       pch = as.numeric(metadata$region) + 15,
       cex = 1.5,
       xlab = paste0("PC1 (", round(variance_explained_samples[1], 2), "%)"),
       ylab = paste0("PC2 (", round(variance_explained_samples[2], 2), "%)"),
       main = "PCA - Treatment & Region")
  legend("topright", 
         legend = c(paste("Treatment:", levels(metadata$treatment)),
                   paste("Region:", levels(metadata$region))),
         col = c(rainbow(length(levels(metadata$treatment))), rep(1, length(levels(metadata$region)))),
         pch = c(rep(16, length(levels(metadata$treatment))), 
                 16:(15+length(levels(metadata$region)))))
  
  # PCA by sex and treatment
  plot(pca_samples$x[,1], pca_samples$x[,2], 
       col = rainbow(length(levels(metadata$treatment)))[metadata$treatment],
       pch = as.numeric(metadata$sex) + 15,
       cex = 1.5,
       xlab = paste0("PC1 (", round(variance_explained_samples[1], 2), "%)"),
       ylab = paste0("PC2 (", round(variance_explained_samples[2], 2), "%)"),
       main = "PCA - Treatment & Sex")
  legend("topright", 
         legend = c(paste("Treatment:", levels(metadata$treatment)),
                   paste("Sex:", levels(metadata$sex))),
         col = c(rainbow(length(levels(metadata$treatment))), rep(1, length(levels(metadata$sex)))),
         pch = c(rep(16, length(levels(metadata$treatment))), 
                 16:(15+length(levels(metadata$sex)))))
  
  # Sample variability by treatment
  sample_vars <- apply(logcpm, 2, var)
  boxplot(sample_vars ~ metadata$treatment, 
          main = "Sample Variability by Treatment",
          xlab = "Treatment", ylab = "Sample Variance",
          col = rainbow(length(levels(metadata$treatment))))
  
  # Sample variability by region
  boxplot(sample_vars ~ metadata$region, 
          main = "Sample Variability by Region",
          xlab = "Region", ylab = "Sample Variance",
          col = rainbow(length(levels(metadata$region))))
  
  dev.off()
  
  # Calculate treatment and region effects on PC1 and PC2
  pc1_treatment_var <- summary(lm(pca_samples$x[,1] ~ metadata$treatment))$r.squared
  pc2_treatment_var <- summary(lm(pca_samples$x[,2] ~ metadata$treatment))$r.squared
  pc1_region_var <- summary(lm(pca_samples$x[,1] ~ metadata$region))$r.squared
  pc2_region_var <- summary(lm(pca_samples$x[,2] ~ metadata$region))$r.squared
  
  cat("Treatment explains", round(pc1_treatment_var * 100, 2), "% of PC1 variance\n")
  cat("Treatment explains", round(pc2_treatment_var * 100, 2), "% of PC2 variance\n")
  cat("Region explains", round(pc1_region_var * 100, 2), "% of PC1 variance\n")
  cat("Region explains", round(pc2_region_var * 100, 2), "% of PC2 variance\n")
  
  return(list(
    pca = pca_samples,
    variance_explained = variance_explained_samples,
    variance_components = list(
      pc1_treatment = pc1_treatment_var,
      pc2_treatment = pc2_treatment_var,
      pc1_region = pc1_region_var,
      pc2_region = pc2_region_var
    )
  ))
}

batch_analysis <- analyze_sample_effects(dge_normalized, metadata, logcpm_normalized, output_dirs$batch)

# ----------------------
# 11. COMPREHENSIVE SAMPLE QC
# ----------------------

cat("\n=== COMPREHENSIVE SAMPLE-LEVEL QC ===\n")

# Enhanced sample-level visualizations - FIXED VERSION
create_comprehensive_sample_plots <- function(dge, metadata, logcpm, output_dir) {
  cat("Creating comprehensive sample QC plots...\n")
  
  # FIXED: Use dge$counts instead of counts(dge)
  lib_sizes <- colSums(dge$counts)
  
  # Main QC plot grid
  png(file.path(output_dir, "comprehensive_sample_qc.png"), width = 2000, height = 1600)
  layout(matrix(c(1,2,3,4,5,6,7,8), nrow = 2, byrow = TRUE))
  
  # 1. Library sizes with treatment coloring
  barplot(lib_sizes/1e6, names.arg = metadata$title, las = 2, cex.names = 0.4,
          main = "Library Sizes by Sample", ylab = "Millions of Reads",
          col = rainbow(length(levels(metadata$treatment)))[metadata$treatment])
  abline(h = mean(lib_sizes/1e6), col = "red", lwd = 2, lty = 2)
  
  # 2. Genes detected
  barplot(metadata$genes_detected, names.arg = metadata$title, las = 2, cex.names = 0.4,
          main = "Genes Detected per Sample", ylab = "Number of Genes",
          col = rainbow(length(levels(metadata$treatment)))[metadata$treatment])
  
  # 3. Mitochondrial percentage
  if(max(metadata$mito_percent) > 0) {
    barplot(metadata$mito_percent, names.arg = metadata$title, las = 2, cex.names = 0.4,
            main = "Mitochondrial Gene %", ylab = "% Mitochondrial",
            col = rainbow(length(levels(metadata$treatment)))[metadata$treatment])
  } else {
    plot.new()
    text(0.5, 0.5, "No mitochondrial genes detected", cex = 1.5)
  }
  
  # 4. Log-CPM boxplot
  boxplot(logcpm, las = 2, col = rainbow(length(levels(metadata$treatment)))[metadata$treatment],
          main = "Log2-CPM Distribution", cex.axis = 0.4)
  
  # 5. Sample correlation heatmap
  sample_cor <- cor(logcpm)
  image(sample_cor, main = "Sample Correlation Heatmap", axes = FALSE)
  
  # 6. PCA plot
  pca_res <- prcomp(t(logcpm))
  variance_explained <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100
  plot(pca_res$x[,1], pca_res$x[,2],
       col = rainbow(length(levels(metadata$treatment)))[metadata$treatment],
       pch = as.numeric(metadata$sex) + 15, cex = 1.5,
       xlab = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
       ylab = paste0("PC2 (", round(variance_explained[2], 2), "%)"),
       main = "PCA - Treatment & Sex")
  
  # 7. MDS plot
  plotMDS(dge, col = rainbow(length(levels(metadata$treatment)))[metadata$treatment],
          main = "MDS Plot")
  
  # 8. Density plot
  plotDensities(logcpm, col = rainbow(length(levels(metadata$treatment)))[metadata$treatment],
                legend = FALSE, main = "Expression Density")
  
  dev.off()
  
  return(list(pca = pca_res, variance_explained = variance_explained))
}

sample_qc_results <- create_comprehensive_sample_plots(dge_normalized, metadata, logcpm_normalized, output_dirs$plots)

# ----------------------
# 12. ENHANCED VISUALIZATIONS
# ----------------------

cat("\n=== CREATING ENHANCED VISUALIZATIONS ===\n")

# Create publication-quality PCA plot (only if ggforce is available)
create_publication_pca <- function(pca_res, metadata, variance_explained, output_file) {
  
  pca_data <- data.frame(
    PC1 = pca_res$x[,1],
    PC2 = pca_res$x[,2],
    PC3 = pca_res$x[,3],
    Treatment = metadata$treatment,
    Sex = metadata$sex,
    Region = metadata$region,
    Sample = metadata$title
  )
  
  if ("ggforce" %in% loadedNamespaces()) {
    p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Treatment, shape = Sex)) +
      geom_point(size = 4, alpha = 0.8) +
      geom_mark_ellipse(aes(group = Treatment), alpha = 0.1, expand = 0.01) +
      labs(x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
           y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
           title = "Principal Component Analysis",
           subtitle = paste("n =", nrow(pca_data), "samples")) +
      theme_bw(base_size = 14) +
      theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        panel.grid.minor = element_blank()
      ) +
      scale_color_brewer(type = "qual", palette = "Set1") +
      guides(color = guide_legend(override.aes = list(size = 3)),
             shape = guide_legend(override.aes = list(size = 3)))
    
    ggsave(output_file, p, width = 12, height = 8, dpi = 300)
  } else {
    # Create basic plot if ggforce not available
    png(output_file, width = 1200, height = 800)
    plot(pca_data$PC1, pca_data$PC2,
         col = rainbow(length(levels(pca_data$Treatment)))[pca_data$Treatment],
         pch = as.numeric(pca_data$Sex) + 15, cex = 1.5,
         xlab = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
         ylab = paste0("PC2 (", round(variance_explained[2], 2), "%)"),
         main = "Principal Component Analysis")
    legend("topright", 
           legend = c(paste("Treatment:", levels(pca_data$Treatment)),
                     paste("Sex:", levels(pca_data$Sex))),
           col = c(rainbow(length(levels(pca_data$Treatment))), rep(1, length(levels(pca_data$Sex)))),
           pch = c(rep(16, length(levels(pca_data$Treatment))), 
                   16:(15+length(levels(pca_data$Sex)))))
    dev.off()
  }
}

create_publication_pca(sample_qc_results$pca, metadata, 
                      sample_qc_results$variance_explained,
                      file.path(output_dirs$plots, "publication_quality_pca.png"))

# Create comprehensive correlation heatmap
create_enhanced_correlation_heatmap <- function(logcpm, metadata, output_file) {
  
  sample_cor <- cor(logcpm)
  
  annotation_col <- data.frame(
    Treatment = metadata$treatment,
    Sex = metadata$sex,
    Region = metadata$region,
    row.names = colnames(logcpm)
  )
  
  # FIXED: Color schemes with proper length handling
  treatment_colors <- RColorBrewer::brewer.pal(max(3, min(9, length(levels(metadata$treatment)))), "Set1")[1:length(levels(metadata$treatment))]
  sex_colors <- c("skyblue", "pink")[1:length(levels(metadata$sex))]
  region_colors <- RColorBrewer::brewer.pal(max(3, length(levels(metadata$region))), "Set2")[1:length(levels(metadata$region))]
  
  ann_colors <- list(
    Treatment = treatment_colors,
    Sex = sex_colors,
    Region = region_colors
  )
  names(ann_colors$Treatment) <- levels(metadata$treatment)
  names(ann_colors$Sex) <- levels(metadata$sex)
  names(ann_colors$Region) <- levels(metadata$region)
  
  png(output_file, width = 1400, height = 1200, res = 300)
  pheatmap(sample_cor, 
           annotation_col = annotation_col, 
           annotation_colors = ann_colors,
           show_colnames = FALSE,
           show_rownames = FALSE,
           main = "Sample Correlation Heatmap",
           fontsize = 10,
           cellwidth = 8,
           cellheight = 8)
  dev.off()
}

create_enhanced_correlation_heatmap(logcpm_normalized, metadata, 
                                  file.path(output_dirs$plots, "enhanced_correlation_heatmap.png"))

# ----------------------
# 13. DATA EXPORT AND ORGANIZATION
# ----------------------

cat("\n=== SAVING PROCESSED DATA ===\n")

# Save all processed objects
saveRDS(dge_normalized, file.path(output_dirs$data, "dge_normalized_final.rds"))
saveRDS(metadata, file.path(output_dirs$data, "metadata_enhanced.rds"))
saveRDS(logcpm_normalized, file.path(output_dirs$data, "logcpm_normalized.rds"))
saveRDS(sample_qc_results$pca, file.path(output_dirs$data, "pca_results.rds"))

# Save analysis parameters
analysis_params <- list(
  analysis_date = Sys.time(),
  r_version = R.version.string,
  packages_used = tryCatch(names(sessionInfo()$otherPkgs), error = function(e) "Package info not available"),
  outliers_removed = outliers_to_remove,
  original_samples = ncol(counts_original),
  final_samples = ncol(dge_normalized),
  original_genes = nrow(counts_original),
  final_genes = nrow(dge_normalized),
  normalization_method = "TMM",
  filtering_criteria = list(
    min_cpm = 1,
    min_samples = 3,
    edgeR_filter = TRUE,
    variance_filter = TRUE
  )
)
saveRDS(analysis_params, file.path(output_dirs$data, "analysis_parameters.rds"))

# Export key results as CSV
write.csv(metadata, file.path(output_dirs$data, "final_metadata.csv"), row.names = FALSE)
write.csv(outlier_results$summary, file.path(output_dirs$data, "outlier_analysis_summary.csv"), row.names = FALSE)

# ----------------------
# 14. COMPREHENSIVE SUMMARY REPORT
# ----------------------

cat("\n=== GENERATING COMPREHENSIVE SUMMARY REPORT ===\n")

create_comprehensive_report <- function(output_file) {
  sink(output_file)
  
  cat("===============================================\n")
  cat("ENHANCED BULK RNA-SEQ PREPROCESSING REPORT\n")
  cat("===============================================\n\n")
  
  cat("Analysis Information:\n")
  cat("- Date:", as.character(Sys.time()), "\n")
  cat("- R Version:", R.version.string, "\n")
  cat("- Script: Enhanced Processing and QC Pipeline\n\n")
  
  cat("Dataset Overview:\n")
  cat("- Original samples:", ncol(counts_original), "\n")
  cat("- Final samples:", ncol(dge_normalized), "\n")
  cat("- Samples removed:", ncol(counts_original) - ncol(dge_normalized), "\n")
  cat("- Original genes:", nrow(counts_original), "\n")
  cat("- Final genes:", nrow(dge_normalized), "\n")
  cat("- Genes filtered:", nrow(counts_original) - nrow(dge_normalized), "\n\n")
  
  cat("Outlier Analysis:\n")
  if(length(outliers_to_remove) > 0) {
    cat("- Samples removed:", paste(outliers_to_remove, collapse = ", "), "\n")
    cat("- Removal reasons: See outlier_detection_summary.csv\n")
  } else {
    cat("- No outliers detected\n")
  }
  cat("\n")
  
  cat("Experimental Design:\n")
  print(table(metadata$treatment, metadata$sex))
  cat("\n")
  
  if (!is.null(design_assessment$power_analysis)) {
    cat("Power Analysis (minimum group size:", design_assessment$min_group_size, "):\n")
    print(design_assessment$power_analysis)
    cat("\n")
  }
  
  cat("Quality Control Metrics:\n")
  cat("- Library size range:", round(range(metadata$lib_size)/1e6, 2), "million reads\n")
  cat("- Mean library size:", round(mean(metadata$lib_size)/1e6, 2), "million reads\n")
  cat("- Genes detected range:", range(metadata$genes_detected), "\n")
  cat("- Mean genes detected:", round(mean(metadata$genes_detected)), "\n")
  if(max(metadata$mito_percent) > 0) {
    cat("- Mitochondrial % range:", round(range(metadata$mito_percent), 2), "\n")
  }
  cat("\n")
  
  cat("Principal Component Analysis:\n")
  cat("- PC1 variance explained:", round(sample_qc_results$variance_explained[1], 2), "%\n")
  cat("- PC2 variance explained:", round(sample_qc_results$variance_explained[2], 2), "%\n")
  cat("- PC1+PC2 cumulative:", round(sum(sample_qc_results$variance_explained[1:2]), 2), "%\n\n")
  
  cat("Normalization:\n")
  cat("- Method: TMM (Trimmed Mean of M-values)\n")
  cat("- Alternative methods compared: RLE, Upper Quartile\n\n")
  
  cat("Sample-Level Analysis:\n")
  if(!is.null(batch_analysis)) {
    cat("- Treatment explains PC1 variance:", round(batch_analysis$variance_components$pc1_treatment * 100, 2), "%\n")
    cat("- Region explains PC1 variance:", round(batch_analysis$variance_components$pc1_region * 100, 2), "%\n")
  }
  cat("- Each sample treated as unique batch (geo_accession)\n\n")
  
  cat("Files Generated:\n")
  cat("Data Objects:\n")
  cat("- dge_normalized_final.rds: Final DGEList object\n")
  cat("- metadata_enhanced.rds: Complete sample metadata\n")
  cat("- logcpm_normalized.rds: Log-CPM normalized expression\n")
  cat("- pca_results.rds: PCA analysis results\n\n")
  
  cat("Statistical Reports:\n")
  cat("- comprehensive_gene_statistics.csv: Gene-level metrics\n")
  cat("- final_metadata.csv: Sample metadata\n")
  cat("- outlier_analysis_summary.csv: Outlier detection results\n")
  cat("- experimental_design_assessment.txt: Design and power analysis\n\n")
  
  cat("Quality Control Plots:\n")
  cat("- comprehensive_sample_qc.png: Multi-panel sample QC\n")
  cat("- publication_quality_pca.png: Publication-ready PCA\n")
  cat("- enhanced_correlation_heatmap.png: Sample correlation matrix\n")
  cat("- outlier_detection_plots.png: Outlier analysis visualization\n")
  cat("- normalization_comparison.png: Normalization method comparison\n")
  cat("- sample_level_analysis.png: Sample variability assessment\n")
  cat("\n")
  
  cat("Recommendations:\n")
  if(!is.null(design_assessment$power_analysis) && min(design_assessment$power_analysis$power) < 0.8) {
    cat("- Consider increasing sample size for better power\n")
  }
  if(length(outliers_to_remove) > 0) {
    cat("- Investigate reasons for sample outliers\n")
  }
  cat("- Data is ready for differential expression analysis\n")
  
  cat("\n===============================================\n")
  cat("ANALYSIS COMPLETE\n")
  cat("===============================================\n")
  
  sink()
}

create_comprehensive_report(file.path(output_dirs$reports, "comprehensive_analysis_report.txt"))

# Save session information
writeLines(capture.output(sessionInfo()), 
           file.path(output_dirs$reports, "session_info.txt"))

# ----------------------
# 15. FINAL STATUS REPORT
# ----------------------

cat("\n")
cat("===============================================\n")
cat("ENHANCED PREPROCESSING PIPELINE COMPLETE\n")
cat("===============================================\n")
cat("Analysis completed at:", as.character(Sys.time()), "\n")

cat("Final Dataset Summary:\n")
cat("- Samples:", ncol(dge_normalized), "\n")
cat("- Genes:", nrow(dge_normalized), "\n")
cat("- Treatment groups:", length(levels(metadata$treatment)), "\n")
cat("- Output directory:", output_dirs$base, "\n\n")

cat("Key Output Files:\n")
cat("- Main data object: dge_normalized_final.rds\n")
cat("- Enhanced metadata: metadata_enhanced.rds\n")
cat("- Comprehensive report: comprehensive_analysis_report.txt\n")
cat("- Publication plots: publication_quality_pca.png\n\n")

cat("Data is ready for:\n")
cat("- Differential expression analysis\n")
cat("- Pathway analysis\n")
cat("- Advanced statistical modeling\n")
cat("- Publication-quality visualization\n\n")

cat("Next recommended steps:\n")
cat("1. Load dge_normalized_final.rds for DE analysis\n")
cat("2. Review comprehensive_analysis_report.txt\n")
cat("3. Examine QC plots for any concerns\n")
cat("4. Proceed with statistical analysis\n")

cat("===============================================\n")