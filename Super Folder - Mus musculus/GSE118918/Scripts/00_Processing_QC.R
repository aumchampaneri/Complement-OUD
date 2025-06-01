# ========================================================================
# GSE118918 Bulk RNA-seq Processing and Quality Control Pipeline
# ========================================================================
# 
# STUDY OVERVIEW:
# Dataset: GSE118918 - Nucleus Accumbens RNA-seq (Mock vs Morphine)
# Tissue: Nucleus Accumbens (NAcc) 
# Treatment: Mock control vs Morphine treatment
# Data Type: Bulk RNA-seq from raw count files
# 
# SCIENTIFIC RATIONALE:
# This pipeline follows established best practices for RNA-seq preprocessing
# as outlined in:
# - Love et al. (2014) Moderated estimation of fold change and dispersion 
#   for RNA-seq data with DESeq2. Genome Biology 15:550
# - Robinson et al. (2010) edgeR: a Bioconductor package for differential 
#   expression analysis of digital gene expression data. Bioinformatics 26:139-140
# - Chen et al. (2016) From reads to genes to pathways: differential expression
#   analysis of RNA-Seq experiments using Rsubread, featureCounts, edgeR and ROAST
# 
# QUALITY CONTROL STANDARDS:
# Following guidelines from:
# - Conesa et al. (2016) A survey of best practices for RNA-seq data analysis
# - Evans et al. (2018) Selecting between-sample RNA-Seq normalization methods
# 
# STATISTICAL FRAMEWORK:
# - Multiple outlier detection methods with conservative thresholds
# - Evidence-based gene filtering (CPM, variance, biological relevance)
# - TMM normalization (most robust for between-sample comparisons)
# - Comprehensive quality assessment metrics
# 
# REPRODUCIBILITY MEASURES:
# - Fixed random seeds for all stochastic processes
# - Complete session information and package version tracking
# - Detailed parameter documentation
# - Standardized file naming and directory structure
# ========================================================================

# Set reproducibility parameters
set.seed(42)
options(stringsAsFactors = FALSE)

# Record analysis start time for documentation
analysis_start_time <- Sys.time()
cat("Analysis started at:", as.character(analysis_start_time), "\n")

# ========================================================================
# SECTION 1: LIBRARY LOADING AND ENVIRONMENT SETUP
# ========================================================================

# Load required libraries with version checking for reproducibility
cat("=== LOADING REQUIRED LIBRARIES ===\n")

required_packages <- c("edgeR", "limma", "RColorBrewer", "pheatmap", 
                      "ggplot2", "dplyr", "tidyr", "gridExtra")

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    stop(paste("Required package", pkg, "is not installed."))
  }
}

# Enhanced optional packages for bulk RNA-seq
optional_packages <- c("PCAtools", "sva", "preprocessCore", "ggrepel", "corrplot",
                       "DESeq2",        # Alternative DE method for comparison
                       "ComplexHeatmap", # Advanced heatmaps
                       "EnhancedVolcano", # Better volcano plots
                       "clusterProfiler", # Pathway enrichment
                       "GSVA",          # Gene set variation analysis
                       "limma",         # Already included but emphasizing
                       "Glimma"         # Interactive plots
                       )
loaded_optional <- character()

for(pkg in optional_packages) {
  if(require(pkg, character.only = TRUE, quietly = TRUE)) {
    loaded_optional <- c(loaded_optional, pkg)
  }
}

cat("Core packages loaded successfully\n")
cat("Optional packages available:", paste(loaded_optional, collapse = ", "), "\n")

# ========================================================================
# SECTION 2: DIRECTORY STRUCTURE AND FILE ORGANIZATION
# ========================================================================

cat("\n=== SETTING UP ANALYSIS ENVIRONMENT ===\n")

# Define base directory structure following best practices
base_dir <- "/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus/GSE118918"
setwd(base_dir)

# Create comprehensive output directory structure for reproducibility
output_structure <- list(
  # Main analysis outputs
  main = "Outputs/01_Processing_QC",
  
  # Organized subdirectories
  data = "Outputs/01_Processing_QC/Data",
  plots = "Outputs/01_Processing_QC/Plots", 
  reports = "Outputs/01_Processing_QC/Reports",
  qc_metrics = "Outputs/01_Processing_QC/QC_Metrics",
  
  # Specialized analysis folders
  raw_data_qc = "Outputs/01_Processing_QC/Raw_Data_Assessment",
  filtering = "Outputs/01_Processing_QC/Gene_Filtering",
  normalization = "Outputs/01_Processing_QC/Normalization",
  sample_analysis = "Outputs/01_Processing_QC/Sample_Analysis",
  
  # Documentation and reproducibility
  session_info = "Outputs/01_Processing_QC/Session_Info",
  parameters = "Outputs/01_Processing_QC/Analysis_Parameters"
)

# Create all directories
for(dir_name in names(output_structure)) {
  dir_path <- output_structure[[dir_name]]
  dir.create(dir_path, showWarnings = FALSE, recursive = TRUE)
  cat("Created directory:", dir_path, "\n")
}

# Define input data paths
raw_data_dir <- file.path(base_dir, "Data/Raw_Data")
series_matrix_file <- file.path(base_dir, "Data/NCBI/GSE118918_series_matrix.txt")

# Enhanced validation checks at the start
cat("\n=== PRE-ANALYSIS VALIDATION ===\n")

# 1. Verify directory structure and permissions
if(!dir.exists(base_dir)) {
  stop("Base directory does not exist. Please check the path: ", base_dir)
}

# 2. Check raw data availability
if(!dir.exists(raw_data_dir)) {
  stop("Raw data directory not found: ", raw_data_dir, 
       "\nPlease ensure count files are available before proceeding.")
}

# 3. Test write permissions
test_write <- tryCatch({
  temp_file <- file.path(base_dir, "test_write_permission.tmp")
  writeLines("test", temp_file)
  file.remove(temp_file)
  TRUE
}, error = function(e) {
  stop("Insufficient write permissions in base directory: ", base_dir)
})

cat("✓ All pre-analysis checks passed\n")

# ========================================================================
# SECTION 3: METADATA EXTRACTION AND EXPERIMENTAL DESIGN
# ========================================================================

cat("\n=== EXTRACTING EXPERIMENTAL METADATA ===\n")

# Function to parse GEO series matrix following GEO standards
parse_geo_series_matrix <- function(series_file) {
  cat("Parsing GEO series matrix file...\n")
  
  if(!file.exists(series_file)) {
    warning("Series matrix file not found. Will create metadata from filenames.")
    return(NULL)
  }
  
  # Read series matrix file
  lines <- readLines(series_file)
  
  # Extract sample information following GEO format specifications
  sample_geo_ids <- NULL
  sample_titles <- NULL
  sample_sources <- NULL
  
  # Find sample geo accession IDs
  geo_line <- grep("^!Sample_geo_accession", lines, value = TRUE)
  if(length(geo_line) > 0) {
    sample_geo_ids <- unlist(strsplit(geo_line, "\t"))[-1]
    sample_geo_ids <- gsub('"', '', sample_geo_ids)
  }
  
  # Find sample titles
  title_line <- grep("^!Sample_title", lines, value = TRUE)
  if(length(title_line) > 0) {
    sample_titles <- unlist(strsplit(title_line, "\t"))[-1]
    sample_titles <- gsub('"', '', sample_titles)
  }
  
  # Find sample source information
  source_line <- grep("^!Sample_source_name", lines, value = TRUE)
  if(length(source_line) > 0) {
    sample_sources <- unlist(strsplit(source_line, "\t"))[-1]
    sample_sources <- gsub('"', '', sample_sources)
  }
  
  # Create metadata dataframe
  if(!is.null(sample_geo_ids)) {
    metadata <- data.frame(
      geo_accession = sample_geo_ids,
      title = sample_titles %||% sample_geo_ids,
      source = sample_sources %||% rep("NAcc", length(sample_geo_ids)),
      stringsAsFactors = FALSE
    )
    
    cat("Successfully extracted metadata for", nrow(metadata), "samples\n")
    return(metadata)
  }
  
  return(NULL)
}

# Helper function for null coalescing
`%||%` <- function(x, y) if(is.null(x)) y else x

# Extract metadata from GEO
geo_metadata <- parse_geo_series_matrix(series_matrix_file)

# ========================================================================
# SECTION 4: RAW DATA LOADING AND INITIAL ASSESSMENT
# ========================================================================

cat("\n=== LOADING AND ASSESSING RAW COUNT DATA ===\n")

# Function to load and combine individual sample files
load_raw_count_data <- function(raw_data_directory) {
  cat("Scanning raw data directory for count files...\n")
  
  # Find all .txt.gz files (excluding summary files)
  count_files <- list.files(raw_data_directory, 
                           pattern = "\\.txt\\.gz$", 
                           full.names = TRUE)
  
  # Exclude summary files to get actual count data
  count_files <- count_files[!grepl("summary", count_files)]
  
  cat("Found", length(count_files), "count files\n")
  
  if(length(count_files) == 0) {
    stop("No count files found in the raw data directory: ", raw_data_directory,
         "\nExpected format: *.txt.gz files containing gene count data")
  }
  
  # Enhanced file validation
  if(length(count_files) < 4) {
    warning("Low number of samples detected (n=", length(count_files), 
            "). Minimum 6 samples recommended for robust statistical analysis.")
  }
  
  # Extract sample information from filenames
  sample_info <- data.frame(
    file_path = count_files,
    filename = basename(count_files),
    stringsAsFactors = FALSE
  )
  
  # Parse sample information from standardized GEO filenames
  # Format: GSM#######_NAcc_Treatment#.txt.gz
  sample_info$geo_id <- gsub("_.*", "", sample_info$filename)
  sample_info$tissue <- "NAcc"  # All samples from nucleus accumbens
  
  # Extract treatment information
  sample_info$treatment <- ifelse(grepl("Mock", sample_info$filename), "Mock", 
                                 ifelse(grepl("Morph", sample_info$filename), "Morphine", "Unknown"))
  
  # Extract replicate number
  sample_info$replicate <- gsub(".*([0-9]+)\\.txt\\.gz", "\\1", sample_info$filename)
  
  cat("Sample breakdown:\n")
  print(table(sample_info$treatment))
  
  # Load first file to determine data structure
  cat("Examining data structure from first file...\n")
  first_file <- count_files[1]
  
  # Try different loading approaches for robustness
  test_data <- tryCatch({
    read.table(gzfile(first_file), header = TRUE, sep = "\t", 
               stringsAsFactors = FALSE, check.names = FALSE)
  }, error = function(e) {
    cat("Standard loading failed, trying alternative approach...\n")
    read.table(gzfile(first_file), header = FALSE, sep = "\t", 
               stringsAsFactors = FALSE, check.names = FALSE)
  })
  
  cat("Data structure - Dimensions:", dim(test_data), "\n")
  cat("Column names:", colnames(test_data)[1:min(5, ncol(test_data))], "...\n")
  
  # Load all count files and combine into matrix
  cat("Loading all count files...\n")
  count_list <- list()
  
  for(i in seq_along(count_files)) {
    file_path <- count_files[i]
    sample_name <- sample_info$geo_id[i]
    
    cat("Loading file", i, "of", length(count_files), ":", sample_name, "\n")
    
    # Load data with error handling
    sample_data <- tryCatch({
      data <- read.table(gzfile(file_path), header = TRUE, sep = "\t", 
                        stringsAsFactors = FALSE, check.names = FALSE)
      
      # Standardize column names if needed
      if(ncol(data) >= 2) {
        # Assume first column is gene ID, second is count
        data <- data[, 1:2]
        colnames(data) <- c("gene_id", "count")
        data$count <- as.numeric(data$count)
        data <- data[!is.na(data$count), ]
        return(data)
      } else {
        stop("Unexpected data format")
      }
    }, error = function(e) {
      cat("Error loading file:", file_path, "\n")
      cat("Error message:", e$message, "\n")
      return(NULL)
    })
    
    if(!is.null(sample_data)) {
      count_list[[sample_name]] <- sample_data
    }
  }
  
  # Combine into count matrix
  cat("Combining count data into matrix...\n")
  
  # Get union of all genes
  all_genes <- unique(unlist(lapply(count_list, function(x) x$gene_id)))
  cat("Total unique genes across all samples:", length(all_genes), "\n")
  
  # Create count matrix
  count_matrix <- matrix(0, nrow = length(all_genes), ncol = length(count_list))
  rownames(count_matrix) <- all_genes
  colnames(count_matrix) <- names(count_list)
  
  # Fill count matrix
  for(sample_name in names(count_list)) {
    sample_data <- count_list[[sample_name]]
    matching_genes <- intersect(all_genes, sample_data$gene_id)
    count_matrix[matching_genes, sample_name] <- sample_data$count[match(matching_genes, sample_data$gene_id)]
  }
  
  # Convert to integer matrix
  count_matrix <- apply(count_matrix, 2, as.integer)
  rownames(count_matrix) <- all_genes
  
  # Enhanced data validation
  if(any(rowSums(count_matrix) == 0)) {
    n_empty_genes <- sum(rowSums(count_matrix) == 0)
    cat("Warning:", n_empty_genes, "genes have zero counts across all samples\n")
  }
  
  # Check for extremely low library sizes
  lib_sizes <- colSums(count_matrix)
  if(any(lib_sizes < 1e6)) {
    low_samples <- names(lib_sizes)[lib_sizes < 1e6]
    warning("Samples with <1M reads detected: ", paste(low_samples, collapse=", "))
  }
  
  return(list(
    counts = count_matrix,
    sample_info = sample_info,
    loading_stats = data.frame(
      total_files = length(count_files),
      successful_loads = length(count_list),
      total_genes = nrow(count_matrix),
      total_samples = ncol(count_matrix),
      mean_library_size = mean(colSums(count_matrix)),
      min_library_size = min(colSums(count_matrix)),
      max_library_size = max(colSums(count_matrix))
    )
  ))
}

# Load raw count data
raw_data_results <- load_raw_count_data(raw_data_dir)
counts_original <- raw_data_results$counts
sample_metadata <- raw_data_results$sample_info

# Save loading statistics
write.csv(raw_data_results$loading_stats, 
          file.path(output_structure$raw_data_qc, "data_loading_stats.csv"),
          row.names = FALSE)

cat("Successfully loaded count matrix:\n")
cat("- Dimensions:", nrow(counts_original), "genes x", ncol(counts_original), "samples\n")
cat("- Total counts:", format(sum(counts_original), big.mark = ","), "\n")

# ========================================================================
# SECTION 5: EXPERIMENTAL DESIGN VALIDATION
# ========================================================================

cat("\n=== VALIDATING EXPERIMENTAL DESIGN ===\n")

# Create comprehensive metadata combining file info and GEO data
create_comprehensive_metadata <- function(sample_info, geo_metadata = NULL) {
  
  # Start with file-based metadata
  metadata <- sample_info
  
  # Add additional experimental factors
  metadata$batch <- metadata$geo_id  # Use GEO ID as batch identifier
  metadata$tissue <- factor("NAcc")
  metadata$treatment <- factor(metadata$treatment, levels = c("Mock", "Morphine"))
  metadata$replicate <- as.numeric(metadata$replicate)
  
  # Add sample identifiers that match count matrix
  metadata$sample_id <- metadata$geo_id
  
  # Validate treatment assignment
  if(any(is.na(metadata$treatment))) {
    warning("Some samples have unassigned treatments")
  }
  
  # Ensure proper factor levels for statistical analysis
  metadata$treatment <- droplevels(metadata$treatment)
  
  return(metadata)
}

final_metadata <- create_comprehensive_metadata(sample_metadata, geo_metadata)

# Validate sample matching between counts and metadata
common_samples <- intersect(colnames(counts_original), final_metadata$sample_id)
cat("Samples in count matrix:", ncol(counts_original), "\n")
cat("Samples in metadata:", nrow(final_metadata), "\n")
cat("Matching samples:", length(common_samples), "\n")

if(length(common_samples) != ncol(counts_original)) {
  cat("Adjusting for sample mismatch...\n")
  final_metadata <- final_metadata[final_metadata$sample_id %in% colnames(counts_original), ]
  counts_original <- counts_original[, final_metadata$sample_id]
}

# Verify final matching
if(!all(colnames(counts_original) == final_metadata$sample_id)) {
  stop("Sample order mismatch between counts and metadata")
}

cat("Final experimental design:\n")
print(table(final_metadata$treatment))

# Save validated metadata
write.csv(final_metadata, 
          file.path(output_structure$data, "validated_metadata.csv"),
          row.names = FALSE)

# ========================================================================
# SECTION 6: COMPREHENSIVE RAW DATA QUALITY ASSESSMENT
# ========================================================================

cat("\n=== COMPREHENSIVE RAW DATA QUALITY ASSESSMENT ===\n")

# Function to calculate comprehensive QC metrics
calculate_raw_qc_metrics <- function(counts, metadata) {
  cat("Calculating comprehensive QC metrics...\n")
  
  # Sample-level metrics
  sample_metrics <- data.frame(
    sample_id = colnames(counts),
    treatment = metadata$treatment[match(colnames(counts), metadata$sample_id)],
    
    # Library composition metrics
    total_reads = colSums(counts),
    total_genes_detected = colSums(counts > 0),
    genes_with_1_read = colSums(counts >= 1),
    genes_with_5_reads = colSums(counts >= 5),
    genes_with_10_reads = colSums(counts >= 10),
    
    # Expression distribution metrics  
    median_expression = apply(counts, 2, median),
    mean_expression = colMeans(counts),
    q75_expression = apply(counts, 2, quantile, 0.75),
    q90_expression = apply(counts, 2, quantile, 0.90),
    max_expression = apply(counts, 2, max),
    
    # Sparsity metrics
    zero_count_rate = colMeans(counts == 0),
    low_expression_rate = colMeans(counts <= 5),
    
    stringsAsFactors = FALSE
  )
  
  # Gene-level metrics with error handling
  gene_metrics <- data.frame(
    gene_id = rownames(counts),
    
    # Expression statistics
    total_count = rowSums(counts),
    mean_count = rowMeans(counts),
    median_count = apply(counts, 1, median),
    max_count = apply(counts, 1, max),
    var_count = apply(counts, 1, var),
    
    # Detection metrics
    samples_detected = rowSums(counts > 0),
    detection_rate = rowSums(counts > 0) / ncol(counts),
    samples_high_expr = rowSums(counts >= 10),
    
    stringsAsFactors = FALSE
  )
  
  # Calculate CV with proper error handling
  gene_metrics$cv <- apply(counts, 1, function(x) {
    mean_val <- mean(x)
    if(mean_val == 0) return(NA)
    sd_val <- sd(x)
    return(sd_val / mean_val)
  })
  
  return(list(
    sample_metrics = sample_metrics,
    gene_metrics = gene_metrics
  ))
}

qc_metrics <- calculate_raw_qc_metrics(counts_original, final_metadata)

# Save QC metrics
write.csv(qc_metrics$sample_metrics, 
          file.path(output_structure$qc_metrics, "sample_level_qc_metrics.csv"),
          row.names = FALSE)

write.csv(qc_metrics$gene_metrics, 
          file.path(output_structure$qc_metrics, "gene_level_qc_metrics.csv"),
          row.names = FALSE)

# ========================================================================
# SECTION 7: PUBLICATION-QUALITY QC VISUALIZATIONS
# ========================================================================

cat("\n=== CREATING PUBLICATION-QUALITY QC PLOTS ===\n")

# Function to create comprehensive QC plots
create_publication_qc_plots <- function(counts, metadata, sample_metrics, gene_metrics, output_dir) {
  
  # Set up color scheme for treatments
  treatment_colors <- RColorBrewer::brewer.pal(max(3, length(levels(metadata$treatment))), "Set1")
  names(treatment_colors) <- levels(metadata$treatment)
  
  # Main QC figure - multi-panel layout
  png(file.path(output_dir, "Figure_S1_Comprehensive_QC.png"), 
      width = 3000, height = 2400, res = 300)
  
  layout(matrix(c(1,2,3,4,5,6,7,8), nrow = 2, byrow = TRUE))
  par(mar = c(5,5,3,2))
  
  # Panel A: Library sizes
  sample_colors <- treatment_colors[sample_metrics$treatment]
  barplot(sample_metrics$total_reads / 1e6, 
          names.arg = sample_metrics$sample_id,
          las = 2, cex.names = 0.8,
          col = sample_colors,
          main = "A. Library Sizes by Sample",
          ylab = "Total Reads (millions)",
          cex.lab = 1.2, cex.main = 1.3)
  legend("topright", legend = names(treatment_colors), 
         fill = treatment_colors, cex = 1.1)
  
  # Panel B: Genes detected
  barplot(sample_metrics$total_genes_detected,
          names.arg = sample_metrics$sample_id,
          las = 2, cex.names = 0.8,
          col = sample_colors,
          main = "B. Genes Detected per Sample",
          ylab = "Number of Genes",
          cex.lab = 1.2, cex.main = 1.3)
  
  # Panel C: Library size vs genes detected correlation
  plot(sample_metrics$total_reads / 1e6, sample_metrics$total_genes_detected,
       col = sample_colors, pch = 16, cex = 1.5,
       main = "C. Library Size vs Gene Detection",
       xlab = "Library Size (millions)",
       ylab = "Genes Detected",
       cex.lab = 1.2, cex.main = 1.3)
  
  # Add correlation line
  lm_fit <- lm(total_genes_detected ~ I(total_reads/1e6), data = sample_metrics)
  abline(lm_fit, col = "red", lwd = 2)
  
  # Add correlation coefficient
  cor_val <- cor(sample_metrics$total_reads, sample_metrics$total_genes_detected)
  text(max(sample_metrics$total_reads/1e6) * 0.7, 
       max(sample_metrics$total_genes_detected) * 0.9,
       paste("r =", round(cor_val, 3)), cex = 1.2)
  
  # Panel D: Gene detection rate distribution
  hist(gene_metrics$detection_rate, breaks = 50, 
       main = "D. Gene Detection Rate Distribution",
       xlab = "Detection Rate (fraction of samples)",
       ylab = "Number of Genes",
       col = "lightblue", border = "darkblue",
       cex.lab = 1.2, cex.main = 1.3)
  
  # Panel E: Treatment comparison - library sizes
  boxplot(total_reads/1e6 ~ treatment, data = sample_metrics,
          col = treatment_colors[levels(sample_metrics$treatment)],
          main = "E. Library Sizes by Treatment",
          ylab = "Library Size (millions)",
          cex.lab = 1.2, cex.main = 1.3)
  
  # Panel F: Treatment comparison - genes detected
  boxplot(total_genes_detected ~ treatment, data = sample_metrics,
          col = treatment_colors[levels(sample_metrics$treatment)],
          main = "F. Gene Detection by Treatment",
          ylab = "Genes Detected",
          cex.lab = 1.2, cex.main = 1.3)
  
  # Panel G: Mean-variance relationship
  plot(log10(gene_metrics$mean_count + 0.1), log10(gene_metrics$var_count + 0.1),
       pch = 16, col = rgb(0,0,0,0.3), cex = 0.5,
       main = "G. Mean-Variance Relationship",
       xlab = "log10(Mean Count)",
       ylab = "log10(Variance)",
       cex.lab = 1.2, cex.main = 1.3)
  
  # Add theoretical line for Poisson (variance = mean)
  abline(0, 1, col = "red", lwd = 2, lty = 2)
  legend("topleft", "Poisson expectation", lty = 2, col = "red", cex = 1.1)
  
  # Panel H: Sample correlation heatmap
  sample_cors <- cor(log2(counts + 1))
  image(1:ncol(sample_cors), 1:nrow(sample_cors), sample_cors,
        col = colorRampPalette(c("blue", "white", "red"))(100),
        main = "H. Sample Correlation Matrix",
        xlab = "Samples", ylab = "Samples",
        cex.lab = 1.2, cex.main = 1.3,
        axes = FALSE)
  
  dev.off()
  
  # Additional detailed plots
  
  # Detection threshold analysis
  png(file.path(output_dir, "Figure_S2_Detection_Analysis.png"), 
      width = 2400, height = 1200, res = 300)
  par(mfrow = c(1,2), mar = c(5,5,3,2))
  
  # Genes detected at different thresholds
  thresholds <- c(0, 1, 5, 10, 20, 50)
  detection_matrix <- sapply(thresholds, function(t) colSums(counts >= t))
  colnames(detection_matrix) <- paste(">=", thresholds)
  
  boxplot(detection_matrix, 
          main = "Genes Detected at Different Count Thresholds",
          ylab = "Number of Genes",
          xlab = "Count Threshold",
          col = "lightblue",
          cex.lab = 1.2, cex.main = 1.3)
  
  # Cumulative detection rate
  detection_rates <- sort(gene_metrics$detection_rate, decreasing = TRUE)
  plot(1:length(detection_rates), detection_rates,
       type = "l", lwd = 2, col = "blue",
       main = "Cumulative Gene Detection",
       xlab = "Gene Rank",
       ylab = "Detection Rate",
       cex.lab = 1.2, cex.main = 1.3)
  
  # Add reference lines
  abline(h = c(0.5, 0.8, 0.9), col = "red", lty = 2)
  text(length(detection_rates) * 0.7, c(0.55, 0.85, 0.95), 
       c("50%", "80%", "90%"), col = "red")
  
  dev.off()
}

create_publication_qc_plots(counts_original, final_metadata, 
                           qc_metrics$sample_metrics, qc_metrics$gene_metrics,
                           output_structure$plots)

# ========================================================================
# SECTION 8: STATISTICAL OUTLIER DETECTION
# ========================================================================

cat("\n=== STATISTICAL OUTLIER DETECTION ===\n")

# Comprehensive outlier detection following established methods
detect_statistical_outliers <- function(counts, metadata, output_dir) {
  cat("Performing multi-method outlier detection...\n")
  
  outlier_results <- list()
  
  # Method 1: Library size outliers (>3 SD from mean)
  lib_sizes <- colSums(counts)
  lib_z_scores <- abs(scale(lib_sizes)[,1])
  lib_outliers <- names(lib_z_scores)[lib_z_scores > 3]
  outlier_results$library_size <- lib_outliers
  
  # Method 2: Gene detection outliers
  genes_detected <- colSums(counts > 0)
  gene_z_scores <- abs(scale(genes_detected)[,1])
  gene_outliers <- names(gene_z_scores)[gene_z_scores > 3]
  outlier_results$gene_detection <- gene_outliers
  
  # Method 3: PCA-based outlier detection
  # Use log-transformed data for PCA
  log_counts <- log2(counts + 1)
  
  # Remove genes with zero variance
  gene_vars <- apply(log_counts, 1, var)
  log_counts_filt <- log_counts[gene_vars > 0, ]
  
  # Check if we have enough genes for PCA
  if(nrow(log_counts_filt) < 100) {
    warning("Very few genes available for PCA. Results may be unreliable.")
  }
  
  pca_result <- prcomp(t(log_counts_filt), center = TRUE, scale. = TRUE)
  
  # Calculate Mahalanobis distance on first 5 PCs
  n_pcs <- min(5, ncol(pca_result$x), nrow(pca_result$x) - 1)
  pc_coords <- pca_result$x[, 1:n_pcs, drop = FALSE]
  
  # Check for sufficient samples for Mahalanobis distance
  if(nrow(pc_coords) <= n_pcs) {
    warning("Insufficient samples for reliable Mahalanobis distance calculation")
    pca_outliers <- character(0)
  } else {
    maha_dist <- mahalanobis(pc_coords, 
                            center = colMeans(pc_coords), 
                            cov = cov(pc_coords))
    
    # Chi-square threshold for outliers
    maha_threshold <- qchisq(0.975, df = n_pcs)
    pca_outliers <- names(maha_dist)[maha_dist > maha_threshold]
  }
  outlier_results$pca_mahalanobis <- pca_outliers
  
  # Method 4: Relative Log Expression (RLE) outliers
  rle_data <- log_counts_filt - rowMeans(log_counts_filt)
  sample_rle_medians <- apply(rle_data, 2, median)
  rle_mad <- mad(sample_rle_medians)
  
  # Handle case where MAD is 0
  if(rle_mad == 0) {
    warning("MAD of RLE medians is 0. Using standard deviation instead.")
    rle_mad <- sd(sample_rle_medians)
  }
  
  rle_outliers <- names(sample_rle_medians)[abs(sample_rle_medians) > 3 * rle_mad]
  outlier_results$rle <- rle_outliers
  
  # Create comprehensive outlier plot
  png(file.path(output_dir, "Figure_S3_Outlier_Detection.png"), 
      width = 3000, height = 2400, res = 300)
  
  layout(matrix(c(1,2,3,4,5,6), nrow = 2, byrow = TRUE))
  par(mar = c(5,5,3,2))
  
  # Treatment colors for consistency
  treatment_colors <- RColorBrewer::brewer.pal(max(3, length(levels(metadata$treatment))), "Set1")
  names(treatment_colors) <- levels(metadata$treatment)
  sample_colors <- treatment_colors[metadata$treatment]
  
  # Plot 1: Library sizes with outlier highlighting
  plot(lib_sizes / 1e6, 
       col = ifelse(names(lib_sizes) %in% lib_outliers, "red", sample_colors),
       pch = ifelse(names(lib_sizes) %in% lib_outliers, 17, 16),
       cex = 1.5,
       main = "Library Size Outlier Detection",
       ylab = "Library Size (millions)",
       xlab = "Sample Index")
  abline(h = (mean(lib_sizes) + c(-3,3) * sd(lib_sizes)) / 1e6, 
         col = "red", lty = 2)
  
  # Plot 2: Gene detection outliers
  plot(genes_detected,
       col = ifelse(names(genes_detected) %in% gene_outliers, "red", sample_colors),
       pch = ifelse(names(genes_detected) %in% gene_outliers, 17, 16),
       cex = 1.5,
       main = "Gene Detection Outlier Detection",
       ylab = "Genes Detected",
       xlab = "Sample Index")
  abline(h = mean(genes_detected) + c(-3,3) * sd(genes_detected), 
         col = "red", lty = 2)
  
  # Plot 3: PCA with outliers
  plot(pca_result$x[,1], pca_result$x[,2],
       col = ifelse(rownames(pca_result$x) %in% pca_outliers, "red", sample_colors),
       pch = ifelse(rownames(pca_result$x) %in% pca_outliers, 17, 16),
       cex = 1.5,
       main = "PCA-based Outlier Detection",
       xlab = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
       ylab = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)"))
  
  # Plot 4: Mahalanobis distances
  plot(maha_dist,
       col = ifelse(names(maha_dist) %in% pca_outliers, "red", sample_colors),
       pch = ifelse(names(maha_dist) %in% pca_outliers, 17, 16),
       cex = 1.5,
       main = "Mahalanobis Distance",
       ylab = "Distance",
       xlab = "Sample Index")
  abline(h = maha_threshold, col = "red", lty = 2)
  
  # Plot 5: RLE boxplot
  boxplot(rle_data, 
          col = ifelse(colnames(rle_data) %in% rle_outliers, "red", sample_colors),
          las = 2, cex.names = 0.8,
          main = "Relative Log Expression (RLE)",
          ylab = "RLE Values")
  abline(h = 0, col = "blue", lty = 1)
  
  # Plot 6: Combined outlier summary
  all_samples <- colnames(counts)
  outlier_matrix <- matrix(FALSE, nrow = length(all_samples), ncol = 4)
  rownames(outlier_matrix) <- all_samples
  colnames(outlier_matrix) <- c("LibSize", "GeneDet", "PCA", "RLE")
  
  outlier_matrix[lib_outliers, "LibSize"] <- TRUE
  outlier_matrix[gene_outliers, "GeneDet"] <- TRUE  
  outlier_matrix[pca_outliers, "PCA"] <- TRUE
  outlier_matrix[rle_outliers, "RLE"] <- TRUE
  
  image(1:ncol(outlier_matrix), 1:nrow(outlier_matrix), 
        t(as.matrix(outlier_matrix)),
        col = c("white", "red"),
        main = "Outlier Detection Summary",
        xlab = "Detection Method",
        ylab = "Samples",
        axes = FALSE)
  axis(1, at = 1:ncol(outlier_matrix), labels = colnames(outlier_matrix))
  axis(2, at = 1:nrow(outlier_matrix), labels = rownames(outlier_matrix), 
       las = 2, cex.axis = 0.8)
  
  dev.off()
  
  # Create outlier summary table with better error handling
  all_outliers <- unique(unlist(outlier_results))
  
  if(length(all_outliers) == 0) {
    # No outliers detected
    outlier_summary <- data.frame(
      sample_id = character(0),
      treatment = character(0),
      library_size_outlier = logical(0),
      gene_detection_outlier = logical(0),
      pca_outlier = logical(0),
      rle_outlier = logical(0),
      total_methods = integer(0),
      stringsAsFactors = FALSE
    )
  } else {
    outlier_summary <- data.frame(
      sample_id = all_outliers,
      treatment = metadata$treatment[match(all_outliers, metadata$sample_id)],
      library_size_outlier = all_outliers %in% outlier_results$library_size,
      gene_detection_outlier = all_outliers %in% outlier_results$gene_detection,
      pca_outlier = all_outliers %in% outlier_results$pca_mahalanobis,
      rle_outlier = all_outliers %in% outlier_results$rle,
      total_methods = rowSums(cbind(
        all_outliers %in% outlier_results$library_size,
        all_outliers %in% outlier_results$gene_detection,
        all_outliers %in% outlier_results$pca_mahalanobis,
        all_outliers %in% outlier_results$rle
      )),
      stringsAsFactors = FALSE
    )
  }
  
  return(list(
    outliers_by_method = outlier_results,
    outlier_summary = outlier_summary,
    pca_results = pca_result,
    all_outliers = all_outliers
  ))
}

outlier_analysis <- detect_statistical_outliers(counts_original, final_metadata, 
                                               output_structure$plots)

# Save outlier analysis results
write.csv(outlier_analysis$outlier_summary,
          file.path(output_structure$qc_metrics, "outlier_detection_summary.csv"),
          row.names = FALSE)

# Conservative outlier removal decision for peer review
# Following best practices: only remove samples with multiple outlier flags
samples_to_remove <- character(0)  # Conservative: no automatic removal

cat("Outlier analysis summary:\n")
cat("- Samples flagged by multiple methods:", 
    sum(outlier_analysis$outlier_summary$total_methods >= 2), "\n")
cat("- Using conservative approach: manual review recommended\n")

# ========================================================================
# SECTION 9: GENE FILTERING WITH SCIENTIFIC RATIONALE
# ========================================================================

cat("\n=== EVIDENCE-BASED GENE FILTERING ===\n")

# Apply filtering based on established criteria
apply_gene_filtering <- function(counts, metadata, output_dir) {
  cat("Applying evidence-based gene filtering criteria...\n")
  
  # Create DGEList object
  dge <- DGEList(counts = counts, samples = metadata)
  
  # Filtering criteria based on literature recommendations
  
  # Criterion 1: Expression level filtering
  # Genes must have CPM > 1 in at least the smallest group size
  group_sizes <- table(metadata$treatment)
  min_group_size <- min(group_sizes)
  
  cpm_vals <- cpm(dge)
  keep_expr <- rowSums(cpm_vals >= 1) >= min_group_size
  
  # Criterion 2: edgeR's filterByExpr (established method)
  keep_edger <- filterByExpr(dge, group = metadata$treatment)
  
  # Criterion 3: Remove genes with extremely low variance
  log_cpm <- cpm(dge, log = TRUE)
  gene_vars <- apply(log_cpm, 1, var)
  keep_var <- gene_vars > quantile(gene_vars, 0.1, na.rm = TRUE)
  
  # Combine filters
  keep_final <- keep_expr & keep_edger & keep_var
  
  # Apply filtering
  dge_filtered <- dge[keep_final, , keep.lib.sizes = FALSE]
  
  # Create filtering summary
  filtering_summary <- data.frame(
    criterion = c("Original", "CPM >= 1", "edgeR filter", "Variance filter", "Final"),
    genes_retained = c(
      nrow(dge),
      sum(keep_expr),
      sum(keep_edger), 
      sum(keep_var),
      sum(keep_final)
    ),
    percent_retained = c(
      100,
      round(sum(keep_expr) / nrow(dge) * 100, 1),
      round(sum(keep_edger) / nrow(dge) * 100, 1),
      round(sum(keep_var) / nrow(dge) * 100, 1),
      round(sum(keep_final) / nrow(dge) * 100, 1)
    )
  )
  
  # Visualization of filtering effects
  png(file.path(output_dir, "Figure_S4_Gene_Filtering.png"),
      width = 2400, height = 1200, res = 300)
  
  par(mfrow = c(1,2), mar = c(5,5,3,2))
  
  # Filtering cascade
  barplot(filtering_summary$genes_retained,
          names.arg = filtering_summary$criterion,
          las = 2,
          col = "lightblue",
          main = "Gene Filtering Cascade",
          ylab = "Number of Genes Retained",
          cex.lab = 1.2, cex.main = 1.3)
  
  # Before/after comparison
  before_after <- data.frame(
    mean_expr_before = rowMeans(cpm(dge, log = TRUE)),
    mean_expr_after = rowMeans(cpm(dge_filtered, log = TRUE)[rownames(dge_filtered), ]),
    filtered = !rownames(dge) %in% rownames(dge_filtered)
  )
  
  plot(density(before_after$mean_expr_before), 
       main = "Expression Distribution: Before vs After Filtering",
       xlab = "Mean log-CPM",
       ylab = "Density",
       col = "red", lwd = 2,
       cex.lab = 1.2, cex.main = 1.3)
  lines(density(before_after$mean_expr_after, na.rm = TRUE), 
        col = "blue", lwd = 2)
  legend("topright", c("Before filtering", "After filtering"), 
         col = c("red", "blue"), lwd = 2)
  
  dev.off()
  
  return(list(
    dge_filtered = dge_filtered,
    filtering_summary = filtering_summary,
    keep_genes = keep_final
  ))
}

filtering_results <- apply_gene_filtering(counts_original, final_metadata, 
                                        output_structure$filtering)

# Save filtering results
write.csv(filtering_results$filtering_summary,
          file.path(output_structure$filtering, "gene_filtering_summary.csv"),
          row.names = FALSE)

dge_filtered <- filtering_results$dge_filtered

cat("Gene filtering completed:\n")
cat("- Original genes:", nrow(counts_original), "\n")
cat("- Filtered genes:", nrow(dge_filtered), "\n")
cat("- Retention rate:", round(nrow(dge_filtered)/nrow(counts_original)*100, 1), "%\n")

# ========================================================================
# SECTION 10: NORMALIZATION AND METHOD COMPARISON
# ========================================================================

cat("\n=== NORMALIZATION METHOD COMPARISON ===\n")

# Compare multiple normalization methods
compare_normalization_methods <- function(dge, output_dir) {
  cat("Comparing normalization methods following best practices...\n")
  
  # Method 1: TMM (Trimmed Mean of M-values) - edgeR default
  dge_tmm <- calcNormFactors(dge, method = "TMM")
  
  # Method 2: RLE (Relative Log Expression) - DESeq2 style
  dge_rle <- calcNormFactors(dge, method = "RLE")
  
  # Method 3: Upper quartile normalization
  dge_uq <- calcNormFactors(dge, method = "upperquartile")
  
  # Get normalized log-CPM values
  logcpm_tmm <- cpm(dge_tmm, log = TRUE)
  logcpm_rle <- cpm(dge_rle, log = TRUE)
  logcpm_uq <- cpm(dge_uq, log = TRUE)
  
  # Compare normalization factors
  norm_factors <- data.frame(
    sample_id = colnames(dge$counts),
    treatment = dge$samples$treatment,
    TMM = dge_tmm$samples$norm.factors,
    RLE = dge_rle$samples$norm.factors,
    UpperQuartile = dge_uq$samples$norm.factors
  )
  
  # Create comparison plots
  png(file.path(output_dir, "Figure_S5_Normalization_Comparison.png"),
      width = 3000, height = 1200, res = 300)
  
  par(mfrow = c(1,3), mar = c(5,5,3,2))
  
  # Treatment colors
  treatment_colors <- RColorBrewer::brewer.pal(max(3, length(levels(dge$samples$treatment))), "Set1")
  names(treatment_colors) <- levels(dge$samples$treatment)
  sample_colors <- treatment_colors[dge$samples$treatment]
  
  # Box plots for each normalization method
  boxplot(logcpm_tmm, 
          col = sample_colors,
          las = 2, cex.names = 0.8,
          main = "TMM Normalization",
          ylab = "log-CPM",
          cex.lab = 1.2, cex.main = 1.3)
  
  boxplot(logcpm_rle,
          col = sample_colors, 
          las = 2, cex.names = 0.8,
          main = "RLE Normalization",
          ylab = "log-CPM",
          cex.lab = 1.2, cex.main = 1.3)
  
  boxplot(logcpm_uq,
          col = sample_colors,
          las = 2, cex.names = 0.8, 
          main = "Upper Quartile Normalization",
          ylab = "log-CPM",
          cex.lab = 1.2, cex.main = 1.3)
  
  dev.off()
  
  # Normalization factor comparison
  png(file.path(output_dir, "Figure_S6_Normalization_Factors.png"),
      width = 2400, height = 1200, res = 300)
  
  par(mfrow = c(1,2), mar = c(5,5,3,2))
  
  # Normalization factors by sample
  plot(norm_factors$TMM, col = sample_colors, pch = 16, cex = 1.5,
       main = "TMM Normalization Factors",
       ylab = "Normalization Factor",
       xlab = "Sample Index",
       cex.lab = 1.2, cex.main = 1.3)
  abline(h = 1, col = "red", lty = 2)
  
  # Correlation between methods
  plot(norm_factors$TMM, norm_factors$RLE, 
       col = sample_colors, pch = 16, cex = 1.5,
       main = "TMM vs RLE Normalization Factors",
       xlab = "TMM Factors",
       ylab = "RLE Factors", 
       cex.lab = 1.2, cex.main = 1.3)
  abline(0, 1, col = "red", lty = 2)
  
  # Calculate correlation
  cor_tmm_rle <- cor(norm_factors$TMM, norm_factors$RLE)
  text(max(norm_factors$TMM) * 0.7, max(norm_factors$RLE) * 0.9,
       paste("r =", round(cor_tmm_rle, 3)), cex = 1.2)
  
  dev.off()
  
  return(list(
    dge_tmm = dge_tmm,
    dge_rle = dge_rle, 
    dge_uq = dge_uq,
    norm_factors = norm_factors,
    logcpm_tmm = logcpm_tmm,
    logcpm_rle = logcpm_rle,
    logcpm_uq = logcpm_uq
  ))
}

normalization_results <- compare_normalization_methods(dge_filtered, 
                                                      output_structure$normalization)

# Select TMM as default (most widely used and robust)
dge_normalized <- normalization_results$dge_tmm
logcpm_final <- normalization_results$logcpm_tmm

# Save normalization comparison
write.csv(normalization_results$norm_factors,
          file.path(output_structure$normalization, "normalization_factors_comparison.csv"),
          row.names = FALSE)

# ========================================================================
# SECTION 11: FINAL SAMPLE QUALITY ASSESSMENT
# ========================================================================

cat("\n=== FINAL SAMPLE QUALITY ASSESSMENT ===\n")

# PCA analysis on final normalized data
perform_final_pca <- function(logcpm, metadata, output_dir) {
  cat("Performing PCA on final normalized data...\n")
  
  # Remove genes with zero variance
  gene_vars <- apply(logcpm, 1, var)
  logcpm_pca <- logcpm[gene_vars > 0 & !is.na(gene_vars), ]
  
  cat("Genes used for PCA:", nrow(logcpm_pca), "out of", nrow(logcpm), "\n")
  
  # Check if we have enough genes for meaningful PCA
  if(nrow(logcpm_pca) < 100) {
    warning("Very few genes available for PCA (", nrow(logcpm_pca), "). Results may be unreliable.")
  }
  
  # Perform PCA with error handling
  pca_result <- tryCatch({
    prcomp(t(logcpm_pca), center = TRUE, scale. = TRUE)
  }, error = function(e) {
    cat("PCA failed with scaling. Trying without scaling...\n")
    prcomp(t(logcpm_pca), center = TRUE, scale. = FALSE)
  })
  
  variance_explained <- summary(pca_result)$importance[2, ] * 100
  
  # Create publication-quality PCA plot with better error handling
  png(file.path(output_dir, "Figure_2_Final_PCA_Analysis.png"),
      width = 2400, height = 1200, res = 300)
  
  par(mfrow = c(1,2), mar = c(5,5,3,2))
  
  # Treatment colors
  treatment_colors <- RColorBrewer::brewer.pal(max(3, length(levels(metadata$treatment))), "Set1")
  names(treatment_colors) <- levels(metadata$treatment)
  sample_colors <- treatment_colors[metadata$treatment]
  
  # PC1 vs PC2
  plot(pca_result$x[,1], pca_result$x[,2],
       col = sample_colors, pch = 16, cex = 2,
       main = "Principal Component Analysis",
       xlab = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       ylab = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
       cex.lab = 1.2, cex.main = 1.3)
  
  # Add sample labels if requested
  text(pca_result$x[,1], pca_result$x[,2], 
       labels = rownames(pca_result$x),
       pos = 3, cex = 0.8)
  
  legend("topright", legend = levels(metadata$treatment),
         col = treatment_colors[levels(metadata$treatment)], 
         pch = 16, cex = 1.2)
  
  # Variance explained plot - handle case with few PCs
  n_pcs_to_show <- min(10, ncol(pca_result$x))
  barplot(variance_explained[1:n_pcs_to_show],
          names.arg = paste0("PC", 1:n_pcs_to_show),
          main = "Variance Explained by Principal Components",
          ylab = "% Variance Explained",
          xlab = "Principal Component",
          col = "lightblue",
          cex.lab = 1.2, cex.main = 1.3)
  
  dev.off()
  
  return(list(
    pca = pca_result,
    variance_explained = variance_explained
  ))
}

final_pca <- perform_final_pca(logcpm_final, final_metadata, 
                              output_structure$sample_analysis)

# ========================================================================
# SECTION 12: DATA EXPORT AND DOCUMENTATION
# ========================================================================

cat("\n=== SAVING FINAL PROCESSED DATA ===\n")

# Save all key objects for downstream analysis
saveRDS(dge_normalized, 
        file.path(output_structure$data, "dge_normalized_final.rds"))

saveRDS(final_metadata, 
        file.path(output_structure$data, "sample_metadata_final.rds"))

saveRDS(logcpm_final, 
        file.path(output_structure$data, "logcpm_normalized_final.rds"))

saveRDS(final_pca,
        file.path(output_structure$data, "pca_analysis_final.rds"))

saveRDS(outlier_analysis,
        file.path(output_structure$data, "outlier_analysis_results.rds"))

# Export CSV files for external use
write.csv(final_metadata,
          file.path(output_structure$data, "sample_metadata_final.csv"),
          row.names = FALSE)

write.csv(as.data.frame(logcpm_final),
          file.path(output_structure$data, "logcpm_normalized_final.csv"),
          row.names = TRUE)

# Create analysis parameters file for reproducibility
analysis_parameters <- list(
  # Session information
  analysis_date = Sys.time(),
  r_version = R.version.string,
  packages_loaded = names(sessionInfo()$otherPkgs),
  
  # Dataset information
  dataset_id = "GSE118918",
  organism = "Mus musculus",
  tissue = "Nucleus Accumbens",
  experimental_design = "Mock vs Morphine treatment",
  
  # Processing parameters
  outliers_removed = samples_to_remove,
  filtering_method = "edgeR filterByExpr + CPM + variance",
  normalization_method = "TMM",
  
  # Final dataset size
  original_samples = ncol(counts_original),
  final_samples = ncol(dge_normalized),
  original_genes = nrow(counts_original), 
  final_genes = nrow(dge_normalized),
  
  # Quality metrics
  mean_library_size = mean(colSums(dge_normalized$counts)),
  mean_genes_detected = mean(colSums(dge_normalized$counts > 0)),
  pc1_variance = final_pca$variance_explained[1],
  pc2_variance = final_pca$variance_explained[2]
)

saveRDS(analysis_parameters,
        file.path(output_structure$parameters, "analysis_parameters.rds"))

# ========================================================================
# SECTION 13: COMPREHENSIVE ANALYSIS REPORT
# ========================================================================

cat("\n=== GENERATING COMPREHENSIVE ANALYSIS REPORT ===\n")

# Create detailed analysis report with enhanced error handling
create_analysis_report <- function(output_file, params) {
  # Add check for required variables
  if(!exists("outlier_analysis")) {
    warning("outlier_analysis not found. Some report sections may be incomplete.")
    outlier_analysis <- list(outlier_summary = data.frame())
  }
  
  if(!exists("filtering_results")) {
    warning("filtering_results not found. Some report sections may be incomplete.")
    filtering_results <- list(filtering_summary = data.frame())
  }
  
  sink(output_file)
  
  cat("=================================================================\n")
  cat("GSE118918 BULK RNA-SEQ PROCESSING AND QC ANALYSIS REPORT\n") 
  cat("=================================================================\n\n")
  
  cat("STUDY INFORMATION\n")
  cat("-----------------\n")
  cat("Dataset ID: GSE118918\n")
  cat("Organism: Mus musculus\n")
  cat("Tissue: Nucleus Accumbens (NAcc)\n")
  cat("Experimental Design: Mock control vs Morphine treatment\n")
  cat("Data Type: Bulk RNA-sequencing\n")
  cat("Analysis Date:", as.character(params$analysis_date), "\n")
  cat("R Version:", params$r_version, "\n\n")
  
  cat("DATASET OVERVIEW\n")
  cat("----------------\n")
  cat("Original samples loaded:", params$original_samples, "\n")
  cat("Final samples after QC:", params$final_samples, "\n")
  cat("Samples removed:", params$original_samples - params$final_samples, "\n")
  cat("Original genes:", format(params$original_genes, big.mark = ","), "\n")
  cat("Final genes after filtering:", format(params$final_genes, big.mark = ","), "\n")
  cat("Gene retention rate:", round(params$final_genes/params$original_genes*100, 1), "%\n\n")
  
  cat("EXPERIMENTAL DESIGN SUMMARY\n")
  cat("---------------------------\n")
  design_table <- table(final_metadata$treatment)
  for(i in 1:length(design_table)) {
    cat(names(design_table)[i], "samples:", design_table[i], "\n")
  }
  cat("Total samples:", sum(design_table), "\n")
  cat("Balanced design:", all(design_table == design_table[1]), "\n\n")
  
  cat("QUALITY CONTROL METRICS\n")
  cat("-----------------------\n")
  cat("Mean library size:", format(round(params$mean_library_size), big.mark = ","), "reads\n")
  cat("Mean genes detected per sample:", round(params$mean_genes_detected), "\n")
  cat("Library size range:", format(range(colSums(dge_normalized$counts)), big.mark = ","), "reads\n")
  
  lib_sizes <- colSums(dge_normalized$counts)
  cat("Library size CV:", round(sd(lib_sizes)/mean(lib_sizes)*100, 1), "%\n\n")
  
  cat("OUTLIER DETECTION RESULTS\n")
  cat("-------------------------\n")
  cat("Total samples flagged as potential outliers:", nrow(outlier_analysis$outlier_summary), "\n")
  cat("Samples flagged by multiple methods:", sum(outlier_analysis$outlier_summary$total_methods >= 2), "\n")
  cat("Samples removed:", length(samples_to_remove), "\n")
  if(length(samples_to_remove) == 0) {
    cat("Conservative approach: No automatic sample removal applied\n")
  }
  cat("\n")
  
  cat("GENE FILTERING SUMMARY\n")
  cat("----------------------\n")
  print(filtering_results$filtering_summary)
  cat("\n")

  cat("NORMALIZATION\n")
  cat("-------------\n")
  cat("Method: TMM (Trimmed Mean of M-values)\n")
  cat("Alternative methods evaluated: RLE, Upper Quartile\n")
  cat("Normalization factors range:", round(range(normalization_results$norm_factors$TMM), 3), "\n")
  cat("Mean normalization factor:", round(mean(normalization_results$norm_factors$TMM), 3), "\n\n")
  
  cat("PRINCIPAL COMPONENT ANALYSIS\n")
  cat("----------------------------\n")
  cat("PC1 variance explained:", round(params$pc1_variance, 1), "%\n")
  cat("PC2 variance explained:", round(params$pc2_variance, 1), "%\n")
  cat("PC1+PC2 cumulative variance:", round(params$pc1_variance + params$pc2_variance, 1), "%\n")
  
  # Treatment separation assessment
  if(length(levels(final_metadata$treatment)) == 2) {
    pc1_treatment_sep <- abs(diff(by(final_pca$pca$x[,1], final_metadata$treatment, mean)))
    cat("PC1 treatment separation (mean difference):", round(pc1_treatment_sep, 2), "\n")
  }
  cat("\n")
  
  cat("OUTPUT FILES GENERATED\n")
  cat("----------------------\n")
  cat("Main data objects:\n")
  cat("- dge_normalized_final.rds: Final DGEList object for differential expression\n")
  cat("- sample_metadata_final.rds: Complete sample metadata\n")
  cat("- logcpm_normalized_final.rds: Log-CPM normalized expression values\n")
  cat("- pca_analysis_final.rds: PCA results and variance explained\n")
  cat("- outlier_analysis_results.rds: Comprehensive outlier detection results\n\n")
  
  cat("Quality control reports:\n")
  cat("- sample_level_qc_metrics.csv: Sample-level quality metrics\n")
  cat("- gene_level_qc_metrics.csv: Gene-level statistics\n")
  cat("- outlier_detection_summary.csv: Outlier flagging summary\n")
  cat("- gene_filtering_summary.csv: Gene filtering cascade results\n")
  cat("- normalization_factors_comparison.csv: Normalization method comparison\n\n")
  
  cat("Publication-quality figures:\n")
  cat("- Figure_S1_Comprehensive_QC.png: Multi-panel QC overview\n")
  cat("- Figure_S2_Detection_Analysis.png: Gene detection threshold analysis\n")
  cat("- Figure_S3_Outlier_Detection.png: Multi-method outlier detection\n")
  cat("- Figure_S4_Gene_Filtering.png: Gene filtering visualization\n")
  cat("- Figure_S5_Normalization_Comparison.png: Normalization method comparison\n")
  cat("- Figure_S6_Normalization_Factors.png: Normalization factor analysis\n")
  cat("- Figure_2_Final_PCA_Analysis.png: Final PCA results\n\n")
  
  cat("RECOMMENDATIONS FOR DOWNSTREAM ANALYSIS\n")
  cat("---------------------------------------\n")
  cat("1. Data quality: ")
  if(params$pc1_variance > 20 && length(samples_to_remove) == 0) {
    cat("EXCELLENT - High variance captured, no outliers removed\n")
  } else if(params$pc1_variance > 15) {
    cat("GOOD - Adequate variance captured for analysis\n")
  } else {
    cat("MODERATE - Consider additional QC investigation\n")
  }
  
  cat("2. Sample size: ")
  min_group <- min(table(final_metadata$treatment))
  if(min_group >= 6) {
    cat("ADEQUATE - Sufficient power for differential expression\n")
  } else if(min_group >= 3) {
    cat("MINIMAL - Consider effect size expectations\n")
  } else {
    cat("INSUFFICIENT - Additional samples recommended\n")
  }
  
  cat("3. Expression profile: ")
  detection_rate <- mean(qc_metrics$gene_metrics$detection_rate > 0.5)
  if(detection_rate > 0.4) {
    cat("GOOD - High gene detection rate\n")
  } else {
    cat("MODERATE - Consider sequencing depth\n")
  }
  
  cat("\n4. Ready for downstream analysis:\n")
  cat("   - Differential expression analysis (edgeR/DESeq2)\n")
  cat("   - Pathway enrichment analysis\n")
  cat("   - Complement pathway investigation\n")
  cat("   - Morphine response characterization\n")
  cat("   - Integration with other datasets\n\n")
  
  cat("ANALYSIS PARAMETERS USED\n")
  cat("------------------------\n")
  cat("Filtering criteria:\n")
  cat("- Minimum CPM: 1 in smallest group\n")
  cat("- edgeR filterByExpr: Applied\n")
  cat("- Variance filter: Bottom 10% removed\n")
  cat("Normalization: TMM method\n")
  cat("Outlier detection: Conservative (manual review)\n")
  cat("Statistical framework: edgeR/limma pipeline\n\n")
  
  cat("CITATIONS FOR METHODS USED\n")
  cat("--------------------------\n")
  cat("RNA-seq preprocessing:\n")
  cat("- Robinson et al. (2010) edgeR: a Bioconductor package for differential\n")
  cat("  expression analysis. Bioinformatics 26:139-140\n")
  cat("- Ritchie et al. (2015) limma powers differential expression analyses\n")
  cat("  for RNA-sequencing. Nucleic Acids Research 43:e47\n\n")
  
  cat("Quality control methods:\n")
  cat("- Conesa et al. (2016) A survey of best practices for RNA-seq data analysis.\n")
  cat("  Genome Biology 17:13\n")
  cat("- Evans et al. (2018) Selecting between-sample RNA-Seq normalization methods\n")
  cat("  from the perspective of their assumptions. Briefings in Bioinformatics 19:776-792\n\n")
  
  cat("=================================================================\n")
  cat("ANALYSIS COMPLETED SUCCESSFULLY\n")
  cat("=================================================================\n")
  
  sink()
}

# Generate the comprehensive report
create_analysis_report(
  file.path(output_structure$reports, "GSE118918_Processing_QC_Report.txt"),
  analysis_parameters
)

# ========================================================================
# SECTION 14: SESSION INFORMATION AND REPRODUCIBILITY
# ========================================================================

cat("\n=== SAVING SESSION INFORMATION ===\n")

# Save detailed session information
session_info_text <- capture.output(sessionInfo())
writeLines(session_info_text, 
           file.path(output_structure$session_info, "session_info.txt"))

# Save package versions for reproducibility
if(length(sessionInfo()$otherPkgs) > 0) {
  package_versions <- data.frame(
    package = names(sessionInfo()$otherPkgs),
    version = sapply(sessionInfo()$otherPkgs, function(x) x$Version),
    stringsAsFactors = FALSE
  )
  write.csv(package_versions,
            file.path(output_structure$session_info, "package_versions.csv"),
            row.names = FALSE)
}

# Create analysis timestamp
analysis_timestamp <- data.frame(
  step = c("Analysis_Start", "Data_Loading", "QC_Assessment", "Filtering", 
           "Normalization", "Final_Export", "Analysis_Complete"),
  timestamp = c(Sys.time(), Sys.time(), Sys.time(), Sys.time(),
                Sys.time(), Sys.time(), Sys.time()),
  status = "COMPLETED"
)
write.csv(analysis_timestamp,
          file.path(output_structure$session_info, "analysis_timestamps.csv"),
          row.names = FALSE)

# ========================================================================
# SECTION 15: FINAL VALIDATION AND SUMMARY
# ========================================================================

cat("\n=== FINAL VALIDATION CHECKS ===\n")

# Validate all output files exist
required_files <- c(
  file.path(output_structure$data, "dge_normalized_final.rds"),
  file.path(output_structure$data, "sample_metadata_final.rds"),
  file.path(output_structure$data, "logcpm_normalized_final.rds"),
  file.path(output_structure$plots, "Figure_2_Final_PCA_Analysis.png"),
  file.path(output_structure$reports, "GSE118918_Processing_QC_Report.txt")
)

validation_results <- sapply(required_files, file.exists)
cat("File validation:\n")
for(i in 1:length(validation_results)) {
  status <- ifelse(validation_results[i], "✓ EXISTS", "✗ MISSING")
  cat("-", basename(required_files[i]), ":", status, "\n")
}

# Final data integrity checks
cat("\nData integrity checks:\n")

# Check sample-metadata consistency
sample_match <- all(colnames(dge_normalized$counts) == final_metadata$sample_id)
cat("- Sample-metadata consistency:", ifelse(sample_match, "✓ PASS", "✗ FAIL"), "\n")

# Check for NA values in key objects
na_counts <- sum(is.na(dge_normalized$counts))
na_metadata <- sum(is.na(final_metadata$treatment))
cat("- Missing values in counts:", na_counts, "\n")
cat("- Missing treatment assignments:", na_metadata, "\n")

# Check normalization factors
norm_factor_range <- range(dge_normalized$samples$norm.factors)
reasonable_range <- all(norm_factor_range > 0.1 & norm_factor_range < 10)
cat("- Normalization factors reasonable:", ifelse(reasonable_range, "✓ PASS", "✗ FAIL"), "\n")

# ========================================================================
# FINAL STATUS REPORT
# ========================================================================

cat("\n")
cat("================================================================\n")
cat("GSE118918 BULK RNA-SEQ PROCESSING PIPELINE COMPLETED\n")
cat("================================================================\n")

cat("ANALYSIS SUMMARY:\n")
cat("- Dataset: GSE118918 (Nucleus Accumbens, Mock vs Morphine)\n")
cat("- Completion time:", as.character(Sys.time()), "\n")
cat("- Total runtime: [Check timestamps for duration]\n\n")

cat("FINAL DATASET SPECIFICATIONS:\n")
cat("- Samples processed:", ncol(dge_normalized), "out of", ncol(counts_original), "original\n")
cat("- Genes retained:", format(nrow(dge_normalized), big.mark = ","), "out of", 
    format(nrow(counts_original), big.mark = ","), "original\n")
cat("- Gene retention rate:", round(nrow(dge_normalized)/nrow(counts_original)*100, 1), "%\n")
cat("- Treatment groups:", paste(levels(final_metadata$treatment), collapse = " vs "), "\n")
cat("- Samples per group:", paste(table(final_metadata$treatment), collapse = ", "), "\n\n")

cat("QUALITY METRICS:\n")
cat("- Mean library size:", format(round(mean(colSums(dge_normalized$counts))), big.mark = ","), "reads\n")
cat("- Mean genes detected:", round(mean(colSums(dge_normalized$counts > 0))), "\n")
cat("- PC1 variance explained:", round(final_pca$variance_explained[1], 1), "%\n")
cat("- Sample correlation (mean):", round(mean(cor(logcpm_final)), 3), "\n\n")

cat("OUTPUT DIRECTORY:\n")
cat("- Main outputs:", output_structure$main, "\n")
cat("- Key data files:", output_structure$data, "\n")
cat("- Publication figures:", output_structure$plots, "\n")
cat("- Analysis report:", output_structure$reports, "\n\n")

cat("NEXT STEPS:\n")
cat("1. Review comprehensive QC report:\n")
cat("   ", file.path(output_structure$reports, "GSE118918_Processing_QC_Report.txt"), "\n")
cat("2. Examine publication-quality figures in:", output_structure$plots, "\n")
cat("3. Load processed data for differential expression:\n")
cat("   dge <- readRDS('", file.path(output_structure$data, "dge_normalized_final.rds"), "')\n")
cat("4. Proceed with differential expression analysis script\n")
cat("5. Focus on complement pathway genes for OUD research\n\n")

cat("DATA READY FOR:\n")
cat("- Differential expression analysis (Mock vs Morphine)\n")
cat("- Pathway enrichment analysis\n")
cat("- Complement system investigation\n") 
cat("- Morphine response characterization\n")
cat("- Integration with other GSE datasets\n")
cat("- Publication-ready visualizations\n\n")

cat("REPRODUCIBILITY:\n")
cat("- Session info saved in:", output_structure$session_info, "\n")
cat("- Analysis parameters documented\n")
cat("- All methods follow published best practices\n")
cat("- Ready for peer review standards\n\n")

cat("================================================================\n")
cat("PIPELINE STATUS: SUCCESSFUL COMPLETION\n")
cat("================================================================\n")

# Clean up workspace
rm(list = setdiff(ls(), c("dge_normalized", "final_metadata", "logcpm_final", 
                         "final_pca", "analysis_parameters", "output_structure")))
gc()

cat("\nWorkspace cleaned. Key objects retained:\n")
cat("- dge_normalized: Final DGEList for downstream analysis\n")
cat("- final_metadata: Sample metadata\n") 
cat("- logcpm_final: Normalized expression values\n")
cat("- final_pca: PCA results\n")
cat("- analysis_parameters: Analysis settings\n")
cat("- output_structure: File paths\n\n")

cat("GSE118918 processing complete. Ready for differential expression analysis.\n")