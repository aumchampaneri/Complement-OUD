# ========================================================================
# GSE239387 Nucleus Accumbens RNA-seq Processing and Quality Control Pipeline
# ========================================================================
# 
# STUDY OVERVIEW:
# Dataset: GSE239387 - Nucleus Accumbens RNA-seq (Control vs Morphine)
# Tissue: Nucleus Accumbens (NAcc) 
# Treatment: Control vs Morphine treatment
# Data Type: Pre-processed FPKM values with comprehensive annotations
# 
# DATA STRUCTURE NOTES:
# - Pre-calculated FPKM values (not raw counts)
# - Extensive functional annotations included
# - KEGG pathway, GO terms, and transcription factor data available
# - Different processing approach needed compared to raw count datasets
# 
# SCIENTIFIC RATIONALE:
# Since data is pre-processed as FPKM, we'll focus on:
# - Quality assessment of pre-processed data
# - Log transformation for downstream analysis
# - Leverage rich annotation data for pathway analysis
# - Prepare data for cross-dataset comparison with GSE118918
# ========================================================================

# Set reproducibility parameters
set.seed(42)
options(stringsAsFactors = FALSE)

# Record analysis start time
analysis_start_time <- Sys.time()
cat("GSE239387 Analysis started at:", as.character(analysis_start_time), "\n")

# ========================================================================
# SECTION 1: LIBRARY LOADING AND ENVIRONMENT SETUP
# ========================================================================

cat("=== LOADING REQUIRED LIBRARIES ===\n")

required_packages <- c("dplyr", "tidyr", "ggplot2", "pheatmap", "RColorBrewer", 
                      "corrplot", "VennDiagram", "ggrepel", "gridExtra", "tibble")

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    stop(paste("Required package", pkg, "is not installed."))
  }
}

# Additional packages for FPKM data analysis
optional_packages <- c("limma", "PCAtools", "ComplexHeatmap", "EnhancedVolcano")
loaded_optional <- character()

for(pkg in optional_packages) {
  if(require(pkg, character.only = TRUE, quietly = TRUE)) {
    loaded_optional <- c(loaded_optional, pkg)
  }
}

cat("Core packages loaded successfully\n")
cat("Optional packages available:", paste(loaded_optional, collapse = ", "), "\n")

# ========================================================================
# SECTION 2: DIRECTORY STRUCTURE AND DATA LOADING
# ========================================================================

cat("\n=== SETTING UP ANALYSIS ENVIRONMENT ===\n")

# Define base directory
base_dir <- "/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus/GSE239387"
setwd(base_dir)

# Create comprehensive output directory structure
output_dirs <- c(
  "Outputs/01_Processing_QC",
  "Outputs/01_Processing_QC/Data",
  "Outputs/01_Processing_QC/Figures", 
  "Outputs/01_Processing_QC/Reports",
  "Outputs/01_Processing_QC/Session_Info",
  "Outputs/01_Processing_QC/Tables"
)

for(dir in output_dirs) {
  if(!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat("Created directory:", dir, "\n")
  }
}

# Load the pre-processed data
cat("\n=== LOADING GSE239387 DATA ===\n")
data_file <- "Data/GSE239387_all.genes.expression.annot.txt"

if(!file.exists(data_file)) {
  stop("Data file not found: ", data_file)
}

# Read the comprehensive annotation file with proper handling
raw_data <- read.delim(data_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("Data loaded successfully\n")
cat("Dimensions:", nrow(raw_data), "genes x", ncol(raw_data), "columns\n")

# ========================================================================
# SECTION 3: DATA STRUCTURE ANALYSIS AND PREPARATION
# ========================================================================

cat("\n=== ANALYZING DATA STRUCTURE ===\n")

# Extract expression data (FPKM columns)
expression_cols <- grep("_fpkm$", colnames(raw_data), value = TRUE)
cat("Expression columns found:", length(expression_cols), "\n")
cat("Sample columns:", paste(expression_cols, collapse = ", "), "\n")

# Extract expression matrix and convert to numeric
expression_data <- raw_data[, expression_cols]
rownames(expression_data) <- raw_data$id

# Ensure all expression data is numeric
for(col in 1:ncol(expression_data)) {
  expression_data[, col] <- as.numeric(expression_data[, col])
}

# Handle any NA values
expression_data[is.na(expression_data)] <- 0

# Create sample metadata
sample_names <- gsub("_NAc_fpkm", "", expression_cols)
sample_metadata <- data.frame(
  sample_id = sample_names,
  treatment = ifelse(grepl("^Ctrl", sample_names), "Control", "Morphine"),
  replicate = as.numeric(gsub(".*([0-9])$", "\\1", sample_names)),
  stringsAsFactors = FALSE
)

rownames(sample_metadata) <- sample_names
colnames(expression_data) <- sample_names

cat("Sample metadata created:\n")
print(sample_metadata)

# Extract annotation data
annotation_data <- raw_data[, !colnames(raw_data) %in% expression_cols]
cat("Annotation columns:", ncol(annotation_data), "\n")

# ========================================================================
# SECTION 4: ENHANCED DATA QUALITY ASSESSMENT
# ========================================================================

cat("\n=== COMPREHENSIVE QUALITY ASSESSMENT ===\n")

# Basic statistics
cat("Expression data summary:\n")
cat("Min value:", min(expression_data, na.rm = TRUE), "\n")
cat("Max value:", max(expression_data, na.rm = TRUE), "\n") 
cat("Number of zeros:", sum(expression_data == 0, na.rm = TRUE), "\n")
cat("Percentage zeros:", round(sum(expression_data == 0, na.rm = TRUE) / length(as.matrix(expression_data)) * 100, 2), "%\n")

# Check for missing values
missing_data <- sum(is.na(expression_data))
cat("Missing values in expression data:", missing_data, "\n")

# Distribution analysis
cat("\n=== EXPRESSION DISTRIBUTION ANALYSIS ===\n")

# Create log2(FPKM + 1) transformed data for analysis
log_fpkm <- log2(expression_data + 1)

# Enhanced sample-level statistics - FIX: Use correct column name
sample_stats <- data.frame(
  sample = colnames(expression_data),
  treatment = sample_metadata$treatment,  # Fixed: was $Treatment
  total_fpkm = colSums(expression_data, na.rm = TRUE),
  mean_fpkm = colMeans(expression_data, na.rm = TRUE),
  median_fpkm = apply(expression_data, 2, median, na.rm = TRUE),
  detected_genes = colSums(expression_data > 0.1, na.rm = TRUE),
  high_expr_genes = colSums(expression_data > 10, na.rm = TRUE),
  very_high_expr_genes = colSums(expression_data > 100, na.rm = TRUE)
)

print(sample_stats)

# Comprehensive gene-level statistics
gene_stats <- data.frame(
  ensembl_id = rownames(expression_data),
  symbol = annotation_data$Symbol[match(rownames(expression_data), annotation_data$id)],
  mean_fpkm = rowMeans(expression_data, na.rm = TRUE),
  median_fpkm = apply(expression_data, 1, median, na.rm = TRUE),
  max_fpkm = apply(expression_data, 1, max, na.rm = TRUE),
  min_fpkm = apply(expression_data, 1, min, na.rm = TRUE),
  sd_fpkm = apply(expression_data, 1, sd, na.rm = TRUE),
  cv = apply(expression_data, 1, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)),
  detected_samples = rowSums(expression_data > 0.1, na.rm = TRUE),
  ctrl_mean = rowMeans(expression_data[, grep("Ctrl", colnames(expression_data))], na.rm = TRUE),
  mor_mean = rowMeans(expression_data[, grep("MOR", colnames(expression_data))], na.rm = TRUE),
  stringsAsFactors = FALSE
)

# Calculate fold change
gene_stats$log2fc_mor_vs_ctrl <- log2((gene_stats$mor_mean + 0.001) / (gene_stats$ctrl_mean + 0.001))

# Filter out genes with infinite CV values
gene_stats$cv[is.infinite(gene_stats$cv)] <- NA

cat("Gene detection summary:\n")
cat("Total genes:", nrow(gene_stats), "\n")
cat("Genes detected in all samples:", sum(gene_stats$detected_samples == ncol(expression_data)), "\n")
cat("Genes detected in >50% samples:", sum(gene_stats$detected_samples > ncol(expression_data)/2), "\n")
cat("Genes with valid symbols:", sum(!is.na(gene_stats$symbol) & gene_stats$symbol != ""), "\n")

# Enhanced filtering criteria
cat("\n=== APPLYING ENHANCED FILTERING CRITERIA ===\n")
pre_filter_genes <- nrow(gene_stats)

# Multi-criteria filtering
filtered_genes <- gene_stats %>%
  filter(
    mean_fpkm >= 0.5,           # Minimum expression
    detected_samples >= 2,       # Detected in at least 2 samples
    max_fpkm >= 1.0,            # At least 1 FPKM in one sample
    !is.na(cv) & !is.infinite(cv), # Valid CV - Fixed: added !is.na check
    cv < 10,                    # Not extremely variable
    !is.na(symbol) & symbol != "" # Has gene symbol
  )

cat("Genes before filtering:", pre_filter_genes, "\n")
cat("Genes after filtering:", nrow(filtered_genes), "\n")
cat("Retention rate:", round(nrow(filtered_genes)/pre_filter_genes*100, 2), "%\n")

# ========================================================================
# SECTION 5: ENHANCED COMPLEMENT GENE IDENTIFICATION
# ========================================================================

cat("\n=== ENHANCED COMPLEMENT PATHWAY ANALYSIS ===\n")

# Comprehensive complement gene patterns with better organization
complement_patterns <- list(
  classical = c("^C1Q[ABC]", "^C1R$", "^C1S$", "^C2$", "^C4[AB]$"),
  alternative = c("^CFB$", "^CFD$", "^CFP$", "^CFH$", "^CFI$", "PROPERDIN"),
  lectin = c("^MBL[12]", "^MASP[12]", "^FCN[123]", "^COLEC"),
  membrane_attack = c("^C[3-9]$", "^C5AR[12]", "^C3AR[12]"),
  regulators = c("^CD[45][569]$", "^C4BP[AB]$", "^SERPING1$", "^DAF$", "^MCP$"),
  receptors = c("^CR[1-4]$", "^C1QR", "^ITGAM$", "^ITGAX$", "^ITGB2$"),
  associated = c("CLUSTERIN", "^CLU$", "VITRONECTIN", "^VTN$", "^APOE$", "^APOH$")
)

# Check available column names first
available_go_cols <- colnames(annotation_data)[grep("GO", colnames(annotation_data), ignore.case = TRUE)]
cat("Available GO columns:", paste(available_go_cols, collapse = ", "), "\n")

# Search in multiple annotation fields with improved logic and proper column name handling
complement_genes <- annotation_data %>%
  filter(
    !is.na(Symbol) & Symbol != "" & (
      # Pattern matching in symbol
      grepl(paste(unlist(complement_patterns), collapse = "|"), Symbol, ignore.case = TRUE) |
      # Description search
      grepl("complement", Description, ignore.case = TRUE) |
      # Pathway search (including KEGG complement pathway)
      grepl("complement|04610", Pathway, ignore.case = TRUE) |
      # GO term search - use proper column references with periods
      {if("GO.Component" %in% colnames(annotation_data)) 
        grepl("complement", annotation_data[["GO.Component"]], ignore.case = TRUE) else FALSE} |
      {if("GO.Function" %in% colnames(annotation_data)) 
        grepl("complement", annotation_data[["GO.Function"]], ignore.case = TRUE) else FALSE} |
      {if("GO.Process" %in% colnames(annotation_data)) 
        grepl("complement", annotation_data[["GO.Process"]], ignore.case = TRUE) else FALSE}
    )
  ) %>%
  distinct(Symbol, .keep_all = TRUE) %>%
  arrange(Symbol)

cat("Found", nrow(complement_genes), "potential complement genes\n")

if(nrow(complement_genes) > 0) {
  cat("Complement genes detected:\n")
  
  # Enhanced complement gene analysis with expression data and error handling
  valid_comp_ids <- complement_genes$id[complement_genes$id %in% rownames(expression_data)]
  
  if(length(valid_comp_ids) > 0) {
    complement_with_expr <- complement_genes %>%
      filter(id %in% valid_comp_ids) %>%
      select(id, Symbol, Description) %>%
      mutate(
        mean_ctrl_fpkm = sapply(id, function(x) {
          if(x %in% rownames(expression_data)) {
            ctrl_cols <- grepl("Ctrl", colnames(expression_data))
            mean(as.numeric(expression_data[x, ctrl_cols]), na.rm = TRUE)
          } else { NA }
        }),
        mean_mor_fpkm = sapply(id, function(x) {
          if(x %in% rownames(expression_data)) {
            mor_cols <- grepl("MOR", colnames(expression_data))
            mean(as.numeric(expression_data[x, mor_cols]), na.rm = TRUE)
          } else { NA }
        }),
        detected_samples = sapply(id, function(x) {
          if(x %in% rownames(expression_data)) {
            sum(as.numeric(expression_data[x, ]) > 0.1, na.rm = TRUE)
          } else { 0 }
        })
      ) %>%
      filter(!is.na(mean_ctrl_fpkm) & !is.na(mean_mor_fpkm)) %>%
      mutate(
        fold_change = (mean_mor_fpkm + 0.001) / (mean_ctrl_fpkm + 0.001),
        log2_fc = log2(fold_change),
        # Classify complement pathway
        pathway_class = case_when(
          grepl(paste(complement_patterns$classical, collapse = "|"), Symbol, ignore.case = TRUE) ~ "Classical",
          grepl(paste(complement_patterns$alternative, collapse = "|"), Symbol, ignore.case = TRUE) ~ "Alternative",
          grepl(paste(complement_patterns$lectin, collapse = "|"), Symbol, ignore.case = TRUE) ~ "Lectin",
          grepl(paste(complement_patterns$membrane_attack, collapse = "|"), Symbol, ignore.case = TRUE) ~ "Membrane Attack",
          grepl(paste(complement_patterns$regulators, collapse = "|"), Symbol, ignore.case = TRUE) ~ "Regulator",
          grepl(paste(complement_patterns$receptors, collapse = "|"), Symbol, ignore.case = TRUE) ~ "Receptor",
          TRUE ~ "Associated"
        ),
        # Statistical significance (simple t-test)
        p_value = sapply(id, function(x) {
          if(x %in% rownames(expression_data)) {
            ctrl_cols <- grepl("Ctrl", colnames(expression_data))
            mor_cols <- grepl("MOR", colnames(expression_data))
            ctrl_vals <- as.numeric(expression_data[x, ctrl_cols])
            mor_vals <- as.numeric(expression_data[x, mor_cols])
            if(all(ctrl_vals == 0) & all(mor_vals == 0)) return(1)
            tryCatch(t.test(mor_vals, ctrl_vals)$p.value, error = function(e) 1)
          } else { 1 }
        })
      ) %>%
      arrange(desc(mean_ctrl_fpkm))
    
    print(complement_with_expr)
    
    # Pathway summary
    pathway_summary <- complement_with_expr %>%
      group_by(pathway_class) %>%
      summarise(
        count = n(),
        mean_expression = round(mean(mean_ctrl_fpkm + mean_mor_fpkm)/2, 2),
        upregulated = sum(log2_fc > 0.5),
        downregulated = sum(log2_fc < -0.5),
        .groups = 'drop'
      )
    
    cat("\nComplement pathway summary:\n")
    print(pathway_summary)
    
    # Save complement gene analysis
    write.csv(complement_with_expr, "Outputs/01_Processing_QC/Data/complement_genes_with_expression.csv", row.names = FALSE)
    write.csv(pathway_summary, "Outputs/01_Processing_QC/Data/complement_pathway_summary.csv", row.names = FALSE)
  }
} else {
  cat("No complement genes found with current search criteria\n")
}

# ========================================================================
# SECTION 6: COMPREHENSIVE VISUALIZATION AND QC PLOTS
# ========================================================================

cat("\n=== GENERATING COMPREHENSIVE QC VISUALIZATIONS ===\n")

# Filter genes with zero variance for PCA
gene_variances <- apply(log_fpkm, 1, var, na.rm = TRUE)
variable_genes <- log_fpkm[gene_variances > 0 & !is.na(gene_variances), ]

cat("Genes retained for PCA (non-zero variance):", nrow(variable_genes), "out of", nrow(log_fpkm), "\n")

# Function to safely create plots
safe_plot <- function(plot_func, filename_base, ...) {
  tryCatch({
    # PDF version
    pdf(paste0("Outputs/01_Processing_QC/Figures/", filename_base, ".pdf"), ...)
    result <- plot_func()
    dev.off()
    
    # PNG version
    png(paste0("Outputs/01_Processing_QC/Figures/", filename_base, ".png"), 
        width = 800, height = 600, res = 100)
    plot_func()
    dev.off()
    
    cat("Successfully created:", filename_base, "\n")
  }, error = function(e) {
    cat("Error creating", filename_base, ":", e$message, "\n")
    dev.off()  # Ensure device is closed
  })
}

# Enhanced sample correlation heatmap
cor_matrix <- cor(log_fpkm, use = "complete.obs")

safe_plot(function() {
  pheatmap::pheatmap(cor_matrix,
           main = "Sample Correlation (log2 FPKM)",
           annotation_col = sample_metadata[, "treatment", drop = FALSE],
           show_rownames = TRUE,
           show_colnames = TRUE,
           color = colorRampPalette(c("blue", "white", "red"))(100))
}, "01_Sample_Correlation_Heatmap", width = 10, height = 8)

# PCA analysis with filtered genes
pca_result <- prcomp(t(variable_genes), center = TRUE, scale. = TRUE)
pca_data <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2], 
  sample = rownames(pca_result$x)
)
pca_data$treatment <- sample_metadata$treatment[match(pca_data$sample, rownames(sample_metadata))]

# Enhanced PCA plot - PDF version
pdf("Outputs/01_Processing_QC/Figures/02_PCA_Plot.pdf", width = 8, height = 6)
p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = treatment)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = sample)) +
  labs(title = "PCA of GSE239387 Samples",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2,1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2,2] * 100, 1), "%)")) +
  theme_bw()
print(p_pca)
dev.off()

# Enhanced PCA plot - PNG version
png("Outputs/01_Processing_QC/Figures/02_PCA_Plot.png", width = 800, height = 600, res = 100)
print(p_pca)
dev.off()

# Enhanced expression distribution plots with better density visualization
safe_plot(function() {
  par(mfrow = c(2, 2))
  
  # Overall distribution
  boxplot(log_fpkm, main = "Log2(FPKM + 1) Distributions", 
          xlab = "Samples", ylab = "Log2(FPKM + 1)", las = 2,
          col = ifelse(grepl("Ctrl", colnames(log_fpkm)), "lightblue", "lightcoral"))
  
  # Treatment comparison
  ctrl_data <- log_fpkm[, grepl("Ctrl", colnames(log_fpkm))]
  mor_data <- log_fpkm[, grepl("MOR", colnames(log_fpkm))]
  
  hist(rowMeans(ctrl_data), breaks = 50, col = rgb(0,0,1,0.5), 
       main = "Mean Expression Distribution", xlab = "Log2(FPKM + 1)", 
       xlim = range(c(rowMeans(ctrl_data), rowMeans(mor_data))))
  hist(rowMeans(mor_data), breaks = 50, col = rgb(1,0,0,0.5), add = TRUE)
  legend("topright", c("Control", "Morphine"), fill = c(rgb(0,0,1,0.5), rgb(1,0,0,0.5)))
  
  # Detection rates
  detection_ctrl <- rowSums(ctrl_data > 0.1) / ncol(ctrl_data)
  detection_mor <- rowSums(mor_data > 0.1) / ncol(mor_data)
  plot(detection_ctrl, detection_mor, pch = 16, col = rgb(0,0,0,0.6),
       xlab = "Control Detection Rate", ylab = "Morphine Detection Rate",
       main = "Gene Detection Comparison")
  abline(0, 1, col = "red", lty = 2)
  
  # Sample-wise statistics
  sample_medians <- apply(log_fpkm, 2, median, na.rm = TRUE)
  barplot(sample_medians, las = 2, main = "Median Expression per Sample",
          col = ifelse(grepl("Ctrl", names(sample_medians)), "lightblue", "lightcoral"))
  
}, "03_Expression_Distributions", width = 12, height = 10)

# Additional QC plot - sample distances - PDF version
pdf("Outputs/01_Processing_QC/Figures/04_Sample_Distances.pdf", width = 10, height = 8)
dist_matrix <- dist(t(variable_genes))
hc <- hclust(dist_matrix)
plot(hc, main = "Sample Clustering (Euclidean Distance)", 
     xlab = "Samples", ylab = "Distance")
dev.off()

# Additional QC plot - sample distances - PNG version
png("Outputs/01_Processing_QC/Figures/04_Sample_Distances.png", width = 1000, height = 800, res = 100)
plot(hc, main = "Sample Clustering (Euclidean Distance)", 
     xlab = "Samples", ylab = "Distance")
dev.off()

# Complement gene analysis and visualization
if(nrow(complement_genes) > 0) {
  cat("Analyzing complement genes...\n")
  
  # Get complement gene expression data
  comp_expr_data <- log_fpkm[complement_genes$id[complement_genes$id %in% rownames(log_fpkm)], ]
  
  if(nrow(comp_expr_data) > 0) {
    # Complement gene heatmap - PDF version
    pdf("Outputs/01_Processing_QC/Figures/05_Complement_Genes_Expression.pdf", width = 12, height = 8)
    
    # Prepare data for heatmap
    comp_matrix <- as.matrix(comp_expr_data)
    rownames(comp_matrix) <- complement_genes$Symbol[match(rownames(comp_matrix), complement_genes$id)]
    
    # Remove rows with all zeros
    comp_matrix <- comp_matrix[rowSums(comp_matrix) > 0, ]
    
    if(nrow(comp_matrix) > 1) {
      # Create annotation for samples
      col_annotation <- data.frame(
        Treatment = ifelse(grepl("Ctrl", colnames(comp_matrix)), "Control", "Morphine")
      )
      rownames(col_annotation) <- colnames(comp_matrix)
      
      # Create heatmap
      pheatmap::pheatmap(comp_matrix,  # Use explicit namespace
               main = "Complement Genes Expression Heatmap",
               annotation_col = col_annotation,
               scale = "row",
               clustering_distance_rows = "correlation",
               clustering_distance_cols = "correlation",
               fontsize_row = 8,
               fontsize_col = 8)
    }
    
    dev.off()
    
    # Complement gene heatmap - PNG version
    if(nrow(comp_matrix) > 1) {
      png("Outputs/01_Processing_QC/Figures/05_Complement_Genes_Expression.png", width = 1200, height = 800, res = 100)
      pheatmap::pheatmap(comp_matrix,
               main = "Complement Genes Expression Heatmap",
               annotation_col = col_annotation,
               scale = "row",
               clustering_distance_rows = "correlation",
               clustering_distance_cols = "correlation",
               fontsize_row = 8,
               fontsize_col = 8)
      dev.off()
    }
    
    # Complement gene statistics
    comp_stats <- filtered_genes[filtered_genes$ensembl_id %in% complement_genes$id, ]
    comp_stats <- comp_stats[order(comp_stats$mean_fpkm, decreasing = TRUE), ]
    
    cat("\nTop complement genes by expression:\n")
    print(head(comp_stats[, c("symbol", "mean_fpkm", "log2fc_mor_vs_ctrl", "detected_samples")], 10))
  }
}

# Enhanced QC visualizations - PDF version
pdf("Outputs/01_Processing_QC/Figures/06_Enhanced_QC_Plots.pdf", width = 14, height = 10)

# Multi-panel plot
par(mfrow = c(2, 3))

# 1. Expression distribution by treatment
ctrl_means <- rowMeans(expression_data[, grepl("Ctrl", colnames(expression_data))])
mor_means <- rowMeans(expression_data[, grepl("MOR", colnames(expression_data))])

plot(log2(ctrl_means + 1), log2(mor_means + 1), 
     xlab = "Control log2(FPKM + 1)", ylab = "Morphine log2(FPKM + 1)",
     main = "Expression Correlation: Control vs Morphine", pch = 16, col = rgb(0, 0, 0, 0.6))
abline(0, 1, col = "red", lty = 2)

# 2. CV distribution
hist(gene_stats$cv[gene_stats$cv < 5 & !is.na(gene_stats$cv) & !is.infinite(gene_stats$cv)], 
     breaks = 50, main = "Coefficient of Variation Distribution", 
     xlab = "CV", col = "lightblue")

# 3. Detection rate by expression level
detection_bins <- cut(gene_stats$mean_fpkm, breaks = c(0, 0.1, 1, 10, 100, Inf))
detection_by_expr <- tapply(gene_stats$detected_samples, detection_bins, mean, na.rm = TRUE)
barplot(detection_by_expr, main = "Detection Rate by Expression Level",
        ylab = "Mean Samples Detected", xlab = "FPKM Bins")

# 4. Sample-wise detection
barplot(sample_stats$detected_genes, names.arg = sample_stats$sample,
        main = "Genes Detected per Sample", ylab = "Number of Genes",
        las = 2, col = ifelse(grepl("Ctrl", sample_stats$sample), "blue", "red"))

# 5. Variance vs Mean relationship
plot(log2(gene_stats$mean_fpkm + 1), log2(gene_variances + 1),
     xlab = "Mean log2(FPKM + 1)", ylab = "Variance log2(FPKM + 1)",
     main = "Mean-Variance Relationship", pch = 16, col = rgb(0, 0, 0, 0.6))

# 6. Top variable genes
top_var_genes <- head(gene_stats[order(gene_stats$cv, decreasing = TRUE, na.last = TRUE), ], 20)
top_var_genes <- top_var_genes[!is.na(top_var_genes$cv), ]  # Remove NAs
barplot(top_var_genes$cv[1:min(15, nrow(top_var_genes))], 
        names.arg = top_var_genes$symbol[1:min(15, nrow(top_var_genes))],
        main = "Top Variable Genes", ylab = "CV", las = 2)

dev.off()

# Enhanced QC visualizations - PNG version
png("Outputs/01_Processing_QC/Figures/06_Enhanced_QC_Plots.png", width = 1400, height = 1000, res = 100)
# Multi-panel plot
par(mfrow = c(2, 3))

# Recreate all 6 plots for PNG
plot(log2(ctrl_means + 1), log2(mor_means + 1), 
     xlab = "Control log2(FPKM + 1)", ylab = "Morphine log2(FPKM + 1)",
     main = "Expression Correlation: Control vs Morphine", pch = 16, col = rgb(0, 0, 0, 0.6))
abline(0, 1, col = "red", lty = 2)

hist(gene_stats$cv[gene_stats$cv < 5 & !is.na(gene_stats$cv) & !is.infinite(gene_stats$cv)], 
     breaks = 50, main = "Coefficient of Variation Distribution", 
     xlab = "CV", col = "lightblue")

barplot(detection_by_expr, main = "Detection Rate by Expression Level",
        ylab = "Mean Samples Detected", xlab = "FPKM Bins")

barplot(sample_stats$detected_genes, names.arg = sample_stats$sample,
        main = "Genes Detected per Sample", ylab = "Number of Genes",
        las = 2, col = ifelse(grepl("Ctrl", sample_stats$sample), "blue", "red"))

plot(log2(gene_stats$mean_fpkm + 1), log2(gene_variances + 1),
     xlab = "Mean log2(FPKM + 1)", ylab = "Variance log2(FPKM + 1)",
     main = "Mean-Variance Relationship", pch = 16, col = rgb(0, 0, 0, 0.6))

barplot(top_var_genes$cv[1:min(15, nrow(top_var_genes))], 
        names.arg = top_var_genes$symbol[1:min(15, nrow(top_var_genes))],
        main = "Top Variable Genes", ylab = "CV", las = 2)

dev.off()

# Complement gene analysis and visualization
if(nrow(complement_genes) > 0) {
  cat("Analyzing complement genes...\n")
  
  # Get complement gene expression data
  comp_expr_data <- log_fpkm[complement_genes$id[complement_genes$id %in% rownames(log_fpkm)], ]
  
  if(nrow(comp_expr_data) > 0) {
    # Complement gene heatmap - PDF version
    pdf("Outputs/01_Processing_QC/Figures/05_Complement_Genes_Expression.pdf", width = 12, height = 8)
    
    # Prepare data for heatmap
    comp_matrix <- as.matrix(comp_expr_data)
    rownames(comp_matrix) <- complement_genes$Symbol[match(rownames(comp_matrix), complement_genes$id)]
    
    # Remove rows with all zeros
    comp_matrix <- comp_matrix[rowSums(comp_matrix) > 0, ]
    
    if(nrow(comp_matrix) > 1) {
      # Create annotation for samples
      col_annotation <- data.frame(
        Treatment = ifelse(grepl("Ctrl", colnames(comp_matrix)), "Control", "Morphine")
      )
      rownames(col_annotation) <- colnames(comp_matrix)
      
      # Create heatmap
      pheatmap::pheatmap(comp_matrix,  # Use explicit namespace
               main = "Complement Genes Expression Heatmap",
               annotation_col = col_annotation,
               scale = "row",
               clustering_distance_rows = "correlation",
               clustering_distance_cols = "correlation",
               fontsize_row = 8,
               fontsize_col = 8)
    }
    
    dev.off()
    
    # Complement gene heatmap - PNG version
    if(nrow(comp_matrix) > 1) {
      png("Outputs/01_Processing_QC/Figures/05_Complement_Genes_Expression.png", width = 1200, height = 800, res = 100)
      pheatmap::pheatmap(comp_matrix,
               main = "Complement Genes Expression Heatmap",
               annotation_col = col_annotation,
               scale = "row",
               clustering_distance_rows = "correlation",
               clustering_distance_cols = "correlation",
               fontsize_row = 8,
               fontsize_col = 8)
      dev.off()
    }
    
    # Complement gene statistics
    comp_stats <- filtered_genes[filtered_genes$ensembl_id %in% complement_genes$id, ]
    comp_stats <- comp_stats[order(comp_stats$mean_fpkm, decreasing = TRUE), ]
    
    cat("\nTop complement genes by expression:\n")
    print(head(comp_stats[, c("symbol", "mean_fpkm", "log2fc_mor_vs_ctrl", "detected_samples")], 10))
  }
}

# Enhanced QC visualizations - PDF version
pdf("Outputs/01_Processing_QC/Figures/06_Enhanced_QC_Plots.pdf", width = 14, height = 10)

# Multi-panel plot
par(mfrow = c(2, 3))

# 1. Expression distribution by treatment
ctrl_means <- rowMeans(expression_data[, grepl("Ctrl", colnames(expression_data))])
mor_means <- rowMeans(expression_data[, grepl("MOR", colnames(expression_data))])

plot(log2(ctrl_means + 1), log2(mor_means + 1), 
     xlab = "Control log2(FPKM + 1)", ylab = "Morphine log2(FPKM + 1)",
     main = "Expression Correlation: Control vs Morphine", pch = 16, col = rgb(0, 0, 0, 0.6))
abline(0, 1, col = "red", lty = 2)

# 2. CV distribution
hist(gene_stats$cv[gene_stats$cv < 5 & !is.na(gene_stats$cv) & !is.infinite(gene_stats$cv)], 
     breaks = 50, main = "Coefficient of Variation Distribution", 
     xlab = "CV", col = "lightblue")

# 3. Detection rate by expression level
detection_bins <- cut(gene_stats$mean_fpkm, breaks = c(0, 0.1, 1, 10, 100, Inf))
detection_by_expr <- tapply(gene_stats$detected_samples, detection_bins, mean, na.rm = TRUE)
barplot(detection_by_expr, main = "Detection Rate by Expression Level",
        ylab = "Mean Samples Detected", xlab = "FPKM Bins")

# 4. Sample-wise detection
barplot(sample_stats$detected_genes, names.arg = sample_stats$sample,
        main = "Genes Detected per Sample", ylab = "Number of Genes",
        las = 2, col = ifelse(grepl("Ctrl", sample_stats$sample), "blue", "red"))

# 5. Variance vs Mean relationship
plot(log2(gene_stats$mean_fpkm + 1), log2(gene_variances + 1),
     xlab = "Mean log2(FPKM + 1)", ylab = "Variance log2(FPKM + 1)",
     main = "Mean-Variance Relationship", pch = 16, col = rgb(0, 0, 0, 0.6))

# 6. Top variable genes
top_var_genes <- head(gene_stats[order(gene_stats$cv, decreasing = TRUE, na.last = TRUE), ], 20)
top_var_genes <- top_var_genes[!is.na(top_var_genes$cv), ]  # Remove NAs
barplot(top_var_genes$cv[1:min(15, nrow(top_var_genes))], 
        names.arg = top_var_genes$symbol[1:min(15, nrow(top_var_genes))],
        main = "Top Variable Genes", ylab = "CV", las = 2)

dev.off()

# Enhanced QC visualizations - PNG version
png("Outputs/01_Processing_QC/Figures/06_Enhanced_QC_Plots.png", width = 1400, height = 1000, res = 100)
# Multi-panel plot
par(mfrow = c(2, 3))

# Recreate all 6 plots for PNG
plot(log2(ctrl_means + 1), log2(mor_means + 1), 
     xlab = "Control log2(FPKM + 1)", ylab = "Morphine log2(FPKM + 1)",
     main = "Expression Correlation: Control vs Morphine", pch = 16, col = rgb(0, 0, 0, 0.6))
abline(0, 1, col = "red", lty = 2)

hist(gene_stats$cv[gene_stats$cv < 5 & !is.na(gene_stats$cv) & !is.infinite(gene_stats$cv)], 
     breaks = 50, main = "Coefficient of Variation Distribution", 
     xlab = "CV", col = "lightblue")

barplot(detection_by_expr, main = "Detection Rate by Expression Level",
        ylab = "Mean Samples Detected", xlab = "FPKM Bins")

barplot(sample_stats$detected_genes, names.arg = sample_stats$sample,
        main = "Genes Detected per Sample", ylab = "Number of Genes",
        las = 2, col = ifelse(grepl("Ctrl", sample_stats$sample), "blue", "red"))

plot(log2(gene_stats$mean_fpkm + 1), log2(gene_variances + 1),
     xlab = "Mean log2(FPKM + 1)", ylab = "Variance log2(FPKM + 1)",
     main = "Mean-Variance Relationship", pch = 16, col = rgb(0, 0, 0, 0.6))

barplot(top_var_genes$cv[1:min(15, nrow(top_var_genes))], 
        names.arg = top_var_genes$symbol[1:min(15, nrow(top_var_genes))],
        main = "Top Variable Genes", ylab = "CV", las = 2)

dev.off()

# Complement gene analysis and visualization
if(nrow(complement_genes) > 0) {
  cat("Analyzing complement genes...\n")
  
  # Get complement gene expression data
  comp_expr_data <- log_fpkm[complement_genes$id[complement_genes$id %in% rownames(log_fpkm)], ]
  
  if(nrow(comp_expr_data) > 0) {
    # Complement gene heatmap - PDF version
    pdf("Outputs/01_Processing_QC/Figures/05_Complement_Genes_Expression.pdf", width = 12, height = 8)
    
    # Prepare data for heatmap
    comp_matrix <- as.matrix(comp_expr_data)
    rownames(comp_matrix) <- complement_genes$Symbol[match(rownames(comp_matrix), complement_genes$id)]
    
    # Remove rows with all zeros
    comp_matrix <- comp_matrix[rowSums(comp_matrix) > 0, ]
    
    if(nrow(comp_matrix) > 1) {
      # Create annotation for samples
      col_annotation <- data.frame(
        Treatment = ifelse(grepl("Ctrl", colnames(comp_matrix)), "Control", "Morphine")
      )
      rownames(col_annotation) <- colnames(comp_matrix)
      
      # Create heatmap
      pheatmap::pheatmap(comp_matrix,  # Use explicit namespace
               main = "Complement Genes Expression Heatmap",
               annotation_col = col_annotation,
               scale = "row",
               clustering_distance_rows = "correlation",
               clustering_distance_cols = "correlation",
               fontsize_row = 8,
               fontsize_col = 8)
    }
    
    dev.off()
    
    # Complement gene heatmap - PNG version
    if(nrow(comp_matrix) > 1) {
      png("Outputs/01_Processing_QC/Figures/05_Complement_Genes_Expression.png", width = 1200, height = 800, res = 100)
      pheatmap::pheatmap(comp_matrix,
               main = "Complement Genes Expression Heatmap",
               annotation_col = col_annotation,
               scale = "row",
               clustering_distance_rows = "correlation",
               clustering_distance_cols = "correlation",
               fontsize_row = 8,
               fontsize_col = 8)
      dev.off()
    }
    
    # Complement gene statistics
    comp_stats <- filtered_genes[filtered_genes$ensembl_id %in% complement_genes$id, ]
    comp_stats <- comp_stats[order(comp_stats$mean_fpkm, decreasing = TRUE), ]
    
    cat("\nTop complement genes by expression:\n")
    print(head(comp_stats[, c("symbol", "mean_fpkm", "log2fc_mor_vs_ctrl", "detected_samples")], 10))
  }
}

# Enhanced QC visualizations - PDF version
pdf("Outputs/01_Processing_QC/Figures/06_Enhanced_QC_Plots.pdf", width = 14, height = 10)

# Multi-panel plot
par(mfrow = c(2, 3))

# 1. Expression distribution by treatment
ctrl_means <- rowMeans(expression_data[, grepl("Ctrl", colnames(expression_data))])
mor_means <- rowMeans(expression_data[, grepl("MOR", colnames(expression_data))])

plot(log2(ctrl_means + 1), log2(mor_means + 1), 
     xlab = "Control log2(FPKM + 1)", ylab = "Morphine log2(FPKM + 1)",
     main = "Expression Correlation: Control vs Morphine", pch = 16, col = rgb(0, 0, 0, 0.6))
abline(0, 1, col = "red", lty = 2)

# 2. CV distribution
hist(gene_stats$cv[gene_stats$cv < 5 & !is.na(gene_stats$cv) & !is.infinite(gene_stats$cv)], 
     breaks = 50, main = "Coefficient of Variation Distribution", 
     xlab = "CV", col = "lightblue")

# 3. Detection rate by expression level
detection_bins <- cut(gene_stats$mean_fpkm, breaks = c(0, 0.1, 1, 10, 100, Inf))
detection_by_expr <- tapply(gene_stats$detected_samples, detection_bins, mean, na.rm = TRUE)
barplot(detection_by_expr, main = "Detection Rate by Expression Level",
        ylab = "Mean Samples Detected", xlab = "FPKM Bins")

# 4. Sample-wise detection
barplot(sample_stats$detected_genes, names.arg = sample_stats$sample,
        main = "Genes Detected per Sample", ylab = "Number of Genes",
        las = 2, col = ifelse(grepl("Ctrl", sample_stats$sample), "blue", "red"))

# 5. Variance vs Mean relationship
plot(log2(gene_stats$mean_fpkm + 1), log2(gene_variances + 1),
     xlab = "Mean log2(FPKM + 1)", ylab = "Variance log2(FPKM + 1)",
     main = "Mean-Variance Relationship", pch = 16, col = rgb(0, 0, 0, 0.6))

# 6. Top variable genes
top_var_genes <- head(gene_stats[order(gene_stats$cv, decreasing = TRUE, na.last = TRUE), ], 20)
top_var_genes <- top_var_genes[!is.na(top_var_genes$cv), ]  # Remove NAs
barplot(top_var_genes$cv[1:min(15, nrow(top_var_genes))], 
        names.arg = top_var_genes$symbol[1:min(15, nrow(top_var_genes))],
        main = "Top Variable Genes", ylab = "CV", las = 2)

dev.off()

# Enhanced QC visualizations - PNG version
png("Outputs/01_Processing_QC/Figures/06_Enhanced_QC_Plots.png", width = 1400, height = 1000, res = 100)
# Multi-panel plot
par(mfrow = c(2, 3))

# Recreate all 6 plots for PNG
plot(log2(ctrl_means + 1), log2(mor_means + 1), 
     xlab = "Control log2(FPKM + 1)", ylab = "Morphine log2(FPKM + 1)",
     main = "Expression Correlation: Control vs Morphine", pch = 16, col = rgb(0, 0, 0, 0.6))
abline(0, 1, col = "red", lty = 2)

hist(gene_stats$cv[gene_stats$cv < 5 & !is.na(gene_stats$cv) & !is.infinite(gene_stats$cv)], 
     breaks = 50, main = "Coefficient of Variation Distribution", 
     xlab = "CV", col = "lightblue")

barplot(detection_by_expr, main = "Detection Rate by Expression Level",
        ylab = "Mean Samples Detected", xlab = "FPKM Bins")

barplot(sample_stats$detected_genes, names.arg = sample_stats$sample,
        main = "Genes Detected per Sample", ylab = "Number of Genes",
        las = 2, col = ifelse(grepl("Ctrl", sample_stats$sample), "blue", "red"))

plot(log2(gene_stats$mean_fpkm + 1), log2(gene_variances + 1),
     xlab = "Mean log2(FPKM + 1)", ylab = "Variance log2(FPKM + 1)",
     main = "Mean-Variance Relationship", pch = 16, col = rgb(0, 0, 0, 0.6))

barplot(top_var_genes$cv[1:min(15, nrow(top_var_genes))], 
        names.arg = top_var_genes$symbol[1:min(15, nrow(top_var_genes))],
        main = "Top Variable Genes", ylab = "CV", las = 2)

dev.off()

# Complement gene analysis and visualization
if(nrow(complement_genes) > 0) {
  cat("Analyzing complement genes...\n")
  
  # Get complement gene expression data
  comp_expr_data <- log_fpkm[complement_genes$id[complement_genes$id %in% rownames(log_fpkm)], ]
  
  if(nrow(comp_expr_data) > 0) {
    # Complement gene heatmap - PDF version
    pdf("Outputs/01_Processing_QC/Figures/05_Complement_Genes_Expression.pdf", width = 12, height = 8)
    
    # Prepare data for heatmap
    comp_matrix <- as.matrix(comp_expr_data)
    rownames(comp_matrix) <- complement_genes$Symbol[match(rownames(comp_matrix), complement_genes$id)]
    
    # Remove rows with all zeros
    comp_matrix <- comp_matrix[rowSums(comp_matrix) > 0, ]
    
    if(nrow(comp_matrix) > 1) {
      # Create annotation for samples
      col_annotation <- data.frame(
        Treatment = ifelse(grepl("Ctrl", colnames(comp_matrix)), "Control", "Morphine")
      )
      rownames(col_annotation) <- colnames(comp_matrix)
      
      # Create heatmap
      pheatmap::pheatmap(comp_matrix,  # Use explicit namespace
               main = "Complement Genes Expression Heatmap",
               annotation_col = col_annotation,
               scale = "row",
               clustering_distance_rows = "correlation",
               clustering_distance_cols = "correlation",
               fontsize_row = 8,
               fontsize_col = 8)
    }
    
    dev.off()
    
    # Complement gene heatmap - PNG version
    if(nrow(comp_matrix) > 1) {
      png("Outputs/01_Processing_QC/Figures/05_Complement_Genes_Expression.png", width = 1200, height = 800, res = 100)
      pheatmap::pheatmap(comp_matrix,
               main = "Complement Genes Expression Heatmap",
               annotation_col = col_annotation,
               scale = "row",
               clustering_distance_rows = "correlation",
               clustering_distance_cols = "correlation",
               fontsize_row = 8,
               fontsize_col = 8)
      dev.off()
    }
    
    # Complement gene statistics
    comp_stats <- filtered_genes[filtered_genes$ensembl_id %in% complement_genes$id, ]
    comp_stats <- comp_stats[order(comp_stats$mean_fpkm, decreasing = TRUE), ]
    
    cat("\nTop complement genes by expression:\n")
    print(head(comp_stats[, c("symbol", "mean_fpkm", "log2fc_mor_vs_ctrl", "detected_samples")], 10))
  }
}

# Enhanced QC visualizations - PDF version
pdf("Outputs/01_Processing_QC/Figures/06_Enhanced_QC_Plots.pdf", width = 14, height = 10)

# Multi-panel plot
par(mfrow = c(2, 3))

# 1. Expression distribution by treatment
ctrl_means <- rowMeans(expression_data[, grepl("Ctrl", colnames(expression_data))])
mor_means <- rowMeans(expression_data[, grepl("MOR", colnames(expression_data))])

plot(log2(ctrl_means + 1), log2(mor_means + 1), 
     xlab = "Control log2(FPKM + 1)", ylab = "Morphine log2(FPKM + 1)",
     main = "Expression Correlation: Control vs Morphine", pch = 16, col = rgb(0, 0, 0, 0.6))
abline(0, 1, col = "red", lty = 2)

# 2. CV distribution
hist(gene_stats$cv[gene_stats$cv < 5 & !is.na(gene_stats$cv) & !is.infinite(gene_stats$cv)], 
     breaks = 50, main = "Coefficient of Variation Distribution", 
     xlab = "CV", col = "lightblue")

# 3. Detection rate by expression level
detection_bins <- cut(gene_stats$mean_fpkm, breaks = c(0, 0.1, 1, 10, 100, Inf))
detection_by_expr <- tapply(gene_stats$detected_samples, detection_bins, mean, na.rm = TRUE)
barplot(detection_by_expr, main = "Detection Rate by Expression Level",
        ylab = "Mean Samples Detected", xlab = "FPKM Bins")

# 4. Sample-wise detection
barplot(sample_stats$detected_genes, names.arg = sample_stats$sample,
        main = "Genes Detected per Sample", ylab = "Number of Genes",
        las = 2, col = ifelse(grepl("Ctrl", sample_stats$sample), "blue", "red"))

# 5. Variance vs Mean relationship
plot(log2(gene_stats$mean_fpkm + 1), log2(gene_variances + 1),
     xlab = "Mean log2(FPKM + 1)", ylab = "Variance log2(FPKM + 1)",
     main = "Mean-Variance Relationship", pch = 16, col = rgb(0, 0, 0, 0.6))

# 6. Top variable genes
top_var_genes <- head(gene_stats[order(gene_stats$cv, decreasing = TRUE, na.last = TRUE), ], 20)
top_var_genes <- top_var_genes[!is.na(top_var_genes$cv), ]  # Remove NAs
barplot(top_var_genes$cv[1:min(15, nrow(top_var_genes))], 
        names.arg = top_var_genes$symbol[1:min(15, nrow(top_var_genes))],
        main = "Top Variable Genes", ylab = "CV", las = 2)

dev.off()

# Enhanced QC visualizations - PNG version
png("Outputs/01_Processing_QC/Figures/06_Enhanced_QC_Plots.png", width = 1400, height = 1000, res = 100)
# Multi-panel plot
par(mfrow = c(2, 3))

# Recreate all 6 plots for PNG
plot(log2(ctrl_means + 1), log2(mor_means + 1), 
     xlab = "Control log2(FPKM + 1)", ylab = "Morphine log2(FPKM + 1)",
     main = "Expression Correlation: Control vs Morphine", pch = 16, col = rgb(0, 0, 0, 0.6))
abline(0, 1, col = "red", lty = 2)

hist(gene_stats$cv[gene_stats$cv < 5 & !is.na(gene_stats$cv) & !is.infinite(gene_stats$cv)], 
     breaks = 50, main = "Coefficient of Variation Distribution", 
     xlab = "CV", col = "lightblue")

barplot(detection_by_expr, main = "Detection Rate by Expression Level",
        ylab = "Mean Samples Detected", xlab = "FPKM Bins")

barplot(sample_stats$detected_genes, names.arg = sample_stats$sample,
        main = "Genes Detected per Sample", ylab = "Number of Genes",
        las = 2, col = ifelse(grepl("Ctrl", sample_stats$sample), "blue", "red"))

plot(log2(gene_stats$mean_fpkm + 1), log2(gene_variances + 1),
     xlab = "Mean log2(FPKM + 1)", ylab = "Variance log2(FPKM + 1)",
     main = "Mean-Variance Relationship", pch = 16, col = rgb(0, 0, 0, 0.6))

barplot(top_var_genes$cv[1:min(15, nrow(top_var_genes))], 
        names.arg = top_var_genes$symbol[1:min(15, nrow(top_var_genes))],
        main = "Top Variable Genes", ylab = "CV", las = 2)

dev.off()

# Complement gene analysis and visualization
if(nrow(complement_genes) > 0) {
  cat("Analyzing complement genes...\n")
  
  # Get complement gene expression data
  comp_expr_data <- log_fpkm[complement_genes$id[complement_genes$id %in% rownames(log_fpkm)], ]
  
  if(nrow(comp_expr_data) > 0) {
    # Complement gene heatmap - PDF version
    pdf("Outputs/01_Processing_QC/Figures/05_Complement_Genes_Expression.pdf", width = 12, height = 8)
    
    # Prepare data for heatmap
    comp_matrix <- as.matrix(comp_expr_data)
    rownames(comp_matrix) <- complement_genes$Symbol[match(rownames(comp_matrix), complement_genes$id)]
    
    # Remove rows with all zeros
    comp_matrix <- comp_matrix[rowSums(comp_matrix) > 0, ]
    
    if(nrow(comp_matrix) > 1) {
      # Create annotation for samples
      col_annotation <- data.frame(
        Treatment = ifelse(grepl("Ctrl", colnames(comp_matrix)), "Control", "Morphine")
      )
      rownames(col_annotation) <- colnames(comp_matrix)
      
      # Create heatmap
      pheatmap::pheatmap(comp_matrix,  # Use explicit namespace
               main = "Complement Genes Expression Heatmap",
               annotation_col = col_annotation,
               scale = "row",
               clustering_distance_rows = "correlation",
               clustering_distance_cols = "correlation",
               fontsize_row = 8,
               fontsize_col = 8)
    }
    
    dev.off()
    
    # Complement gene heatmap - PNG version
    if(nrow(comp_matrix) > 1) {
      png("Outputs/01_Processing_QC/Figures/05_Complement_Genes_Expression.png", width = 1200, height = 800, res = 100)
      pheatmap::pheatmap(comp_matrix,
               main = "Complement Genes Expression Heatmap",
               annotation_col = col_annotation,
               scale = "row",
               clustering_distance_rows = "correlation",
               clustering_distance_cols = "correlation",
               fontsize_row = 8,
               fontsize_col = 8)
      dev.off()
    }
    
    # Complement gene statistics
    comp_stats <- filtered_genes[filtered_genes$ensembl_id %in% complement_genes$id, ]
    comp_stats <- comp_stats[order(comp_stats$mean_fpkm, decreasing = TRUE), ]
    
    cat("\nTop complement genes by expression:\n")
    print(head(comp_stats[, c("symbol", "mean_fpkm", "log2fc_mor_vs_ctrl", "detected_samples")], 10))
  }
}

# Enhanced QC visualizations - PDF version
pdf("Outputs/01_Processing_QC/Figures/06_Enhanced_QC_Plots.pdf", width = 14, height = 10)

# Multi-panel plot
par(mfrow = c(2, 3))

# 1. Expression distribution by treatment
ctrl_means <- rowMeans(expression_data[, grepl("Ctrl", colnames(expression_data))])
mor_means <- rowMeans(expression_data[, grepl("MOR", colnames(expression_data))])

plot(log2(ctrl_means + 1), log2(mor_means + 1), 
     xlab = "Control log2(FPKM + 1)", ylab = "Morphine log2(FPKM + 1)",
     main = "Expression Correlation: Control vs Morphine", pch = 16, col = rgb(0, 0, 0, 0.6))
abline(0, 1, col = "red", lty = 2)

# 2. CV distribution
hist(gene_stats$cv[gene_stats$cv < 5 & !is.na(gene_stats$cv) & !is.infinite(gene_stats$cv)], 
     breaks = 50, main = "Coefficient of Variation Distribution", 
     xlab = "CV", col = "lightblue")

# 3. Detection rate by expression level
detection_bins <- cut(gene_stats$mean_fpkm, breaks = c(0, 0.1, 1, 10, 100, Inf))
detection_by_expr <- tapply(gene_stats$detected_samples, detection_bins, mean, na.rm = TRUE)
barplot(detection_by_expr, main = "Detection Rate by Expression Level",
        ylab = "Mean Samples Detected", xlab = "FPKM Bins")

# 4. Sample-wise detection
barplot(sample_stats$detected_genes, names.arg = sample_stats$sample,
        main = "Genes Detected per Sample", ylab = "Number of Genes",
        las = 2, col = ifelse(grepl("Ctrl", sample_stats$sample), "blue", "red"))

# 5. Variance vs Mean relationship
plot(log2(gene_stats$mean_fpkm + 1), log2(gene_variances + 1),
     xlab = "Mean log2(FPKM + 1)", ylab = "Variance log2(FPKM + 1)",
     main = "Mean-Variance Relationship", pch = 16, col = rgb(0, 0, 0, 0.6))

# 6. Top variable genes
top_var_genes <- head(gene_stats[order(gene_stats$cv, decreasing = TRUE, na.last = TRUE), ], 20)
top_var_genes <- top_var_genes[!is.na(top_var_genes$cv), ]  # Remove NAs
barplot(top_var_genes$cv[1:min(15, nrow(top_var_genes))], 
        names.arg = top_var_genes$symbol[1:min(15, nrow(top_var_genes))],
        main = "Top Variable Genes", ylab = "CV", las = 2)

dev.off()

# Enhanced QC visualizations - PNG version
png("Outputs/01_Processing_QC/Figures/06_Enhanced_QC_Plots.png", width = 1400, height = 1000, res = 100)
# Multi-panel plot
par(mfrow = c(2, 3))

# Recreate all 6 plots for PNG
plot(log2(ctrl_means + 1), log2(mor_means + 1), 
     xlab = "Control log2(FPKM + 1)", ylab = "Morphine log2(FPKM + 1)",
     main = "Expression Correlation: Control vs Morphine", pch = 16, col = rgb(0, 0, 0, 0.6))
abline(0, 1, col = "red", lty = 2)

hist(gene_stats$cv[gene_stats$cv < 5 & !is.na(gene_stats$cv) & !is.infinite(gene_stats$cv)], 
     breaks = 50, main = "Coefficient of Variation Distribution", 
     xlab = "CV", col = "lightblue")

barplot(detection_by_expr, main = "Detection Rate by Expression Level",
        ylab = "Mean Samples Detected", xlab = "FPKM Bins")

barplot(sample_stats$detected_genes, names.arg = sample_stats$sample,
        main = "Genes Detected per Sample", ylab = "Number of Genes",
        las = 2, col = ifelse(grepl("Ctrl", sample_stats$sample), "blue", "red"))

plot(log2(gene_stats$mean_fpkm + 1), log2(gene_variances + 1),
     xlab = "Mean log2(FPKM + 1)", ylab = "Variance log2(FPKM + 1)",
     main = "Mean-Variance Relationship", pch = 16, col = rgb(0, 0, 0, 0.6))

barplot(top_var_genes$cv[1:min(15, nrow(top_var_genes))], 
        names.arg = top_var_genes$symbol[1:min(15, nrow(top_var_genes))],
        main = "Top Variable Genes", ylab = "CV", las = 2)

dev.off()

# Complement gene analysis and visualization
if(nrow(complement_genes) > 0) {
  cat("Analyzing complement genes...\n")
  
  # Get complement gene expression data
  comp_expr_data <- log_fpkm[complement_genes$id[complement_genes$id %in% rownames(log_fpkm)], ]
  
  if(nrow(comp_expr_data) > 0) {
    # Complement gene heatmap - PDF version
    pdf("Outputs/01_Processing_QC/Figures/05_Complement_Genes_Expression.pdf", width = 12, height = 8)
    
    # Prepare data for heatmap
    comp_matrix <- as.matrix(comp_expr_data)
    rownames(comp_matrix) <- complement_genes$Symbol[match(rownames(comp_matrix), complement_genes$id)]
    
    # Remove rows with all zeros
    comp_matrix <- comp_matrix[rowSums(comp_matrix) > 0, ]
    
    if(nrow(comp_matrix) > 1) {
      # Create annotation for samples
      col_annotation <- data.frame(
        Treatment = ifelse(grepl("Ctrl", colnames(comp_matrix)), "Control", "Morphine")
      )
      rownames(col_annotation) <- colnames(comp_matrix)
      
      # Create heatmap
      pheatmap::pheatmap(comp_matrix,  # Use explicit namespace
               main = "Complement Genes Expression Heatmap",
               annotation_col = col_annotation,
               scale = "row",
               clustering_distance_rows = "correlation",
               clustering_distance_cols = "correlation",
               fontsize_row = 8,
               fontsize_col = 8)
    }
    
    dev.off()
    
    # Complement gene heatmap - PNG version
    if(nrow(comp_matrix) > 1) {
      png("Outputs/01_Processing_QC/Figures/05_Complement_Genes_Expression.png", width = 1200, height = 800, res = 100)
      pheatmap::pheatmap(comp_matrix,
               main = "Complement Genes Expression Heatmap",
               annotation_col = col_annotation,
               scale = "row",
               clustering_distance_rows = "correlation",
               clustering_distance_cols = "correlation",
               fontsize_row = 8,
               fontsize_col = 8)
      dev.off()
    }
    
    # Complement gene statistics
    comp_stats <- filtered_genes[filtered_genes$ensembl_id %in% complement_genes$id, ]
    comp_stats <- comp_stats[order(comp_stats$mean_fpkm, decreasing = TRUE), ]
    
    cat("\nTop complement genes by expression:\n")
    print(head(comp_stats[, c("symbol", "mean_fpkm", "log2fc_mor_vs_ctrl", "detected_samples")], 10))
  }
}

# Enhanced QC visualizations - PDF version
pdf("Outputs/01_Processing_QC/Figures/06_Enhanced_QC_Plots.pdf", width = 14, height = 10)

# Multi-panel plot
par(mfrow = c(2, 3))

# 1. Expression distribution by treatment
ctrl_means <- rowMeans(expression_data[, grepl("Ctrl", colnames(expression_data))])
mor_means <- rowMeans(expression_data[, grepl("MOR", colnames(expression_data))])

plot(log2(ctrl_means + 1), log2(mor_means + 1), 
     xlab = "Control log2(FPKM + 1)", ylab = "Morphine log2(FPKM + 1)",
     main = "Expression Correlation: Control vs Morphine", pch = 16, col = rgb(0, 0, 0, 0.6))
abline(0, 1, col = "red", lty = 2)

# 2. CV distribution
hist(gene_stats$cv[gene_stats$cv < 5 & !is.na(gene_stats$cv) & !is.infinite(gene_stats$cv)], 
     breaks = 50, main = "Coefficient of Variation Distribution", 
     xlab = "CV", col = "lightblue")

# 3. Detection rate by expression level
detection_bins <- cut(gene_stats$mean_fpkm, breaks = c(0, 0.1, 1, 10, 100, Inf))
detection_by_expr <- tapply(gene_stats$detected_samples, detection_bins, mean, na.rm = TRUE)
barplot(detection_by_expr, main = "Detection Rate by Expression Level",
        ylab = "Mean Samples Detected", xlab = "FPKM Bins")

# 4. Sample-wise detection
barplot(sample_stats$detected_genes, names.arg = sample_stats$sample,
        main = "Genes Detected per Sample", ylab = "Number of Genes",
        las = 2, col = ifelse(grepl("Ctrl", sample_stats$sample), "blue", "red"))

# 5. Variance vs Mean relationship
plot(log2(gene_stats$mean_fpkm + 1), log2(gene_variances + 1),
     xlab = "Mean log2(FPKM + 1)", ylab = "Variance log2(FPKM + 1)",
     main = "Mean-Variance Relationship", pch = 16, col = rgb(0, 0, 0, 0.6))

# 6. Top variable genes
top_var_genes <- head(gene_stats[order(gene_stats$cv, decreasing = TRUE, na.last = TRUE), ], 20)
top_var_genes <- top_var_genes[!is.na(top_var_genes$cv), ]  # Remove NAs
barplot(top_var_genes$cv[1:min(15, nrow(top_var_genes))], 
        names.arg = top_var_genes$symbol[1:min(15, nrow(top_var_genes))],
        main = "Top Variable Genes", ylab = "CV", las = 2)

dev.off()

# Enhanced QC visualizations - PNG version
png("Outputs/01_Processing_QC/Figures/06_Enhanced_QC_Plots.png", width = 1400, height = 1000, res = 100)
# Multi-panel plot
par(mfrow = c(2, 3))

# Recreate all 6 plots for PNG
plot(log2(ctrl_means + 1), log2(mor_means + 1), 
     xlab = "Control log2(FPKM + 1)", ylab = "Morphine log2(FPKM + 1)",
     main = "Expression Correlation: Control vs Morphine", pch = 16, col = rgb(0, 0, 0, 0.6))
abline(0, 1, col = "red", lty = 2)

hist(gene_stats$cv[gene_stats$cv < 5 & !is.na(gene_stats$cv) & !is.infinite(gene_stats$cv)], 
     breaks = 50, main = "Coefficient of Variation Distribution", 
     xlab = "CV", col = "lightblue")

barplot(detection_by_expr, main = "Detection Rate by Expression Level",
        ylab = "Mean Samples Detected", xlab = "FPKM Bins")

barplot(sample_stats$detected_genes, names.arg = sample_stats$sample,
        main = "Genes Detected per Sample", ylab = "Number of Genes",
        las = 2, col = ifelse(grepl("Ctrl", sample_stats$sample), "blue", "red"))

plot(log2(gene_stats$mean_fpkm + 1), log2(gene_variances + 1),
     xlab = "Mean log2(FPKM + 1)", ylab = "Variance log2(FPKM + 1)",
     main = "Mean-Variance Relationship", pch = 16, col = rgb(0, 0, 0, 0.6))

barplot(top_var_genes$cv[1:min(15, nrow(top_var_genes))], 
        names.arg = top_var_genes$symbol[1:min(15, nrow(top_var_genes))],
        main = "Top Variable Genes", ylab = "CV", las = 2)

dev.off()

# Complement gene analysis and visualization
if(nrow(complement_genes) > 0) {
  cat("Analyzing complement genes...\n")
  
  # Get complement gene expression data
  comp_expr_data <- log_fpkm[complement_genes$id[complement_genes$id %in% rownames(log_fpkm)], ]
  
  if(nrow(comp_expr_data) > 0) {
    # Complement gene heatmap - PDF version
    pdf("Outputs/01_Processing_QC/Figures/05_Complement_Genes_Expression.pdf", width = 12, height = 8)
    
    # Prepare data for heatmap
    comp_matrix <- as.matrix(comp_expr_data)
    rownames(comp_matrix) <- complement_genes$Symbol[match(rownames(comp_matrix), complement_genes$id)]
    
    # Remove rows with all zeros
    comp_matrix <- comp_matrix[rowSums(comp_matrix) > 0, ]
    
    if(nrow(comp_matrix) > 1) {
      # Create annotation for samples
      col_annotation <- data.frame(
        Treatment = ifelse(grepl("Ctrl", colnames(comp_matrix)), "Control", "Morphine")
      )
      rownames(col_annotation) <- colnames(comp_matrix)
      
      # Create heatmap
      pheatmap::pheatmap(comp_matrix,  # Use explicit namespace
               main = "Complement Genes Expression Heatmap",
               annotation_col = col_annotation,
               scale = "row",
               clustering_distance_rows = "correlation",
               clustering_distance_cols = "correlation",
               fontsize_row = 8,
               fontsize_col = 8)
    }
    
    dev.off()
    
    # Complement gene heatmap - PNG version
    if(nrow(comp_matrix) > 1) {
      png("Outputs/01_Processing_QC/Figures/05_Complement_Genes_Expression.png", width = 1200, height = 800, res = 100)
      pheatmap::pheatmap(comp_matrix,
               main = "Complement Genes Expression Heatmap",
               annotation_col = col_annotation,
               scale = "row",
               clustering_distance_rows = "correlation",
               clustering_distance_cols = "correlation",
               fontsize_row = 8,
               fontsize_col = 8)
      dev.off()
    }
    
    # Complement gene statistics
    comp_stats <- filtered_genes[filtered_genes$ensembl_id %in% complement_genes$id, ]
    comp_stats <- comp_stats[order(comp_stats$mean_fpkm, decreasing = TRUE), ]
    
    cat("\nTop complement genes by expression:\n")
    print(head(comp_stats[, c("symbol", "mean_fpkm", "log2fc_mor_vs_ctrl", "detected_samples")], 10))
  }
}

# Enhanced QC visualizations - PDF version
pdf("Outputs/01_Processing_QC/Figures/06_Enhanced_QC_Plots.pdf", width = 14, height = 10)

# Multi-panel plot
par(mfrow = c(2, 3))

# 1. Expression distribution by treatment
ctrl_means <- rowMeans(expression_data[, grepl("Ctrl", colnames(expression_data))])
mor_means <- rowMeans(expression_data[, grepl("MOR", colnames(expression_data))])

plot(log2(ctrl_means + 1), log2(mor_means + 1), 
     xlab = "Control log2(FPKM + 1)", ylab = "Morphine log2(FPKM + 1)",
     main = "Expression Correlation: Control vs Morphine", pch = 16, col = rgb(0, 0, 0, 0.6))
abline(0, 1, col = "red", lty = 2)

# 2. CV distribution
hist(gene_stats$cv[gene_stats$cv < 5 & !is.na(gene_stats$cv) & !is.infinite(gene_stats$cv)], 
     breaks = 50, main = "Coefficient of Variation Distribution", 
     xlab = "CV", col = "lightblue")

# 3. Detection rate by expression level
detection_bins <- cut(gene_stats$mean_fpkm, breaks = c(0, 0.1, 1, 10, 100, Inf))
detection_by_expr <- tapply(gene_stats$detected_samples, detection_bins, mean, na.rm = TRUE)
barplot(detection_by_expr, main = "Detection Rate by Expression Level",
        ylab = "Mean Samples Detected", xlab = "FPKM Bins")

# 4. Sample-wise detection
barplot(sample_stats$detected_genes, names.arg = sample_stats$sample,
        main = "Genes Detected per Sample", ylab = "Number of Genes",
        las = 2, col = ifelse(grepl("Ctrl", sample_stats$sample), "blue", "red"))

# 5. Variance vs Mean relationship
plot(log2(gene_stats$mean_fpkm + 1), log2(gene_variances + 1),
     xlab = "Mean log2(FPKM + 1)", ylab = "Variance log2(FPKM + 1)",
     main = "Mean-Variance Relationship", pch = 16, col = rgb(0, 0, 0, 0.6))

# 6. Top variable genes
top_var_genes <- head(gene_stats[order(gene_stats$cv, decreasing = TRUE, na.last = TRUE), ], 20)
top_var_genes <- top_var_genes[!is.na(top_var_genes$cv), ]  # Remove NAs
barplot(top_var_genes$cv[1:min(15, nrow(top_var_genes))], 
        names.arg = top_var_genes$symbol[1:min(15, nrow(top_var_genes))],
        main = "Top Variable Genes", ylab = "CV", las = 2)

dev.off()

# Enhanced QC visualizations - PNG version
png("Outputs/01_Processing_QC/Figures/06_Enhanced_QC_Plots.png", width = 1400, height = 1000, res = 100)
# Multi-panel plot
par(mfrow = c(2, 3))

# Recreate all 6 plots for PNG
plot(log2(ctrl_means + 1), log2(mor_means + 1), 
     xlab = "Control log2(FPKM + 1)", ylab = "Morphine log2(FPKM + 1)",
     main = "Expression Correlation: Control vs Morphine", pch = 16, col = rgb(0, 0, 0, 0.6))
abline(0, 1, col = "red", lty = 2)

hist(gene_stats$cv[gene_stats$cv < 5 & !is.na(gene_stats$cv) & !is.infinite(gene_stats$cv)], 
     breaks = 50, main = "Coefficient of Variation Distribution", 
     xlab = "CV", col = "lightblue")

barplot(detection_by_expr, main = "Detection Rate by Expression Level",
        ylab = "Mean Samples Detected", xlab = "FPKM Bins")

barplot(sample_stats$detected_genes, names.arg = sample_stats$sample,
        main = "Genes Detected per Sample", ylab = "Number of Genes",
        las = 2, col = ifelse(grepl("Ctrl", sample_stats$sample), "blue", "red"))

plot(log2(gene_stats$mean_fpkm + 1), log2(gene_variances + 1),
     xlab = "Mean log2(FPKM + 1)", ylab = "Variance log2(FPKM + 1)",
     main = "Mean-Variance Relationship", pch = 16, col = rgb(0, 0, 0, 0.6))

barplot(top_var_genes$cv[1:min(15, nrow(top_var_genes))], 
        names.arg = top_var_genes$symbol[1:min(15, nrow(top_var_genes))],
        main = "Top Variable Genes", ylab = "CV", las = 2)

dev.off()

# Complement gene analysis and visualization
if(nrow(complement_genes) > 0) {
  cat("Analyzing complement genes...\n")
  
  # Get complement gene expression data
  comp_expr_data <- log_fpkm[complement_genes$id[complement_genes$id %in% rownames(log_fpkm)], ]
  
  if(nrow(comp_expr_data) > 0) {
    # Complement gene heatmap - PDF version
    pdf("Outputs/01_Processing_QC/Figures/05_Complement_Genes_Expression.pdf", width = 12, height = 8)
    
    # Prepare data for heatmap
    comp_matrix <- as.matrix(comp_expr_data)
    rownames(comp_matrix) <- complement_genes$Symbol[match(rownames(comp_matrix), complement_genes$id)]
    
    # Remove rows with all zeros
    comp_matrix <- comp_matrix[rowSums(comp_matrix) > 0, ]
    
    if(nrow(comp_matrix) > 1) {
      # Create annotation for samples
      col_annotation <- data.frame(
        Treatment = ifelse(grepl("Ctrl", colnames(comp_matrix)), "Control", "Morphine")
      )
      rownames(col_annotation) <- colnames(comp_matrix)
      
      # Create heatmap
      pheatmap::pheatmap(comp_matrix,  # Use explicit namespace
               main = "Complement Genes Expression Heatmap",
               annotation_col = col_annotation,
               scale = "row",
               clustering_distance_rows = "correlation",
               clustering_distance_cols = "correlation",
               fontsize_row = 8,
               fontsize_col = 8)
    }
    
    dev.off()
    
    # Complement gene heatmap - PNG version
    if(nrow(comp_matrix) > 1) {
      png("Outputs/01_Processing_QC/Figures/05_Complement_Genes_Expression.png", width = 1200, height = 800, res = 100)
      pheatmap::pheatmap(comp_matrix,
               main = "Complement Genes Expression Heatmap",
               annotation_col = col_annotation,
               scale = "row",
               clustering_distance_rows = "correlation",
               clustering_distance_cols = "correlation",
               fontsize_row = 8,
               fontsize_col = 8)
      dev.off()
    }
    
    # Complement gene statistics
    comp_stats <- filtered_genes[filtered_genes$ensembl_id %in% complement_genes$id, ]
    comp_stats <- comp_stats[order(comp_stats$mean_fpkm, decreasing = TRUE), ]
    
    cat("\nTop complement genes by expression:\n")
    print(head(comp_stats[, c("symbol", "mean_fpkm", "log2fc_mor_vs_ctrl", "detected_samples")], 10))
  }
}

# ========================================================================
# SECTION 7: ENHANCED ANNOTATION ANALYSIS
# ========================================================================

cat("\n=== ENHANCED ANNOTATION ANALYSIS ===\n")

# Fix annotation column names - need to check actual column names
available_cols <- colnames(annotation_data)
cat("Available annotation columns:", paste(available_cols, collapse = ", "), "\n")

# Build annotation summary dynamically based on available columns
annotation_types <- character()
genes_annotated <- numeric()
percentages <- numeric()

# Check Symbol
if("Symbol" %in% available_cols) {
  annotation_types <- c(annotation_types, "Symbol")
  genes_annotated <- c(genes_annotated, sum(!is.na(annotation_data$Symbol) & annotation_data$Symbol != ""))
  percentages <- c(percentages, sum(!is.na(annotation_data$Symbol) & annotation_data$Symbol != "") / nrow(annotation_data) * 100)
}

# Check Description
if("Description" %in% available_cols) {
  annotation_types <- c(annotation_types, "Description")
  genes_annotated <- c(genes_annotated, sum(!is.na(annotation_data$Description) & annotation_data$Description != ""))
  percentages <- c(percentages, sum(!is.na(annotation_data$Description) & annotation_data$Description != "") / nrow(annotation_data) * 100)
}

# Check other columns
for(col_check in c("KEGG_A_class", "KEGG_B_class", "Pathway", "GO.Component", "GO Component", 
                   "GO.Function", "GO Function", "GO.Process", "GO Process", "TF_family")) {
  if(col_check %in% available_cols) {
    annotation_types <- c(annotation_types, col_check)
    empty_vals <- c("-", "--", "")
    genes_annotated <- c(genes_annotated, sum(!is.na(annotation_data[[col_check]]) & 
                                              !annotation_data[[col_check]] %in% empty_vals))
    percentages <- c(percentages, sum(!is.na(annotation_data[[col_check]]) & 
                                      !annotation_data[[col_check]] %in% empty_vals) / nrow(annotation_data) * 100)
  }
}

annotation_summary <- data.frame(
  annotation_type = annotation_types,
  genes_annotated = genes_annotated,
  percentage = round(percentages, 2)
)

print(annotation_summary)
write.csv(annotation_summary, "Outputs/01_Processing_QC/Tables/annotation_summary.csv", row.names = FALSE)

# ========================================================================
# SECTION 8: SAVE PROCESSED DATA
# ========================================================================

cat("\n=== SAVING PROCESSED DATA ===\n")

write.csv(filtered_genes, "Outputs/01_Processing_QC/Data/filtered_gene_stats.csv", row.names = FALSE)
write.csv(sample_stats, "Outputs/01_Processing_QC/Data/sample_statistics.csv", row.names = FALSE)
write.csv(expression_data, "Outputs/01_Processing_QC/Data/expression_matrix_fpkm.csv")
write.csv(log_fpkm, "Outputs/01_Processing_QC/Data/expression_matrix_log2fpkm.csv")

# Create filtered expression matrix for downstream analysis
filtered_log_fpkm <- log_fpkm[filtered_genes$ensembl_id, ]
write.csv(filtered_log_fpkm, "Outputs/01_Processing_QC/Data/filtered_expression_matrix_log2fpkm.csv")
saveRDS(filtered_log_fpkm, "Outputs/01_Processing_QC/Data/filtered_log2fpkm.rds")

if(nrow(complement_genes) > 0) {
  write.csv(complement_genes, "Outputs/01_Processing_QC/Data/complement_genes_detected.csv", row.names = FALSE)
  if(exists("comp_stats")) {
    write.csv(comp_stats, "Outputs/01_Processing_QC/Data/complement_gene_statistics.csv", row.names = FALSE)
  }
}

# ========================================================================
# SECTION 9: GENERATE COMPREHENSIVE PROCESSING REPORT
# ========================================================================

cat("\n=== GENERATING COMPREHENSIVE PROCESSING REPORT ===\n")

# Enhanced report with additional metrics
report_text <- paste0(
  "=================================================================\n",
  "GSE239387 RNA-SEQ PROCESSING AND QC ANALYSIS REPORT\n",
  "=================================================================\n\n",
  "STUDY INFORMATION\n",
  "-----------------\n",
  "Dataset ID: GSE239387\n",
  "Organism: Mus musculus\n",
  "Tissue: Nucleus Accumbens (NAcc)\n",
  "Experimental Design: Control vs Morphine treatment\n",
  "Data Type: Pre-processed FPKM values with annotations\n",
  "Analysis Date: ", as.character(Sys.time()), "\n",
  "R Version: ", R.version.string, "\n\n",
  
  "DATASET OVERVIEW\n",
  "----------------\n",
  "Total genes: ", nrow(expression_data), "\n",
  "Total samples: ", ncol(expression_data), "\n", 
  "Control samples: ", sum(sample_metadata$treatment == "Control"), "\n",
  "Morphine samples: ", sum(sample_metadata$treatment == "Morphine"), "\n",
  "Genes with symbols: ", sum(!is.na(annotation_data$Symbol) & annotation_data$Symbol != ""), "\n",
  "Complement genes detected: ", nrow(complement_genes), "\n",
  "Expressed genes (filtered): ", nrow(filtered_genes), "\n\n",
  
  "DATA QUALITY METRICS\n",
  "--------------------\n",
  "Mean FPKM across samples: ", round(mean(sample_stats$mean_fpkm), 2), "\n",
  "Median FPKM across samples: ", round(median(sample_stats$median_fpkm), 2), "\n",
  "Mean genes detected per sample: ", round(mean(sample_stats$detected_genes)), "\n",
  "Range of genes detected: ", min(sample_stats$detected_genes), " - ", max(sample_stats$detected_genes), "\n",
  "Percentage of zero values: ", round(sum(expression_data == 0) / length(as.matrix(expression_data)) * 100, 2), "%\n",
  "Variable genes for PCA: ", nrow(variable_genes), " out of ", nrow(log_fpkm), "\n",
  "Gene retention after filtering: ", round(nrow(filtered_genes) / nrow(gene_stats) * 100, 2), "%\n",
  "PC1 variance explained: ", round(summary(pca_result)$importance[2,1] * 100, 1), "%\n",
  "PC2 variance explained: ", round(summary(pca_result)$importance[2,2] * 100, 1), "%\n\n",
  
  "COMPLEMENT SYSTEM ANALYSIS\n",
  "--------------------------\n",
  "Complement genes found: ", nrow(complement_genes), "\n",
  if(nrow(complement_genes) > 0) {
    paste0("Mean complement expression (Control): ", 
           round(mean(complement_with_expr$mean_ctrl_fpkm), 2), "\n",
           "Mean complement expression (Morphine): ", 
           round(mean(complement_with_expr$mean_mor_fpkm), 2), "\n")
  } else {
    "No complement genes detected with current criteria\n"
  },
  "\n",
  
  "READY FOR DOWNSTREAM ANALYSIS\n",
  "------------------------------\n",
  "- Differential expression analysis (Control vs Morphine)\n",
  "- Pathway enrichment analysis using extensive annotations\n", 
  "- Complement system investigation\n",
  "- Cross-dataset comparison with GSE118918\n",
  "- Publication-ready visualizations\n\n",
  
  "OUTPUT FILES GENERATED\n",
  "----------------------\n",
  "- Expression matrices (FPKM and log2-transformed)\n",
  "- Filtered expression data for DE analysis\n",
  "- Sample and gene statistics\n",
  "- Annotation summaries\n",
  "- Complement gene analysis\n",
  "- Quality control visualizations\n",
  "- PCA results and plots\n\n",
  
  "=================================================================\n",
  "ANALYSIS COMPLETED SUCCESSFULLY\n",
  "================================================================="
)

writeLines(report_text, "Outputs/01_Processing_QC/Reports/GSE239387_Processing_QC_Report.txt")

# Save session info
writeLines(capture.output(sessionInfo()), 
          "Outputs/01_Processing_QC/Session_Info/session_info.txt")

# Summary output - Fix string concatenation
cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("ENHANCED GSE239387 PROCESSING SUMMARY\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("Total genes in dataset:", nrow(raw_data), "\n")
cat("Genes after filtering:", nrow(filtered_genes), "\n")
cat("Retention rate:", round(nrow(filtered_genes)/nrow(raw_data)*100, 2), "%\n")
cat("Complement genes detected:", nrow(complement_genes), "\n")
cat("Complement genes in filtered set:", 
    sum(complement_genes$id %in% filtered_genes$ensembl_id), "\n")
cat("Samples analyzed:", ncol(expression_data), "\n")
cat("Control samples:", sum(sample_stats$treatment == "Control"), "\n")
cat("Morphine samples:", sum(sample_stats$treatment == "Morphine"), "\n")
cat("\nProcessing completed successfully!\n")
cat("Results saved in: Outputs/01_Processing_QC/\n")
cat("GSE239387 processing complete. Ready for differential expression analysis.\n")
