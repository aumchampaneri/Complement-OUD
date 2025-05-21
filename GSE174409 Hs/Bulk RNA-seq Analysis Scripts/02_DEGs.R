library(limma)
library(edgeR)
library(readr)
library(dplyr)
library(tibble)
library(stringr)

cat("==== 1. Load expression matrix and metadata ====\n")
expr <- readRDS('/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE174409/QC/logcpm_filtered_normalized.rds')
cat("Expression matrix dimensions:", dim(expr), "\n")

meta <- readRDS('/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE174409/QC/metadata.rds')
cat("Metadata dimensions:", dim(meta), "\n")
cat("Metadata columns:", paste(colnames(meta), collapse=", "), "\n")

# Identify sample ID column in metadata
potential_id_cols <- c("title", "sample_name", "sample_id", "sample", "id", "geo_accession", "run")
found_cols <- intersect(potential_id_cols, colnames(meta))

if (length(found_cols) > 0) {
  id_col <- found_cols[1]
  cat("Using", id_col, "as sample identifier\n")
  meta$sample_id <- meta[[id_col]]
} else {
  # Find column that matches expression column names
  char_cols <- sapply(meta, function(x) is.character(x) || is.factor(x))
  char_col_names <- names(char_cols)[char_cols]
  match_counts <- sapply(char_col_names, function(col) {
    sum(as.character(meta[[col]]) %in% colnames(expr))
  })

  if (max(match_counts) > 0) {
    id_col <- char_col_names[which.max(match_counts)]
    cat("Using", id_col, "as sample identifier (", max(match_counts), "matches)\n")
    meta$sample_id <- as.character(meta[[id_col]])
  } else {
    stop("Could not identify a suitable sample ID column")
  }
}

# Check for required condition columns
req_cols <- c("diagnosis", "sex", "region")
missing_cols <- req_cols[!req_cols %in% tolower(colnames(meta))]

if (length(missing_cols) > 0) {
  cat("Warning: Missing expected columns:", paste(missing_cols, collapse=", "), "\n")
  cat("Available columns:", paste(colnames(meta), collapse=", "), "\n")
  stop("Required condition columns not found in metadata")
}

# Standardize column names (case-insensitive match)
for (col in req_cols) {
  idx <- which(tolower(colnames(meta)) == col)
  if (length(idx) > 0) {
    colnames(meta)[idx] <- col
  }
}

# Match samples between expression data and metadata
common_samples <- intersect(colnames(expr), meta$sample_id)
cat("Matched samples:", length(common_samples), "\n")

if (length(common_samples) == 0) {
  # Try to find common pattern in names
  expr_names <- colnames(expr)
  meta_names <- meta$sample_id

  cat("First few expression column names:", paste(head(expr_names), collapse=", "), "\n")
  cat("First few metadata sample IDs:", paste(head(meta_names), collapse=", "), "\n")

  stop("No matching samples between expression data and metadata")
}

# Subset data to matched samples
expr_matched <- expr[, common_samples]
meta_matched <- meta[meta$sample_id %in% common_samples, ]
# Ensure metadata is in same order as expression columns
meta_matched <- meta_matched[match(colnames(expr_matched), meta_matched$sample_id), ]

cat("==== 2. Create design matrix ====\n")
# Create model design based on diagnosis
# Adjust this section based on your specific comparison of interest
meta_matched$diagnosis <- factor(meta_matched$diagnosis)
cat("Diagnosis levels:", paste(levels(meta_matched$diagnosis), collapse=", "), "\n")

if (length(levels(meta_matched$diagnosis)) < 2) {
  stop("Need at least 2 diagnosis groups for comparison")
}

# Create a design matrix with diagnosis as the primary factor
# Also account for sex and region as covariates
design <- model.matrix(~ diagnosis + sex + region, data=meta_matched)
cat("Design matrix created with dimensions:", dim(design), "\n")
cat("Design factors:", colnames(design), "\n")

cat("==== 3. Fit linear model with limma ====\n")
fit <- lmFit(expr_matched, design)
fit <- eBayes(fit)

cat("==== 4. Extract results for all comparisons ====\n")
# Get all possible pairwise comparisons of diagnosis
diagnosis_levels <- levels(meta_matched$diagnosis)
all_comparisons <- list()

# Extract multiple comparisons if there are more than 2 diagnosis groups
if (length(diagnosis_levels) > 2) {
  for (i in 2:length(diagnosis_levels)) {
    coef_name <- paste0("diagnosis", diagnosis_levels[i])
    if (coef_name %in% colnames(design)) {
      all_comparisons[[paste0(diagnosis_levels[i], "_vs_", diagnosis_levels[1])]] <-
        topTable(fit, coef=coef_name, number=Inf, sort.by="P") %>%
        rownames_to_column("gene_id")
    }
  }

  # Combine all results
  all_results <- bind_rows(all_comparisons, .id="comparison")
  write_csv(all_results, "GSE174409/NeuroinflammationResults/DEG_all_comparisons.csv")
  cat("Saved all pairwise comparisons to GSE174409/NeuroinflammationResults/DEG_all_comparisons.csv\n")

  # Default comparison for volcano plot
  default_comparison <- names(all_comparisons)[1]
  deg_results <- all_comparisons[[default_comparison]]
  cat("Using", default_comparison, "as default comparison for DEG_results.csv\n")
} else {
  # Just one comparison
  deg_results <- topTable(fit, coef=2, number=Inf, sort.by="P") %>%
    rownames_to_column("gene_id")
}

# Create directory if it doesn't exist
dir.create("GSE174409/NeuroinflammationResults", showWarnings=FALSE, recursive=TRUE)

# Save results for volcano plot
write_csv(deg_results, "GSE174409/NeuroinflammationResults/DEG_results.csv")
cat("Saved primary DEG results to GSE174409/NeuroinflammationResults/DEG_results.csv\n")

# Show summary of results
sig_genes <- sum(deg_results$adj.P.Val < 0.05, na.rm=TRUE)
up_genes <- sum(deg_results$adj.P.Val < 0.05 & deg_results$logFC > 0, na.rm=TRUE)
down_genes <- sum(deg_results$adj.P.Val < 0.05 & deg_results$logFC < 0, na.rm=TRUE)

cat("\n==== 5. Summary of differential expression ====\n")
cat("Total genes tested:", nrow(deg_results), "\n")
cat("Significantly differential (FDR < 0.05):", sig_genes, "\n")
cat("Upregulated:", up_genes, "\n")
cat("Downregulated:", down_genes, "\n")

cat("\nDEG analysis complete. Results are ready for visualization.\n")