# 00_Metadata.R - Fixed version to properly preserve Ensembl gene IDs

# Load libraries
library(dplyr)
library(stringr)
library(readr)
library(tibble)
library(data.table)
library(edgeR)  # Added for DGE object creation

# === 1. Load gene expression matrix ===
expr_file <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE174409/GSE174409_raw_counts_02102020.csv"
expr_df <- read_csv(expr_file)

# Check the format of the first column (should contain Ensembl IDs)
cat("First few row identifiers:\n")
print(head(expr_df[[1]]))

# Explicitly set rownames to be the Ensembl IDs in first column
ensembl_ids <- expr_df[[1]]
expr_df <- expr_df[, -1]  # Remove the first column
rownames(expr_df) <- ensembl_ids  # Set rownames to Ensembl IDs

cat("Expression matrix dimensions:\n")
print(dim(expr_df))
print(head(expr_df[, 1:5]))

# Verify that rownames are properly set to Ensembl IDs
cat("\nVerifying row names contain Ensembl IDs:\n")
print(head(rownames(expr_df)))

# === 2. Extract sample information ===
sample_names <- colnames(expr_df)
metadata_df <- data.frame(
  sample_id = sample_names,
  stringsAsFactors = FALSE
)

# Extract region from sample names
metadata_df$region <- case_when(
  str_detect(sample_names, "^HN") ~ "NAC",
  str_detect(sample_names, "^HBA") ~ "DLPFC",
  TRUE ~ NA_character_
)

# Extract sample numbers
metadata_df$sample_num <- as.integer(str_extract(sample_names, "\\d+"))

# Assign diagnosis based on sample numbers (assuming 1-20 are controls, 21-40 are OUD)
metadata_df$diagnosis <- case_when(
  metadata_df$sample_num <= 20 ~ "CONT",
  metadata_df$sample_num > 20 ~ "OUD",
  TRUE ~ NA_character_
)

# === 3-6. Process additional metadata (unchanged) ===
# [All your existing metadata processing code here]
series_file <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE174409/GSE174409_series_matrix.txt"
cat("Reading series matrix file for additional metadata...\n")

# Function to extract metadata from series matrix file
extract_series_metadata <- function(file_path) {
  lines <- readLines(file_path)

  # Get GSM IDs
  gsm_line <- grep("!Sample_geo_accession", lines, value = TRUE)
  gsm_ids <- unlist(strsplit(gsm_line, "\t"))[-1]
  gsm_ids <- gsub("\"", "", gsm_ids)

  # Create dataframe with GSM IDs
  metadata <- data.frame(gsm_id = gsm_ids, stringsAsFactors = FALSE)

  # Extract title and source
  for(pattern in c("!Sample_title", "!Sample_source_name_ch1")) {
    line <- grep(pattern, lines, value = TRUE)
    if(length(line) > 0) {
      parts <- unlist(strsplit(line, "\t"))
      col_name <- tolower(gsub("!Sample_|_ch1", "", pattern))
      values <- gsub("\"", "", parts[-1])
      metadata[[col_name]] <- values
    }
  }

  # Process characteristic lines
  char_lines <- lines[grep("!Sample_characteristics_ch1", lines)]

  for(line in char_lines) {
    parts <- unlist(strsplit(line, "\t"))
    values <- gsub("\"", "", parts[-1])

    # Check if the first value has a "label: value" format
    first_value <- values[1]
    if(length(first_value) > 0 && grepl(": ", first_value)) {
      # Extract the label from the first value
      meta_type <- trimws(sub(":.*$", "", first_value))

      # Clean column name
      col_name <- gsub("[^a-zA-Z0-9_]", "_", tolower(meta_type))
      col_name <- gsub("^_+|_+$", "", col_name)

      # Make column name unique if already exists
      if(col_name %in% names(metadata)) {
        col_name <- paste0(col_name, "_", sum(grepl(paste0("^", col_name), names(metadata))))
      }

      # Extract all values (removing the type prefix)
      clean_values <- sapply(values, function(x) {
        if(grepl(": ", x)) {
          return(trimws(sub("^.*?: ", "", x)))
        } else {
          return(NA)
        }
      })

      # Add to metadata
      metadata[[col_name]] <- clean_values
    } else {
      # For lines without "label: value" format
      label <- gsub("!Sample_characteristics_ch1: ", "", parts[1])
      label <- gsub("\"", "", label)

      # If label is missing, use a generic name
      if(label == "") {
        label <- paste0("characteristic_", sum(grepl("^characteristic_", names(metadata))) + 1)
      }

      # Clean column name
      col_name <- gsub("[^a-zA-Z0-9_]", "_", tolower(label))
      col_name <- gsub("^_+|_+$", "", col_name)

      # Add to metadata
      metadata[[col_name]] <- values
    }
  }

  # Print available metadata columns
  cat("\nExtracted metadata columns from series matrix:\n")
  print(names(metadata))

  return(metadata)
}

# Get metadata from series matrix
series_metadata <- tryCatch({
  extract_series_metadata(series_file)
}, error = function(e) {
  cat("Error reading series matrix file:", e$message, "\n")
  return(NULL)
})

# === 4. Map GSM IDs to sample names ===
if(!is.null(series_metadata)) {
  # Create a direct mapping between GSM IDs and sample names based on order
  mapping <- data.frame(
    gsm_id = series_metadata$gsm_id,
    sample_id = sample_names,
    stringsAsFactors = FALSE
  )

  # Add this mapping to series_metadata
  series_metadata <- merge(series_metadata, mapping, by="gsm_id")

  # Now extract any useful metadata columns from series_metadata
  exclude_cols <- c("gsm_id", "sample_id", "title", "source")
  useful_cols <- setdiff(colnames(series_metadata), exclude_cols)

  if(length(useful_cols) > 0) {
    # Create a lookup table
    lookup <- series_metadata[, c("sample_id", useful_cols)]
    rownames(lookup) <- lookup$sample_id
    lookup$sample_id <- NULL

    # Add these columns to metadata_df
    for(col in useful_cols) {
      metadata_df[[col]] <- lookup[metadata_df$sample_id, col]
    }

    cat("Added", length(useful_cols), "metadata columns from series matrix file\n")

    # Optional: Convert common metadata fields to appropriate types
    numeric_cols <- c("age", "pmi", "rin", "ph")
    for(col in intersect(numeric_cols, colnames(metadata_df))) {
      metadata_df[[col]] <- as.numeric(metadata_df[[col]])
    }

    factor_cols <- c("sex", "race", "manner_of_death")
    for(col in intersect(factor_cols, colnames(metadata_df))) {
      metadata_df[[col]] <- as.factor(metadata_df[[col]])
    }
  } else {
    cat("No additional metadata columns found in series matrix file\n")
  }
}

# === 5. Check metadata columns ===
cat("\nFinal metadata columns:\n")
print(colnames(metadata_df))

# === 6. Check result ===
cat("\nMetadata preview:\n")
print(head(metadata_df))

# Count samples by group
cat("\nSample counts by group:\n")
print(table(metadata_df$diagnosis, metadata_df$region))

# === 7. Create and save processed data for downstream analysis ===
# Create directories for saving processed data
base_dir <- "/GSE174409 Hs"
raw_dir <- file.path(base_dir, "Data")
qc_dir <- file.path(base_dir, "QC")
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

# Save raw expression matrix with proper Ensembl IDs
saveRDS(expr_df, file.path(raw_dir, "rawcounts.rds"))
write.csv(expr_df, file.path(raw_dir, "rawcounts.csv"))

# Save metadata
saveRDS(metadata_df, file.path(qc_dir, "metadata.rds"))

# Create DGEList object
cat("\nCreating DGEList object...\n")
dge <- DGEList(counts = as.matrix(expr_df), samples = metadata_df)

# Basic filtering (remove genes with very low counts)
keep <- filterByExpr(dge, group = metadata_df$diagnosis)
dge_filtered <- dge[keep, , keep.lib.sizes=FALSE]
cat("Filtered genes:", nrow(dge), "->", nrow(dge_filtered), "\n")

# Normalize
dge_filtered <- calcNormFactors(dge_filtered, method = "TMM")

# Calculate log-CPM values
logcpm_filtered_norm <- cpm(dge_filtered, log = TRUE, prior.count = 2)

# Verify gene IDs are preserved
cat("\nDGE object gene IDs check (first 5):\n")
print(head(rownames(dge_filtered), 5))
cat("\nlogCPM matrix gene IDs check (first 5):\n")
print(head(rownames(logcpm_filtered_norm), 5))

# Save processed objects
saveRDS(dge_filtered, file.path(qc_dir, "dge_filtered_normalized.rds"))
saveRDS(logcpm_filtered_norm, file.path(qc_dir, "logcpm_filtered_normalized.rds"))

cat("\nPreprocessing complete. All files saved with proper Ensembl gene IDs.\n")
cat("Files saved:\n")
cat("- Raw counts:", file.path(raw_dir, "rawcounts.rds"), "\n")
cat("- Metadata:", file.path(qc_dir, "metadata.rds"), "\n")
cat("- Filtered DGE object:", file.path(qc_dir, "dge_filtered_normalized.rds"), "\n")
cat("- Log-CPM values:", file.path(qc_dir, "logcpm_filtered_normalized.rds"), "\n")