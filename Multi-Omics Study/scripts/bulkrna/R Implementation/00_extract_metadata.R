#!/usr/bin/env Rscript

# =============================================================================
# Extract Metadata from GEO Series Matrix for GSE174409
# =============================================================================
# 
# Description: Extract and format sample metadata from GSE174409_series_matrix.txt
#              for use in bulk RNA-seq analysis
# 
# Author: Multi-Omics Study Team
# Date: 2024
# 
# =============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

# Set file paths
series_matrix_file <- "../../data/raw/bulkrna/GSE174409_series_matrix.txt"
output_file <- "../../data/raw/bulkrna/sample_metadata.csv"

cat("Extracting metadata from GEO series matrix...\n")

# Read the series matrix file
lines <- readLines(series_matrix_file)

# Find lines with sample information
sample_title_line <- grep("^!Sample_title", lines)[1]
sample_source_line <- grep("^!Sample_source_name_ch1", lines)[1]
sample_characteristics_lines <- grep("^!Sample_characteristics_ch1", lines)

# Extract sample titles (sample IDs)
sample_titles <- lines[sample_title_line]
sample_titles <- str_split(sample_titles, "\t")[[1]][-1]  # Remove first element
sample_titles <- str_remove_all(sample_titles, '"')

# Extract brain regions
brain_regions <- lines[sample_source_line]
brain_regions <- str_split(brain_regions, "\t")[[1]][-1]
brain_regions <- str_remove_all(brain_regions, '"')

# Extract characteristics
characteristics <- list()
for (i in seq_along(sample_characteristics_lines)) {
  char_line <- lines[sample_characteristics_lines[i]]
  char_values <- str_split(char_line, "\t")[[1]][-1]
  char_values <- str_remove_all(char_values, '"')
  
  # Extract the characteristic name from the first value
  char_name <- str_extract(char_values[1], "^[^:]+")
  
  # Extract values (everything after the colon)
  char_data <- str_remove(char_values, "^[^:]+: ")
  
  characteristics[[char_name]] <- char_data
}

# Create metadata data frame
metadata <- data.frame(
  sample_id = sample_titles,
  region = brain_regions,
  stringsAsFactors = FALSE
)

# Add characteristics to metadata
for (char_name in names(characteristics)) {
  if (char_name == "diagnosis") {
    # Map diagnosis to condition
    metadata$condition <- ifelse(characteristics[[char_name]] == "OUD", "OUD", "Control")
  } else if (char_name == "Sex") {
    metadata$sex <- characteristics[[char_name]]
  } else if (char_name == "age") {
    metadata$age <- as.numeric(characteristics[[char_name]])
  } else if (char_name == "race") {
    metadata$race <- characteristics[[char_name]]
  } else if (char_name == "pmi") {
    metadata$pmi <- as.numeric(characteristics[[char_name]])
  } else if (char_name == "avgph") {
    metadata$ph <- as.numeric(characteristics[[char_name]])
  } else if (char_name == "rpf_rin") {
    metadata$rin <- as.numeric(characteristics[[char_name]])
  } else if (char_name == "tissuestoragetime_month") {
    metadata$storage_time_months <- as.numeric(characteristics[[char_name]])
  }
}

# Create batch variable based on region (simple batching strategy)
metadata$batch <- ifelse(metadata$region == "NAC", "Batch_1", "Batch_2")

# Reorder columns for clarity
metadata <- metadata %>%
  select(sample_id, condition, sex, age, region, race, pmi, ph, rin, 
         storage_time_months, batch) %>%
  arrange(sample_id)

# Display summary
cat("Metadata extraction summary:\n")
cat("Total samples:", nrow(metadata), "\n")
cat("Conditions:", table(metadata$condition), "\n")
cat("Regions:", table(metadata$region), "\n")
cat("Sex distribution:", table(metadata$sex), "\n")

# Save metadata
write_csv(metadata, output_file)
cat("Metadata saved to:", output_file, "\n")

# Display first few rows
cat("\nFirst 10 rows of extracted metadata:\n")
print(head(metadata, 10))

cat("\n✅ Metadata extraction completed successfully!\n")