# Load required libraries
library(edgeR)
library(limma)
library(readr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tidyr)  # Required for pivot_wider function

# Install EnhancedVolcano if not available
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  cat("Installing EnhancedVolcano package...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("EnhancedVolcano")
}
library(EnhancedVolcano)

# Load preprocessed data from QC step
qc_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE289002/QC"
dge_filtered <- readRDS(file.path(qc_dir, "dge_filtered_normalized.rds"))
metadata <- readRDS(file.path(qc_dir, "metadata.rds"))

# Create output directory for results
results_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE289002/DE_results"
dir.create(results_dir, showWarnings = FALSE)

# Create output directory for sex-treatment focused analysis
sex_focus_dir <- file.path(results_dir, "sex_treatment_analysis")
dir.create(sex_focus_dir, showWarnings = FALSE)

# Try to get gene mapping, but continue if Ensembl is unavailable
gene_mapping_success <- FALSE
id_to_symbol <- NULL
id_to_description <- NULL

# Try using gene mapping if biomaRt is available
tryCatch({
  library(biomaRt)
  cat("Attempting to map Ensembl IDs to gene symbols...\n")

  # Try using a mirror site instead of the main Ensembl site
  ensembl_ids <- rownames(dge_filtered)

  # Try different mirrors if one fails
  mirrors <- c("useast", "uswest", "asia")
  for (mirror in mirrors) {
    tryCatch({
      cat("Trying Ensembl mirror:", mirror, "\n")
      mart <- useEnsembl(biomart = "ensembl",
                         dataset = "mmusculus_gene_ensembl",
                         mirror = mirror)

      # Get gene symbols and descriptions
      gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),
                         filters = "ensembl_gene_id",
                         values = ensembl_ids,
                         mart = mart)

      # Create lookup tables
      id_to_symbol <- setNames(gene_info$external_gene_name, gene_info$ensembl_gene_id)
      id_to_description <- setNames(gene_info$description, gene_info$ensembl_gene_id)

      # Handle missing mappings
      missing_ids <- ensembl_ids[!ensembl_ids %in% gene_info$ensembl_gene_id]
      if (length(missing_ids) > 0) {
        id_to_symbol[missing_ids] <- missing_ids
        id_to_description[missing_ids] <- ""
      }

      cat("Successfully connected to Ensembl mirror:", mirror, "\n")
      gene_mapping_success <- TRUE
    }, error = function(e) {
      cat("Failed to connect to Ensembl mirror:", mirror, "\n")
    })

    if (gene_mapping_success) break
  }

  if (!gene_mapping_success) {
    cat("All Ensembl mirrors failed. Continuing without gene symbol mapping.\n")
    # Create empty mapping to avoid errors
    id_to_symbol <- setNames(rownames(dge_filtered), rownames(dge_filtered))
    id_to_description <- setNames(rep("", length(rownames(dge_filtered))), rownames(dge_filtered))
  }

}, error = function(e) {
  cat("Error in gene mapping:", conditionMessage(e), "\n")
  cat("Continuing without gene symbol mapping.\n")
  # Create empty mapping to avoid errors
  id_to_symbol <- setNames(rownames(dge_filtered), rownames(dge_filtered))
  id_to_description <- setNames(rep("", length(rownames(dge_filtered))), rownames(dge_filtered))
})

# Function to perform differential expression analysis for a region
analyze_region <- function(region) {
  # Subset data for the current region
  region_meta <- metadata[metadata$region == region, ]
  region_samples <- region_meta$title
  region_dge <- dge_filtered[, region_samples]

  # Verify samples match
  cat("\nProcessing region:", region, "\n")
  cat("Number of samples:", ncol(region_dge), "\n")

  # Create region-specific results directory
  region_dir <- file.path(results_dir, paste0("region_", region))
  dir.create(region_dir, showWarnings = FALSE)

  # Make sure treatment and sex are factors
  region_meta$treatment <- factor(region_meta$treatment,
                               levels = c("Sal", "Mor + 24h", "Mor + 2W", "Chronic mor"))
  region_meta$sex <- factor(region_meta$sex)

  # Create design matrix
  design <- model.matrix(~ treatment + sex + treatment:sex, data = region_meta)
  colnames(design) <- make.names(colnames(design))

  # Display design matrix structure
  cat("Design matrix column names:\n")
  print(colnames(design))

  # Estimate dispersion
  region_dge <- estimateDisp(region_dge, design)

  # Fit model
  fit <- glmQLFit(region_dge, design)

  # Create contrasts with more descriptive names
  contrasts <- list(
    "Morphine_24hr_Withdrawal_vs_Saline" = makeContrasts(treatmentMor...24h, levels = design),
    "Morphine_2Week_Withdrawal_vs_Saline" = makeContrasts(treatmentMor...2W, levels = design),
    "Chronic_Morphine_vs_Saline" = makeContrasts(treatmentChronic.mor, levels = design),
    "Female_vs_Male" = makeContrasts(-sexmale, levels = design)  # Use negative to flip direction
  )

  # Interaction contrasts with more descriptive names
  interactions <- list(
    "Sex_Difference_Morphine_24hr_Withdrawal" = makeContrasts(treatmentMor...24h.sexmale, levels = design),
    "Sex_Difference_Morphine_2Week_Withdrawal" = makeContrasts(treatmentMor...2W.sexmale, levels = design),
    "Sex_Difference_Chronic_Morphine" = makeContrasts(treatmentChronic.mor.sexmale, levels = design)
  )

  # Combine all contrasts
  all_contrasts <- c(contrasts, interactions)

  # Perform differential expression for each contrast
  for (contrast_name in names(all_contrasts)) {
    contrast <- all_contrasts[[contrast_name]]

    # Test for differential expression
    qlf <- glmQLFTest(fit, contrast = contrast)

    # Extract results
    res <- topTags(qlf, n = Inf, sort.by = "PValue")
    res_df <- as.data.frame(res)

    # Add gene information
    res_df$ensembl_id <- rownames(res_df)
    if (gene_mapping_success) {
      res_df$gene_symbol <- id_to_symbol[rownames(res_df)]
      res_df$description <- id_to_description[rownames(res_df)]

      # Move columns to beginning for better readability
      res_df <- res_df[, c("ensembl_id", "gene_symbol", "description",
                          setdiff(colnames(res_df), c("ensembl_id", "gene_symbol", "description")))]
    }

    # Save results to CSV
    write.csv(res_df, file.path(region_dir, paste0(contrast_name, ".csv")), row.names = FALSE)

    # Count significant genes
    sig_genes <- sum(res_df$FDR < 0.05)
    sig_genes_fc <- sum(res_df$FDR < 0.05 & abs(res_df$logFC) > 1)
    cat(contrast_name, "- Significant genes (FDR<0.05):", sig_genes, "\n")
    cat(contrast_name, "- Significant genes (FDR<0.05 & |logFC|>1):", sig_genes_fc, "\n")
  }

  return(all_contrasts)
}

# Modified function to handle duplicate gene symbols and average samples by condition
analyze_gene_set <- function(gene_set_name, gene_symbols) {
  cat("\nAnalyzing", gene_set_name, "genes across treatment groups and sexes...\n")

  # Create directory for this gene set
  gene_set_dir <- file.path(sex_focus_dir, paste0(gene_set_name, "_genes"))
  dir.create(gene_set_dir, showWarnings = FALSE)

  # Remove duplicates from gene symbols
  unique_gene_symbols <- unique(gene_symbols)
  if (length(unique_gene_symbols) < length(gene_symbols)) {
    cat("  Note: Removed", length(gene_symbols) - length(unique_gene_symbols),
        "duplicate gene symbols\n")
    gene_symbols <- unique_gene_symbols
  }

  # Map gene symbols to Ensembl IDs if mapping is available
  if (gene_mapping_success) {
    # Find matching Ensembl IDs for gene set
    reverse_mapping <- setNames(names(id_to_symbol), id_to_symbol)
    gene_ensembl <- reverse_mapping[gene_symbols]
    gene_ensembl <- gene_ensembl[!is.na(gene_ensembl)]

    if (length(gene_ensembl) == 0) {
      cat("  No matching Ensembl IDs found for", gene_set_name, "genes\n")
      return(FALSE)
    }

    cat("  Found", length(gene_ensembl), "matching Ensembl IDs\n")

    # Get expression data for these genes
    gene_expr <- cpm(dge_filtered, log=TRUE)[gene_ensembl, ]

    # Set row names to gene symbols, ensuring uniqueness
    row_labels <- id_to_symbol[rownames(gene_expr)]
    row_labels[is.na(row_labels)] <- rownames(gene_expr)[is.na(row_labels)]

    # Make row names unique if necessary
    if (anyDuplicated(row_labels) > 0) {
      cat("  Making duplicate gene symbols unique\n")
      row_labels <- make.unique(row_labels)
    }

    rownames(gene_expr) <- row_labels

    # Calculate average expression for each sex-treatment combination
    avg_expr <- matrix(0, nrow = nrow(gene_expr),
                      ncol = 2 * length(unique(metadata$treatment)))
    rownames(avg_expr) <- rownames(gene_expr)

    # Create column names for the averaged data
    col_names <- c()
    col_idx <- 1

    # First create columns for females then males to keep same ordering
    for (sex in c("female", "male")) {
      for (treatment in unique(metadata$treatment)) {
        # Find samples for this sex-treatment combination
        samples <- metadata$title[metadata$sex == sex & metadata$treatment == treatment]

        if (length(samples) > 0) {
          # Calculate average expression for this combination
          avg_expr[, col_idx] <- rowMeans(gene_expr[, samples, drop=FALSE])
          col_names <- c(col_names, paste0(sex, "_", treatment))
          col_idx <- col_idx + 1
        }
      }
    }

    # Remove any unused columns and assign column names
    avg_expr <- avg_expr[, 1:(col_idx-1), drop = FALSE]
    colnames(avg_expr) <- col_names

    # Create annotation dataframe for the averaged data
    sex_treatment <- colnames(avg_expr)
    anno <- data.frame(
      Sex = sapply(strsplit(sex_treatment, "_"), `[`, 1),
      Treatment = sapply(strsplit(sex_treatment, "_"), function(x) paste(x[-1], collapse="_")),
      row.names = sex_treatment
    )

    # Define colors
    anno_colors <- list(
      Treatment = c("Sal" = "lightblue", "Mor + 24h" = "orange",
                   "Mor + 2W" = "darkgreen", "Chronic mor" = "red"),
      Sex = c("male" = "blue", "female" = "pink")
    )

    # Calculate gap position between sexes
    sex_change_idx <- which(diff(as.numeric(factor(anno$Sex))) != 0)

    # Generate expression heatmap with smaller color scale (-3 to 3)
    png(file.path(gene_set_dir, paste0(gene_set_name, "_expression_heatmap.png")),
        width = 900, height = 1000, res = 150)
    pheatmap(avg_expr, scale = "row",
            annotation_col = anno,
            annotation_colors = anno_colors,
            cluster_cols = FALSE,
            cluster_rows = TRUE,
            gaps_col = sex_change_idx,
            main = paste0(gene_set_name, " genes - Average expression by sex and treatment"),
            fontsize_row = 10,
            breaks = seq(-3, 3, length.out = 101))  # Set scale from -3 to 3
    dev.off()

    # Create a sex difference heatmap
    sex_diff_by_treatment <- data.frame(matrix(0, nrow=nrow(gene_expr),
                                              ncol=length(unique(metadata$treatment))))
    rownames(sex_diff_by_treatment) <- rownames(gene_expr)
    colnames(sex_diff_by_treatment) <- unique(metadata$treatment)

    # Calculate average difference between sexes for each treatment
    for (treatment in unique(metadata$treatment)) {
      male_samples <- metadata$title[metadata$sex == "male" & metadata$treatment == treatment]
      female_samples <- metadata$title[metadata$sex == "female" & metadata$treatment == treatment]

      if (length(male_samples) > 0 && length(female_samples) > 0) {
        male_avg <- rowMeans(gene_expr[, male_samples, drop=FALSE])
        female_avg <- rowMeans(gene_expr[, female_samples, drop=FALSE])
        sex_diff_by_treatment[, treatment] <- female_avg - male_avg
      }
    }

    # Generate sex difference heatmap with smaller color scale (-3 to 3)
    png(file.path(gene_set_dir, paste0(gene_set_name, "_sex_difference_heatmap.png")),
        width = 700, height = 1000, res = 150)
    pheatmap(sex_diff_by_treatment,
            color = colorRampPalette(c("blue", "white", "red"))(100),
            cluster_rows = TRUE,
            cluster_cols = FALSE,
            main = paste0(gene_set_name, " genes - Sex differences by treatment (female - male)"),
            fontsize_row = 10,
            breaks = seq(-3, 3, length.out = 101))  # Set scale from -3 to 3
    dev.off()

    # Save sex difference data
    write.csv(as.data.frame(sex_diff_by_treatment),
              file.path(gene_set_dir, paste0(gene_set_name, "_sex_differences.csv")))

    return(TRUE)
  } else {
    cat("Gene symbol mapping unavailable, skipping", gene_set_name, "analysis\n")
    return(FALSE)
  }
}

# Modified function to create sex-treatment heatmaps with averaged data
create_sex_treatment_heatmap <- function(contrast_name) {
  # Create a dataframe to store combined data across regions
  combined_samples <- c()
  combined_metadata <- data.frame()

  # Combine data across brain regions
  for (region in unique(metadata$region)) {
    # Get region-specific data
    region_meta <- metadata[metadata$region == region, ]
    region_samples <- region_meta$title

    # Load contrast results
    region_dir <- file.path(results_dir, paste0("region_", region))
    contrast_file <- file.path(region_dir, paste0(contrast_name, ".csv"))

    if (file.exists(contrast_file)) {
      # Add to combined data
      combined_samples <- c(combined_samples, region_samples)
      combined_metadata <- rbind(combined_metadata, region_meta)
    }
  }

  if (length(combined_samples) == 0) {
    cat("No data found for contrast:", contrast_name, "\n")
    return(FALSE)
  }

  # Use a single contrast file to get top genes
  sample_region <- unique(metadata$region)[1]
  sample_region_dir <- file.path(results_dir, paste0("region_", sample_region))
  sample_contrast_file <- file.path(sample_region_dir, paste0(contrast_name, ".csv"))

  if (file.exists(sample_contrast_file)) {
    de_results <- read.csv(sample_contrast_file)

    # Get top 50 genes from this contrast
    top_genes <- de_results$ensembl_id[order(de_results$FDR)[1:min(50, nrow(de_results))]]

    # Extract expression data across all regions
    all_expr <- cpm(dge_filtered[, combined_samples], log=TRUE)[top_genes, ]

    # Use gene symbols for row names if available
    if (gene_mapping_success) {
      row_labels <- id_to_symbol[rownames(all_expr)]
      row_labels[is.na(row_labels)] <- rownames(all_expr)[is.na(row_labels)]
      rownames(all_expr) <- row_labels
    }

    # Calculate average expression for each region-sex-treatment combination
    unique_regions <- unique(combined_metadata$region)
    unique_treatments <- unique(combined_metadata$treatment)
    unique_sexes <- unique(combined_metadata$sex)

    # Create empty matrix for averaged data
    avg_expr <- matrix(0,
                     nrow = nrow(all_expr),
                     ncol = length(unique_regions) * length(unique_treatments) * length(unique_sexes))
    rownames(avg_expr) <- rownames(all_expr)

    # Create column names and calculate averages
    col_names <- c()
    col_idx <- 1

    # Create mapping from condition to column
    condition_map <- data.frame(
      region = character(),
      sex = character(),
      treatment = character(),
      col_idx = integer(),
      stringsAsFactors = FALSE
    )

    # First create female, then male columns to maintain order
    for (region in unique_regions) {
      for (sex in c("female", "male")) {
        for (treatment in levels(factor(unique_treatments))) {
          # Find samples for this region-sex-treatment combination
          samples <- combined_metadata$title[
            combined_metadata$region == region &
            combined_metadata$sex == sex &
            combined_metadata$treatment == treatment
          ]

          if (length(samples) > 0) {
            # Calculate average expression for this combination
            avg_expr[, col_idx] <- rowMeans(all_expr[, samples, drop=FALSE])
            col_name <- paste(region, sex, treatment, sep="_")
            col_names <- c(col_names, col_name)

            # Add to condition map
            condition_map <- rbind(condition_map, data.frame(
              region = region,
              sex = sex,
              treatment = treatment,
              col_idx = col_idx,
              stringsAsFactors = FALSE
            ))

            col_idx <- col_idx + 1
          }
        }
      }
    }

    # Remove unused columns and set column names
    avg_expr <- avg_expr[, 1:(col_idx-1), drop = FALSE]
    colnames(avg_expr) <- col_names

    # Create annotation dataframe
    anno <- data.frame(
      Region = sapply(strsplit(colnames(avg_expr), "_"), `[`, 1),
      Sex = sapply(strsplit(colnames(avg_expr), "_"), `[`, 2),
      Treatment = sapply(strsplit(colnames(avg_expr), "_"), function(x) paste(x[-(1:2)], collapse = "_")),
      row.names = colnames(avg_expr),
      stringsAsFactors = FALSE
    )

    # Define colors
    anno_colors <- list(
      Region = setNames(rainbow(length(unique_regions)), unique_regions),
      Treatment = c("Sal" = "lightblue", "Mor + 24h" = "orange",
                   "Mor + 2W" = "darkgreen", "Chronic mor" = "red"),
      Sex = c("male" = "blue", "female" = "pink")
    )

    # Order columns: first by region, then by sex, then by treatment
    ordered_cols <- order(anno$Region, anno$Sex, anno$Treatment)
    ordered_expr <- avg_expr[, ordered_cols]
    ordered_anno <- anno[ordered_cols, ]

    # Calculate gap positions between regions and sexes
    region_change_idx <- which(diff(as.numeric(factor(ordered_anno$Region))) != 0)
    sex_change_idx <- which(diff(as.numeric(factor(ordered_anno$Sex))) != 0)
    all_gaps <- sort(unique(c(region_change_idx, sex_change_idx)))

    # Generate heatmap with averaged data and smaller color scale (-3 to 3)
    png(file.path(sex_focus_dir, paste0(contrast_name, "_sex_treatment_heatmap.png")),
        width = 900, height = 1000, res = 150)
    pheatmap(ordered_expr, scale = "row",
            annotation_col = ordered_anno,
            annotation_colors = anno_colors,
            cluster_cols = FALSE,  # Don't cluster columns to preserve ordering
            cluster_rows = TRUE,   # Still cluster rows by expression pattern
            gaps_col = all_gaps,   # Add visual gaps between groups
            main = paste0(contrast_name, " - Average gene expression by region, sex & treatment"),
            fontsize_row = 8,
            show_colnames = FALSE,
            breaks = seq(-3, 3, length.out = 101))  # Set scale from -3 to 3
    dev.off()

    # Create a separate file just for the sex differences
    # Fix: Ensure row names don't contain missing values before creating data frame
    row_names <- rownames(all_expr)
    is_invalid <- is.na(row_names) | row_names == ""
    if (any(is_invalid)) {
      # Replace invalid row names with a placeholder
      row_names[is_invalid] <- paste0("Gene_", which(is_invalid))
      rownames(all_expr) <- row_names
    }
    sex_diff_summary <- data.frame(Gene = row_names, stringsAsFactors = FALSE)

    # Calculate and add sex differences for each region and treatment
    for (region in unique_regions) {
      for (treatment in levels(factor(unique_treatments))) {
        # Find the column indices for this region/treatment combination
        female_col <- which(condition_map$region == region &
                           condition_map$sex == "female" &
                           condition_map$treatment == treatment)
        male_col <- which(condition_map$region == region &
                         condition_map$sex == "male" &
                         condition_map$treatment == treatment)

        if (length(female_col) > 0 && length(male_col) > 0) {
          # Calculate sex difference (female - male)
          diff_col_name <- paste0(region, "_", treatment, "_sexDiff")
          sex_diff_summary[[diff_col_name]] <- avg_expr[, female_col] - avg_expr[, male_col]
        }
      }
    }

    # Save the summary
    write.csv(sex_diff_summary,
              file.path(sex_focus_dir, paste0(contrast_name, "_sex_diff_summary.csv")),
              row.names = FALSE)

    return(TRUE)
  }

  return(FALSE)
}

# Process each brain region to generate basic results files
all_contrasts <- NULL
for (region in unique(metadata$region)) {
  # Process this region and get the contrasts
  region_contrasts <- analyze_region(region)

  # Save contrasts from first region
  if (is.null(all_contrasts)) {
    all_contrasts <- region_contrasts
  }
}

# Create sex-focused heatmaps for each contrast
cat("\n\nGenerating sex-focused visualizations across all brain regions...\n")
for (contrast_name in names(all_contrasts)) {
  cat("Creating sex-treatment heatmap for", contrast_name, "\n")
  create_sex_treatment_heatmap(contrast_name)
}

# Create overall MDS plot grouped by sex and treatment
cat("\nCreating overall MDS plot grouped by sex and treatment...\n")
png(file.path(sex_focus_dir, "overall_MDS_plot_by_sex_treatment.png"), width = 800, height = 700, res = 150)
plotMDS(dge_filtered,
        col = as.numeric(metadata$treatment) + (as.numeric(metadata$sex) - 1) * 4,
        main = "MDS Plot Grouped by Sex and Treatment")
legend("topright",
       legend = c(paste(rep(c("Female", "Male"), each = 4),
                        rep(levels(factor(metadata$treatment)), 2), sep="-")),
       col = 1:8,
       pch = 16)
dev.off()

# Set path to KEGG pathway files
kegg_files_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE289002/kegg_outputs"
cat("\nReading KEGG pathway files from:", kegg_files_dir, "\n")

# Get all CSV files in the KEGG directory
kegg_files <- list.files(kegg_files_dir, pattern = "\\.csv$", full.names = TRUE)

if (length(kegg_files) > 0) {
  cat("Found", length(kegg_files), "KEGG pathway files\n")

  for (kegg_file in kegg_files) {
    # Extract pathway name from filename
    pathway_name <- gsub("\\.csv$", "", basename(kegg_file))
    cat("Processing pathway:", pathway_name, "\n")

    # Read gene list from CSV file
    tryCatch({
      pathway_genes_df <- read.csv(kegg_file, stringsAsFactors = FALSE)

      # Determine which column contains gene symbols
      possible_gene_cols <- c("gene_symbol", "symbol", "gene", "Gene", "Symbol",
                             "GeneSymbol", "SYMBOL", "gene_name", "Name")
      gene_col <- NULL
      for (col in possible_gene_cols) {
        if (col %in% colnames(pathway_genes_df)) {
          gene_col <- col
          break
        }
      }

      if (is.null(gene_col) && ncol(pathway_genes_df) > 0) {
        # If no recognized column name, use the first column
        gene_col <- colnames(pathway_genes_df)[1]
      }

      if (is.null(gene_col)) {
        cat("  Warning: Could not identify gene column in", pathway_name, "\n")
        next
      }

      # Extract gene symbols
      pathway_genes <- pathway_genes_df[[gene_col]]
      pathway_genes <- pathway_genes[!is.na(pathway_genes) & pathway_genes != ""]

      if (length(pathway_genes) == 0) {
        cat("  Warning: No genes found in", pathway_name, "\n")
        next
      }

      cat("  Found", length(pathway_genes), "genes in pathway\n")

      # Analyze this pathway's genes
      analyze_gene_set(pathway_name, pathway_genes)
    }, error = function(e) {
      cat("  Error processing", pathway_name, ":", conditionMessage(e), "\n")
    })
  }
} else {
  cat("No KEGG pathway files found in:", kegg_files_dir, "\n")
  cat("Falling back to pre-defined gene sets...\n")

  # Fall back to predefined gene sets if no KEGG files are found
  complement_genes <- c("C1qa", "C1qb", "C1qc", "C3", "C4b", "Cfb", "Cfd",
                       "C1ra", "C1rb", "C1s", "C2", "C4", "C5", "C6", "C7", "C8a",
                       "C8b", "C8g", "C9", "Cfh", "Cfi", "Cd55", "Cd59a")

  tlr_genes <- c("Tlr1", "Tlr2", "Tlr3", "Tlr4", "Tlr5", "Tlr6", "Tlr7", "Tlr8", "Tlr9",
                 "Tlr11", "Tlr12", "Tlr13", "Myd88", "Irak1", "Irak4", "Traf6")

  # Analyze specific gene sets
  analyze_gene_set("Complement", complement_genes)
  analyze_gene_set("TLR", tlr_genes)
}

cat("\nAnalysis complete. Results saved in:", results_dir, "\n")
cat("Sex-specific analysis results saved in:", sex_focus_dir, "\n")