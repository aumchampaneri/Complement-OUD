# Analysis Outputs_0 - Pathview visualization of DESeq2 results
# This script visualizes gene expression data on KEGG pathways

# Load required libraries
library(clusterProfiler)
library(org.Mm.eg.db)  # Mouse database
library(pathview)
library(dplyr)
library(readr)
library(stringr)

# Set base directory and create output directory
base_dir <- "/"
deseq2_dir <- file.path(base_dir, "GSE207128/DESeq2_outputs")
pathview_dir <- file.path(base_dir, "GSE207128/Pathview_results")
dir.create(pathview_dir, showWarnings = FALSE, recursive = TRUE)

# Find all DESeq2 result files
deseq2_files <- list.files(deseq2_dir, pattern = "deseq2_results_.*\\.csv$", full.names = TRUE)
cat("Found", length(deseq2_files), "DESeq2 result files\n")

# Function to convert mouse gene IDs to KEGG IDs (necessary for pathview)
convert_to_kegg_ids <- function(gene_ids, gene_symbols, fold_changes) {
  # Create a data frame for results
  result <- data.frame(
    original_id = gene_ids,
    gene_symbol = gene_symbols,
    fold_change = fold_changes,
    kegg_id = NA_character_,
    stringsAsFactors = FALSE
  )

  # Clean Ensembl IDs (remove version numbers)
  clean_ids <- gsub("\\..*$", "", gene_ids)

  # First get ENTREZID (intermediate step)
  entrez_ids <- mapIds(org.Mm.eg.db,
                       keys = clean_ids,
                       column = "ENTREZID",
                       keytype = "ENSEMBL",
                       multiVals = "first")

  valid_entrez <- !is.na(entrez_ids)

  # Then map ENTREZID to KEGG
  if(sum(valid_entrez) > 0) {
    valid_entrez_ids <- entrez_ids[valid_entrez]
    kegg_ids <- mapIds(org.Mm.eg.db,
                       keys = valid_entrez_ids,
                       column = "PATH",
                       keytype = "ENTREZID",
                       multiVals = "first")

    # Use direct KEGG ID format for pathview (remove "path:" prefix)
    kegg_ids <- gsub("path:", "", kegg_ids)

    # Update the result data frame
    result$kegg_id[valid_entrez] <- kegg_ids
  }

  # Try using symbols for unmapped genes
  if(sum(is.na(result$kegg_id)) > 0) {
    unmapped_idx <- which(is.na(result$kegg_id))
    symbols_to_try <- gene_symbols[unmapped_idx]

    # Get ENTREZID from symbols
    symbol_entrez <- mapIds(org.Mm.eg.db,
                           keys = symbols_to_try,
                           column = "ENTREZID",
                           keytype = "SYMBOL",
                           multiVals = "first")

    valid_symbols <- !is.na(symbol_entrez)

    if(sum(valid_symbols) > 0) {
      # Then map to KEGG
      valid_symbol_entrez <- symbol_entrez[valid_symbols]
      symbol_kegg <- mapIds(org.Mm.eg.db,
                           keys = valid_symbol_entrez,
                           column = "PATH",
                           keytype = "ENTREZID",
                           multiVals = "first")

      symbol_kegg <- gsub("path:", "", symbol_kegg)
      result$kegg_id[unmapped_idx[valid_symbols]] <- symbol_kegg
    }
  }

  # Return only rows with valid KEGG IDs
  result <- result[!is.na(result$kegg_id), ]

  if(nrow(result) > 0) {
    cat("Successfully mapped", nrow(result), "genes to KEGG IDs\n")
    return(result)
  } else {
    cat("Failed to map genes to KEGG IDs\n")
    return(NULL)
  }
}

# Pathways of interest - complement and inflammation
pathways_of_interest <- c(
  "mmu04610", # Complement and coagulation cascades
  "mmu04620", # Toll-like receptor signaling
  "mmu04621", # NOD-like receptor signaling
  "mmu04622", # RIG-I-like receptor signaling
  "mmu04623", # Cytosolic DNA-sensing
  "mmu04625", # C-type lectin receptor signaling
  "mmu04060", # Cytokine-cytokine receptor interaction
  "mmu04066", # HIF-1 signaling
  "mmu04668", # TNF signaling
  "mmu04657", # IL-17 signaling
  "mmu04217", # Necroptosis
  "mmu04064", # NF-kappa B signaling
  "mmu04151"  # PI3K-Akt signaling
)

# Process each DESeq2 result file
for (file in deseq2_files) {
  # Extract cell type from filename
  cell_type <- str_extract(basename(file), "(?<=deseq2_results_).*(?=_OUD_vs_Normal)")
  cat("\nAnalyzing:", cell_type, "\n")

  # Read DESeq2 results
  deseq_results <- read_csv(file)

  # Check if we have at least some genes
  if (nrow(deseq_results) < 50) {
    cat("Too few genes for pathway visualization, skipping...\n")
    next
  }

  # Map to KEGG IDs
  kegg_mapping <- convert_to_kegg_ids(
    deseq_results$gene_id,
    deseq_results$gene_symbol,
    deseq_results$log2FoldChange
  )

  if (is.null(kegg_mapping) || nrow(kegg_mapping) < 20) {
    cat("Failed to map enough genes to KEGG IDs, skipping...\n")
    next
  }

  # Create named vector for pathview
  gene_data <- kegg_mapping$fold_change
  names(gene_data) <- kegg_mapping$kegg_id

  # Create directory for this cell type
  cell_type_dir <- file.path(pathview_dir, gsub(" ", "_", cell_type))
  dir.create(cell_type_dir, showWarnings = FALSE)

  # Generate pathway visualizations
  for (pathway in pathways_of_interest) {
    tryCatch({
      # Set KEGG species (mouse)
      path_out <- file.path(cell_type_dir, paste0(pathway, "_", gsub(" ", "_", cell_type)))

      # Create the pathway visualization
      pathview(
        gene.data = gene_data,
        pathway.id = pathway,
        species = "mmu",  # Mouse
        out.suffix = gsub(" ", "_", cell_type),
        kegg.dir = cell_type_dir,
        limit = list(gene = c(-2, 2)),  # Color limits for fold change
        low = "blue",  # Downregulated genes
        mid = "white", # No change
        high = "red"   # Upregulated genes
      )

      cat("Generated pathway visualization for", pathway, "\n")
    }, error = function(e) {
      cat("Error generating pathway visualization for", pathway, ":", e$message, "\n")
    })
  }
}

cat("\nPathview analysis complete. Results saved in:", pathview_dir, "\n")