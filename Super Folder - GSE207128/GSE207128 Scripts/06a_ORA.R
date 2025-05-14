# GSE207128 - Over-Representation Analysis (ORA) of DESeq2 results
# This script performs pathway enrichment analysis on DEGs identified by DESeq2

# Load required libraries
library(clusterProfiler)
library(org.Mm.eg.db)  # Changed to mouse database instead of human
library(enrichplot)
library(DOSE)
library(ReactomePA)
library(dplyr)
library(readr)
library(ggplot2)
library(stringr)

# Set base directory and create output directory
base_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD"
deseq2_dir <- file.path(base_dir, "GSE207128/DESeq2_outputs")
ora_dir <- file.path(base_dir, "GSE207128/ORA_results")
dir.create(ora_dir, showWarnings = FALSE, recursive = TRUE)

# Find all DESeq2 result files
deseq2_files <- list.files(deseq2_dir, pattern = "deseq2_results_.*\\.csv$", full.names = TRUE)
cat("Found", length(deseq2_files), "DESeq2 result files\n")

# Function to convert mouse gene IDs to Entrez IDs
convert_to_entrez <- function(gene_ids, gene_symbols) {
  # Examine some IDs to understand format
  cat("Sample gene IDs:", head(gene_ids, 3), "\n")
  cat("Sample gene symbols:", head(gene_symbols, 3), "\n")

  # Create mapping result dataframe
  result <- data.frame(
    original_id = gene_ids,
    gene_symbol = gene_symbols,
    entrezid = NA_character_,
    stringsAsFactors = FALSE
  )

  # Try using clean ENS IDs first
  clean_ids <- gsub("\\..*$", "", gene_ids)  # Remove version numbers

  # Map from Ensembl ID to Entrez ID
  entrez_ids <- mapIds(org.Mm.eg.db,
                       keys = clean_ids,
                       column = "ENTREZID",
                       keytype = "ENSEMBL",
                       multiVals = "first")

  valid_ids <- !is.na(entrez_ids)
  if(sum(valid_ids) > 0) {
    cat("Converted", sum(valid_ids), "genes using Ensembl IDs\n")
    result$entrezid[valid_ids] <- entrez_ids[valid_ids]
  }

  # Try using symbols for any remaining unmapped genes
  if(sum(is.na(result$entrezid)) > 0 && !all(is.na(gene_symbols))) {
    unmapped_idx <- which(is.na(result$entrezid))
    symbols_to_try <- gene_symbols[unmapped_idx]
    symbol_entrez <- mapIds(org.Mm.eg.db,
                           keys = symbols_to_try,
                           column = "ENTREZID",
                           keytype = "SYMBOL",
                           multiVals = "first")

    valid_symbols <- !is.na(symbol_entrez)
    if(sum(valid_symbols) > 0) {
      cat("Converted", sum(valid_symbols), "additional genes using symbols\n")
      result$entrezid[unmapped_idx[valid_symbols]] <- symbol_entrez[valid_symbols]
    }
  }

  # Return only rows with valid entrez IDs
  result <- result[!is.na(result$entrezid), ]

  if(nrow(result) > 0) {
    cat("Successfully converted", nrow(result), "out of", length(gene_ids), "genes\n")
    return(result)
  } else {
    cat("Failed to convert any genes to Entrez IDs\n")
    return(NULL)
  }
}

# Function to perform ORA analysis
run_ora_analysis <- function(gene_list, universe, output_prefix) {
  results <- list()
  min_genes <- 3

  # GO Biological Process enrichment
  tryCatch({
    go_bp <- enrichGO(gene = gene_list,
                      universe = universe,
                      OrgDb = org.Mm.eg.db,  # Changed to mouse
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2,
                      keyType = "ENTREZID")

    if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
      results[["GO_BP"]] <- go_bp
      pdf(file.path(ora_dir, paste0(output_prefix, "_GO_BP.pdf")), width = 10, height = c(8, 12))
      print(dotplot(go_bp, showCategory = 20))
      dev.off()

      # Save results table
      write.csv(go_bp@result, file.path(ora_dir, paste0(output_prefix, "_GO_BP.csv")), row.names = FALSE)
    }
  }, error = function(e) {
    message("Error in GO BP analysis: ", e$message)
  })

  # GO Molecular Function enrichment
  tryCatch({
    go_mf <- enrichGO(gene = gene_list,
                      universe = universe,
                      OrgDb = org.Mm.eg.db,  # Changed to mouse
                      ont = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2,
                      keyType = "ENTREZID")

    if (!is.null(go_mf) && nrow(go_mf@result) > 0) {
      results[["GO_MF"]] <- go_mf
      pdf(file.path(ora_dir, paste0(output_prefix, "_GO_MF.pdf")), width = 10, height = c(8, 12))
      print(dotplot(go_mf, showCategory = 20))
      dev.off()

      # Save results table
      write.csv(go_mf@result, file.path(ora_dir, paste0(output_prefix, "_GO_MF.csv")), row.names = FALSE)
    }
  }, error = function(e) {
    message("Error in GO MF analysis: ", e$message)
  })

  # GO Cellular Component enrichment
  tryCatch({
    go_cc <- enrichGO(gene = gene_list,
                      universe = universe,
                      OrgDb = org.Mm.eg.db,  # Changed to mouse
                      ont = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2,
                      keyType = "ENTREZID")

    if (!is.null(go_cc) && nrow(go_cc@result) > 0) {
      results[["GO_CC"]] <- go_cc
      pdf(file.path(ora_dir, paste0(output_prefix, "_GO_CC.pdf")), width = 10, height = c(8, 12))
      print(dotplot(go_cc, showCategory = 20))
      dev.off()

      # Save results table
      write.csv(go_cc@result, file.path(ora_dir, paste0(output_prefix, "_GO_CC.csv")), row.names = FALSE)
    }
  }, error = function(e) {
    message("Error in GO CC analysis: ", e$message)
  })

  # KEGG pathway enrichment
  tryCatch({
    kegg <- enrichKEGG(gene = gene_list,
                       universe = universe,
                       organism = "mmu",  # Changed to mouse
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

    if (!is.null(kegg) && nrow(kegg@result) > 0) {
      results[["KEGG"]] <- kegg
      pdf(file.path(ora_dir, paste0(output_prefix, "_KEGG.pdf")), width = 10, height = c(8, 12))
      print(dotplot(kegg, showCategory = 20))
      dev.off()

      # Save results table
      write.csv(kegg@result, file.path(ora_dir, paste0(output_prefix, "_KEGG.csv")), row.names = FALSE)
    }
  }, error = function(e) {
    message("Error in KEGG analysis: ", e$message)
  })

  # Skip Reactome (less reliable for mouse)

  return(results)
}

# Process each DESeq2 result file
for (file in deseq2_files) {
  # Extract cell type from filename
  cell_type <- str_extract(basename(file), "(?<=deseq2_results_).*(?=_OUD_vs_Normal)")
  cat("\nAnalyzing:", cell_type, "\n")

  # Read DESeq2 results
  deseq_results <- read_csv(file)

  # Filter for significant DEGs
  sig_genes <- deseq_results %>% filter(padj < 0.05)
  cat("Significant DEGs:", nrow(sig_genes), "\n")

  if (nrow(sig_genes) < 3) {
    cat("No significant DEGs found (less than 3), skipping...\n")
    next
  }

  # Save significant genes for reference
  sig_genes_file <- file.path(ora_dir, paste0(cell_type, "_significant_genes.csv"))
  write.csv(sig_genes, sig_genes_file, row.names = FALSE)
  cat("Saved significant gene list to:", sig_genes_file, "\n")

  # Convert gene IDs to Entrez IDs
  entrez_mapping <- convert_to_entrez(sig_genes$gene_id, sig_genes$gene_symbol)

  if (is.null(entrez_mapping) || nrow(entrez_mapping) < 3) {
    cat("Failed to convert enough gene IDs to Entrez IDs, skipping...\n")
    next
  }

  # Get background universe
  universe_mapping <- convert_to_entrez(deseq_results$gene_id, deseq_results$gene_symbol)
  if (is.null(universe_mapping) || nrow(universe_mapping) < 10) {
    cat("WARNING: Failed to convert enough universe genes, using default background\n")
    universe <- NULL
  } else {
    universe <- unique(universe_mapping$entrezid)
    cat("Using", length(universe), "genes as background universe\n")
  }

  # Run ORA analysis
  gene_list <- unique(entrez_mapping$entrezid)
  cat("Running ORA with", length(gene_list), "genes\n")
  output_prefix <- paste0(gsub(" ", "_", cell_type))
  ora_results <- run_ora_analysis(gene_list, universe, output_prefix)

  # Create enrichment map if multiple enrichment results are available
  if (length(ora_results) >= 2) {
    tryCatch({
      # Get names of result categories that have data
      result_names <- names(ora_results)

      # Create the map
      pdf(file.path(ora_dir, paste0(output_prefix, "_enrichment_map.pdf")), width = 12, height = 10)
      print(emapplot(ora_results, showCategory = 10))
      dev.off()

      # Create cnetplot for the first enrichment result
      pdf(file.path(ora_dir, paste0(output_prefix, "_cnetplot.pdf")), width = 12, height = 10)
      print(cnetplot(ora_results[[result_names[1]]], categorySize="pvalue", foldChange=NULL))
      dev.off()
    }, error = function(e) {
      message("Error creating enrichment visualizations: ", e$message)
    })
  }
}

cat("\nORA analysis complete. Results saved in:", ora_dir, "\n")