# Analysis Outputs_0 - Gene Set Enrichment Analysis (GSEA) of DESeq2 results
# This script performs pathway enrichment analysis using GSEA methodology

# Load required libraries
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(DOSE)
library(dplyr)
library(readr)
library(ggplot2)
library(stringr)

# Set base directory and create output directory
base_dir <- "/"
deseq2_dir <- file.path(base_dir, "GSE207128/DESeq2_outputs")
gsea_dir <- file.path(base_dir, "GSE207128/GSEA_results")
dir.create(gsea_dir, showWarnings = FALSE, recursive = TRUE)

# Find all DESeq2 result files
deseq2_files <- list.files(deseq2_dir, pattern = "deseq2_results_.*\\.csv$", full.names = TRUE)
cat("Found", length(deseq2_files), "DESeq2 result files\n")

# Function to convert mouse gene IDs to Entrez IDs
convert_to_entrez <- function(gene_ids, gene_symbols, fold_changes) {
  # Create mapping result dataframe
  result <- data.frame(
    original_id = gene_ids,
    gene_symbol = gene_symbols,
    fold_change = fold_changes,
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

# Function to prepare gene ranking for GSEA
prepare_gene_list <- function(entrez_mapping) {
  # Get unique entrez IDs (some might be duplicated)
  gene_list <- entrez_mapping %>%
    group_by(entrezid) %>%
    summarize(fold_change = mean(fold_change)) %>%
    arrange(desc(fold_change))

  # Create named vector for GSEA
  gene_ranks <- gene_list$fold_change
  names(gene_ranks) <- gene_list$entrezid

  return(gene_ranks)
}

# Function to perform GSEA analysis
run_gsea_analysis <- function(gene_list, output_prefix) {
  results <- list()
  min_genes <- 15  # GSEA typically needs more genes than ORA

  if (length(gene_list) < min_genes) {
    cat("Too few genes for GSEA analysis\n")
    return(results)
  }

  # GO Biological Process enrichment
  tryCatch({
    go_bp <- gseGO(geneList = gene_list,
                   OrgDb = org.Mm.eg.db,
                   ont = "BP",
                   minGSSize = 10,
                   maxGSSize = 500,
                   pvalueCutoff = 0.05,
                   verbose = FALSE,
                   keyType = "ENTREZID")

    if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
      results[["GO_BP"]] <- go_bp
      pdf(file.path(gsea_dir, paste0(output_prefix, "_GO_BP_gsea.pdf")), width = 10, height = 8)
      print(dotplot(go_bp, showCategory = 20, title = "GO Biological Process"))
      dev.off()

      # GSEA plot for top pathways
      if(nrow(go_bp@result) >= 1) {
        top_pathways <- head(go_bp@result$ID, 5)
        for(pathway in top_pathways) {
          tryCatch({
            pdf(file.path(gsea_dir, paste0(output_prefix, "_GO_BP_", pathway, "_gseaplot.pdf")), width = 10, height = 6)
            print(gseaplot2(go_bp, geneSetID = pathway, title = go_bp@result$Description[go_bp@result$ID == pathway]))
            dev.off()
          }, error = function(e) {
            message("Error plotting GSEA for pathway ", pathway, ": ", e$message)
          })
        }
      }

      # Save results table
      write.csv(go_bp@result, file.path(gsea_dir, paste0(output_prefix, "_GO_BP_gsea.csv")), row.names = FALSE)
    }
  }, error = function(e) {
    message("Error in GO BP GSEA analysis: ", e$message)
  })

  # GO Molecular Function enrichment
  tryCatch({
    go_mf <- gseGO(geneList = gene_list,
                   OrgDb = org.Mm.eg.db,
                   ont = "MF",
                   minGSSize = 10,
                   maxGSSize = 500,
                   pvalueCutoff = 0.05,
                   verbose = FALSE,
                   keyType = "ENTREZID")

    if (!is.null(go_mf) && nrow(go_mf@result) > 0) {
      results[["GO_MF"]] <- go_mf
      pdf(file.path(gsea_dir, paste0(output_prefix, "_GO_MF_gsea.pdf")), width = 10, height = 8)
      print(dotplot(go_mf, showCategory = 20, title = "GO Molecular Function"))
      dev.off()

      # Save results table
      write.csv(go_mf@result, file.path(gsea_dir, paste0(output_prefix, "_GO_MF_gsea.csv")), row.names = FALSE)
    }
  }, error = function(e) {
    message("Error in GO MF GSEA analysis: ", e$message)
  })

  # GO Cellular Component enrichment
  tryCatch({
    go_cc <- gseGO(geneList = gene_list,
                   OrgDb = org.Mm.eg.db,
                   ont = "CC",
                   minGSSize = 10,
                   maxGSSize = 500,
                   pvalueCutoff = 0.05,
                   verbose = FALSE,
                   keyType = "ENTREZID")

    if (!is.null(go_cc) && nrow(go_cc@result) > 0) {
      results[["GO_CC"]] <- go_cc
      pdf(file.path(gsea_dir, paste0(output_prefix, "_GO_CC_gsea.pdf")), width = 10, height = 8)
      print(dotplot(go_cc, showCategory = 20, title = "GO Cellular Component"))
      dev.off()

      # Save results table
      write.csv(go_cc@result, file.path(gsea_dir, paste0(output_prefix, "_GO_CC_gsea.csv")), row.names = FALSE)
    }
  }, error = function(e) {
    message("Error in GO CC GSEA analysis: ", e$message)
  })

  # KEGG pathway enrichment
  tryCatch({
    kegg <- gseKEGG(geneList = gene_list,
                    organism = "mmu",
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 0.05,
                    verbose = FALSE)

    if (!is.null(kegg) && nrow(kegg@result) > 0) {
      results[["KEGG"]] <- kegg
      pdf(file.path(gsea_dir, paste0(output_prefix, "_KEGG_gsea.pdf")), width = 10, height = 8)
      print(dotplot(kegg, showCategory = 20, title = "KEGG Pathways"))
      dev.off()

      # Save results table
      write.csv(kegg@result, file.path(gsea_dir, paste0(output_prefix, "_KEGG_gsea.csv")), row.names = FALSE)
    }
  }, error = function(e) {
    message("Error in KEGG GSEA analysis: ", e$message)
  })

  return(results)
}

# Process each DESeq2 result file
for (file in deseq2_files) {
  # Extract cell type from filename
  cell_type <- str_extract(basename(file), "(?<=deseq2_results_).*(?=_OUD_vs_Normal)")
  cat("\nAnalyzing:", cell_type, "\n")

  # Read DESeq2 results
  deseq_results <- read_csv(file)

  # Filter out low expression genes for better ranking
  filtered_results <- deseq_results %>%
    filter(!is.na(log2FoldChange)) %>%
    arrange(desc(abs(log2FoldChange)))

  cat("Total genes for GSEA:", nrow(filtered_results), "\n")

  if (nrow(filtered_results) < 100) {
    cat("Too few genes for meaningful GSEA, skipping...\n")
    next
  }

  # Convert gene IDs to Entrez IDs
  entrez_mapping <- convert_to_entrez(
    filtered_results$gene_id,
    filtered_results$gene_symbol,
    filtered_results$log2FoldChange
  )

  if (is.null(entrez_mapping) || nrow(entrez_mapping) < 50) {
    cat("Failed to convert enough gene IDs to Entrez IDs, skipping...\n")
    next
  }

  # Prepare gene list for GSEA
  gene_list <- prepare_gene_list(entrez_mapping)
  cat("Prepared ranked gene list with", length(gene_list), "genes\n")

  # Run GSEA analysis
  output_prefix <- paste0(gsub(" ", "_", cell_type))
  gsea_results <- run_gsea_analysis(gene_list, output_prefix)

  # Create enrichment map if multiple enrichment results are available
  if (length(gsea_results) >= 2) {
    tryCatch({
      pdf(file.path(gsea_dir, paste0(output_prefix, "_gsea_enrichment_map.pdf")), width = 12, height = 10)
      print(emapplot(gsea_results, showCategory = 10))
      dev.off()

      # Save a combined top results summary
      all_results <- data.frame()
      for (category in names(gsea_results)) {
        if (nrow(gsea_results[[category]]@result) > 0) {
          top_results <- head(gsea_results[[category]]@result, 10) %>%
            mutate(Category = category)
          all_results <- rbind(all_results, top_results)
        }
      }

      if (nrow(all_results) > 0) {
        write.csv(all_results, file.path(gsea_dir, paste0(output_prefix, "_top_gsea_results.csv")), row.names = FALSE)
      }
    }, error = function(e) {
      message("Error creating enrichment visualizations: ", e$message)
    })
  }
}

cat("\nGSEA analysis complete. Results saved in:", gsea_dir, "\n")