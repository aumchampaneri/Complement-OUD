# R

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(org.Mm.eg.db)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(pheatmap)
  library(RColorBrewer)
})

project_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE289002"
de_results_dir <- file.path(project_dir, "DE_results")
enrichment_dir <- file.path(project_dir, "Enrichment_results")
dir.create(enrichment_dir, showWarnings = FALSE, recursive = TRUE)

#--- Helper: Load gene mapping and create lookup tables
load_gene_mapping <- function(mapping_file) {
  message("[DEBUG] Loading gene mapping from: ", mapping_file)
  gene_mapping <- tryCatch(read.csv(mapping_file, stringsAsFactors = FALSE), error = function(e) stop("[ERROR] Could not read gene mapping file: ", e$message))
  message("[DEBUG] Mapping file loaded. Rows: ", nrow(gene_mapping))
  ensembl_to_symbol <- setNames(gene_mapping$gene_symbol, gene_mapping$ensembl_id)
  ensembl_ids <- gene_mapping$ensembl_id
  ensembl2entrez <- bitr(ensembl_ids, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  message("[DEBUG] bitr mapping done. Mapped: ", nrow(ensembl2entrez), " of ", length(ensembl_ids))
  ensembl_to_entrez <- setNames(ensembl2entrez$ENTREZID, ensembl2entrez$ENSEMBL)
  message("[DEBUG] Loaded gene mappings for ", length(ensembl_to_entrez), " genes")
  list(symbol = ensembl_to_symbol, entrez = ensembl_to_entrez)
}

#--- Helper: Robust dotplot wrapper
safe_dotplot <- function(enrich_obj, out_path, title) {
  res <- enrich_obj@result
  # Check for at least 2 rows and at least 2 unique non-NA Description values
  if (is.null(res) || nrow(res) < 2 || length(unique(na.omit(res$Description))) < 2) {
    message("[DEBUG] Not enough terms for dotplot: ", title)
    return()
  }
  # Check that the plotting variable (GeneRatio, Count, or p.adjust) is not constant or all NA
  plot_cols <- c("GeneRatio", "Count", "p.adjust")
  plot_col <- plot_cols[plot_cols %in% colnames(res)][1]
  if (!is.null(plot_col)) {
    vals <- res[[plot_col]]
    if (all(is.na(vals)) || length(unique(na.omit(vals))) < 2) {
      message("[DEBUG] Not enough variation in ", plot_col, " for dotplot: ", title)
      return()
    }
  }
  # Try-catch for dotplot to avoid fatal errors
  tryCatch({
    png(out_path, width = 1000, height = 800, res = 100)
    print(dotplot(enrich_obj, showCategory = 20, title = title))
    dev.off()
  }, error = function(e) {
    message("[DEBUG] dotplot failed for ", title, ": ", e$message)
    try(dev.off(), silent = TRUE)
  })
}

#--- Helper: Perform enrichment for a single contrast
perform_enrichment <- function(de_file, out_dir, contrast, region, ensembl_to_entrez) {
  message("[DEBUG] Reading DE file: ", de_file)
  if (!file.exists(de_file)) {
    message("[ERROR] DE file does not exist: ", de_file)
    return(NULL)
  }
  de_results <- tryCatch(read.csv(de_file, row.names = 1, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(de_results) || nrow(de_results) == 0) {
    message("[ERROR] Failed to read or empty DE file: ", de_file)
    return(NULL)
  }
  de_results$ensembl_id <- rownames(de_results)
  required_cols <- c("FDR", "ensembl_id", "logFC")
  if (!all(required_cols %in% colnames(de_results))) {
    message("[ERROR] Required columns missing in ", de_file, ": ", paste(setdiff(required_cols, colnames(de_results)), collapse = ", "))
    return(NULL)
  }
  sig_genes <- de_results[de_results$FDR < 0.05 & !is.na(de_results$FDR), ]
  message("[DEBUG] Significant genes (FDR < 0.05): ", nrow(sig_genes))
  if (nrow(sig_genes) == 0) {
    message("[DEBUG] No significant genes for ", contrast, " in ", region)
    return(NULL)
  }
  sig_entrez <- ensembl_to_entrez[sig_genes$ensembl_id]
  sig_entrez <- sig_entrez[!is.na(sig_entrez)]
  message("[DEBUG] Non-NA Entrez IDs: ", length(sig_entrez))
  if (length(sig_entrez) == 0) {
    message("[DEBUG] No Entrez IDs mapped for significant genes in ", contrast, " in ", region)
    return(NULL)
  }
  all_genes <- de_results$ensembl_id
  gene_list <- de_results$logFC
  names(gene_list) <- ensembl_to_entrez[all_genes]
  gene_list <- gene_list[!is.na(names(gene_list))]
  gene_list <- sort(gene_list, decreasing = TRUE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  # GO BP
  message("[DEBUG] Running enrichGO BP...")
  go_bp <- tryCatch(
    enrichGO(gene = sig_entrez, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2),
    error = function(e) NULL
  )
  if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
    message("[DEBUG] GO BP enrichment found ", nrow(go_bp@result), " terms")
    write.csv(go_bp@result, file.path(out_dir, "GO_BP_results.csv"), row.names = FALSE)
    safe_dotplot(go_bp, file.path(out_dir, "GO_BP_dotplot.png"), paste0("GO BP - ", region, " - ", contrast))
  } else {
    message("[DEBUG] No GO BP enrichment results for ", contrast, " in ", region)
  }
  # GO MF
  message("[DEBUG] Running enrichGO MF...")
  go_mf <- tryCatch(
    enrichGO(gene = sig_entrez, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2),
    error = function(e) NULL
  )
  if (!is.null(go_mf) && nrow(go_mf@result) > 0) {
    message("[DEBUG] GO MF enrichment found ", nrow(go_mf@result), " terms")
    write.csv(go_mf@result, file.path(out_dir, "GO_MF_results.csv"), row.names = FALSE)
    safe_dotplot(go_mf, file.path(out_dir, "GO_MF_dotplot.png"), paste0("GO MF - ", region, " - ", contrast))
  } else {
    message("[DEBUG] No GO MF enrichment results for ", contrast, " in ", region)
  }
  # GO CC
  message("[DEBUG] Running enrichGO CC...")
  go_cc <- tryCatch(
    enrichGO(gene = sig_entrez, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2),
    error = function(e) NULL
  )
  if (!is.null(go_cc) && nrow(go_cc@result) > 0) {
    message("[DEBUG] GO CC enrichment found ", nrow(go_cc@result), " terms")
    write.csv(go_cc@result, file.path(out_dir, "GO_CC_results.csv"), row.names = FALSE)
    safe_dotplot(go_cc, file.path(out_dir, "GO_CC_dotplot.png"), paste0("GO CC - ", region, " - ", contrast))
  } else {
    message("[DEBUG] No GO CC enrichment results for ", contrast, " in ", region)
  }
  # KEGG
  message("[DEBUG] Running enrichKEGG...")
  kegg <- tryCatch(
    enrichKEGG(gene = sig_entrez, organism = "mmu", keyType = "ncbi-geneid", pAdjustMethod = "BH", pvalueCutoff = 0.05),
    error = function(e) NULL
  )
  if (!is.null(kegg) && nrow(kegg@result) > 0) {
    message("[DEBUG] KEGG enrichment found ", nrow(kegg@result), " terms")
    write.csv(kegg@result, file.path(out_dir, "KEGG_results.csv"), row.names = FALSE)
    safe_dotplot(kegg, file.path(out_dir, "KEGG_dotplot.png"), paste0("KEGG - ", region, " - ", contrast))
  } else {
    message("[DEBUG] No KEGG enrichment results for ", contrast, " in ", region)
  }
  invisible(list(go_bp = go_bp, go_mf = go_mf, go_cc = go_cc, kegg = kegg))
}

#--- Main: Run enrichment for all regions/contrasts
main <- function() {
  mapping_file <- file.path(project_dir, "ensembl_to_symbol.csv")
  if (!file.exists(mapping_file)) stop("[ERROR] Gene mapping file not found.")
  mapping <- load_gene_mapping(mapping_file)
  regions <- list.dirs(de_results_dir, full.names = FALSE, recursive = FALSE)
  regions <- regions[regions != ""]
  message("[DEBUG] Found region folders: ", paste(regions, collapse = ", "))
  for (region_folder in regions) {
    region_path <- file.path(de_results_dir, region_folder)
    if (!dir.exists(region_path)) next
    region <- region_folder
    region_output <- file.path(enrichment_dir, region)
    dir.create(region_output, showWarnings = FALSE, recursive = TRUE)
    de_files <- list.files(region_path, pattern = "\\.csv$", full.names = TRUE)
    message("[DEBUG] Region: ", region, " DE files found: ", paste(basename(de_files), collapse = ", "))
    for (de_file in de_files) {
      contrast <- gsub("\\.csv$", "", basename(de_file))
      contrast_dir <- file.path(region_output, paste0(region, "_", contrast))
      message("[DEBUG] Processing contrast: ", contrast, " in region: ", region)
      perform_enrichment(de_file, contrast_dir, contrast, region, mapping$entrez)
    }
  }
  message("Enrichment analysis complete. Results saved in: ", enrichment_dir)
}

main()