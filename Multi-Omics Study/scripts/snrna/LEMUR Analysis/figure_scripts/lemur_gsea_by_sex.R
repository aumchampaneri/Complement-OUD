# Complement-OUD/Multi-Omics Study/scripts/snrna/LEMUR Analysis/figure_scripts/lemur_gsea_by_brain_area.R
# Robust LEMUR + GSEA/ORA analysis by brain area
# Defensive against segfaults, malformed input, and package errors

library(zellkonverter)
library(SingleCellExperiment)
library(scran)
library(scater)
library(lemur)
library(tidyverse)
library(ggplot2)
library(viridis)
library(patchwork)
library(Matrix)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
# Optional: for full crash isolation
# library(callr)

# Optional: CLI args for input/output/species
args <- commandArgs(trailingOnly = TRUE)
BASE_DIR <- "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study"
INPUT_H5AD <- ifelse(length(args) > 0, args[1], file.path(BASE_DIR, "data/raw/snrna/GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5ad"))
OUTPUT_DIR <- ifelse(length(args) > 1, args[2], file.path(BASE_DIR, "scripts/snrna/LEMUR Analysis/outputs/panelA_by_brain_area"))
SPECIES <- ifelse(length(args) > 2, args[3], "Homo sapiens")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

PARAMS <- list(
  min_genes = 200,
  max_genes = 6000,
  max_mito_pct = 20,
  min_counts = 500,
  max_counts = 50000,
  n_hvg = 4000,
  n_embedding = 20,
  design_formula = "~ condition",
  batch_key = "donor_id",
  condition_key = "condition",
  fdr_threshold = 0.05,
  effect_size_threshold = 0.1,
  figure_dpi = 300,
  figure_width = 10,
  figure_height = 8
)

cat("📊 Loading H5AD data...\n")
sce <- zellkonverter::readH5AD(INPUT_H5AD, use_hdf5 = TRUE)
cat("✅ Data loaded: ", nrow(sce), "genes ×", ncol(sce), "cells\n")

meta_data <- colData(sce)
colnames(meta_data)[colnames(meta_data) == "ID"] <- "donor_id"
colnames(meta_data)[colnames(meta_data) == "Dx_OUD"] <- "condition"
if (!"region" %in% colnames(meta_data)) {
  possible_region_cols <- c("area", "brain_area", "brain_region", "region_name", "Region", "Area")
  found_region <- possible_region_cols[possible_region_cols %in% colnames(meta_data)]
  if (length(found_region) > 0) {
    message("Renaming column '", found_region[1], "' to 'region' for analysis.")
    colnames(meta_data)[colnames(meta_data) == found_region[1]] <- "region"
  } else {
    stop(
      "❌ Missing required column: region\n",
      "Available columns: ", paste(colnames(meta_data), collapse = ", "), "\n"
    )
  }
}
colData(sce) <- meta_data

cat("Available assays in SCE after loading:\n")
print(names(assays(sce)))
if (!"counts" %in% names(assays(sce))) {
  possible_counts <- c("X", "data", "logcounts")
  found_counts <- possible_counts[possible_counts %in% names(assays(sce))]
  if (length(found_counts) > 0) {
    message("Assigning assay '", found_counts[1], "' to 'counts' for downstream analysis.")
    assays(sce)$counts <- assays(sce)[[found_counts[1]]]
  } else {
    stop(
      "❌ No 'counts' assay found in SingleCellExperiment object.\n",
      "Available assays: ", paste(names(assays(sce)), collapse = ", "), "\n"
    )
  }
}

required_cols <- c("donor_id", "condition", "region")
missing_cols <- setdiff(required_cols, colnames(colData(sce)))
if (length(missing_cols) > 0) {
  stop("❌ Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# ---- Detect Sex column ----
sex_col <- NULL
possible_sex_cols <- c("Sex", "sex", "gender", "Gender")
for (col in possible_sex_cols) {
  if (col %in% colnames(colData(sce))) {
    sex_col <- col
    break
  }
}
if (is.null(sex_col)) {
  stop("❌ No Sex column found in metadata. Checked: ", paste(possible_sex_cols, collapse = ", "))
}
sexes <- unique(as.character(colData(sce)[[sex_col]]))
cat("🚻 Sex categories found:", paste(sexes, collapse = ", "), "\n")

# Save function (choose CSV or RDS)
save_results <- function(obj, path, use_rds = FALSE) {
  if (use_rds) {
    saveRDS(obj, file = sub("\\.csv$", ".rds", path))
  } else {
    write.csv(obj, file = path, row.names = FALSE)
  }
}

# Defensive, segmentation-safe enrichment runner
run_safe_gsea <- function(fun, ..., label = "") {
  cat("  [", label, "] Running enrichment...\n")
  result <- tryCatch({
    res <- fun(...)
    if (is.null(res) || nrow(as.data.frame(res)) == 0) {
      cat("  [", label, "] No significant results.\n")
      return(NULL)
    }
    cat("  [", label, "] Enrichment complete.\n")
    return(res)
  }, error = function(e) {
    cat("  [", label, "] Enrichment failed: ", e$message, "\n")
    return(NULL)
  }, finally = {
    gc()
  })
  result
}

cat("\n============================\n")
cat("Processing sex-stratified analysis (full dataset)\n")
cat("============================\n")

# QC and filtering on full dataset
is_mito <- grepl("^MT-|^mt-", rownames(sce))
sce <- addPerCellQCMetrics(sce, subsets = list(Mito = is_mito))
is_ribo <- grepl("^RP[SL]|^Rp[sl]", rownames(sce))
sce <- addPerCellQCMetrics(sce, subsets = list(Ribo = is_ribo))
colData(sce)$n_genes <- colSums(counts(sce) > 0)
colData(sce)$total_counts <- colSums(counts(sce))

cells_to_keep <- (
  sce$n_genes >= PARAMS$min_genes &
    sce$n_genes <= PARAMS$max_genes &
    sce$total_counts >= PARAMS$min_counts &
    sce$total_counts <= PARAMS$max_counts &
    sce$subsets_Mito_percent <= PARAMS$max_mito_pct
)
sce <- sce[, cells_to_keep]
genes_to_keep <- rowSums(counts(sce) > 0) >= 10
sce <- sce[genes_to_keep, ]
cat("After filtering:", nrow(sce), "genes ×", ncol(sce), "cells\n")

clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters = clusters)
sce <- logNormCounts(sce)
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, n = PARAMS$n_hvg)
sce_hvg <- sce[hvg, ]

colnames(colData(sce_hvg)) <- make.names(colnames(colData(sce_hvg)), unique = TRUE)
sce_hvg$condition <- droplevels(as.factor(sce_hvg$condition))
sce_hvg$sex <- droplevels(as.factor(sce_hvg[[sex_col]]))

# Always use ~ sex + condition as design formula
if (length(unique(sce_hvg$condition)) < 2 || length(unique(sce_hvg$sex)) < 2) {
  cat("⚠️ Skipping analysis: only one condition or sex present after filtering\n")
} else {
  design_formula <- ~ sex + condition
  cat("Using design formula:", deparse(design_formula), "\n")
  fit <- lemur(sce_hvg,
    design = design_formula,
    n_embedding = PARAMS$n_embedding,
    verbose = TRUE
  )
  cat("LEMUR model fitted for sex-stratified analysis\n")

  de_res <- test_de(fit, contrast = cond(condition = "OUD") - cond(condition = "None"))
  de_results_df <- data.frame(
    gene = rownames(de_res),
    effect_size = rowData(de_res)$effect_size,
    pval = rowData(de_res)$pval,
    adj_pval = p.adjust(rowData(de_res)$pval, method = "BH"),
    stringsAsFactors = FALSE
  )
  sig_results <- de_results_df[
    de_results_df$adj_pval < PARAMS$fdr_threshold &
    abs(de_results_df$effect_size) > PARAMS$effect_size_threshold,
  ]
  cat("Significant DE genes:", nrow(sig_results), "\n")

  output_dir_sex <- file.path(OUTPUT_DIR, "sex_condition")
  dir.create(output_dir_sex, recursive = TRUE, showWarnings = FALSE)
  save_results(de_results_df, file.path(output_dir_sex, "lemur_de_results_sex_condition.csv"))
  save_results(sig_results, file.path(output_dir_sex, "lemur_sig_de_genes_sex_condition.csv"))

  # Prepare gene lists and validate
  gene_list <- de_results_df$effect_size
  names(gene_list) <- de_results_df$gene
  gene_list <- sort(gene_list, decreasing = TRUE)
  gene_df <- bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  gene_df <- gene_df[!duplicated(gene_df$SYMBOL) & !duplicated(gene_df$ENTREZID), ] # 1:1 mapping
  gene_list_entrez <- gene_list[gene_df$SYMBOL]
  names(gene_list_entrez) <- gene_df$ENTREZID
  gene_list_entrez <- gene_list_entrez[!is.na(names(gene_list_entrez))]
  sig_gene_symbols <- sig_results$gene
  sig_gene_entrez <- gene_df$ENTREZID[gene_df$SYMBOL %in% sig_gene_symbols]

  # Pre-validation for GSEA input
  if (length(gene_list_entrez) < 10 || is.null(names(gene_list_entrez)) ||
    any(is.na(names(gene_list_entrez))) || !is.numeric(gene_list_entrez)) {
    warning("⚠️ Skipping GSEA/ORA: Invalid gene_list_entrez (too few genes or bad formatting).")
  } else if (length(sig_gene_entrez) < 5) {
    warning("⚠️ Skipping ORA: Too few significant genes.")
  } else {
    # Hallmark
    if (requireNamespace("msigdbr", quietly = TRUE)) {
      library(msigdbr)
      hallmark <- msigdbr(species = SPECIES, category = "H")
      hallmark_df <- hallmark[, c("gs_name", "entrez_gene")]
      gsea_hallmark <- run_safe_gsea(GSEA, geneList = gene_list_entrez, TERM2GENE = hallmark_df, pvalueCutoff = 0.05, verbose = FALSE, label = "Hallmark GSEA")
      save_results(as.data.frame(gsea_hallmark), file.path(output_dir_sex, "hallmark_gsea_results_sex_condition.csv"))
      ora_hallmark <- run_safe_gsea(enricher, gene = sig_gene_entrez, TERM2GENE = hallmark_df, pvalueCutoff = 0.05, label = "Hallmark ORA")
      save_results(as.data.frame(ora_hallmark), file.path(output_dir_sex, "hallmark_ora_results_sex_condition.csv"))
    }

    # GO
    for (ont in c("BP", "CC", "MF")) {
    gsea_go <- run_safe_gsea(gseGO, geneList = gene_list_entrez, OrgDb = org.Hs.eg.db, ont = ont, keyType = "ENTREZID", minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE, label = paste0("GO_", ont, " GSEA"))
    save_results(as.data.frame(gsea_go), file.path(OUTPUT_DIR, paste0("GO_", ont, "_gsea_results_", area_tag, ".csv")))
    ora_go <- run_safe_gsea(enrichGO, gene = sig_gene_entrez, OrgDb = org.Hs.eg.db, ont = ont, keyType = "ENTREZID", pvalueCutoff = 0.05, label = paste0("GO_", ont, " ORA"))
    save_results(as.data.frame(ora_go), file.path(OUTPUT_DIR, paste0("GO_", ont, "_ora_results_", area_tag, ".csv")))
  }

  # KEGG
  if (requireNamespace("KEGGREST", quietly = TRUE)) {
    kegg_org <- "hsa"
    gsea_kegg <- run_safe_gsea(gseKEGG, geneList = gene_list_entrez, organism = kegg_org, pvalueCutoff = 0.05, verbose = FALSE, label = "KEGG GSEA")
    save_results(as.data.frame(gsea_kegg), file.path(OUTPUT_DIR, paste0("KEGG_gsea_results_", area_tag, ".csv")))
    ora_kegg <- run_safe_gsea(enrichKEGG, gene = sig_gene_entrez, organism = kegg_org, pvalueCutoff = 0.05, label = "KEGG ORA")
    save_results(as.data.frame(ora_kegg), file.path(OUTPUT_DIR, paste0("KEGG_ora_results_", area_tag, ".csv")))
  }

  # Reactome
  if (requireNamespace("ReactomePA", quietly = TRUE)) {
    library(ReactomePA)
    gsea_reactome <- run_safe_gsea(gsePathway, geneList = gene_list_entrez, organism = "human", pvalueCutoff = 0.05, verbose = FALSE, label = "Reactome GSEA")
    save_results(as.data.frame(gsea_reactome), file.path(OUTPUT_DIR, paste0("Reactome_gsea_results_", area_tag, ".csv")))
    ora_reactome <- run_safe_gsea(enrichPathway, gene = sig_gene_entrez, organism = "human", pvalueCutoff = 0.05, readable = TRUE, label = "Reactome ORA")
    save_results(as.data.frame(ora_reactome), file.path(OUTPUT_DIR, paste0("Reactome_ora_results_", area_tag, ".csv")))
  }

  # TFs
  if (requireNamespace("msigdbr", quietly = TRUE)) {
    tfsets <- msigdbr(species = SPECIES, category = "C3")
    tfsets <- tfsets[grepl("TFT", tfsets$gs_subcat), ]
    tf_df <- tfsets[, c("gs_name", "entrez_gene")]
    gsea_tf <- run_safe_gsea(GSEA, geneList = gene_list_entrez, TERM2GENE = tf_df, pvalueCutoff = 0.05, verbose = FALSE, label = "TF GSEA")
    save_results(as.data.frame(gsea_tf), file.path(OUTPUT_DIR, paste0("TF_gsea_results_", area_tag, ".csv")))
    ora_tf <- run_safe_gsea(enricher, gene = sig_gene_entrez, TERM2GENE = tf_df, pvalueCutoff = 0.05, label = "TF ORA")
    save_results(as.data.frame(ora_tf), file.path(OUTPUT_DIR, paste0("TF_ora_results_", area_tag, ".csv")))
  }

  # Explicit garbage collection
  gc()
  cat("Completed brain area:", area, "\n")
}

cat("\n✅ Panel A analysis complete. Results saved to:", OUTPUT_DIR, "\n")

# ---- Optional: Unit test block for GSEA segfault debugging ----
if (interactive()) {
  cat("\n[Unit Test] Testing gseGO on small gene set...\n")
  test_gene_list <- head(gene_list_entrez, 100)
  test_result <- tryCatch(
    {
      gseGO(
        geneList = test_gene_list, OrgDb = org.Hs.eg.db, ont = "BP",
        keyType = "ENTREZID", verbose = FALSE
      )
    },
    error = function(e) e$message
  )
  print(test_result)
}

# ---- Optional: callr isolation for GSEA/ORA (uncomment to use) ----
# run_callr_gsea <- function(fun, ..., label = "") {
#   cat("  [", label, "] Running enrichment in isolated process...\n")
#   result <- tryCatch({
#     res <- callr::r(function(fun, ...) fun(...), args = list(fun, ...))
#     if (is.null(res) || nrow(as.data.frame(res)) == 0) {
#       cat("  [", label, "] No significant results.\n")
#       return(NULL)
#     }
#     cat("  [", label, "] Enrichment complete.\n")
#     return(res)
#   }, error = function(e) {
#     cat("  [", label, "] Enrichment failed: ", e$message, "\n")
#     return(NULL)
#   }, finally = {
#     gc()
#   })
#   result
# }
