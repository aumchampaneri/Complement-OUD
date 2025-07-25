Complement-OUD/Multi-Omics Study/scripts/snrna/LEMUR Pipeline/de_enrichment/de_enrichment_workflow.R
```

```R
# DE & Enrichment Workflow for LEMUR snRNA-seq Analysis
# ----------------------------------------------------
# Loads LEMUR fits, extracts pseudobulk per cell type, runs LEMUR test_de for DE,
# and performs ORA, GSEA, RRHO, SCENIC, WGCNA, and GWAS integration.
# Modular and ready for expansion.

# ==============================
# 0. Setup & Libraries
# ==============================
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(clusterProfiler)
  library(ReactomePA)
  library(msigdbr)
  library(rrho)
  # library(SCENIC) # Uncomment if SCENIC is installed/configured
  # library(WGCNA)  # Uncomment if WGCNA is installed/configured
  library(Matrix)
  library(readr)
  library(lemur)
  library(org.Hs.eg.db)
})

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path(dirname(getwd()), "outputs", paste0("de_enrichment_", timestamp))
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

log_message <- function(message) {
  cat("[", format(Sys.time(), "%H:%M:%S"), "] ", message, "\n")
}

# ==============================
# 1. Load LEMUR Fits & DE Results
# ==============================
lemur_obj_dir <- file.path("..", "lemur_analysis", "objects")
sex_levels <- c("M", "F") # Update as needed
region_levels <- c("Caudate", "Putamen") # Update as needed

lemur_fits <- list()
lemur_de_results <- list()

for (sex in sex_levels) {
  fit_path <- file.path(lemur_obj_dir, paste0("lemur_fit_sex_", sex, "_", timestamp, ".rds"))
  de_path  <- file.path(lemur_obj_dir, paste0("lemur_de_sex_", sex, "_", timestamp, ".rds"))
  if (file.exists(fit_path) && file.exists(de_path)) {
    lemur_fits[[paste0("sex_", sex)]] <- readRDS(fit_path)
    lemur_de_results[[paste0("sex_", sex)]] <- readRDS(de_path)
    log_message(paste("Loaded LEMUR fit and DE for sex:", sex))
  }
}
for (region in region_levels) {
  fit_path <- file.path(lemur_obj_dir, paste0("lemur_fit_region_", region, "_", timestamp, ".rds"))
  de_path  <- file.path(lemur_obj_dir, paste0("lemur_de_region_", region, "_", timestamp, ".rds"))
  if (file.exists(fit_path) && file.exists(de_path)) {
    lemur_fits[[paste0("region_", region)]] <- readRDS(fit_path)
    lemur_de_results[[paste0("region_", region)]] <- readRDS(de_path)
    log_message(paste("Loaded LEMUR fit and DE for region:", region))
  }
}

# ==============================
# 2. Extract Pseudobulk per Cell Type
# ==============================
# (Assume SCE object is available as 'sce' or load from checkpoint)
sce_path <- file.path(lemur_obj_dir, "sce_filtered.rds")
if (file.exists(sce_path)) {
  sce <- readRDS(sce_path)
} else {
  stop("Filtered SCE object not found. Please check pipeline outputs.")
}

cell_types <- unique(colData(sce)$celltype)
pseudobulk_list <- list()
for (ct in cell_types) {
  cells_ct <- which(colData(sce)$celltype == ct)
  # Aggregate counts by donor/sample for pseudobulk
  donors <- colData(sce)$donor_id[cells_ct]
  pb_mat <- t(sapply(unique(donors), function(d) {
    if (sum(donors == d) > 1) {
      Matrix::rowSums(counts(sce)[, cells_ct[donors == d], drop = FALSE])
    } else {
      counts(sce)[, cells_ct[donors == d], drop = FALSE]
    }
  }))
  colnames(pb_mat) <- unique(donors)
  pseudobulk_list[[ct]] <- pb_mat
  log_message(paste("Extracted pseudobulk for cell type:", ct))
}

# ==============================
# 3. Differential Expression (LEMUR test_de)
# ==============================
de_results_ct <- list()
for (name in names(lemur_fits)) {
  fit <- lemur_fits[[name]]
  de <- lemur_de_results[[name]]
  # Optionally, run additional contrasts or cell-type specific DE
  # Example: test_de(fit, contrast = ...)
  de_results_ct[[name]] <- de
  # Save DE table
  saveRDS(de, file.path(output_dir, paste0("DE_", name, "_", timestamp, ".rds")))
  log_message(paste("Saved DE results for:", name))
}

# ==============================
# 4. Enrichment Analyses
# ==============================
enrichment_results <- list()
for (name in names(de_results_ct)) {
  de_table <- de_results_ct[[name]]
  sig_genes <- rownames(de_table)[de_table$adj_pval < 0.05]
  ranked_genes <- setNames(de_table$logFC, rownames(de_table))
  # ORA: GO, Reactome, MSigDB
  enrich_go <- tryCatch({
    enrichGO(gene = sig_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
  }, error = function(e) NULL)
  enrich_reactome <- tryCatch({
    enrichPathway(gene = sig_genes, organism = "human", pAdjustMethod = "BH", qvalueCutoff = 0.05)
  }, error = function(e) NULL)
  msigdb_sets <- msigdbr(species = "Homo sapiens", category = "H")
  enrich_msig <- tryCatch({
    enricher(gene = sig_genes, TERM2GENE = msigdb_sets[,c("gs_name", "gene_symbol")])
  }, error = function(e) NULL)
  enrichment_results[[paste0(name, "_GO")]] <- enrich_go
  enrichment_results[[paste0(name, "_Reactome")]] <- enrich_reactome
  enrichment_results[[paste0(name, "_MSigDB")]] <- enrich_msig
  saveRDS(enrich_go, file.path(output_dir, paste0("enrichment_", name, "_GO_", timestamp, ".rds")))
  saveRDS(enrich_reactome, file.path(output_dir, paste0("enrichment_", name, "_Reactome_", timestamp, ".rds")))
  saveRDS(enrich_msig, file.path(output_dir, paste0("enrichment_", name, "_MSigDB_", timestamp, ".rds")))
  log_message(paste("Saved enrichment results for:", name))
  # GSEA
  gsea_res <- tryCatch({
    GSEA(ranked_genes, TERM2GENE = msigdb_sets[,c("gs_name", "gene_symbol")])
  }, error = function(e) NULL)
  enrichment_results[[paste0(name, "_GSEA")]] <- gsea_res
  saveRDS(gsea_res, file.path(output_dir, paste0("gsea_", name, "_", timestamp, ".rds")))
  log_message(paste("Saved GSEA results for:", name))
}

# ==============================
# 5. RRHO Analysis
# ==============================
# Example: Compare ranked lists between regions or sexes
# (Requires two ranked lists, e.g., logFC from two DE tables)
if (length(de_results_ct) >= 2) {
  names_list <- names(de_results_ct)
  for (i in 1:(length(names_list)-1)) {
    for (j in (i+1):length(names_list)) {
      ranked1 <- setNames(de_results_ct[[names_list[i]]]$logFC, rownames(de_results_ct[[names_list[i]]]))
      ranked2 <- setNames(de_results_ct[[names_list[j]]]$logFC, rownames(de_results_ct[[names_list[j]]]))
      rrho_res <- tryCatch({
        RRHO(ranked1, ranked2, alternative = "enrichment")
      }, error = function(e) NULL)
      saveRDS(rrho_res, file.path(output_dir, paste0("rrho_", names_list[i], "_vs_", names_list[j], "_", timestamp, ".rds")))
      log_message(paste("Saved RRHO results for:", names_list[i], "vs", names_list[j]))
    }
  }
}

# ==============================
# 6. SCENIC & WGCNA (Network Analysis)
# ==============================
# Placeholder: Uncomment and configure as needed
# for (ct in cell_types) {
#   # SCENIC: Run on pseudobulk or single-cell expression for each cell type
#   # scenic_res <- runSCENIC(pseudobulk_list[[ct]], ...)
#   # saveRDS(scenic_res, file.path(output_dir, paste0("scenic_", ct, "_", timestamp, ".rds")))
#   # WGCNA: Run on pseudobulk or single-cell expression for each cell type
#   # wgcna_res <- blockwiseModules(pseudobulk_list[[ct]], ...)
#   # saveRDS(wgcna_res, file.path(output_dir, paste0("wgcna_", ct, "_", timestamp, ".rds")))
# }

# ==============================
# 7. GWAS Integration (LDSC/MAGMA)
# ==============================
# Placeholder: If significant DEGs, run GWAS enrichment
# if (length(sig_genes) > 0) {
#   # ldsc_res <- runLDSC(sig_genes, gwas_data, ...)
#   # saveRDS(ldsc_res, file.path(output_dir, paste0("ldsc_", name, "_", timestamp, ".rds")))
# }

# ==============================
# 8. Summary & Session Info
# ==============================
log_message("DE/Enrichment workflow complete.")
log_message(capture.output(sessionInfo()))
