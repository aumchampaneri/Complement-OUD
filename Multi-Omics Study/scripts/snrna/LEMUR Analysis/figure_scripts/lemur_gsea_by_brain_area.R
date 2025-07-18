# Script: lemur_gsea_by_brain_area.R
# Purpose: Run LEMUR and GSEA stratified by brain area for Figure Panel A (Pathway enrichment by brain area → OUD vs Control)
# Author: [Your Name]
# Date: [Today's Date]
# Description:
#   For each brain area, this script subsets the snRNA-seq data, runs LEMUR for OUD vs Control,
#   performs GSEA on the DE results, and saves area-specific outputs and plots.
#   This script is designed to be run after the main LEMUR pipeline and assumes the same environment.

# =============================================================================
# SETUP
# =============================================================================

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
# Add any other packages your GSEA pipeline requires

# ---- File paths ----
BASE_DIR <- "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study"
INPUT_H5AD <- file.path(BASE_DIR, "data/raw/snrna/GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5ad")
OUTPUT_DIR <- file.path(BASE_DIR, "scripts/snrna/LEMUR Analysis/outputs/panelA_by_brain_area")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Analysis parameters ----
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

# =============================================================================
# LOAD DATA
# =============================================================================

cat("📊 Loading H5AD data...\n")
sce <- zellkonverter::readH5AD(INPUT_H5AD, use_hdf5 = TRUE)
cat("✅ Data loaded: ", nrow(sce), "genes ×", ncol(sce), "cells\n")

# ---- Ensure metadata columns ----
meta_data <- colData(sce)
colnames(meta_data)[colnames(meta_data) == "ID"] <- "donor_id"
colnames(meta_data)[colnames(meta_data) == "Dx_OUD"] <- "condition"

# ---- Robust region column detection ----
if (!"region" %in% colnames(meta_data)) {
  # Try common alternatives
  possible_region_cols <- c("area", "brain_area", "brain_region", "region_name", "Region", "Area")
  found_region <- possible_region_cols[possible_region_cols %in% colnames(meta_data)]
  if (length(found_region) > 0) {
    message("Renaming column '", found_region[1], "' to 'region' for analysis.")
    colnames(meta_data)[colnames(meta_data) == found_region[1]] <- "region"
  } else {
    stop(
      "❌ Missing required column: region\n",
      "Available columns: ", paste(colnames(meta_data), collapse = ", "), "\n",
      "If your region column has a different name, please rename it to 'region' or add logic to detect it."
    )
  }
}
colData(sce) <- meta_data

# ---- Ensure 'counts' assay exists ----
cat("Available assays in SCE after loading:\n")
print(names(assays(sce)))
if (!"counts" %in% names(assays(sce))) {
  # Try common alternatives
  possible_counts <- c("X", "data", "logcounts")
  found_counts <- possible_counts[possible_counts %in% names(assays(sce))]
  if (length(found_counts) > 0) {
    message("Assigning assay '", found_counts[1], "' to 'counts' for downstream analysis.")
    assays(sce)$counts <- assays(sce)[[found_counts[1]]]
  } else {
    stop(
      "❌ No 'counts' assay found in SingleCellExperiment object.\n",
      "Available assays: ", paste(names(assays(sce)), collapse = ", "), "\n",
      "Please ensure your H5AD contains a raw count matrix or update the script to use the correct assay."
    )
  }
}

# ---- Check for required columns ----
required_cols <- c("donor_id", "condition", "region")
missing_cols <- setdiff(required_cols, colnames(colData(sce)))
if (length(missing_cols) > 0) {
  stop("❌ Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# =============================================================================
# STRATIFY BY BRAIN AREA
# =============================================================================

areas <- unique(as.character(sce$region))
cat("🧠 Brain areas found:", paste(areas, collapse = ", "), "\n")

for (area in areas) {
  cat("\n============================\n")
  cat("Processing brain area:", area, "\n")
  cat("============================\n")

  # Subset SCE for this brain area
  sce_area <- sce[, sce$region == area]
  cat("Cells in area:", ncol(sce_area), "\n")

  # Skip if too few cells
  if (ncol(sce_area) < 100) {
    cat("⚠️  Skipping area (too few cells)\n")
    next
  }

  # ---- Quality control ----
  is_mito <- grepl("^MT-|^mt-", rownames(sce_area))
  sce_area <- addPerCellQCMetrics(sce_area, subsets = list(Mito = is_mito))
  is_ribo <- grepl("^RP[SL]|^Rp[sl]", rownames(sce_area))
  sce_area <- addPerCellQCMetrics(sce_area, subsets = list(Ribo = is_ribo))
  colData(sce_area)$n_genes <- colSums(counts(sce_area) > 0)
  colData(sce_area)$total_counts <- colSums(counts(sce_area))

  cells_to_keep <- (
    sce_area$n_genes >= PARAMS$min_genes &
      sce_area$n_genes <= PARAMS$max_genes &
      sce_area$total_counts >= PARAMS$min_counts &
      sce_area$total_counts <= PARAMS$max_counts &
      sce_area$subsets_Mito_percent <= PARAMS$max_mito_pct
  )
  sce_area <- sce_area[, cells_to_keep]
  genes_to_keep <- rowSums(counts(sce_area) > 0) >= 10
  sce_area <- sce_area[genes_to_keep, ]
  cat("After filtering:", nrow(sce_area), "genes ×", ncol(sce_area), "cells\n")

  # ---- Normalization and HVG selection ----
  clusters <- quickCluster(sce_area)
  sce_area <- computeSumFactors(sce_area, clusters = clusters)
  sce_area <- logNormCounts(sce_area)
  dec <- modelGeneVar(sce_area)
  hvg <- getTopHVGs(dec, n = PARAMS$n_hvg)
  sce_area_hvg <- sce_area[hvg, ]

  # ---- LEMUR model ----
  colnames(colData(sce_area_hvg)) <- make.names(colnames(colData(sce_area_hvg)), unique = TRUE)
  sce_area_hvg$condition <- droplevels(as.factor(sce_area_hvg$condition))

  # Always use ~condition as design formula (no donor control)
  if (length(unique(sce_area_hvg$condition)) < 2) {
    cat("⚠️ Skipping area:", area, "- only one condition present after filtering\n")
    next
  }
  design_formula <- ~condition
  cat("Using design formula:", deparse(design_formula), "\n")

  # Print design matrix rank for debugging
  X <- model.matrix(design_formula, data = colData(sce_area_hvg))
  cat("✅ Design matrix rank:", Matrix::rankMatrix(X)[1], "of", ncol(X), "columns\n")

  fit <- lemur(sce_area_hvg,
    design = design_formula,
    n_embedding = PARAMS$n_embedding,
    verbose = TRUE
  )
  cat("LEMUR model fitted for", area, "\n")

  # Run LEMUR DE with default contrast (OUD vs None)
  de_res <- test_de(fit, contrast = cond(condition = "OUD") - cond(condition = "None"))
  # Check for empty DE results
  if (nrow(de_res) == 0 || length(rowData(de_res)$effect_size) == 0) {
    cat("⚠️ No DE results for area:", area, "- skipping downstream analysis\n")
    next
  }

  # Extract results directly from rowData
  de_results_df <- data.frame(
    gene = rownames(de_res),
    effect_size = rowData(de_res)$effect_size,
    pval = rowData(de_res)$pval,
    adj_pval = p.adjust(rowData(de_res)$pval, method = "BH"),
    stringsAsFactors = FALSE
  )

  # Apply filtering
  sig_results <- de_results_df[
    de_results_df$adj_pval < PARAMS$fdr_threshold &
      abs(de_results_df$effect_size) > PARAMS$effect_size_threshold,
  ]
  cat("Significant DE genes:", nrow(sig_results), "\n")

  # ---- Save DE results ----
  area_tag <- gsub("[^A-Za-z0-9]", "_", area)
  write.csv(de_results_df,
    file = file.path(OUTPUT_DIR, paste0("lemur_de_results_", area_tag, ".csv")),
    row.names = FALSE
  )
  write.csv(sig_results,
    file = file.path(OUTPUT_DIR, paste0("lemur_sig_de_genes_", area_tag, ".csv")),
    row.names = FALSE
  )

  # ---- Enhanced GSEA & ORA: Hallmark, GO, KEGG, Reactome, TFs ----
  # This section performs both GSEA and ORA for multiple gene set collections.
  # It is robust to missing packages and will skip unavailable analyses.
  # Results are saved as CSVs and PNG plots, labeled by brain area and gene set type.

  # Helper: Save and plot top pathways
  save_and_plot_gsea <- function(gsea_obj, prefix, area, area_tag, top_n = 10) {
    if (!is.null(gsea_obj) && nrow(as.data.frame(gsea_obj)) > 0) {
      gsea_df <- as.data.frame(gsea_obj)
      write.csv(gsea_df,
        file = file.path(OUTPUT_DIR, paste0(prefix, "_gsea_results_", area_tag, ".csv")),
        row.names = FALSE
      )
      # Plot top pathways
      top_gsea <- head(gsea_df, top_n)
      p <- ggplot(top_gsea, aes(x = reorder(Description, NES), y = NES, fill = p.adjust < 0.05)) +
        geom_col() +
        coord_flip() +
        labs(
          title = paste("Top", prefix, "GSEA Pathways -", area),
          x = "Pathway", y = "Normalized Enrichment Score (NES)"
        ) +
        scale_fill_manual(values = c("grey70", "tomato")) +
        theme_minimal()
      ggsave(
        filename = file.path(OUTPUT_DIR, paste0(prefix, "_gsea_top_", area_tag, ".png")),
        plot = p, width = 8, height = 5, dpi = PARAMS$figure_dpi
      )
      cat(prefix, "GSEA results and plot saved for", area, "\n")
    } else {
      cat("⚠️  ", prefix, "GSEA failed or no significant results for", area, "\n")
    }
  }

  save_and_plot_ora <- function(ora_obj, prefix, area, area_tag, top_n = 10) {
    if (!is.null(ora_obj) && nrow(as.data.frame(ora_obj)) > 0) {
      ora_df <- as.data.frame(ora_obj)
      write.csv(ora_df,
        file = file.path(OUTPUT_DIR, paste0(prefix, "_ora_results_", area_tag, ".csv")),
        row.names = FALSE
      )
      # Plot top pathways
      top_ora <- head(ora_df, top_n)
      p <- ggplot(top_ora, aes(x = reorder(Description, Count), y = Count, fill = p.adjust < 0.05)) +
        geom_col() +
        coord_flip() +
        labs(
          title = paste("Top", prefix, "ORA Pathways -", area),
          x = "Pathway", y = "Gene Count"
        ) +
        scale_fill_manual(values = c("grey70", "tomato")) +
        theme_minimal()
      ggsave(
        filename = file.path(OUTPUT_DIR, paste0(prefix, "_ora_top_", area_tag, ".png")),
        plot = p, width = 8, height = 5, dpi = PARAMS$figure_dpi
      )
      cat(prefix, "ORA results and plot saved for", area, "\n")
    } else {
      cat("⚠️  ", prefix, "ORA failed or no significant results for", area, "\n")
    }
  }

  # ---- Prepare gene lists ----
  if (requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    library(clusterProfiler)
    library(org.Hs.eg.db)
    # Prepare ranked gene list for GSEA
    gene_list <- de_results_df$effect_size
    names(gene_list) <- de_results_df$gene
    gene_list <- sort(gene_list, decreasing = TRUE)
    # Map gene symbols to Entrez IDs
    gene_df <- bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    gene_list_entrez <- gene_list[gene_df$SYMBOL]
    names(gene_list_entrez) <- gene_df$ENTREZID
    gene_list_entrez <- gene_list_entrez[!is.na(names(gene_list_entrez))]
    sig_gene_symbols <- sig_results$gene
    sig_gene_entrez <- gene_df$ENTREZID[gene_df$SYMBOL %in% sig_gene_symbols]

    # ---- Hallmark gene sets (MSigDB) ----
    if (requireNamespace("msigdbr", quietly = TRUE)) {
      library(msigdbr)
      hallmark <- msigdbr(species = "Homo sapiens", category = "H")
      hallmark_list <- split(hallmark$entrez_gene, hallmark$gs_name)
      # GSEA
      gsea_hallmark <- tryCatch(
        {
          GSEA(geneList = gene_list_entrez, TERM2GENE = hallmark_list, pvalueCutoff = 0.05, verbose = FALSE)
        },
        error = function(e) NULL
      )
      save_and_plot_gsea(gsea_hallmark, "hallmark", area, area_tag)
      # ORA
      ora_hallmark <- tryCatch(
        {
          enricher(gene = sig_gene_entrez, TERM2GENE = hallmark_list, pvalueCutoff = 0.05)
        },
        error = function(e) NULL
      )
      save_and_plot_ora(ora_hallmark, "hallmark", area, area_tag)
    } else {
      cat("⚠️  msigdbr not available, skipping Hallmark for", area, "\n")
    }

    # ---- GO: BP, CC, MF ----
    for (ont in c("BP", "CC", "MF")) {
      # GSEA
      gsea_go <- tryCatch(
        {
          gseGO(
            geneList = gene_list_entrez,
            OrgDb = org.Hs.eg.db,
            ont = ont,
            keyType = "ENTREZID",
            minGSSize = 10,
            maxGSSize = 500,
            pvalueCutoff = 0.05,
            verbose = FALSE
          )
        },
        error = function(e) NULL
      )
      save_and_plot_gsea(gsea_go, paste0("GO_", ont), area, area_tag)
      # ORA
      ora_go <- tryCatch(
        {
          enrichGO(
            gene = sig_gene_entrez,
            OrgDb = org.Hs.eg.db,
            ont = ont,
            keyType = "ENTREZID",
            pvalueCutoff = 0.05
          )
        },
        error = function(e) NULL
      )
      save_and_plot_ora(ora_go, paste0("GO_", ont), area, area_tag)
    }

    # ---- KEGG ----
    if (requireNamespace("KEGGREST", quietly = TRUE)) {
      kegg_org <- "hsa"
      # GSEA
      gsea_kegg <- tryCatch(
        {
          gseKEGG(
            geneList = gene_list_entrez,
            organism = kegg_org,
            pvalueCutoff = 0.05,
            verbose = FALSE
          )
        },
        error = function(e) NULL
      )
      save_and_plot_gsea(gsea_kegg, "KEGG", area, area_tag)
      # ORA
      ora_kegg <- tryCatch(
        {
          enrichKEGG(
            gene = sig_gene_entrez,
            organism = kegg_org,
            pvalueCutoff = 0.05
          )
        },
        error = function(e) NULL
      )
      save_and_plot_ora(ora_kegg, "KEGG", area, area_tag)
    } else {
      cat("⚠️  KEGGREST not available, skipping KEGG for", area, "\n")
    }

    # ---- Reactome ----
    if (requireNamespace("ReactomePA", quietly = TRUE)) {
      library(ReactomePA)
      # GSEA
      gsea_reactome <- tryCatch(
        {
          gsePathway(
            geneList = gene_list_entrez,
            organism = "human",
            pvalueCutoff = 0.05,
            verbose = FALSE
          )
        },
        error = function(e) NULL
      )
      save_and_plot_gsea(gsea_reactome, "Reactome", area, area_tag)
      # ORA
      ora_reactome <- tryCatch(
        {
          enrichPathway(
            gene = sig_gene_entrez,
            organism = "human",
            pvalueCutoff = 0.05,
            readable = TRUE
          )
        },
        error = function(e) NULL
      )
      save_and_plot_ora(ora_reactome, "Reactome", area, area_tag)
    } else {
      cat("⚠️  ReactomePA not available, skipping Reactome for", area, "\n")
    }

    # ---- Transcription Factor (TF) targets (MSigDB C3:TFT) ----
    if (requireNamespace("msigdbr", quietly = TRUE)) {
      tfsets <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT")
      tf_list <- split(tfsets$entrez_gene, tfsets$gs_name)
      # GSEA
      gsea_tf <- tryCatch(
        {
          GSEA(geneList = gene_list_entrez, TERM2GENE = tf_list, pvalueCutoff = 0.05, verbose = FALSE)
        },
        error = function(e) NULL
      )
      save_and_plot_gsea(gsea_tf, "TF", area, area_tag)
      # ORA
      ora_tf <- tryCatch(
        {
          enricher(gene = sig_gene_entrez, TERM2GENE = tf_list, pvalueCutoff = 0.05)
        },
        error = function(e) NULL
      )
      save_and_plot_ora(ora_tf, "TF", area, area_tag)
    } else {
      cat("⚠️  msigdbr not available, skipping TFs for", area, "\n")
    }

    # ---- Flexible Filtering: Inflammatory Pathways ----
    # Example: Save filtered results for pathways containing "inflamm", "cytokine", "immune"
    filter_keywords <- c("inflamm", "cytokine", "immune", "complement")
    filter_and_save <- function(df, prefix, area_tag) {
      if (!is.null(df) && nrow(df) > 0) {
        hits <- df[Reduce(`|`, lapply(filter_keywords, function(k) grepl(k, df$Description, ignore.case = TRUE))), ]
        if (nrow(hits) > 0) {
          write.csv(hits, file = file.path(OUTPUT_DIR, paste0(prefix, "_filtered_inflammatory_", area_tag, ".csv")), row.names = FALSE)
        }
      }
    }
    # Apply to all GSEA/ORA results above (example for Hallmark GSEA)
    if (exists("gsea_hallmark")) filter_and_save(as.data.frame(gsea_hallmark), "hallmark_gsea", area_tag)
    if (exists("ora_hallmark")) filter_and_save(as.data.frame(ora_hallmark), "hallmark_ora", area_tag)
    if (exists("gsea_go")) filter_and_save(as.data.frame(gsea_go), "GO_gsea", area_tag)
    if (exists("ora_go")) filter_and_save(as.data.frame(ora_go), "GO_ora", area_tag)
    if (exists("gsea_kegg")) filter_and_save(as.data.frame(gsea_kegg), "KEGG_gsea", area_tag)
    if (exists("ora_kegg")) filter_and_save(as.data.frame(ora_kegg), "KEGG_ora", area_tag)
    if (exists("gsea_reactome")) filter_and_save(as.data.frame(gsea_reactome), "Reactome_gsea", area_tag)
    if (exists("ora_reactome")) filter_and_save(as.data.frame(ora_reactome), "Reactome_ora", area_tag)
    if (exists("gsea_tf")) filter_and_save(as.data.frame(gsea_tf), "TF_gsea", area_tag)
    if (exists("ora_tf")) filter_and_save(as.data.frame(ora_tf), "TF_ora", area_tag)
  } else {
    cat("⚠️  clusterProfiler/org.Hs.eg.db not available, skipping all GSEA/ORA for", area, "\n")
  }
}

cat("\n✅ Panel A analysis complete. Results saved to:", OUTPUT_DIR, "\n")
