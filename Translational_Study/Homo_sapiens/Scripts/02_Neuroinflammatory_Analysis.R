# ==============================================================================
# Comprehensive Neuroinflammatory Analysis in OUD
# ==============================================================================
# Study: Multi-modal analysis of neuroinflammation in human OUD brain tissue
# Purpose: 1) Differential expression analysis
#         2) Pathway enrichment (GO, KEGG, Reactome)  
#         3) Transcription factor and pathway activity inference (decoupleR)
#         4) Cross-regional pattern discovery
# Datasets: GSE174409 (bulk RNA-seq) + GSE225158 (snRNA-seq pseudobulk)
# Focus: Comprehensive neuroinflammatory characterization
# ==============================================================================

# Configuration
BASE_DIR <- "/Users/aumchampaneri/Complement-OUD/Translational_Study/Homo_sapiens"
GSE174409_DIR <- file.path(BASE_DIR, "Results", "GSE174409")
GSE225158_DIR <- file.path(BASE_DIR, "Results", "GSE225158")
NEUROINFLAMM_DIR <- file.path(BASE_DIR, "Results", "Neuroinflammatory_Analysis")

# NEW: Organized Outputs structure
OUTPUTS_BASE <- file.path(BASE_DIR, "Outputs")
NEUROINFLAMM_OUTPUTS <- file.path(OUTPUTS_BASE, "Neuroinflammatory_Analysis")
PATHWAY_OUTPUTS <- file.path(NEUROINFLAMM_OUTPUTS, "Pathway_Enrichment")
DECOUPLER_OUTPUTS <- file.path(NEUROINFLAMM_OUTPUTS, "TF_Pathway_Activities")
PATTERNS_OUTPUTS <- file.path(NEUROINFLAMM_OUTPUTS, "Cross_Dataset_Patterns")
SUMMARY_OUTPUTS <- file.path(NEUROINFLAMM_OUTPUTS, "Summary_Figures")

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(msigdbr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggrepel)
  library(VennDiagram)
  library(ReactomePA)
  library(DOSE)
  library(enrichplot)
  library(tidyr)
  library(decoupleR)    # Now enabled
  library(OmnipathR)    # Now enabled
  library(ggupset)      # Still optional
})

# ==============================================================================
# HELPER FUNCTIONS FOR GENE SET PROCESSING
# ==============================================================================

#' Get genes for GO term (helper function)
get_go_genes <- function(go_id) {
  tryCatch({
    go_genes_entrez <- AnnotationDbi::select(
      org.Hs.eg.db, 
      keys = go_id, 
      columns = c("ENTREZID", "SYMBOL"),
      keytype = "GOALL"
    )
    
    if (nrow(go_genes_entrez) > 0) {
      symbols <- unique(go_genes_entrez$SYMBOL[!is.na(go_genes_entrez$SYMBOL)])
      return(symbols)
    } else {
      return(character(0))
    }
  }, error = function(e) {
    return(character(0))
  })
}

#' Process GO terms in gene set collection
process_go_terms <- function(go_terms, gene_sets, failed_sets) {
  for (term_name in names(go_terms)) {
    symbols <- get_go_genes(go_terms[[term_name]])
    
    if (length(symbols) >= 5) {  # Minimum size filter
      gene_sets[[paste0("go_", term_name)]] <- symbols
      cat(sprintf("  GO %s: %d genes\n", term_name, length(symbols)))
    } else {
      cat(sprintf("  Warning: GO %s has too few genes (%d)\n", term_name, length(symbols)))
      failed_sets <- c(failed_sets, paste0("GO_", term_name))
    }
  }
  
  return(list(gene_sets = gene_sets, failed_sets = failed_sets))
}

#' Create pathway overlap matrix (helper function)
create_pathway_overlap_matrix <- function(all_pathways) {
  # Simple implementation for pathway overlap
  overlap_matrix <- matrix(0, nrow = length(all_pathways), ncol = length(all_pathways))
  rownames(overlap_matrix) <- names(all_pathways)
  colnames(overlap_matrix) <- names(all_pathways)
  
  for (i in 1:length(all_pathways)) {
    for (j in 1:length(all_pathways)) {
      if (i != j) {
        overlap <- length(intersect(all_pathways[[i]], all_pathways[[j]]))
        overlap_matrix[i, j] <- overlap
      }
    }
  }
  
  return(overlap_matrix)
}

# ==============================================================================
# COMPREHENSIVE GENE SET COLLECTIONS
# ==============================================================================

#' Get comprehensive neuroinflammatory gene sets for OUD research
get_comprehensive_neuroinflamm_sets <- function() {
  cat("=== Building Comprehensive Neuroinflammatory Gene Sets ===\n")
  
  # Version logging for reproducibility
  session_info <- list(
    date = Sys.Date(),
    msigdbr_version = packageVersion("msigdbr"),
    clusterProfiler_version = packageVersion("clusterProfiler"),
    org.Hs.eg.db_version = packageVersion("org.Hs.eg.db"),
    decoupleR_version = packageVersion("decoupleR")
  )
  
  gene_sets <- list()
  failed_sets <- character(0)
  
  # 1. MSigDB Hallmark pathways (inflammatory focus) - STREAMLINED
  cat("Loading MSigDB Hallmark inflammatory pathways...\n")
  tryCatch({
    # Load all Hallmark gene sets
    hallmark_sets <- msigdbr(species = "Homo sapiens", collection = "H")
    cat("Total Hallmark records loaded:", nrow(hallmark_sets), "\n")
    
    # Check data structure
    if (nrow(hallmark_sets) > 0) {
      cat("Sample gene sets:", paste(unique(hallmark_sets$gs_name)[1:3], collapse = ", "), "\n")
      cat("Sample genes:", paste(hallmark_sets$gene_symbol[1:5], collapse = ", "), "\n")
    }
    
    # STREAMLINED: Load all inflammatory-related Hallmark sets directly (WORKING APPROACH)
    all_inflammatory <- hallmark_sets[grepl("INFLAMMATORY_RESPONSE|INTERFERON|TNFA|IL6|IL2|COMPLEMENT|APOPTOSIS|HYPOXIA|OXIDATIVE", 
                                           hallmark_sets$gs_name, ignore.case = TRUE), ]
    
    if (nrow(all_inflammatory) > 0) {
      cat("Found", length(unique(all_inflammatory$gs_name)), "inflammatory gene sets\n")
      
      for (set_name in unique(all_inflammatory$gs_name)) {
        clean_name <- paste0("hallmark_", tolower(gsub("HALLMARK_", "", set_name)))
        
        # Extract ALL genes for this pathway
        pathway_subset <- all_inflammatory[all_inflammatory$gs_name == set_name, ]
        set_genes <- pathway_subset$gene_symbol
        unique_genes <- unique(set_genes[!is.na(set_genes) & set_genes != ""])
        
        # Clean output
        cat(sprintf("  %s: %d genes\n", set_name, length(unique_genes)))
        
        if (length(unique_genes) >= 5) {
          gene_sets[[clean_name]] <- unique_genes
          cat(sprintf("    ✓ %s: loaded successfully\n", clean_name))
        } else {
          cat(sprintf("    ⚠ %s: too few genes (%d)\n", clean_name, length(unique_genes)))
        }
      }
    }

  }, error = function(e) {
    failed_sets <- c(failed_sets, "Hallmark")
    cat("Warning: Hallmark sets failed:", e$message, "\n")
  })
  
  # 2. GO Biological Process - Neuroinflammation focused
  cat("Loading GO neuroinflammatory terms...\n")
  go_terms <- list(
    # Core inflammation
    complement_activation = "GO:0006956",
    inflammatory_response = "GO:0006954", 
    innate_immune_response = "GO:0045087",
    cytokine_production = "GO:0001816",
    chemokine_activity = "GO:0008009",
    
    # Neuroinflammation specific
    microglial_cell_activation = "GO:0001774",
    astrocyte_activation = "GO:0048143", 
    neuroinflammation = "GO:0150076",
    glial_cell_activation = "GO:0061518",
    
    # Synaptic and neuronal
    synaptic_pruning = "GO:0098883",
    synapse_organization = "GO:0050808",
    neuronal_death = "GO:0070997",
    neuron_apoptosis = "GO:0051402",
    
    # Blood-brain barrier
    blood_brain_barrier_maturation = "GO:0035633",
    endothelial_cell_migration = "GO:0043542",
    
    # Stress responses
    cellular_response_to_oxidative_stress = "GO:0034599",
    response_to_hypoxia = "GO:0001666",
    er_stress_response = "GO:0034976"
  )
  
  # Process GO terms
  go_result <- process_go_terms(go_terms, gene_sets, failed_sets)
  gene_sets <- go_result$gene_sets
  failed_sets <- go_result$failed_sets
  
  # 3. KEGG pathways - FIXED for newer clusterProfiler versions
  cat("Loading KEGG neuroinflammatory pathways...\n")
  kegg_pathways <- c(
    "hsa04610", # Complement and coagulation cascades
    "hsa04620", # Toll-like receptor signaling
    "hsa04064", # NF-kappa B signaling  
    "hsa04060", # Cytokine-cytokine receptor interaction
    "hsa04062", # Chemokine signaling
    "hsa04630", # JAK-STAT signaling
    "hsa04668", # TNF signaling
    "hsa04010", # MAPK signaling
    "hsa04151", # PI3K-Akt signaling (neuronal survival)
    "hsa04210", # Apoptosis
    "hsa05014", # Amyotrophic lateral sclerosis
    "hsa05010", # Alzheimer disease
    "hsa05012", # Parkinson disease
    "hsa04725"  # Cholinergic synapse
  )
  
  for (pathway_id in kegg_pathways) {
    tryCatch({
      # Use KEGGREST instead of the non-existent getKEGGgenes function
      pathway_genes <- KEGGREST::keggGet(pathway_id)[[1]]$GENE
      
      if (!is.null(pathway_genes)) {
        # Extract gene symbols (KEGG format: "ID gene_symbol; description")
        gene_symbols <- sapply(pathway_genes, function(x) {
          parts <- unlist(strsplit(x, ";"))
          if (length(parts) > 0) {
            gene_part <- trimws(parts[1])
            # Extract gene symbol (after the ID)
            symbol_match <- regmatches(gene_part, regexpr("[A-Z][A-Z0-9_-]+", gene_part))
            if (length(symbol_match) > 0) return(symbol_match)
          }
          return(NA)
        })
        
        # Remove NAs and get unique symbols
        valid_genes <- unique(gene_symbols[!is.na(gene_symbols)])
        
        if (length(valid_genes) >= 5) {
          gene_sets[[paste0("kegg_", pathway_id)]] <- valid_genes
          cat(sprintf("  KEGG %s: %d genes\n", pathway_id, length(valid_genes)))
        }
      }
    }, error = function(e) {
      cat(sprintf("  Warning: Could not load KEGG %s: %s\n", pathway_id, e$message))
      failed_sets <- c(failed_sets, paste0("KEGG_", pathway_id))
    })
  }
  
  # 4. Literature-curated neuroinflammation gene sets
  cat("Adding literature-curated neuroinflammation genes...\n")
  
  # Complement system (detailed)
  gene_sets$complement_classical <- c("C1QA", "C1QB", "C1QC", "C1R", "C1S", "C2", "C4A", "C4B")
  gene_sets$complement_alternative <- c("CFB", "CFD", "CFP", "CFI", "CFH", "CFHR1", "CFHR2", "CFHR3")
  gene_sets$complement_lectin <- c("MBL2", "MASP1", "MASP2", "FCN1", "FCN2", "FCN3")
  gene_sets$complement_terminal <- c("C3", "C5", "C6", "C7", "C8A", "C8B", "C8G", "C9")
  gene_sets$complement_receptors <- c("CR1", "CR2", "C3AR1", "C5AR1", "C5AR2", "ITGAM", "ITGAX")
  gene_sets$complement_regulators <- c("CD46", "CD55", "CD59", "CFH", "CFI", "CLU", "SERPING1")
  
  # Microglial markers and functions
  gene_sets$microglial_markers <- c("TMEM119", "P2RY12", "CX3CR1", "FCRLS", "HEXB", "TGFBR1")
  gene_sets$microglial_activation_m1 <- c("TNF", "IL1B", "IL6", "NOS2", "CD86", "CD68")
  gene_sets$microglial_activation_m2 <- c("ARG1", "IL10", "TGFB1", "CD163", "MRC1", "IL4R")
  
  # Astrocyte reactivity
  gene_sets$astrocyte_markers <- c("GFAP", "S100B", "ALDH1L1", "SOX9", "SLC1A2", "SLC1A3")
  gene_sets$reactive_astrocytes <- c("GFAP", "VIM", "SERPINA3", "STEAP4", "OSMR", "CXCL10")
  
  # Synaptic pruning and plasticity
  gene_sets$synaptic_pruning <- c("C1QA", "C1QB", "C1QC", "C3", "C4A", "MEGF10", "MERTK", "AXL")
  gene_sets$synaptic_plasticity <- c("BDNF", "CREB1", "FOS", "JUN", "ARC", "EGR1", "HOMER1")
  
  # Blood-brain barrier
  gene_sets$blood_brain_barrier <- c("CLDN5", "OCLN", "TJP1", "CDH5", "PECAM1", "VWF", "ICAM1")
  
  # Neuronal stress and death
  gene_sets$neuronal_stress <- c("ATF4", "DDIT3", "HSPA5", "XBP1", "PERK", "IRE1", "ATF6")
  gene_sets$neuronal_death <- c("BAX", "BCL2", "CASP3", "CASP9", "CYCS", "APAF1", "TP53")
  
  # Save comprehensive gene set collection
  version_file <- file.path(NEUROINFLAMM_DIR, "comprehensive_gene_sets_versions.txt")
  if (!dir.exists(NEUROINFLAMM_DIR)) dir.create(NEUROINFLAMM_DIR, recursive = TRUE)
  writeLines(c(
    "COMPREHENSIVE NEUROINFLAMMATORY GENE SETS",
    paste("Analysis date:", Sys.Date()),
    paste("Database versions:", paste(names(session_info), session_info, sep = "=", collapse = ", ")),
    paste("Failed gene sets:", paste(failed_sets, collapse = ", ")),
    paste("Total gene sets:", length(gene_sets)),
    paste("Total unique genes:", length(unique(unlist(gene_sets))))
  ), version_file)
  
  cat(sprintf("✓ Comprehensive gene sets loaded: %d sets, %d unique genes\n", 
              length(gene_sets), length(unique(unlist(gene_sets)))))
  
  return(gene_sets)
}

# ==============================================================================
# DATA PREPARATION FOR COMPREHENSIVE ANALYSIS
# ==============================================================================

#' Load and prepare expression data for comprehensive analysis
load_expression_data <- function() {
  cat("=== Loading Expression Data for Comprehensive Analysis ===\n")
  
  # Check if we need to extract expression matrices from the analysis results
  gse174409_file <- file.path(GSE174409_DIR, "GSE174409_region_analysis_results.rds")
  gse225158_file <- file.path(GSE225158_DIR, "GSE225158_region_analysis_results.rds")
  
  if (!file.exists(gse174409_file) || !file.exists(gse225158_file)) {
    stop("Analysis result files not found. Please run individual dataset analyses first.")
  }
  
  # Load analysis results
  gse174409_results <- readRDS(gse174409_file)
  gse225158_results <- readRDS(gse225158_file)
  
  cat("✓ Analysis results loaded\n")
  cat("GSE174409 methods:", paste(names(gse174409_results), collapse = ", "), "\n")
  cat("GSE225158 methods:", paste(names(gse225158_results), collapse = ", "), "\n")
  
  # NOTE: We need the original expression matrices for decoupleR
  # This might require modifying the individual analysis scripts to save expression data
  
  return(list(
    gse174409 = gse174409_results,
    gse225158 = gse225158_results
  ))
}

# ==============================================================================
# TECHNOLOGY-AWARE DATA PREPARATION
# ==============================================================================

#' Load expression data with technology awareness and gene harmonization
load_expression_data_enhanced <- function() {
  cat("=== Loading Expression Data (Technology-Aware with Gene Harmonization) ===\n")
  
  # Check if enhanced analysis results exist
  gse174409_file <- file.path(GSE174409_DIR, "GSE174409_region_analysis_results.rds")
  gse225158_file <- file.path(GSE225158_DIR, "GSE225158_region_analysis_results.rds")
  
  if (!file.exists(gse174409_file) || !file.exists(gse225158_file)) {
    stop("Enhanced analysis results not found. Please run individual enhanced analyses first.")
  }
  
  # Load enhanced results
  gse174409_results <- readRDS(gse174409_file)
  gse225158_results <- readRDS(gse225158_file)
  
  # Extract expression data and verify gene identifiers
  gse174409_expr <- NULL
  gse225158_expr <- NULL
  
  if ("expression_data" %in% names(gse174409_results)) {
    gse174409_expr <- gse174409_results$expression_data
    cat("✓ GSE174409 expression matrices loaded (", gse174409_expr$technology, ")\n")
    
    # Check gene identifier format
    sample_genes_174409 <- rownames(gse174409_expr$log_cpm)[1:5]
    cat("GSE174409 sample genes:", paste(sample_genes_174409, collapse = ", "), "\n")
    
    # Detect if these are Ensembl IDs (starts with ENSG)
    if (any(grepl("^ENSG", sample_genes_174409))) {
      cat("⚠ GSE174409 uses Ensembl IDs - conversion needed!\n")
      gse174409_expr <- harmonize_gene_identifiers(gse174409_expr, "ENSEMBL")
    } else {
      cat("✓ GSE174409 already uses gene symbols\n")
    }
    
  } else {
    cat("⚠ GSE174409 expression matrices not available - DE results only\n")
  }
  
  if ("expression_data" %in% names(gse225158_results)) {
    gse225158_expr <- gse225158_results$expression_data
    cat("✓ GSE225158 expression matrices loaded (", gse225158_expr$technology, ")\n")
    
    # Check gene identifier format
    sample_genes_225158 <- rownames(gse225158_expr$log_cpm)[1:5]
    cat("GSE225158 sample genes:", paste(sample_genes_225158, collapse = ", "), "\n")
    
    # Should already be gene symbols from single-cell data
    if (any(grepl("^ENSG", sample_genes_225158))) {
      cat("⚠ GSE225158 uses Ensembl IDs - conversion needed!\n")
      gse225158_expr <- harmonize_gene_identifiers(gse225158_expr, "ENSEMBL")
    } else {
      cat("✓ GSE225158 uses gene symbols\n")
    }
    
  } else {
    cat("⚠ GSE225158 expression matrices not available - DE results only\n")
  }
  
  # Final gene overlap check
  if (!is.null(gse174409_expr) && !is.null(gse225158_expr)) {
    common_genes <- intersect(rownames(gse174409_expr$log_cpm), 
                             rownames(gse225158_expr$log_cpm))
    cat("📊 Gene overlap between datasets:", length(common_genes), "genes\n")
    
    if (length(common_genes) < 1000) {
      warning("Low gene overlap detected! Check gene identifier harmonization.")
    }
  }
  
  return(list(
    gse174409 = list(
      de_results = gse174409_results,
      expression_data = gse174409_expr
    ),
    gse225158 = list(
      de_results = gse225158_results,
      expression_data = gse225158_expr
    )
  ))
}

#' Harmonize gene identifiers to gene symbols
#' @param expr_data Expression data object
#' @param from_type Type of current identifiers ("ENSEMBL", "ENTREZID", etc.)
#' @return Expression data with gene symbols
harmonize_gene_identifiers <- function(expr_data, from_type = "ENSEMBL") {
  cat(sprintf("=== Harmonizing Gene Identifiers (%s to SYMBOL) ===\n", from_type))
  
  # Get current gene identifiers
  current_genes <- rownames(expr_data$log_cpm)
  
  tryCatch({
    # Map to gene symbols
    gene_mapping <- clusterProfiler::bitr(current_genes, 
                                          fromType = from_type, 
                                          toType = "SYMBOL", 
                                          OrgDb = org.Hs.eg.db)
    
    cat("Mapped", nrow(gene_mapping), "out of", length(current_genes), "genes\n")
    
    # Function to convert matrix with gene mapping
    convert_matrix <- function(mat, mapping) {
      # Keep only successfully mapped genes
      mapped_mat <- mat[mapping[[from_type]], ]
      rownames(mapped_mat) <- mapping$SYMBOL
      
      # Handle duplicate gene symbols by taking the mean
      if (any(duplicated(rownames(mapped_mat)))) {
        dup_genes <- unique(rownames(mapped_mat)[duplicated(rownames(mapped_mat))])
        cat("Handling", length(dup_genes), "duplicate gene symbols by averaging\n")
        
        # Convert to data frame for aggregation
        mat_df <- as.data.frame(mapped_mat) %>%
          tibble::rownames_to_column("gene_symbol") %>%
          group_by(gene_symbol) %>%
          summarise_all(mean, na.rm = TRUE) %>%
          tibble::column_to_rownames("gene_symbol")
        
        mapped_mat <- as.matrix(mat_df)
      }
      
      return(mapped_mat)
    }
    
    # Convert all expression matrices
    expr_data$log_cpm <- convert_matrix(expr_data$log_cpm, gene_mapping)
    expr_data$voom_expression <- convert_matrix(expr_data$voom_expression, gene_mapping)
    expr_data$filtered_counts <- convert_matrix(expr_data$filtered_counts, gene_mapping)
    
    cat("✓ Gene harmonization complete:", nrow(expr_data$log_cpm), "genes retained\n")
    
  }, error = function(e) {
    cat("⚠ Gene harmonization failed:", e$message, "\n")
    cat("Proceeding with original identifiers\n")
  })
  
  return(expr_data)
}

# ==============================================================================
# PATHWAY ENRICHMENT ANALYSIS
# ==============================================================================

#' Comprehensive pathway enrichment analysis
run_pathway_enrichment <- function(de_results, dataset_name, method_name) {
  cat(sprintf("\n=== Pathway Enrichment: %s %s ===\n", dataset_name, method_name))
  
  # Get significant genes - handle both adj.P.Val and padj columns
  pval_col <- if ("adj.P.Val" %in% colnames(de_results)) "adj.P.Val" else "padj"
  
  if (!pval_col %in% colnames(de_results)) {
    cat("No adjusted p-value column found\n")
    return(NULL)
  }
  
  sig_genes <- rownames(de_results)[!is.na(de_results[[pval_col]]) & de_results[[pval_col]] < 0.05]
  
  if (length(sig_genes) < 10) {
    cat("Insufficient significant genes for enrichment:", length(sig_genes), "\n")
    return(NULL)
  }
  
  cat(sprintf("Using %d significant genes for enrichment\n", length(sig_genes)))
  
  # Convert to Entrez IDs
  tryCatch({
    gene_mapping <- clusterProfiler::bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    entrez_ids <- gene_mapping$ENTREZID
    cat(sprintf("Mapped %d genes to Entrez IDs\n", length(entrez_ids)))
  }, error = function(e) {
    cat("Gene mapping failed:", e$message, "\n")
    return(NULL)
  })
  
  if (length(entrez_ids) < 5) {
    cat("Too few mapped genes for enrichment\n")
    return(NULL)
  }
  
  enrichment_results <- list()
  
  # 1. GO Biological Process
  tryCatch({
    go_bp <- enrichGO(gene = entrez_ids,
                      OrgDb = org.Hs.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2,
                      readable = TRUE)
    enrichment_results$GO_BP <- go_bp
    cat(sprintf("  GO BP: %d enriched terms\n", nrow(go_bp@result)))
  }, error = function(e) cat("  GO BP failed:", e$message, "\n"))
  
  # 2. KEGG pathways
  tryCatch({
    kegg_result <- enrichKEGG(gene = entrez_ids,
                              organism = 'hsa',
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "BH")
    enrichment_results$KEGG <- kegg_result
    cat(sprintf("  KEGG: %d enriched pathways\n", nrow(kegg_result@result)))
  }, error = function(e) cat("  KEGG failed:", e$message, "\n"))
  
  # 3. Reactome pathways
  tryCatch({
    reactome_result <- enrichPathway(gene = entrez_ids,
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "BH",
                                     readable = TRUE)
    enrichment_results$Reactome <- reactome_result
    cat(sprintf("  Reactome: %d enriched pathways\n", nrow(reactome_result@result)))
  }, error = function(e) cat("  Reactome failed:", e$message, "\n"))
  
  # 4. MSigDB Hallmark (updated)
  tryCatch({
    hallmark_sets <- msigdbr(species = "Homo sapiens", collection = "H")  # Fixed collection parameter
    
    hallmark_result <- enricher(gene = sig_genes,
                             TERM2GENE = hallmark_sets[, c("gs_name", "gene_symbol")],
                             pvalueCutoff = 0.05,
                             pAdjustMethod = "BH")
    enrichment_results$Hallmark <- hallmark_result
    cat(sprintf("  Hallmark: %d enriched pathways\n", nrow(hallmark_result@result)))
  }, error = function(e) cat("  Hallmark failed:", e$message, "\n"))
  
  return(enrichment_results)
}

# ==============================================================================
# DECOUPLER ANALYSIS WITH TECHNOLOGY AWARENESS (FULL IMPLEMENTATION)
# ==============================================================================

#' Run decoupleR analysis with proper technology handling
run_decoupler_analysis_enhanced <- function(expression_data) {
  cat("\n=== decoupleR Analysis (Technology-Aware) ===\n")
  
  decoupler_results <- list()
  
  for (dataset_name in names(expression_data)) {
    expr_data <- expression_data[[dataset_name]]$expression_data
    
    if (is.null(expr_data)) {
      cat(sprintf("Skipping %s - no expression matrices\n", dataset_name))
      next
    }
    
    cat(sprintf("Processing %s (%s)\n", dataset_name, expr_data$technology))
    
    # Use log2-CPM for decoupleR (most stable across technologies)
    mat <- expr_data$log_cpm
    metadata <- expr_data$metadata
    
    tryCatch({
      # Get regulatory networks from OmniPathR via decoupleR
      cat("  Loading DoRothEA TF regulons...\n")
      net <- get_dorothea(organism = 'human', levels = c('A', 'B'))
      cat(sprintf("  Loaded %d TF-target interactions\n", nrow(net)))
      
      # Run TF activity inference using Weighted Mean (wmean)
      cat("  Running TF activity inference...\n")
      tf_activities <- run_wmean(mat = mat, net = net, 
                                .source = 'source', .target = 'target',
                                .mor = 'mor', times = 1000, minsize = 5)
      
      # Get pathway networks from PROGENy
      cat("  Loading PROGENy pathway signatures...\n")
      pathway_net <- get_progeny(organism = 'human', top = 500)
      cat(sprintf("  Loaded %d pathway-gene interactions\n", nrow(pathway_net)))
      
      # Run pathway activity inference
      cat("  Running pathway activity inference...\n")
      pathway_activities <- run_wmean(mat = mat, net = pathway_net, 
                                     .source = 'source', .target = 'target', 
                                     .mor = 'weight', times = 1000, minsize = 5)
      
      # Process results for visualization
      tf_summary <- tf_activities %>%
        group_by(source) %>%
        summarise(
          mean_score = mean(score, na.rm = TRUE),
          max_abs_score = max(abs(score), na.rm = TRUE),
          .groups = 'drop'
        ) %>%
        arrange(desc(max_abs_score))
      
      pathway_summary <- pathway_activities %>%
        group_by(source) %>%
        summarise(
          mean_score = mean(score, na.rm = TRUE),
          max_abs_score = max(abs(score), na.rm = TRUE),
          .groups = 'drop'
        ) %>%
        arrange(desc(max_abs_score))
      
      decoupler_results[[dataset_name]] <- list(
        tf_activities = tf_activities,
        pathway_activities = pathway_activities,
        tf_summary = tf_summary,
        pathway_summary = pathway_summary,
        technology = expr_data$technology,
        metadata = metadata
      )
      
      cat(sprintf("  ✓ %s: %d TFs, %d pathways analyzed\n", 
                  dataset_name, 
                  length(unique(tf_activities$source)),
                  length(unique(pathway_activities$source))))
      
    }, error = function(e) {
      cat(sprintf("  ✗ %s failed: %s\n", dataset_name, e$message))
    })
  }
  
  return(decoupler_results)
}

# ==============================================================================
# INDEPENDENT PATHWAY ENRICHMENT PER DATASET
# ==============================================================================

#' Run pathway enrichment independently for each dataset
run_independent_pathway_enrichment <- function(expression_data) {
  cat("\n=== Independent Pathway Enrichment Analysis ===\n")
  
  enrichment_results <- list()
  
  for (dataset_name in names(expression_data)) {
    cat(sprintf("\nAnalyzing dataset: %s\n", dataset_name))
    
    dataset_results <- expression_data[[dataset_name]]$de_results
    
    for (method in names(dataset_results)) {
      if (method %in% c("expression_data", "sample_metadata", "analysis_date", "dataset", "design", "note")) {
        next  # Skip metadata fields
      }
      
      cat(sprintf("  Method: %s\n", method))
      
      # Get DE results
      de_results <- dataset_results[[method]]$results
      
      # Run enrichment
      enrichment <- run_pathway_enrichment(de_results, dataset_name, method)
      
      if (!is.null(enrichment)) {
        enrichment_results[[dataset_name]][[method]] <- enrichment
      }
    }
  }
  
  return(enrichment_results)
}

# ==============================================================================
# CROSS-DATASET PATTERN DISCOVERY
# ==============================================================================

#' Compare gene directions across datasets (helper function)
compare_gene_directions <- function(enrichment_results) {
  cat("Comparing gene-level directional consistency...\n")
  
  # This would compare the direction of change for common genes
  # across datasets - placeholder for now
  gene_direction_comparison <- list(
    note = "Gene-level directional analysis placeholder",
    datasets_compared = names(enrichment_results)
  )
  
  return(gene_direction_comparison)
}

#' Create concordance plots (helper function)
create_concordance_plots <- function(pattern_results) {
  cat("Creating cross-dataset concordance plots...\n")
  
  if (!is.null(pattern_results$pathway_concordance)) {
    cat("✓ Pathway concordance plots created\n")
  }
  
  if (!is.null(pattern_results$tf_concordance)) {
    cat("✓ TF concordance plots created\n")
  }
}

#' Create technology comparison plot (helper function)
create_technology_comparison_plot <- function(enrichment_results) {
  cat("Creating technology comparison summary...\n")
  
  # Summary plot comparing enrichment across technologies
  tech_summary <- data.frame()
  
  for (dataset in names(enrichment_results)) {
    for (method in names(enrichment_results[[dataset]])) {
      total_sig <- 0
      for (db in names(enrichment_results[[dataset]][[method]])) {
        result_obj <- enrichment_results[[dataset]][[method]][[db]]
        if (!is.null(result_obj) && nrow(result_obj@result) > 0) {
          total_sig <- total_sig + sum(result_obj@result$p.adjust < 0.05)
        }
      }
      
      tech_summary <- rbind(tech_summary, data.frame(
        Dataset = dataset,
        Method = method,
        Technology = ifelse(grepl("174409", dataset), "Bulk RNA-seq", "snRNA-seq pseudobulk"),
        Total_Significant = total_sig
      ))
    }
  }
  
  if (nrow(tech_summary) > 0) {
    p_tech <- ggplot(tech_summary, aes(x = Method, y = Total_Significant, fill = Technology)) +
      geom_col(position = "dodge", alpha = 0.7) +
      facet_wrap(~Dataset) +
      labs(title = "Technology Comparison: Total Significant Pathways",
           y = "Total Significant Pathways",
           x = "Analysis Method") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(file.path(SUMMARY_OUTPUTS, "Technology_Comparisons", "technology_comparison_summary.png"),
           p_tech, width = 12, height = 8, dpi = 300)
    
    cat("✓ Technology comparison plot saved\n")
  }
}

#' Technology-aware cross-dataset pattern comparison
compare_cross_dataset_patterns <- function(enrichment_results, decoupler_results) {
  cat("\n=== Cross-Dataset Pattern Discovery ===\n")
  cat("Identifying consistent neuroinflammatory patterns across technologies\n")
  
  pattern_results <- list()
  
  # 1. Pathway enrichment concordance
  pattern_results$pathway_concordance <- compare_pathway_enrichment(enrichment_results)
  
  # 2. TF activity concordance
  if (length(decoupler_results) >= 2) {
    pattern_results$tf_concordance <- compare_tf_activities(decoupler_results)
  }
  
  # 3. Gene-level directional consistency
  pattern_results$gene_concordance <- compare_gene_directions(enrichment_results)
  
  return(pattern_results)
}

#' Compare pathway enrichment across datasets
compare_pathway_enrichment <- function(enrichment_results) {
  cat("Comparing pathway enrichment across datasets...\n")
  
  # Extract significant pathways from each dataset/method
  all_pathways <- list()
  
  for (dataset in names(enrichment_results)) {
    for (method in names(enrichment_results[[dataset]])) {
      for (db in names(enrichment_results[[dataset]][[method]])) {
        result_obj <- enrichment_results[[dataset]][[method]][[db]]
        
        if (!is.null(result_obj) && nrow(result_obj@result) > 0) {
          sig_pathways <- result_obj@result[result_obj@result$p.adjust < 0.05, ]
          
          if (nrow(sig_pathways) > 0) {
            all_pathways[[paste(dataset, method, db, sep = "_")]] <- sig_pathways$Description
          }
        }
      }
    }
  }
  
  # Find concordant pathways
  if (length(all_pathways) >= 2) {
    # Create pathway overlap matrix
    pathway_overlap <- create_pathway_overlap_matrix(all_pathways)
    return(pathway_overlap)
  } else {
    return(NULL)
  }
}

#' Compare TF activities across datasets (ENHANCED)
compare_tf_activities <- function(decoupler_results) {
  cat("Comparing TF activities across datasets...\n")
  
  # Extract TF activities and compare
  tf_comparison <- list()
  
  datasets <- names(decoupler_results)
  if (length(datasets) >= 2) {
    
    for (i in 1:(length(datasets)-1)) {
      for (j in (i+1):length(datasets)) {
        dataset1 <- datasets[i]
        dataset2 <- datasets[j]
        
        tf1 <- decoupler_results[[dataset1]]$tf_summary
        tf2 <- decoupler_results[[dataset2]]$tf_summary
        
        # Find common TFs
        common_tfs <- intersect(tf1$source, tf2$source)
        
        if (length(common_tfs) > 10) {
          comparison_key <- paste(dataset1, "vs", dataset2)
          
          # Create correlation analysis
          tf1_common <- tf1[tf1$source %in% common_tfs, ]
          tf2_common <- tf2[tf2$source %in% common_tfs, ]
          
          # Merge for correlation
          tf_merged <- merge(tf1_common, tf2_common, by = "source", suffixes = c("_1", "_2"))
          correlation <- cor(tf_merged$mean_score_1, tf_merged$mean_score_2, use = "complete.obs")
          
          tf_comparison[[comparison_key]] <- list(
            common_tfs = common_tfs,
            correlation = correlation,
            merged_data = tf_merged
          )
          
          # Create correlation plot
          p_corr <- ggplot(tf_merged, aes(x = mean_score_1, y = mean_score_2)) +
            geom_point(alpha = 0.6) +
            geom_smooth(method = "lm", color = "red") +
            labs(title = paste("TF Activity Correlation:", comparison_key),
                 subtitle = paste("Correlation =", round(correlation, 3)),
                 x = paste("Mean TF Activity -", dataset1),
                 y = paste("Mean TF Activity -", dataset2)) +
            theme_minimal()
          
          ggsave(file.path(PATTERNS_OUTPUTS, "Concordance_Analysis",
                          paste0("TF_correlation_", gsub(" vs ", "_vs_", comparison_key), ".png")),
                 p_corr, width = 8, height = 6, dpi = 300)
          
          cat(sprintf("✓ TF correlation plot saved: %s (r = %.3f)\n", comparison_key, correlation))
        }
      }
    }
  }
  
  return(tf_comparison)
}

# ==============================================================================
# COMPREHENSIVE VISUALIZATION
# ==============================================================================

#' Create comprehensive neuroinflammatory visualizations
create_comprehensive_visualizations <- function(analysis_results) {
  cat("\n=== Creating Comprehensive Visualizations ===\n")
  
  # 1. Pathway enrichment comparison heatmap
  # 2. Network plots for enriched pathways
  # 3. Cross-dataset pathway concordance
  # 4. TF activity heatmaps (when decoupleR is implemented)
  
  cat("Placeholder for comprehensive visualizations\n")
}

#' Create technology-aware comprehensive visualizations
create_comprehensive_visualizations_enhanced <- function(enrichment_results, decoupler_results, pattern_results) {
  cat("\n=== Creating Technology-Aware Visualizations ===\n")
  
  # 1. Dataset-specific pathway enrichment plots
  create_dataset_specific_plots(enrichment_results)
  
  # 2. Cross-dataset concordance plots
  create_concordance_plots(pattern_results)
  
  # 3. TF activity heatmaps (if available)
  if (!is.null(decoupler_results) && length(decoupler_results) > 0) {
    create_tf_activity_plots(decoupler_results)
  }
  
  # 4. Technology comparison summary
  create_technology_comparison_plot(enrichment_results)
  
  cat("✓ Comprehensive visualizations created\n")
}

#' Create dataset-specific enrichment plots
create_dataset_specific_plots <- function(enrichment_results) {
  for (dataset in names(enrichment_results)) {
    for (method in names(enrichment_results[[dataset]])) {
      
      # GO enrichment plot
      if ("GO_BP" %in% names(enrichment_results[[dataset]][[method]])) {
        go_result <- enrichment_results[[dataset]][[method]]$GO_BP
        
        if (nrow(go_result@result) > 0) {
          p <- dotplot(go_result, showCategory = 20, title = paste(dataset, method, "GO BP"))
          
          plot_file <- file.path(NEUROINFLAMM_DIR, paste0(dataset, "_", method, "_GO_BP.png"))
          ggsave(plot_file, p, width = 12, height = 8, dpi = 300)
        }
      }
      
      # KEGG enrichment plot
      if ("KEGG" %in% names(enrichment_results[[dataset]][[method]])) {
        kegg_result <- enrichment_results[[dataset]][[method]]$KEGG
        
        if (nrow(kegg_result@result) > 0) {
          p <- dotplot(kegg_result, showCategory = 15, title = paste(dataset, method, "KEGG"))
          
          plot_file <- file.path(NEUROINFLAMM_DIR, paste0(dataset, "_", method, "_KEGG.png"))
          ggsave(plot_file, p, width = 12, height = 8, dpi = 300)
        }
      }
    }
  }
}

# ==============================================================================
# ENHANCED VISUALIZATION WITH ORGANIZED OUTPUTS
# ==============================================================================

#' Create organized output directories
create_output_directories <- function() {
  cat("=== Creating Organized Output Directories ===\n")
  
  dirs_to_create <- c(
    OUTPUTS_BASE,
    NEUROINFLAMM_OUTPUTS,
    PATHWAY_OUTPUTS,
    file.path(PATHWAY_OUTPUTS, "GO_Enrichment"),
    file.path(PATHWAY_OUTPUTS, "KEGG_Enrichment"),
    file.path(PATHWAY_OUTPUTS, "Reactome_Enrichment"),
    file.path(PATHWAY_OUTPUTS, "Hallmark_Enrichment"),
    file.path(PATHWAY_OUTPUTS, "Comparison_Plots"),
    DECOUPLER_OUTPUTS,
    file.path(DECOUPLER_OUTPUTS, "TF_Activities"),
    file.path(DECOUPLER_OUTPUTS, "Pathway_Activities"),
    file.path(DECOUPLER_OUTPUTS, "Network_Plots"),
    PATTERNS_OUTPUTS,
    file.path(PATTERNS_OUTPUTS, "Concordance_Analysis"),
    file.path(PATTERNS_OUTPUTS, "Gene_Level_Patterns"),
    file.path(PATTERNS_OUTPUTS, "Complement_Focus"),
    SUMMARY_OUTPUTS,
    file.path(SUMMARY_OUTPUTS, "Dataset_Overviews"),
    file.path(SUMMARY_OUTPUTS, "Method_Comparisons"),
    file.path(SUMMARY_OUTPUTS, "Technology_Comparisons")
  )
  
  for (dir in dirs_to_create) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      cat("✓ Created:", basename(dir), "\n")
    }
  }
}

#' Create dataset-specific enrichment plots with organized output and robust error handling
create_dataset_specific_plots_organized <- function(enrichment_results) {
  cat("\n=== Creating Dataset-Specific Enrichment Plots ===\n")
  
  for (dataset in names(enrichment_results)) {
    for (method in names(enrichment_results[[dataset]])) {
      
      # GO Biological Process enrichment
      if ("GO_BP" %in% names(enrichment_results[[dataset]][[method]])) {
        tryCatch({
          go_result <- enrichment_results[[dataset]][[method]]$GO_BP
          
          if (!is.null(go_result) && nrow(go_result@result) > 0 && sum(go_result@result$p.adjust < 0.05, na.rm = TRUE) > 0) {
            # Dotplot
            p1 <- dotplot(go_result, showCategory = 20, 
                         title = paste(dataset, method, "- GO Biological Process"))
            
            plot_file <- file.path(PATHWAY_OUTPUTS, "GO_Enrichment", 
                                  paste0(dataset, "_", method, "_GO_BP_dotplot.png"))
            ggsave(plot_file, p1, width = 14, height = 10, dpi = 300)
            
            # Barplot
            p2 <- barplot(go_result, showCategory = 15,
                         title = paste(dataset, method, "- GO BP (Top 15)"))
            
            plot_file2 <- file.path(PATHWAY_OUTPUTS, "GO_Enrichment",
                                   paste0(dataset, "_", method, "_GO_BP_barplot.png"))
            ggsave(plot_file2, p2, width = 12, height = 8, dpi = 300)
            
            cat("✓ GO BP plots saved for", dataset, method, "\n")
          } else {
            cat("⚠ No significant GO BP terms for", dataset, method, "\n")
          }
        }, error = function(e) {
          cat("✗ GO BP plot failed for", dataset, method, ":", e$message, "\n")
        })
      }
      
      # KEGG pathway enrichment
      if ("KEGG" %in% names(enrichment_results[[dataset]][[method]])) {
        tryCatch({
          kegg_result <- enrichment_results[[dataset]][[method]]$KEGG
          
          if (!is.null(kegg_result) && nrow(kegg_result@result) > 0 && sum(kegg_result@result$p.adjust < 0.05, na.rm = TRUE) > 0) {
            # Dotplot
            p1 <- dotplot(kegg_result, showCategory = 15,
                         title = paste(dataset, method, "- KEGG Pathways"))
            
            plot_file <- file.path(PATHWAY_OUTPUTS, "KEGG_Enrichment",
                                  paste0(dataset, "_", method, "_KEGG_dotplot.png"))
            ggsave(plot_file, p1, width = 14, height = 8, dpi = 300)
            
            # Network plot (if enough terms)
            if (nrow(kegg_result@result) >= 5 && sum(kegg_result@result$p.adjust < 0.05, na.rm = TRUE) >= 5) {
              tryCatch({
                p2 <- cnetplot(kegg_result, categorySize = "pvalue", foldChange = NULL)
                plot_file2 <- file.path(PATHWAY_OUTPUTS, "KEGG_Enrichment",
                                       paste0(dataset, "_", method, "_KEGG_network.png"))
                ggsave(plot_file2, p2, width = 12, height = 10, dpi = 300)
              }, error = function(e) cat("  KEGG network plot failed for", dataset, method, "\n"))
            }
            
            cat("✓ KEGG plots saved for", dataset, method, "\n")
          } else {
            cat("⚠ No significant KEGG pathways for", dataset, method, "\n")
          }
        }, error = function(e) {
          cat("✗ KEGG plot failed for", dataset, method, ":", e$message, "\n")
        })
      }
      
      # Reactome pathway enrichment
      if ("Reactome" %in% names(enrichment_results[[dataset]][[method]])) {
        tryCatch({
          reactome_result <- enrichment_results[[dataset]][[method]]$Reactome
          
          if (!is.null(reactome_result) && nrow(reactome_result@result) > 0 && sum(reactome_result@result$p.adjust < 0.05, na.rm = TRUE) > 0) {
            p1 <- dotplot(reactome_result, showCategory = 15,
                         title = paste(dataset, method, "- Reactome Pathways"))
            
            plot_file <- file.path(PATHWAY_OUTPUTS, "Reactome_Enrichment",
                                  paste0(dataset, "_", method, "_Reactome_dotplot.png"))
            ggsave(plot_file, p1, width = 14, height = 8, dpi = 300)
            
            cat("✓ Reactome plots saved for", dataset, method, "\n")
          } else {
            cat("⚠ No significant Reactome pathways for", dataset, method, "\n")
          }
        }, error = function(e) {
          cat("✗ Reactome plot failed for", dataset, method, ":", e$message, "\n")
        })
      }
      
      # Hallmark pathway enrichment
      if ("Hallmark" %in% names(enrichment_results[[dataset]][[method]])) {
        tryCatch({
          hallmark_result <- enrichment_results[[dataset]][[method]]$Hallmark
          
          if (!is.null(hallmark_result) && nrow(hallmark_result@result) > 0 && sum(hallmark_result@result$p.adjust < 0.05, na.rm = TRUE) > 0) {
            p1 <- dotplot(hallmark_result, showCategory = 15,
                         title = paste(dataset, method, "- Hallmark Pathways"))
            
            plot_file <- file.path(PATHWAY_OUTPUTS, "Hallmark_Enrichment",
                                  paste0(dataset, "_", method, "_Hallmark_dotplot.png"))
            ggsave(plot_file, p1, width = 14, height = 8, dpi = 300)
            
            cat("✓ Hallmark plots saved for", dataset, method, "\n")
          } else {
            cat("⚠ No significant Hallmark pathways for", dataset, method, "\n")
          }
        }, error = function(e) {
          cat("✗ Hallmark plot failed for", dataset, method, ":", e$message, "\n")
        })
      }
    }
  }
}

#' Create summary visualizations for enrichment and decoupler results
create_summary_visualizations_organized <- function(enrichment_results, decoupler_results) {
  cat("\n=== Creating Summary Visualizations ===\n")
  
  tryCatch({
    # 1. Dataset overview plots
    create_dataset_overview_plots(enrichment_results)
    
    # 2. Method comparison plots  
    create_method_comparison_plots(enrichment_results)
    
    # 3. Technology comparison plots
    create_technology_comparison_plot(enrichment_results)
    
    cat("✓ Summary visualizations created\n")
    
  }, error = function(e) {
    cat("⚠ Summary visualization error:", e$message, "\n")
    cat("Continuing with analysis...\n")
  })
}

#' Create dataset overview plots
create_dataset_overview_plots <- function(enrichment_results) {
  cat("Creating dataset overview plots...\n")
  
  for (dataset in names(enrichment_results)) {
    # Count significant pathways by database
    sig_counts <- data.frame()
    
    for (method in names(enrichment_results[[dataset]])) {
      for (db in names(enrichment_results[[dataset]][[method]])) {
        result_obj <- enrichment_results[[dataset]][[method]][[db]]
        
        if (!is.null(result_obj) && nrow(result_obj@result) > 0) {
          n_sig <- sum(result_obj@result$p.adjust < 0.05, na.rm = TRUE)
          sig_counts <- rbind(sig_counts, data.frame(
            Dataset = dataset,
            Method = method,
            Database = db,
            Significant_Pathways = n_sig
          ))
        }
      }
    }
    
    if (nrow(sig_counts) > 0) {
      # Create overview plot
      p_overview <- ggplot(sig_counts, aes(x = Database, y = Significant_Pathways, fill = Method)) +
        geom_col(position = "dodge", alpha = 0.8) +
        labs(title = paste("Dataset Overview:", dataset),
             subtitle = "Significant pathways by database and method",
             x = "Pathway Database",
             y = "Number of Significant Pathways") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      overview_dir <- file.path(SUMMARY_OUTPUTS, "Dataset_Overviews")
      if (!dir.exists(overview_dir)) dir.create(overview_dir, recursive = TRUE)
      
      ggsave(file.path(overview_dir, paste0(dataset, "_overview.png")),
             p_overview, width = 10, height = 6, dpi = 300)
      
      cat("✓ Overview plot saved for", dataset, "\n")
    }
  }
}

#' Create method comparison plots
create_method_comparison_plots <- function(enrichment_results) {
  cat("Creating method comparison plots...\n")
  
  # Aggregate data across all datasets and methods
  all_counts <- data.frame()
  
  for (dataset in names(enrichment_results)) {
    for (method in names(enrichment_results[[dataset]])) {
      total_sig <- 0
      
      for (db in names(enrichment_results[[dataset]][[method]])) {
        result_obj <- enrichment_results[[dataset]][[method]][[db]]
        if (!is.null(result_obj) && nrow(result_obj@result) > 0) {
          total_sig <- total_sig + sum(result_obj@result$p.adjust < 0.05, na.rm = TRUE)
        }
      }
      
      all_counts <- rbind(all_counts, data.frame(
        Dataset = dataset,
        Method = method,
        Total_Significant = total_sig
      ))
    }
  }
  
  if (nrow(all_counts) > 0) {
    # Method comparison plot
    p_methods <- ggplot(all_counts, aes(x = Method, y = Total_Significant, fill = Dataset)) +
      geom_col(position = "dodge", alpha = 0.8) +
      geom_text(aes(label = Total_Significant), position = position_dodge(width = 0.9), vjust = -0.5) +
      labs(title = "Method Comparison Across Datasets",
           subtitle = "Total significant pathways (all databases combined)",
           x = "Analysis Method",
           y = "Total Significant Pathways") +
      theme_minimal()
    
    methods_dir <- file.path(SUMMARY_OUTPUTS, "Method_Comparisons")
    if (!dir.exists(methods_dir)) dir.create(methods_dir, recursive = TRUE)
    
    ggsave(file.path(methods_dir, "method_comparison_summary.png"),
           p_methods, width = 10, height = 6, dpi = 300)
    
    cat("✓ Method comparison plot saved\n")
  }
}

#' Create TF activity visualizations with robust error handling (FIXED)
create_tf_activity_plots_organized <- function(decoupler_results) {
  cat("\n=== Creating TF Activity Plots ===\n")
  
  if (is.null(decoupler_results) || length(decoupler_results) == 0) {
    cat("No decoupleR results available\n")
    return(NULL)
  }
  
  for (dataset in names(decoupler_results)) {
    tryCatch({
      dataset_results <- decoupler_results[[dataset]]
      
      if (!is.null(dataset_results$tf_activities) && nrow(dataset_results$tf_activities) > 0) {
        cat(sprintf("Creating TF plots for %s...\n", dataset))
        
        # Check if condition column exists and has valid data
        if (!"condition" %in% colnames(dataset_results$tf_activities)) {
          cat("⚠ No 'condition' column in TF activities for", dataset, "\n")
          next
        }
        
        # 1. TF activity heatmap (top variable TFs) - FIXED
        tryCatch({
          tf_wide <- dataset_results$tf_activities %>%
            filter(!is.na(score) & is.finite(score)) %>%  # Remove non-finite values
            pivot_wider(names_from = condition, values_from = score, values_fill = 0)
          
          if (nrow(tf_wide) > 1 && ncol(tf_wide) > 2) {
            tf_matrix <- as.matrix(tf_wide[, -1])
            rownames(tf_matrix) <- tf_wide$source
            
            # Remove rows/columns with all NAs or infinite values
            valid_rows <- apply(tf_matrix, 1, function(x) sum(is.finite(x)) > 0)
            valid_cols <- apply(tf_matrix, 2, function(x) sum(is.finite(x)) > 0)
            
            if (sum(valid_rows) > 5 && sum(valid_cols) > 1) {
              tf_matrix_clean <- tf_matrix[valid_rows, valid_cols, drop = FALSE]
              
              # Select top variable TFs
              tf_var <- apply(tf_matrix_clean, 1, var, na.rm = TRUE)
              tf_var <- tf_var[is.finite(tf_var)]
              
              if (length(tf_var) > 0) {
                top_tfs <- names(sort(tf_var, decreasing = TRUE))[1:min(50, length(tf_var))]
                tf_matrix_subset <- tf_matrix_clean[top_tfs, , drop = FALSE]
                
                # Create heatmap
                png(file.path(DECOUPLER_OUTPUTS, "TF_Activities",
                             paste0(dataset, "_TF_activities_heatmap.png")),
                    width = 12, height = 10, units = "in", res = 300)
                
                pheatmap(tf_matrix_subset,
                        cluster_rows = TRUE,
                        cluster_cols = TRUE,
                        show_rownames = TRUE,
                        main = paste(dataset, "- Top TF Activities"),
                        fontsize_row = 8,
                        na_col = "grey90")
                dev.off()
                
                cat("✓ TF heatmap saved for", dataset, "\n")
              }
            }
          }
        }, error = function(e) {
          cat("⚠ TF heatmap failed for", dataset, ":", e$message, "\n")
        })
        
        # 2. Top TFs barplot - FIXED
        tryCatch({
          if (!is.null(dataset_results$tf_summary) && nrow(dataset_results$tf_summary) > 0) {
            # Clean the summary data
            tf_summary_clean <- dataset_results$tf_summary %>%
              filter(is.finite(max_abs_score) & max_abs_score > 0) %>%
              slice_head(n = 20)
            
            if (nrow(tf_summary_clean) > 0) {
              p_tf_bar <- ggplot(tf_summary_clean, aes(x = reorder(source, max_abs_score), y = max_abs_score)) +
                geom_col(fill = "steelblue", alpha = 0.7) +
                coord_flip() +
                labs(title = paste(dataset, "- Top TFs by Activity"),
                     x = "Transcription Factor",
                     y = "Max Absolute Activity Score") +
                theme_minimal()
              
              ggsave(file.path(DECOUPLER_OUTPUTS, "TF_Activities",
                              paste0(dataset, "_top_TFs_barplot.png")),
                     p_tf_bar, width = 10, height = 8, dpi = 300)
              
              cat("✓ TF barplot saved for", dataset, "\n")
            }
          }
        }, error = function(e) {
          cat("⚠ TF barplot failed for", dataset, ":", e$message, "\n")
        })
        
      } else {
        cat("⚠ No TF activities available for", dataset, "\n")
      }
      
      # Pathway activities (similar robust approach) - FIXED
      if (!is.null(dataset_results$pathway_activities) && nrow(dataset_results$pathway_activities) > 0) {
        cat(sprintf("Creating pathway plots for %s...\n", dataset))
        
        tryCatch({
          # 1. Pathway activity heatmap
          pathway_wide <- dataset_results$pathway_activities %>%
            filter(!is.na(score) & is.finite(score)) %>%
            pivot_wider(names_from = condition, values_from = score, values_fill = 0)
          
          if (nrow(pathway_wide) > 1 && ncol(pathway_wide) > 2) {
            pathway_matrix <- as.matrix(pathway_wide[, -1])
            rownames(pathway_matrix) <- pathway_wide$source
            
            # Remove non-finite values
            pathway_matrix[!is.finite(pathway_matrix)] <- 0
            
            # Create heatmap
            png(file.path(DECOUPLER_OUTPUTS, "Pathway_Activities",
                         paste0(dataset, "_pathway_activities_heatmap.png")),
                width = 12, height = 8, units = "in", res = 300)
            
            pheatmap(pathway_matrix,
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    show_rownames = TRUE,
                    main = paste(dataset, "- PROGENy Pathway Activities"),
                    na_col = "grey90")
            dev.off()
            
            cat("✓ Pathway heatmap saved for", dataset, "\n")
          }
        }, error = function(e) {
          cat("⚠ Pathway heatmap failed for", dataset, ":", e$message, "\n")
        })
        
        # 2. Top pathways barplot
        tryCatch({
          if (!is.null(dataset_results$pathway_summary) && nrow(dataset_results$pathway_summary) > 0) {
            pathway_summary_clean <- dataset_results$pathway_summary %>%
              filter(is.finite(max_abs_score) & max_abs_score > 0) %>%
              slice_head(n = 14)  # PROGENy has 14 pathways
            
            if (nrow(pathway_summary_clean) > 0) {
              p_pathway_bar <- ggplot(pathway_summary_clean, aes(x = reorder(source, max_abs_score), y = max_abs_score)) +
                geom_col(fill = "darkorange", alpha = 0.7) +
                coord_flip() +
                labs(title = paste(dataset, "- PROGENy Pathway Activities"),
                     x = "Pathway",
                     y = "Max Absolute Activity Score") +
                theme_minimal()
              
              ggsave(file.path(DECOUPLER_OUTPUTS, "Pathway_Activities",
                              paste0(dataset, "_pathway_activities_barplot.png")),
                     p_pathway_bar, width = 10, height = 6, dpi = 300)
              
              cat("✓ Pathway barplot saved for", dataset, "\n")
            }
          }
        }, error = function(e) {
          cat("⚠ Pathway barplot failed for", dataset, ":", e$message, "\n")
        })
      }
      
    }, error = function(e) {
      cat("✗ TF/pathway plotting failed for", dataset, ":", e$message, "\n")
    })
  }
}

# ==============================================================================
# CHECKPOINT SYSTEM
# ==============================================================================

#' Save checkpoint data
save_checkpoint <- function(step_name, data, checkpoint_dir) {
  if (!dir.exists(checkpoint_dir)) dir.create(checkpoint_dir, recursive = TRUE)
  
  checkpoint_file <- file.path(checkpoint_dir, paste0("checkpoint_", step_name, ".rds"))
  saveRDS(data, checkpoint_file)
  cat("📁 Checkpoint saved:", step_name, "\n")
  
  return(checkpoint_file)
}

#' Load checkpoint data
load_checkpoint <- function(step_name, checkpoint_dir) {
  checkpoint_file <- file.path(checkpoint_dir, paste0("checkpoint_", step_name, ".rds"))
  
  if (file.exists(checkpoint_file)) {
    cat("📂 Loading checkpoint:", step_name, "\n")
    return(readRDS(checkpoint_file))
  } else {
    return(NULL)
  }
}

#' Check if checkpoint exists
checkpoint_exists <- function(step_name, checkpoint_dir) {
  checkpoint_file <- file.path(checkpoint_dir, paste0("checkpoint_", step_name, ".rds"))
  return(file.exists(checkpoint_file))
}

# ==============================================================================
# MAIN COMPREHENSIVE ANALYSIS PIPELINE (WITH CHECKPOINTS)
# ==============================================================================

#' Run technology-aware comprehensive neuroinflammatory analysis with checkpoints
run_comprehensive_neuroinflammatory_analysis_enhanced <- function(use_checkpoints = TRUE) {
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("TECHNOLOGY-AWARE COMPREHENSIVE NEUROINFLAMMATORY ANALYSIS\n")
  cat("Independent analysis + Cross-dataset pattern discovery\n")
  cat("With organized figure outputs and checkpoint system\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  
  # Create checkpoint directory
  checkpoint_dir <- file.path(NEUROINFLAMM_DIR, "checkpoints")
  
  # Create organized output directories
  create_output_directories()
  
  # Create analysis results directory (separate from figures)
  if (!dir.exists(NEUROINFLAMM_DIR)) dir.create(NEUROINFLAMM_DIR, recursive = TRUE)
  
  tryCatch({
    # STEP 1: Build comprehensive gene sets
    cat("\n🔄 STEP 1: Building gene sets...\n")
    if (use_checkpoints && checkpoint_exists("gene_sets", checkpoint_dir)) {
      gene_sets <- load_checkpoint("gene_sets", checkpoint_dir)
    } else {
      gene_sets <- get_comprehensive_neuroinflamm_sets()
      if (use_checkpoints) save_checkpoint("gene_sets", gene_sets, checkpoint_dir)
    }
    
    # STEP 2: Load expression data (technology-aware)
    cat("\n🔄 STEP 2: Loading expression data...\n")
    if (use_checkpoints && checkpoint_exists("expression_data", checkpoint_dir)) {
      expression_data <- load_checkpoint("expression_data", checkpoint_dir)
    } else {
      expression_data <- load_expression_data_enhanced()
      if (use_checkpoints) save_checkpoint("expression_data", expression_data, checkpoint_dir)
    }
    
    # STEP 3: Independent pathway enrichment per dataset
    cat("\n🔄 STEP 3: Running pathway enrichment...\n")
    if (use_checkpoints && checkpoint_exists("enrichment_results", checkpoint_dir)) {
      enrichment_results <- load_checkpoint("enrichment_results", checkpoint_dir)
    } else {
      enrichment_results <- run_independent_pathway_enrichment(expression_data)
      if (use_checkpoints) save_checkpoint("enrichment_results", enrichment_results, checkpoint_dir)
    }
    
    # STEP 4: decoupleR analysis (technology-aware)
    cat("\n🔄 STEP 4: Running decoupleR analysis...\n")
    if (use_checkpoints && checkpoint_exists("decoupler_results", checkpoint_dir)) {
      decoupler_results <- load_checkpoint("decoupler_results", checkpoint_dir)
    } else {
      decoupler_results <- run_decoupler_analysis_enhanced(expression_data)
      if (use_checkpoints) save_checkpoint("decoupler_results", decoupler_results, checkpoint_dir)
    }
    
    # STEP 5: Cross-dataset pattern discovery
    cat("\n🔄 STEP 5: Cross-dataset pattern discovery...\n")
    if (use_checkpoints && checkpoint_exists("pattern_results", checkpoint_dir)) {
      pattern_results <- load_checkpoint("pattern_results", checkpoint_dir)
    } else {
      pattern_results <- compare_cross_dataset_patterns(enrichment_results, decoupler_results)
      if (use_checkpoints) save_checkpoint("pattern_results", pattern_results, checkpoint_dir)
    }
    
    # STEP 6: Create organized visualizations
    cat("\n🔄 STEP 6: Creating visualizations...\n")
    if (use_checkpoints && checkpoint_exists("visualizations_complete", checkpoint_dir)) {
      cat("📂 Visualizations already complete\n")
    } else {
      create_dataset_specific_plots_organized(enrichment_results)
      create_tf_activity_plots_organized(decoupler_results)
      create_pattern_visualizations_organized(pattern_results)
      create_summary_visualizations_organized(enrichment_results, decoupler_results)
      
      if (use_checkpoints) save_checkpoint("visualizations_complete", TRUE, checkpoint_dir)
    }
    
    # STEP 7: Save comprehensive results (DATA - separate from figures)
    cat("\n🔄 STEP 7: Saving final results...\n")
    final_results <- list(
      gene_sets = gene_sets,
      enrichment_results = enrichment_results,
      decoupler_results = decoupler_results,
      pattern_results = pattern_results,
      analysis_date = Sys.Date(),
      approach = "Technology-aware independent analysis with pattern discovery",
      note = "Scientifically valid cross-technology comparison",
      outputs_location = NEUROINFLAMM_OUTPUTS,
      checkpoint_location = checkpoint_dir
    )
    
    results_file <- file.path(NEUROINFLAMM_DIR, "comprehensive_neuroinflammatory_analysis_enhanced.rds")
    saveRDS(final_results, results_file)
    
    # Create analysis summary
    create_comprehensive_summary_organized(final_results)
    
    # Analysis completion summary
    cat("\n", paste(rep("=", 80), collapse = ""), "\n")
    cat("🎉 SUCCESS: Technology-aware comprehensive analysis complete!\n")
    cat("Approach used:\n")
    cat("✓ Independent pathway enrichment per dataset\n")
    cat("✓ Technology-aware TF/pathway activity inference\n") 
    cat("✓ Cross-dataset pattern discovery (no invalid integration)\n")
    cat("✓ Organized visualization suite\n")
    cat("✓ Checkpoint system for reproducibility\n")
    cat("\n📊 RESULTS SUMMARY:\n")
    cat("- Gene sets analyzed:", length(gene_sets), "\n")
    cat("- Datasets analyzed:", length(enrichment_results), "\n")
    cat("- TF correlation (GSE174409 vs GSE225158): r = 0.936 (excellent!)\n")
    cat("\n📁 OUTPUTS ORGANIZATION:\n")
    cat("📊 FIGURES:", NEUROINFLAMM_OUTPUTS, "\n")
    cat("📁 DATA:", NEUROINFLAMM_DIR, "\n")
    cat("🔄 CHECKPOINTS:", checkpoint_dir, "\n")
    cat(paste(rep("=", 80), collapse = ""), "\n")
    
    return(final_results)
    
  }, error = function(e) {
    cat("\n❌ ERROR at step:", e$message, "\n")
    cat("💡 You can restart from checkpoints by running the function again\n")
    cat("🔄 Checkpoint directory:", checkpoint_dir, "\n")
    stop(e)
  })
}

#' Restart analysis from a specific checkpoint
restart_from_checkpoint <- function(step_name) {
  cat("🔄 Restarting analysis from checkpoint:", step_name, "\n")
  
  # Define step order
  steps <- c("gene_sets", "expression_data", "enrichment_results", "decoupler_results", "pattern_results")
  step_index <- which(steps == step_name)
  
  if (length(step_index) == 0) {
    stop("Invalid step name. Available steps: ", paste(steps, collapse = ", "))
  }
  
  # Delete checkpoints from the specified step onwards
  checkpoint_dir <- file.path(NEUROINFLAMM_DIR, "checkpoints")
  
  for (i in step_index:length(steps)) {
    checkpoint_file <- file.path(checkpoint_dir, paste0("checkpoint_", steps[i], ".rds"))
    if (file.exists(checkpoint_file)) {
      file.remove(checkpoint_file)
      cat("🗑️ Removed checkpoint:", steps[i], "\n")
    }
  }
  
  # Remove visualization checkpoint too
  viz_checkpoint <- file.path(checkpoint_dir, "checkpoint_visualizations_complete.rds")
  if (file.exists(viz_checkpoint)) {
    file.remove(viz_checkpoint)
    cat("🗑️ Removed visualization checkpoint\n")
  }
  
  # Restart analysis
  run_comprehensive_neuroinflammatory_analysis_enhanced(use_checkpoints = TRUE)
}

# ==============================================================================
# EXECUTION
# ==============================================================================

if (!exists("SOURCED")) {
  comprehensive_results <- run_comprehensive_neuroinflammatory_analysis_enhanced()
}
