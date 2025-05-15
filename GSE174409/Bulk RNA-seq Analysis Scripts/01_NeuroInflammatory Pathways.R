# 01_NeuroInflammatory Pathways.R
# Analysis of neuroinflammation-related genes in the OUD RNA-seq dataset using KEGG pathways

# ----------------------
# 1. SETUP AND DATA LOADING
# ----------------------

# Load required libraries
suppressMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(pheatmap)
  library(dplyr)
  library(RColorBrewer)
  library(biomaRt)    # For Ensembl ID mapping
  library(reshape2)   # For data reshaping
  library(tidyr)      # For separate function
  library(KEGGREST)   # For KEGG pathway queries
  library(plyr)       # For rbind.fill
  library(clusterProfiler) # For pathway enrichment
})

# Set directories
base_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE174409"
qc_dir <- file.path(base_dir, "QC")
output_dir <- file.path(base_dir, "NeuroinflammationResults")
dir.create(output_dir, showWarnings = FALSE)

# Load preprocessed data
dge_filtered <- readRDS(file.path(qc_dir, "dge_filtered_normalized.rds"))
metadata_df <- readRDS(file.path(qc_dir, "metadata.rds"))
logcpm_filtered_norm <- readRDS(file.path(qc_dir, "logcpm_filtered_normalized.rds"))

# ----------------------
# 2. MAP ENSEMBL IDS TO GENE SYMBOLS
# ----------------------

# Get Ensembl IDs from data
ensembl_ids <- rownames(logcpm_filtered_norm)
cat("First few gene IDs in dataset:", head(ensembl_ids), "\n")

# Set up biomaRt to convert Ensembl IDs to gene symbols
cat("Converting Ensembl IDs to gene symbols using biomaRt...\n")
ensembl <- tryCatch({
  useMart("ensembl", dataset = "hsapiens_gene_ensembl")
}, error = function(e) {
  cat("Error connecting to Ensembl: ", e$message, "\n")
  cat("Using archived version from April 2023...\n")
  useMart("ENSEMBL_MART_ENSEMBL", host="https://apr2023.archive.ensembl.org",
          dataset="hsapiens_gene_ensembl")
})

mapping_result <- tryCatch({
  getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id", "description"),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = ensembl
  )
}, error = function(e) {
  cat("Error retrieving data from biomaRt: ", e$message, "\n")
  # Return empty data frame with correct structure
  data.frame(
    ensembl_gene_id = character(0),
    external_gene_name = character(0),
    entrezgene_id = character(0),
    description = character(0),
    stringsAsFactors = FALSE
  )
})

# Create mapping dictionary
gene_symbols <- setNames(
  mapping_result$external_gene_name,
  mapping_result$ensembl_gene_id
)

# Count how many IDs were mapped
mapped_count <- sum(ensembl_ids %in% names(gene_symbols))
cat("Successfully mapped", mapped_count, "out of", length(ensembl_ids), "Ensembl IDs to gene symbols\n")

# Save the mapping for reference
mapping_df <- data.frame(
  EnsemblID = ensembl_ids,
  Symbol = gene_symbols[ensembl_ids],
  stringsAsFactors = FALSE
)
write.csv(mapping_df, file.path(output_dir, "ensembl_to_symbol_mapping.csv"), row.names = FALSE)

# ----------------------
# 3. IDENTIFY NEUROINFLAMMATION GENES USING KEGG PATHWAYS
# ----------------------

# Define relevant KEGG pathways for neuroinflammation
neuro_pathways <- c(
  "hsa04060", # Cytokine-cytokine receptor interaction
  "hsa04061", # Viral protein interaction with cytokine and cytokine receptor
  "hsa04062", # Chemokine signaling pathway
  "hsa04620", # Toll-like receptor signaling pathway
  "hsa04621", # NOD-like receptor signaling pathway
  "hsa04622", # RIG-I-like receptor signaling pathway
  "hsa04623", # Cytosolic DNA-sensing pathway
  "hsa04625", # C-type lectin receptor signaling pathway
  "hsa04650", # Natural killer cell mediated cytotoxicity
  "hsa04657", # IL-17 signaling pathway
  "hsa04658", # Th1 and Th2 cell differentiation
  "hsa04659", # Th17 cell differentiation
  "hsa04660", # T cell receptor signaling pathway
  "hsa04662", # B cell receptor signaling pathway
  "hsa04668", # TNF signaling pathway
  "hsa04672", # Intestinal immune network for IgA production
  "hsa04610", # Complement and coagulation cascades
  "hsa04640", # Hematopoietic cell lineage
  "hsa04380"  # Osteoclast differentiation (includes inflammation mediators)
)

# Retrieve genes from KEGG pathways
cat("Retrieving genes from KEGG pathways...\n")
neuro_genes_kegg <- list()

for (pathway in neuro_pathways) {
  tryCatch({
    # Get pathway information with timeout
    pathway_info <- NULL
    attempts <- 0
    max_attempts <- 3

    while(is.null(pathway_info) && attempts < max_attempts) {
      attempts <- attempts + 1
      pathway_info <- tryCatch({
        keggGet(pathway)
      }, error = function(e) {
        cat("Attempt", attempts, "- Error retrieving pathway", pathway, ":", e$message, "\n")
        if(attempts < max_attempts) Sys.sleep(2)  # Wait before retrying
        return(NULL)
      })
    }

    if(is.null(pathway_info)) {
      cat("Failed to retrieve pathway", pathway, "after", max_attempts, "attempts. Skipping.\n")
      next
    }

    if(length(pathway_info) > 0 && "GENE" %in% names(pathway_info[[1]])) {
      # Extract gene information
      genes <- pathway_info[[1]]$GENE
      # KEGG returns a vector where even indices are gene IDs and odd indices are descriptions
      gene_ids <- genes[seq(1, length(genes), 2)]
      gene_symbols <- sapply(strsplit(genes[seq(2, length(genes), 2)], ";"), `[`, 1)
      gene_symbols <- trimws(gene_symbols)

      cat("Retrieved", length(gene_ids), "genes from pathway", pathway, "\n")

      # Store in list
      neuro_genes_kegg[[pathway]] <- data.frame(
        EntrezID = gene_ids,
        Symbol = gene_symbols,
        Pathway = pathway,
        stringsAsFactors = FALSE
      )
    } else {
      cat("No genes found in pathway", pathway, "\n")
    }
  }, error = function(e) {
    cat("Error processing pathway", pathway, ":", e$message, "\n")
  })
}

# Check if we retrieved any pathways
if(length(neuro_genes_kegg) == 0) {
  cat("WARNING: Failed to retrieve any pathways from KEGG. Using backup list of neuroinflammation genes.\n")

  # Fallback to a defined list of neuroinflammation genes
  neuro_symbols <- c(
    # Cytokines and inflammation mediators
    "IL1B", "IL1A", "IL6", "TNF", "IL10", "TGFB1", "IL18", "IFNG", "IL4", "IL13",
    # Receptors and signaling
    "TLR2", "TLR4", "NLRP3", "CASP1", "TNFRSF1A", "TNFRSF1B",
    "STAT1", "STAT3", "JAK1", "JAK2", "RELA", "NFKB1",
    # Glia markers
    "GFAP", "S100B", "AIF1", "P2RY12", "TMEM119", "CX3CR1", "CX3CL1",
    "CD68", "CD14", "CD163", "MRC1", "APOE",
    # Complement components
    "C1QA", "C1QB", "C1QC", "C3", "C4A", "C4B", "C1R", "C1S",
    "CFB", "CFH", "CFI", "C5AR1", "CR1", "VSIG4"
  )

  # Map symbols to Ensembl IDs using our existing mapping
  symbol_to_ensembl <- setNames(
    mapping_df$EnsemblID,
    mapping_df$Symbol
  )

  # Find genes in our dataset
  neuro_ensembl <- symbol_to_ensembl[neuro_symbols]
  neuro_ensembl <- neuro_ensembl[!is.na(neuro_ensembl)]

  # Create mapping
  neuro_mapping <- data.frame(
    EnsemblID = neuro_ensembl,
    Symbol = names(neuro_ensembl),
    Source = "Manual list",
    stringsAsFactors = FALSE
  )
} else {
  # Combine all pathway genes
  all_pathway_genes <- rbind.fill(neuro_genes_kegg)
  all_pathway_genes <- all_pathway_genes[!duplicated(all_pathway_genes$Symbol), ]

  # Get Entrez to Ensembl mapping
  entrez_to_ensembl <- tryCatch({
    getBM(
      attributes = c("entrezgene_id", "ensembl_gene_id", "external_gene_name"),
      filters = "entrezgene_id",
      values = all_pathway_genes$EntrezID,
      mart = ensembl
    )
  }, error = function(e) {
    cat("Error retrieving Entrez to Ensembl mapping: ", e$message, "\n")
    # Try direct symbol mapping as fallback
    getBM(
      attributes = c("external_gene_name", "ensembl_gene_id"),
      filters = "external_gene_name",
      values = all_pathway_genes$Symbol,
      mart = ensembl
    )
  })

  # Find neuroinflammation genes in our dataset
  if("entrezgene_id" %in% colnames(entrez_to_ensembl)) {
    # Use Entrez ID mapping
    neuro_ensembl <- entrez_to_ensembl$ensembl_gene_id[entrez_to_ensembl$ensembl_gene_id %in% ensembl_ids]
  } else {
    # Use symbol mapping
    neuro_ensembl <- entrez_to_ensembl$ensembl_gene_id[entrez_to_ensembl$ensembl_gene_id %in% ensembl_ids]
  }

  neuro_ensembl <- unique(neuro_ensembl)

  # Create a mapping for analysis with pathway information
  neuro_mapping <- data.frame(
    EnsemblID = neuro_ensembl,
    Symbol = mapping_df$Symbol[match(neuro_ensembl, mapping_df$EnsemblID)],
    Source = "KEGG pathways",
    stringsAsFactors = FALSE
  )

  # Save pathway details for reference
  pathway_details <- data.frame(
    PathwayID = rep(names(neuro_genes_kegg), sapply(neuro_genes_kegg, nrow)),
    rbind.fill(neuro_genes_kegg)
  )
  write.csv(pathway_details, file.path(output_dir, "kegg_pathway_details.csv"), row.names = FALSE)
}

cat("Found", length(neuro_ensembl), "neuroinflammation-related genes in dataset\n")

# Save neuroinflammation gene list
write.csv(neuro_mapping, file.path(output_dir, "neuroinflammation_gene_list.csv"), row.names = FALSE)

# ----------------------
# 4. EXPLORATORY ANALYSIS
# ----------------------

if(length(neuro_ensembl) > 0) {
  # Extract expression data for neuroinflammation genes
  neuro_expr <- logcpm_filtered_norm[neuro_ensembl, ]

  # Create a version with gene symbols as row names for visualization
  neuro_expr_symbols <- neuro_expr
  rownames(neuro_expr_symbols) <- neuro_mapping$Symbol

  # Save expression data
  write.csv(neuro_expr_symbols, file.path(output_dir, "neuroinflammation_expression.csv"))

  # Create annotation for heatmap
  anno <- data.frame(
    Diagnosis = metadata_df$diagnosis,
    Region = metadata_df$region,
    row.names = colnames(neuro_expr)
  )

  # Color scheme
  anno_colors <- list(
    Diagnosis = c("CONT" = "blue", "OUD" = "red"),
    Region = c("NAC" = "purple", "DLPFC" = "green")
  )

  # Generate heatmap with gene symbols
  tryCatch({
    png(file.path(output_dir, "neuroinflammation_heatmap.png"), width = 1200, height = 1600, res = 120)
    pheatmap(
      neuro_expr_symbols,
      annotation_col = anno,
      annotation_colors = anno_colors,
      scale = "row",
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_colnames = FALSE,
      fontsize_row = 10,
      main = "Neuroinflammation-Related Gene Expression"
    )
    dev.off()
  }, error = function(e) {
    cat("Error generating heatmap:", e$message, "\n")
  })

  # PCA analysis
  tryCatch({
    pca_neuro <- prcomp(t(neuro_expr))
    var_explained <- (pca_neuro$sdev^2) / sum(pca_neuro$sdev^2) * 100

    # PCA plot data
    pca_data <- data.frame(
      PC1 = pca_neuro$x[,1],
      PC2 = pca_neuro$x[,2],
      Diagnosis = metadata_df$diagnosis,
      Region = metadata_df$region
    )

    # PCA plot
    png(file.path(output_dir, "neuroinflammation_pca.png"), width = 900, height = 700, res = 100)
    ggplot(pca_data, aes(x = PC1, y = PC2, color = Diagnosis, shape = Region)) +
      geom_point(size = 3) +
      labs(
        title = "PCA of Neuroinflammation Genes",
        x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
        y = paste0("PC2 (", round(var_explained[2], 1), "%)")
      ) +
      theme_bw() +
      scale_color_manual(values = c("CONT" = "blue", "OUD" = "red")) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        legend.position = "right"
      )
    dev.off()
  }, error = function(e) {
    cat("Error generating PCA plot:", e$message, "\n")
  })

  # ----------------------
  # 5. DIFFERENTIAL EXPRESSION ANALYSIS
  # ----------------------

  # Design matrix with interaction term
  design <- model.matrix(~diagnosis * region, data = metadata_df)
  print("Available coefficients in design matrix:")
  print(colnames(design))

  # Fit model for neuroinflammation genes
  dge_neuro <- dge_filtered[neuro_ensembl, ]
  v_neuro <- voom(dge_neuro, design, plot = FALSE)
  fit <- lmFit(v_neuro, design)
  fit <- eBayes(fit)

  # Add gene symbols to all results
  add_symbols <- function(res_df) {
    res_df$Symbol <- neuro_mapping$Symbol[match(rownames(res_df), neuro_mapping$EnsemblID)]
    return(res_df)
  }

  # Extract results for main effects and interaction
  results_diag <- topTable(fit, coef = "diagnosisOUD", n = Inf)
  results_region <- topTable(fit, coef = "regionNAC", n = Inf)
  results_interaction <- topTable(fit, coef = "diagnosisOUD:regionNAC", n = Inf)

  # Apply symbol mapping
  results_diag <- add_symbols(results_diag)
  results_region <- add_symbols(results_region)
  results_interaction <- add_symbols(results_interaction)

  # Save all results
  write.csv(results_diag, file.path(output_dir, "DE_diagnosis.csv"))
  write.csv(results_region, file.path(output_dir, "DE_region.csv"))
  write.csv(results_interaction, file.path(output_dir, "DE_interaction.csv"))

  # Function to create volcano plots with proper error handling
  create_volcano <- function(results, title, filename) {
    tryCatch({
      png(file.path(output_dir, filename), width = 800, height = 700)

      # Main plot
      plot(results$logFC, -log10(results$P.Value),
           pch = 20,
           main = title,
           xlab = "log2 fold change",
           ylab = "-log10 p-value",
           col = ifelse(results$adj.P.Val < 0.05, "red", "black"))
      abline(h = -log10(0.05), col = "blue", lty = 2)
      abline(v = c(-1, 1), col = "blue", lty = 2)

      # Label significant genes if they exist
      sig_genes <- subset(results, adj.P.Val < 0.05)
      if(nrow(sig_genes) > 0 && "Symbol" %in% colnames(sig_genes) &&
         sum(!is.na(sig_genes$Symbol)) > 0) {
        # Filter out NA symbols
        sig_genes <- sig_genes[!is.na(sig_genes$Symbol),]
        if(nrow(sig_genes) > 0) {
          text(x = sig_genes$logFC,
               y = -log10(sig_genes$P.Value),
               labels = sig_genes$Symbol,
               pos = 3,
               cex = 0.8)
        }
      }

      dev.off()
    }, error = function(e) {
      cat("Error creating volcano plot", filename, ":", e$message, "\n")
    })
  }

  # Create volcano plots for all comparisons
  create_volcano(results_diag, "Volcano Plot - OUD vs Control", "volcano_diagnosis.png")
  create_volcano(results_region, "Volcano Plot - DLPFC vs NAC", "volcano_region.png")
  create_volcano(results_interaction, "Volcano Plot - Interaction (OUD:DLPFC)", "volcano_interaction.png")

  # ----------------------
  # 6. REGION-SPECIFIC EFFECTS
  # ----------------------

  # Create the contrast matrix manually to avoid syntax issues with the interaction term
  contrast_matrix <- matrix(0, ncol = 2, nrow = ncol(design))
  rownames(contrast_matrix) <- colnames(design)
  colnames(contrast_matrix) <- c("OUD_vs_CONT_in_DLPFC", "OUD_vs_CONT_in_NAC")

  # Set up the contrasts
  # For DLPFC (reference level): only the main effect is needed
  contrast_matrix["diagnosisOUD", "OUD_vs_CONT_in_DLPFC"] <- 1

  # For NAC: main effect + interaction
  contrast_matrix["diagnosisOUD", "OUD_vs_CONT_in_NAC"] <- 1
  contrast_matrix["diagnosisOUD:regionNAC", "OUD_vs_CONT_in_NAC"] <- 1

  # Print the contrast matrix to verify
  print(contrast_matrix)

  # Fit contrast
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)

  # Extract results for each region
  results_DLPFC <- topTable(fit2, coef = "OUD_vs_CONT_in_DLPFC", n = Inf)
  results_NAC <- topTable(fit2, coef = "OUD_vs_CONT_in_NAC", n = Inf)

  # Add gene symbols
  results_DLPFC <- add_symbols(results_DLPFC)
  results_NAC <- add_symbols(results_NAC)

  # Save region-specific results
  write.csv(results_DLPFC, file.path(output_dir, "DE_OUD_in_DLPFC.csv"))
  write.csv(results_NAC, file.path(output_dir, "DE_OUD_in_NAC.csv"))

  # Create volcano plots for region-specific effects
  create_volcano(results_NAC, "Volcano Plot - OUD vs Control in NAC", "volcano_OUD_in_NAC.png")
  create_volcano(results_DLPFC, "Volcano Plot - OUD vs Control in DLPFC", "volcano_OUD_in_DLPFC.png")

  # ----------------------
  # 7. PATHWAY ENRICHMENT ANALYSIS
  # ----------------------

  # Function to perform pathway enrichment on DE genes
  perform_enrichment <- function(results_df, title, output_prefix) {
    tryCatch({
      # Get significant genes (FDR < 0.1)
      sig_genes <- subset(results_df, adj.P.Val < 0.1)

      if(nrow(sig_genes) >= 10) { # Need minimum genes for meaningful enrichment
        # Get entrez IDs for significant genes
        sig_genes_entrez <- mapping_result$entrezgene_id[
          match(rownames(sig_genes), mapping_result$ensembl_gene_id)]
        sig_genes_entrez <- sig_genes_entrez[!is.na(sig_genes_entrez)]

        # Background set (all neuroinflammation genes)
        background_entrez <- mapping_result$entrezgene_id[
          match(neuro_ensembl, mapping_result$ensembl_gene_id)]
        background_entrez <- background_entrez[!is.na(background_entrez)]

        # KEGG pathway enrichment
        if(length(sig_genes_entrez) >= 10) {
          tryCatch({
            kk <- enrichKEGG(
              gene = sig_genes_entrez,
              universe = background_entrez,
              organism = "hsa",
              pvalueCutoff = 0.1
            )

            if(nrow(kk) > 0) {
              # Save results
              kk_df <- as.data.frame(kk)
              write.csv(kk_df, file.path(output_dir, paste0(output_prefix, "_KEGG_enrichment.csv")))

              # Plot
              png(file.path(output_dir, paste0(output_prefix, "_KEGG_dotplot.png")),
                  width = 900, height = 700, res = 100)
              print(dotplot(kk, showCategory = 15, title = paste0("KEGG Enrichment: ", title)))
              dev.off()
            }
          }, error = function(e) {
            cat("Error in KEGG enrichment for", title, ":", e$message, "\n")
          })
        }

        # GO Biological Process enrichment
        if(length(sig_genes_entrez) >= 10) {
          tryCatch({
            ego <- enrichGO(
              gene = sig_genes_entrez,
              universe = background_entrez,
              OrgDb = "org.Hs.eg.db",
              ont = "BP",
              pAdjustMethod = "BH",
              pvalueCutoff = 0.1,
              readable = TRUE
            )

            if(nrow(ego) > 0) {
              # Save results
              ego_df <- as.data.frame(ego)
              write.csv(ego_df, file.path(output_dir, paste0(output_prefix, "_GO_enrichment.csv")))

              # Plot
              png(file.path(output_dir, paste0(output_prefix, "_GO_dotplot.png")),
                  width = 900, height = 700, res = 100)
              print(dotplot(ego, showCategory = 15, title = paste0("GO Enrichment: ", title)))
              dev.off()
            }
          }, error = function(e) {
            cat("Error in GO enrichment for", title, ":", e$message, "\n")
          })
        }
      }
    }, error = function(e) {
      cat("Error in enrichment analysis for", title, ":", e$message, "\n")
    })
  }

  # Check if clusterProfiler and org.Hs.eg.db are available
  if("package:clusterProfiler" %in% search() || requireNamespace("clusterProfiler", quietly = TRUE)) {
    # Load org.Hs.eg.db if needed
    if(!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      cat("Warning: org.Hs.eg.db package not available. Skipping GO enrichment.\n")
    }

    # Perform enrichment for each comparison
    perform_enrichment(results_diag, "OUD vs Control", "diag")
    perform_enrichment(results_region, "NAC vs DLPFC", "region")
    perform_enrichment(results_interaction, "OUD:Region Interaction", "interaction")
    perform_enrichment(results_NAC, "OUD vs Control in NAC", "NAC")
    perform_enrichment(results_DLPFC, "OUD vs Control in DLPFC", "DLPFC")
  } else {
    cat("Warning: clusterProfiler package not available. Skipping pathway enrichment analysis.\n")
  }

  # ----------------------
  # 8. SUMMARY VISUALIZATION OF REGION-SPECIFIC EFFECTS
  # ----------------------

  # Identify top genes with region-specific effects
  sig_genes_any <- unique(c(
    rownames(results_NAC[results_NAC$adj.P.Val < 0.1, ]),
    rownames(results_DLPFC[results_DLPFC$adj.P.Val < 0.1, ]),
    rownames(results_interaction[results_interaction$adj.P.Val < 0.1, ])
  ))

  if(length(sig_genes_any) > 0) {
    # Extract expression for these genes
    sig_data <- neuro_expr[sig_genes_any, ]

    # Calculate average expression by group
    tryCatch({
      avg_expr <- sapply(split(1:ncol(sig_data), paste0(metadata_df$diagnosis, "_", metadata_df$region)),
                        function(idx) rowMeans(sig_data[, idx, drop = FALSE]))

      # Convert to data frame for plotting
      expr_long <- reshape2::melt(avg_expr)
      colnames(expr_long) <- c("GeneID", "Group", "Expression")

      # Add gene symbols
      expr_long$Symbol <- neuro_mapping$Symbol[match(expr_long$GeneID, neuro_mapping$EnsemblID)]

      # Split group into diagnosis and region
      expr_long <- expr_long %>%
        tidyr::separate(Group, into = c("Diagnosis", "Region"), sep = "_")

      # Plot expression patterns
      png(file.path(output_dir, "region_specific_expression.png"), width = 1000, height = 800)
      ggplot(expr_long, aes(x = Region, y = Expression, color = Diagnosis, group = Diagnosis)) +
        geom_point(size = 2) +
        geom_line() +
        facet_wrap(~Symbol, scales = "free_y") +
        theme_bw() +
        scale_color_manual(values = c("CONT" = "blue", "OUD" = "red")) +
        labs(title = "Region-Specific Expression Patterns of Key Neuroinflammation Genes",
             y = "Log2 CPM Expression") +
        theme(strip.background = element_rect(fill = "lightgrey"),
              strip.text = element_text(face = "bold"),
              axis.text.x = element_text(angle = 45, hjust = 1))
      dev.off()
    }, error = function(e) {
      cat("Error generating region-specific expression plot:", e$message, "\n")
    })
  }

  # ----------------------
  # 9. FINAL REPORT
  # ----------------------

  # Generate summary report
  sink(file.path(output_dir, "summary.txt"))
  cat("=== NEUROINFLAMMATION GENE ANALYSIS SUMMARY ===\n\n")
  cat("Neuroinflammation genes found in dataset:", length(neuro_ensembl), "\n")

  if(exists("pathway_details")) {
    pathway_counts <- table(pathway_details$PathwayID)
    cat("\nPathways used:\n")
    for(pw in names(pathway_counts)) {
      cat("  ", pw, ": ", pathway_counts[pw], " genes\n", sep="")
    }
  }

  cat("\n=== DIFFERENTIAL EXPRESSION RESULTS ===\n")
  cat("Main effect of OUD (vs Control):\n")
  cat("  Significant genes (FDR < 0.05):", sum(results_diag$adj.P.Val < 0.05), "\n")
  cat("  Significant genes (FDR < 0.1):", sum(results_diag$adj.P.Val < 0.1), "\n")

  cat("\nMain effect of Region (NAC vs DLPFC):\n")
  cat("  Significant genes (FDR < 0.05):", sum(results_region$adj.P.Val < 0.05), "\n")
  cat("  Significant genes (FDR < 0.1):", sum(results_region$adj.P.Val < 0.1), "\n")

  cat("\nInteraction effect (OUD × Region):\n")
  cat("  Significant genes (FDR < 0.05):", sum(results_interaction$adj.P.Val < 0.05), "\n")
  cat("  Significant genes (FDR < 0.1):", sum(results_interaction$adj.P.Val < 0.1), "\n")

  cat("\nRegion-specific OUD effects:\n")
  cat("  NAC: Significant genes (FDR < 0.05):", sum(results_NAC$adj.P.Val < 0.05), "\n")
  cat("  DLPFC: Significant genes (FDR < 0.05):", sum(results_DLPFC$adj.P.Val < 0.05), "\n")

  cat("\n\n=== TOP GENES BY CONDITION ===\n")

  # Function to print top genes
  print_top_genes <- function(results_df, label, n = 5) {
    cat("\n", label, ":\n", sep="")
    if(nrow(results_df) > 0 && "Symbol" %in% colnames(results_df)) {
      top_genes <- results_df[order(results_df$P.Value), c("Symbol", "logFC", "P.Value", "adj.P.Val")]
      print(head(top_genes, n))
    } else {
      cat("No genes found or symbol information missing.\n")
    }
  }

  print_top_genes(results_diag, "Top genes by OUD effect")
  print_top_genes(results_region, "Top genes by region effect")
  print_top_genes(results_interaction, "Top genes with interaction effect")
  print_top_genes(results_NAC, "Top genes with OUD effect in NAC")
  print_top_genes(results_DLPFC, "Top genes with OUD effect in DLPFC")

  cat("\n\nAnalysis complete. Results saved to:", output_dir, "\n")
  sink()
} else {
  # If no genes were found
  sink(file.path(output_dir, "summary.txt"))
  cat("=== NEUROINFLAMMATION GENE ANALYSIS SUMMARY ===\n\n")
  cat("No neuroinflammation-related genes were found in the dataset.\n")
  cat("This could be due to:\n")
  cat("1. Mapping issues between Ensembl IDs and gene symbols\n")
  cat("2. Neuroinflammation genes not present in the filtered dataset\n")
  cat("3. Failed connection to KEGG database\n\n")

  cat("First few gene IDs in the dataset:", head(ensembl_ids), "\n\n")
  cat("Check the ensembl_to_symbol_mapping.csv file to see if the mapping was successful\n")
  sink()
}

cat("✓ Neuroinflammation analysis complete. Output saved to:", output_dir, "\n")