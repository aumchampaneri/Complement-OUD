# Enhanced GSEA Visualization Script with improved error handling and organization

# ==== Library Imports ====
library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(enrichplot)
library(clusterProfiler)
library(ComplexHeatmap)
library(circlize)  # Added for colorRamp2 function
library(treemapify)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggridges)
library(tibble)
library(grid)

# Resolve namespace conflicts
select <- dplyr::select  # Explicitly use dplyr's select

# ==== Path Configuration ====
base_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD"
gsea_output_dir <- file.path(base_dir, "GSE225158/GSEA outputs")
plot_output_dir <- file.path(gsea_output_dir, "plots")
dir.create(plot_output_dir, showWarnings = FALSE, recursive = TRUE)

# ==== Data Loading ====
# Read and preprocess the GSEA data
# gsea_file_path <- file.path(gsea_output_dir, "gsea_results_F_OUD_vs_F_None.csv")
gsea_file_path <- file.path(gsea_output_dir, "gsea_F_inflammatory_complement_pathways.csv")
if(!file.exists(gsea_file_path)) {
  stop(paste("GSEA results file not found:", gsea_file_path))
}
df <- read.csv(gsea_file_path)

# Print column names and row count to troubleshoot
cat("Available columns in GSEA results:", paste(colnames(df), collapse=", "), "\n")
cat("Number of GSEA results:", nrow(df), "\n")

# Load original DESeq2 results for heatmap
deseq_results_path <- file.path(base_dir, "GSE225158/DESeq2 outputs/deseq2_results_F_OUD_vs_F_None.csv")
if(!file.exists(deseq_results_path)) {
  cat("WARNING: DESeq2 results file not found:", deseq_results_path, "\n")
  res <- NULL
} else {
  res <- read.csv(deseq_results_path)
  cat("DESeq2 results loaded with", nrow(res), "rows\n")
}

# ==== Data Preprocessing ====
# Ensure we have entrez IDs in the results if DESeq2 data exists
if (!is.null(res)) {
  if (!"entrez" %in% colnames(res)) {
    res$entrez <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = res$gene,
      column = "ENTREZID",
      keytype = "SYMBOL",
      multiVals = "first"
    )
    cat("Mapped", sum(!is.na(res$entrez)), "out of", nrow(res), "genes to Entrez IDs\n")
  }
}

# Determine which adjusted p-value column to use
p_col <- if("p.adjust" %in% colnames(df)) {
  "p.adjust"
} else if("adj.P.Val" %in% colnames(df)) {
  "adj.P.Val"
} else if("padj" %in% colnames(df)) {
  "padj"
} else {
  cat("WARNING: Could not find adjusted p-value column, using 'pvalue'\n")
  "pvalue"
}

# Filter for significant pathways
df_filtered <- df %>%
  filter(!is.na(!!sym(p_col)), !is.na(NES), !!sym(p_col) < 0.05) %>%
  arrange(desc(abs(NES)))

# Check if we have enough pathways
if(nrow(df_filtered) == 0) {
  cat("WARNING: No significant pathways found! Check significance threshold.\n")
} else {
  cat("Found", nrow(df_filtered), "significant pathways\n")
  # Limit to top 25 by absolute NES
  df_filtered <- df_filtered %>% slice_head(n = 25)
  cat("Using top 25 pathways for visualization\n")
}

# ==== GSEA Object Loading ====
gsea_obj_paths <- c(
  file.path(gsea_output_dir, "gsea_obj_F_OUD_vs_F_None.rds")
)

gsea_data <- NULL
gsea_obj_path <- NULL

for (path in gsea_obj_paths) {
  if(file.exists(path)) {
    cat("Found GSEA object at:", path, "\n")
    tryCatch({
      gsea_data <- readRDS(path)
      gsea_obj_path <- path
      cat("Successfully loaded GSEA object\n")
      break
    }, error = function(e) {
      cat("ERROR loading GSEA object from", path, ":", e$message, "\n")
    })
  }
}

if(is.null(gsea_data)) {
  cat("WARNING: Could not find or load GSEA object. Some plots will be skipped.\n")
}

# ==== Plot Generation ====
# 1. CUSTOM DOT PLOT
if(nrow(df_filtered) > 0) {
  tryCatch({
    cat("Creating custom dot plot...\n")
    df_filtered <- df_filtered %>% arrange(NES)
    p <- ggplot(df_filtered, aes(x = NES, y = reorder(Description, NES),
                               color = NES, size = -log10(!!sym(p_col)))) +
      geom_point() +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      theme_bw(base_size = 12) +
      labs(
        x = "Normalized Enrichment Score (NES)",
        y = "Pathway",
        color = "NES",
        size = paste0("-log10(", p_col, ")"),
        title = "GSEA Pathway Enrichment"
      ) +
      theme(axis.text.y = element_text(size = 9),
            plot.title = element_text(face = "bold", size = 14))
    dot_plot_file <- file.path(plot_output_dir, "gsea_dotplot.png")
    ggsave(dot_plot_file,
         p + theme(plot.margin = margin(5, 5, 5, 20)),
         width = 20, height = 10, dpi = 300)
    cat("Custom dot plot saved:", dot_plot_file, "\n")
  }, error = function(e) {
    cat("ERROR creating custom dot plot:", e$message, "\n")
  })
} else {
  cat("Skipping custom dot plot - no significant pathways\n")
}

# 2. ENRICHPLOT DOTPLOT
if(!is.null(gsea_data) && nrow(df_filtered) > 0) {
  tryCatch({
    cat("Creating enrichplot dotplot...\n")
    top_pathways <- head(df_filtered, 20)
    matched_idx <- match(top_pathways$ID, gsea_data@result$ID)
    matched_idx <- matched_idx[!is.na(matched_idx)]

    if(length(matched_idx) > 0) {
      gsea_subset <- gsea_data
      gsea_subset@result <- gsea_data@result[matched_idx, ]

      enrich_dotplot_file <- file.path(plot_output_dir, "gsea_enrichment_dotplot.png")
      png(enrich_dotplot_file, width = 1000, height = 800, res = 100)
      p <- dotplot(gsea_subset, showCategory=20,
                   title="Enriched Pathways Dotplot")
      print(p)
      dev.off()
      cat("Enrichment dotplot saved:", enrich_dotplot_file, "\n")
    } else {
      cat("WARNING: No matching pathway IDs found for enrichment dotplot\n")
    }
  }, error = function(e) {
    cat("ERROR creating enrichment dotplot:", e$message, "\n")
  })
} else {
  if(is.null(gsea_data)) {
    cat("Skipping enrichment dotplot - GSEA object not loaded\n")
  } else {
    cat("Skipping enrichment dotplot - no significant pathways\n")
  }
}

# 3. ENRICHMENT PLOTS - OVERLAID VERSION
cat("Creating enrichment plots...\n")
if(!is.null(gsea_data) && nrow(df_filtered) > 0) {
  tryCatch({
    # Get top pathways from our filtered dataset
    top_pathway_ids <- match(df_filtered$ID[1:min(5, nrow(df_filtered))], gsea_data@result$ID)
    top_pathway_ids <- top_pathway_ids[!is.na(top_pathway_ids)]

    if(length(top_pathway_ids) > 0) {
      # Create one overlaid plot with all top pathways
      cat("Creating overlaid enrichment plot for top pathways\n")
      overlaid_plot_file <- file.path(plot_output_dir, "gsea_overlaid_enrichment.png")

      png(overlaid_plot_file, width = 1500, height = 1000, res = 100)
      p <- gseaplot2(gsea_data,
                   geneSetID = top_pathway_ids,
                   title = "Top Enriched Pathways",
                   subplots = 1:3,  # Include all three subplots
                   pvalue_table = TRUE)
      print(p)
      dev.off()
      cat("Saved overlaid plot:", overlaid_plot_file, "\n")

      # Also create individual plots for reference
      for (i in seq_along(top_pathway_ids)) {
        pathway_id <- top_pathway_ids[i]
        pathway_name <- gsea_data@result$Description[pathway_id]
        pathway_name_clean <- gsub("[^[:alnum:]]", "_", pathway_name)

        cat("Creating individual enrichment plot for", pathway_name, "\n")
        gseaplot_file <- file.path(plot_output_dir, paste0("gsea_", i, "_", pathway_name_clean, ".png"))

        png(gseaplot_file, width = 1500, height = 1000, res = 100)
        p <- gseaplot2(gsea_data, geneSetID = pathway_id, title = pathway_name)
        print(p)
        dev.off()
      }
    } else {
      cat("WARNING: No matching pathway IDs found for enrichment plots\n")
    }
  }, error = function(e) {
    cat("ERROR creating enrichment plots:", e$message, "\n")
  })
} else {
  if(is.null(gsea_data)) {
    cat("Skipping enrichment plots - GSEA object not loaded\n")
  } else {
    cat("Skipping enrichment plots - no significant pathways\n")
  }
}

# 4. ENHANCED BAR PLOT
if(nrow(df_filtered) > 0) {
  tryCatch({
    cat("Creating bar plot...\n")
    p <- ggplot(df_filtered, aes(x = reorder(Description, NES), y = NES,
                               fill = NES > 0)) +
      geom_col(width = 0.7) +
      geom_text(aes(label = sprintf("%.2f", NES),
                  hjust = ifelse(NES > 0, -0.1, 1.1)),
              size = 3) +
      scale_fill_manual(values = c("steelblue", "firebrick"),
                      labels = c("Down", "Up"),
                      name = "Regulation") +
      coord_flip() +
      theme_bw() +
      theme(legend.position = "right",
            axis.text.y = element_text(size = 9)) +
      labs(title = "GSEA Pathway Enrichment Scores",
           x = "Pathway",
           y = "Normalized Enrichment Score")
    bar_plot_file <- file.path(plot_output_dir, "pathway_barplot.png")
    ggsave(bar_plot_file, p, width = 14, height = 10, dpi = 300)
    cat("Bar plot saved:", bar_plot_file, "\n")
  }, error = function(e) {
    cat("ERROR creating bar plot:", e$message, "\n")
  })
} else {
  cat("Skipping bar plot - no significant pathways\n")
}

# 5. HEATMAP AND NETWORK OF LEADING EDGE GENES - Fixed version
cat("Creating leading edge visualizations...\n")
if(nrow(df_filtered) > 0 && "core_enrichment" %in% colnames(df_filtered) && !is.null(res)) {
  tryCatch({
    # Get top pathways (limit to 8 for visibility)
    top_pathways <- df_filtered %>%
      slice_head(n = 8) %>%
      select(ID, Description, NES, core_enrichment)

    # Parse core enrichment genes
    genes_list <- lapply(top_pathways$core_enrichment, function(x) {
      if(is.na(x) || x == "") return(character(0))
      strsplit(x, "/")[[1]]
    })
    names(genes_list) <- top_pathways$Description

    # Check if we have genes
    if(all(lengths(genes_list) == 0)) {
      cat("WARNING: No genes found in core_enrichment strings\n")
    } else {
      # 1. IMPROVED HEATMAP WITH ANNOTATION
      all_genes <- unique(unlist(genes_list))
      cat("Found", length(all_genes), "unique core enrichment genes\n")

      # Get gene symbols for better readability
      gene_symbols <- bitr(all_genes,
                         fromType = "ENTREZID",
                         toType = "SYMBOL",
                         OrgDb = org.Hs.eg.db)

      # Create expression matrix for genes found in DESeq results
      expr_data <- res %>%
        filter(entrez %in% all_genes) %>%
        select(gene, entrez, log2FoldChange, pvalue) %>%
        arrange(desc(abs(log2FoldChange)))

      # Create pathway annotation - which gene belongs to which pathway
      pathway_membership <- matrix(0,
                                 nrow = nrow(expr_data),
                                 ncol = length(genes_list),
                                 dimnames = list(expr_data$gene, names(genes_list)))

      for(i in seq_along(genes_list)) {
        pathway_membership[expr_data$gene[expr_data$entrez %in% genes_list[[i]]], i] <- 1
      }

      # Create enhanced heatmap
      heatmap_file <- file.path(plot_output_dir, "leading_edge_heatmap.png")
      png(heatmap_file, width = 1200, height = 1700, res = 100)

      # Create color scales
      col_fc <- circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
      col_pval <- circlize::colorRamp2(c(0, 0.05, 0.1), c("red", "white", "lightgray"))
      col_membership <- c("0" = "white", "1" = "darkgreen")

      # Create annotations
      ha_column <- HeatmapAnnotation(
        NES = top_pathways$NES,
        col = list(NES = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
      )

      ha_row <- rowAnnotation(
        p_value = expr_data$pvalue,
        col = list(p_value = col_pval)
      )

      # Draw the heatmaps
      ht_list <- Heatmap(as.matrix(expr_data$log2FoldChange),
                       name = "log2FC",
                       col = col_fc,
                       row_names_gp = gpar(fontsize = 9),
                       show_row_names = length(all_genes) < 50, # Show names only if not too many
                       column_title = "Log2 Fold Change") +
                Heatmap(pathway_membership,
                       name = "In Pathway",
                       col = col_membership,
                       show_row_names = FALSE,
                       top_annotation = ha_column,
                       right_annotation = ha_row,
                       column_title = "Pathway Membership",
                       column_names_gp = gpar(fontsize = 9))

      draw(ht_list, padding = unit(c(2, 10, 2, 2), "mm"))
      dev.off()
      cat("Leading edge heatmap saved:", heatmap_file, "\n")

      # 2. GENE-CONCEPT NETWORK - Fixed to avoid gene_list error
      if(!is.null(gsea_data)) {
        concept_net_file <- file.path(plot_output_dir, "gene_concept_network.png")
        png(concept_net_file, width = 1400, height = 1000, res = 100)

        # Match pathway IDs to those in the GSEA object
        matched_idx <- match(top_pathways$ID, gsea_data@result$ID)
        matched_idx <- matched_idx[!is.na(matched_idx)]

        if(length(matched_idx) > 0) {
          # Create fold change vector for coloring genes
          gene_fc <- res$log2FoldChange
          names(gene_fc) <- res$entrez
          gene_fc <- gene_fc[!is.na(names(gene_fc))]

          gsea_subset <- gsea_data
          gsea_subset@result <- gsea_data@result[matched_idx, ]

          tryCatch({
            p <- cnetplot(gsea_subset,
                       showCategory = min(5, length(matched_idx)),
                       categorySize = "pvalue",
                       foldChange = gene_fc,
                       colorEdge = TRUE)
            print(p)
          }, error = function(e) {
            cat("WARNING: Error in cnetplot:", e$message,
                "\nTrying simpler network plot...\n")
            # Fallback to simpler network
            p <- cnetplot(gsea_subset,
                       showCategory = min(5, length(matched_idx)),
                       colorEdge = TRUE)
            print(p)
          })
        }
        dev.off()
        cat("Gene-concept network saved:", concept_net_file, "\n")
      }

      # 3. RIDGEPLOT OF GENE RANKS
      if(!is.null(gsea_data)) {
        ridge_file <- file.path(plot_output_dir, "leading_edge_ridgeplot.png")
        png(ridge_file, width = 1000, height = 800, res = 100)

        # Match pathway IDs
        matched_idx <- match(top_pathways$ID, gsea_data@result$ID)
        matched_idx <- matched_idx[!is.na(matched_idx)]

        if(length(matched_idx) > 0) {
          gsea_subset <- gsea_data
          gsea_subset@result <- gsea_data@result[matched_idx, ]

          p <- ridgeplot(gsea_subset, showCategory = min(8, length(matched_idx)))
          print(p)
        }
        dev.off()
        cat("Ridge plot saved:", ridge_file, "\n")
      }
    }
  }, error = function(e) {
    cat("ERROR creating leading edge visualizations:", e$message, "\n")
  })
} else {
  cat("Skipping leading edge visualizations - missing required data\n")
}

# 6. NETWORK PLOT
cat("Creating network plot...\n")
if(!is.null(gsea_data) && nrow(df_filtered) > 0) {
  tryCatch({
    # Create emapplot for network visualization
    network_plot_file <- file.path(plot_output_dir, "pathway_network.png")

    # Use only top 20 pathways to avoid overcrowded network
    top_ids <- head(df_filtered$ID, 20)
    matched_idx <- match(top_ids, gsea_data@result$ID)
    matched_idx <- matched_idx[!is.na(matched_idx)]

    if(length(matched_idx) > 0) {
      gsea_subset <- gsea_data
      gsea_subset@result <- gsea_data@result[matched_idx, ]

      # Create and save the network plot
      png(network_plot_file, width = 1500, height = 1000, res = 100)
      p <- emapplot(pairwise_termsim(gsea_subset),
                  showCategory = min(20, length(matched_idx)))
      print(p)
      dev.off()
      cat("Network plot saved:", network_plot_file, "\n")
    } else {
      cat("WARNING: No matching pathway IDs found for network plot\n")
    }
  }, error = function(e) {
    cat("ERROR creating network plot:", e$message, "\n")
  })
} else {
  if(is.null(gsea_data)) {
    cat("Skipping network plot - GSEA object not loaded\n")
  } else if(nrow(df_filtered) == 0) {
    cat("Skipping network plot - no significant pathways\n")
  }
}

# 7. Save session info for reproducibility
info_file <- file.path(plot_output_dir, "session_info.txt")
writeLines(capture.output(sessionInfo()), info_file)
cat("Session info saved:", info_file, "\n")

cat("All visualizations completed! Files saved to:", plot_output_dir, "\n")