# Enhanced GSEA Visualization Script
library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(enrichplot)
library(clusterProfiler)
library(ComplexHeatmap)
library(treemapify)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggridges)
library(tibble)

# Define paths
base_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD"
gsea_output_dir <- file.path(base_dir, "GSE225158/GSEA outputs")
plot_output_dir <- file.path(gsea_output_dir, "plots")
dir.create(plot_output_dir, showWarnings = FALSE, recursive = TRUE)

# Read and preprocess the GSEA data
gsea_file_path <- file.path(gsea_output_dir, "gsea_results_M_OUD_vs_M_None.csv")
df <- read.csv(gsea_file_path)

# Print column names to troubleshoot
cat("Available columns in GSEA results:", paste(colnames(df), collapse=", "), "\n")

# Load original DESeq2 results for heatmap
deseq_results_path <- file.path(base_dir, "GSE225158/DESeq2 outputs/deseq2_results_M_OUD_vs_M_None.csv")
res <- read.csv(deseq_results_path)

# Ensure we have entrez IDs in the results
if (!"entrez" %in% colnames(res)) {
  res$entrez <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = res$gene,
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )
}

# Determine which adjusted p-value column to use - MORE ROBUST
p_col <- if("p.adjust" %in% colnames(df)) {
  "p.adjust"
} else if("adj.P.Val" %in% colnames(df)) {
  "adj.P.Val"
} else if("padj" %in% colnames(df)) {
  "padj"
} else {
  message("Could not find adjusted p-value column, using 'pvalue'")
  "pvalue"
}

# Filter for significant pathways - ADDED NES CHECK
df_filtered <- df %>%
  filter(!is.na(!!sym(p_col)), !is.na(NES), !!sym(p_col) < 0.05) %>%
  arrange(desc(abs(NES))) %>%
  slice_head(n = 25)  # Limit to top 25 by absolute NES

# 1. IMPROVED DOT PLOT
message("Creating dot plot...")
df_filtered <- df_filtered %>% arrange(NES)
p <- ggplot(df_filtered, aes(x = NES, y = reorder(Description, NES),
                             color = NES, size = -log10(!!sym(p_col)))) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_bw(base_size = 12) + # PUBLICATION-READY THEME
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "Pathway",
    color = "NES",
    size = paste0("-log10(", p_col, ")"),
    title = "GSEA Pathway Enrichment"
  ) +
  theme(axis.text.y = element_text(size = 9),
        plot.title = element_text(face = "bold", size = 14))
ggsave(file.path(plot_output_dir, "gsea_dotplot.png"),
       p + theme(plot.margin = margin(5, 5, 5, 20)),
       width = 20, height = 10, dpi = 300)

# 2. ENRICHMENT PLOTS
message("Creating enrichment plots...")
# Load the saved GSEA object
gsea_obj_path <- file.path(gsea_output_dir, "gsea_obj.rds")
if(file.exists(gsea_obj_path)) {
  gsea_data <- readRDS(gsea_obj_path)

  # Get top pathways from our filtered dataset to ensure proper ID matching
  top_pathway_ids <- match(df_filtered$ID[1:min(5, nrow(df_filtered))], gsea_data@result$ID)
  top_pathway_ids <- top_pathway_ids[!is.na(top_pathway_ids)]

  if(length(top_pathway_ids) > 0) {
    pdf(file.path(plot_output_dir, "enrichment_plots.pdf"), width=15, height=10)
    print(gseaplot2(gsea_data, geneSetID=top_pathway_ids, pvalue_table=TRUE))
    dev.off()
  }
}

# 3. ENHANCED BAR PLOT with value labels and better coloring
message("Creating bar plot...")
p <- ggplot(df_filtered, aes(x = reorder(Description, NES), y = NES,
                             fill = NES > 0)) + # CATEGORICAL COLORING
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", NES),
                hjust = ifelse(NES > 0, -0.1, 1.1)),
            size = 3) + # ADD VALUE LABELS
  scale_fill_manual(values = c("steelblue", "firebrick"),
                    labels = c("Down", "Up"),
                    name = "Regulation") + # CLEAR LABELS
  coord_flip() +
  theme_bw() + # PUBLICATION-READY
  theme(legend.position = "right",
        axis.text.y = element_text(size = 9),
        plot.title = element_text(face = "bold", size = 14)) +
  labs(title = "GSEA Pathway Enrichment Scores",
       x = "",
       y = "Normalized Enrichment Score (NES)")
ggsave(file.path(plot_output_dir, "pathway_barplot.png"), p,
       width = 16, height = 8, dpi = 300)

# 4. HEATMAP OF LEADING EDGE GENES
message("Creating heatmap...")
if("core_enrichment" %in% colnames(df_filtered)) {
  top_pathways <- df_filtered %>% head(10)
  all_genes <- c()

  for(i in 1:nrow(top_pathways)) {
    if(!is.na(top_pathways$core_enrichment[i])) {
      genes <- unlist(strsplit(as.character(top_pathways$core_enrichment[i]), "/"))
      all_genes <- c(all_genes, genes)
    }
  }

  unique_genes <- unique(all_genes)

  # Only proceed if we have genes
  if(length(unique_genes) > 0) {
    expr_df <- res %>%
      filter(entrez %in% unique_genes) %>%
      arrange(desc(abs(log2FoldChange)))

    gene_symbols <- mapIds(
      org.Hs.eg.db,
      keys = expr_df$entrez,
      column = "SYMBOL",
      keytype = "ENTREZID",
      multiVals = "first"
    )

    expr_df$gene_symbol <- gene_symbols
    if(nrow(expr_df) > 50) expr_df <- expr_df %>% head(50)

    # Use dplyr::select explicitly to avoid namespace conflicts
    expr_mat <- expr_df %>%
      dplyr::select(gene_symbol, log2FoldChange) %>%
      tibble::column_to_rownames("gene_symbol") %>%
      as.matrix()

    png(file.path(plot_output_dir, "leading_edge_heatmap.png"), width=800, height=1000, res=120)
    print(Heatmap(expr_mat, name="log2FC",
                  cluster_rows=TRUE,
                  show_row_names=TRUE,
                  row_names_gp=gpar(fontsize=8)))
    dev.off()
  }
}

# 5. NETWORK PLOT with CONSISTENT PATHWAY SELECTION
message("Creating network plot...")
if(file.exists(gsea_obj_path)) {
  gsea_data <- readRDS(gsea_obj_path)

  # Create properly formatted gene list
  gene_list_for_viz <- res$log2FoldChange
  names(gene_list_for_viz) <- as.character(res$entrez)
  gene_list_for_viz <- gene_list_for_viz[!is.na(names(gene_list_for_viz))]
  gene_list_for_viz <- gene_list_for_viz[!duplicated(names(gene_list_for_viz))]

  tryCatch({
    # Try dotplot first (most reliable)
    png(file.path(plot_output_dir, "pathway_dotplot.png"), width=1000, height=800, res=120)
    print(dotplot(gsea_data, showCategory=15, title="Top Enriched Pathways"))
    dev.off()

    # Use SAME TOP PATHWAYS as in other plots for consistency
    png(file.path(plot_output_dir, "pathway_gene_network.png"), width=1200, height=1000, res=120)
    top_pathways_ids <- df_filtered$ID[1:min(5, nrow(df_filtered))]
    print(cnetplot(gsea_data,
                  categorySize="pvalue",
                  showCategory=top_pathways_ids,
                  foldChange=gene_list_for_viz,
                  colorEdge=FALSE))
    dev.off()
  }, error = function(e) {
    message("Network plot failed: ", e$message)
  })
}

# 6. Save session info for reproducibility
writeLines(capture.output(sessionInfo()),
           file.path(plot_output_dir, "session_info.txt"))

message("All visualizations completed! Files saved to:", plot_output_dir)