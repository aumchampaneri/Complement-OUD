# Enhanced GSEA Visualization Script
library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(enrichplot)     # For advanced GSEA plots
library(clusterProfiler) # For handling GSEA data
library(ComplexHeatmap)  # For heatmap
library(treemapify)     # For treemap
library(AnnotationDbi)  # For gene ID mapping
library(org.Hs.eg.db)   # For human gene annotations
library(ggridges)       # For ridge plots

# Define paths
base_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD"
gsea_output_dir <- file.path(base_dir, "GSE225158/GSEA outputs")
plot_output_dir <- file.path(gsea_output_dir, "plots")
dir.create(plot_output_dir, showWarnings = FALSE, recursive = TRUE)

# Read and preprocess the GSEA data
gsea_file_path <- file.path(gsea_output_dir, "gsea_inflammatory_complement_pathways.csv")
df <- read.csv(gsea_file_path)

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

# Filter for significant pathways
df_filtered <- df %>%
  filter(!is.na(padj), padj < 0.05) %>%  # Exclude NA p-values
  arrange(desc(abs(NES))) %>%
  slice_head(n = 25)  # Limit to top 25 by absolute NES

# 1. IMPROVED DOT PLOT
message("Creating dot plot...")
df_filtered <- df_filtered %>% arrange(NES)
p <- ggplot(df_filtered, aes(x = NES, y = reorder(Description, NES), color = NES, size = -log10(padj))) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "Pathway",
    color = "NES",
    size = "-log10(padj)",
    title = "GSEA Pathway Enrichment"
  ) +
  theme(axis.text.y = element_text(size = 8))
ggsave(file.path(plot_output_dir, "gsea_dotplot.png"), p + theme(plot.margin = margin(5, 5, 5, 20)),
       width = 20, height = 10)

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
    pdf(file.path(plot_output_dir, "enrichment_plots.pdf"), width=12, height=10)
    print(gseaplot2(gsea_data, geneSetID=top_pathway_ids, pvalue_table=TRUE))
    dev.off()
  }
}

# 3. RIDGE PLOT
message("Creating ridge plot...")
p <- ggplot(df_filtered, aes(x = NES, y = reorder(Description, NES), fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01, quantile_lines = TRUE, quantiles = 2) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.y = element_text(size = 8)) +
  labs(title = "Distribution of Enrichment Scores",
       x = "Normalized Enrichment Score (NES)",
       y = "")
ggsave(file.path(plot_output_dir, "ridge_plot.png"), p, width = 10, height = 8)

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

    expr_mat <- expr_df %>%
      select(gene_symbol, log2FoldChange) %>%
      column_to_rownames("gene_symbol") %>%
      as.matrix()

    png(file.path(plot_output_dir, "leading_edge_heatmap.png"), width=800, height=1000)
    print(Heatmap(expr_mat, name="log2FC",
                  cluster_rows=TRUE,
                  show_row_names=TRUE,
                  row_names_gp=gpar(fontsize=8)))
    dev.off()
  }
}

# 5. CATEGORY TREEMAP
message("Creating treemap...")
df_categorized <- df_filtered %>%
  mutate(Category = case_when(
    grepl("INFLAMMATION|CYTOKINE|IMMUNE|IL|TNF|INTERFERON", Description, ignore.case=TRUE) ~ "Inflammation",
    grepl("COMPLEMENT", Description, ignore.case=TRUE) ~ "Complement",
    grepl("APOPTOSIS|DEATH", Description, ignore.case=TRUE) ~ "Cell Death",
    grepl("SIGNAL|PATHWAY", Description, ignore.case=TRUE) ~ "Signaling",
    TRUE ~ "Other"
  ))

p <- ggplot(df_categorized, aes(area = abs(NES), fill = NES,
                              label = Description, subgroup = Category)) +
  geom_treemap() +
  geom_treemap_subgroup_border() +
  geom_treemap_subgroup_text(place = "centre", grow = TRUE, alpha = 0.5, colour = "black") +
  geom_treemap_text(aes(label = Description), place = "center", grow = FALSE,
                  size = 10, colour = "black") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme(legend.position = "right") +
  labs(title = "GSEA Pathway Enrichment by Category", fill = "NES")
ggsave(file.path(plot_output_dir, "pathway_treemap.png"), p, width = 15, height = 12)

# 6. PATHWAY NETWORK PLOT
message("Creating network plot...")
if(file.exists(gsea_obj_path)) {
  gsea_data <- readRDS(gsea_obj_path)

  # Try to create network plot with filtered significant pathways
  tryCatch({
    sig_pathway_ids <- match(df_filtered$ID, gsea_data@result$ID)
    sig_pathway_ids <- sig_pathway_ids[!is.na(sig_pathway_ids)]

    if(length(sig_pathway_ids) > 0) {
      png(file.path(plot_output_dir, "pathway_network.png"), width=1200, height=1000)
      print(emapplot(gsea_data, showCategory=sig_pathway_ids[1:min(15, length(sig_pathway_ids))], color="NES"))
      dev.off()
    }
  }, error = function(e) {
    # Fallback to default behavior if the above fails
    message("Using default network plot...")
    png(file.path(plot_output_dir, "pathway_network.png"), width=1200, height=1000)
    print(emapplot(gsea_data, showCategory=15, color="NES"))
    dev.off()
  })
}

message("All visualizations completed! Files saved to:", plot_output_dir)