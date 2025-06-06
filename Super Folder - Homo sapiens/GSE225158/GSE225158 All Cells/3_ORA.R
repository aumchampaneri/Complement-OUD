# Consolidated ORA and Visualization Script with Error Handling and Organized Output

library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
library(enrichplot)
library(ggplot2)
library(readr)
library(dplyr)

# Create main output directories
dir.create("GSE225158/ORA outputs", recursive = TRUE, showWarnings = FALSE)
dir.create("GSE225158/ORA outputs/visualizations", recursive = TRUE, showWarnings = FALSE)

# Create a directory for KEGG data
kegg_dir <- "GSE225158/ORA outputs/visualizations/kegg_data"
dir.create(kegg_dir, recursive = TRUE, showWarnings = FALSE)

# Function to set up high-quality PNG output
create_high_res_png <- function(filename, width_in = 8, height_in = 6, dpi = 300) {
  png(filename, width = width_in, height = height_in, units = "in", res = dpi)
}

# Define a safe pathview function with error handling and improved quality
pathview_safe <- function(gene_fc, pathway_id, pathway_name, contrast, kegg_dir, out_dir) {
  tryCatch({
    cat(sprintf("Processing pathway %s (%s)...\n", pathway_id, pathway_name))
    pathview(
      gene.data = gene_fc,
      pathway.id = pathway_id,
      species = "hsa",
      out.suffix = paste0(contrast, "_", pathway_name),
      kegg.dir = kegg_dir,
      out.dir = out_dir,
      node.sum = "mean",
      limit = list(gene = 2, cpd = 1),
      low = list(gene = "blue", cpd = "blue"),
      mid = list(gene = "white", cpd = "white"),
      high = list(gene = "red", cpd = "red"),
      same.layer = TRUE,
      new.signature = FALSE,
      pdf.size = c(10, 8)  # Creates larger PDF output
    )
    cat("Success!\n")
  }, error = function(e) {
    cat(sprintf("Error processing pathway %s: %s\n", pathway_id, e$message))
  })
}

# Function to create MA plot with improved resolution
create_ma_plot <- function(res, contrast, out_dir) {
  tryCatch({
    cat("Creating MA plot...\n")
    res_ma <- res[res$baseMean > 0, ]

    # High-res PNG
    plot_file <- file.path(out_dir, paste0('MA_plot_', contrast, '.png'))
    create_high_res_png(plot_file, width_in = 8, height_in = 6, dpi = 300)
    with(res_ma, plot(baseMean, log2FoldChange, pch=20, cex=0.5, log='x',
         main=paste('MA Plot:', contrast), xlab='Mean Expression', ylab='log2 Fold Change',
         cex.axis = 1.2, cex.lab = 1.2, cex.main = 1.4))
    abline(h=0, col='red')
    dev.off()

    # PDF version
    pdf_file <- file.path(out_dir, paste0('MA_plot_', contrast, '.pdf'))
    pdf(pdf_file, width = 8, height = 6)
    with(res_ma, plot(baseMean, log2FoldChange, pch=20, cex=0.5, log='x',
         main=paste('MA Plot:', contrast), xlab='Mean Expression', ylab='log2 Fold Change',
         cex.axis = 1.2, cex.lab = 1.2, cex.main = 1.4))
    abline(h=0, col='red')
    dev.off()

    cat("Success!\n")
  }, error = function(e) {
    cat(sprintf("Error creating MA plot: %s\n", e$message))
  })
}

# Function to create Volcano plot with improved resolution
create_volcano_plot <- function(res, contrast, out_dir) {
  tryCatch({
    cat("Creating Volcano plot...\n")
    res$padj[is.na(res$padj)] <- 1
    res$significant <- res$padj < 0.05

    # Base plot with improved aesthetics
    p <- ggplot(res, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
      geom_point(alpha=0.6, size=1.2) +
      scale_color_manual(values=c('grey','red')) +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = "right",
        plot.title = element_text(size=16, face="bold"),
        axis.title = element_text(size=14, face="bold"),
        axis.text = element_text(size=12)
      ) +
      labs(title=paste('Volcano Plot:', contrast), x='log2 Fold Change', y='-log10 Adjusted p-value')

    # High-res PNG
    plot_file <- file.path(out_dir, paste0('Volcano_plot_', contrast, '.png'))
    ggsave(plot_file, p, width=8, height=7, dpi=300)

    # PDF version
    pdf_file <- file.path(out_dir, paste0('Volcano_plot_', contrast, '.pdf'))
    ggsave(pdf_file, p, width=8, height=7)

    cat("Success!\n")
  }, error = function(e) {
    cat(sprintf("Error creating Volcano plot: %s\n", e$message))
  })
}

# Function to create enrichment barplot with improved resolution
create_barplot <- function(enrich_obj, term_type, contrast, out_dir) {
  tryCatch({
    cat(sprintf("Creating barplot for %s...\n", term_type))
    if(is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) {
      cat("No enrichment results to plot.\n")
      return()
    }

    # Create enhanced barplot with better styling
    p <- barplot(enrich_obj, showCategory=20) +
         ggtitle(paste0("Enriched ", term_type, " Terms: ", contrast)) +
         theme_minimal(base_size = 14) +
         theme(
           plot.title = element_text(size=16, face="bold"),
           axis.title = element_text(size=14, face="bold"),
           axis.text.y = element_text(size=10),
           axis.text.x = element_text(size=12)
         )

    # High-res PNG
    plot_file <- file.path(out_dir, paste0('Barplot_', term_type, '_', contrast, '.png'))
    ggsave(plot_file, p, width=10, height=10, dpi=300)

    # PDF version
    pdf_file <- file.path(out_dir, paste0('Barplot_', term_type, '_', contrast, '.pdf'))
    ggsave(pdf_file, p, width=10, height=10)

    cat("Success!\n")
  }, error = function(e) {
    cat(sprintf("Error creating barplot: %s\n", e$message))
  })
}

# Function to create enrichment dotplot with improved resolution
create_dotplot <- function(enrich_obj, term_type, contrast, out_dir) {
  tryCatch({
    cat(sprintf("Creating dotplot for %s...\n", term_type))
    if(is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) {
      cat("No enrichment results to plot.\n")
      return()
    }

    p <- dotplot(enrich_obj, showCategory=20) +
         ggtitle(paste0("Enriched ", term_type, " Terms: ", contrast)) +
         theme_minimal(base_size = 14) +
         theme(
           plot.title = element_text(size=16, face="bold"),
           axis.title = element_text(size=14, face="bold"),
           axis.text.y = element_text(size=10),
           axis.text.x = element_text(size=12),
           legend.title = element_text(size=12, face="bold")
         )

    # High-res PNG
    plot_file <- file.path(out_dir, paste0('Dotplot_', term_type, '_', contrast, '.png'))
    ggsave(plot_file, p, width=10, height=10, dpi=300)

    # PDF version
    pdf_file <- file.path(out_dir, paste0('Dotplot_', term_type, '_', contrast, '.pdf'))
    ggsave(pdf_file, p, width=10, height=10)

    cat("Success!\n")
  }, error = function(e) {
    cat(sprintf("Error creating dotplot: %s\n", e$message))
  })
}

# Function to create enrichment map with improved resolution
create_emapplot <- function(enrich_obj, term_type, contrast, out_dir) {
  tryCatch({
    cat(sprintf("Creating emapplot for %s...\n", term_type))
    if(is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) {
      cat("No enrichment results to plot.\n")
      return()
    }

    # Only use top 50 terms to avoid overcrowded plots
    if(length(enrich_obj$ID) > 50) {
      enrich_obj_subset <- enrich_obj[1:50]
    } else {
      enrich_obj_subset <- enrich_obj
    }

    emap <- emapplot(pairwise_termsim(enrich_obj_subset), showCategory=30, node_label_size=4,
                    edge_width=0.5, label_format=30) +
            theme(legend.title=element_text(size=14, face="bold"),
                  legend.text=element_text(size=12))

    # High-res PNG
    plot_file <- file.path(out_dir, paste0('Emapplot_', term_type, '_', contrast, '.png'))
    ggsave(plot_file, emap, width=12, height=10, dpi=300)

    # PDF version
    pdf_file <- file.path(out_dir, paste0('Emapplot_', term_type, '_', contrast, '.pdf'))
    ggsave(pdf_file, emap, width=12, height=10)

    cat("Success!\n")
  }, error = function(e) {
    cat(sprintf("Error creating emapplot: %s\n", e$message))
  })
}

# Function to create concept network plot with improved resolution
create_cnetplot <- function(enrich_obj, term_type, contrast, out_dir) {
  tryCatch({
    cat(sprintf("Creating cnetplot for %s...\n", term_type))
    if(is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) {
      cat("No enrichment results to plot.\n")
      return()
    }

    cnet <- cnetplot(enrich_obj, showCategory=10, foldChange=NULL,
                    node_label_size=4, node_label_fontface="bold") +
            theme(legend.title=element_text(size=14, face="bold"),
                 legend.text=element_text(size=12))

    # High-res PNG
    plot_file <- file.path(out_dir, paste0('Cnetplot_', term_type, '_', contrast, '.png'))
    ggsave(plot_file, cnet, width=12, height=10, dpi=300)

    # PDF version
    pdf_file <- file.path(out_dir, paste0('Cnetplot_', term_type, '_', contrast, '.pdf'))
    ggsave(pdf_file, cnet, width=12, height=10)

    cat("Success!\n")
  }, error = function(e) {
    cat(sprintf("Error creating cnetplot: %s\n", e$message))
  })
}

# List of contrasts to analyze
contrast_names <- c(
  "F_OUD_vs_F_None",
  "M_OUD_vs_M_None"
  # "F_OUD_vs_M_OUD",
  # "F_None_vs_M_None"
)

# Process each contrast
for (contrast in contrast_names) {
  cat(sprintf("\n\n======= Processing contrast: %s =======\n", contrast))

  # Create contrast-specific visualization directory
  contrast_viz_dir <- paste0("GSE225158/ORA outputs/visualizations/", contrast)
  dir.create(contrast_viz_dir, recursive = TRUE, showWarnings = FALSE)

  # Load and filter for significant DE genes
  deseq_results_path <- paste0('GSE225158/DESeq2 outputs/deseq2_results_', contrast, '.csv')
  deseq_results <- read_csv(deseq_results_path)

  # Generate MA and Volcano plots
  create_ma_plot(deseq_results, contrast, contrast_viz_dir)
  create_volcano_plot(deseq_results, contrast, contrast_viz_dir)

  # Filter for significant genes
  de_genes <- deseq_results %>%
    filter(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 0.5) %>%
    pull(gene)

  # Remove any NA or empty entries
  de_genes <- de_genes[!is.na(de_genes) & de_genes != ""]
  cat(sprintf("Number of significant DE genes: %d\n", length(de_genes)))

  if (length(de_genes) < 10) {
    cat("Too few DE genes for meaningful ORA analysis. Skipping...\n")
    next
  }

  # Convert gene symbols to Entrez IDs with error handling
  gene_map <- bitr(de_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  cat(sprintf("Original genes: %d, Successfully mapped: %d, Lost: %d (%.1f%%)\n",
            length(de_genes),
            length(unique(gene_map$SYMBOL)),
            length(de_genes) - length(unique(gene_map$SYMBOL)),
            (length(de_genes) - length(unique(gene_map$SYMBOL)))/length(de_genes)*100))

  entrez_ids <- unique(gene_map$ENTREZID)

  # Get all genes as background (universe)
  all_genes <- read_csv('GSE225158/counts.csv')$gene
  all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]

  # Map all genes to Entrez IDs
  all_gene_map <- bitr(all_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  universe <- unique(all_gene_map$ENTREZID)
  cat(sprintf("Background universe size: %d genes\n", length(universe)))

  # Run ORA (GO Biological Process)
  ego_BP <- enrichGO(
    gene = entrez_ids,
    universe = universe,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )

  # Save results
  write.csv(as.data.frame(ego_BP), file = paste0("GSE225158/ORA outputs/ora_results_BP_", contrast, ".csv"), row.names = FALSE)

  # Create enrichment plots for BP
  create_barplot(ego_BP, "BP", contrast, contrast_viz_dir)
  create_dotplot(ego_BP, "BP", contrast, contrast_viz_dir)
  create_emapplot(ego_BP, "BP", contrast, contrast_viz_dir)
  create_cnetplot(ego_BP, "BP", contrast, contrast_viz_dir)

  # Molecular Function
  ego_MF <- enrichGO(
    gene = entrez_ids,
    universe = universe,
    OrgDb = org.Hs.eg.db,
    ont = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  write.csv(as.data.frame(ego_MF), file = paste0("GSE225158/ORA outputs/ora_results_MF_", contrast, ".csv"), row.names = FALSE)

  # Create enrichment plots for MF
  create_barplot(ego_MF, "MF", contrast, contrast_viz_dir)
  create_dotplot(ego_MF, "MF", contrast, contrast_viz_dir)
  create_emapplot(ego_MF, "MF", contrast, contrast_viz_dir)
  create_cnetplot(ego_MF, "MF", contrast, contrast_viz_dir)

  # Cellular Component
  ego_CC <- enrichGO(
    gene = entrez_ids,
    universe = universe,
    OrgDb = org.Hs.eg.db,
    ont = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  write.csv(as.data.frame(ego_CC), file = paste0("GSE225158/ORA outputs/ora_results_CC_", contrast, ".csv"), row.names = FALSE)

  # Create enrichment plots for CC
  create_barplot(ego_CC, "CC", contrast, contrast_viz_dir)
  create_dotplot(ego_CC, "CC", contrast, contrast_viz_dir)
  create_emapplot(ego_CC, "CC", contrast, contrast_viz_dir)
  create_cnetplot(ego_CC, "CC", contrast, contrast_viz_dir)

  # Pathway visualization setup
  gene_map <- bitr(deseq_results$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  de_data <- merge(deseq_results, gene_map, by.x = "gene", by.y = "SYMBOL")
  gene_fc <- de_data$log2FoldChange
  names(gene_fc) <- de_data$ENTREZID

  # Visualize pathways with error handling - using contrast-specific folder
  # Original pathways
  pathview_safe(gene_fc, "hsa04610", "complement", contrast, kegg_dir, contrast_viz_dir)
  pathview_safe(gene_fc, "hsa04080", "neuroligand", contrast, kegg_dir, contrast_viz_dir)

  # Neurotransmitter systems
  pathview_safe(gene_fc, "hsa04728", "dopamine", contrast, kegg_dir, contrast_viz_dir)
  pathview_safe(gene_fc, "hsa04726", "serotonin", contrast, kegg_dir, contrast_viz_dir)
  pathview_safe(gene_fc, "hsa04727", "gaba", contrast, kegg_dir, contrast_viz_dir)
  pathview_safe(gene_fc, "hsa04724", "glutamate", contrast, kegg_dir, contrast_viz_dir)

  # Inflammatory and immune pathways
  pathview_safe(gene_fc, "hsa04668", "TNF", contrast, kegg_dir, contrast_viz_dir)
  pathview_safe(gene_fc, "hsa04060", "cytokine", contrast, kegg_dir, contrast_viz_dir)

  # Signaling pathways
  pathview_safe(gene_fc, "hsa04630", "JAKSTAT", contrast, kegg_dir, contrast_viz_dir)
  pathview_safe(gene_fc, "hsa04010", "MAPK", contrast, kegg_dir, contrast_viz_dir)
  pathview_safe(gene_fc, "hsa04151", "PI3K", contrast, kegg_dir, contrast_viz_dir)

  # Other relevant pathways
  pathview_safe(gene_fc, "hsa04723", "endocannabinoid", contrast, kegg_dir, contrast_viz_dir)
  pathview_safe(gene_fc, "hsa04713", "circadian", contrast, kegg_dir, contrast_viz_dir)
  pathview_safe(gene_fc, "hsa04210", "apoptosis", contrast, kegg_dir, contrast_viz_dir)
  pathview_safe(gene_fc, "hsa04620", "TLR", contrast, kegg_dir, contrast_viz_dir)
}

cat("\nAnalysis complete!\n")