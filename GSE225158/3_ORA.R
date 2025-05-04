# Consolidated ORA and Visualization Script

library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
library(enrichplot)
library(ggplot2)
library(readr)
library(dplyr)

# Create output directories
dir.create("GSE225158/ORA outputs", recursive = TRUE, showWarnings = FALSE)
dir.create("GSE225158/ORA outputs/visualizations", recursive = TRUE, showWarnings = FALSE)

# List of contrasts to analyze
contrast_names <- c(
  "F_OUD_vs_F_None",
  "M_OUD_vs_M_None",
  "F_OUD_vs_M_OUD",
  "F_None_vs_M_None"
)

# Process each contrast
for (contrast in contrast_names) {
  cat(sprintf("\n\n======= Processing contrast: %s =======\n", contrast))
  
  # Load and filter for significant DE genes
  deseq_results_path <- paste0('GSE225158/DESeq2 outputs/deseq2_results_', contrast, '.csv')
  deseq_results <- read_csv(deseq_results_path)
  
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
  
  # Visualizations - with null checks
  if (!is.null(ego_BP) && nrow(ego_BP) > 0) {
    dotplot_BP <- dotplot(ego_BP, showCategory = 15, title = paste0("Biological Process - ", contrast))
    ggsave(paste0("GSE225158/ORA outputs/visualizations/dotplot_BP_", contrast, ".png"), dotplot_BP, width = 10, height = 8)
    
    emap_BP <- emapplot(pairwise_termsim(ego_BP), showCategory = min(15, nrow(ego_BP)))
    ggsave(paste0("GSE225158/ORA outputs/visualizations/emap_BP_", contrast, ".png"), emap_BP, width = 10, height = 8)
    
    cnet_BP <- cnetplot(ego_BP, showCategory = min(5, nrow(ego_BP)), foldChange = NULL)
    ggsave(paste0("GSE225158/ORA outputs/visualizations/cnet_BP_", contrast, ".png"), cnet_BP, width = 12, height = 10)
    
    barplot_BP <- barplot(ego_BP, showCategory = 15)
    ggsave(paste0("GSE225158/ORA outputs/visualizations/barplot_BP_", contrast, ".png"), barplot_BP, width = 10, height = 8)
  }
  
  if (!is.null(ego_MF) && nrow(ego_MF) > 0) {
    dotplot_MF <- dotplot(ego_MF, showCategory = 15, title = paste0("Molecular Function - ", contrast))
    ggsave(paste0("GSE225158/ORA outputs/visualizations/dotplot_MF_", contrast, ".png"), dotplot_MF, width = 10, height = 8)
    
    if (nrow(ego_MF) > 1) {
      emap_MF <- emapplot(pairwise_termsim(ego_MF), showCategory = min(15, nrow(ego_MF)))
      ggsave(paste0("GSE225158/ORA outputs/visualizations/emap_MF_", contrast, ".png"), emap_MF, width = 10, height = 8)
    }
    
    cnet_MF <- cnetplot(ego_MF, showCategory = min(5, nrow(ego_MF)), foldChange = NULL)
    ggsave(paste0("GSE225158/ORA outputs/visualizations/cnet_MF_", contrast, ".png"), cnet_MF, width = 12, height = 10)
  }
  
  if (!is.null(ego_CC) && nrow(ego_CC) > 0) {
    dotplot_CC <- dotplot(ego_CC, showCategory = 15, title = paste0("Cellular Component - ", contrast))
    ggsave(paste0("GSE225158/ORA outputs/visualizations/dotplot_CC_", contrast, ".png"), dotplot_CC, width = 10, height = 8)
  }
  
  # MA and Volcano Plots
  # MA Plot
  res_ma <- deseq_results[deseq_results$baseMean > 0, ]
  png(paste0("GSE225158/ORA outputs/visualizations/MA_plot_", contrast, ".png"), width=800, height=600)
  with(res_ma, plot(baseMean, log2FoldChange, pch=20, cex=0.5, log='x',
       main=paste('MA Plot:', contrast), xlab='Mean Expression', ylab='log2 Fold Change'))
  abline(h=0, col='red')
  dev.off()
  
  # Volcano Plot
  deseq_results$padj[is.na(deseq_results$padj)] <- 1
  deseq_results$significant <- deseq_results$padj < 0.05 & abs(deseq_results$log2FoldChange) > 0.5
  png(paste0("GSE225158/ORA outputs/visualizations/Volcano_plot_", contrast, ".png"), width=800, height=600)
  p <- ggplot(deseq_results, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
    geom_point(alpha=0.6) +
    scale_color_manual(values=c('grey','red')) +
    theme_minimal() +
    labs(title=paste('Volcano Plot:', contrast), x='log2 Fold Change', y='-log10 Adjusted p-value')
  print(p)
  dev.off()

  # Pathway visualization - focus on complement pathway
  gene_map <- bitr(deseq_results$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  de_data <- merge(deseq_results, gene_map, by.x = "gene", by.y = "SYMBOL")
  gene_fc <- de_data$log2FoldChange
  names(gene_fc) <- de_data$ENTREZID

  # Create a directory for KEGG data if it doesn't exist
  kegg_dir <- "GSE225158/ORA outputs/visualizations/kegg_data"
  dir.create(kegg_dir, recursive = TRUE, showWarnings = FALSE)

  # Visualize complement pathway
  pathview(
      gene.data = gene_fc,
      pathway.id = "hsa04610",  # Complement and coagulation cascades
      species = "hsa",
      out.suffix = paste0(contrast, "_complement"),
      kegg.dir = kegg_dir,
      out.dir = "GSE225158/ORA outputs/visualizations/"
  )

  # Add neuroactive ligand pathway
  pathview(
      gene.data = gene_fc,
      pathway.id = "hsa04080",  # Neuroactive ligand-receptor interaction
      species = "hsa",
      out.suffix = paste0(contrast, "_neuroligand"),
      kegg.dir = kegg_dir,
      out.dir = "GSE225158/ORA outputs/visualizations/"
  )
}