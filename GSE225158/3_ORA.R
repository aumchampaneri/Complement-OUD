library(clusterProfiler)
library(org.Hs.eg.db)
library(readr)
library(dplyr)

# Create output directory if it doesn't exist
dir.create("GSE225158/ORA outputs", recursive = TRUE, showWarnings = FALSE)

# Load and filter for significant DE genes
deseq_results <- read_csv('GSE225158/DESeq2 outputs/deseq2_results_M_OUD_vs_M_None.csv')

# Filter for significant genes (typically padj < 0.05 and abs(log2FoldChange) > some threshold)
de_genes <- deseq_results %>%
  filter(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 0.5) %>%
  pull(gene)

# Remove any NA or empty entries
de_genes <- de_genes[!is.na(de_genes) & de_genes != ""]
cat(sprintf("Number of significant DE genes: %d\n", length(de_genes)))

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

# Run ORA (GO Biological Process) with explicit parameters
ego <- enrichGO(
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
write.csv(as.data.frame(ego), file = "GSE225158/ORA outputs/ora_results_BP.csv", row.names = FALSE)

# Run additional ontologies
# Molecular Function
egoMF <- enrichGO(
  gene = entrez_ids,
  universe = universe,
  OrgDb = org.Hs.eg.db,
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)
write.csv(as.data.frame(egoMF), file = "GSE225158/ORA outputs/ora_results_MF.csv", row.names = FALSE)

# Cellular Component
egoCC <- enrichGO(
  gene = entrez_ids,
  universe = universe,
  OrgDb = org.Hs.eg.db,
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)
write.csv(as.data.frame(egoCC), file = "GSE225158/ORA outputs/ora_results_CC.csv", row.names = FALSE)