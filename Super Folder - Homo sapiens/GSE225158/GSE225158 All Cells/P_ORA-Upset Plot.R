library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# Load gene list
de_genes <- read.csv('GSE225158/DESeq2 outputs/deseq2_results_F_OUD_vs_F_None.csv')$gene
de_genes <- de_genes[!is.na(de_genes) & de_genes != ""]

# Map to Entrez IDs
gene_map <- bitr(de_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
entrez_ids <- unique(gene_map$ENTREZID)

# Run enrichment
ego <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

# Upset plot
png("ora_upsetplot.png", width = 1400, height = 1000)
upsetplot(ego)
dev.off()