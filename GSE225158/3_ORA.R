library(clusterProfiler)
library(org.Hs.eg.db)
library(readr)
library(dplyr)

# Load DE genes (assumes a 'gene' column)
de_genes <- read_csv('/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/DESeq2 outputs/deseq2_results_M_OUD_vs_M_None.csv')$gene

# Remove any NA or empty entries
de_genes <- de_genes[!is.na(de_genes) & de_genes != ""]

# Convert gene symbols to Entrez IDs
gene_map <- bitr(de_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
entrez_ids <- unique(gene_map$ENTREZID)

# Run ORA (GO Biological Process)
ego <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

# Save results
write.csv(as.data.frame(ego), file = "GSE225158/ORA outputs/ora_results.csv", row.names = FALSE)