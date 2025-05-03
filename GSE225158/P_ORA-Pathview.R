library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)

# Prepare gene data: named vector of log2FC
de_data <- read.csv("GSE225158/DESeq2 outputs/deseq2_results_M_OUD_vs_M_None.csv")
gene_map <- bitr(de_data$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
de_data <- merge(de_data, gene_map, by.x = "gene", by.y = "SYMBOL")
gene_fc <- de_data$log2FoldChange
names(gene_fc) <- de_data$ENTREZID

# Visualize the pathway
pathview(
    gene.data = gene_fc,
    pathway.id = "hsa04610",
    species = "hsa",
    out.suffix = "M_OUD_vs_M_None"
)