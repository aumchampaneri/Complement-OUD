library(dplyr)
library(stringr)
library(readr)

# Read genes of interest and filtered GSEA results
gene_df <- read_csv('GSE225158 Gene Library/gene_names.csv')
genes_of_interest <- gene_df$`Gene Name`
gsea_res <- read_csv('GSE225158/gsea_inflammatory_complement_pathways.csv')

# Map gene symbols to Entrez IDs (if needed)
library(org.Hs.eg.db)
library(AnnotationDbi)
entrez_ids <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = genes_of_interest,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)
entrez_ids <- na.omit(entrez_ids) %>% as.character()

# Filter pathways containing your genes of interest
gsea_res$contains_gene <- sapply(
  gsea_res$core_enrichment,
  function(x) any(unlist(str_split(as.character(x), "/")) %in% entrez_ids)
)
focus_pathways <- gsea_res %>% filter(contains_gene)

# Save or view the results
write_csv(focus_pathways, 'GSE225158/pathways_with_genes_of_interest.csv')