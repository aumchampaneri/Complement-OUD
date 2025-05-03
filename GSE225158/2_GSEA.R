# Load libraries
library(clusterProfiler)
library(msigdbr)
library(dplyr)
library(readr)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Load DESeq2 results
res <- read_csv('DESeq2 outputs/deseq2_results_M_OUD_vs_M_None.csv')

# Map gene symbols to Entrez IDs
res$entrez <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = res$gene,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

# Filter for valid Entrez IDs and stats
res <- res %>%
  filter(!is.na(entrez), !duplicated(entrez), !is.na(stat))

# Prepare ranked gene list (Entrez IDs)
gene_list <- res$stat
names(gene_list) <- as.character(res$entrez)
gene_list <- sort(gene_list, decreasing = TRUE)

# Get all MSigDB gene sets (all collections, Entrez IDs)
msig_all <- msigdbr(species = "Homo sapiens") %>%
  dplyr::select(gs_name, entrez_gene = ncbi_gene) %>%
  filter(!is.na(entrez_gene))

# Run GSEA
gsea_res <- GSEA(
  geneList = gene_list,
  TERM2GENE = msig_all,
  pvalueCutoff = 0.05,
  minGSSize = 5,
  maxGSSize = 5000,
  nPermSimple = 10000
)

# Save results
write.csv(gsea_res@result, 'GSE225158/GSEA outputs/gsea_results_M_OUD_vs_M_None.csv', row.names = FALSE)