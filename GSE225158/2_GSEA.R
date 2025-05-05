# Load libraries
library(clusterProfiler)
library(msigdbr)
library(dplyr)
library(readr)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Load DESeq2 results
res <- read_csv('GSE225158/DESeq2 outputs/deseq2_results_M_OUD_vs_M_None.csv')

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

# Break ties in gene statistics by adding minor random variation
set.seed(42)  # For reproducibility
gene_list_ties_broken <- gene_list + rnorm(length(gene_list), 0, 1e-10)
names(gene_list_ties_broken) <- names(gene_list)
gene_list_ties_broken <- sort(gene_list_ties_broken, decreasing = TRUE)

# Get all MSigDB gene sets (all collections, Entrez IDs)
msig_all <- msigdbr(species = "Homo sapiens") %>%
  dplyr::select(gs_name, entrez_gene = ncbi_gene) %>%
  filter(!is.na(entrez_gene))

# Run GSEA with improved parameters
gsea_res <- GSEA(
  geneList = gene_list_ties_broken,  # Use the version with ties broken
  TERM2GENE = msig_all,
  pvalueCutoff = 0.05,
  minGSSize = 5,
  maxGSSize = 5000,
  nPermSimple = 100000,  # Increased for better p-value estimation
  eps = 0  # Better estimate very small p-values
)

# Save results
write.csv(gsea_res@result, 'GSE225158/GSEA outputs/gsea_results_M_OUD_vs_M_None.csv', row.names = FALSE)

# Save GSEA object for visualization
saveRDS(gsea_res, file.path("GSE225158/GSEA outputs", "gsea_obj_M_OUD_vs_M_None.rds"))

# Check column names to find adjusted p-value column
gsea_cols <- colnames(gsea_res@result)
print(paste("Available columns:", paste(gsea_cols, collapse=", ")))

# Filter for inflammatory and complement pathways
# Read complement genes from CSV file
complement_genes_df <- read.csv("GSE225158/KEGG outputs/kegg_complement_unique_genes.csv", stringsAsFactors = FALSE)
complement_genes <- complement_genes_df$gene  # Adjust column name if needed

# Get Entrez IDs for complement genes
complement_entrez <- mapIds(
  org.Hs.eg.db,
  keys = complement_genes,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)
complement_entrez <- complement_entrez[!is.na(complement_entrez)]
complement_entrez <- as.character(complement_entrez)

# Filter pathways by keywords or complement genes in core enrichment
inflammatory_complement_pathways <- gsea_res@result %>%
  filter(
    # Original keyword filter
    grepl("INFLAMMATORY|INFLAMMATION|COMPLEMENT|CYTOKINE|IMMUNE|IL|TNF|INTERFERON",
          Description, ignore.case = TRUE) |
    # New filter: check if any complement genes are in core enrichment
    sapply(core_enrichment, function(genes) {
      if(is.na(genes)) return(FALSE)
      gene_list <- unlist(strsplit(genes, "/"))
      return(any(gene_list %in% complement_entrez))
    })
  )

# Use p.adjust column instead of padj (based on standard clusterProfiler output)
if("p.adjust" %in% colnames(gsea_res@result)) {
  inflammatory_complement_pathways <- inflammatory_complement_pathways %>%
    arrange(p.adjust)
} else if("adj.P.Val" %in% colnames(gsea_res@result)) {
  inflammatory_complement_pathways <- inflammatory_complement_pathways %>%
    arrange(adj.P.Val)
} else {
  # If we can't find adjusted p-value column, just keep the order as is
  print("Could not find adjusted p-value column. Results not sorted by significance.")
}

# Save inflammatory and complement pathways
write.csv(inflammatory_complement_pathways,
          'GSE225158/GSEA outputs/gsea_M_inflammatory_complement_pathways.csv',
          row.names = FALSE)