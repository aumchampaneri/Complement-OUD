library(dplyr)
library(stringr)
library(readr)

# Read genes of interest and filtered GSEA results
gene_df <- read_csv('GSE225158/kegg_unique_genes.csv')
genes_of_interest <- gene_df$`gene`
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

# Count overlap for each pathway
gsea_res$num_genes_of_interest <- sapply(
  gsea_res$core_enrichment,
  function(x) sum(unlist(str_split(as.character(x), "/")) %in% entrez_ids)
)

# Find the pathway(s) with the highest count
max_count <- max(gsea_res$num_genes_of_interest)
max_rows <- gsea_res %>% filter(num_genes_of_interest == max_count)

# For each such pathway, list the genes of interest present
max_rows$genes_in_pathway <- sapply(
  max_rows$core_enrichment,
  function(x) {
    entrez_in_pathway <- unlist(str_split(as.character(x), "/"))
    matched_entrez <- intersect(entrez_in_pathway, entrez_ids)
    # Map back to gene symbols
    if (length(matched_entrez) > 0) {
      symbols <- AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys = matched_entrez,
        column = "SYMBOL",
        keytype = "ENTREZID",
        multiVals = "first"
      )
      paste(na.omit(symbols), collapse = ", ")
    } else {
      NA
    }
  }
)

# Print pathway description, count, and gene symbols
print(as.data.frame(max_rows %>% dplyr::select(Description, num_genes_of_interest, genes_in_pathway)), right = FALSE, row.names = FALSE)