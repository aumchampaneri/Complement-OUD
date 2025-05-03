library(dplyr)
library(stringr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(readr)

# Read GSEA results
gsea_res <- read_csv('GSE225158/gsea_results_M_OUD_vs_M_None.csv')

# Read gene names and check column names
gene_df <- read_csv('GSE225158 Gene Library/gene_names.csv')
print(colnames(gene_df))  # Check actual column names

# Use the correct column name (adjust as needed)
gene_list <- gene_df$`Gene Name`

# Continue as before
entrez_ids <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = gene_list,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)
entrez_ids <- na.omit(entrez_ids) %>% as.character()

# Define keywords
keywords <- c(
  "INFLAM", "INTERLEUKIN", "CYTOKINE", "TNF", "NF_KB", "IL", "TOLL_LIKE", "INTERFERON"
)

# Function to check gene overlap
has_gene_overlap <- function(core_enrichment, entrez_ids) {
  pathway_ids <- unlist(str_split(as.character(core_enrichment), "/"))
  any(pathway_ids %in% entrez_ids)
}

# Filter for pathways matching keywords or gene overlap
focus_res <- gsea_res %>%
  filter(
    grepl(paste(keywords, collapse = "|"), ID, ignore.case = TRUE) | # Uncoment this line to filter by keywords
      mapply(has_gene_overlap, core_enrichment, MoreArgs = list(entrez_ids = entrez_ids))
  )

# Save results
write_csv(focus_res, 'GSE225158/gsea_inflammatory_complement_pathways.csv')