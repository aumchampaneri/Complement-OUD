# R
library(readr)
library(dplyr)
library(tibble)
library(stringr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)

# 1. Read GSEA results
gsea <- read_csv('GSEA outputs/gsea_results_M_OUD_vs_M_None.csv')

# 2. Choose pathway of interest
pathway_id <- "LAKE_ADULT_KIDNEY_C19_COLLECTING_DUCT_INTERCALATED_CELLS_TYPE_A_MEDULLA"

# 3. Extract core_enrichment string for the pathway
core_enrichment_str <- gsea %>%
  filter(ID == pathway_id) %>%
  pull(core_enrichment)

if (length(core_enrichment_str) == 0 || is.na(core_enrichment_str)) {
  stop("Pathway not found or core_enrichment is missing.")
}

# 4. Extract and clean Entrez IDs
entrez_ids <- unlist(str_split(as.character(core_enrichment_str), "/"))
entrez_ids <- str_trim(entrez_ids)
entrez_ids <- entrez_ids[entrez_ids != ""]
entrez_ids <- entrez_ids[grepl("^[0-9]+$", entrez_ids)]

if (length(entrez_ids) == 0) {
  stop("No valid Entrez IDs found in core_enrichment.")
}

# 5. Map Entrez IDs to gene symbols
gene_symbols <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = entrez_ids,
  columns = c("SYMBOL"),
  keytype = "ENTREZID"
) %>% pull(SYMBOL) %>% unique()
gene_symbols <- gene_symbols[!is.na(gene_symbols)]

if (length(gene_symbols) == 0) {
  stop("No gene symbols mapped from Entrez IDs.")
}

# 6. Read DESeq2 results
expr <- read_csv('DESeq2 outputs/deseq2_results_M_OUD_vs_M_None.csv')

# 7. Filter for leading edge genes (now using gene symbols)
barplot_data <- expr %>%
  filter(gene %in% gene_symbols) %>%
  dplyr::select(gene, log2FoldChange)

if (nrow(barplot_data) == 0) {
  stop("No matching leading edge genes found in DESeq2 results after mapping.")
}

# 8. Barplot
p <- ggplot(barplot_data, aes(x = reorder(gene, log2FoldChange), y = log2FoldChange)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = paste("log2FoldChange of Leading Edge Genes for", pathway_id),
    x = "Gene",
    y = "log2FoldChange"
  ) +
  theme_minimal()

ggsave(
  '/Users/aumchampaneri/PycharmProjects/Complement-OUD/leading_edge_genes_barplot.png',
  p + theme(plot.margin = margin(5, 5, 5, 20)), # increase left margin
  width = 20, height = 13 # increase width
)