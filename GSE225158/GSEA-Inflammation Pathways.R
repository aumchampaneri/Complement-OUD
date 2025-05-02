# Load required library
library(dplyr)

# Read GSEA results
gsea_res <- read.csv('GSE225158/gsea_results_M_OUD_vs_M_None.csv', stringsAsFactors = FALSE)

# Read gene names from the CSV file
gene_list <- read.csv('GSE225158 Gene Library/gene_names.csv', stringsAsFactors = FALSE)$Gene.Name

# Define keywords for inflammatory and complement pathways
keywords <- c(
  "INFLAM", "INTERLEUKIN", "CYTOKINE", "TNF", "NF_KB", "IL", "TOLL_LIKE", "INTERFERON",
  "COMPLEMENT", "C3", "C4", "C5", "C1Q", "C1R", "C1S", "C2", "C6", "C7", "C8", "C9"
)

# Filter for pathways containing any keyword (case-insensitive)
focus_res <- gsea_res %>%
  dplyr::filter(grepl(paste(keywords, collapse = "|"), ID, ignore.case = TRUE))

# View top results
print(head(focus_res))

# Save filtered results to a new CSV
write.csv(focus_res, 'GSE225158/gsea_inflammatory_complement_pathways.csv', row.names = FALSE)