# R

if (!requireNamespace("pathview", quietly = TRUE)) BiocManager::install("pathview")
if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) BiocManager::install("org.Mm.eg.db")
library(pathview)
library(readr)
library(dplyr)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(tibble)
library(stringr)

de_files <- c(
  "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE289002/DE_results/sex_treatment/Chronic_vs_Sal_female.csv",
  "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE289002/DE_results/sex_treatment/Chronic_vs_Sal_male.csv",
  "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE289002/DE_results/sex_treatment/Mor2W_vs_Sal_female.csv",
  "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE289002/DE_results/sex_treatment/Mor2W_vs_Sal_male.csv",
  "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE289002/DE_results/sex_treatment/Mor24h_vs_Sal_female.csv",
  "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE289002/DE_results/sex_treatment/Mor24h_vs_Sal_male.csv"
)

for (de_file in de_files) {
  cat("Processing:", de_file, "\n")
  de_results <- read_csv(de_file, show_col_types = FALSE)
  de_results <- de_results %>% column_to_rownames(var = colnames(de_results)[1])
  ensembl_ids <- rownames(de_results)
  entrez_map <- AnnotationDbi::select(org.Mm.eg.db, keys = ensembl_ids, keytype = "ENSEMBL", columns = "ENTREZID")
  de_results$ENTREZID <- entrez_map$ENTREZID[match(ensembl_ids, entrez_map$ENSEMBL)]
  gene_data <- de_results$logFC
  names(gene_data) <- de_results$ENTREZID
  gene_data <- gene_data[!is.na(names(gene_data))]
  suffix <- str_replace(basename(de_file), "\\.csv$", "")
  pv_out <- pathview(
    gene.data = gene_data,
    pathway.id = "mmu04610",
    species = "mmu",
    out.suffix = suffix
  )
}
cat("All Pathview plots generated.\n")