# R
library(edgeR)
library(pheatmap)
library(org.Mm.eg.db)
library(AnnotationDbi)

gene_list_path <- '/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE289002/kegg_outputs/GSE207128_mouse_complement_cascade_genes.csv'
expr_path <- '/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE289002/QC/dge_filtered_normalized.rds'
metadata_path <- '/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE289002/QC/metadata.rds'
heatmap_path <- '/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE289002/DE_results/sex_treatment/plots/complement_genes_heatmap.png'

# 1. Read complement gene list
complement_genes <- toupper(read.csv(gene_list_path)$gene)

# 2. Load expression matrix and metadata
dge <- readRDS(expr_path)
metadata <- readRDS(metadata_path)
dge <- dge[, metadata$title]

# 3. Map Ensembl IDs to gene symbols
ensembl_ids <- rownames(dge)
gene_symbols <- mapIds(org.Mm.eg.db, keys=ensembl_ids, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
expr_matrix <- cpm(dge, log=TRUE)
rownames(expr_matrix) <- toupper(gene_symbols)

# 4. Subset to complement genes
expr_subset <- expr_matrix[rownames(expr_matrix) %in% complement_genes & !is.na(rownames(expr_matrix)), ]

# 5. Z-score scaling (row-wise)
expr_z <- t(scale(t(expr_subset)))

# 6. Aggregate by treatment and sex
metadata$group <- paste(metadata$sex, metadata$treatment, sep="_")
group_means <- aggregate(t(expr_z), by=list(metadata$group), FUN=mean)
rownames(group_means) <- group_means$Group.1
group_means <- t(group_means[, -1])

# 7. Build descriptive group labels
treatment_desc <- c(
  "Sal" = "Saline control",
  "Mor + 24h" = "Morphine 24h",
  "Mor + 2W" = "Morphine 2W",
  "Chronic mor" = "Chronic morphine"
)
group_counts <- table(metadata$group)
group_labels <- sapply(colnames(group_means), function(g) {
  parts <- unlist(strsplit(g, "_"))
  sex <- parts[1]
  treatment <- paste(parts[-1], collapse="_")
  desc <- treatment_desc[treatment]
  n <- group_counts[g]
  paste0(sex, ": ", desc, " (n=", n, ")")
})
colnames(group_means) <- group_labels

# 8. Plot heatmap
pheatmap(
  group_means,
  color = colorRampPalette(c("navy", "white", "firebrick"))(50),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  main = "Complement Gene Expression (z-score)",
  filename = heatmap_path
)