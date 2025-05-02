# R

library(DESeq2)

# Load counts and metadata
counts <- read.csv('/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/counts.csv', row.names=1, check.names=FALSE)
meta <- read.csv('/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/meta.csv', row.names=1, check.names=FALSE)

# Reorder meta to match counts columns
meta <- meta[colnames(counts), ]

# Create combined group
meta$Sex_Dx_OUD <- factor(paste(meta$Sex, meta$Dx_OUD, sep="_"))

# Set up DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData=counts, colData=meta, design=~Sex_Dx_OUD+Age+PMI)
dds <- DESeq(dds)

# List of contrasts: c("group", "level1", "level2")
contrasts <- list(
  c("Sex_Dx_OUD", "F_OUD", "F_None"),
  c("Sex_Dx_OUD", "M_OUD", "M_None"),
  c("Sex_Dx_OUD", "F_OUD", "M_OUD"),
  c("Sex_Dx_OUD", "F_None", "M_None")
)
contrast_names <- c(
  "F_OUD_vs_F_None",
  "M_OUD_vs_M_None",
  "F_OUD_vs_M_OUD",
  "F_None_vs_M_None"
)

# Loop through contrasts and save results
for (i in seq_along(contrasts)) {
  res <- results(dds, contrast=contrasts[[i]])
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res)  # <-- Use rownames from DESeq2 results
  out_path <- paste0('/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/deseq2_results_', contrast_names[i], '.csv')
  write.csv(res_df, file=out_path, row.names=FALSE)
}