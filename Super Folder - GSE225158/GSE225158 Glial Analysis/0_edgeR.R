library(edgeR)

# Load counts and metadata
counts <- read.csv('/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/counts.csv', row.names=1, check.names=FALSE)
meta <- read.csv('/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/meta.csv', row.names=1, check.names=FALSE)

# Filter only glial cell types
# Identify which entries in meta$celltype3 are glial cells
# Assuming glial cells are marked as "Astro", "Oligo", "OPC", "Microglia" etc.
glial_types <- c("Astrocytes", "Oligos_Pre", "Oligos", "OPCs", "Microglia") # Adjust these names based on your data
is_glial <- meta$celltype3 %in% glial_types

# Filter metadata to keep only glial cells
meta_glial <- meta[is_glial, ]

# Filter counts to keep only columns corresponding to glial cells
counts_glial <- counts[, rownames(meta_glial)]

# Now continue with the standard edgeR workflow using the filtered data
# Reorder meta to match counts columns (should already be aligned, but just to be safe)
meta_glial <- meta_glial[colnames(counts_glial), ]

# Create combined group
meta_glial$Sex_Dx_OUD <- factor(paste(meta_glial$Sex, meta_glial$Dx_OUD, sep="_"))

# Create DGEList object
dge <- DGEList(counts = counts_glial)

# Filter low count genes
keep <- filterByExpr(dge, group=meta_glial$Sex_Dx_OUD)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Calculate normalization factors using TMM method
dge <- calcNormFactors(dge, method="TMM")

# Create design matrix without intercept to get coefficients for each group
design <- model.matrix(~0+Sex_Dx_OUD+Age+PMI, data=meta_glial)

# Rename columns for easier handling
colnames(design) <- gsub("Sex_Dx_OUD", "", colnames(design))

# Estimate dispersion parameters
dge <- estimateDisp(dge, design)

# Fit quasi-likelihood model
fit <- glmQLFit(dge, design)

# Create edgeR outputs directory
dir.create('/Users/aumchampaneri/PycharmProjects/Complement-OUD/Glial Analysis/edgeR Glial outputs/', showWarnings = FALSE, recursive = TRUE)

# Define contrasts
contrast_names <- c(
  "F_OUD_vs_F_None",
  "M_OUD_vs_M_None",
  "F_OUD_vs_M_OUD",
  "F_None_vs_M_None"
)

# Create contrast matrix
my_contrasts <- makeContrasts(
  F_OUD_vs_F_None = F_OUD - F_None,
  M_OUD_vs_M_None = M_OUD - M_None,
  F_OUD_vs_M_OUD = F_OUD - M_OUD,
  F_None_vs_M_None = F_None - M_None,
  levels = design
)

# Loop through contrasts and save results
for (i in seq_along(contrast_names)) {
  # Get contrast name
  contrast_name <- contrast_names[i]

  # Perform test
  qlf <- glmQLFTest(fit, contrast=my_contrasts[, contrast_name])

  # Extract results
  res_df <- topTags(qlf, n=Inf)$table
  res_df$gene <- rownames(res_df)

  # Save results
  out_path <- paste0('/Users/aumchampaneri/PycharmProjects/Complement-OUD/Glial Analysis/edgeR Glial outputs/edgeR_results_', contrast_name, '.csv')
  write.csv(res_df, file=out_path, row.names=FALSE)
}

# Export normalized counts (log CPM values)
norm_counts <- cpm(dge, log=TRUE)
write.csv(norm_counts, file='/Users/aumchampaneri/PycharmProjects/Complement-OUD/Glial Analysis/edgeR Glial outputs/log_cpm_counts.csv')