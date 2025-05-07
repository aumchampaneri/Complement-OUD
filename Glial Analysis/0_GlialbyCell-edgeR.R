library(edgeR)

# Load counts and metadata
counts <- read.csv('/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/counts.csv', row.names=1, check.names=FALSE)
meta <- read.csv('/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/meta.csv', row.names=1, check.names=FALSE)

# Glial subtypes to process
glial_types <- c("Microglia", "Astrocytes", "OPCs", "Oligos", "Oligos_Pre")

# Output folder
base_out_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/Glial Analysis/edge R Per_Celltype"

for (glial in glial_types) {
  cat("Processing:", glial, "\n")

  # Subset metadata and counts
  meta_sub <- meta[meta$celltype3 == glial, ]
  counts_sub <- counts[, rownames(meta_sub)]
  meta_sub <- meta_sub[colnames(counts_sub), ]

  # Create group variable
  meta_sub$Sex_Dx_OUD <- factor(paste(meta_sub$Sex, meta_sub$Dx_OUD, sep = "_"))

  # Create DGEList object
  dge <- DGEList(counts = counts_sub)

  # Filter low-expression genes
  keep <- filterByExpr(dge, group = meta_sub$Sex_Dx_OUD)
  dge <- dge[keep, , keep.lib.sizes = FALSE]

  # TMM normalization
  dge <- calcNormFactors(dge, method = "TMM")

  # Design matrix
  design <- model.matrix(~0 + Sex_Dx_OUD + Age + PMI, data = meta_sub)
  colnames(design) <- gsub("Sex_Dx_OUD", "", colnames(design))

  # Estimate dispersion
  dge <- estimateDisp(dge, design)

  # Fit model
  fit <- glmQLFit(dge, design)

  # Output directory
  out_dir <- file.path(base_out_dir, glial)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # Define and run only within-sex OUD contrasts
  contrast_names <- c("F_OUD_vs_F_None", "M_OUD_vs_M_None")
  my_contrasts <- makeContrasts(
    F_OUD_vs_F_None = F_OUD - F_None,
    M_OUD_vs_M_None = M_OUD - M_None,
    levels = design
  )

  for (i in seq_along(contrast_names)) {
    contrast_name <- contrast_names[i]
    qlf <- glmQLFTest(fit, contrast = my_contrasts[, contrast_name])
    res_df <- topTags(qlf, n = Inf)$table
    res_df$gene <- rownames(res_df)
    write.csv(res_df, file = file.path(out_dir, paste0("edgeR_results_", contrast_name, ".csv")), row.names = FALSE)
  }

  # Export normalized logCPM
  log_cpm <- cpm(dge, log = TRUE)
  write.csv(log_cpm, file = file.path(out_dir, "log_cpm_counts.csv"))
}
