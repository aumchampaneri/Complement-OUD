# Load required libraries
library(edgeR)
library(limma)

# Load preprocessed data
qc_dir <- '/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE289002/QC'
dge <- readRDS(file.path(qc_dir, 'dge_filtered_normalized.rds'))
metadata <- readRDS(file.path(qc_dir, 'metadata.rds'))

# Remove region as a factor (ignore region)
metadata$sex <- factor(metadata$sex)
metadata$treatment <- factor(metadata$treatment, levels = c('Sal', 'Mor + 24h', 'Mor + 2W', 'Chronic mor'))

# Ensure sample order matches
dge <- dge[, metadata$title]

# Design matrix: sex, treatment, and interaction
design <- model.matrix(~ sex * treatment, data = metadata)
colnames(design) <- make.names(colnames(design))

# Estimate dispersion and fit model
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)

# Example contrasts:
# 1. Treatment effect within each sex
contrasts <- list(
  'Mor24h_vs_Sal_female' = makeContrasts(treatmentMor...24h, levels=design),
  'Mor24h_vs_Sal_male'   = makeContrasts(treatmentMor...24h + sexmale.treatmentMor...24h, levels=design),
  'Mor2W_vs_Sal_female'  = makeContrasts(treatmentMor...2W, levels=design),
  'Mor2W_vs_Sal_male'    = makeContrasts(treatmentMor...2W + sexmale.treatmentMor...2W, levels=design),
  'Chronic_vs_Sal_female' = makeContrasts(treatmentChronic.mor, levels=design),
  'Chronic_vs_Sal_male'   = makeContrasts(treatmentChronic.mor + sexmale.treatmentChronic.mor, levels=design)
)

# 2. Sex difference within each treatment
sexdiff_contrasts <- list(
  'SexDiff_Sal'         = makeContrasts(sexmale, levels=design),
  'SexDiff_Mor24h'      = makeContrasts(sexmale + sexmale.treatmentMor...24h, levels=design),
  'SexDiff_Mor2W'       = makeContrasts(sexmale + sexmale.treatmentMor...2W, levels=design),
  'SexDiff_Chronic'     = makeContrasts(sexmale + sexmale.treatmentChronic.mor, levels=design)
)

# Output directory
results_dir <- '/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE289002/DE_results/sex_treatment'
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Run all contrasts and save results
for (name in names(contrasts)) {
  qlf <- glmQLFTest(fit, contrast=contrasts[[name]])
  res <- topTags(qlf, n=Inf, sort.by='PValue')
  write.csv(as.data.frame(res), file.path(results_dir, paste0(name, '.csv')), row.names=TRUE)
}

for (name in names(sexdiff_contrasts)) {
  qlf <- glmQLFTest(fit, contrast=sexdiff_contrasts[[name]])
  res <- topTags(qlf, n=Inf, sort.by='PValue')
  write.csv(as.data.frame(res), file.path(results_dir, paste0(name, '.csv')), row.names=TRUE)
}

cat('edgeR analysis complete. Results saved in:', results_dir, '\n')