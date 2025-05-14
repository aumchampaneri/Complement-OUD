# Load required libraries
library(edgeR)
library(limma)
library(readr)
library(dplyr)
library(ggplot2)

# File paths
count_file <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE289002/GSE289002_mouse_raw_counts.csv"
metadata_file <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE289002/mouse_metadata.csv"

# Read data
counts <- read.csv(count_file, row.names = 1)
metadata <- read.csv(metadata_file)

# Ensure sample order matches between counts and metadata
counts <- counts[, metadata$geo_accession]
all(colnames(counts) == metadata$geo_accession) # Should return TRUE

# Create separate analyses for each brain region
for (region in unique(metadata$region)) {
  # Subset data for the current region
  region_meta <- metadata[metadata$region == region, ]
  region_counts <- counts[, region_meta$geo_accession]

  # Create DGEList object
  dge <- DGEList(counts = region_counts)

  # Filter low expressed genes
  keep <- filterByExpr(dge)
  dge <- dge[keep, , keep.lib.sizes=FALSE]

  # Normalize
  dge <- calcNormFactors(dge, method = "TMM")

  # Create design matrix with treatment and sex as factors
  region_meta$treatment <- factor(region_meta$treatment,
                                levels = c("Sal", "Mor + 24h", "Mor + 2W", "Chronic mor"))
  region_meta$sex <- factor(region_meta$sex)

  # Create design matrix: ~ treatment + sex + treatment:sex
  design <- model.matrix(~ treatment + sex + treatment:sex, data = region_meta)

  # Estimate dispersion
  dge <- estimateDisp(dge, design)

  # Fit model
  fit <- glmQLFit(dge, design)

  # Perform comparisons
  # 1. Each treatment vs. Saline
  contrasts <- list(
    "Mor24h_vs_Saline" = makeContrasts(treatmentMor + 24h, levels = design),
    "Mor2W_vs_Saline" = makeContrasts(treatmentMor + 2W, levels = design),
    "ChronicMor_vs_Saline" = makeContrasts(treatmentChronic mor, levels = design),
    "Female_vs_Male" = makeContrasts(sexfemale, levels = design)
  )

  # Interaction contrasts
  interactions <- list(
    "Mor24h_Sex_Interaction" = makeContrasts(`treatmentMor + 24h:sexfemale`, levels = design),
    "Mor2W_Sex_Interaction" = makeContrasts(`treatmentMor + 2W:sexfemale`, levels = design),
    "ChronicMor_Sex_Interaction" = makeContrasts(`treatmentChronic mor:sexfemale`, levels = design)
  )

  # Combine all contrasts
  all_contrasts <- c(contrasts, interactions)

  # Results directory
  results_dir <- paste0("/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE289002/DE_results_", region)
  dir.create(results_dir, showWarnings = FALSE)

  # Perform differential expression for each contrast
  for (contrast_name in names(all_contrasts)) {
    contrast <- all_contrasts[[contrast_name]]

    # Test for differential expression
    qlf <- glmQLFTest(fit, contrast = contrast)

    # Extract results
    res <- topTags(qlf, n = Inf, sort.by = "PValue")
    res_df <- as.data.frame(res)

    # Add gene names as a column
    res_df$gene_id <- rownames(res_df)

    # Save results to CSV
    write.csv(res_df, file.path(results_dir, paste0(contrast_name, ".csv")))

    # Volcano plot
    png(file.path(results_dir, paste0(contrast_name, "_volcano.png")), width = 800, height = 600)
    with(res_df, plot(logFC, -log10(FDR),
                     pch = 20, main = contrast_name,
                     xlab = "Log2 Fold Change", ylab = "-log10(FDR)"))
    with(subset(res_df, FDR < 0.05 & abs(logFC) > 1),
         points(logFC, -log10(FDR), pch = 20, col = "red"))
    abline(h = -log10(0.05), col = "blue", lty = 2)
    abline(v = c(-1, 1), col = "blue", lty = 2)
    dev.off()
  }

  # Generate MDS plot for sample clustering
  png(file.path(results_dir, "MDS_plot.png"), width = 1000, height = 800)
  plotMDS(dge, col = as.numeric(region_meta$treatment) + 4 * as.numeric(region_meta$sex),
          labels = paste(region_meta$treatment, region_meta$sex, sep = "_"))
  legend("topright", legend = c(levels(region_meta$treatment), levels(region_meta$sex)),
         col = c(1:4, 5:8), pch = 19)
  dev.off()
}

# Look specifically at complement-related genes (if interested)
complement_genes <- c("C1qa", "C1qb", "C1qc", "C3", "C4b", "Cfb", "Cfd")
# You'll need to map these to the actual IDs in your dataset, possibly using a gene annotation database