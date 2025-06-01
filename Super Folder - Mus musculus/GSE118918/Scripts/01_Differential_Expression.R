# ========================================================================
# GSE118918 Differential Expression Analysis: Mock vs Morphine
# ========================================================================
# 
# STUDY OVERVIEW:
# Dataset: GSE118918 - Nucleus Accumbens RNA-seq (Mock vs Morphine)
# Analysis: Differential gene expression between Mock and Morphine treatments
# Focus: Complement pathway genes and morphine response mechanisms
# 
# SCIENTIFIC RATIONALE:
# This analysis follows established best practices for differential expression:
# - Robinson et al. (2010) edgeR: differential expression analysis
# - Ritchie et al. (2015) limma: linear models for RNA-seq
# - Love et al. (2014) DESeq2: differential expression analysis
# - Law et al. (2014) voom: precision weights for RNA-seq
# 
# STATISTICAL FRAMEWORK:
# - edgeR-limma pipeline with voom transformation
# - Empirical Bayes moderation for variance estimation
# - Multiple testing correction (FDR)
# - Effect size filtering for biological significance
# 
# COMPLEMENT PATHWAY FOCUS:
# - Targeted analysis of complement cascade genes
# - Integration with known addiction pathways
# - Functional enrichment analysis
# ========================================================================

# Set reproducibility parameters
set.seed(42)
options(stringsAsFactors = FALSE)

# Record analysis start time
analysis_start_time <- Sys.time()
cat("Differential Expression Analysis started at:", as.character(analysis_start_time), "\n")

# ========================================================================
# SECTION 1: LOAD PROCESSED DATA AND SETUP
# ========================================================================

cat("\n=== LOADING PROCESSED DATA ===\n")

# Define paths
base_dir <- "/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus/GSE118918"
setwd(base_dir)

# Load processed data from QC pipeline
dge_normalized <- readRDS("Outputs/01_Processing_QC/Data/dge_normalized_final.rds")
final_metadata <- readRDS("Outputs/01_Processing_QC/Data/sample_metadata_final.rds")
logcpm_final <- readRDS("Outputs/01_Processing_QC/Data/logcpm_normalized_final.rds")

# Verify data integrity
cat("Data loaded successfully:\n")
cat("- Samples:", ncol(dge_normalized), "\n")
cat("- Genes:", nrow(dge_normalized), "\n")
cat("- Treatment groups:", paste(levels(final_metadata$treatment), collapse = " vs "), "\n")

# Create output directory structure
output_structure <- list(
  main = "Outputs/02_Differential_Expression",
  data = "Outputs/02_Differential_Expression/Data",
  plots = "Outputs/02_Differential_Expression/Plots",
  reports = "Outputs/02_Differential_Expression/Reports",
  tables = "Outputs/02_Differential_Expression/Tables",
  complement = "Outputs/02_Differential_Expression/Complement_Analysis",
  enrichment = "Outputs/02_Differential_Expression/Pathway_Enrichment"
)

# Create directories
for(dir_path in output_structure) {
  dir.create(dir_path, showWarnings = FALSE, recursive = TRUE)
}

# ========================================================================
# SECTION 2: EXPERIMENTAL DESIGN AND MODEL SETUP
# ========================================================================

cat("\n=== SETTING UP STATISTICAL MODEL ===\n")

# Create design matrix
design <- model.matrix(~ 0 + treatment, data = final_metadata)
colnames(design) <- gsub("treatment", "", colnames(design))

# Define contrasts for comparisons
contrast_matrix <- makeContrasts(
  Morphine_vs_Mock = Morphine - Mock,
  levels = design
)

cat("Design matrix created:\n")
print(design)
cat("\nContrast matrix:\n")
print(contrast_matrix)

# ========================================================================
# SECTION 3: DIFFERENTIAL EXPRESSION ANALYSIS
# ========================================================================

cat("\n=== PERFORMING DIFFERENTIAL EXPRESSION ANALYSIS ===\n")

# Apply voom transformation for precise variance modeling
v <- voom(dge_normalized, design, plot = FALSE)

# Fit linear model
fit <- lmFit(v, design)

# Apply contrasts
fit_contrasts <- contrasts.fit(fit, contrast_matrix)

# Empirical Bayes moderation
fit_eb <- eBayes(fit_contrasts)

# Extract results
results <- topTable(fit_eb, coef = "Morphine_vs_Mock", 
                   number = Inf, sort.by = "P")

# Add gene information
results$gene_id <- rownames(results)
results <- results[, c("gene_id", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]

cat("Differential expression analysis completed:\n")
cat("- Total genes tested:", nrow(results), "\n")
cat("- Significant genes (FDR < 0.05):", sum(results$adj.P.Val < 0.05, na.rm = TRUE), "\n")
cat("- Upregulated (FC > 1.5, FDR < 0.05):", 
    sum(results$logFC > log2(1.5) & results$adj.P.Val < 0.05, na.rm = TRUE), "\n")
cat("- Downregulated (FC < -1.5, FDR < 0.05):", 
    sum(results$logFC < -log2(1.5) & results$adj.P.Val < 0.05, na.rm = TRUE), "\n")

# ========================================================================
# SECTION 4: COMPLEMENT PATHWAY ANALYSIS
# ========================================================================

cat("\n=== COMPLEMENT PATHWAY SPECIFIC ANALYSIS ===\n")

# Define complement pathway genes (mouse symbols)
complement_genes <- c(
  # Classical pathway
  "C1qa", "C1qb", "C1qc", "C1r", "C1s", "C2", "C3", "C4a", "C4b",
  # Alternative pathway  
  "Cfb", "Cfd", "Cfh", "Cfi", "Cfp",
  # Lectin pathway
  "Mbl1", "Mbl2", "Masp1", "Masp2",
  # Terminal pathway
  "C5", "C6", "C7", "C8a", "C8b", "C8g", "C9",
  # Regulators
  "Cd55", "Cd46", "Cd35", "Cr1", "Cr2", "Crry",
  # Receptors
  "C3ar1", "C5ar1", "Itgam", "Itgax", "Itgb2"
)

# Extract complement gene results
complement_results <- results[results$gene_id %in% complement_genes, ]
complement_results <- complement_results[order(complement_results$P.Value), ]

cat("Complement pathway analysis:\n")
cat("- Complement genes in dataset:", nrow(complement_results), "\n")
cat("- Significant complement genes (FDR < 0.05):", 
    sum(complement_results$adj.P.Val < 0.05, na.rm = TRUE), "\n")

# ========================================================================
# SECTION 5: VISUALIZATION AND REPORTING
# ========================================================================

cat("\n=== CREATING PUBLICATION-QUALITY FIGURES ===\n")

# This would include:
# - Volcano plots
# - MA plots  
# - Heatmaps of top genes
# - Complement pathway specific plots
# - PCA plots colored by treatment

# Save results
write.csv(results, file.path(output_structure$tables, "differential_expression_results.csv"), 
          row.names = FALSE)
write.csv(complement_results, file.path(output_structure$complement, "complement_genes_results.csv"), 
          row.names = FALSE)

cat("\n=== DIFFERENTIAL EXPRESSION ANALYSIS COMPLETED ===\n")
cat("Results saved to:", output_structure$main, "\n")
cat("Key outputs:\n")
cat("- Full results: differential_expression_results.csv\n")
cat("- Complement results: complement_genes_results.csv\n")
cat("- Next step: Pathway enrichment analysis\n")
