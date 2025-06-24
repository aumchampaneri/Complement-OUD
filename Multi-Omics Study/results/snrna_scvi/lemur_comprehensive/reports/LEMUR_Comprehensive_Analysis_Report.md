# 🌊 LEMUR Comprehensive Analysis Report

**Generated**: 2025-06-23 23:55:17
**Dataset**: GSE225158 - Human Striatal Single-Cell RNA-seq
**Analysis**: OUD vs Control with Sex Effects
**Framework**: LEMUR with Harmony Batch Correction

## 📊 Dataset Overview

- **Total Cells**: 98,848
- **Genes Analyzed**: 3,000 (highly variable)
- **Samples**: 22 brain samples
- **Design Formula**: `~ level1 + Sex`
- **Embedding Dimensions**: 15

## 🧬 Differential Expression Results

### Main Contrasts

#### Oud Vs Control

- **Significant Genes**: 46
- **Coefficient Threshold**: 0.000987 (adaptive)
- **FDR Threshold**: 0.05

**Top 10 Significant Genes**:

| Gene | Coefficient | P-value | FDR |
|------|-------------|---------|-----|
| CAMKMT | -0.0285 | 0.00e+00 | 0.00e+00 |
| ANKRD44 | -0.0230 | 0.00e+00 | 0.00e+00 |
| LINC00882 | 0.0265 | 0.00e+00 | 0.00e+00 |
| ESRRG | -0.0273 | 0.00e+00 | 0.00e+00 |
| DPYSL5 | -0.0220 | 0.00e+00 | 0.00e+00 |
| AC079352.1 | 0.0438 | 0.00e+00 | 0.00e+00 |
| CTNNA2 | 0.0205 | 0.00e+00 | 0.00e+00 |
| GLUL | -0.0195 | 1.11e-15 | 4.16e-13 |
| PDE1A | 0.0194 | 1.78e-15 | 5.92e-13 |
| PLCL1 | -0.0189 | 8.22e-15 | 2.46e-12 |

#### Male Vs Female

- **Significant Genes**: 35
- **Coefficient Threshold**: 0.001149 (adaptive)
- **FDR Threshold**: 0.05

**Top 10 Significant Genes**:

| Gene | Coefficient | P-value | FDR |
|------|-------------|---------|-----|
| AC079352.1 | -0.0466 | 0.00e+00 | 0.00e+00 |
| LMCD1-AS1 | 0.0800 | 0.00e+00 | 0.00e+00 |
| GLUL | 0.0328 | 0.00e+00 | 0.00e+00 |
| PDE1A | 0.0612 | 0.00e+00 | 0.00e+00 |
| FOXP1 | -0.0489 | 0.00e+00 | 0.00e+00 |
| AC012494.1 | -0.0366 | 0.00e+00 | 0.00e+00 |
| LRRC7 | 0.0475 | 0.00e+00 | 0.00e+00 |
| MTRNR2L12 | -0.0271 | 1.05e-13 | 3.96e-11 |
| TMEFF2 | 0.0254 | 3.13e-12 | 1.04e-09 |
| SMYD3 | -0.0233 | 1.46e-10 | 4.37e-08 |

### Sex-Stratified Analysis


#### F Subset
- **Significant Genes**: 62
- **Coefficient Threshold**: 0.001603

#### M Subset
- **Significant Genes**: 44
- **Coefficient Threshold**: 0.001239

## 🔬 Technical Details

### Statistical Methods
- **Differential Expression**: LEMUR (Latent Expression Model for Unified Regression)
- **Batch Correction**: Harmony integration across samples
- **Multiple Testing**: Benjamini-Hochberg FDR correction
- **Significance Criteria**: FDR < 0.05 AND |coefficient| > adaptive threshold

### Adaptive Thresholds
The analysis uses data-adaptive significance thresholds based on the 75th percentile
of effect sizes, rather than arbitrary fixed cutoffs. This approach is more appropriate for
single-cell data where effect sizes tend to be smaller.

### Memory Usage
- **RAM Available**: 64 GB
- **Estimated Usage**: ~2.2 GB
- **Efficiency**: Full dataset analysis uses less memory than expected

## 📋 Files Generated

### Tables
- `*_results.csv`: Full differential expression results
- `*_significant.csv`: Significant genes only

### Plots
- `lemur_umap_overview.png`: UMAP visualization of cell embeddings
- `volcano_*.png`: Volcano plots for each contrast
- `top_genes_heatmap.png`: Heatmap of top differentially expressed genes
- `coefficient_distributions.png`: Distribution of effect sizes

### Reports
- `LEMUR_Comprehensive_Analysis_Report.md`: This comprehensive report

## 🎯 Key Findings

1. **Statistical Power**: Full dataset (98K cells) provides excellent statistical power
2. **Effect Detection**: Adaptive thresholds reveal biologically meaningful differences
3. **Batch Correction**: Harmony successfully integrates across 22 samples
4. **Sex Differences**: 2 sex-specific analyses completed

## 🔬 Next Steps

1. **Pathway Enrichment**: Analyze significant genes for biological pathways
2. **Cell Type Specificity**: Determine which cell types drive observed effects
3. **Validation**: Compare results with bulk RNA-seq or other datasets
4. **Functional Analysis**: Literature review and experimental validation

---
*Analysis completed using LEMUR v2.0 - Comprehensive Pipeline*
*For questions or support, refer to the analysis documentation*
