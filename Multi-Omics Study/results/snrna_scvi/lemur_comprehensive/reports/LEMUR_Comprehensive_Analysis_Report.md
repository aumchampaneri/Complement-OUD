# 🌊 LEMUR Comprehensive Analysis Report

**Generated**: 2025-06-11 17:57:10  
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

- **Significant Genes**: 43
- **Coefficient Threshold**: 0.000826 (adaptive)
- **FDR Threshold**: 0.05

**Top 10 Significant Genes**:

| Gene | Coefficient | P-value | FDR |
|------|-------------|---------|-----|
| AC012494.1 | -0.0372 | 0.00e+00 | 0.00e+00 |
| PARD3B | 0.0271 | 0.00e+00 | 0.00e+00 |
| MTRNR2L12 | -0.0288 | 0.00e+00 | 0.00e+00 |
| DPP10 | 0.0378 | 0.00e+00 | 0.00e+00 |
| ST6GALNAC3 | 0.0168 | 2.22e-16 | 1.33e-13 |
| LINC00882 | -0.0165 | 8.88e-16 | 4.44e-13 |
| MACF1 | -0.0146 | 1.28e-12 | 5.49e-10 |
| SMYD3 | 0.0140 | 8.98e-12 | 3.37e-09 |
| PLCL1 | 0.0133 | 8.39e-11 | 2.80e-08 |
| RGS7 | -0.0129 | 3.02e-10 | 9.05e-08 |

#### Male Vs Female

- **Significant Genes**: 46
- **Coefficient Threshold**: 0.000956 (adaptive)
- **FDR Threshold**: 0.05

**Top 10 Significant Genes**:

| Gene | Coefficient | P-value | FDR |
|------|-------------|---------|-----|
| AC006059.1 | 0.0224 | 0.00e+00 | 0.00e+00 |
| MTRNR2L12 | 0.0199 | 0.00e+00 | 0.00e+00 |
| AC012494.1 | 0.0340 | 0.00e+00 | 0.00e+00 |
| PLD5 | -0.0212 | 0.00e+00 | 0.00e+00 |
| PDE4B | -0.0191 | 0.00e+00 | 0.00e+00 |
| DPP10 | 0.0290 | 0.00e+00 | 0.00e+00 |
| ST6GALNAC3 | -0.0324 | 0.00e+00 | 0.00e+00 |
| LINC00882 | 0.0182 | 4.44e-16 | 1.67e-13 |
| ANKRD44 | 0.0179 | 1.33e-15 | 4.44e-13 |
| AC073050.1 | -0.0161 | 5.74e-13 | 1.72e-10 |

### Sex-Stratified Analysis


#### F Subset
- **Significant Genes**: 55
- **Coefficient Threshold**: 0.001429

#### M Subset
- **Significant Genes**: 39
- **Coefficient Threshold**: 0.001420

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
