# LEMUR Analysis Outputs Directory

## 📁 Directory Purpose
This directory contains all outputs generated by the LEMUR analysis pipeline for the GSE225158 snRNA-seq dataset (OUD vs Control comparison).

## 📂 Directory Structure

```
outputs/
├── tables/           # Statistical results and data tables
├── plots/            # Visualizations and figures
├── reports/          # Summary reports and logs
└── README.md         # This file
```

## 📊 Expected Output Files

### Tables Directory (`tables/`)
- **`lemur_de_results.csv`** - Complete differential expression results
  - Columns: gene, log_fc, pvalue, adj_pvalue, significant
- **`top_de_genes.csv`** - Top 100 significant DE genes
- **`de_neighborhoods.csv`** - Summary of differential neighborhoods
- **`qc_metrics.csv`** - Quality control metrics summary

### Plots Directory (`plots/`)
- **`umap_overview.png`** - Multi-panel UMAP visualization
  - Condition, batch, UMI counts, gene counts
- **`volcano_plot.png`** - Volcano plot of DE results
- **`de_neighborhoods_umap.png`** - DE neighborhoods on UMAP (if found)
- **`qc_metrics.png`** - Quality control metric distributions
- **`lemur_diagnostics.png`** - LEMUR model diagnostic plots

### Reports Directory (`reports/`)
- **`analysis_report.txt`** - Human-readable analysis summary
- **`analysis_summary.rds`** - Structured R summary object
- **`session_info.txt`** - R session information
- **`lemur_analysis.log`** - Detailed analysis log

## 📈 Additional Files
- **`lemur_workspace.RData`** - Complete R workspace for reproducibility
- **`lemur_fit.rds`** - LEMUR model object
- **`lemur_embedding.csv`** - LEMUR latent embedding
- **`harmony_embedding.csv`** - Harmony-corrected embedding

## 🔍 Understanding the Results

### Differential Expression Results
- **log_fc**: Log2 fold change (OUD vs Control)
- **pvalue**: Raw p-value from LEMUR test
- **adj_pvalue**: FDR-adjusted p-value (Benjamini-Hochberg)
- **significant**: Boolean flag for FDR < 0.05

### Visualizations
- **UMAP plots**: Check for batch effects and biological signal
- **Volcano plot**: Overview of effect sizes and significance
- **DE neighborhoods**: Spatial organization of differential effects

### Quality Metrics
- Cells retained after filtering
- Genes used in analysis
- Batch correction effectiveness
- Model convergence statistics

## 📋 File Naming Convention
- CSV files: Comma-separated values, Excel-compatible
- PNG files: High-resolution figures (300 DPI)
- RDS files: R binary objects for reproducibility
- TXT files: Human-readable reports

## 🗂️ File Sizes (Typical)
- Tables: 1-50 MB (depending on gene count)
- Plots: 1-10 MB per figure
- Workspace: 100-500 MB (complete analysis)
- Total: ~200-1000 MB

## 📅 Timestamp Convention
Files are generated with timestamps in the analysis log. Check `reports/analysis_report.txt` for run details.

## 🔒 Data Integrity
- All results are version-controlled through R session info
- Random seed (42) ensures reproducibility
- Input data checksums logged in reports

## 📞 Troubleshooting
If output files are missing or incomplete:
1. Check `reports/lemur_analysis.log` for errors
2. Verify input data integrity
3. Check disk space availability
4. Review R session warnings in reports

---
**Generated by**: LEMUR Analysis Pipeline  
**Dataset**: GSE225158 OUD vs Control snRNA-seq  
**Analysis Date**: Automatically timestamped in reports  
