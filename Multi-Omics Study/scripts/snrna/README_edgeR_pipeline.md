# edgeR Analysis Pipeline - GSE225158 OUD Study

## Overview

This directory contains a comprehensive pipeline for differential expression analysis using edgeR followed by pathway enrichment analysis. The pipeline analyzes single-cell RNA-seq data from GSE225158 (OUD vs Control) with a focus on complement system pathways.

## Pipeline Components

### 1. 📊 Differential Expression (`03b_differential_expression_edgeR.py`)
- **Purpose**: Robust DE analysis using edgeR with pseudo-bulk aggregation
- **Method**: Generalized linear models with proper dispersion estimation
- **Contrasts**: 7 different comparisons including subgroup and interaction analyses
- **Output**: CSV files with DE results, plots, and statistical summaries

### 2. 🛤️ Pathway Enrichment (`04_pathway_enrichment_edgeR.py`)
- **Purpose**: Comprehensive pathway enrichment from edgeR results
- **Methods**: Gene Ontology (GO), KEGG, Reactome pathway analysis
- **Features**: Complement system focus, visualization, cross-contrast comparison
- **Output**: Enrichment tables, plots, and biological interpretation

### 3. 🚀 Complete Pipeline (`run_edgeR_pipeline.py`)
- **Purpose**: Automated execution of the full analysis workflow
- **Features**: Prerequisite checking, error handling, comprehensive reporting
- **Output**: Combined results and pipeline execution report

## Quick Start

### Option 1: Run Complete Pipeline (Recommended)
```bash
cd "Multi-Omics Study/scripts/snrna"
python3 run_edgeR_pipeline.py
```

### Option 2: Run Individual Steps
```bash
# Step 1: Differential Expression
python3 03b_differential_expression_edgeR.py

# Step 2: Pathway Enrichment (requires Step 1 results)
python3 04_pathway_enrichment_edgeR.py
```

## Requirements

### Python Dependencies
- `scanpy` - Single-cell analysis
- `pandas` - Data manipulation
- `numpy` - Numerical computing
- `matplotlib` - Plotting
- `scipy` - Statistical functions
- `rpy2` - R integration

### R Dependencies
- `edgeR` - Differential expression analysis
- `clusterProfiler` - Pathway enrichment
- `org.Hs.eg.db` - Human gene annotations
- `enrichplot` - Enrichment visualization
- `ReactomePA` - Reactome pathway analysis
- `DOSE` - Disease ontology analysis

### Data Requirements
- Input file: `GSE225158_annotated_scvi.h5ad`
- Location: `data/processed/snrna_scvi/`
- Format: Scanpy AnnData object with cell type annotations

## Analysis Details

### Contrasts Analyzed
1. **OUD_vs_Control_Putamen** - OUD vs Control in Putamen region
2. **OUD_vs_Control_Caudate** - OUD vs Control in Caudate region
3. **OUD_vs_Control_Male** - OUD vs Control in Male subjects
4. **OUD_vs_Control_Female** - OUD vs Control in Female subjects
5. **OUD_Effect_Male_vs_Female** - Interaction: OUD effect differences by sex
6. **OUD_Effect_Putamen_vs_Caudate** - Interaction: OUD effect differences by region
7. **Pooled_OUD_vs_Control** - Overall OUD vs Control (all samples)

### Statistical Thresholds
- **FDR**: < 0.05 (multiple testing correction)
- **Log Fold Change**: |logFC| > 0.5
- **Minimum genes per pathway**: 5
- **Maximum genes per pathway**: 500

### Output Structure
```
results/snrna_scvi/
├── differential_expression_edgeR/
│   ├── edgeR_*_results.csv           # Full DE results per contrast
│   ├── edgeR_*_significant.csv       # Significant genes only
│   ├── edgeR_summary.csv             # Analysis summary
│   ├── edgeR_all_contrasts.csv       # Combined results
│   ├── pseudobulk_counts.csv         # Aggregated count data
│   ├── pseudobulk_metadata.csv       # Sample metadata
│   └── plots/                        # Visualization files
│
├── pathway_enrichment_edgeR/
│   ├── GO_BP_*_enrichment.csv        # Gene Ontology Biological Process
│   ├── GO_MF_*_enrichment.csv        # Gene Ontology Molecular Function
│   ├── GO_CC_*_enrichment.csv        # Gene Ontology Cellular Component
│   ├── KEGG_*_enrichment.csv         # KEGG pathway enrichment
│   ├── Reactome_*_enrichment.csv     # Reactome pathway enrichment
│   ├── all_pathway_enrichment_results.csv  # Combined enrichment
│   ├── complement_pathway_enrichment.csv   # Complement-focused results
│   ├── pathway_enrichment_summary.txt      # Analysis summary
│   └── plots/                        # Enrichment visualizations
│
└── edgeR_pipeline_report.txt         # Complete pipeline summary
```

## Key Features

### 🧬 Complement System Focus
- Automated identification of complement-related pathways
- Keywords: complement activation, classical/alternative/lectin pathways, MAC
- Separate analysis and visualization of complement results

### 📈 Comprehensive Visualizations
- **Volcano plots** - DE gene distributions
- **MA plots** - Mean-difference relationships
- **Heatmaps** - Pathway enrichment across contrasts
- **Dot plots** - Enrichment significance and gene counts
- **Network plots** - Pathway relationships
- **Comparison plots** - Cross-contrast analysis

### 🔄 Interaction Analysis
- Sex-specific OUD effects (Male vs Female interaction)
- Region-specific OUD effects (Putamen vs Caudate interaction)
- Proper statistical modeling using GLM interaction terms

### 🎯 Quality Control
- Minimum sample size checks
- Gene filtering (low expression removal)
- Dispersion estimation validation
- Multiple testing correction (Benjamini-Hochberg)

## Troubleshooting

### Common Issues

1. **R packages not found**
   ```r
   # Install missing packages
   if (!require("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   BiocManager::install(c("edgeR", "clusterProfiler", "org.Hs.eg.db"))
   ```

2. **rpy2 installation issues**
   ```bash
   # Ensure R is in PATH, then:
   pip install rpy2
   ```

3. **Memory issues with large datasets**
   - Increase available memory
   - Consider subsetting data for testing

4. **No significant pathways found**
   - Check gene ID conversion (symbols to Entrez)
   - Verify significance thresholds
   - Ensure minimum gene count requirements

### Expected Runtime
- **edgeR analysis**: 5-15 minutes (depending on sample size)
- **Pathway enrichment**: 10-30 minutes (depending on number of contrasts)
- **Total pipeline**: 15-45 minutes

## Interpretation Guide

### Differential Expression Results
- **logFC > 0**: Upregulated in OUD vs Control
- **logFC < 0**: Downregulated in OUD vs Control
- **FDR**: Adjusted p-value (controls false discovery rate)
- **Count**: Number of samples contributing to the analysis

### Pathway Enrichment Results
- **pvalue**: Raw enrichment p-value
- **qvalue**: FDR-adjusted p-value
- **Count**: Number of DE genes in pathway
- **geneID**: Contributing genes (gene symbols)

### Biological Interpretation
1. **Review top significant pathways** in each contrast
2. **Focus on complement pathways** for study relevance
3. **Compare patterns across** brain regions and sexes
4. **Validate key findings** with literature and experiments
5. **Consider pathway interactions** and crosstalk

## Citation

If you use this pipeline, please cite:
- **edgeR**: Robinson MD, McCarthy DJ, Smyth GK (2010). Bioinformatics
- **clusterProfiler**: Yu G, Wang LG, Han Y, He QY (2012). OMICS
- **ReactomePA**: Yu G, He QY (2016). Molecular BioSystems

## Support

For questions or issues:
1. Check this README and troubleshooting section
2. Review the console output for specific error messages
3. Verify all dependencies are properly installed
4. Ensure input data format matches requirements

## Version History

- **v1.0**: Initial implementation with comprehensive DE and pathway analysis
- **Features**: edgeR GLM, multi-contrast analysis, complement focus, automated pipeline