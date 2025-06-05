# Comprehensive edgeR Analysis Pipeline - GSE225158 OUD Study

## Overview

This directory contains a comprehensive pipeline for differential expression analysis using edgeR followed by **advanced pathway enrichment analysis**. The pipeline analyzes single-cell RNA-seq data from GSE225158 (OUD vs Control) with specialized focus on complement system and neuroinflammation pathways.

## Pipeline Components

### 1. 📊 Differential Expression (`03b_differential_expression_edgeR.py`)
- **Purpose**: Robust DE analysis using edgeR with pseudo-bulk aggregation
- **Method**: Generalized linear models with proper dispersion estimation
- **Contrasts**: 7 different comparisons including subgroup and interaction analyses
- **Output**: CSV files with DE results, plots, and statistical summaries

### 2. 🛤️ **Advanced Pathway Enrichment** (`04_pathway_enrichment_edgeR.py`)
- **Purpose**: Comprehensive multi-level pathway enrichment analysis
- **Basic Methods**: Gene Ontology (GO), KEGG, Reactome pathway analysis
- **Advanced Methods**: GSEA, MSigDB enrichment, leading edge analysis
- **Specialized Analysis**: Complement system, neuroinflammation pathways
- **Comparative Analysis**: Pathway crosstalk, functional modules, meta-analysis
- **Visualizations**: Heatmaps, volcano plots, network maps, UpSet plots
- **Output**: Organized multi-level results with comprehensive biological insights

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
- **Pathway Enrichment**: p < 0.05 (nominal) for pathway inclusion
- **GSEA**: NES (Normalized Enrichment Score) with p < 0.25
- **Meta-Analysis**: Fisher's method for combining p-values
- **Minimum genes per pathway**: 5 (ORA), 15 (GSEA)
- **Maximum genes per pathway**: 500

### Analysis Methods
- **Over-Representation Analysis (ORA)**: Significant genes vs background
- **Gene Set Enrichment Analysis (GSEA)**: Ranked gene list analysis
- **Leading Edge Analysis**: Core genes driving enrichments
- **Pathway Crosstalk**: Jaccard similarity between pathway gene sets
- **Functional Modules**: Hierarchical clustering of pathways

### Comprehensive Output Structure
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
├── pathway_enrichment_edgeR/         # 🆕 ORGANIZED HIERARCHICAL STRUCTURE
│   ├── 01_basic_enrichment/          # Basic pathway enrichment
│   │   ├── GO_BP_*_enrichment.csv    # Gene Ontology Biological Process
│   │   ├── GO_MF_*_enrichment.csv    # Gene Ontology Molecular Function
│   │   ├── GO_CC_*_enrichment.csv    # Gene Ontology Cellular Component
│   │   ├── KEGG_*_enrichment.csv     # KEGG pathway enrichment
│   │   ├── Reactome_*_enrichment.csv # Reactome pathway enrichment
│   │   └── all_basic_enrichment_results.csv
│   │
│   ├── 02_advanced_analysis/         # 🆕 Advanced enrichment methods
│   │   ├── gsea_analysis/            # Gene Set Enrichment Analysis
│   │   │   ├── GSEA_GO_BP_*.csv      # GSEA GO Biological Process
│   │   │   └── GSEA_KEGG_*.csv       # GSEA KEGG pathways
│   │   ├── msigdb_enrichment/        # MSigDB pathway collections
│   │   │   ├── MSigDB_Hallmark_*.csv # Hallmark gene sets
│   │   │   └── MSigDB_C2_Canonical_*.csv # Canonical pathways
│   │   └── leading_edge/             # Leading edge gene analysis
│   │       └── leading_edge_*.csv    # Core genes driving enrichment
│   │
│   ├── 03_specialized_analysis/      # 🆕 Disease-specific analyses
│   │   ├── complement_system/        # Complement pathway focus
│   │   │   └── complement_pathway_enrichment.csv
│   │   └── neuroinflammation/        # Neuroinflammation pathways
│   │       └── neuroinflammation_pathways.csv
│   │
│   ├── 04_visualizations/            # 🆕 Comprehensive plots
│   │   ├── plots/                    # Basic visualization plots
│   │   ├── network_maps/             # Pathway network visualizations
│   │   ├── upset_plots/              # Pathway overlap analysis
│   │   └── volcano_plots_*.png       # Pathway significance plots
│   │
│   ├── 05_comparative_analysis/      # 🆕 Cross-contrast comparisons
│   │   ├── pathway_crosstalk/        # Pathway interaction analysis
│   │   ├── functional_modules/       # Clustered pathway modules
│   │   ├── pathway_occurrence_matrix.csv # Pattern analysis
│   │   └── meta_pathway_analysis.csv # Meta-analysis results
│   │
│   └── 06_summary_reports/           # Analysis summaries
│       └── pathway_enrichment_summary.txt
│
└── edgeR_pipeline_report.txt         # Complete pipeline summary
```

## Key Features

### 🔬 **Multi-Level Enrichment Analysis**
- **Basic Enrichment**: GO, KEGG, Reactome pathway analysis
- **Advanced Methods**: GSEA with ranked gene lists, MSigDB collections
- **Leading Edge Analysis**: Core genes driving pathway enrichments
- **Meta-Analysis**: Combined evidence across contrasts

### 🧬 **Specialized Disease Focus**
- **Complement System**: Automated identification of complement pathways
  - Keywords: complement activation, classical/alternative/lectin pathways, MAC
- **Neuroinflammation**: Addiction and neuroinflammation pathway focus
  - Keywords: microglia, cytokines, dopamine, opioid, synaptic transmission

### 📊 **Advanced Comparative Analysis**
- **Pathway Crosstalk**: Interaction analysis between pathways
- **Functional Modules**: Clustering of related pathways
- **Pattern Analysis**: Pathway occurrence patterns across contrasts
- **UpSet Plots**: Sophisticated overlap visualization

### 📈 **Comprehensive Visualizations**
- **Basic Plots**: Heatmaps, dot plots, volcano plots
- **Network Visualizations**: Enrichment maps, pathway networks
- **Comparative Plots**: Cross-contrast analysis, overlap matrices
- **Advanced Graphics**: UpSet plots, module visualizations

### 🔄 **Robust Statistical Framework**
- **Multiple Methods**: Over-representation + GSEA approaches
- **Proper Corrections**: FDR adjustment across all analyses
- **Quality Controls**: Minimum gene set sizes, filtering
- **Meta-Analysis**: Fisher's method for combining evidence

### 🎯 **Interaction Analysis**
- Sex-specific OUD effects (Male vs Female interaction)
- Region-specific OUD effects (Putamen vs Caudate interaction)
- Proper statistical modeling using GLM interaction terms

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
- **Basic pathway enrichment**: 10-20 minutes
- **Advanced analyses (GSEA, MSigDB)**: 15-30 minutes
- **Specialized & comparative analyses**: 5-10 minutes
- **Total comprehensive pipeline**: 35-75 minutes

## Interpretation Guide

### Differential Expression Results
- **logFC > 0**: Upregulated in OUD vs Control
- **logFC < 0**: Downregulated in OUD vs Control
- **FDR**: Adjusted p-value (controls false discovery rate)
- **Count**: Number of samples contributing to the analysis

### Basic Pathway Enrichment Results
- **pvalue**: Raw enrichment p-value
- **qvalue**: FDR-adjusted p-value
- **Count**: Number of DE genes in pathway
- **geneID**: Contributing genes (gene symbols)

### GSEA Results
- **NES**: Normalized Enrichment Score (positive = upregulated, negative = downregulated)
- **core_enrichment**: Leading edge genes driving the enrichment
- **setSize**: Total genes in the pathway

### Advanced Analysis Results
- **Leading Edge**: Core genes consistently driving pathway enrichments
- **Crosstalk Matrix**: Similarity scores between pathway pairs
- **Functional Modules**: Pathway clusters based on biological similarity
- **Meta-Analysis**: Combined statistical evidence across contrasts

### Biological Interpretation Workflow
1. **Start with basic enrichment** - identify key biological processes
2. **Examine GSEA results** - understand directional changes
3. **Focus on specialized analyses** - complement system, neuroinflammation
4. **Explore comparative results** - pathway patterns, crosstalk, modules
5. **Validate leading edge genes** - prioritize for experimental follow-up
6. **Consider meta-analysis** - pathways with strongest overall evidence

## Citation

If you use this pipeline, please cite:
- **edgeR**: Robinson MD, McCarthy DJ, Smyth GK (2010). Bioinformatics
- **clusterProfiler**: Yu G, Wang LG, Han Y, He QY (2012). OMICS
- **ReactomePA**: Yu G, He QY (2016). Molecular BioSystems
- **MSigDB**: Liberzon A, et al. (2011). Bioinformatics
- **GSEA**: Subramanian A, et al. (2005). PNAS

## Support

For questions or issues:
1. Check this README and troubleshooting section
2. Review the console output for specific error messages
3. Verify all dependencies are properly installed
4. Ensure input data format matches requirements
5. Check organized output directories for expected results

## Version History

- **v1.0**: Initial implementation with basic DE and pathway analysis
- **v2.0**: 🆕 **Comprehensive Advanced Analysis** - GSEA, MSigDB, specialized analyses
- **Features**: 
  - Multi-level enrichment analysis (ORA + GSEA)
  - Specialized disease focus (complement, neuroinflammation)
  - Advanced comparative analysis (crosstalk, modules, meta-analysis)
  - Organized hierarchical output structure
  - Comprehensive visualization suite