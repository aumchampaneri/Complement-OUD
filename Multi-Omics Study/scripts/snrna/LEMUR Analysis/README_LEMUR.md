# LEMUR Analysis Pipeline for GSE225158 snRNA-seq Data

## 🧬 Overview

This directory contains a comprehensive R-based pipeline for performing LEMUR (Latent Embedding Multivariate Regression) analysis on single-cell RNA-seq data from the GSE225158 dataset, comparing OUD (Opioid Use Disorder) vs Control samples in striatal tissue.

## 📁 Files Description

### Main Analysis Scripts
- **`run_LEMUR_GSE225158.R`** - Main analysis pipeline script
- **`lemur_utils.R`** - Utility functions for LEMUR analysis
- **`run_lemur.sh`** - Shell script for easy execution

### Input Data
- **Input H5AD**: `/data/raw/snrna/GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5ad`
- **Metadata columns**: 
  - `donor_id` - Batch identifier
  - `condition` - Treatment condition ("Control", "OUD")

## 🚀 Quick Start

### Prerequisites
- R (≥ 4.0.0)
- Required R packages (automatically installed)

### Option 1: Using the Shell Script (Recommended)
```bash
# Navigate to the script directory
cd "Multi-Omics Study/scripts/snrna/LEMUR Analysis"

# Run complete analysis with package installation
./run_lemur.sh --install-packages

# Or run analysis only (if packages already installed)
./run_lemur.sh
```

### Option 2: Running R Script Directly
```bash
# In R console or RStudio
source("run_LEMUR_GSE225158.R")
```

## 📊 Analysis Pipeline Steps

### 1. **Data Loading & QC**
- Load H5AD file using `zellkonverter`
- Compute per-cell QC metrics (n_genes, total_counts, pct_mito)
- Apply quality control filters:
  - Genes: 200-6000 per cell
  - UMI counts: 500-50000 per cell
  - Mitochondrial %: <20%

### 2. **Normalization & Feature Selection**
- Size factor normalization using `scran`
- Log-transformation
- Select 3000 highly variable genes (HVGs)

### 3. **Dimensionality Reduction**
- PCA on HVGs (50 components)
- UMAP for visualization (30 PC dimensions)

### 4. **LEMUR Model Fitting**
- Design formula: `~ donor_id + condition`
- Embedding dimensions: 20
- Batch-aware latent space modeling

### 5. **Batch Correction**
- Harmony integration on LEMUR embedding
- Correction for donor-specific effects

### 6. **Differential Testing**
- LEMUR-native differential testing
- Contrast: OUD vs Control
- Multiple testing correction (FDR)

### 7. **Neighborhood Analysis**
- Identify differentially expressed neighborhoods
- Spatial context of DE effects

### 8. **Results & Visualization**
- Comprehensive result tables
- UMAP plots (condition, batch, QC metrics)
- Volcano plots
- DE neighborhood visualization

## 📂 Output Structure

```
scripts/snrna/LEMUR Analysis/outputs/
├── tables/
│   ├── lemur_de_results.csv          # Complete DE results
│   ├── top_de_genes.csv              # Top significant genes
│   └── de_neighborhoods.csv          # Neighborhood summary
├── plots/
│   ├── umap_overview.png             # Multi-panel UMAP
│   ├── volcano_plot.png              # Volcano plot
│   └── de_neighborhoods_umap.png     # DE neighborhoods
├── reports/
│   ├── analysis_report.txt           # Human-readable summary
│   └── analysis_summary.rds          # Structured summary
└── lemur_workspace.RData             # Complete R workspace
```

## 🔧 Configuration Parameters

Key parameters can be modified in the `PARAMS` list:

```r
PARAMS <- list(
  # Quality control
  min_genes = 200,
  max_genes = 6000,
  max_mito_pct = 20,
  
  # LEMUR settings
  n_hvg = 3000,
  n_embedding = 20,
  design_formula = "~ donor_id + condition",
  
  # Statistics
  fdr_threshold = 0.05,
  effect_size_threshold = 0.1
)
```

## 📋 Required R Packages

### Core Bioconductor
- `SingleCellExperiment` - Single-cell data containers
- `scran` - Normalization and QC
- `scater` - Single-cell analysis utilities
- `lemur` - LEMUR analysis

### Data Processing
- `tidyverse` - Data manipulation
- `zellkonverter` - H5AD file reading
- `HDF5Array` - Efficient data storage

### Visualization
- `ggplot2` - Plotting
- `viridis` - Color palettes
- `patchwork` - Plot composition

### Dimensionality Reduction
- `uwot` - UMAP implementation
- `harmony` - Batch correction

## 🔍 Expected Results

### Typical Output Summary
### Expected Output Summary
- **Input**: ~50,000-100,000 cells, ~20,000-30,000 genes
- **After QC**: ~40,000-80,000 cells, ~15,000-20,000 genes
- **HVGs**: 3,000 most variable genes
- **DE genes**: 500-2,000 significant genes (FDR < 0.05)
- **Runtime**: 30-60 minutes (depending on data size)
- **Output Location**: `outputs/` directory within LEMUR Analysis folder

### Key Result Files
1. **`outputs/tables/lemur_de_results.csv`** - Complete differential expression results
2. **`outputs/tables/top_de_genes.csv`** - Top 100 significant genes
3. **`outputs/plots/umap_overview.png`** - Quality control and batch visualization
4. **`outputs/plots/volcano_plot.png`** - DE significance overview

## 🐛 Troubleshooting

### Common Issues

1. **Package Installation Errors**
   ```bash
   # Navigate to LEMUR Analysis directory first
   cd "Multi-Omics Study/scripts/snrna/LEMUR Analysis"
   # Install Bioconductor packages manually
   ./run_lemur.sh --install-packages --skip-analysis
   ```

2. **Memory Issues**
   ```r
   # Reduce parameters for large datasets
   PARAMS$n_hvg <- 2000  # Reduce HVGs
   PARAMS$n_embedding <- 15  # Reduce embedding dimensions
   ```

3. **H5AD Loading Issues**
   ```r
   # Check file path and permissions
   file.exists(INPUT_H5AD)
   # Ensure zellkonverter is properly installed
   ```

### Performance Optimization
- Use HDF5 backend for large datasets
- Reduce embedding dimensions for faster computation
- Filter more stringently for memory constraints

## 📈 Interpreting Results

### Differential Expression Results
- **log_fc**: Log2 fold change (OUD vs Control)
- **pvalue**: Raw p-value from LEMUR test
- **adj_pvalue**: FDR-adjusted p-value
- **significant**: FDR < 0.05 flag

### Visualization Guide
- **UMAP plots**: Check for batch effects and condition separation
- **Volcano plot**: Overview of DE magnitude and significance
- **DE neighborhoods**: Spatial organization of differential effects

## 🔗 References

1. **LEMUR**: Fischer, D.S., et al. "Modeling intercellular communication in tissues using spatial graphs of cells." *Nature Biotechnology* (2022).

2. **GSE225158**: Seney, M.L., et al. "Transcriptional Alterations in the Striatum of Subjects with Opioid Use Disorder." *Biological Psychiatry* (2023).

3. **Analysis Framework**: Amemiya, H.M., Kundaje, A. & Boyle, A.P. "The ENCODE Blacklist: Identification of Problematic Regions of the Genome." *Scientific Reports* (2019).

## 📞 Support

For questions or issues:
1. Check troubleshooting section above
2. Review log files in the output directory
3. Ensure all dependencies are properly installed
4. Verify input data format and metadata columns

---

**Generated by**: LEMUR Analysis Pipeline v1.0  
**Date**: 2024  
**Compatible with**: R ≥ 4.0.0, Bioconductor ≥ 3.14