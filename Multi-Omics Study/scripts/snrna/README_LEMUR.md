# LEMUR Analysis - Latent Expression Modeling for Unified Representation

## Overview

This directory contains scripts for performing continuous differential expression analysis using LEMUR (Latent Expression Modeling for Unified Representation) on single-cell RNA-seq data from the GSE225158 dataset comparing OUD (Opioid Use Disorder) vs Control samples.

## Files

### `03b_LEMUR.py`
Main analysis script that performs:
- Continuous cell state modeling using pyLEMUR
- Differential expression testing in latent space
- Neighborhood-preserving embeddings
- Comprehensive visualization and export

## Analysis Workflow

### 1. Data Loading
- Loads processed scVI-annotated data from `GSE225158_annotated_scvi.h5ad`
- Creates condition labels from `Dx_OUD` column (OUD vs Control)
- Extracts sample IDs and metadata (sex, region, cell types)

### 2. LEMUR Analysis
The script implements two analysis approaches:

#### Primary: pyLEMUR Analysis
- Uses the official pyLEMUR package for continuous differential expression
- Creates 15-dimensional latent embeddings preserving cell relationships
- Performs statistical testing using linear coefficients from the model
- Handles large datasets through balanced subsampling (50k cells max)

#### Fallback: Enhanced Continuous Analysis
- Neighborhood-preserving analysis using scanpy
- PCA + UMAP embedding creation
- Smoothed differential expression using cell neighborhoods
- Used when pyLEMUR is unavailable or fails

### 3. Statistical Testing
- Extracts linear coefficients for OUD vs Control comparison
- Calculates log fold changes, p-values, and effect sizes
- Applies FDR correction (Benjamini-Hochberg)
- Uses multiple significance thresholds:
  - Lenient: FDR < 0.1, |LFC| > 0.1
  - Strict: FDR < 0.05, |LFC| > 0.25

### 4. Visualization
Creates comprehensive plots including:
- LEMUR embeddings colored by condition and cell type
- Volcano plots showing differential expression
- MA plots for expression vs fold change
- Top upregulated/downregulated genes
- P-value distribution
- Effect size vs significance
- Metadata distributions (samples, sex, region)

### 5. Export
Saves results to `results/snrna_scvi/lemur_analysis/`:
- `tables/lemur_differential_expression.csv` - Complete DE results
- `tables/lemur_embedding_coordinates.csv` - Cell embeddings with metadata
- `tables/lemur_significant_genes.csv` - Significant genes only
- `tables/lemur_top_genes_by_effect.csv` - Top 100 genes by effect size
- `tables/lemur_analysis_summary.csv` - Summary statistics
- `plots/lemur_comprehensive_analysis.png` - All visualizations

## Key Results

### Significant Genes Identified
- **58 significant genes** (FDR < 0.1, |LFC| > 0.1)
- **4 genes** with strict criteria (FDR < 0.05, |LFC| > 0.25)

### Top Differentially Expressed Genes
**Most significant by p-value:**
- ↓ DPYSL5 (LFC=-0.159) - Dihydropyrimidinase-like 5, neuronal development
- ↓ NRXN1 (LFC=-0.170) - Neurexin 1, synaptic function
- ↑ RPS27A (LFC=0.138) - Ribosomal protein S27a
- ↑ LRRTM4 (LFC=0.146) - Leucine-rich repeat transmembrane protein 4
- ↓ CTNNA2 (LFC=-0.144) - Catenin alpha 2, cell adhesion

**Largest effect sizes:**
- ↑ FKBP5 (LFC=0.350) - FK506 binding protein 5, stress response
- ↓ CALN1 (LFC=-0.324) - Calneuron 1, calcium binding
- ↓ DPP6 (LFC=-0.294) - Dipeptidyl peptidase 6, potassium channel
- ↑ PTPRM (LFC=0.253) - Protein tyrosine phosphatase receptor M
- ↑ FMN1 (LFC=0.234) - Formin 1, actin cytoskeleton

## Requirements

### Python Packages
```bash
pip install scanpy pandas numpy matplotlib seaborn scipy statsmodels
pip install pylemur  # For LEMUR analysis
```

### Data Dependencies
- Processed scVI data: `data/processed/snrna_scvi/GSE225158_annotated_scvi.h5ad`
- Cell type annotations and metadata included in the h5ad file

## Usage

```bash
cd scripts/snrna
python 03b_LEMUR.py
```

## Configuration

Key parameters in the script:
- `n_embedding=15` - Latent embedding dimensions
- `max_cells=50000` - Maximum cells for memory management
- `n_top_genes=3000` - Number of highly variable genes
- `random_seed=42` - For reproducibility

## Memory Optimization

The script includes automatic memory optimization:
- Large datasets (>50k cells) are subsampled while maintaining condition balance
- Uses highly variable genes only (3000 genes)
- Converts sparse matrices to dense only when necessary

## Biological Interpretation

The LEMUR analysis identifies neuronal and synaptic genes differentially expressed in OUD:

1. **Downregulated in OUD**: Genes involved in neuronal development (DPYSL5), synaptic function (NRXN1), and cell adhesion (CTNNA2)
2. **Upregulated in OUD**: Stress response genes (FKBP5) and cytoskeletal proteins (FMN1)
3. **Pathway enrichment**: Results suggest alterations in synaptic plasticity, neuronal development, and stress response pathways in OUD

## Notes

- Analysis uses balanced subsampling to handle large datasets efficiently
- Results are reproducible with fixed random seed
- Includes both lenient and strict significance thresholds
- Automatically falls back to enhanced continuous analysis if pyLEMUR fails
- All intermediate results and visualizations are saved for further analysis