# LEMUR Analysis Directory Index

## 📁 Directory Overview
This directory contains all scripts and configuration files for LEMUR (Latent Embedding Multivariate Regression) analysis of the GSE225158 snRNA-seq dataset comparing OUD vs Control samples.

## 📋 File Inventory

### 🎯 Core Analysis Scripts
- **`run_LEMUR_GSE225158.R`** - Main end-to-end LEMUR analysis pipeline
- **`00_LEMUR.R`** - Original LEMUR script (legacy)

### 🛠️ Utility Scripts
- **`lemur_utils.R`** - Helper functions for advanced visualizations and analysis
- **`test_lemur_setup.R`** - Environment validation and dependency check
- **`validate_h5ad.R`** - H5AD file structure and metadata validation

### 🚀 Execution Scripts
- **`run_lemur.sh`** - Shell script for one-command execution with package management

### ⚙️ Configuration
- **`lemur_config.yaml`** - Comprehensive configuration file for analysis parameters

### 📚 Documentation
- **`README_LEMUR.md`** - Complete documentation and user guide
- **`INDEX.md`** - This file - directory navigation guide

## 🚀 Quick Start Guide

### 1. First-Time Setup
```bash
# Navigate to this directory
cd "Multi-Omics Study/scripts/snrna/LEMUR Analysis"

# Test your environment
Rscript test_lemur_setup.R

# Validate your data
Rscript validate_h5ad.R
```

### 2. Run Analysis
```bash
# Option A: Full automated run with package installation
./run_lemur.sh --install-packages

# Option B: Run analysis only (if packages already installed)
./run_lemur.sh

# Option C: Direct R execution
Rscript run_LEMUR_GSE225158.R
```

### 3. Troubleshooting
```bash
# Install packages only
./run_lemur.sh --install-packages --skip-analysis

# Get help
./run_lemur.sh --help
```

## 📊 Expected Workflow

1. **Validation** → `test_lemur_setup.R` & `validate_h5ad.R`
2. **Configuration** → Edit `lemur_config.yaml` if needed
3. **Execution** → `run_lemur.sh` or `run_LEMUR_GSE225158.R`
4. **Results** → Check `../../results/snrna/lemur_analysis/`

## 📁 Input Data
- **H5AD File**: `../../data/raw/snrna/GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5ad`
- **Required Metadata**: `donor_id` (batch), `condition` (OUD/Control)

## 📂 Output Location
Results will be saved to: `../../results/snrna/lemur_analysis/`

## 🔧 Customization
- Modify parameters in `lemur_config.yaml`
- Edit `PARAMS` list in `run_LEMUR_GSE225158.R`
- Use functions from `lemur_utils.R` for custom analyses

## 📞 Support
- Read `README_LEMUR.md` for detailed documentation
- Run validation scripts for troubleshooting
- Check log files in output directory

---
**Last Updated**: 2024  
**Compatible With**: R ≥ 4.0.0, Bioconductor ≥ 3.14