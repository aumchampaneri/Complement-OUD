# 🧬 GSEA Analysis Pipeline - Consolidated & Streamlined

This directory contains a **consolidated, production-ready** Gene Set Enrichment Analysis (GSEA) pipeline for analyzing the 49 significant genes identified through LEMUR (Latent Embedding Multivariate Regression) analysis.

## 🎯 **Quick Start**

### **Single Command Execution:**
```bash
cd "Multi-Omics Study/scripts/snrna/LEMUR Analysis/LEMUR GSEA"
Rscript gsea_master_pipeline.R
```

**That's it!** The pipeline handles everything automatically: validation, analysis, visualization, and reporting.

---

## 📁 **Clean Directory Structure**

```
LEMUR GSEA/
├── 📜 gsea_master_pipeline.R    # 🚀 MAIN SCRIPT - Run this!
├── ⚙️  gsea_config.R            # Configuration parameters
├── 🔧 gsea_functions.R          # Utility functions
├── 📖 README.md                 # This documentation
├── 📂 outputs/                  # All analysis results
│   ├── 📊 plots/               # Publication-ready visualizations
│   ├── 📋 tables/              # CSV result files
│   ├── 📄 reports/             # Analysis reports
│   └── 📜 logs/                # Execution logs
└── 📂 backup/                   # Original scripts (archived)
    ├── run_GSEA_GSE225158.R
    ├── run_complete_GSEA_pipeline.R
    ├── create_improved_plots.R
    └── validate_setup.R
```

---

## 🚀 **What the Pipeline Does**

### **🔧 Automated Setup**
- ✅ **Package Management** - Installs missing packages automatically
- ✅ **Data Validation** - Checks input data integrity
- ✅ **System Requirements** - Verifies R version and memory
- ✅ **Directory Creation** - Sets up output folders

### **🧬 Comprehensive Analysis**
- **📊 Over-Representation Analysis (ORA)** - Primary method (works reliably)
- **🔬 Gene Set Enrichment Analysis (GSEA)** - Advanced ranked analysis  
- **🗃️ Multiple Databases**:
  - MSigDB (Hallmark, GO terms, Canonical pathways)
  - KEGG (Metabolic and signaling pathways)
  - Reactome (Curated biological pathways)

### **🎨 Publication-Ready Visualizations**
- **Enhanced Dot Plots** - Gene counts with significance coloring
- **Enhanced Bar Plots** - Fold enrichment visualization
- **Volcano Plots** - Statistical vs biological significance
- **Database Summary** - Cross-platform comparison
- **High Resolution** - 300 DPI PNG format

### **📋 Comprehensive Reporting**
- **Markdown Reports** - Human-readable analysis summaries
- **CSV Tables** - Machine-readable results for further analysis
- **Execution Logs** - Complete audit trail
- **R Workspace** - Full reproducibility

---

## 📊 **Expected Results**

Based on the 49 significant genes from LEMUR analysis, you can expect:

### **🔥 Key Pathway Categories**
1. **🧠 Synaptic Function** (15+ pathways)
   - GRIA4 (AMPA receptor), NRXN3 (Neurexin), EPB41L2
2. **📈 Transcription/Signaling** (12+ pathways)  
   - STAT3 (IL6-JAK-STAT3), LIFR, IGF1R
3. **🔗 Cell Adhesion** (8+ pathways)
   - CDH19, CDH20 (Cadherins), CLDN11
4. **⚡ Ion Channels** (6+ pathways)
   - KCNMB4 (BK channels), PDE1A (cGMP)

### **📈 Analysis Statistics**
- **~400+ enriched pathways** across all databases
- **86% gene ID conversion rate** (43/50 genes)
- **Multiple significance levels** (p < 0.05, FDR < 0.05)
- **Cross-database validation** for robust findings

---

## ⚙️ **Configuration**

### **Key Parameters** (in `gsea_config.R`)
```r
# Statistical thresholds
STATS$pvalue_cutoff = 0.05      # P-value cutoff
STATS$qvalue_cutoff = 0.2       # FDR cutoff for GSEA
STATS$min_gs_size = 10          # Minimum gene set size
STATS$max_gs_size = 500         # Maximum gene set size

# Plot settings
PLOT_PARAMS$max_pathways_plot = 20    # Show top 20 pathways
PLOT_PARAMS$width = 12               # Plot width (inches)
PLOT_PARAMS$height = 8               # Plot height (inches)
PLOT_PARAMS$dpi = 300               # Resolution
```

### **Workflow Control** (enable/disable features)
```r
WORKFLOW$run_gsea = TRUE        # Run GSEA analysis
WORKFLOW$run_ora = TRUE         # Run ORA analysis  
WORKFLOW$run_kegg = TRUE        # Include KEGG pathways
WORKFLOW$run_reactome = TRUE    # Include Reactome pathways
WORKFLOW$create_enhanced_plots = TRUE  # Generate visualizations
```

---

## 🔧 **Customization**

### **Adding New Databases**
```r
# In gsea_config.R, add to PATHWAY_DATABASES:
my_database = list(
  enabled = TRUE,
  name = "My Custom Database",
  description = "Custom pathway collection"
)
```

### **Modifying Plot Types**
```r
# In gsea_config.R, control plot generation:
PLOT_TYPES$enhanced_dotplot = TRUE
PLOT_TYPES$volcano_plot = TRUE  
PLOT_TYPES$network_plot = FALSE  # Disable if not needed
```

### **Adjusting Filtering**
```r
# Modify neuroscience keywords in gsea_config.R:
NEURO_KEYWORDS <- c(
  "synaptic", "neuron", "addiction", "reward",
  "dopamine", "serotonin", "GABA", "glutamate",
  # Add your custom keywords here
)
```

---

## 📦 **Requirements**

### **R Version**
- **Minimum:** R 4.0.0
- **Recommended:** R 4.3.0+

### **Required Packages** (auto-installed)
```r
# Bioconductor packages
clusterProfiler, org.Hs.eg.db, msigdbr, enrichplot,
DOSE, ReactomePA, pathview

# CRAN packages  
ggplot2, dplyr, readr, stringr, RColorBrewer,
viridis, ggrepel, scales
```

### **System Requirements**
- **Memory:** 4+ GB RAM recommended
- **Storage:** ~500 MB for results
- **Internet:** Optional (for some databases)

---

## 🕐 **Runtime Expectations**

- **Total Runtime:** ~5-15 minutes
- **Package Installation:** ~2-5 minutes (first run only)
- **Data Loading:** ~30 seconds  
- **Analysis:** ~2-5 minutes
- **Visualization:** ~2-5 minutes
- **Reporting:** ~1 minute

---

## 📊 **Output Guide**

### **📈 Key Files to Review**
1. **`reports/gsea_analysis_report.md`** - Start here! Comprehensive summary
2. **`plots/database_summary.png`** - Overview of all results
3. **`tables/ora_hallmark_results.csv`** - Top-level pathway findings
4. **`plots/*_enhanced_dotplot.png`** - Best publication figures

### **📁 Complete Output Structure**
```
outputs/
├── plots/
│   ├── database_summary.png              # Cross-database comparison
│   ├── hallmark_enhanced_dotplot.png     # Hallmark pathways
│   ├── go_biological_process_enhanced_dotplot.png
│   ├── kegg_enhanced_dotplot.png         # KEGG pathways
│   ├── reactome_enhanced_dotplot.png     # Reactome pathways
│   └── [15+ additional plots]
├── tables/
│   ├── analysis_summary.csv              # Key metrics
│   ├── ora_hallmark_results.csv          # Hallmark results
│   ├── ora_go_bp_results.csv            # GO Biological Process
│   ├── kegg_results.csv                 # KEGG pathways
│   ├── reactome_results.csv             # Reactome pathways
│   └── input_gene_list.csv              # Input gene information
├── reports/
│   └── gsea_analysis_report.md           # Comprehensive report
└── logs/
    └── gsea_pipeline_[timestamp].log     # Execution log
```

---

## 🆘 **Troubleshooting**

### **Common Issues & Solutions**

#### **"Input file not found"**
```bash
# Ensure LEMUR analysis completed first
ls "../outputs/tables/top_de_genes.csv"
```

#### **Package installation errors**
```r
# Install BiocManager first
install.packages("BiocManager")
BiocManager::install("clusterProfiler")
```

#### **Memory issues**
```r
# Check available memory
memory.limit()
# Increase if needed (Windows)
memory.limit(size = 8000)
```

#### **Internet connectivity issues**
```r
# Disable web-based analyses in config
WORKFLOW$run_enrichr = FALSE
```

### **Getting Help**
1. **Check logs:** `outputs/logs/` directory
2. **Verify input:** Ensure LEMUR analysis completed
3. **Update packages:** Run `update.packages()`
4. **Check R version:** `R.version.string`

---

## 🔬 **Scientific Interpretation**

### **Understanding Results**
- **P-value:** Statistical significance of enrichment
- **FDR/Q-value:** Multiple testing corrected p-value  
- **Fold Enrichment:** Observed vs expected gene ratio
- **Gene Count:** Number of input genes in pathway
- **Gene Ratio:** Proportion of input genes in pathway

### **Key Pathways for OUD Research**
- **IL6-JAK-STAT3:** Neuroinflammation and addiction
- **Synaptic Function:** Direct addiction mechanisms
- **cGMP Signaling:** Second messenger dysregulation
- **Cell Adhesion:** Neural connectivity changes

---

## 📚 **References & Citations**

### **Methods Papers**
- **clusterProfiler:** Yu et al. (2012) OMICS
- **MSigDB:** Liberzon et al. (2011) Bioinformatics
- **Reactome:** Fabregat et al. (2018) Nucleic Acids Research

### **Biological Context**
- **STAT3 in addiction:** Multiple publications on neuroplasticity
- **Synaptic function in OUD:** Extensive literature on glutamate signaling
- **cGMP in reward:** Second messenger system research

---

## 📝 **Version History**

- **v2.0** (Current) - Consolidated pipeline with enhanced plotting
- **v1.5** - Added improved visualization functions
- **v1.0** - Original multi-script implementation

---

## 🤝 **Contributing**

### **Extending the Pipeline**
1. **Add new databases** in `gsea_config.R`
2. **Create custom plot types** in `gsea_functions.R`
3. **Modify analysis parameters** as needed
4. **Test thoroughly** before production use

### **Reporting Issues**
- Provide **full error messages**
- Include **R session info** (`sessionInfo()`)
- Share **input data structure** (first few rows)
- Attach **log files** from `outputs/logs/`

---

**🎉 Ready to discover the pathways underlying opioid use disorder!**

*This consolidated pipeline provides a robust, reproducible framework for pathway analysis that can be easily adapted for other multi-omics studies.*