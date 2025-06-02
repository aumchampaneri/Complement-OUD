# Comprehensive Cross-Species OUD Analysis

A 12-week bioinformatics project implementing a complete cross-species meta-analysis comparing mouse opioid use disorder (OUD) models to human OUD patients, with focus on complement system and neuroinflammatory responses.

## Project Overview

### Objective
Identify conserved molecular mechanisms of opioid use disorder across species to inform therapeutic target discovery, with particular emphasis on the complement cascade and neuroinflammatory pathways.

### Datasets

#### Mouse Datasets (3)
- **GSE118918**: Bulk RNA-seq, NAcc, males only, acute morphine (20mg/kg, 4h post-injection)
- **GSE207128**: Single-cell RNA-seq, Amygdala, males only, chronic dependence/withdrawal (15mg/kg x14 doses)  
- **GSE289002**: Bulk RNA-seq, NAc+PFC, both sexes, temporal progression (Control→Intoxication→1d withdrawal→14d withdrawal)

#### Human Datasets (2)
- **GSE174409**: Bulk RNA-seq, DLPFC+NAcc, both sexes, 20 OUD vs 20 controls
- **GSE225158**: Single-cell RNA-seq, Caudate+Putamen, both sexes, OUD vs controls

## Directory Structure

```
Translational Study/
├── Scripts/
│   ├── 00_Setup/                    # Project initialization
│   ├── 01_Data_Processing/          # Data download and harmonization
│   ├── 02_Mouse_Analysis/           # Within-species mouse analysis
│   ├── 03_Human_Analysis/           # Within-species human analysis
│   ├── 04_Cross_Species_Integration/# Cross-species comparison
│   ├── 05_Sex_Analysis/             # Sex-specific analysis
│   ├── 06_Complement_Analysis/      # Complement-focused analysis
│   ├── 07_Visualization/            # Publication figures
│   └── 08_Reporting/                # Final integration and reports
├── Data/
│   ├── Raw/                         # Downloaded raw data
│   ├── Processed/                   # Quality-controlled processed data
│   ├── Integrated/                  # Cross-species integrated data
│   └── Orthologs/                   # Ortholog mapping files
├── Results/
│   ├── Mouse_Analysis/              # Mouse-specific results
│   ├── Human_Analysis/              # Human-specific results
│   ├── Cross_Species/               # Cross-species integration results
│   ├── Sex_Specific/                # Sex-stratified analysis results
│   ├── Complement_Focus/            # Complement pathway analysis
│   └── Final_Integration/           # Final integrated results
└── Figures/
    ├── QC/                          # Quality control plots
    ├── Cross_Species/               # Cross-species comparison plots
    ├── Complement/                  # Complement-specific visualizations
    └── Publication/                 # Publication-ready figures
```

## Implementation Phases

### Phase 1: Setup and Data Processing (Weeks 1-2)
- **Task 1.1**: Create comprehensive directory structure ✅
- **Task 1.2**: Install and load required packages ✅
- **Task 1.3**: Process human datasets ✅  
- **Task 1.4**: Harmonize mouse datasets ✅
- **Task 1.5**: Create ortholog mapping ✅

### Phase 2: Within-Species Analysis (Weeks 3-4)
- Mouse meta-analysis across datasets
- Mouse regional and temporal analysis
- Mouse sex-specific analysis
- Human bulk RNA-seq analysis
- Human single-cell analysis

### Phase 3: Cross-Species Integration (Weeks 5-7)
- Direct cross-species comparison
- Regional conservation analysis  
- Pathway-level integration
- Cell-type integration

### Phase 4: Complement-Focused Analysis (Weeks 6-8)
- Define complement gene sets
- Cross-species complement analysis
- Cell-type specific complement responses

### Phase 5: Sex-Specific Analysis (Weeks 8-9)
- Cross-species sex differences
- Sex-stratified meta-analysis

### Phase 6: Advanced Integration & Network Analysis (Weeks 9-10)
- Co-expression network analysis
- Therapeutic target prioritization

### Phase 7: Visualization & Reporting (Weeks 11-12)
- Publication figures
- Statistical reporting
- Final integration report

## Getting Started

### Prerequisites
- R ≥ 4.0.0
- Minimum 16GB RAM (32GB+ recommended)
- 500GB+ available disk space
- Stable internet connection for data downloads

### Installation

1. **Navigate to the project directory:**
   ```bash
   cd "/Users/aumchampaneri/Complement-OUD/Translational Study"
   ```

2. **Run the master setup script:**
   ```r
   source("Scripts/00_Setup/00_Master_Setup.R")
   ```

3. **Install required packages:**
   ```r
   source("Scripts/00_Setup/01_Package_Setup.R")
   ```

### Phase 1 Execution

Run the Phase 1 scripts in order:

1. **Human data processing:**
   ```r
   source("Scripts/01_Data_Processing/01_Human_Data_Processing.R")
   ```

2. **Mouse data harmonization:**
   ```r
   source("Scripts/01_Data_Processing/02_Mouse_Data_Harmonization.R")
   ```

3. **Ortholog mapping:**
   ```r
   source("Scripts/01_Data_Processing/03_Ortholog_Mapping.R")
   ```

## Key Features

### Statistical Methods
- **Differential Expression**: DESeq2 for bulk, Seurat/edgeR for single-cell
- **Meta-analysis**: Random effects models, rank-based methods
- **Multiple Testing**: FDR correction using Benjamini-Hochberg
- **Cross-species**: Ortholog mapping, effect size correlation
- **Network Analysis**: WGCNA for co-expression, igraph for visualization

### Quality Control Standards
- **Gene Filtering**: >10 reads in >50% of samples for bulk, >3 cells for single-cell
- **Sample QC**: Remove outliers >3 SD from mean on key metrics
- **Batch Effects**: ComBat-seq or similar for technical batch correction
- **Ortholog Confidence**: High-confidence 1:1 orthologs for primary analysis

### Complement Gene Sets
- **Classical pathway**: C1qa, C1qb, C1qc, C1r, C1s, C2, C3, C4a, C4b
- **Alternative pathway**: Cfb, Cfd, Cfh, Cfi, Cfp
- **Lectin pathway**: Masp1, Masp2, Mbl1, Mbl2
- **Terminal pathway**: C5, C6, C7, C8a, C8b, C8g, C9
- **Regulators**: Cd55, Cd46, Cd35, Cr1, Cr2, C4bp
- **Receptors**: C3ar1, C5ar1, Cr1, Cr2

## Expected Outcomes

### Success Metrics
- Identify >100 conserved OUD genes across species (FDR < 0.05)
- Demonstrate complement pathway conservation (pathway p < 0.001)
- Generate >20 high-confidence therapeutic targets
- Complete regional conservation analysis for NAcc and PFC/DLPFC
- Characterize sex-specific patterns in both species

### Deliverables
- Comprehensive cross-species gene expression database
- High-confidence ortholog mapping
- Conserved OUD gene signatures
- Complement pathway conservation analysis
- Sex-specific molecular signatures
- Therapeutic target prioritization
- Publication-ready figures and manuscripts

## Technical Specifications

### Performance Requirements
- **Expected Runtime**: ~200 CPU hours total
- **Peak Memory**: ~64GB for single-cell integration
- **Storage**: ~500GB for all results

### Output Formats
- All results in CSV/Excel format with comprehensive annotations
- High-resolution figures (300 DPI) suitable for publication
- Reproducible R scripts with session info and package versions
- Comprehensive documentation of all methods and parameters

## Critical Implementation Notes

1. **Study Heterogeneity**: Different OUD models require careful interpretation
2. **Sex Differences**: Use male-only for primary cross-species, sex-stratified for within-species
3. **Literature Validation**: Use known addiction genes as positive controls
4. **Complement Focus**: Focus on but don't limit analysis to only complement genes
5. **Reproducibility**: All analyses must be fully documented and repeatable

## Support and Documentation

### Generated Reports
- Quality control reports for all datasets
- Ortholog mapping validation reports
- Statistical summaries for all comparisons
- Method documentation with parameters
- Session information and package versions

### Log Files
- All scripts generate detailed log outputs
- Error messages and warnings are captured
- Progress tracking for long-running analyses

## License

This project is part of academic research on opioid use disorder. Please cite appropriately if using any components.

## Contact

For questions about the analysis pipeline or implementation, please refer to the comprehensive documentation generated with each analysis phase.

---

**Last Updated**: June 2025  
**Analysis Version**: 1.0  
**Phase Status**: Phase 1 Implementation Complete ✅
