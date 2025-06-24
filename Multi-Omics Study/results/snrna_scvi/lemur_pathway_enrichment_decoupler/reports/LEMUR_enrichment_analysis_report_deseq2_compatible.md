# LEMUR Comprehensive Pathway Enrichment Analysis Report

## Overview
Analysis of 5 LEMUR contrasts using decoupleR
- **Main effects**: 2
- **Sex-stratified**: 2
- **Interaction effects**: 1

## Summary Statistics

### Main Effects
| Analysis Type | Total Significant Features |
|---|---|
| Transcription Factors | 36 |
| Pathways (PROGENy) | 2 |
| Hallmark Gene Sets | 20 |

### Sex-Stratified Effects
| Analysis Type | Total Significant Features |
|---|---|
| Transcription Factors | 36 |
| Pathways (PROGENy) | 9 |
| Hallmark Gene Sets | 20 |

### Interaction Effects
| Analysis Type | Total Significant Features |
|---|---|
| Transcription Factors | 18 |
| Pathways (PROGENy) | 6 |
| Hallmark Gene Sets | 10 |

## Main Effects Results

### Oud Vs Control

**Transcription Factors:** 18 significant

- ATF6: -5.273
- CTNNB1: 1.585
- FOXO1: -0.123
- GATA3: 2.224
- GLI1: 0.077

**Hallmarks:** 10 significant

- ADIPOGENESIS: -1.040
- APICAL_JUNCTION: -1.048
- COMPLEMENT: 2.394
- ESTROGEN_RESPONSE_EARLY: -1.095
- FATTY_ACID_METABOLISM: -1.051

### Male Vs Female

**Transcription Factors:** 18 significant

- AP1: 0.299
- ATF6: 0.163
- ESR1: 0.558
- FOXO3: 0.178
- GATA3: 0.254

**Pathways:** 2 significant

- MAPK: 1.366
- PI3K: 1.329

**Hallmarks:** 10 significant

- ANDROGEN_RESPONSE: -1.342
- APICAL_JUNCTION: 1.211
- COMPLEMENT: 1.299
- ESTROGEN_RESPONSE_LATE: -1.266
- FATTY_ACID_METABOLISM: -1.311


## Sex-Stratified Results

### Oud Vs Control Female

**Transcription Factors:** 18 significant

- AR: 1.078
- ATF6: -1.213
- ESR1: -1.630
- GATA1: 0.866
- GATA3: -1.371

**Pathways:** 2 significant

- Hypoxia: 1.349
- NFkB: 1.236

**Hallmarks:** 10 significant

- APICAL_JUNCTION: -1.242
- APICAL_SURFACE: 1.222
- COMPLEMENT: -1.206
- EPITHELIAL_MESENCHYMAL_TRANSITION: 1.251
- ESTROGEN_RESPONSE_EARLY: 1.294

### Oud Vs Control Male

**Transcription Factors:** 18 significant

- AR: -0.235
- ATF1: 0.379
- ATF3: 0.580
- DLX4: -0.374
- DNMT1: -0.200

**Pathways:** 7 significant

- EGFR: -1.234
- Estrogen: -1.446
- MAPK: -1.260
- NFkB: -1.269
- PI3K: -1.292

**Hallmarks:** 10 significant

- ALLOGRAFT_REJECTION: -1.351
- ANDROGEN_RESPONSE: -1.189
- BILE_ACID_METABOLISM: -1.379
- INTERFERON_GAMMA_RESPONSE: -1.332
- KRAS_SIGNALING_DN: -1.411


## Interaction Effects Results

### Sex Oud Interaction

**Transcription Factors:** 18 significant

- AR: -2.883
- ESR1: 1.127
- GATA3: 1.142
- GLI1: 0.867
- ID2: -7.742

**Pathways:** 6 significant

- Androgen: -1.375
- EGFR: -1.258
- Estrogen: -1.290
- Hypoxia: -1.285
- NFkB: -1.253

**Hallmarks:** 10 significant

- ANDROGEN_RESPONSE: -1.397
- APICAL_JUNCTION: 1.415
- APICAL_SURFACE: -1.336
- ESTROGEN_RESPONSE_LATE: -1.347
- KRAS_SIGNALING_DN: -1.473


## Latent Space Analysis Integration

### LEMUR Embedding Analysis
- **Variance Explained:** 0.906 of total variance
- **Significant Correlations:** 241
- **Coordinate Systems:** scVI UMAP, LEMUR UMAP, LEMUR Components

### Cross-Analysis Interpretation
The pathway enrichment results can be interpreted alongside latent space analysis:
- **Variance components** show which biological axes LEMUR captures
- **Embedding correlations** reveal which covariates drive latent structure
- **Pathway activities** explain the functional meaning of latent components
- **TF activities** identify regulatory drivers of observed patterns

### Recommended Integration Steps
1. Compare pathway activities with LEMUR component correlations
2. Identify TFs that correlate with high-variance components
3. Validate pathway patterns across different coordinate systems
4. Cross-reference with latent space metadata overlays


## Cross-Method Comparison Guide

### Direct DESeq2 Equivalents for Comparison:

### LEMUR-Unique Insights:
- **Oud_Vs_Control_Female**: Sex-specific OUD vulnerability patterns not captured in standard DESeq2
- **Oud_Vs_Control_Male**: Sex-specific OUD vulnerability patterns not captured in standard DESeq2
- **Sex_Oud_Interaction**: Sex-specific OUD vulnerability patterns not captured in standard DESeq2
