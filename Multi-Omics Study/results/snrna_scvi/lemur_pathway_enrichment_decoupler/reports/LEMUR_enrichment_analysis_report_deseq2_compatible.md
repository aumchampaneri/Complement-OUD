# LEMUR Pathway Enrichment Analysis Report (DESeq2-Compatible)

## Overview
Analysis of 6 LEMUR contrasts using decoupleR
- **DESeq2-Compatible contrasts**: 4
- **LEMUR-Specific contrasts**: 2

## Summary Statistics

### DESeq2-Compatible Contrasts
| Analysis Type | Total Significant Features |
|---|---|
| Transcription Factors | 53 |
| Pathways (PROGENy) | 23 |
| Hallmark Gene Sets | 25 |

### LEMUR-Specific Contrasts
| Analysis Type | Total Significant Features |
|---|---|
| Transcription Factors | 71 |
| Pathways (PROGENy) | 8 |
| Hallmark Gene Sets | 22 |

## DESeq2-Compatible Results

### OUD vs Control Caudate

**Transcription Factors:** 21 significant

- ETV6: -3.517
- FOXD3: -3.295
- FOXP2: -6.646
- HAND2: 4.541
- HOXC8: -3.371

**Pathways:** 7 significant

- EGFR: 3.831
- Hypoxia: -2.412
- JAK-STAT: 5.579
- MAPK: 4.176
- NFkB: 4.607

**Hallmarks:** 18 significant

- ANDROGEN_RESPONSE: -4.240
- APICAL_JUNCTION: -2.755
- CHOLESTEROL_HOMEOSTASIS: -6.145
- EPITHELIAL_MESENCHYMAL_TRANSITION: -3.182
- FATTY_ACID_METABOLISM: -3.012

### OUD vs Control Putamen

**Transcription Factors:** 12 significant

- ARID4B: 3.680
- CREB1: -3.397
- CREBZF: -6.461
- ESRRG: -4.807
- GABPA: 3.354

**Pathways:** 8 significant

- Androgen: -3.071
- EGFR: -3.269
- JAK-STAT: -2.591
- MAPK: -3.930
- NFkB: -3.559

**Hallmarks:** 1 significant

- ANDROGEN_RESPONSE: -5.072

### OUD Effect Male vs Female Caudate

**Transcription Factors:** 16 significant

- ARID4B: 3.343
- DLX1: -4.409
- ESRRG: -3.728
- FOXP2: 4.952
- GLIS3: -4.173

**Pathways:** 4 significant

- Androgen: -4.856
- Hypoxia: 9.221
- TGFb: 4.602
- p53: 2.813

**Hallmarks:** 1 significant

- HYPOXIA: 4.248

### OUD Effect Male vs Female Putamen

**Transcription Factors:** 4 significant

- ID4: 4.262
- MYT1: 4.605
- OLIG1: 6.610
- ZBTB4: -4.900

**Pathways:** 4 significant

- Hypoxia: -4.737
- JAK-STAT: 2.599
- VEGF: 3.386
- p53: 3.129

**Hallmarks:** 5 significant

- ALLOGRAFT_REJECTION: 3.837
- ANDROGEN_RESPONSE: 3.909
- EPITHELIAL_MESENCHYMAL_TRANSITION: -4.217
- INTERFERON_GAMMA_RESPONSE: 2.912
- NOTCH_SIGNALING: 2.990


## LEMUR-Specific Results

### LEMUR Sex OUD Vulnerability Caudate

**Transcription Factors:** 17 significant

- EBF3: -3.260
- ETV7: 3.539
- FOSB: -4.643
- FOXG1: -3.347
- MECP2: -3.311

**Pathways:** 3 significant

- Hypoxia: -8.439
- JAK-STAT: 4.452
- TGFb: 3.094

**Hallmarks:** 6 significant

- EPITHELIAL_MESENCHYMAL_TRANSITION: -3.057
- GLYCOLYSIS: -4.150
- HYPOXIA: -5.572
- INTERFERON_GAMMA_RESPONSE: 3.267
- REACTIVE_OXYGEN_SPECIES_PATHWAY: -2.777

### LEMUR Sex OUD Vulnerability Putamen

**Transcription Factors:** 54 significant

- ATF1: 3.319
- ATF4: 3.582
- BCLAF1: 3.285
- CREB1: 4.251
- DLX1: 7.194

**Pathways:** 5 significant

- Hypoxia: 5.509
- NFkB: 6.063
- TGFb: 3.187
- TNFa: 5.043
- VEGF: 2.639

**Hallmarks:** 16 significant

- ANDROGEN_RESPONSE: 2.938
- APICAL_JUNCTION: 3.877
- APOPTOSIS: 3.030
- COMPLEMENT: 2.507
- ESTROGEN_RESPONSE_EARLY: 3.153


## Cross-Method Comparison Guide

### Direct DESeq2 Equivalents for Comparison:
- **OUD_vs_Control_Caudate** ↔ **03_OUD_vs_Control_Caudate_results.csv**
- **OUD_vs_Control_Putamen** ↔ **02_OUD_vs_Control_Putamen_results.csv**
- **OUD_Effect_Male_vs_Female_Caudate** ↔ **08_OUD_Effect_Male_vs_Female_Caudate_results.csv**
- **OUD_Effect_Male_vs_Female_Putamen** ↔ **07_OUD_Effect_Male_vs_Female_Putamen_results.csv**

### LEMUR-Unique Insights:
- **LEMUR_Sex_OUD_Vulnerability_Caudate**: Sex-specific OUD vulnerability patterns not captured in DESeq2
- **LEMUR_Sex_OUD_Vulnerability_Putamen**: Sex-specific OUD vulnerability patterns not captured in DESeq2
