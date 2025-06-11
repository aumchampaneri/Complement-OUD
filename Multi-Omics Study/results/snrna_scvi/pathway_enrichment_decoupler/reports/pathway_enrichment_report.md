# Pathway Enrichment Analysis Report
## GSE225158 - OUD vs Control Analysis

### Analysis Overview
- **Total Contrasts Analyzed**: 10
- **Enrichment Methods Used**:
  - Transcription Factor Activity (CollecTRI)
  - Pathway Activity (PROGENy)
  - Hallmark Gene Set Enrichment (MSigDB)

### Transcription Factor Analysis Summary
- **Total Significant TFs**: 626
- **Significant TFs by Contrast**:
  - OUD_vs_Control_Pooled: 3 TFs
  - OUD_vs_Control_Putamen: 1 TFs
  - OUD_vs_Control_Caudate: 10 TFs
  - OUD_vs_Control_Male: 51 TFs
  - OUD_vs_Control_Female: 26 TFs
  - OUD_Effect_Male_vs_Female: 161 TFs
  - OUD_Effect_Male_vs_Female_Putamen: 197 TFs
  - OUD_Effect_Male_vs_Female_Caudate: 170 TFs
  - Control_Effect_Male_vs_Female: 3 TFs
  - Control_Effect_Male_vs_Female_Putamen: 4 TFs

### Pathway Analysis Summary (PROGENy)
- **Total Significant Pathways**: 87
- **Significant Pathways by Contrast**:
  - OUD_vs_Control_Pooled: 7 pathways
  - OUD_vs_Control_Putamen: 9 pathways
  - OUD_vs_Control_Caudate: 7 pathways
  - OUD_vs_Control_Male: 7 pathways
  - OUD_vs_Control_Female: 10 pathways
  - OUD_Effect_Male_vs_Female: 11 pathways
  - OUD_Effect_Male_vs_Female_Putamen: 11 pathways
  - OUD_Effect_Male_vs_Female_Caudate: 13 pathways
  - Control_Effect_Male_vs_Female: 8 pathways
  - Control_Effect_Male_vs_Female_Putamen: 4 pathways

### Hallmark Gene Set Analysis Summary
- **Total Significant Hallmarks**: 239
- **Significant Hallmarks by Contrast**:
  - OUD_vs_Control_Pooled: 15 hallmarks
  - OUD_vs_Control_Putamen: 13 hallmarks
  - OUD_vs_Control_Caudate: 15 hallmarks
  - OUD_vs_Control_Male: 41 hallmarks
  - OUD_vs_Control_Female: 26 hallmarks
  - OUD_Effect_Male_vs_Female: 35 hallmarks
  - OUD_Effect_Male_vs_Female_Putamen: 37 hallmarks
  - OUD_Effect_Male_vs_Female_Caudate: 34 hallmarks
  - Control_Effect_Male_vs_Female: 15 hallmarks
  - Control_Effect_Male_vs_Female_Putamen: 8 hallmarks

### Key Findings
#### Most Activated Across Contrasts
**TF Analysis:**
- KDM5D: activated in 5 contrasts
- REL: activated in 4 contrasts
- MYC: activated in 3 contrasts
- EBF3: activated in 3 contrasts
- PROP1: activated in 3 contrasts

**Pathway Analysis:**
- VEGF: activated in 6 contrasts
- EGFR: activated in 5 contrasts
- Estrogen: activated in 4 contrasts
- MAPK: activated in 4 contrasts
- NFkB: activated in 4 contrasts

**Hallmark Analysis:**
- MYC_TARGETS_V1: activated in 6 contrasts
- E2F_TARGETS: activated in 5 contrasts
- MTORC1_SIGNALING: activated in 5 contrasts
- MYC_TARGETS_V2: activated in 5 contrasts
- OXIDATIVE_PHOSPHORYLATION: activated in 4 contrasts

### Files Generated
- Individual contrast plots in `/plots/`
- Comprehensive heatmaps comparing all contrasts
- Detailed results tables in `/tables/`
- Summary statistics and this report in `/reports/`

### Analysis Parameters
- **Significance Threshold**: p-adj < 0.05
- **Enrichment Method**: Univariate Linear Model (ULM)
- **Input Statistics**: t-statistics from DESeq2
- **Networks Used**: CollecTRI, PROGENy, MSigDB Hallmarks
