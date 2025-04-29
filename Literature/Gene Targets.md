# Complement & Inflammatory Pathway Gene List for snRNA-seq Analysis
_Focus: Opiate Dependence in Human Brain (GSE225158)_

## Classical Complement Pathway
| **Gene** | **Function** | **Rationale** |
|---------|--------------|---------------|
| `C1QA`, `C1QB`, `C1QC` | Subunits of C1q; bind immune complexes to initiate classical pathway | Measure activation potential of classical complement cascade |
| `C1R`, `C1S` | Serine proteases in the C1 complex; cleave C4 and C2 | Essential for classical-pathway signal propagation |
| `C2` | Combines with C4b to form C3 convertase (C4b2a) | Required for C3 activation; pathway strength indicator |
| `C4A`, `C4B` | Cleaved by C1s to C4b (binds surfaces) and C4a (anaphylatoxin) | C4 is central to classical pathway amplification |
| `C3` | Central node of all pathways; cleaved to C3a/C3b | Key hub of complement system; used in multiple pathways |

## Alternative Complement Pathway
| **Gene** | **Function** | **Rationale** |
|---------|--------------|---------------|
| `CFB` | Combines with C3b to form C3 convertase (C3bBb) | Essential for amplification loop in AP |
| `CFD` | Cleaves factor B to activate the AP C3 convertase | Rate-limiting enzyme of AP |
| `CFP` *(Properdin)* | Stabilizes AP C3/C5 convertases | Enhances AP activity and initiation |
| `CFH` | Regulates AP by promoting degradation of C3b | Prevents overactivation; modulates self/non-self |
| `CFI` | Protease that deactivates C3b (and C4b) | Downregulates AP and classical pathways |
| `CR1` *(CD35)* | Receptor/cofactor promoting C3b/C4b inactivation | Inhibits complement activation on cell surfaces |
| `CR2` *(CD21)* | Receptor for C3d; links complement to adaptive immunity | Monitors complement-adaptive crosstalk |

## Complosome (Intracellular Complement)
| **Gene** | **Function** | **Rationale** |
|---------|--------------|---------------|
| `C3`, `C5` | Intracellular complement stores cleaved to C3a/C5a | Drive metabolic signaling (via C3aR/C5aR) inside cells |
| `CTSL` | Lysosomal protease that cleaves C3 intracellularly | Required for noncanonical C3 activation in T cells |
| `C3AR1`, `C5AR1` | Receptors for intracellular C3a/C5a | Mediate autocrine signaling, influence metabolism |
| `C5AR2` *(C5L2)* | Modulates C5aR1 signaling | Balances inflammatory vs. regulatory responses |
| `CD46` *(MCP)* | Cofactor for C3b degradation; also a complement receptor | Central to complosome-driven Th1 signaling and metabolic reprogramming |

## Inflammatory Crosstalk Pathways
### NF-κB Signaling
| **Gene** | **Function** | **Rationale** |
|---------|--------------|---------------|
| `NFKB1`, `RELA` | Canonical NF-κB subunits (p50/p65) | Induced by C5a; drive cytokine expression |
| `NFKB2`, `RELB` | Non-canonical NF-κB subunits | Broader pathway readout |
| `NFKBIA` *(IκBα)* | NF-κB inhibitor; feedback regulator | Indicates pathway activation and repression cycle |

### NLRP3 Inflammasome & IL-1 Axis
| **Gene** | **Function** | **Rationale** |
|---------|--------------|---------------|
| `NLRP3`, `PYCARD`, `CASP1` | Inflammasome components | Activated by MAC/C5a; promote IL-1β release |
| `IL1B`, `IL18` | Pro-inflammatory cytokines | Final outputs of inflammasome pathway |
| `IL1R1`, `IL1R2`, `IL1RN` | IL-1 receptors/antagonists | Measure responsiveness and regulatory balance |

### TNF Pathway
| **Gene** | **Function** | **Rationale** |
|---------|--------------|---------------|
| `TNF`, `TNFRSF1A`, `TNFRSF1B` | TNF and its receptors | Complement signaling can induce TNF |
| `TNFAIP3` *(A20)* | Negative feedback regulator | Terminates NF-κB/TNF signaling |
| `MAPK8`, `MAPK14` | JNK/p38 MAPKs downstream of TNF/IL-1 | Amplify inflammatory gene expression |

## Interferon Response
| **Gene** | **Function** | **Rationale** |
|---------|--------------|---------------|
| `IFNG`, `IFNGR1`, `IFNGR2` | Th1 cytokine and receptors | Induced via CD46 and key in brain inflammation. Modulates immune responses and neuronal inflammation. |
| `STAT1`, `IRF1`, `IRF7` | Transcription factors in IFN signaling | Involved in IFN-induced gene transcription. Complement–IFN cross-regulation suspected. |
| `IFNA2`, `IFNA6`, `IFNA13`, `IFNA5`, `IFNA21`, `IFNW1` | Type I IFNs (alpha subtypes and omega) | Alternative type I interferons to replace missing IFNA1 and IFNB1. Reflect viral/inflammatory response overlaps. |
| `IL6` | Acute-phase cytokine induced by C5a, TNF | Reliable marker of systemic inflammation and immune activation. Important for inflammatory responses in the brain. |
| `IL10` | Regulatory cytokine; induced via CD46 | Marker of anti-inflammatory feedback during Th1 contraction. Modulates immune responses in inflammation and tissue repair. |