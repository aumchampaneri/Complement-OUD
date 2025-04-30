import pandas as pd

# Dictionary mapping pathway names to gene lists
pathway_dict = {
    "classical_pathway": ["C1QA", "C1QB", "C1QC", "C1R", "C1S", "C2", "C4A", "C4B"],
    "alternative_pathway": ["CFB", "CFD", "CFP", "CFH", "CFI", "CR1", "CR2"],
    "complosome": ["C3", "C5", "CTSL", "C3AR1", "C5AR1", "C5AR2", "CD46"],
    "nfkb_pathway": ["NFKB1", "RELA", "NFKB2", "RELB", "NFKBIA"],
    "inflammasome": ["NLRP3", "PYCARD", "CASP1", "IL1B", "IL18", "IL1R1", "IL1R2", "IL1RN"],
    "tnf_pathway": ["TNF", "TNFRSF1A", "TNFRSF1B", "TNFAIP3", "MAPK8", "MAPK14"],
    "interferon_pathway": ["IFNG", "IFNGR1", "IFNGR2", "STAT1", "IRF1", "IRF7",
                           "IFNA2", "IFNA6", "IFNA13", "IFNA5", "IFNA21", "IFNW1", "IL6", "IL10"]
}

# Flatten into decoupleR gene_set format
gene_set = pd.DataFrame([
    {"source": pathway, "target": gene, "weight": 1.0}
    for pathway, genes in pathway_dict.items()
    for gene in genes
])

# Preview
print(gene_set.head())

# Optionally save to CSV
gene_set.to_csv("multi_pathway_gene_set.csv", index=False)