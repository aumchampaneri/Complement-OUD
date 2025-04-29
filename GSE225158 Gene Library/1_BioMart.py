import scanpy as sc
'''
This script uses `ScanPy` to query BioMart for Ensembl IDs based on a passed Gene list
- "/Literature/Gene Targets.md" -> Explains what genes and why

Inputs:
    None
Outputs:
    - .csv file with the following columns:
        - Gene Name
        - Ensembl ID
'''
# Query BioMart annotations for Homo sapiens
    # Returns a DataFrame with all Ensembl gene IDs and external gene names in the Homo sapiens genome
biomart_query = sc.queries.biomart_annotations("hsapiens", # Query for Homo sapiens
                                              ["ensembl_gene_id", "external_gene_name"],
                                               host='www.ensembl.org')

complement_genes = {
    "classical_pathway": ["C1QA", "C1QB", "C1QC", "C1R", "C1S", "C2", "C4A", "C4B", "C3"],
    "alternative_pathway": ["CFB", "CFD", "CFP", "CFH", "CFI", "CR1", "CR2"],
    "complosome": ["C3", "C5", "CTSL", "C3AR1", "C5AR1", "C5AR2", "CD46"],
    "nfkb_pathway": ["NFKB1", "RELA", "NFKB2", "RELB", "NFKBIA"],
    "inflammasome": ["NLRP3", "PYCARD", "CASP1", "IL1B", "IL18", "IL1R1", "IL1R2", "IL1RN"],
    "tnf_pathway": ["TNF", "TNFRSF1A", "TNFRSF1B", "TNFAIP3", "MAPK8", "MAPK14"],
    "interferon_pathway": ["IFNG", "IFNGR1", "IFNGR2", "STAT1", "IRF1", "IRF7", "IFNA2", "IFNA6", "IFNA13", "IFNA5", "IFNA5", "IFNA21", "IFNW1", "IL6", "IL10"]
}

# Flatten to a unique list of genes
all_genes = sorted(set(gene for genes in complement_genes.values() for gene in genes))

# Filter BioMart query to match your genes
gene_df = biomart_query[biomart_query["external_gene_name"].isin(all_genes)].copy()

# Export result to CSV
gene_df.rename(columns={"external_gene_name": "Gene Name", "ensembl_gene_id": "Ensembl ID"}, inplace=True)
gene_df.sort_values(by="Gene Name", inplace=True)
gene_df.to_csv("init_gene_dictionary.csv", index=False)


# Report missing genes
found_genes = set(gene_df["Gene Name"])
missing_genes = sorted(set(all_genes) - found_genes)

if missing_genes:
    print(f"{len(missing_genes)} genes were not found in the BioMart query:")
    for gene in missing_genes:
        print(f"  - {gene}")

    # Optional: save missing genes to a file
    with open("missing_genes.txt", "w") as f:
        for gene in missing_genes:
            f.write(gene + "\n")
else:
    print("All genes were successfully found in the BioMart query.")