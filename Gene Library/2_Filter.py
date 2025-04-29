import pandas as pd
import scanpy as sc
'''
This script filters a gene list to keep only those genes that are present in an AnnData object.

Inputs:
    adata: AnnData object containing gene expression data
Outputs:
    - .csv file with the following columns:
        - Gene Name
        - Ensembl ID
'''
# Load the Ensembl-mapped gene list
gene_df = pd.read_csv("init_gene_dictionary.csv")

# Load AnnData object
adata = sc.read_h5ad("/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE233279/GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5ad")

# adata.var_names are gene names (e.g., C1QA, C3, etc.)
adata_gene_names = adata.var_names

# Filter gene_df to keep only those gene names in adata
filtered_gene_df = gene_df[gene_df["Gene Name"].isin(adata_gene_names)].copy()

# Save filtered list
filtered_gene_df.to_csv("gene_library.csv", index=False)

# log genes that weren't found
all_input_genes = gene_df["Gene Name"].tolist()
missing_genes = sorted(set(all_input_genes) - set(adata_gene_names))

# Print missing genes
if missing_genes:
    print("Missing genes:")
    for gene in missing_genes:
        print(gene)