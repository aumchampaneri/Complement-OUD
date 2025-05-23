# BioMart CSV to find Ensembl IDs from KEGG Query
# KEGG = /Users/aumchampaneri/PycharmProjects/Complement-OUD/gget census UMAP/Accessory Data/GSE207128_mouse_complement_cascade_genes.csv
# BioMart = /Users/aumchampaneri/PycharmProjects/Complement-OUD/gget census UMAP/Accessory Data/biomart_mmusculus.csv
import pandas as pd
import os
import csv

# Load the KEGG genes CSV
kegg_file = '/Super Folder - Mus musculus/gget census UMAP/Accessory Data/GSE207128_mouse_complement_cascade_genes.csv'
kegg_df = pd.read_csv(kegg_file)

# Load the BioMart CSV
biomart_file = '/Super Folder - Mus musculus/gget census UMAP/Accessory Data/biomart_mmusculus.csv'
biomart_df = pd.read_csv(biomart_file)

# Merge the KEGG and BioMart data on the gene name
merged_df = pd.merge(kegg_df, biomart_df, left_on='gene', right_on='external_gene_name', how='inner')

# Save the merged data to a new CSV
output_file = '/Super Folder - Mus musculus/gget census UMAP/Accessory Data/kegg_biomart_alignment.csv'
merged_df.to_csv(output_file, index=False)

print(f"✅ Merged data saved to: {output_file}")