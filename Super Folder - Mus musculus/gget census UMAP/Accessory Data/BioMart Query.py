import scanpy as sc
import pandas as pd

# Returns a DataFrame with all Ensembl gene IDs and external gene names in the Mus musculus genome
biomart_query = sc.queries.biomart_annotations("mmusculus", # Query for Mus musculus
                                              ["external_gene_name", "ensembl_gene_id"],
                                               host='www.ensembl.org')

# Convert the query result to a DataFrame
biomart_df = pd.DataFrame(biomart_query)

# Save the DataFrame to a CSV file
biomart_df.to_csv("/Users/aumchampaneri/PycharmProjects/Complement-OUD/gget census UMAP/Accessory Data/biomart_mmusculus.csv", index=False)