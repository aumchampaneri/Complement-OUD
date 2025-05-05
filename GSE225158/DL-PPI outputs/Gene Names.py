import pandas as pd
import requests
import time
import os

# Load data
male_df = pd.read_csv(
    "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/DL-PPI outputs/dl_ppi_viz_data_male_0.4.csv",
    names=["protein1", "protein2", "prediction", "method"])
female_df = pd.read_csv(
    "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/DL-PPI outputs/dl_ppi_viz_data_female_0.4.csv",
    names=["protein1", "protein2", "prediction", "method"])

# Get all unique accession numbers
all_accessions = set(pd.concat([
    male_df["protein1"], male_df["protein2"],
    female_df["protein1"], female_df["protein2"]
]).unique())

print(f"Found {len(all_accessions)} unique UniProt IDs")


# Function to convert UniProt accession to gene name
def get_gene_name(uniprot_id):
    try:
        url = f"https://www.uniprot.org/uniprot/{uniprot_id}.xml"
        response = requests.get(url)
        if response.status_code == 200:
            # Extract gene name from XML - looking for name in first gene tag
            xml = response.text
            gene_start = xml.find("<gene")
            if gene_start > 0:
                name_start = xml.find("<name", gene_start)
                if name_start > 0:
                    name_content_start = xml.find(">", name_start) + 1
                    name_content_end = xml.find("</name>", name_content_start)
                    if name_content_start > 0 and name_content_end > 0:
                        return xml[name_content_start:name_content_end]
        return uniprot_id  # Return original ID if extraction fails
    except:
        return uniprot_id


# Create mapping dictionary
print("Converting UniProt IDs to gene names...")
acc_to_gene = {}
for i, acc in enumerate(all_accessions):
    if i % 10 == 0:  # Progress update
        print(f"Processed {i}/{len(all_accessions)} proteins")

    gene_name = get_gene_name(acc)
    acc_to_gene[acc] = gene_name
    time.sleep(0.5)  # Avoid overloading the server

# Save mapping to a file for future reference
mapping_df = pd.DataFrame(list(acc_to_gene.items()), columns=["UniProt_ID", "Gene_Name"])
mapping_df.to_csv(
    "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/DL-PPI outputs/uniprot_gene_mapping.csv",
    index=False)
print("UniProt to gene name mapping saved.")

# Create new dataframes with gene names instead of accession numbers
male_gene_df = male_df.copy()
female_gene_df = female_df.copy()

# Replace UniProt IDs with gene names
male_gene_df["protein1"] = male_gene_df["protein1"].map(acc_to_gene)
male_gene_df["protein2"] = male_gene_df["protein2"].map(acc_to_gene)
female_gene_df["protein1"] = female_gene_df["protein1"].map(acc_to_gene)
female_gene_df["protein2"] = female_gene_df["protein2"].map(acc_to_gene)

# Save new CSV files with gene names
male_gene_df.to_csv(
    "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/DL-PPI outputs/dl_ppi_viz_data_male_0.4_genes.csv",
    index=False)
female_gene_df.to_csv(
    "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/DL-PPI outputs/dl_ppi_viz_data_female_0.4_genes.csv",
    index=False)

print("Created new CSV files with gene names:")
print("- dl_ppi_viz_data_male_0.4_genes.csv")
print("- dl_ppi_viz_data_female_0.4_genes.csv")