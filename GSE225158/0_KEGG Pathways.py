from bioservices import KEGG
import csv

# Initialize KEGG client
kegg = KEGG()
kegg.organism = "hsa"  # Homo sapiens

# Define KEGG pathways of interest -> comment out some pathways for debugging
brain_inflammation_pathways = {
    "complement_cascade": "hsa04610",
    "nfkb_pathway": "hsa04064",
    "toll_like_receptor": "hsa04620",
    "nod_like_receptor": "hsa04621",
    "jak_stat_pathway": "hsa04630",
    "cytokine_signaling": "hsa04060",
    "tnf_signaling": "hsa04668",             # Optional
    "chemokine_signaling": "hsa04062",       # Optional
    "il17_pathway": "hsa04657"               # Optional
}

rows = []
unique_genes = set()

# Fetch genes from each pathway
for name, kegg_id in brain_inflammation_pathways.items():
    try:
        print(f"Fetching {name} ({kegg_id})...")
        raw = kegg.get(kegg_id)
        parsed = kegg.parse(raw)

        print(f"First 10 rows of {name}:")
        print(dict(list(parsed.get("GENE", {}).items())[:10]))

        # Extract the gene entries
        if "GENE" in parsed:
            for gene_id, desc in parsed["GENE"].items():
                # desc example: 'AKT3; AKT serine/threonine kinase 3 [KO:K04456] [EC:2.7.11.1]'
                parts = desc.split(';')
                gene = parts[0].strip()
                rows.append({"pathway": name, "gene": gene})
                unique_genes.add(gene)
        else:
            print(f"No GENE data found for {name}")

    except Exception as e:
        print(f"Error fetching {name}: {e}")

# Export full pathway-gene mapping to CSV
with open("KEGG exports/kegg_inflammatory_pathways.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["pathway", "gene"])
    writer.writeheader()
    writer.writerows(rows)

# Export unique gene list to CSV
with open("KEGG exports/kegg_unique_genes.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["gene"])
    for gene in sorted(unique_genes):
        writer.writerow([gene])

print("\n✅ Exported:")
print(" - kegg_inflammatory_pathways.csv (pathway-gene pairs)")
print(" - kegg_unique_genes.csv (unique gene symbols)")