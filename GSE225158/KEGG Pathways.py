from bioservices import KEGG
import csv

# Initialize KEGG client
kegg = KEGG()
kegg.organism = "hsa"  # Homo sapiens

# Define KEGG pathways of interest
brain_inflammation_pathways = {
    "complement_cascade": "hsa04610",
    "nfkb_pathway": "hsa04064",
    "toll_like_receptor": "hsa04620",
    "nod_like_receptor": "hsa04621",
    "jak_stat_pathway": "hsa04630",
    "cytokine_signaling": "hsa04060",        # Cytokine–cytokine receptor interaction
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

        for entry in parsed.get("GENE", {}).values():
            parts = entry.split(';')[0].split()
            if len(parts) >= 2:
                gene = parts[1]
                rows.append({"pathway": name, "gene": gene})
                unique_genes.add(gene)
    except Exception as e:
        print(f"Error fetching {name}: {e}")

# Export full pathway-gene mapping
with open("kegg_inflammatory_pathways.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["pathway", "gene"])
    writer.writeheader()
    writer.writerows(rows)

# Export unique gene list
with open("kegg_unique_genes.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["gene"])
    for gene in sorted(unique_genes):
        writer.writerow([gene])

print("\n✅ Exported:")
print(" - kegg_inflammatory_pathways.csv (pathway-gene pairs)")
print(" - kegg_unique_genes.csv (unique gene symbols)")