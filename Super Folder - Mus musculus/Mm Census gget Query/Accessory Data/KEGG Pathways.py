from bioservices import KEGG
import csv
import os
import datetime

# Initialize KEGG client
kegg = KEGG()
kegg.organism = "mmu"  # Mus musculus (mouse)

# Define KEGG pathways of interest for mice
mouse_inflammation_pathways = {
    "complement_cascade": "mmu04610",
    # "nfkb_pathway": "mmu04064",
    # "toll_like_receptor": "mmu04620",
    # "nod_like_receptor": "mmu04621",
    # "jak_stat_pathway": "mmu04630",
    # "cytokine_signaling": "mmu04060",
    # "tnf_signaling": "mmu04668",
    # "chemokine_signaling": "mmu04062",
    # "il17_pathway": "mmu04657"
}

rows = []
unique_genes = set()

# Fetch genes from each pathway
for name, kegg_id in mouse_inflammation_pathways.items():
    try:
        print(f"Fetching {name} ({kegg_id})...")
        raw = kegg.get(kegg_id)
        parsed = kegg.parse(raw)

        print(f"First 10 rows of {name}:")
        print(dict(list(parsed.get("GENE", {}).items())[:10]))

        # Extract the gene entries
        if "GENE" in parsed:
            for gene_id, desc in parsed["GENE"].items():
                # desc example: 'Akt3; thymoma viral proto-oncogene 3 [KO:K04456] [EC:2.7.11.1]'
                parts = desc.split(';')
                gene = parts[0].strip()
                rows.append({"pathway": name, "gene": gene})
                unique_genes.add(gene)
        else:
            print(f"No GENE data found for {name}")

    except Exception as e:
        print(f"Error fetching {name}: {e}")

# Create a unique timestamp-based directory
timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
output_dir = f"/Mm Census gget Query/Accessory Data"

# Create the directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Export full pathway-gene mapping to CSV
pathway_file = os.path.join(output_dir, "GSE207128_mouse_inflammatory_pathways_gene_mapping.csv")
with open(pathway_file, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["pathway", "gene"])
    writer.writeheader()
    writer.writerows(rows)

# Export unique gene list to CSV
unique_file = os.path.join(output_dir, "GSE207128_mouse_inflammatory_unique_genes.csv")
with open(unique_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["gene"])
    for gene in sorted(unique_genes):
        writer.writerow([gene])

# Export pathway-specific gene lists (one file per pathway)
for pathway in mouse_inflammation_pathways:
    pathway_genes = [row["gene"] for row in rows if row["pathway"] == pathway]
    pathway_file = os.path.join(output_dir, f"GSE207128_mouse_{pathway}_genes.csv")
    with open(pathway_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["gene"])
        for gene in sorted(pathway_genes):
            writer.writerow([gene])

print(f"\n✅ Results exported to: {output_dir}")
print(f" - GSE207128_mouse_inflammatory_pathways_gene_mapping.csv (all pathway-gene pairs)")
print(f" - GSE207128_mouse_inflammatory_unique_genes.csv (all unique genes across pathways)")
print(" - Individual pathway gene lists (one file per pathway)")