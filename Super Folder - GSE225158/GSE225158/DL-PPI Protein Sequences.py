# protein_downloader.py
import os
import pickle
import time
from mygene import MyGeneInfo
import requests
from Bio import Entrez
from tqdm import tqdm


def filter_gene_list(gene_list):
    """Filter out likely non-coding RNAs and non-protein entries"""
    filtered_gene_list = []
    for gene in gene_list:
        # Skip genes with period in name (likely non-coding RNAs or genomic contigs)
        if isinstance(gene, str) and '.' in gene:
            continue
        # Skip genes with common non-coding RNA prefixes
        if isinstance(gene, str) and any(gene.startswith(prefix) for prefix in
                                         ['LINC', 'MIR', 'SNORD', 'AC', 'AL', 'AP']):
            continue
        filtered_gene_list.append(gene)

    return filtered_gene_list


def download_protein_sequences(gene_list, cache_dir="protein_cache"):
    """Download protein sequences for a list of genes and cache them"""
    # Create cache directory
    os.makedirs(cache_dir, exist_ok=True)

    # Initialize MyGene API
    mg = MyGeneInfo()

    # Set email for NCBI API
    Entrez.email = "aum.champaneri@outlook.com"
    Entrez.api_key = "bee8c07e73f4e5935bdacaa9ecc653e9dc09"

    # Filter gene list
    print("Filtering gene list...")
    filtered_genes = filter_gene_list(gene_list)
    print(f"Filtered out {len(gene_list) - len(filtered_genes)} likely non-coding genes")
    print(f"Proceeding with {len(filtered_genes)} potential protein-coding genes")

    # Load existing cache if available
    gene_to_uniprot_file = os.path.join(cache_dir, "gene_to_uniprot.pkl")
    protein_sequences_file = os.path.join(cache_dir, "protein_sequences.pkl")

    gene_to_uniprot = {}
    protein_sequences_dict = {}

    if os.path.exists(gene_to_uniprot_file):
        with open(gene_to_uniprot_file, 'rb') as f:
            gene_to_uniprot = pickle.load(f)
        print(f"Loaded {len(gene_to_uniprot)} existing gene-to-UniProt mappings")

    if os.path.exists(protein_sequences_file):
        with open(protein_sequences_file, 'rb') as f:
            protein_sequences_dict = pickle.load(f)
        print(f"Loaded {len(protein_sequences_dict)} existing protein sequences")

    # Find genes not yet mapped
    genes_to_query = [g for g in filtered_genes if g not in gene_to_uniprot]
    print(f"Need to map {len(genes_to_query)} new genes")

    # Query in batches
    batch_size = 100
    for i in range(0, len(genes_to_query), batch_size):
        batch = genes_to_query[i:i + batch_size]
        try:
            # Query gene info
            gene_info = mg.querymany(batch, scopes='symbol', fields='uniprot', species='human', returnall=True)

            # Process results - ensure we extract string IDs only
            for entry in gene_info['out']:
                gene_symbol = entry.get('query', '')
                if 'uniprot' in entry:
                    uniprot_data = entry['uniprot']

                    # Handle different return formats from API
                    if isinstance(uniprot_data, dict) and 'Swiss-Prot' in uniprot_data:
                        # If it's a dict with 'Swiss-Prot' key, use that value
                        uniprot_id = uniprot_data['Swiss-Prot']
                        if isinstance(uniprot_id, list) and uniprot_id:
                            gene_to_uniprot[gene_symbol] = uniprot_id[0]
                        elif isinstance(uniprot_id, str):
                            gene_to_uniprot[gene_symbol] = uniprot_id

                    elif isinstance(uniprot_data, list) and uniprot_data:
                        # If it's a list, use the first value
                        if isinstance(uniprot_data[0], dict) and 'Swiss-Prot' in uniprot_data[0]:
                            uniprot_id = uniprot_data[0]['Swiss-Prot']
                            if isinstance(uniprot_id, str):
                                gene_to_uniprot[gene_symbol] = uniprot_id
                        else:
                            gene_to_uniprot[gene_symbol] = uniprot_data[0]

                    elif isinstance(uniprot_data, str):
                        # If it's a string, use it directly
                        gene_to_uniprot[gene_symbol] = uniprot_data

            print(f"Batch {i // batch_size + 1}: Mapped {len(gene_info['out'])} genes")

            # Save progress after each batch
            with open(gene_to_uniprot_file, 'wb') as f:
                pickle.dump(gene_to_uniprot, f)

        except Exception as e:
            print(f"Error in batch query {i // batch_size + 1}: {e}")

        # Add delay between batches to avoid rate limiting
        time.sleep(1)

    # Clean up non-string values
    for gene, uniprot_id in list(gene_to_uniprot.items()):
        if not isinstance(uniprot_id, str):
            print(f"Removing non-string UniProt ID for gene {gene}: {uniprot_id}")
            gene_to_uniprot.pop(gene)

    # Download sequences for mapped UniProt IDs
    uniprot_ids_to_query = [u for u in gene_to_uniprot.values() if u not in protein_sequences_dict]
    print(f"Need to download {len(uniprot_ids_to_query)} new protein sequences")

    for uniprot_id in tqdm(uniprot_ids_to_query, desc="Downloading protein sequences"):
        try:
            # Get protein sequence from UniProt
            url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
            response = requests.get(url)

            if response.status_code == 200:
                # Parse FASTA
                sequence_lines = response.text.strip().split('\n')
                if len(sequence_lines) > 1:
                    # Extract sequence from FASTA (skip the header line)
                    sequence = ''.join(sequence_lines[1:])
                    protein_sequences_dict[uniprot_id] = sequence

            # Save progress every 10 proteins
            if len(protein_sequences_dict) % 10 == 0:
                with open(protein_sequences_file, 'wb') as f:
                    pickle.dump(protein_sequences_dict, f)

        except Exception as e:
            print(f"Error downloading sequence for {uniprot_id}: {e}")

        # Short delay to avoid rate limiting
        time.sleep(0.5)

    # Final save
    with open(gene_to_uniprot_file, 'wb') as f:
        pickle.dump(gene_to_uniprot, f)

    with open(protein_sequences_file, 'wb') as f:
        pickle.dump(protein_sequences_dict, f)

    print(f"Successfully mapped {len(gene_to_uniprot)} genes to UniProt IDs")
    print(f"Downloaded {len(protein_sequences_dict)} protein sequences")
    print(f"Data saved to {cache_dir}")

    return gene_to_uniprot, protein_sequences_dict


if __name__ == "__main__":
    # Example usage - load all genes from your project
    import pandas as pd

    # Load pathway data
    try:
        pathway_file = "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/KEGG outputs/kegg_inflammatory_pathways.csv"
        pathway_df = pd.read_csv(pathway_file)
        all_pathway_genes = pathway_df.gene.unique().tolist()
    except Exception as e:
        print(f"Error loading pathway data: {e}")
        all_pathway_genes = []

    # Load sex-specific genes
    try:
        m_file = "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/deseq2_results/male_specific_degs.csv"
        f_file = "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/deseq2_results/female_specific_degs.csv"
        male_genes = pd.read_csv(m_file)['gene'].tolist()
        female_genes = pd.read_csv(f_file)['gene'].tolist()
    except Exception as e:
        print(f"Error loading sex-specific genes: {e}")
        male_genes = []
        female_genes = []

    # Combine all genes
    all_genes = list(set(all_pathway_genes + male_genes + female_genes))

    # Download and cache protein sequences
    download_protein_sequences(all_genes)