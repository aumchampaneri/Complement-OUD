import os
import requests
import pandas as pd
import time
import numpy as np
from Bio import SeqIO
from io import StringIO
import mygene
import json
import sys
import random
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
from Bio.Align import PairwiseAligner


def map_genes_to_uniprot_via_mygene(gene_ids, batch_size=50):
    """Map gene IDs to UniProt using mygene.info instead of UniProt API"""
    print(f"Mapping {len(gene_ids)} genes to UniProt via mygene.info...")
    mg = mygene.MyGeneInfo()
    mappings = []

    # Convert to list if needed
    gene_ids_list = list(gene_ids) if isinstance(gene_ids, set) else gene_ids

    # Process in batches
    for i in range(0, len(gene_ids_list), batch_size):
        batch = gene_ids_list[i:i+batch_size]
        print(f"Processing batch {i//batch_size + 1}/{(len(gene_ids_list)-1)//batch_size + 1}")

        try:
            # Request both Entrez ID and UniProt ID
            results = mg.querymany(batch, scopes='entrezgene',
                                  fields=['uniprot', 'symbol'],
                                  species='human')

            for result in results:
                entrez_id = result.get('query')
                if not result.get('notfound', False) and 'uniprot' in result:
                    # Try Swiss-Prot first (reviewed entries)
                    uniprot_id = None
                    if 'Swiss-Prot' in result['uniprot']:
                        uniprot_id = result['uniprot']['Swiss-Prot']
                    # Fall back to TrEMBL if needed
                    elif 'TrEMBL' in result['uniprot']:
                        uniprot_id = result['uniprot']['TrEMBL']

                    if uniprot_id:
                        # Handle both string and list responses
                        if isinstance(uniprot_id, list):
                            uniprot_id = uniprot_id[0]

                        mappings.append({
                            "gene_symbol": entrez_id,
                            "uniprot_id": uniprot_id
                        })
        except Exception as e:
            print(f"Error in batch {i//batch_size + 1}: {e}")

        # Be nice to the API
        time.sleep(0.5)

    df = pd.DataFrame(mappings)
    print(f"Successfully mapped {len(df)}/{len(gene_ids)} genes using mygene.info")
    return df


def get_entrez_ids(gene_symbols, batch_size=100):
    """Map gene symbols to Entrez IDs using mygene with batching"""
    print(f"Converting {len(gene_symbols)} gene symbols to Entrez IDs...")
    mg = mygene.MyGeneInfo()
    symbol_to_entrez = {}

    # Convert set to list if needed
    gene_symbols_list = list(gene_symbols) if isinstance(gene_symbols, set) else gene_symbols

    # Process in batches to avoid timeout issues
    for i in range(0, len(gene_symbols_list), batch_size):
        batch = gene_symbols_list[i:i + batch_size]
        try:
            results = mg.querymany(batch, scopes='symbol', fields='entrezgene', species='human')

            for result in results:
                if 'entrezgene' in result and not result.get('notfound', False):
                    symbol_to_entrez[result.get('query', '')] = str(result.get('entrezgene'))
        except Exception as e:
            print(f"Error in batch {i // batch_size + 1}: {e}")

        # Add slight delay to be nice to the API
        time.sleep(0.5)

    print(f"Successfully mapped {len(symbol_to_entrez)}/{len(gene_symbols)} genes to Entrez IDs")
    return symbol_to_entrez


def fetch_sequence_ncbi(protein_id, retries=3):
    """Fetch protein sequence from NCBI as a fallback"""
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

    for attempt in range(retries):
        try:
            # First get the GI number for the protein
            esearch_url = f"{base_url}esearch.fcgi?db=protein&term={protein_id}&retmode=json"
            search_response = requests.get(esearch_url, timeout=15)
            search_response.raise_for_status()
            search_data = search_response.json()

            if not search_data.get('esearchresult', {}).get('idlist', []):
                print(f"No NCBI record found for {protein_id}")
                time.sleep(1)
                continue

            ncbi_id = search_data['esearchresult']['idlist'][0]

            # Now fetch the sequence
            efetch_url = f"{base_url}efetch.fcgi?db=protein&id={ncbi_id}&rettype=fasta&retmode=text"
            fetch_response = requests.get(efetch_url, timeout=15)
            fetch_response.raise_for_status()

            # Parse the FASTA format
            fasta_io = StringIO(fetch_response.text)
            seq = str(next(SeqIO.parse(fasta_io, "fasta")).seq)
            return seq

        except Exception as e:
            print(f"NCBI error for {protein_id}: {e}, attempt {attempt+1}/{retries}")
            time.sleep(2)

    return None


def fetch_sequence(uniprot_id, retries=3):
    """Fetches protein sequence from UniProt or NCBI"""
    # Try UniProt first
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"

    for attempt in range(retries):
        try:
            r = requests.get(url, timeout=15)
            if r.status_code == 200:
                fasta_io = StringIO(r.text)
                seq = str(next(SeqIO.parse(fasta_io, "fasta")).seq)
                print(f"Fetched sequence for {uniprot_id}: {len(seq)} aa")
                return seq
            else:
                print(f"UniProt error for {uniprot_id} (HTTP {r.status_code})")
                break  # Move to NCBI faster if UniProt returns error
        except Exception as e:
            print(f"UniProt error for {uniprot_id}: {e}")
            time.sleep(1)

    # Try NCBI as fallback
    print(f"Trying NCBI for {uniprot_id}")
    seq = fetch_sequence_ncbi(uniprot_id, retries)
    if seq:
        print(f"Fetched sequence from NCBI for {uniprot_id}: {len(seq)} aa")
        return seq

    print(f"All attempts failed for {uniprot_id}")
    return None


def advanced_local_prediction(seq1, seq2):
    """Advanced local PPI prediction using multiple sequence features with updated Biopython methods"""
    if not seq1 or not seq2:
        return {"prediction": 0.0, "method": "local_failure", "reason": "missing_sequence"}

    try:
        # Use PairwiseAligner instead of pairwise2
        aligner = PairwiseAligner()
        aligner.mode = 'local'
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -0.5

        # Run alignment
        alignments = aligner.align(seq1, seq2)

        # Get the best alignment
        if len(alignments) == 0:
            return {"prediction": 0.2, "method": "local_default", "reason": "alignment_failed"}

        best_alignment = alignments[0]
        alignment_score = best_alignment.score

        # Calculate length-normalized score
        norm_score = alignment_score / (min(len(seq1), len(seq2)))

        # Calculate sequence length ratio as another feature
        length_ratio = min(len(seq1), len(seq2)) / max(len(seq1), len(seq2))

        # Calculate sequence composition similarity (simple approach)
        aa_counts1 = {aa: seq1.count(aa) / len(seq1) for aa in set(seq1)}
        aa_counts2 = {aa: seq2.count(aa) / len(seq2) for aa in set(seq2)}

        # Find common amino acids
        common_aa = set(aa_counts1.keys()) & set(aa_counts2.keys())
        composition_similarity = 0

        for aa in common_aa:
            composition_similarity += min(aa_counts1[aa], aa_counts2[aa])

        # Combine features with empirical weights
        combined_score = (0.5 * norm_score + 0.3 * length_ratio + 0.2 * composition_similarity)

        # Scale to [0-1] range using sigmoid function
        final_score = 1 / (1 + np.exp(-5 * (combined_score - 0.3)))

        return {
            "prediction": round(final_score, 4),
            "method": "local_advanced",
            "align_score": round(alignment_score, 2),
            "norm_score": round(norm_score, 4),
            "length_ratio": round(length_ratio, 4),
            "composition_sim": round(composition_similarity, 4)
        }

    except Exception as e:
        print(f"Error in advanced local prediction: {e}")
        return {"prediction": 0.1, "method": "local_error", "reason": str(e)}


def predict_ppi(protein_pairs, batch_size=5, output_dir=None):
    """Predicts protein-protein interactions using local methods since APIs fail"""
    results = []
    total_batches = (len(protein_pairs) + batch_size - 1) // batch_size

    # Pre-fetch all sequences to save time
    print("Pre-fetching protein sequences...")
    unique_proteins = set()
    for acc1, acc2 in protein_pairs:
        unique_proteins.add(acc1)
        unique_proteins.add(acc2)

    protein_sequences = {}
    with ThreadPoolExecutor(max_workers=5) as executor:
        futures = {executor.submit(fetch_sequence, p): p for p in unique_proteins}
        for future in tqdm(futures, total=len(unique_proteins)):
            protein = futures[future]
            try:
                seq = future.result()
                if seq:
                    protein_sequences[protein] = seq
            except Exception as e:
                print(f"Failed to fetch {protein}: {e}")

    print(f"Successfully fetched {len(protein_sequences)}/{len(unique_proteins)} protein sequences")

    # Now process protein pairs using the cached sequences
    for i in range(0, len(protein_pairs), batch_size):
        batch = protein_pairs[i:i+batch_size]
        print(f"\nProcessing batch {i//batch_size + 1}/{total_batches} ({len(batch)} pairs)...")

        batch_results = []
        for acc1, acc2 in batch:
            seq1 = protein_sequences.get(acc1)
            seq2 = protein_sequences.get(acc2)

            result = {
                "acc1": acc1,
                "acc2": acc2,
                "prediction": 0,
                "confidence": 0,
                "method": "unknown"
            }

            if not seq1 or not seq2:
                result["method"] = "failed_missing_sequence"
                batch_results.append(result)
                continue

            # Use advanced local prediction
            prediction_result = advanced_local_prediction(seq1, seq2)

            result["prediction"] = prediction_result.get("prediction", 0)
            result["confidence"] = 0.5  # Local method has limited confidence
            result["method"] = prediction_result.get("method", "local")
            result["details"] = prediction_result

            batch_results.append(result)

        results.extend(batch_results)

        # Save intermediate results to avoid losing progress
        if output_dir and i % (batch_size * 5) == 0:
            temp_df = pd.DataFrame(results)
            temp_df.to_csv(f"{output_dir}/ppi_predictions_batch_{i//batch_size}.csv", index=False)

        print(f"Batch {i//batch_size + 1} completed")

    return pd.DataFrame(results)


def main():
    print("Starting protein-protein interaction prediction pipeline...")

    # Create output directory
    out_dir = "DL-PPI outputs"
    os.makedirs(out_dir, exist_ok=True)

    # Load gene lists
    print("\nLoading gene lists...")
    try:
        comp_df = pd.read_csv("/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/KEGG outputs/kegg_complement_unique_genes.csv")
        deseq_df = pd.read_csv("/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/DESeq2 outputs/deseq2_results_M_OUD_vs_M_None.csv")

        # Get significant genes and complement genes
        sig = deseq_df[(deseq_df.padj < 0.05) & (abs(deseq_df.log2FoldChange) > 1)].gene.tolist()
        comp = comp_df.gene.tolist()

        # Filter out entries that aren't gene symbols
        comp = [c for c in comp if isinstance(c, str) and not c.startswith('[')]

        print(f"Found {len(sig)} significant genes and {len(comp)} complement genes")
    except Exception as e:
        print(f"Error loading gene lists: {e}")
        return

    # Convert gene symbols to Entrez IDs
    symbol_to_entrez = get_entrez_ids(set(sig + comp))

    # Convert genes to Entrez IDs
    sig_entrez = [symbol_to_entrez.get(s) for s in sig if s in symbol_to_entrez]
    comp_entrez = [symbol_to_entrez.get(c) for c in comp if c in symbol_to_entrez]

    # Remove None values
    sig_entrez = [s for s in sig_entrez if s]
    comp_entrez = [c for c in comp_entrez if c]

    print(f"Converted {len(sig_entrez)}/{len(sig)} significant genes to Entrez IDs")
    print(f"Converted {len(comp_entrez)}/{len(comp)} complement genes to Entrez IDs")

    # Create reverse mapping for later use
    entrez_to_symbol = {v: k for k, v in symbol_to_entrez.items()}

    # Map genes to UniProt - UPDATED TO USE MYGENE
    print("\nMapping genes to UniProt via mygene.info...")
    mapping = map_genes_to_uniprot_via_mygene(set(sig_entrez + comp_entrez))

    if mapping.empty:
        print("ERROR: Failed to map any genes to UniProt. Exiting.")
        return

    print(f"Successfully mapped {len(mapping)} genes to UniProt")
    print(f"Sample of mapping results:\n{mapping.head()}")

    # Save mapping results
    mapping.to_csv(f"{out_dir}/gene_to_uniprot_mapping.csv", index=False)

    # Add original gene symbols to mapping
    try:
        mapping['original_gene'] = mapping['gene_symbol'].apply(
            lambda x: entrez_to_symbol.get(x, "Unknown")
        )
        print("Added original gene symbols to mapping")
    except Exception as e:
        print(f"Error adding original gene symbols: {e}")

    # Build protein accession pairs
    print("\nBuilding protein pairs...")
    try:
        df_map = dict(zip(mapping.gene_symbol, mapping.uniprot_id))

        # Get which complement and significant genes have UniProt mappings
        comp_entrez_mapped = [c for c in comp_entrez if c in df_map]
        sig_entrez_mapped = [s for s in sig_entrez if s in df_map]

        print(f"Mapped complement genes: {len(comp_entrez_mapped)}/{len(comp_entrez)}")
        print(f"Mapped significant genes: {len(sig_entrez_mapped)}/{len(sig_entrez)}")

        # Build all possible pairs
        pairs = []
        for c in comp_entrez_mapped:
            for s in sig_entrez_mapped:
                pairs.append((df_map[c], df_map[s]))

        print(f"Created {len(pairs)} protein pairs for PPI prediction")

        if not pairs:
            print("No valid protein pairs created. Exiting.")
            return

        # For initial testing, use a smaller subset
        test_pairs = random.sample(pairs, min(10, len(pairs)))
        print(f"\nTesting PPI prediction with {len(test_pairs)} random pairs...")

        # Run predictions on test pairs
        df_results = predict_ppi(test_pairs, output_dir=out_dir)

        # Add gene symbols to results for readability
        reverse_map = {}
        for eid, uid in df_map.items():
            gene_symbol = entrez_to_symbol.get(eid, eid)
            reverse_map[uid] = gene_symbol

        df_results['gene1'] = df_results['acc1'].map(reverse_map)
        df_results['gene2'] = df_results['acc2'].map(reverse_map)

        # Clean and format for output
        if 'details' in df_results.columns:
            df_results = df_results.drop(columns=['details'])

        print(f"\nPPI prediction results:\n{df_results}")

        # Save results
        df_results.to_csv(f"{out_dir}/ppi_predictions.csv", index=False)
        print(f"Results saved to {out_dir}/ppi_predictions.csv")

        # Ask if user wants to run full prediction
        if len(pairs) > len(test_pairs):
            user_input = input(f"\nDo you want to run predictions for all {len(pairs)} pairs? (y/n): ")
            if user_input.lower() == 'y':
                print("\nRunning full PPI prediction...")
                full_results = predict_ppi(pairs, output_dir=out_dir)
                full_results['gene1'] = full_results['acc1'].map(reverse_map)
                full_results['gene2'] = full_results['acc2'].map(reverse_map)

                if 'details' in full_results.columns:
                    full_results = full_results.drop(columns=['details'])

                full_results.to_csv(f"{out_dir}/full_ppi_predictions.csv", index=False)
                print(f"Full results saved to {out_dir}/full_ppi_predictions.csv")

    except Exception as e:
        print(f"Error during protein pair processing: {e}")
        import traceback
        traceback.print_exc()

    print("\nPPI prediction pipeline completed")


if __name__ == "__main__":
    main()