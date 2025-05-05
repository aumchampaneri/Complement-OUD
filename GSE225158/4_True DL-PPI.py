import torch
from torch.utils.data import Dataset, DataLoader
from transformers import AutoTokenizer, AutoModel
import pandas as pd
import numpy as np
import os
import random
import requests
import time
from tqdm import tqdm
from tqdm.auto import tqdm as tqdm_auto
from Bio import SeqIO
from io import StringIO
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score


# 1. Define the Deep Learning PPI model using pre-trained ESM-2
class DeepLearningPPI:
    def __init__(self, model_name="facebook/esm2_t33_650M_UR50D"):
        """Initialize with a pre-trained protein language model"""
        print(f"Loading {model_name} model...")
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.tokenizer = AutoTokenizer.from_pretrained(model_name)
        self.model = AutoModel.from_pretrained(model_name).to(self.device)

        # Define a classifier head on top of the pre-trained model
        self.classifier = torch.nn.Sequential(
            torch.nn.Linear(2 * 1280, 512),  # 1280 is ESM-2's hidden size
            torch.nn.ReLU(),
            torch.nn.Dropout(0.1),
            torch.nn.Linear(512, 128),
            torch.nn.ReLU(),
            torch.nn.Dropout(0.1),
            torch.nn.Linear(128, 1),
            torch.nn.Sigmoid()
        ).to(self.device)

        print(f"Model loaded on {self.device}")

    def get_embeddings(self, sequence, max_length=1022):
        """Get protein embedding from pre-trained model"""
        # Truncate sequence if needed
        if len(sequence) > max_length:
            sequence = sequence[:max_length]

        # Tokenize and encode
        inputs = self.tokenizer(sequence, return_tensors="pt").to(self.device)

        # Get embeddings from the model
        with torch.no_grad():
            outputs = self.model(**inputs)

        # Use the [CLS] token embedding as protein representation
        embedding = outputs.last_hidden_state[:, 0, :].cpu().numpy()
        return embedding[0]

    def predict_interaction(self, seq1, seq2):
        """Predict interaction between two proteins"""
        # Get embeddings for both proteins
        emb1 = torch.tensor(self.get_embeddings(seq1)).unsqueeze(0).to(self.device)
        emb2 = torch.tensor(self.get_embeddings(seq2)).unsqueeze(0).to(self.device)

        # Concatenate embeddings
        concat_emb = torch.cat([emb1, emb2], dim=1)

        # Get prediction
        with torch.no_grad():
            prediction = self.classifier(concat_emb).item()

        return {
            "prediction": round(prediction, 4),
            "method": "deep_learning_esm2"
        }


# 2. Dataset class for protein pairs
class PPIDataset(Dataset):
    def __init__(self, protein_pairs, labels, protein_sequences, dl_model):
        self.protein_pairs = protein_pairs
        self.labels = labels
        self.protein_sequences = protein_sequences
        self.dl_model = dl_model

    def __len__(self):
        return len(self.protein_pairs)

    def __getitem__(self, idx):
        acc1, acc2 = self.protein_pairs[idx]
        seq1 = self.protein_sequences.get(acc1, "")
        seq2 = self.protein_sequences.get(acc2, "")

        # Get embeddings
        emb1 = torch.tensor(self.dl_model.get_embeddings(seq1))
        emb2 = torch.tensor(self.dl_model.get_embeddings(seq2))

        # Concatenate embeddings
        concat_emb = torch.cat([emb1, emb2], dim=0)

        return concat_emb, torch.tensor([self.labels[idx]], dtype=torch.float32)


# 3. Functions to get training data from public databases
def get_complement_proteins(csv_path="KEGG outputs/kegg_complement_unique_genes.csv"):
    """Load complement cascade proteins from CSV file"""
    try:
        # Check for absolute or relative path
        if not os.path.exists(csv_path):
            csv_path = f"/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/{csv_path}"
            if not os.path.exists(csv_path):
                raise FileNotFoundError(f"Could not find {csv_path}")

        comp_df = pd.read_csv(csv_path)
        complement_genes = comp_df['gene'].tolist()

        # Filter out entries that aren't gene symbols
        complement_genes = [c for c in complement_genes if isinstance(c, str) and not c.startswith('[')]
        print(f"Loaded {len(complement_genes)} complement cascade genes from {csv_path}")
        return complement_genes
    except Exception as e:
        print(f"Error loading complement genes: {e}")
        # Fallback to a minimal list if file loading fails
        fallback_genes = ["C3", "C5", "CFB", "CFD", "CFH"]
        print(f"Using fallback list of {len(fallback_genes)} genes")
        return fallback_genes


def get_string_interactions(gene_list, score_threshold=700, species=9606):
    """Get protein interactions from STRING database"""
    string_api_url = "https://string-db.org/api/tsv/network"

    # Get interactions
    params = {
        "identifiers": "%0d".join(gene_list),
        "species": species,
        "network_type": "physical",
        "required_score": score_threshold
    }

    response = requests.get(string_api_url, params=params)
    interactions_df = pd.read_csv(StringIO(response.text), sep='\t')

    print(f"Downloaded {len(interactions_df)} interactions from STRING")
    return interactions_df


def get_uniprot_sequences(gene_list, species="human"):
    """Get protein sequences from UniProt for a list of genes"""
    sequences = {}
    uniprot_url = "https://rest.uniprot.org/uniprotkb/search?query="

    for gene in tqdm(gene_list, desc="Getting UniProt sequences"):
        try:
            query = f"{gene}+AND+organism_id:9606+AND+reviewed:true"
            response = requests.get(f"{uniprot_url}{query}&format=fasta")

            if response.status_code == 200 and response.text:
                fasta_data = StringIO(response.text)
                for record in SeqIO.parse(fasta_data, "fasta"):
                    acc = record.id.split('|')[1] if '|' in record.id else record.id
                    sequences[gene] = {
                        "uniprot_id": acc,
                        "sequence": str(record.seq)
                    }
                    break  # Take first reviewed entry

            time.sleep(0.5)  # Be gentle with the API
        except Exception as e:
            print(f"Error getting sequence for {gene}: {e}")

    return sequences


def create_training_data(interactions_df, sequences, complement_genes):
    """Create positive and negative examples for training"""
    # Create mapping from gene name to uniprot id
    gene_to_uniprot = {gene: data["uniprot_id"] for gene, data in sequences.items()
                       if gene in sequences}

    # Create positive examples (known interactions)
    positive_pairs = []
    for _, row in interactions_df.iterrows():
        gene1, gene2 = row['preferredName_A'], row['preferredName_B']

        if gene1 in gene_to_uniprot and gene2 in gene_to_uniprot:
            positive_pairs.append((gene_to_uniprot[gene1], gene_to_uniprot[gene2]))

    # Create negative examples (random pairs not in positive set)
    negative_pairs = []
    uniprot_ids = list(gene_to_uniprot.values())
    positive_set = set(tuple(sorted(pair)) for pair in positive_pairs)

    while len(negative_pairs) < len(positive_pairs):
        gene1 = random.choice(uniprot_ids)
        gene2 = random.choice(uniprot_ids)
        if gene1 != gene2 and tuple(sorted((gene1, gene2))) not in positive_set:
            negative_pairs.append((gene1, gene2))
            positive_set.add(tuple(sorted((gene1, gene2))))

    # Create labels
    pairs = positive_pairs + negative_pairs
    labels = [1] * len(positive_pairs) + [0] * len(negative_pairs)

    # Verify that sequences exist for all pairs
    protein_sequences = {data["uniprot_id"]: data["sequence"]
                         for gene, data in sequences.items()}

    valid_pairs = []
    valid_labels = []
    for pair, label in zip(pairs, labels):
        if pair[0] in protein_sequences and pair[1] in protein_sequences:
            valid_pairs.append(pair)
            valid_labels.append(label)

    return valid_pairs, valid_labels, protein_sequences


def create_complement_training_data():
    """Create a complement cascade training dataset"""
    # Get complement proteins from CSV
    complement_genes = get_complement_proteins()
    print(f"Working with {len(complement_genes)} complement cascade proteins")

    # Get interactions from STRING
    interactions_df = get_string_interactions(complement_genes)

    # Expand gene set to include interactors
    expanded_genes = list(set(
        list(interactions_df['preferredName_A']) +
        list(interactions_df['preferredName_B'])
    ))
    print(f"Expanded to {len(expanded_genes)} genes including interactors")

    # Get protein sequences
    sequences = get_uniprot_sequences(expanded_genes)

    # Create training data
    pairs, labels, protein_sequences = create_training_data(
        interactions_df, sequences, complement_genes)

    print(f"Created {len(pairs)} training pairs ({sum(labels)} positive, {len(labels) - sum(labels)} negative)")

    # Save to files
    os.makedirs("DL-PPI outputs/data", exist_ok=True)

    # Save interactions
    interactions_df.to_csv("DL-PPI outputs/data/complement_interactions.csv", index=False)

    # Save training data
    training_df = pd.DataFrame({
        'protein1': [p[0] for p in pairs],
        'protein2': [p[1] for p in pairs],
        'interaction_label': labels
    })
    training_df.to_csv("DL-PPI outputs/data/known_interactions.csv", index=False)

    # Save sequences as FASTA
    with open("DL-PPI outputs/data/protein_sequences.fasta", 'w') as f:
        for acc, seq in protein_sequences.items():
            f.write(f">{acc}\n{seq}\n")

    return training_df, protein_sequences


# 4. Training function for the model
def train_ppi_model(train_pairs, train_labels, protein_sequences, batch_size=8, epochs=5):
    """Train the deep learning PPI model on your data"""
    # Initialize model
    dl_model = DeepLearningPPI()

    # Create datasets
    train_pairs_val, val_pairs, train_labels_val, val_labels = train_test_split(
        train_pairs, train_labels, test_size=0.2, random_state=42)

    # Create datasets
    train_dataset = PPIDataset(train_pairs_val, train_labels_val, protein_sequences, dl_model)
    val_dataset = PPIDataset(val_pairs, val_labels, protein_sequences, dl_model)

    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size)

    # Define optimizer and loss
    optimizer = torch.optim.Adam(dl_model.classifier.parameters(), lr=1e-4)
    criterion = torch.nn.BCELoss()

    # Training loop
    for epoch in range(epochs):
        dl_model.classifier.train()
        total_loss = 0

        for embeddings, labels in tqdm_auto(train_loader, desc=f"Epoch {epoch + 1}/{epochs}"):
            embeddings = embeddings.to(dl_model.device)
            labels = labels.to(dl_model.device)

            optimizer.zero_grad()
            outputs = dl_model.classifier(embeddings)
            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()

            total_loss += loss.item()

        # Validation
        dl_model.classifier.eval()
        val_loss = 0
        all_preds = []
        all_labels = []

        with torch.no_grad():
            for embeddings, labels in val_loader:
                embeddings = embeddings.to(dl_model.device)
                labels = labels.to(dl_model.device)

                outputs = dl_model.classifier(embeddings)
                loss = criterion(outputs, labels)
                val_loss += loss.item()

                all_preds.extend(outputs.cpu().numpy())
                all_labels.extend(labels.cpu().numpy())

        # Calculate metrics
        auc = roc_auc_score(all_labels, all_preds)
        print(f"Epoch {epoch + 1}: Train Loss: {total_loss / len(train_loader):.4f}, "
              f"Val Loss: {val_loss / len(val_loader):.4f}, AUC: {auc:.4f}")

    # Save the trained model
    os.makedirs("DL-PPI outputs/models", exist_ok=True)
    torch.save(dl_model.classifier.state_dict(), "DL-PPI outputs/models/classifier_head.pt")

    return dl_model


# 5. Prediction function using deep learning model
def predict_ppi_deep_learning(protein_pairs, protein_sequences, batch_size=5, output_dir=None):
    """Predicts PPIs using deep learning model"""
    results = []
    total_batches = (len(protein_pairs) + batch_size - 1) // batch_size

    # Initialize the DL model
    dl_model = DeepLearningPPI()

    # Load trained model if exists
    if os.path.exists("DL-PPI outputs/models/classifier_head.pt"):
        print("Loading pre-trained classifier head...")
        dl_model.classifier.load_state_dict(
            torch.load("DL-PPI outputs/models/classifier_head.pt")
        )
    else:
        print("WARNING: No pre-trained model found. Using untrained model.")

    # Process protein pairs
    for i in range(0, len(protein_pairs), batch_size):
        batch = protein_pairs[i:i + batch_size]
        print(f"\nProcessing batch {i // batch_size + 1}/{total_batches} ({len(batch)} pairs)...")

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

            # Use deep learning model for prediction
            prediction_result = dl_model.predict_interaction(seq1, seq2)

            result["prediction"] = prediction_result["prediction"]
            result["confidence"] = 0.8  # Higher confidence with deep learning
            result["method"] = prediction_result["method"]

            batch_results.append(result)

        results.extend(batch_results)

        # Save intermediate results
        if output_dir and i % (batch_size * 5) == 0:
            temp_df = pd.DataFrame(results)
            temp_df.to_csv(f"{output_dir}/dl_ppi_predictions_batch_{i // batch_size}.csv", index=False)

        print(f"Batch {i // batch_size + 1} completed")

    return pd.DataFrame(results)


# 6. Main function
def main():
    print("Starting deep learning protein-protein interaction prediction pipeline...")

    # Create output directory
    out_dir = "DL-PPI outputs"
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(f"{out_dir}/data", exist_ok=True)

    # Check if we have existing training data
    if not os.path.exists(f"{out_dir}/data/known_interactions.csv"):
        print("\nCreating complement cascade training data...")
        training_df, protein_sequences_dict = create_complement_training_data()
    else:
        print("\nLoading existing training data...")
        training_df = pd.read_csv(f"{out_dir}/data/known_interactions.csv")

        # Load protein sequences from FASTA
        protein_sequences_dict = {}
        with open(f"{out_dir}/data/protein_sequences.fasta", 'r') as f:
            current_id = None
            current_seq = ""
            for line in f:
                if line.startswith('>'):
                    if current_id:
                        protein_sequences_dict[current_id] = current_seq
                    current_id = line.strip()[1:]  # Remove '>'
                    current_seq = ""
                else:
                    current_seq += line.strip()
            if current_id:
                protein_sequences_dict[current_id] = current_seq

    # Create training pairs and labels
    train_pairs = list(zip(training_df['protein1'], training_df['protein2']))
    train_labels = training_df['interaction_label'].tolist()

    print(
        f"Loaded {len(train_pairs)} training pairs, {sum(train_labels)} positive, {len(train_labels) - sum(train_labels)} negative")

    # Check if model exists or needs training
    if os.path.exists(f"{out_dir}/models/classifier_head.pt"):
        print("\nFound existing trained model.")
        # Load model for inference
        dl_model = DeepLearningPPI()
        dl_model.classifier.load_state_dict(torch.load(f"{out_dir}/models/classifier_head.pt"))
    else:
        print("\nTraining new model...")
        dl_model = train_ppi_model(train_pairs, train_labels, protein_sequences_dict)

    # Generate test pairs from complement genes and DEGs
    print("\nGenerating test pairs from complement genes and DEGs...")
    comp_file = "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/KEGG outputs/kegg_complement_unique_genes.csv"
    deg_file = "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/DESeq2 outputs/deseq2_results_F_OUD_vs_F_None.csv"

    # Load gene lists
    try:
        comp_df = pd.read_csv(comp_file)
        deseq_df = pd.read_csv(deg_file)

        # Get significant genes and complement genes
        sig = deseq_df[(deseq_df.padj < 0.05) & (abs(deseq_df.log2FoldChange) > 1)].gene.tolist()
        comp = comp_df.gene.tolist()

        # Filter out entries that aren't gene symbols
        comp = [c for c in comp if isinstance(c, str) and not c.startswith('[')]

        print(f"Found {len(sig)} significant genes and {len(comp)} complement genes")
    except Exception as e:
        print(f"Error loading gene lists: {e}")
        return

    # Map genes to UniProt IDs
    print("\nMapping genes to UniProt IDs...")
    from mygene import MyGeneInfo
    mg = MyGeneInfo()

    # Get mappings for all genes
    all_genes = set(comp + sig)
    gene_info = mg.querymany(all_genes, scopes='symbol', fields='uniprot', species='human', returnall=True)

    # Process results
    gene_to_uniprot = {}
    for entry in gene_info['out']:
        if 'uniprot' in entry and 'Swiss-Prot' in entry['uniprot']:
            gene_symbol = entry['query']
            uniprot_id = entry['uniprot']['Swiss-Prot']
            if isinstance(uniprot_id, list):  # Handle multiple UniProt IDs
                uniprot_id = uniprot_id[0]
            gene_to_uniprot[gene_symbol] = uniprot_id

    # Map genes to UniProt IDs
    comp_uniprot = [gene_to_uniprot.get(g) for g in comp if g in gene_to_uniprot]
    sig_uniprot = [gene_to_uniprot.get(g) for g in sig if g in gene_to_uniprot]

    # Remove None values
    comp_uniprot = [c for c in comp_uniprot if c]
    sig_uniprot = [s for s in sig_uniprot if s]

    print(f"Mapped {len(comp_uniprot)}/{len(comp)} complement genes to UniProt")
    print(f"Mapped {len(sig_uniprot)}/{len(sig)} significant genes to UniProt")

    # Create test pairs
    test_pairs = []
    for c in comp_uniprot:
        for s in sig_uniprot:
            if c in protein_sequences_dict and s in protein_sequences_dict:
                test_pairs.append((c, s))

    print(f"Created {len(test_pairs)} protein pairs for prediction")

    # Run predictions
    print("\nPredicting interactions...")
    results_df = predict_ppi_deep_learning(test_pairs, protein_sequences_dict,
                                           batch_size=5, output_dir=out_dir)

    # Save final results
    results_df.to_csv(f"{out_dir}/dl_ppi_predictions.csv", index=False)

    # Create reverse mapping for visualization
    uniprot_to_gene = {v: k for k, v in gene_to_uniprot.items()}

    # Create visualization-friendly format
    viz_df = pd.DataFrame({
        'gene1': [uniprot_to_gene.get(p[0], p[0]) for p in test_pairs],
        'gene2': [uniprot_to_gene.get(p[1], p[1]) for p in test_pairs],
        'prediction': results_df['prediction'].values,
        'method': results_df['method'].values
    })

    # Keep only high-confidence predictions for visualization
    viz_df = viz_df[viz_df['prediction'] >= 0.4]
    viz_df.to_csv(f"{out_dir}/dl_ppi_viz_data.csv", index=False)

    print(f"\nFound {len(viz_df)} high-confidence interactions.")
    print(f"Results saved to {out_dir}/dl_ppi_predictions.csv")

    # Visualization (if available)
    try:
        from visualization import (interactive_ppi_network, create_3d_ppi_network,
                                   create_sankey_diagram, plot_network_metrics)

        print("\nGenerating visualizations...")
        print("Generating interactive PPI network...")
        interactive_ppi_network(viz_df, min_score=0.5)

        print("Generating 3D PPI network...")
        create_3d_ppi_network(viz_df, min_score=0.5)

        print("Generating Sankey diagram...")
        create_sankey_diagram(viz_df, min_score=0.7)

        print("Analyzing network metrics...")
        plot_network_metrics(viz_df, min_score=0.5)
    except ImportError:
        print("\nVisualization module not found. Skipping visualizations.")

    print("\nDL-PPI prediction pipeline completed.")


if __name__ == "__main__":
    main()