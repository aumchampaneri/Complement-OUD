import os
import torch
import pandas as pd
import numpy as np
import json
import time
import gc
import pickle
from tqdm import tqdm
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from scipy import stats
import requests

os.environ["OMP_NUM_THREADS"] = "1"  # Limit OpenMP threads
torch.set_num_threads(4)  # Reduce parallel threads

# Constants
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
torch.manual_seed(RANDOM_SEED)


# Memory optimization for M1 Mac
def optimize_for_m1_mac():
    """Configure environment for optimal performance on M1 Mac with limited RAM"""
    # Enable MPS (Metal Performance Shaders) for M1 acceleration
    if hasattr(torch.backends, 'mps') and torch.backends.mps.is_available():
        device = torch.device('mps')
        print("Using MPS (Apple Silicon GPU) acceleration")
    else:
        device = torch.device('cpu')
        print("MPS not available, using CPU")

    # Limit memory usage
    os.environ["TOKENIZERS_PARALLELISM"] = "false"
    torch.set_num_threads(4)  # Prevent thread explosion

    return device


# Function to ensure directory exists
def ensure_path(path):
    os.makedirs(path, exist_ok=True)
    return path


# Function to load cached protein data instead of making API calls
def load_cached_protein_data(gene_list, cache_dir="protein_cache"):
    """Load pre-cached protein sequences and gene-to-UniProt mappings"""
    import os
    import pickle

    gene_to_uniprot = {}
    protein_sequences_dict = {}

    # Define cache file paths
    gene_to_uniprot_file = os.path.join(cache_dir, "gene_to_uniprot.pkl")
    protein_sequences_file = os.path.join(cache_dir, "protein_sequences.pkl")

    # Check if cache exists
    if not os.path.exists(gene_to_uniprot_file) or not os.path.exists(protein_sequences_file):
        print(f"Warning: Protein cache not found in {cache_dir}. Run protein_downloader.py first.")
        return {}, {}

    # Load cached gene-to-UniProt mappings
    with open(gene_to_uniprot_file, 'rb') as f:
        full_gene_to_uniprot = pickle.load(f)

    # Load cached protein sequences
    with open(protein_sequences_file, 'rb') as f:
        protein_sequences_dict = pickle.load(f)

    # Filter gene mappings to only include genes in our list
    for gene in gene_list:
        if gene in full_gene_to_uniprot:
            gene_to_uniprot[gene] = full_gene_to_uniprot[gene]

    print(f"Loaded {len(protein_sequences_dict)} protein sequences from cache")
    print(f"Found UniProt mappings for {len(gene_to_uniprot)} out of {len(gene_list)} genes")

    return gene_to_uniprot, protein_sequences_dict


# Enhanced DeepLearningPPI class with memory optimizations
class DeepLearningPPI:
    def __init__(self, model_name="facebook/esm2_t33_650M_UR50D", embedding_cache=None):
        from transformers import AutoTokenizer, EsmModel, AutoConfig
        import torch
        import torch.nn as nn

        # Set device to MPS if available
        self.device = optimize_for_m1_mac()

        # Initialize embedding cache
        if embedding_cache is None:
            self.embedding_cache = {}
        else:
            self.embedding_cache = embedding_cache
            print(f"Loaded {len(self.embedding_cache)} embeddings from cache")

        print(f"Loading {model_name} model with memory optimizations...")

        # Load model with memory optimizations
        config = AutoConfig.from_pretrained(model_name)
        self.tokenizer = AutoTokenizer.from_pretrained(model_name)
        self.model = EsmModel.from_pretrained(model_name, config=config).to(self.device)
        print(f"Model loaded on {self.device}")

        # Create classification head
        embedding_dim = self.model.config.hidden_size
        self.classifier = nn.Sequential(
            nn.Linear(embedding_dim * 2, 512),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(512, 256),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(256, 1),
            nn.Sigmoid()
        ).to(self.device)

    def get_embeddings(self, sequence, max_length=512, use_cache=True):
        if not isinstance(sequence, str):
            print(f"Input sequence provided is not a string: {type(sequence)}. Converting...")
            sequence = str(sequence)

        # Check if we already have this embedding cached
        if use_cache and sequence in self.embedding_cache:
            return self.embedding_cache[sequence]

        # Truncate long sequences to max_length
        if len(sequence) > max_length:
            sequence = sequence[:max_length]

        # Tokenize and get embeddings
        try:
            inputs = self.tokenizer(sequence, return_tensors="pt", padding=True).to(self.device)
            with torch.no_grad():
                outputs = self.model(**inputs)

            # Use the mean of the last hidden state (CLS token)
            embeddings = outputs.last_hidden_state[:, 0, :].cpu().numpy()

            # Cache the embedding
            if use_cache:
                self.embedding_cache[sequence] = embeddings[0]

            return embeddings[0]
        except Exception as e:
            print(f"Error getting embeddings: {e}")
            # Return zero embedding vector if there's an error
            return np.zeros(self.model.config.hidden_size)

    def predict_interaction(self, seq1, seq2):
        # Get embeddings for both proteins
        emb1 = self.get_embeddings(seq1)
        emb2 = self.get_embeddings(seq2)

        # Concatenate embeddings
        combined = np.concatenate([emb1, emb2])
        combined_tensor = torch.tensor(combined, dtype=torch.float32).unsqueeze(0).to(self.device)

        # Predict interaction
        self.classifier.eval()
        with torch.no_grad():
            pred = self.classifier(combined_tensor)

        return pred.item()

    def load_embedding_cache(self, cache_file):
        if os.path.exists(cache_file):
            try:
                with open(cache_file, 'rb') as f:
                    self.embedding_cache = pickle.load(f)
                print(f"Loaded {len(self.embedding_cache)} embeddings from cache")
            except Exception as e:
                print(f"Error loading embedding cache: {e}")
                self.embedding_cache = {}
        else:
            print("No embedding cache found, starting with empty cache")
            self.embedding_cache = {}


# Enhanced function to load inflammatory pathway data
def load_inflammatory_pathways(
        pathway_file="/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/KEGG outputs/kegg_inflammatory_pathways.csv"):
    """Load inflammatory pathway data and handle multiple pathways per gene"""
    try:
        pathway_df = pd.read_csv(pathway_file)
        print(f"Loaded {len(pathway_df)} gene-pathway associations")

        # Count pathways and genes
        unique_pathways = pathway_df['pathway'].unique()
        print(f"Found {len(unique_pathways)} unique pathways")

        # Show top pathways by gene count
        pathway_counts = pathway_df.groupby('pathway')['gene'].count().sort_values(ascending=False)
        print("Top 5 pathways by gene count:")
        for p, count in pathway_counts.head(5).items():
            print(f"  - {p}: {count} genes")

        # Count multi-pathway genes
        gene_pathway_counts = pathway_df.groupby('gene')['pathway'].count()
        multi_pathway_genes = gene_pathway_counts[gene_pathway_counts > 1]
        print(f"Found {len(multi_pathway_genes)} genes that belong to multiple pathways")

        return pathway_df
    except Exception as e:
        print(f"Error loading pathway data: {e}")
        return pd.DataFrame()


# Function to load sex-specific genes
def load_sex_specific_genes(sex="M"):
    """Load differentially expressed genes from DESeq2 results"""
    try:
        if sex.upper() == "M":
            file_path = "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/DESeq2 outputs/deseq2_results_M_OUD_vs_M_None.csv"
        else:
            file_path = "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/DESeq2 outputs/deseq2_results_F_OUD_vs_F_None.csv"

        deg_df = pd.read_csv(file_path)
        gene_list = deg_df['gene'].tolist()
        print(f"Loaded {len(gene_list)} {sex} sex-specific genes")
        return gene_list
    except Exception as e:
        print(f"Error loading {sex} sex-specific genes: {e}")
        return []


# Function to create complement training data
def create_complement_training_data(output_dir="DL-PPI outputs"):
    """Create training data for the model"""
    import requests
    import re
    from io import StringIO

    print("Creating training data for protein-protein interactions...")

    # Get known PPIs from STRING database
    try:
        # Known C5/C3AR1 interaction
        known_interactions = [
            {"protein1": "P01031", "protein2": "Q16581", "interaction": 1},  # C5 and C3AR1
        ]

        # Add more known interactions from complement pathway
        complement_pairs = [
            ("P01024", "P17927"),  # C3 and CR1
            ("P01024", "P15529"),  # C3 and CD46
            ("P01024", "P13987"),  # C3 and CD59
            ("P01031", "P21730"),  # C5 and C5AR1
            ("P02748", "P13987"),  # C9 and CD59
            ("P00736", "P13987"),  # C1R and CD59
            ("P00736", "P05155"),  # C1R and SERPING1
        ]

        for pair in complement_pairs:
            known_interactions.append({"protein1": pair[0], "protein2": pair[1], "interaction": 1})

        # Create equal number of negative samples (non-interacting proteins)
        # Use proteins from different cellular compartments as negative samples
        negative_pairs = [
            ("P01024", "P01375"),  # C3 and TNF
            ("P01031", "P01584"),  # C5 and IL1B
            ("P02748", "P05231"),  # C9 and IL6
            ("P00736", "P10145"),  # C1R and CXCL8
            ("P13987", "P16581"),  # CD59 and ELAM1
            ("P15529", "P29460"),  # CD46 and IL12B
            ("P17927", "P01375"),  # CR1 and TNF
        ]

        for pair in negative_pairs:
            known_interactions.append({"protein1": pair[0], "protein2": pair[1], "interaction": 0})

        # Download protein sequences
        protein_sequences_dict = {}
        for interaction in known_interactions:
            for protein_key in ["protein1", "protein2"]:
                uniprot_id = interaction[protein_key]
                if uniprot_id not in protein_sequences_dict:
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
                    except Exception as e:
                        print(f"Error downloading sequence for {uniprot_id}: {e}")

        # Save data to files
        ensure_path(f"{output_dir}/data")

        # Create dataframe
        training_df = pd.DataFrame(known_interactions)
        training_df.columns = ['protein1', 'protein2', 'interaction_label']

        # Shuffle data
        training_df = training_df.sample(frac=1, random_state=RANDOM_SEED).reset_index(drop=True)

        # Save training data
        training_df.to_csv(f"{output_dir}/data/known_interactions.csv", index=False)

        # Save protein sequences
        with open(f"{output_dir}/data/protein_sequences.pkl", 'wb') as f:
            pickle.dump(protein_sequences_dict, f)

        print(f"Created training data with {len(training_df)} protein interactions")
        print(f"Saved protein sequences for {len(protein_sequences_dict)} proteins")

        return training_df, protein_sequences_dict

    except Exception as e:
        print(f"Error creating training data: {e}")
        return pd.DataFrame(), {}


# Helper to generate pathway-based pairs
def generate_pathway_based_pairs(gene_list, pathway_df, protein_sequences_dict, gene_to_uniprot):
    """Generate protein pairs within pathways"""
    # Create pathway groups
    pathway_groups = pathway_df.groupby('pathway')['gene'].apply(list).to_dict()

    pairs = []
    for pathway, genes in pathway_groups.items():
        # Filter genes to include only those in both our gene list and pathway
        pathway_genes = list(set(genes).intersection(set(gene_list)))

        # Get UniProt IDs with sequences
        uniprot_ids = [gene_to_uniprot.get(g) for g in pathway_genes]
        uniprot_ids_raw = [u for u in uniprot_ids if u]
        uniprot_ids = [u for u in uniprot_ids_raw if u in protein_sequences_dict]
        print(
            f"Pathway {pathway}: Found {len(pathway_genes)} genes, {len(uniprot_ids_raw)} with UniProt IDs, "
            f"{len(uniprot_ids)} with sequences"
        )

        # Create all protein pairs within this pathway
        for i in range(len(uniprot_ids)):
            for j in range(i + 1, len(uniprot_ids)):
                pairs.append((uniprot_ids[i], uniprot_ids[j]))

    print(f"Generated {len(pairs)} pathway-based protein pairs")
    return pairs


# Helper to get pathway information for a gene pair
def get_pathway_info(gene1, gene2, pathway_df):
    """Get pathway information for a gene pair"""
    gene1_pathways = pathway_df[pathway_df['gene'] == gene1]['pathway'].unique().tolist()
    gene2_pathways = pathway_df[pathway_df['gene'] == gene2]['pathway'].unique().tolist()

    # Find common pathways
    common_pathways = list(set(gene1_pathways) & set(gene2_pathways))

    return {
        'gene1_pathways': gene1_pathways,
        'gene2_pathways': gene2_pathways,
        'common_pathways': common_pathways,
        'has_common_pathway': len(common_pathways) > 0
    }


# Dataset class for protein-protein interactions
class PPIDataset(Dataset):
    def __init__(self, pairs, labels, sequences_dict, model):
        self.pairs = pairs
        self.labels = labels
        self.sequences_dict = sequences_dict
        self.model = model

        # Precompute and cache embeddings for all proteins to speed up training
        for pair in pairs:
            for protein in pair:
                if protein in sequences_dict and protein not in model.embedding_cache:
                    model.get_embeddings(sequences_dict[protein], use_cache=True)

    def __len__(self):
        return len(self.pairs)

    def __getitem__(self, idx):
        protein1, protein2 = self.pairs[idx]

        # Get sequences
        seq1 = self.sequences_dict[protein1]
        seq2 = self.sequences_dict[protein2]

        # Get embeddings
        emb1 = self.model.get_embeddings(seq1)
        emb2 = self.model.get_embeddings(seq2)

        # Concatenate embeddings
        combined = np.concatenate([emb1, emb2])

        # Get label
        label = self.labels[idx]

        return torch.tensor(combined, dtype=torch.float32), torch.tensor([label], dtype=torch.float32)


# Training function with memory optimization
def train_ppi_model(train_pairs, train_labels, protein_sequences_dict,
                    batch_size=8, epochs=20, patience=3, lr_scheduler=True,
                    output_dir="DL-PPI outputs/models"):
    """Train model with memory optimizations for M1 Mac"""
    print("Initializing model...")
    model = DeepLearningPPI()

    # Create dataset and dataloader
    dataset = PPIDataset(train_pairs, train_labels, protein_sequences_dict, model)

    # Split data for training and validation
    train_indices, val_indices = train_test_split(
        range(len(dataset)),
        test_size=0.2,
        random_state=RANDOM_SEED,
        stratify=train_labels
    )

    train_sampler = torch.utils.data.SubsetRandomSampler(train_indices)
    val_sampler = torch.utils.data.SubsetRandomSampler(val_indices)

    train_loader = DataLoader(
        dataset,
        batch_size=batch_size,
        sampler=train_sampler,
        num_workers=0  # No multiprocessing for RAM conservation
    )

    val_loader = DataLoader(
        dataset,
        batch_size=batch_size,
        sampler=val_sampler,
        num_workers=0  # No multiprocessing for RAM conservation
    )

    # Initialize optimizer
    optimizer = torch.optim.Adam(model.classifier.parameters(), lr=1e-4, weight_decay=1e-5)

    # Learning rate scheduler
    if lr_scheduler:
        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimizer, mode='min', factor=0.5, patience=2, verbose=True
        )

    # Loss function
    criterion = torch.nn.BCELoss()

    # Early stopping variables
    best_val_loss = float('inf')
    patience_counter = 0

    print("Starting training...")
    for epoch in range(epochs):
        # Training phase
        model.classifier.train()
        train_loss = 0
        train_count = 0

        for inputs, targets in train_loader:
            inputs = inputs.to(model.device)
            targets = targets.to(model.device)

            # Forward pass
            optimizer.zero_grad()
            outputs = model.classifier(inputs)
            loss = criterion(outputs, targets)

            # Backward pass and update
            loss.backward()
            optimizer.step()

            # Accumulate loss
            train_loss += loss.item() * inputs.size(0)
            train_count += inputs.size(0)

            # Memory management
            del inputs, targets, outputs
            torch.cuda.empty_cache() if torch.cuda.is_available() else None

        train_loss = train_loss / train_count

        # Validation phase
        model.classifier.eval()
        val_loss = 0
        val_correct = 0
        val_count = 0

        with torch.no_grad():
            for inputs, targets in val_loader:
                inputs = inputs.to(model.device)
                targets = targets.to(model.device)

                # Forward pass
                outputs = model.classifier(inputs)
                loss = criterion(outputs, targets)

                # Calculate accuracy
                preds = (outputs > 0.5).float()
                val_correct += (preds == targets).sum().item()

                # Accumulate loss
                val_loss += loss.item() * inputs.size(0)
                val_count += inputs.size(0)

                # Memory management
                del inputs, targets, outputs, preds
                torch.cuda.empty_cache() if torch.cuda.is_available() else None

        val_loss = val_loss / val_count
        val_acc = val_correct / val_count

        # Update learning rate if needed
        if lr_scheduler:
            scheduler.step(val_loss)

        # Print progress
        print(
            f"Epoch {epoch + 1}/{epochs} - Train Loss: {train_loss:.4f}, Val Loss: {val_loss:.4f}, Val Acc: {val_acc:.4f}")

        # Check for early stopping
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            patience_counter = 0
            # Save best model
            os.makedirs(output_dir, exist_ok=True)
            torch.save(model.classifier.state_dict(), f"{output_dir}/classifier_head.pt")
            print(f"Saved model checkpoint with validation loss {val_loss:.4f}")
        else:
            patience_counter += 1
            if patience_counter >= patience:
                print(f"Early stopping at epoch {epoch + 1}")
                break

        # Garbage collection
        gc.collect()

    # Load best model
    model.classifier.load_state_dict(torch.load(f"{output_dir}/classifier_head.pt"))
    return model


# Function to predict PPIs
def predict_ppi_deep_learning(protein_pairs, protein_sequences_dict, model=None,
                              batch_size=4, model_path="DL-PPI outputs/models/classifier_head.pt",
                              output_dir="DL-PPI outputs/results"):
    """Predict protein-protein interactions using deep learning"""
    print(f"Predicting interactions for {len(protein_pairs)} protein pairs...")

    # Create model if not provided
    if model is None:
        model = DeepLearningPPI()
        model.classifier.load_state_dict(torch.load(model_path))

    # Ensure model is in evaluation mode
    model.classifier.eval()

    # Process in batches to save memory
    predictions = []

    # Create batches
    batches = [protein_pairs[i:i + batch_size] for i in range(0, len(protein_pairs), batch_size)]

    with torch.no_grad():
        for i, batch in enumerate(tqdm(batches, desc="Predicting")):
            batch_predictions = []

            for pair in batch:
                protein1, protein2 = pair

                # Check if both proteins are in sequences
                if protein1 not in protein_sequences_dict or protein2 not in protein_sequences_dict:
                    print(f"Warning: Missing sequence for {protein1} or {protein2}")
                    continue

                # Get sequences
                seq1 = protein_sequences_dict[protein1]
                seq2 = protein_sequences_dict[protein2]

                # Predict interaction
                pred = model.predict_interaction(seq1, seq2)

                # Store result
                batch_predictions.append({
                    'acc1': protein1,
                    'acc2': protein2,
                    'prediction': pred,
                    'method': 'dl-ppi'
                })

            # Free memory immediately
            torch.cuda.empty_cache() if torch.cuda.is_available() else None

            # Add to overall predictions
            predictions.extend(batch_predictions)

            # Perform garbage collection every few batches
            if i % 5 == 0:
                gc.collect()

    # Create results dataframe
    results_df = pd.DataFrame(predictions)

    # Save results
    os.makedirs(output_dir, exist_ok=True)
    results_df.to_csv(f"{output_dir}/dl_ppi_predictions.csv", index=False)

    print(f"Prediction complete. Results saved to {output_dir}/dl_ppi_predictions.csv")
    return results_df


# Helper to precompute all embeddings
def precompute_all_embeddings(protein_sequences_dict, model, cache_file):
    """Precompute embeddings for all proteins"""
    print("Precomputing embeddings for all proteins...")

    # Check if already cached
    if hasattr(model, 'embedding_cache') and len(model.embedding_cache) >= len(protein_sequences_dict):
        print(f"Embeddings already cached for {len(model.embedding_cache)} proteins")
        return

    # Compute embeddings
    for protein_id, sequence in tqdm(protein_sequences_dict.items()):
        if protein_id not in model.embedding_cache:
            model.get_embeddings(sequence, use_cache=True)

    # Save cache
    print(f"Saving {len(model.embedding_cache)} embeddings to cache...")
    with open(cache_file, 'wb') as f:
        pickle.dump(model.embedding_cache, f)


# Helper to create visualization data
def create_viz_data(results_df, uniprot_to_gene, threshold=0.4):
    """Create visualization-friendly format from results dataframe"""
    # Extract protein pairs
    pairs = list(zip(results_df['acc1'], results_df['acc2']))

    # Create visualization-friendly format
    viz_df = pd.DataFrame({
        'protein1': [uniprot_to_gene.get(p[0], p[0]) for p in pairs],
        'protein2': [uniprot_to_gene.get(p[1], p[1]) for p in pairs],
        'prediction': results_df['prediction'].values,
        'method': results_df['method'].values
    })

    # Keep only high-confidence predictions for visualization
    viz_df = viz_df[viz_df['prediction'] >= threshold]
    return viz_df


# Compare male and female networks
def compare_networks(male_df, female_df, output_dir="DL-PPI outputs/comparison"):
    """Compare male and female networks"""
    print("Comparing male and female networks...")
    os.makedirs(output_dir, exist_ok=True)

    # 1. FIND SHARED AND SEX-SPECIFIC INTERACTIONS
    male_pairs = set(zip(male_df['protein1'], male_df['protein2']))
    female_pairs = set(zip(female_df['protein1'], female_df['protein2']))

    # Find shared and unique interactions
    shared_pairs = male_pairs.intersection(female_pairs)
    male_specific = male_pairs - shared_pairs
    female_specific = female_pairs - shared_pairs

    print(f"Found {len(shared_pairs)} shared interactions")
    print(f"Found {len(male_specific)} male-specific interactions")
    print(f"Found {len(female_specific)} female-specific interactions")

    # Create Venn diagram
    plt.figure(figsize=(10, 7), dpi=300)
    from matplotlib_venn import venn2
    venn = venn2([male_pairs, female_pairs], ('Male', 'Female'))
    venn.get_patch_by_id('10').set_color('blue')
    venn.get_patch_by_id('01').set_color('red')
    venn.get_patch_by_id('11').set_color('purple')
    plt.title('Overlap of Predicted Interactions by Sex', fontsize=16)
    plt.savefig(f"{output_dir}/interaction_overlap_venn.png", dpi=300, bbox_inches='tight')
    plt.close()

    # 2. ANALYZE SHARED INTERACTIONS FOR SCORE DIFFERENCES
    shared_interactions = []

    # For each shared interaction, find scores in both datasets
    for pair in shared_pairs:
        male_score = male_df[
            (male_df['protein1'] == pair[0]) & (male_df['protein2'] == pair[1]) |
            (male_df['protein1'] == pair[1]) & (male_df['protein2'] == pair[0])
            ]['prediction'].values[0]

        female_score = female_df[
            (female_df['protein1'] == pair[0]) & (female_df['protein2'] == pair[1]) |
            (female_df['protein1'] == pair[1]) & (female_df['protein2'] == pair[0])
            ]['prediction'].values[0]

        shared_interactions.append({
            'Protein1': pair[0],
            'Protein2': pair[1],
            'Male Score': male_score,
            'Female Score': female_score,
            'Difference': female_score - male_score
        })

    shared_df = pd.DataFrame(shared_interactions)

    # Get top differential interactions
    if not shared_df.empty:
        top_diff = shared_df.iloc[shared_df['Difference'].abs().argsort()[-15:][::-1]]

        plt.figure(figsize=(12, 8), dpi=300)

        # Prepare data for grouped bar chart
        index = range(len(top_diff))
        bar_width = 0.35

        # Create paired labels
        interaction_labels = [f"{row['Protein1']}-{row['Protein2']}" for _, row in top_diff.iterrows()]

        # Plot bars
        plt.barh([i for i in index], top_diff['Male Score'], bar_width,
                 label='Male', color='blue', alpha=0.7)
        plt.barh([i + bar_width for i in index], top_diff['Female Score'], bar_width,
                 label='Female', color='red', alpha=0.7)

        # Add details
        plt.yticks([i + bar_width / 2 for i in index], interaction_labels)
        plt.xlabel('Interaction Confidence Score')
        plt.title('Top Differential Protein Interactions by Sex', fontsize=16)
        plt.legend()
        plt.grid(axis='x', linestyle='--', alpha=0.5)
        plt.tight_layout()
        plt.savefig(f"{output_dir}/differential_interactions.png", dpi=300)
        plt.close()

    # 3. NETWORK STATISTICS COMPARISON
    print("Creating network statistics comparison...")
    # Create separate graphs for comparison statistics
    G_male = nx.Graph()
    for _, row in male_df.iterrows():
        G_male.add_edge(row['protein1'], row['protein2'], weight=row['prediction'])

    G_female = nx.Graph()
    for _, row in female_df.iterrows():
        G_female.add_edge(row['protein1'], row['protein2'], weight=row['prediction'])

    # Calculate network metrics
    male_stats = {
        'Nodes': G_male.number_of_nodes(),
        'Edges': G_male.number_of_edges(),
        'Density': nx.density(G_male),
        'Transitivity': nx.transitivity(G_male),
        'Avg Clustering': nx.average_clustering(G_male),
        'Avg Path Length': nx.average_shortest_path_length(G_male) if nx.is_connected(G_male) else 'N/A',
        'Diameter': nx.diameter(G_male) if nx.is_connected(G_male) else 'N/A',
    }

    female_stats = {
        'Nodes': G_female.number_of_nodes(),
        'Edges': G_female.number_of_edges(),
        'Density': nx.density(G_female),
        'Transitivity': nx.transitivity(G_female),
        'Avg Clustering': nx.average_clustering(G_female),
        'Avg Path Length': nx.average_shortest_path_length(G_female) if nx.is_connected(G_female) else 'N/A',
        'Diameter': nx.diameter(G_female) if nx.is_connected(G_female) else 'N/A',
    }

    # Create table of network statistics
    stats_df = pd.DataFrame({
        'Metric': list(male_stats.keys()),
        'Male': list(male_stats.values()),
        'Female': list(female_stats.values())
    })

    # Save network statistics
    stats_df.to_csv(f"{output_dir}/network_statistics.csv", index=False)
    print("Network statistics saved to network_statistics.csv")

    # 4. VISUALIZE SEX-SPECIFIC SUB-NETWORKS
    # Create sub-networks with top nodes by degree centrality
    create_subnetwork_viz(
        G_male,
        'Top Male-Specific Protein Interaction Network',
        f"{output_dir}/male_specific_network.png"
    )

    create_subnetwork_viz(
        G_female,
        'Top Female-Specific Protein Interaction Network',
        f"{output_dir}/female_specific_network.png"
    )

    # Return summary results
    return {
        'shared_count': len(shared_pairs),
        'male_specific_count': len(male_specific),
        'female_specific_count': len(female_specific),
        'network_stats': stats_df
    }

# 4. VISUALIZE SEX-SPECIFIC SUB-NETWORKS
# Create sub-networks with top nodes by degree centrality
def create_subnetwork_viz(G, title, filename, top_n=20):
    # Get top nodes by degree centrality
    degree_dict = dict(G.degree())
    top_nodes = sorted(degree_dict.items(), key=lambda x: x[1], reverse=True)[:top_n]
    top_node_names = [n[0] for n in top_nodes]

    # Create subgraph
    sub_G = G.subgraph(top_node_names)

    plt.figure(figsize=(10, 8), dpi=300)
    pos = nx.spring_layout(sub_G, seed=42)

    # Draw nodes
    nx.draw_networkx_nodes(sub_G, pos, node_size=400, node_color='skyblue', alpha=0.8)

    # Draw edges with weights as colors
    edges = sub_G.edges()
    weights = [sub_G[u][v]['weight'] for u, v in edges]
    nx.draw_networkx_edges(sub_G, pos, width=2, alpha=0.5, edge_color=weights, edge_cmap=plt.cm.YlOrRd)

    # Draw labels
    nx.draw_networkx_labels(sub_G, pos, font_size=10, font_family='sans-serif')

    plt.title(title, fontsize=16)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

    # Create visualizations
    create_subnetwork_viz(
        G_male,
        'Top Male-Specific Protein Interaction Network',
        f"{output_dir}/male_specific_network.png"
    )

    create_subnetwork_viz(
        G_female,
        'Top Female-Specific Protein Interaction Network',
        f"{output_dir}/female_specific_network.png"
    )

    # Return summary results
    return {
        'shared_count': len(shared_pairs),
        'male_specific_count': len(male_specific),
        'female_specific_count': len(female_specific),
        'network_stats': stats_df
    }


# Helper to visualize PPI network by pathways
def visualize_pathway_ppi_network(results_df, uniprot_to_gene, pathway_df,
                                 threshold=0.6, output_dir="DL-PPI outputs/networks"):
    """Visualize PPI network colored by pathway"""
    print("Visualizing pathway-based PPI network...")
    os.makedirs(output_dir, exist_ok=True)

    # Filter results by threshold
    high_conf_df = results_df[results_df['prediction'] >= threshold]
    print(f"Using {len(high_conf_df)} high-confidence interactions (score >= {threshold})")

    # Create network
    G = nx.Graph()

    # Add edges with weights
    for _, row in high_conf_df.iterrows():
        protein1 = uniprot_to_gene.get(row['acc1'], row['acc1'])
        protein2 = uniprot_to_gene.get(row['acc2'], row['acc2'])
        G.add_edge(protein1, protein2, weight=row['prediction'])

    # Get pathway info for nodes
    pathway_dict = {}
    for gene in G.nodes():
        pathways = pathway_df[pathway_df['gene'] == gene]['pathway'].unique().tolist()
        pathway_dict[gene] = pathways

    # Assign colors based on dominant pathway
    node_colors = []
    pathway_color_map = {}

    # Get unique pathways for color mapping
    all_pathways = pathway_df['pathway'].unique().tolist()
    cmap = plt.cm.tab20
    for i, pathway in enumerate(all_pathways):
        pathway_color_map[pathway] = cmap(i % 20)

    for node in G.nodes():
        pathways = pathway_dict.get(node, [])
        if pathways:
            # Use first pathway for color (could be enhanced to use most specific pathway)
            node_colors.append(pathway_color_map[pathways[0]])
        else:
            # Grey for nodes without pathway info
            node_colors.append((0.5, 0.5, 0.5, 0.5))

    # Create visualization
    plt.figure(figsize=(12, 10), dpi=300)

    # Position nodes using spring layout
    pos = nx.spring_layout(G, seed=42, k=0.3)

    # Draw network
    nx.draw_networkx_nodes(G, pos, node_size=200, node_color=node_colors, alpha=0.8)

    # Draw edges with weight-based width
    edges = G.edges()
    weights = [G[u][v]['weight'] * 3 for u, v in edges]  # Scale up for visibility
    nx.draw_networkx_edges(G, pos, width=weights, alpha=0.3, edge_color='gray')

    # Draw labels for high-degree nodes
    node_degrees = dict(G.degree())
    top_nodes = [n for n, d in sorted(node_degrees.items(), key=lambda x: x[1], reverse=True)[:30]]
    label_dict = {n: n for n in top_nodes}
    nx.draw_networkx_labels(G, pos, labels=label_dict, font_size=8, font_family='sans-serif')

    # Create legend for pathways
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color,
                                  label=pathway, markersize=10)
                      for pathway, color in pathway_color_map.items()]

    plt.axis('off')
    plt.title('Protein-Protein Interaction Network by Pathway', fontsize=16)
    plt.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5),
               ncol=1, fontsize=8)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/pathway_ppi_network.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Calculate and save network statistics
    stats = {
        'Nodes': G.number_of_nodes(),
        'Edges': G.number_of_edges(),
        'Average Degree': sum(dict(G.degree()).values()) / G.number_of_nodes(),
        'Density': nx.density(G),
        'Connected Components': nx.number_connected_components(G),
        'Average Clustering': nx.average_clustering(G),
    }

    # Find most central nodes
    degree_centrality = nx.degree_centrality(G)
    betweenness_centrality = nx.betweenness_centrality(G)
    eigenvector_centrality = nx.eigenvector_centrality(G, max_iter=1000)

    top_central_nodes = sorted(degree_centrality.items(), key=lambda x: x[1], reverse=True)[:20]

    # Save statistics
    with open(f"{output_dir}/network_stats.txt", 'w') as f:
        f.write("Network Statistics:\n")
        for k, v in stats.items():
            f.write(f"{k}: {v}\n")

        f.write("\nTop 20 Central Nodes (Degree Centrality):\n")
        for node, centrality in top_central_nodes:
            f.write(f"{node}: {centrality:.4f}\n")

    print(f"Network visualization saved to {output_dir}/pathway_ppi_network.png")
    print(f"Network statistics saved to {output_dir}/network_stats.txt")

    return G


# Main function to run the pipeline
def main():
    print("Starting memory-optimized DL-PPI pipeline for both sexes...")

    # Step 1: Load inflammatory pathway data
    pathway_df = load_inflammatory_pathways()

    # Step 2: Load sex-specific genes
    male_genes = load_sex_specific_genes(sex="M")
    female_genes = load_sex_specific_genes(sex="F")
    print(f"Loaded {len(male_genes)} male-specific and {len(female_genes)} female-specific genes")

    # Step 3: Load training data or create it if doesn't exist
    train_file = "DL-PPI outputs/data/known_interactions.csv"
    seq_file = "DL-PPI outputs/data/protein_sequences.pkl"

    if os.path.exists(train_file) and os.path.exists(seq_file):
        print("\nLoading existing training data...")
        training_df = pd.read_csv(train_file)
        with open(seq_file, 'rb') as f:
            protein_sequences_dict = pickle.load(f)
        print(f"Loaded {len(training_df)} training pairs, {training_df['interaction_label'].sum()} positive, {len(training_df) - training_df['interaction_label'].sum()} negative")
    else:
        training_df, protein_sequences_dict = create_complement_training_data()

    # Step 4: Create and train model or load existing
    model_file = "DL-PPI outputs/models/classifier_head.pt"
    embedding_cache_file = "DL-PPI outputs/models/embedding_cache.pkl"

    # Load model
    print("\nLoading facebook/esm2_t33_650M_UR50D model with memory optimizations...")
    dl_model = DeepLearningPPI(embedding_cache=None)

    # Load embedding cache if exists
    if os.path.exists(embedding_cache_file):
        dl_model.load_embedding_cache(embedding_cache_file)

    # Check if trained model exists
    if os.path.exists(model_file):
        print("\nFound existing trained model.")
        dl_model.classifier.load_state_dict(torch.load(model_file))
    else:
        print("\nTraining new model...")
        # Prepare training data
        train_pairs = list(zip(training_df['protein1'], training_df['protein2']))
        train_labels = training_df['interaction_label'].values

        # Train model
        dl_model = train_ppi_model(train_pairs, train_labels, protein_sequences_dict)

    # Step 5: Map genes to UniProt IDs - do this ONCE for all genes
    print("\nLoading cached protein data...")
    all_pathway_genes = pathway_df.gene.unique().tolist()
    all_genes = set(all_pathway_genes + male_genes + female_genes)
    gene_to_uniprot, protein_sequences_dict = load_cached_protein_data(list(all_genes))

    # Create reverse mapping for visualization
    uniprot_to_gene = {v: k for k, v in gene_to_uniprot.items()}

    # Add gene mapping to model for convenient access
    dl_model.uniprot_to_gene = uniprot_to_gene

    # Step 6: Precompute embeddings for all proteins
    precompute_all_embeddings(protein_sequences_dict, dl_model, embedding_cache_file)

    # Step 7: Generate protein pairs for prediction - Male specific
    print("\nGenerating protein pairs for male network...")
    male_genes_in_pathways = list(set(male_genes) & set(all_pathway_genes))
    male_pairs = generate_pathway_based_pairs(male_genes_in_pathways, pathway_df,
                                             protein_sequences_dict, gene_to_uniprot)

    # Step 8: Generate protein pairs for prediction - Female specific
    print("\nGenerating protein pairs for female network...")
    female_genes_in_pathways = list(set(female_genes) & set(all_pathway_genes))
    female_pairs = generate_pathway_based_pairs(female_genes_in_pathways, pathway_df,
                                              protein_sequences_dict, gene_to_uniprot)

    # Step 9: Predict interactions - Male
    print("\nPredicting male-specific interactions...")
    male_results = predict_ppi_deep_learning(male_pairs, protein_sequences_dict, dl_model,
                                           output_dir="DL-PPI outputs/results/male")

    # Step 10: Predict interactions - Female
    print("\nPredicting female-specific interactions...")
    female_results = predict_ppi_deep_learning(female_pairs, protein_sequences_dict, dl_model,
                                             output_dir="DL-PPI outputs/results/female")

    # Step 11: Create visualization data
    male_viz = create_viz_data(male_results, uniprot_to_gene, threshold=0.5)
    female_viz = create_viz_data(female_results, uniprot_to_gene, threshold=0.5)

    # Step 12: Compare networks
    comparison = compare_networks(male_viz, female_viz)

    # Step 13: Visualize pathway-based networks
    print("\nCreating pathway-based network visualizations...")
    male_network = visualize_pathway_ppi_network(male_results, uniprot_to_gene, pathway_df,
                                              output_dir="DL-PPI outputs/networks/male")
    female_network = visualize_pathway_ppi_network(female_results, uniprot_to_gene, pathway_df,
                                                output_dir="DL-PPI outputs/networks/female")

    print("\nPipeline completed successfully!")
    return {
        'dl_model': dl_model,
        'male_results': male_results,
        'female_results': female_results,
        'comparison': comparison
    }


if __name__ == "__main__":
    main()