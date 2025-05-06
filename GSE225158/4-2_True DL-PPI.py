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
import os
import requests

os.environ["OMP_NUM_THREADS"] = "1"  # Limit OpenMP threads
import torch

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
        device = torch.device("mps")
        print("Using MPS (Apple Silicon GPU) acceleration")
    else:
        device = torch.device("cpu")
        print("MPS not available, using CPU")

    # Limit memory usage
    os.environ["TOKENIZERS_PARALLELISM"] = "false"
    torch.set_num_threads(4)  # Prevent thread explosion

    return device


# Function to ensure directory exists
def ensure_path(path):
    os.makedirs(path, exist_ok=True)
    return path


# Enhanced DeepLearningPPI class with memory optimizations
class DeepLearningPPI:
    def __init__(self, model_name="facebook/esm2_t33_650M_UR50D", embedding_cache=None):
        """Initialize with a pre-trained protein language model"""
        print(f"Loading {model_name} model with memory optimizations...")
        # Get optimized device
        self.device = optimize_for_m1_mac()
        self.embedding_cache = embedding_cache or {}

        # Load tokenizer
        from transformers import AutoTokenizer, AutoModel, AutoConfig
        self.tokenizer = AutoTokenizer.from_pretrained(model_name)

        # Load model with memory optimizations
        config = AutoConfig.from_pretrained(model_name)
        with torch.no_grad():
            self.model = AutoModel.from_pretrained(
                model_name,
                config=config,
                torch_dtype=torch.float16,  # Use FP16 instead of FP32
                low_cpu_mem_usage=True
            ).to(self.device)

        # Define a smaller classifier head
        self.classifier = torch.nn.Sequential(
            torch.nn.Linear(2 * 1280, 256),  # Reduced intermediate size
            torch.nn.ReLU(),
            torch.nn.Dropout(0.1),
            torch.nn.Linear(256, 64),  # Reduced intermediate size
            torch.nn.ReLU(),
            torch.nn.Dropout(0.1),
            torch.nn.Linear(64, 1),
            torch.nn.Sigmoid()
        ).to(self.device)

        print(f"Model loaded on {self.device}")

    def get_embeddings(self, sequence, max_length=512, use_cache=True):
        """Get protein embedding with memory optimizations"""
        # Check cache first if enabled
        if use_cache and sequence in self.embedding_cache:
            return self.embedding_cache[sequence]

        # Truncate sequence if needed
        if len(sequence) > max_length:
            sequence = sequence[:max_length]

        # Use context manager to clear cache immediately
        with torch.no_grad():
            inputs = self.tokenizer(sequence, return_tensors="pt").to(self.device)
            outputs = self.model(**inputs)
            # Get embedding and immediately move to CPU and numpy
            embedding = outputs.last_hidden_state[:, 0, :].cpu().numpy()

        # Force garbage collection
        if hasattr(torch.cuda, 'empty_cache'):
            torch.cuda.empty_cache()
        gc.collect()

        # Cache the result if enabled
        if use_cache:
            self.embedding_cache[sequence] = embedding[0]

        return embedding[0]

    def predict_interaction(self, seq1, seq2):
        """Predict interaction between two proteins"""
        # Get embeddings
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

    def load_embedding_cache(self, cache_file):
        """Load embedding cache from file"""
        if os.path.exists(cache_file):
            try:
                with open(cache_file, 'rb') as f:
                    self.embedding_cache = pickle.load(f)
                print(f"Loaded {len(self.embedding_cache)} embeddings from cache")
            except Exception as e:
                print(f"Error loading embedding cache: {e}")
                self.embedding_cache = {}
        else:
            print("No embedding cache found, starting fresh")
            self.embedding_cache = {}


# Enhanced function to load inflammatory pathway data
def load_inflammatory_pathways(
        pathway_file="/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/KEGG outputs/kegg_inflammatory_pathways.csv"):
    """Load inflammatory pathway data and handle multiple pathways per gene"""
    try:
        # Load the pathway file
        path_df = pd.read_csv(pathway_file)

        # Get statistics on pathways
        pathways = path_df['pathway'].unique()
        pathway_counts = path_df.groupby('pathway').count()

        print(f"Loaded {len(path_df)} gene-pathway associations")
        print(f"Found {len(pathways)} unique pathways")
        print(f"Top 5 pathways by gene count:")
        for pathway, count in pathway_counts.sort_values('gene', ascending=False).head(5).iterrows():
            print(f"  - {pathway}: {count['gene']} genes")

        # Check for overlapping genes
        gene_pathway_counts = path_df.groupby('gene').count()
        overlapping_genes = gene_pathway_counts[gene_pathway_counts['pathway'] > 1]
        print(f"Found {len(overlapping_genes)} genes that belong to multiple pathways")

        return path_df
    except Exception as e:
        print(f"Error loading inflammatory pathway data: {e}")
        return pd.DataFrame(columns=["pathway", "gene"])


# Function to load sex-specific genes
def load_sex_specific_genes(sex="M"):
    """Load differentially expressed genes from DESeq2 results"""
    try:
        # Use the appropriate DESeq2 file based on sex
        if sex == "M":
            file_path = "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/DESeq2 outputs/deseq2_results_M_OUD_vs_M_None.csv"
        else:
            file_path = "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/DESeq2 outputs/deseq2_results_F_OUD_vs_F_None.csv"

        # Load the DESeq2 results
        deg_df = pd.read_csv(file_path)

        # Check column names and use appropriate ones
        padj_col = 'padj' if 'padj' in deg_df.columns else ('p.adj' if 'p.adj' in deg_df.columns else 'FDR')
        log2fc_col = 'log2FoldChange' if 'log2FoldChange' in deg_df.columns else 'log2FC'

        # Filter for significant DEGs (padj < 0.05 and |log2FC| > 1)
        filtered_df = deg_df[(deg_df[padj_col] < 0.05) & (abs(deg_df[log2fc_col]) > 1)]

        # Extract gene names
        gene_column = "gene" if "gene" in deg_df.columns else deg_df.columns[0]
        genes = filtered_df[gene_column].unique().tolist()

        print(f"Loaded {len(genes)} {sex} sex-specific genes")
        return genes
    except Exception as e:
        print(f"Error loading sex-specific genes: {e}")
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
        # Use a subset of inflammatory proteins as a seed
        base_url = "https://string-db.org/api"
        output_format = "tsv-no-header"
        method = "network"

        # Download PPI data - FIX: Use proper URL encoding for multiple identifiers
        proteins = ["C3", "C5", "CFB", "CFH", "CFD"]
        request_url = f"{base_url}/{output_format}/{method}?identifiers={'%0D'.join(proteins)}&species=9606&add_nodes=20&network_type=physical"

        print(f"Requesting URL: {request_url}")
        response = requests.get(request_url)

        if response.status_code != 200:
            print(f"API request failed with status code {response.status_code}")
            print(f"Response content: {response.text[:200]}...")
            raise Exception("Failed to retrieve data from STRING database")

        # Parse response and handle different column formats
        response_text = response.text.strip()
        if not response_text:
            raise Exception("Empty response from STRING database")

        # Debug the response
        print(f"First few lines of API response:")
        print('\n'.join(response_text.split('\n')[:5]))

        ppi_df = pd.read_csv(StringIO(response_text), sep='\t', header=None)

        # Handle different column formats
        if len(ppi_df.columns) == 3:
            ppi_df.columns = ['protein1', 'protein2', 'score']
        elif len(ppi_df.columns) == 4:
            ppi_df.columns = ['protein1', 'protein2', 'score', 'additional']
            # Keep only the columns we need
            ppi_df = ppi_df[['protein1', 'protein2', 'score']]
        else:
            print(f"Warning: Unexpected column count in response: {len(ppi_df.columns)}")
            print(f"Columns: {ppi_df.columns}")
            # Create minimal required columns
            if len(ppi_df.columns) >= 2:
                ppi_df = ppi_df.iloc[:, :2]
                ppi_df.columns = ['protein1', 'protein2']
                ppi_df['score'] = 900  # Default high score
            else:
                raise Exception(f"Insufficient columns in API response: {len(ppi_df.columns)}")

        # Filter for high confidence interactions
        if 'score' in ppi_df.columns:
            ppi_df = ppi_df[ppi_df['score'] > 700]

        print(f"Downloaded {len(ppi_df)} high-confidence PPIs")

        # Create backup data if no PPIs found
        if len(ppi_df) == 0:
            print("No PPIs found, creating sample data for testing")
            sample_proteins = ['P01024', 'P01031', 'P00751', 'P08603', 'P00746']
            sample_pairs = []
            for i in range(len(sample_proteins)):
                for j in range(i + 1, len(sample_proteins)):
                    sample_pairs.append((sample_proteins[i], sample_proteins[j], 900))
            ppi_df = pd.DataFrame(sample_pairs, columns=['protein1', 'protein2', 'score'])

        # Create positive and negative examples
        proteins = set(ppi_df['protein1'].tolist() + ppi_df['protein2'].tolist())

        # Generate negative examples
        neg_pairs = []
        num_pos = len(ppi_df)

        while len(neg_pairs) < num_pos:
            p1 = np.random.choice(list(proteins))
            p2 = np.random.choice(list(proteins))

            if p1 != p2 and not ppi_df[(ppi_df['protein1'] == p1) & (ppi_df['protein2'] == p2)].any().any():
                if not ppi_df[(ppi_df['protein1'] == p2) & (ppi_df['protein2'] == p1)].any().any():
                    neg_pairs.append((p1, p2, 0))  # Include label 0 for negative examples

        # Create balanced dataset
        pos_data = [(row['protein1'], row['protein2'], 1) for _, row in ppi_df.iterrows()]

        # Combine and shuffle
        all_data = pos_data + neg_pairs
        np.random.shuffle(all_data)

        # Create DataFrame with explicit column names
        train_df = pd.DataFrame(all_data, columns=['protein1', 'protein2', 'interaction_label'])
        print(f"Created {len(train_df)} training pairs, {len(pos_data)} positive, {len(neg_pairs)} negative")

        # Get protein sequences
        protein_sequences = {}
        for protein_id in tqdm(proteins, desc="Downloading protein sequences"):
            try:
                # Check if it's a STRING ID (has 9606.ENSP format)
                if protein_id.startswith('9606.ENSP'):
                    # Extract the Ensembl ID part
                    ensembl_id = protein_id.split('.')[1]

                    # Use mygene to convert Ensembl ID to UniProt
                    from mygene import MyGeneInfo
                    mg = MyGeneInfo()
                    result = mg.query(ensembl_id, scopes='ensembl.protein', fields='uniprot', species='human')

                    if result['hits'] and 'uniprot' in result['hits'][0]:
                        # Get the UniProt ID
                        if 'Swiss-Prot' in result['hits'][0]['uniprot']:
                            uniprot_id = result['hits'][0]['uniprot']['Swiss-Prot']
                        elif isinstance(result['hits'][0]['uniprot'], list) and len(result['hits'][0]['uniprot']) > 0:
                            uniprot_id = result['hits'][0]['uniprot'][0]
                        else:
                            print(f"No Swiss-Prot entry found for {protein_id}")
                            continue

                        # Download sequence using proper UniProt ID
                        response = requests.get(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta")
                        if response.status_code == 200:
                            lines = response.text.strip().split('\n')
                            sequence = ''.join(lines[1:])
                            protein_sequences[protein_id] = sequence
                        else:
                            print(
                                f"Could not retrieve sequence for mapped ID {uniprot_id}, status: {response.status_code}")
                    else:
                        print(f"Could not map {protein_id} to UniProt")
                else:
                    # Try original method for non-Ensembl IDs (assumed to be UniProt)
                    response = requests.get(f"https://rest.uniprot.org/uniprotkb/{protein_id}.fasta")
                    if response.status_code == 200:
                        lines = response.text.strip().split('\n')
                        sequence = ''.join(lines[1:])
                        protein_sequences[protein_id] = sequence
                    else:
                        print(f"Could not retrieve sequence for {protein_id}, status: {response.status_code}")
            except Exception as e:
                print(f"Error processing {protein_id}: {e}")

        print(f"Downloaded {len(protein_sequences)} protein sequences")

        # Save data
        os.makedirs(f"{output_dir}/data", exist_ok=True)
        train_df.to_csv(f"{output_dir}/data/known_interactions.csv", index=False)

        # Save protein sequences as FASTA
        with open(f"{output_dir}/data/protein_sequences.fasta", 'w') as f:
            for protein_id, sequence in protein_sequences.items():
                f.write(f">{protein_id}\n{sequence}\n")

        return train_df, protein_sequences

    except Exception as e:
        print(f"Error creating training data: {e}")
        # Create fallback minimal dataset so the pipeline doesn't fail
        fallback_proteins = ['P01024', 'P01031', 'P00751', 'P08603', 'P00746']
        fallback_sequences = {
            'P01024': 'MGPTSGPSLLLLLLTHLPLALGSPMYSIITPNILRLESEETMVLEAHDAQGDVPVTVTVHDFPGKKLVLSSEKTVLTPATNHMGNVTFTIPANREFKSEKGRNKFVTVQATFGTQVVEKVVLVSLQSGYLFIQTDKTIYTPGSTVLYRIFTVNHKLLPVGRTVMVNIENPEGIPVKQDSLSSQNQLGVLPLSWDIPELVNMGQWKIRAYYENSPQQVFSTEFEVKEYVLPSFEVIVEPTEKFYYIYNEKGLEVTITARFLYGKKVEGTAFVIFGIQDGEQRISLPESLKRIPIEDGSGEVVLSRKVLLDGVQNPRAEDLVGKSLYVSATVILHSGSDMVQAERSGIPIVTSPYQIHFTKTPKYFKPGMPFDLMVFVTNPDGSPAYRVPVAVQGEDTVQSLTQGDGVAKLSINTHPSQKPLSITVRTKKQELSEAEQATRTMQALPYSTVGNSNNYLHLSVLRTELRPGETLNVNFLLRMDRAHEAKIRYYTYLIMNKGRLLKAGRQVREPGQDLVVLPLSITTDFIPSFRLVAYYTLIGASGQREVVADSVWVDVKDSCVGSLVVKSGQSEDRQPVPGQQMTLKIEGDHGARVVLVAVDKGVFVLNKKNKLTQSKIWDVVEKADIGCTPGSGKDYAGVFSDAGLTFTSSSGQQTAQRAELQCPQPAAR',
            'P01031': 'MGLLGILCFLIFLGKTWGQEQTYVISAPKIFRVGASENIVIQVTGYTEAEISVPANVQFTCPPPQDTSLNTKGICESPGVIPTKTKPLSVFKCDTPPEVNVHVPGLYVIGTSVILQEGHLAFLTPGKPYQAVVFVGLKKDVPVYEIGVPAESTGTWNLGSSLWSHEYDVGHIMFSFGPPTHDVNFPLEITCEQTTKQPWDCPKIKCSRCMPTQGQTVTNTLHYLKGVKVLLENLDVTIETDYAEPEFIIKRPDIAVSTWDPHARTTFKSFIKVLRNSQYSELDFGTVKEVVPEGISVPIEGTVKILMSGSVVSKRKTLYSSVPKEKGQALISPFHEFDIHGDDLAEELDLNISDISGFLGAMYVINFKGDIQVTFTNVVVESGSLVKSNLSSSVVYTKVTLVDYKERLASIQLPELTKMPFTLQVVTACRHGVCSEGDLASIASLKVQEAAISLTSIHPDIIQLIKTNLSLFEHSIAHCFIENVLMTNLELSDRGGSNFDCWTSVRGYQEGKVFYLLEDPTISLSRNSLEISKLKADMVPLNTDFRQFVTLPVLYQLSFDQLALNYIIPQMTNCILYGPDGGSTVMASVSVFLILGAVLLFIMIACVLRRQKRPPDREAPSLKDEFNSLVGFFPVFLSAVFFAIGLVLGVALRHRKRQQKINEQEIAAIGARPQPSLQPPSTTKTPPTKTPSPKPVLVSLPSLLKTDTEQHQVPSQDQVRIQQESPVEGASDTLSPEEEQASGGLAEALAALVVVVLVVVILIALIAYFRQKITPGKCAVNFIQSRQRTSIHIGHANEISVSIHAEDSDSSSSGDHIASSVSVIAREDDELPGDRLYVQLSSPGFTKEKVQINVTQVIPVTLQPLPSFRILMGSILNQLTFNPSELQLVFWEYDPSKMAIPILSIDPLVIHQDWLNMGVWLHAFSHQDGLHSLLVPVSLEDIQGWFNMLGISIAFLLLYFLVRLLLCRRPPRSQGSPQLLKTPDILPTRQPLSSISNRDNRKHSRPSQTSTKSPKLSHNSHMGVVPTASSLNVETQASQGTFSHFFPHQE',
            'P00751': 'MRLLAKIICLMLWAICVAQDCVPKEGSDQNPRGQPNIPGSYPEGPGSYPEGPQIPGSYPEGPGSYPEGPGSYPEGTLSYGPGGGIYGSYPEGPGSYAEGPGSYAEGPGSYAEGPGSYPEGPQIYAEGPGSYAEGPRSYAEGPRIYAEGPGSYAEGPTIYAEGPGSYAEGTGIYAEGPGSYAEGRQSYSEGPQIYPEGPGIYAEGPQIYSEGPGIYPEGPGSYAEGPGIYAEGPTIYAEGPGIYAEGPGIYAEGPGIYAEGPTIYAEGPGIYAEGPGIYPEGPGIYAEGPGIYSEGPRSYAEGPGIYPEGPGIYPEGPGIYAEGPGIYPEGPGIYPEGPGSYDEGPGIYPEGPGAYPEGPGIYPEGPGIYPEGPGIYPEGPGAYTEGPGSYAEGPGSYPEGPGSYPEGPGSYPEGTLSYGPGG',
            'P08603': 'MRLLAKIICLMADSCSPTRLLLYYISYTVLLFALSFIHFSTTQCGYPCEIRRGENAVFQEHGDKVELPAAKPNPTTGPATVVEEVCANRHVYYGAVAGFGGGVKKQVLFDLVPVRSRRVPLHEIYNDGAYHVGQALLTPTEVTYTCNEGYSLIGNHTAACVLIGEINGDYLELDSNTTHMMDFGPCIVPDRKHAHGQSQVLDVQRFGLCRASNMPTTVGRFSSAPNILRLRYHNEEWAEIIKDAYWAKILATPGSVGEHYFTAVVEKLAQLAGKYTYLVLTLRHAKLTGGSPIYTGSSMVFEPPKQGYSILALLDTNLKHVVTAKGKEGAWDADNCCVKKCGKFSACQDSGNDCGDFSDEDDCEPCQHRKCNNPPPQNGGSPCPGIYNEKLNYPGSVKTYSCDPRRLVNGDANQTCTDGEWANLPSFCQREDSFISCKLSPCSQPPQIEHGTINSSRSSQESYAHGTKLSYTCEGGFRISEENETTCYMGKWSSPPQCEGLPCKSPPEISHGVVAHMASPSTSPDTGTTAERKKREVEKTVDCCSLTLGLLIGILVLVFCLYFKKSRKKPSYKPRPTDPLSIPVVTPTLPTSLPTVVETTASTEKCQGAGYVTQSLPVPTAPKLCPPPPQIPNSHNMTTTLNYRDGEKVSVLCQENYLIQEGEEITCKDGRWQSIPLCVEKIPCSQPPQIEHGTINSSRMSLLWDQQKLPSCPQDITVNSSQDWGSVKSYECKPGYRLQGADTTICLENGWSPEPPQCIRKCKS',
            'P00746': 'MGATRSVWGLLGAAAVLWAVSGDSCHWDEDCVEIFTNGKWNDRACGDKSAKCRLIKDPGCEYRNECKEIVEEAEERMEEVDQCVVAHRGGPGSQVWLGRCCPGGECPRKKFTICVPESRYWELKDGCVGFGGGDLKCQGCGNPWKPLPECKGDVKCVPIELRRRQSVLTEIHQGVCLEVDECGAAPGQCLDGSALWEFSCHTSCPVKSQPCGFDACWPEHTWNCSVSSCGGGVQCRRTCTDNNKCEGRWRCDGEACGQQCKTKACDGECCAYVHRGKCVPGFGGGWQCGCHNKCIDKFWCDGDPDCKDGSDEVSCPRVPTCKPPCGHRRCQCHAACQVGWACSQDECLSGLCMGDGICQCDSKLCEGDKFCQRGQCICKDSGRWGCNCERRCEQERLAVRYEACDTDDECSGLEKCKDGVCQCPRRRELRYLGPCKAPGIQCGKGTCVQRSSRTDSGRCICGAGTVGYHCEHSDCPSGYHCECRAGYTGAFCDQDLFCPEAVRCEPVQCPPPKIAAAVLPVAVALAVLIALLVLSAIWYSWKRHSRTANNLDIVELYRIADL'
        }

        fallback_pairs = []
        for i in range(len(fallback_proteins)):
            for j in range(i + 1, len(fallback_proteins)):
                fallback_pairs.append((fallback_proteins[i], fallback_proteins[j], 1 if i % 2 == 0 else 0))

        train_df = pd.DataFrame(fallback_pairs, columns=['protein1', 'protein2', 'interaction_label'])
        print(f"Created {len(train_df)} fallback training pairs")

        return train_df, fallback_sequences


# Helper to generate pathway-based pairs
def generate_pathway_based_pairs(gene_list, pathway_df, protein_sequences_dict, gene_to_uniprot):
    """Generate protein pairs within pathways"""
    # Create pathway groups
    pathway_groups = pathway_df.groupby('pathway')['gene'].apply(list).to_dict()

    pairs = []
    for pathway, genes in pathway_groups.items():
        # Find intersection with our gene list
        pathway_genes = [g for g in genes if g in gene_list]

        # Get UniProt IDs where available
        uniprot_ids = [gene_to_uniprot.get(g) for g in pathway_genes]
        uniprot_ids = [u for u in uniprot_ids if u and u in protein_sequences_dict]

        # Generate all pairwise combinations within pathway
        for i in range(len(uniprot_ids)):
            for j in range(i + 1, len(uniprot_ids)):
                pairs.append((uniprot_ids[i], uniprot_ids[j]))

    print(f"Generated {len(pairs)} pathway-based protein pairs")

    # After getting UniProt IDs
    uniprot_ids = [gene_to_uniprot.get(g) for g in pathway_genes]
    uniprot_ids_raw = [u for u in uniprot_ids if u]
    uniprot_ids = [u for u in uniprot_ids_raw if u in protein_sequences_dict]
    print(
        f"Pathway {pathway}: Found {len(pathway_genes)} genes, {len(uniprot_ids_raw)} with UniProt IDs, {len(uniprot_ids)} with sequences")
    return pairs


# Helper to get pathway information for a gene pair
def get_pathway_info(gene1, gene2, pathway_df):
    """Get pathway information for a gene pair"""
    gene1_pathways = pathway_df[pathway_df['gene'] == gene1]['pathway'].unique().tolist()
    gene2_pathways = pathway_df[pathway_df['gene'] == gene2]['pathway'].unique().tolist()

    # Find common pathways
    common_pathways = list(set(gene1_pathways) & set(gene2_pathways))

    return {
        'shared_pathways': common_pathways,
        'gene1_pathways': gene1_pathways,
        'gene2_pathways': gene2_pathways,
        'is_same_pathway': len(common_pathways) > 0
    }

# Dataset class for protein-protein interactions
class PPIDataset(Dataset):
    def __init__(self, pairs, labels, sequences_dict, model):
        self.pairs = pairs
        self.labels = labels
        self.sequences_dict = sequences_dict
        self.model = model

    def __len__(self):
        return len(self.pairs)

    def __getitem__(self, idx):
        p1, p2 = self.pairs[idx]
        label = self.labels[idx]

        # Get protein sequences with better placeholder
        seq1 = self.sequences_dict.get(p1, "X" * 10)
        seq2 = self.sequences_dict.get(p2, "X" * 10)

        # Get embeddings with error handling
        try:
            emb1 = self.model.get_embeddings(seq1)
            emb2 = self.model.get_embeddings(seq2)
        except Exception as e:
            print(f"Error generating embeddings for {p1}/{p2}: {e}")
            # Return zero embeddings as fallback
            emb1 = np.zeros(1280)
            emb2 = np.zeros(1280)

        return {
            "emb1": torch.tensor(emb1, dtype=torch.float32),
            "emb2": torch.tensor(emb2, dtype=torch.float32),
            "label": torch.tensor(label, dtype=torch.float32)
        }

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
        range(len(dataset)), test_size=0.2, random_state=RANDOM_SEED,
        stratify=train_labels
    )

    train_sampler = torch.utils.data.SubsetRandomSampler(train_indices)
    val_sampler = torch.utils.data.SubsetRandomSampler(val_indices)

    train_loader = DataLoader(
        dataset, batch_size=batch_size, sampler=train_sampler,
        num_workers=0, pin_memory=False
    )

    val_loader = DataLoader(
        dataset, batch_size=batch_size, sampler=val_sampler,
        num_workers=0, pin_memory=False
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
        # Training
        model.classifier.train()
        train_loss = 0.0
        train_pbar = tqdm(train_loader, desc=f"Epoch {epoch + 1}/{epochs} [Train]")

        for batch in train_pbar:
            # Zero gradients
            optimizer.zero_grad()

            # Get embeddings and labels
            emb1 = batch["emb1"].to(model.device)
            emb2 = batch["emb2"].to(model.device)
            labels = batch["label"].unsqueeze(1).to(model.device)

            # Forward pass
            concat_emb = torch.cat([emb1, emb2], dim=1)
            outputs = model.classifier(concat_emb)

            # Calculate loss
            loss = criterion(outputs, labels)

            # Backward pass
            loss.backward()

            # Gradient clipping
            torch.nn.utils.clip_grad_norm_(model.classifier.parameters(), max_norm=1.0)

            # Update weights
            optimizer.step()

            # Update loss
            train_loss += loss.item() * labels.size(0)
            train_pbar.set_postfix({"batch_loss": f"{loss.item():.4f}"})

            # Force garbage collection
            gc.collect()
            if hasattr(torch.cuda, 'empty_cache'):
                torch.cuda.empty_cache()

        # Calculate average training loss
        train_loss = train_loss / len(train_indices)

        # Validation
        model.classifier.eval()
        val_loss = 0.0
        correct = 0
        total = 0

        val_pbar = tqdm(val_loader, desc=f"Epoch {epoch + 1}/{epochs} [Val]")
        with torch.no_grad():
            for batch in val_pbar:
                # Get embeddings and labels
                emb1 = batch["emb1"].to(model.device)
                emb2 = batch["emb2"].to(model.device)
                labels = batch["label"].unsqueeze(1).to(model.device)

                # Forward pass
                concat_emb = torch.cat([emb1, emb2], dim=1)
                outputs = model.classifier(concat_emb)

                # Calculate loss
                loss = criterion(outputs, labels)
                val_loss += loss.item() * labels.size(0)

                # Calculate accuracy
                predicted = (outputs > 0.5).float()
                correct += (predicted == labels).sum().item()
                total += labels.size(0)

                val_pbar.set_postfix({"batch_loss": f"{loss.item():.4f}"})

                # Force garbage collection
                gc.collect()
                if hasattr(torch.cuda, 'empty_cache'):
                    torch.cuda.empty_cache()

        # Calculate average validation loss and accuracy
        val_loss = val_loss / len(val_indices)
        val_acc = correct / total

        print(
            f"Epoch {epoch + 1}/{epochs}: Train Loss: {train_loss:.4f}, Val Loss: {val_loss:.4f}, Val Acc: {val_acc:.4f}")

        # Learning rate scheduler
        if lr_scheduler:
            scheduler.step(val_loss)

        # Early stopping
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            patience_counter = 0

            # Save best model
            os.makedirs(output_dir, exist_ok=True)
            torch.save(model.classifier.state_dict(), f"{output_dir}/classifier_head.pt")
            print(f"Saved best model with val_loss: {val_loss:.4f}")
        else:
            patience_counter += 1

        if patience_counter >= patience:
            print(f"Early stopping after {epoch + 1} epochs")
            break

    # Load best model
    model.classifier.load_state_dict(torch.load(f"{output_dir}/classifier_head.pt"))
    return model


# Function to predict PPIs
def predict_ppi_deep_learning(protein_pairs, protein_sequences_dict, model=None,
                              batch_size=4, model_path="DL-PPI outputs/models/classifier_head.pt",
                              output_dir="DL-PPI outputs/results"):
    """Predict protein-protein interactions using deep learning"""
    print(f"Predicting interactions for {len(protein_pairs)} protein pairs...")
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
    import os
    os.environ["OMP_NUM_THREADS"] = "1"  # Limit OpenMP threads
    import torch
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
            device = torch.device("mps")
            print("Using MPS (Apple Silicon GPU) acceleration")
        else:
            device = torch.device("cpu")
            print("MPS not available, using CPU")

        # Limit memory usage
        os.environ["TOKENIZERS_PARALLELISM"] = "false"
        torch.set_num_threads(4)  # Prevent thread explosion

        return device

    # Function to ensure directory exists
    def ensure_path(path):
        os.makedirs(path, exist_ok=True)
        return path

    # Enhanced DeepLearningPPI class with memory optimizations
    class DeepLearningPPI:
        def __init__(self, model_name="facebook/esm2_t33_650M_UR50D", embedding_cache=None):
            """Initialize with a pre-trained protein language model"""
            print(f"Loading {model_name} model with memory optimizations...")
            # Get optimized device
            self.device = optimize_for_m1_mac()
            self.embedding_cache = embedding_cache or {}

            # Load tokenizer
            from transformers import AutoTokenizer, AutoModel, AutoConfig
            self.tokenizer = AutoTokenizer.from_pretrained(model_name)

            # Load model with memory optimizations
            config = AutoConfig.from_pretrained(model_name)
            with torch.no_grad():
                self.model = AutoModel.from_pretrained(
                    model_name,
                    config=config,
                    torch_dtype=torch.float16,  # Use FP16 instead of FP32
                    low_cpu_mem_usage=True
                ).to(self.device)

            # Define a smaller classifier head
            self.classifier = torch.nn.Sequential(
                torch.nn.Linear(2 * 1280, 256),  # Reduced intermediate size
                torch.nn.ReLU(),
                torch.nn.Dropout(0.1),
                torch.nn.Linear(256, 64),  # Reduced intermediate size
                torch.nn.ReLU(),
                torch.nn.Dropout(0.1),
                torch.nn.Linear(64, 1),
                torch.nn.Sigmoid()
            ).to(self.device)

            print(f"Model loaded on {self.device}")

        def get_embeddings(self, sequence, max_length=512, use_cache=True):
            """Get protein embedding with memory optimizations"""
            # Check cache first if enabled
            if use_cache and sequence in self.embedding_cache:
                return self.embedding_cache[sequence]

            # Truncate sequence if needed
            if len(sequence) > max_length:
                sequence = sequence[:max_length]

            # Use context manager to clear cache immediately
            with torch.no_grad():
                inputs = self.tokenizer(sequence, return_tensors="pt").to(self.device)
                outputs = self.model(**inputs)
                # Get embedding and immediately move to CPU and numpy
                embedding = outputs.last_hidden_state[:, 0, :].cpu().numpy()

            # Force garbage collection
            if hasattr(torch.cuda, 'empty_cache'):
                torch.cuda.empty_cache()
            gc.collect()

            # Cache the result if enabled
            if use_cache:
                self.embedding_cache[sequence] = embedding[0]

            return embedding[0]

        def predict_interaction(self, seq1, seq2):
            """Predict interaction between two proteins"""
            # Get embeddings
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

        def load_embedding_cache(self, cache_file):
            """Load embedding cache from file"""
            if os.path.exists(cache_file):
                try:
                    with open(cache_file, 'rb') as f:
                        self.embedding_cache = pickle.load(f)
                    print(f"Loaded {len(self.embedding_cache)} embeddings from cache")
                except Exception as e:
                    print(f"Error loading embedding cache: {e}")
                    self.embedding_cache = {}
            else:
                print("No embedding cache found, starting fresh")
                self.embedding_cache = {}

    # Enhanced function to load inflammatory pathway data
    def load_inflammatory_pathways(
            pathway_file="/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/KEGG outputs/kegg_inflammatory_pathways.csv"):
        """Load inflammatory pathway data and handle multiple pathways per gene"""
        try:
            # Load the pathway file
            path_df = pd.read_csv(pathway_file)

            # Get statistics on pathways
            pathways = path_df['pathway'].unique()
            pathway_counts = path_df.groupby('pathway').count()

            print(f"Loaded {len(path_df)} gene-pathway associations")
            print(f"Found {len(pathways)} unique pathways")
            print(f"Top 5 pathways by gene count:")
            for pathway, count in pathway_counts.sort_values('gene', ascending=False).head(5).iterrows():
                print(f"  - {pathway}: {count['gene']} genes")

            # Check for overlapping genes
            gene_pathway_counts = path_df.groupby('gene').count()
            overlapping_genes = gene_pathway_counts[gene_pathway_counts['pathway'] > 1]
            print(f"Found {len(overlapping_genes)} genes that belong to multiple pathways")

            return path_df
        except Exception as e:
            print(f"Error loading inflammatory pathway data: {e}")
            return pd.DataFrame(columns=["pathway", "gene"])

    # Function to load sex-specific genes
    def load_sex_specific_genes(sex="M"):
        """Load differentially expressed genes from DESeq2 results"""
        try:
            # Use the appropriate DESeq2 file based on sex
            if sex == "M":
                file_path = "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/DESeq2 outputs/deseq2_results_M_OUD_vs_M_None.csv"
            else:
                file_path = "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/DESeq2 outputs/deseq2_results_F_OUD_vs_F_None.csv"

            # Load the DESeq2 results
            deg_df = pd.read_csv(file_path)

            # Check column names and use appropriate ones
            padj_col = 'padj' if 'padj' in deg_df.columns else ('p.adj' if 'p.adj' in deg_df.columns else 'FDR')
            log2fc_col = 'log2FoldChange' if 'log2FoldChange' in deg_df.columns else 'log2FC'

            # Filter for significant DEGs (padj < 0.05 and |log2FC| > 1)
            filtered_df = deg_df[(deg_df[padj_col] < 0.05) & (abs(deg_df[log2fc_col]) > 1)]

            # Extract gene names
            gene_column = "gene" if "gene" in deg_df.columns else deg_df.columns[0]
            genes = filtered_df[gene_column].unique().tolist()

            print(f"Loaded {len(genes)} {sex} sex-specific genes")
            return genes
        except Exception as e:
            print(f"Error loading sex-specific genes: {e}")
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
            # Use a subset of inflammatory proteins as a seed
            base_url = "https://string-db.org/api"
            output_format = "tsv-no-header"
            method = "network"

            # Download PPI data - FIX: Use proper URL encoding for multiple identifiers
            proteins = ["C3", "C5", "CFB", "CFH", "CFD"]
            request_url = f"{base_url}/{output_format}/{method}?identifiers={'%0D'.join(proteins)}&species=9606&add_nodes=20&network_type=physical"

            print(f"Requesting URL: {request_url}")
            response = requests.get(request_url)

            if response.status_code != 200:
                print(f"API request failed with status code {response.status_code}")
                print(f"Response content: {response.text[:200]}...")
                raise Exception("Failed to retrieve data from STRING database")

            # Parse response and handle different column formats
            response_text = response.text.strip()
            if not response_text:
                raise Exception("Empty response from STRING database")

            # Debug the response
            print(f"First few lines of API response:")
            print('\n'.join(response_text.split('\n')[:5]))

            ppi_df = pd.read_csv(StringIO(response_text), sep='\t', header=None)

            # Handle different column formats
            if len(ppi_df.columns) == 3:
                ppi_df.columns = ['protein1', 'protein2', 'score']
            elif len(ppi_df.columns) == 4:
                ppi_df.columns = ['protein1', 'protein2', 'score', 'additional']
                # Keep only the columns we need
                ppi_df = ppi_df[['protein1', 'protein2', 'score']]
            else:
                print(f"Warning: Unexpected column count in response: {len(ppi_df.columns)}")
                print(f"Columns: {ppi_df.columns}")
                # Create minimal required columns
                if len(ppi_df.columns) >= 2:
                    ppi_df = ppi_df.iloc[:, :2]
                    ppi_df.columns = ['protein1', 'protein2']
                    ppi_df['score'] = 900  # Default high score
                else:
                    raise Exception(f"Insufficient columns in API response: {len(ppi_df.columns)}")

            # Filter for high confidence interactions
            if 'score' in ppi_df.columns:
                ppi_df = ppi_df[ppi_df['score'] > 700]

            print(f"Downloaded {len(ppi_df)} high-confidence PPIs")

            # Create backup data if no PPIs found
            if len(ppi_df) == 0:
                print("No PPIs found, creating sample data for testing")
                sample_proteins = ['P01024', 'P01031', 'P00751', 'P08603', 'P00746']
                sample_pairs = []
                for i in range(len(sample_proteins)):
                    for j in range(i + 1, len(sample_proteins)):
                        sample_pairs.append((sample_proteins[i], sample_proteins[j], 900))
                ppi_df = pd.DataFrame(sample_pairs, columns=['protein1', 'protein2', 'score'])

            # Create positive and negative examples
            proteins = set(ppi_df['protein1'].tolist() + ppi_df['protein2'].tolist())

            # Generate negative examples
            neg_pairs = []
            num_pos = len(ppi_df)

            while len(neg_pairs) < num_pos:
                p1 = np.random.choice(list(proteins))
                p2 = np.random.choice(list(proteins))

                if p1 != p2 and not ppi_df[(ppi_df['protein1'] == p1) & (ppi_df['protein2'] == p2)].any().any():
                    if not ppi_df[(ppi_df['protein1'] == p2) & (ppi_df['protein2'] == p1)].any().any():
                        neg_pairs.append((p1, p2, 0))  # Include label 0 for negative examples

            # Create balanced dataset
            pos_data = [(row['protein1'], row['protein2'], 1) for _, row in ppi_df.iterrows()]

            # Combine and shuffle
            all_data = pos_data + neg_pairs
            np.random.shuffle(all_data)

            # Create DataFrame with explicit column names
            train_df = pd.DataFrame(all_data, columns=['protein1', 'protein2', 'interaction_label'])
            print(f"Created {len(train_df)} training pairs, {len(pos_data)} positive, {len(neg_pairs)} negative")

            # Get protein sequences
            protein_sequences = {}
            for protein_id in tqdm(proteins, desc="Downloading protein sequences"):
                try:
                    # Check if it's a STRING ID (has 9606.ENSP format)
                    if protein_id.startswith('9606.ENSP'):
                        # Extract the Ensembl ID part
                        ensembl_id = protein_id.split('.')[1]

                        # Use mygene to convert Ensembl ID to UniProt
                        from mygene import MyGeneInfo
                        mg = MyGeneInfo()
                        result = mg.query(ensembl_id, scopes='ensembl.protein', fields='uniprot', species='human')

                        if result['hits'] and 'uniprot' in result['hits'][0]:
                            # Get the UniProt ID
                            if 'Swiss-Prot' in result['hits'][0]['uniprot']:
                                uniprot_id = result['hits'][0]['uniprot']['Swiss-Prot']
                            elif isinstance(result['hits'][0]['uniprot'], list) and len(
                                    result['hits'][0]['uniprot']) > 0:
                                uniprot_id = result['hits'][0]['uniprot'][0]
                            else:
                                print(f"No Swiss-Prot entry found for {protein_id}")
                                continue

                            # Download sequence using proper UniProt ID
                            response = requests.get(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta")
                            if response.status_code == 200:
                                lines = response.text.strip().split('\n')
                                sequence = ''.join(lines[1:])
                                protein_sequences[protein_id] = sequence
                            else:
                                print(
                                    f"Could not retrieve sequence for mapped ID {uniprot_id}, status: {response.status_code}")
                        else:
                            print(f"Could not map {protein_id} to UniProt")
                    else:
                        # Try original method for non-Ensembl IDs (assumed to be UniProt)
                        response = requests.get(f"https://rest.uniprot.org/uniprotkb/{protein_id}.fasta")
                        if response.status_code == 200:
                            lines = response.text.strip().split('\n')
                            sequence = ''.join(lines[1:])
                            protein_sequences[protein_id] = sequence
                        else:
                            print(f"Could not retrieve sequence for {protein_id}, status: {response.status_code}")
                except Exception as e:
                    print(f"Error processing {protein_id}: {e}")

            print(f"Downloaded {len(protein_sequences)} protein sequences")

            # Save data
            os.makedirs(f"{output_dir}/data", exist_ok=True)
            train_df.to_csv(f"{output_dir}/data/known_interactions.csv", index=False)

            # Save protein sequences as FASTA
            with open(f"{output_dir}/data/protein_sequences.fasta", 'w') as f:
                for protein_id, sequence in protein_sequences.items():
                    f.write(f">{protein_id}\n{sequence}\n")

            return train_df, protein_sequences

        except Exception as e:
            print(f"Error creating training data: {e}")
            # Create fallback minimal dataset so the pipeline doesn't fail
            fallback_proteins = ['P01024', 'P01031', 'P00751', 'P08603', 'P00746']
            fallback_sequences = {
                'P01024': 'MGPTSGPSLLLLLLTHLPLALGSPMYSIITPNILRLESEETMVLEAHDAQGDVPVTVTVHDFPGKKLVLSSEKTVLTPATNHMGNVTFTIPANREFKSEKGRNKFVTVQATFGTQVVEKVVLVSLQSGYLFIQTDKTIYTPGSTVLYRIFTVNHKLLPVGRTVMVNIENPEGIPVKQDSLSSQNQLGVLPLSWDIPELVNMGQWKIRAYYENSPQQVFSTEFEVKEYVLPSFEVIVEPTEKFYYIYNEKGLEVTITARFLYGKKVEGTAFVIFGIQDGEQRISLPESLKRIPIEDGSGEVVLSRKVLLDGVQNPRAEDLVGKSLYVSATVILHSGSDMVQAERSGIPIVTSPYQIHFTKTPKYFKPGMPFDLMVFVTNPDGSPAYRVPVAVQGEDTVQSLTQGDGVAKLSINTHPSQKPLSITVRTKKQELSEAEQATRTMQALPYSTVGNSNNYLHLSVLRTELRPGETLNVNFLLRMDRAHEAKIRYYTYLIMNKGRLLKAGRQVREPGQDLVVLPLSITTDFIPSFRLVAYYTLIGASGQREVVADSVWVDVKDSCVGSLVVKSGQSEDRQPVPGQQMTLKIEGDHGARVVLVAVDKGVFVLNKKNKLTQSKIWDVVEKADIGCTPGSGKDYAGVFSDAGLTFTSSSGQQTAQRAELQCPQPAAR',
                'P01031': 'MGLLGILCFLIFLGKTWGQEQTYVISAPKIFRVGASENIVIQVTGYTEAEISVPANVQFTCPPPQDTSLNTKGICESPGVIPTKTKPLSVFKCDTPPEVNVHVPGLYVIGTSVILQEGHLAFLTPGKPYQAVVFVGLKKDVPVYEIGVPAESTGTWNLGSSLWSHEYDVGHIMFSFGPPTHDVNFPLEITCEQTTKQPWDCPKIKCSRCMPTQGQTVTNTLHYLKGVKVLLENLDVTIETDYAEPEFIIKRPDIAVSTWDPHARTTFKSFIKVLRNSQYSELDFGTVKEVVPEGISVPIEGTVKILMSGSVVSKRKTLYSSVPKEKGQALISPFHEFDIHGDDLAEELDLNISDISGFLGAMYVINFKGDIQVTFTNVVVESGSLVKSNLSSSVVYTKVTLVDYKERLASIQLPELTKMPFTLQVVTACRHGVCSEGDLASIASLKVQEAAISLTSIHPDIIQLIKTNLSLFEHSIAHCFIENVLMTNLELSDRGGSNFDCWTSVRGYQEGKVFYLLEDPTISLSRNSLEISKLKADMVPLNTDFRQFVTLPVLYQLSFDQLALNYIIPQMTNCILYGPDGGSTVMASVSVFLILGAVLLFIMIACVLRRQKRPPDREAPSLKDEFNSLVGFFPVFLSAVFFAIGLVLGVALRHRKRQQKINEQEIAAIGARPQPSLQPPSTTKTPPTKTPSPKPVLVSLPSLLKTDTEQHQVPSQDQVRIQQESPVEGASDTLSPEEEQASGGLAEALAALVVVVLVVVILIALIAYFRQKITPGKCAVNFIQSRQRTSIHIGHANEISVSIHAEDSDSSSSGDHIASSVSVIAREDDELPGDRLYVQLSSPGFTKEKVQINVTQVIPVTLQPLPSFRILMGSILNQLTFNPSELQLVFWEYDPSKMAIPILSIDPLVIHQDWLNMGVWLHAFSHQDGLHSLLVPVSLEDIQGWFNMLGISIAFLLLYFLVRLLLCRRPPRSQGSPQLLKTPDILPTRQPLSSISNRDNRKHSRPSQTSTKSPKLSHNSHMGVVPTASSLNVETQASQGTFSHFFPHQE',
                'P00751': 'MRLLAKIICLMLWAICVAQDCVPKEGSDQNPRGQPNIPGSYPEGPGSYPEGPQIPGSYPEGPGSYPEGPGSYPEGTLSYGPGGGIYGSYPEGPGSYAEGPGSYAEGPGSYAEGPGSYPEGPQIYAEGPGSYAEGPRSYAEGPRIYAEGPGSYAEGPTIYAEGPGSYAEGTGIYAEGPGSYAEGRQSYSEGPQIYPEGPGIYAEGPQIYSEGPGIYPEGPGSYAEGPGIYAEGPTIYAEGPGIYAEGPGIYAEGPGIYAEGPTIYAEGPGIYAEGPGIYPEGPGIYAEGPGIYSEGPRSYAEGPGIYPEGPGIYPEGPGIYAEGPGIYPEGPGIYPEGPGSYDEGPGIYPEGPGAYPEGPGIYPEGPGIYPEGPGIYPEGPGAYTEGPGSYAEGPGSYPEGPGSYPEGPGSYPEGTLSYGPGG',
                'P08603': 'MRLLAKIICLMADSCSPTRLLLYYISYTVLLFALSFIHFSTTQCGYPCEIRRGENAVFQEHGDKVELPAAKPNPTTGPATVVEEVCANRHVYYGAVAGFGGGVKKQVLFDLVPVRSRRVPLHEIYNDGAYHVGQALLTPTEVTYTCNEGYSLIGNHTAACVLIGEINGDYLELDSNTTHMMDFGPCIVPDRKHAHGQSQVLDVQRFGLCRASNMPTTVGRFSSAPNILRLRYHNEEWAEIIKDAYWAKILATPGSVGEHYFTAVVEKLAQLAGKYTYLVLTLRHAKLTGGSPIYTGSSMVFEPPKQGYSILALLDTNLKHVVTAKGKEGAWDADNCCVKKCGKFSACQDSGNDCGDFSDEDDCEPCQHRKCNNPPPQNGGSPCPGIYNEKLNYPGSVKTYSCDPRRLVNGDANQTCTDGEWANLPSFCQREDSFISCKLSPCSQPPQIEHGTINSSRSSQESYAHGTKLSYTCEGGFRISEENETTCYMGKWSSPPQCEGLPCKSPPEISHGVVAHMASPSTSPDTGTTAERKKREVEKTVDCCSLTLGLLIGILVLVFCLYFKKSRKKPSYKPRPTDPLSIPVVTPTLPTSLPTVVETTASTEKCQGAGYVTQSLPVPTAPKLCPPPPQIPNSHNMTTTLNYRDGEKVSVLCQENYLIQEGEEITCKDGRWQSIPLCVEKIPCSQPPQIEHGTINSSRMSLLWDQQKLPSCPQDITVNSSQDWGSVKSYECKPGYRLQGADTTICLENGWSPEPPQCIRKCKS',
                'P00746': 'MGATRSVWGLLGAAAVLWAVSGDSCHWDEDCVEIFTNGKWNDRACGDKSAKCRLIKDPGCEYRNECKEIVEEAEERMEEVDQCVVAHRGGPGSQVWLGRCCPGGECPRKKFTICVPESRYWELKDGCVGFGGGDLKCQGCGNPWKPLPECKGDVKCVPIELRRRQSVLTEIHQGVCLEVDECGAAPGQCLDGSALWEFSCHTSCPVKSQPCGFDACWPEHTWNCSVSSCGGGVQCRRTCTDNNKCEGRWRCDGEACGQQCKTKACDGECCAYVHRGKCVPGFGGGWQCGCHNKCIDKFWCDGDPDCKDGSDEVSCPRVPTCKPPCGHRRCQCHAACQVGWACSQDECLSGLCMGDGICQCDSKLCEGDKFCQRGQCICKDSGRWGCNCERRCEQERLAVRYEACDTDDECSGLEKCKDGVCQCPRRRELRYLGPCKAPGIQCGKGTCVQRSSRTDSGRCICGAGTVGYHCEHSDCPSGYHCECRAGYTGAFCDQDLFCPEAVRCEPVQCPPPKIAAAVLPVAVALAVLIALLVLSAIWYSWKRHSRTANNLDIVELYRIADL'
            }

            fallback_pairs = []
            for i in range(len(fallback_proteins)):
                for j in range(i + 1, len(fallback_proteins)):
                    fallback_pairs.append((fallback_proteins[i], fallback_proteins[j], 1 if i % 2 == 0 else 0))

            train_df = pd.DataFrame(fallback_pairs, columns=['protein1', 'protein2', 'interaction_label'])
            print(f"Created {len(train_df)} fallback training pairs")

            return train_df, fallback_sequences

    # Helper to generate pathway-based pairs
    def generate_pathway_based_pairs(gene_list, pathway_df, protein_sequences_dict, gene_to_uniprot):
        """Generate protein pairs within pathways"""
        # Create pathway groups
        pathway_groups = pathway_df.groupby('pathway')['gene'].apply(list).to_dict()

        pairs = []
        for pathway, genes in pathway_groups.items():
            # Find intersection with our gene list
            pathway_genes = [g for g in genes if g in gene_list]

            # Get UniProt IDs where available
            uniprot_ids = [gene_to_uniprot.get(g) for g in pathway_genes]
            uniprot_ids = [u for u in uniprot_ids if u and u in protein_sequences_dict]

            # Generate all pairwise combinations within pathway
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
            'shared_pathways': common_pathways,
            'gene1_pathways': gene1_pathways,
            'gene2_pathways': gene2_pathways,
            'is_same_pathway': len(common_pathways) > 0
        }

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
            range(len(dataset)), test_size=0.2, random_state=RANDOM_SEED,
            stratify=train_labels
        )

        train_sampler = torch.utils.data.SubsetRandomSampler(train_indices)
        val_sampler = torch.utils.data.SubsetRandomSampler(val_indices)

        train_loader = DataLoader(
            dataset, batch_size=batch_size, sampler=train_sampler,
            num_workers=0, pin_memory=False
        )

        val_loader = DataLoader(
            dataset, batch_size=batch_size, sampler=val_sampler,
            num_workers=0, pin_memory=False
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
            # Training
            model.classifier.train()
            train_loss = 0.0
            train_pbar = tqdm(train_loader, desc=f"Epoch {epoch + 1}/{epochs} [Train]")

            for batch in train_pbar:
                # Zero gradients
                optimizer.zero_grad()

                # Get embeddings and labels
                emb1 = batch["emb1"].to(model.device)
                emb2 = batch["emb2"].to(model.device)
                labels = batch["label"].unsqueeze(1).to(model.device)

                # Forward pass
                concat_emb = torch.cat([emb1, emb2], dim=1)
                outputs = model.classifier(concat_emb)

                # Calculate loss
                loss = criterion(outputs, labels)

                # Backward pass
                loss.backward()

                # Gradient clipping
                torch.nn.utils.clip_grad_norm_(model.classifier.parameters(), max_norm=1.0)

                # Update weights
                optimizer.step()

                # Update loss
                train_loss += loss.item() * labels.size(0)
                train_pbar.set_postfix({"batch_loss": f"{loss.item():.4f}"})

                # Force garbage collection
                gc.collect()
                if hasattr(torch.cuda, 'empty_cache'):
                    torch.cuda.empty_cache()

            # Calculate average training loss
            train_loss = train_loss / len(train_indices)

            # Validation
            model.classifier.eval()
            val_loss = 0.0
            correct = 0
            total = 0

            val_pbar = tqdm(val_loader, desc=f"Epoch {epoch + 1}/{epochs} [Val]")
            with torch.no_grad():
                for batch in val_pbar:
                    # Get embeddings and labels
                    emb1 = batch["emb1"].to(model.device)
                    emb2 = batch["emb2"].to(model.device)
                    labels = batch["label"].unsqueeze(1).to(model.device)

                    # Forward pass
                    concat_emb = torch.cat([emb1, emb2], dim=1)
                    outputs = model.classifier(concat_emb)

                    # Calculate loss
                    loss = criterion(outputs, labels)
                    val_loss += loss.item() * labels.size(0)

                    # Calculate accuracy
                    predicted = (outputs > 0.5).float()
                    correct += (predicted == labels).sum().item()
                    total += labels.size(0)

                    val_pbar.set_postfix({"batch_loss": f"{loss.item():.4f}"})

                    # Force garbage collection
                    gc.collect()
                    if hasattr(torch.cuda, 'empty_cache'):
                        torch.cuda.empty_cache()

            # Calculate average validation loss and accuracy
            val_loss = val_loss / len(val_indices)
            val_acc = correct / total

            print(
                f"Epoch {epoch + 1}/{epochs}: Train Loss: {train_loss:.4f}, Val Loss: {val_loss:.4f}, Val Acc: {val_acc:.4f}")

            # Learning rate scheduler
            if lr_scheduler:
                scheduler.step(val_loss)

            # Early stopping
            if val_loss < best_val_loss:
                best_val_loss = val_loss
                patience_counter = 0

                # Save best model
                os.makedirs(output_dir, exist_ok=True)
                torch.save(model.classifier.state_dict(), f"{output_dir}/classifier_head.pt")
                print(f"Saved best model with val_loss: {val_loss:.4f}")
            else:
                patience_counter += 1

            if patience_counter >= patience:
                print(f"Early stopping after {epoch + 1} epochs")
                break

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

                for p1, p2 in batch:
                    # Get sequences
                    seq1 = protein_sequences_dict.get(p1, "")
                    seq2 = protein_sequences_dict.get(p2, "")

                    if not seq1 or not seq2:
                        continue

                    # Predict
                    result = model.predict_interaction(seq1, seq2)
                    result["acc1"] = p1
                    result["acc2"] = p2

                    # Add gene names if available
                    if hasattr(model, 'uniprot_to_gene') and p1 in model.uniprot_to_gene:
                        result["gene1"] = model.uniprot_to_gene[p1]
                    if hasattr(model, 'uniprot_to_gene') and p2 in model.uniprot_to_gene:
                        result["gene2"] = model.uniprot_to_gene[p2]

                    batch_predictions.append(result)

                predictions.extend(batch_predictions)

                # Force garbage collection
                gc.collect()
                if hasattr(torch.cuda, 'empty_cache'):
                    torch.cuda.empty_cache()

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
                model.get_embeddings(sequence)

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
                     color='blue', alpha=0.7, label='Male')
            plt.barh([i + bar_width for i in index], top_diff['Female Score'], bar_width,
                     color='red', alpha=0.7, label='Female')

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
        metrics = {
            'Metric': ['Nodes', 'Edges', 'Density', 'Avg. Degree', 'Connected Components'],
            'Male': [
                len(G_male.nodes()),
                len(G_male.edges()),
                nx.density(G_male),
                sum(dict(G_male.degree()).values()) / len(G_male.nodes()),
                nx.number_connected_components(G_male)
            ],
            'Female': [
                len(G_female.nodes()),
                len(G_female.edges()),
                nx.density(G_female),
                sum(dict(G_female.degree()).values()) / len(G_female.nodes()),
                nx.number_connected_components(G_female)
            ]
        }

        # Create a nice comparison table/figure
        metrics_df = pd.DataFrame(metrics)
        metrics_df['Difference %'] = (metrics_df['Female'] - metrics_df['Male']) / metrics_df['Male'] * 100

        # Create a visual table
        plt.figure(figsize=(10, 5), dpi=300)
        plt.axis('off')
        table = plt.table(
            cellText=metrics_df.values,
            colLabels=metrics_df.columns,
            cellLoc='center',
            loc='center',
            cellColours=[['#f2f2f2'] * 4] * 5
        )
        table.auto_set_font_size(False)
        table.set_fontsize(12)
        table.scale(1.2, 1.8)
        plt.title('Network Property Comparison: Male vs Female OUD', fontsize=16, pad=20)
        plt.savefig(f"{output_dir}/network_metrics_comparison.png", dpi=300, bbox_inches='tight')
        plt.close()

        # 4. STATISTICAL COMPARISON OF SCORE DISTRIBUTIONS
        print("Creating statistical comparison of score distributions...")

        # A. Kolmogorov-Smirnov test to statistically compare distributions
        ks_stat, p_value = stats.ks_2samp(male_df['prediction'], female_df['prediction'])
        print(f"KS statistic: {ks_stat}, p-value: {p_value}")

        # B. Quantile-Quantile plot
        plt.figure(figsize=(8, 8), dpi=300)
        male_scores = sorted(male_df['prediction'])
        female_scores = sorted(female_df['prediction'])

        # Make equal length if needed
        min_len = min(len(male_scores), len(female_scores))
        male_quantiles = np.linspace(0, 1, min_len)
        female_quantiles = np.linspace(0, 1, min_len)

        male_sample = np.quantile(male_scores, male_quantiles)
        female_sample = np.quantile(female_scores, female_quantiles)

        plt.scatter(male_sample, female_sample, alpha=0.5)
        plt.plot([0.4, 0.6], [0.4, 0.6], 'r--')  # diagonal line
        plt.xlabel('Male Prediction Scores')
        plt.ylabel('Female Prediction Scores')
        plt.title('Q-Q Plot: Male vs Female Prediction Scores')
        plt.grid(alpha=0.3)

        # Add KS test results to plot
        plt.figtext(0.15, 0.15,
                    f"Kolmogorov-Smirnov test:\nStatistic: {ks_stat:.4f}\np-value: {p_value:.4f}",
                    bbox=dict(facecolor='white', alpha=0.8),
                    fontsize=10)

        plt.savefig(f"{output_dir}/qq_plot.png", dpi=300, bbox_inches='tight')
        plt.close()

        # C. Enhanced histogram/density plot
        plt.figure(figsize=(10, 6), dpi=300)

        # Create histogram/KDE of prediction scores
        sns.histplot(male_df['prediction'], kde=True, color='blue', alpha=0.5, label='Male', bins=20)
        sns.histplot(female_df['prediction'], kde=True, color='red', alpha=0.5, label='Female', bins=20)

        # Add mean lines
        plt.axvline(male_df['prediction'].mean(), color='blue', linestyle='--',
                    label=f'Male mean: {male_df["prediction"].mean():.4f}')
        plt.axvline(female_df['prediction'].mean(), color='red', linestyle='--',
                    label=f'Female mean: {female_df["prediction"].mean():.4f}')

        # Add KS test results to plot
        plt.text(0.42, plt.gca().get_ylim()[1] * 0.9,
                 f"KS test: statistic={ks_stat:.4f}, p-value={p_value:.4f}",
                 bbox=dict(facecolor='white', alpha=0.8))

        plt.xlabel('Prediction Score')
        plt.ylabel('Frequency')
        plt.title('Distribution of Prediction Scores by Sex', fontsize=16)
        plt.legend()
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(f"{output_dir}/prediction_score_distribution.png", dpi=300)
        plt.close()

        print(f"All visualizations saved to {output_dir}")
        return metrics_df

    # Main function to run both sexes in a single pass
    def main():
        print("Starting memory-optimized DL-PPI pipeline for both sexes...")
        run_params = {
            "model": "facebook/esm2_t33_650M_UR50D",
            "batch_size_train": 8,
            "batch_size_predict": 4,
            "epochs": 20,
            "patience": 3,
            "lr": 1e-4,
            "weight_decay": 1e-5,
            "embedding_dim": 1280,
            "timestamp": time.strftime("%Y%m%d-%H%M%S")
        }

        # Create output directories
        out_dir = "DL-PPI outputs"
        os.makedirs(out_dir, exist_ok=True)
        os.makedirs(f"{out_dir}/data", exist_ok=True)
        os.makedirs(f"{out_dir}/models", exist_ok=True)
        os.makedirs(f"{out_dir}/results", exist_ok=True)

        # Create sex-specific output directories
        male_dir = ensure_path(f"{out_dir}/results/male")
        female_dir = ensure_path(f"{out_dir}/results/female")
        comp_dir = ensure_path(f"{out_dir}/results/comparison")

        # Load inflammatory pathway data
        pathway_file = "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/KEGG outputs/kegg_inflammatory_pathways.csv"
        pathway_df = load_inflammatory_pathways(pathway_file)

        if pathway_df.empty:
            print("WARNING: Pathway data unavailable!")
            return

        # Load or create training data
        if not os.path.exists(f"{out_dir}/data/known_interactions.csv"):
            print("\nCreating training data for inflammatory pathways...")
            training_df, protein_sequences_dict = create_complement_training_data(out_dir)
        else:
            print("\nLoading existing training data...")
            training_df = pd.read_csv(f"{out_dir}/data/known_interactions.csv")

            # Load protein sequences from FASTA
            protein_sequences_dict = {}
            with open(f"{out_dir}/data/protein_sequences.fasta", 'r') as f:
                current_id = None
                current_seq = ""
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_id:
                            protein_sequences_dict[current_id] = current_seq
                        current_id = line[1:]
                        current_seq = ""
                    else:
                        current_seq += line
                if current_id:
                    protein_sequences_dict[current_id] = current_seq

        # Create training pairs and labels
        train_pairs = list(zip(training_df['protein1'], training_df['protein2']))
        train_labels = training_df['interaction_label'].tolist()

        print(
            f"Loaded {len(train_pairs)} training pairs, {sum(train_labels)} positive, {len(train_labels) - sum(train_labels)} negative")

        # Initialize deep learning model - ONCE for both sexes
        embedding_cache_file = f"{out_dir}/data/embedding_cache.pkl"
        dl_model = DeepLearningPPI()

        # Load embedding cache if available
        dl_model.load_embedding_cache(embedding_cache_file)

        # Check if model exists or needs training
        if os.path.exists(f"{out_dir}/models/classifier_head.pt"):
            print("\nFound existing trained model.")
            # Load model for inference
            dl_model.classifier.load_state_dict(torch.load(f"{out_dir}/models/classifier_head.pt"))
        else:
            print("\nTraining new model with memory optimizations...")
            # Train model with smaller batch size for M1 Mac
            dl_model = train_ppi_model(train_pairs, train_labels, protein_sequences_dict,
                                       batch_size=8,  # Smaller batch size for M1 Mac
                                       epochs=20,
                                       patience=3,
                                       output_dir=f"{out_dir}/models")

        # Precompute embeddings for all proteins to avoid recomputing
        precompute_all_embeddings(protein_sequences_dict, dl_model, embedding_cache_file)

        # Extract all genes from inflammatory pathways
        all_pathway_genes = pathway_df.gene.unique().tolist()
        # Filter out entries that aren't gene symbols
        all_pathway_genes = [g for g in all_pathway_genes if isinstance(g, str) and not g.startswith('[')]
        print(f"Loaded {len(all_pathway_genes)} genes from inflammatory pathways")

        # Load sex-specific DEGs
        male_genes = load_sex_specific_genes("M")
        female_genes = load_sex_specific_genes("F")

        print(f"Loaded {len(male_genes)} male-specific and {len(female_genes)} female-specific genes")

        # Before generating pairs, map all genes and download sequences
        all_genes = list(set(pathway_df['gene'].tolist() + male_genes + female_genes))
        gene_to_uniprot, protein_sequences_dict = map_genes_and_download_sequences(all_genes, protein_sequences_dict)

        # More comprehensive filtering before mapping
        print("\nMapping genes to UniProt IDs...")
        from mygene import MyGeneInfo
        import re
        mg = MyGeneInfo()

        # Get all genes
        all_genes = set(all_pathway_genes + male_genes + female_genes)

        # More aggressive filtering for non-coding RNAs
        filtered_genes = []
        for gene in all_genes:
            # Skip any gene that matches non-coding RNA patterns
            if (isinstance(gene, str) and
                    not re.match(r'^(AC|AL|AP|RP|LINC|LOC\d+|.*-AS\d+|.*-IT\d+|.*\.\d+)', gene) and
                    not '-AS' in gene and
                    not gene.startswith(('MIR', 'SNORD', 'SNORA', 'TCONS_', 'XLOC_')) and
                    len(gene) > 1):  # Skip single-character entries
                filtered_genes.append(gene)

        print(f"Filtered out {len(all_genes) - len(filtered_genes)} non-protein coding entries")
        print(f"Proceeding with {len(filtered_genes)} potential protein-coding genes")

        # Query only the filtered genes
        gene_info = mg.querymany(filtered_genes, scopes='symbol', fields='uniprot', species='human', returnall=True)

        # Process results
        gene_to_uniprot = {}
        for entry in gene_info['out']:
            if 'uniprot' in entry and 'Swiss-Prot' in entry['uniprot']:
                gene = entry.get('query', '')
                uniprot_id = entry['uniprot']['Swiss-Prot']
                if isinstance(uniprot_id, list):
                    uniprot_id = uniprot_id[0]  # Take the first one
                gene_to_uniprot[gene] = uniprot_id

        # Create reverse mapping for visualization
        uniprot_to_gene = {v: k for k, v in gene_to_uniprot.items()}

        # Add gene mapping to model for convenient access
        dl_model.uniprot_to_gene = uniprot_to_gene

        # Process male data
        male_results = None
        print("\n========== PROCESSING MALE DATA ==========")
        if male_genes:
            # Generate pathway-based pairs for males
            male_pairs = generate_pathway_based_pairs(male_genes, pathway_df, protein_sequences_dict, gene_to_uniprot)
            male_test_pairs = [(p[0], p[1]) for p in male_pairs]

            if male_test_pairs:
                male_results = predict_ppi_deep_learning(
                    male_test_pairs,
                    protein_sequences_dict,
                    model=dl_model,
                    batch_size=run_params["batch_size_predict"],
                    output_dir=male_dir
                )

                # Add metadata to results
                if not male_results.empty:
                    # Add gene names
                    male_results['gene1'] = male_results['acc1'].map(uniprot_to_gene)
                    male_results['gene2'] = male_results['acc2'].map(uniprot_to_gene)

                    # Add pathway information
                    def get_pathway_info_for_row(row):
                        gene1 = row['gene1'] if pd.notna(row['gene1']) else ""
                        gene2 = row['gene2'] if pd.notna(row['gene2']) else ""
                        return get_pathway_info(gene1, gene2, pathway_df)

                    male_results['pathway_info'] = male_results.apply(get_pathway_info_for_row, axis=1)
                    male_results.to_csv(f"{male_dir}/dl_ppi_predictions_with_metadata.csv", index=False)

                    # Save visualization-friendly format
                    male_viz_df = create_viz_data(male_results, uniprot_to_gene, threshold=0.4)
                    male_viz_df.to_csv(f"{out_dir}/results/dl_ppi_viz_data_male_0.4_genes.csv", index=False)

        # Force garbage collection to free memory before processing female data
        gc.collect()
        if hasattr(torch.cuda, 'empty_cache'):
            torch.cuda.empty_cache()

        # Process female data
        female_results = None
        print("\n========== PROCESSING FEMALE DATA ==========")
        if female_genes:
            # Generate pathway-based pairs for females
            female_pairs = generate_pathway_based_pairs(female_genes, pathway_df, protein_sequences_dict,
                                                        gene_to_uniprot)
            female_test_pairs = [(p[0], p[1]) for p in female_pairs]

            if female_test_pairs:
                female_results = predict_ppi_deep_learning(
                    female_test_pairs,
                    protein_sequences_dict,
                    model=dl_model,
                    batch_size=run_params["batch_size_predict"],
                    output_dir=female_dir
                )

                # Add metadata to results
                if not female_results.empty:
                    # Add gene names
                    female_results['gene1'] = female_results['acc1'].map(uniprot_to_gene)
                    female_results['gene2'] = female_results['acc2'].map(uniprot_to_gene)

                    # Add pathway information
                    def get_pathway_info_for_row(row):
                        gene1 = row['gene1'] if pd.notna(row['gene1']) else ""
                        gene2 = row['gene2'] if pd.notna(row['gene2']) else ""
                        return get_pathway_info(gene1, gene2, pathway_df)

                    female_results['pathway_info'] = female_results.apply(get_pathway_info_for_row, axis=1)
                    female_results.to_csv(f"{female_dir}/dl_ppi_predictions_with_metadata.csv", index=False)

                    # Save visualization-friendly format
                    female_viz_df = create_viz_data(female_results, uniprot_to_gene, threshold=0.4)
                    female_viz_df.to_csv(f"{out_dir}/results/dl_ppi_viz_data_female_0.4_genes.csv", index=False)

        # Compare networks if both have results
        if male_results is not None and female_results is not None and not male_results.empty and not female_results.empty:
            print("\n========== COMPARING MALE AND FEMALE NETWORKS ==========")
            # Convert UniProt IDs to gene symbols for visualization
            male_viz_df = pd.DataFrame({
                'protein1': male_results['gene1'],
                'protein2': male_results['gene2'],
                'prediction': male_results['prediction'],
                'method': male_results['method']
            })

            female_viz_df = pd.DataFrame({
                'protein1': female_results['gene1'],
                'protein2': female_results['gene2'],
                'prediction': female_results['prediction'],
                'method': female_results['method']
            })

            # Compare the networks
            network_metrics = compare_networks(male_viz_df, female_viz_df, output_dir=comp_dir)

            # Save metrics
            network_metrics.to_csv(f"{comp_dir}/network_metrics.csv", index=False)

        # Save embedding cache before finishing
        print("\nSaving embedding cache...")
        with open(embedding_cache_file, 'wb') as f:
            pickle.dump(dl_model.embedding_cache, f)

        # Save run parameters for reproducibility
        with open(f"{out_dir}/run_parameters.json", 'w') as f:
            json.dump(run_params, f, indent=2)

        print("\nDL-PPI prediction pipeline completed for both sexes.")
        print(f"Results saved to {out_dir}")
        print(f"Run parameters saved to {out_dir}/run_parameters.json")

    if __name__ == "__main__":
        main()
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

            for p1, p2 in batch:
                # Get sequences
                seq1 = protein_sequences_dict.get(p1, "")
                seq2 = protein_sequences_dict.get(p2, "")

                if not seq1 or not seq2:
                    continue

                # Predict
                result = model.predict_interaction(seq1, seq2)
                result["acc1"] = p1
                result["acc2"] = p2

                # Add gene names if available
                if hasattr(model, 'uniprot_to_gene') and p1 in model.uniprot_to_gene:
                    result["gene1"] = model.uniprot_to_gene[p1]
                if hasattr(model, 'uniprot_to_gene') and p2 in model.uniprot_to_gene:
                    result["gene2"] = model.uniprot_to_gene[p2]

                batch_predictions.append(result)

            predictions.extend(batch_predictions)

            # Force garbage collection
            gc.collect()
            if hasattr(torch.cuda, 'empty_cache'):
                torch.cuda.empty_cache()

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
            model.get_embeddings(sequence)

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
                 color='blue', alpha=0.7, label='Male')
        plt.barh([i + bar_width for i in index], top_diff['Female Score'], bar_width,
                 color='red', alpha=0.7, label='Female')

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
    metrics = {
        'Metric': ['Nodes', 'Edges', 'Density', 'Avg. Degree', 'Connected Components'],
        'Male': [
            len(G_male.nodes()),
            len(G_male.edges()),
            nx.density(G_male),
            sum(dict(G_male.degree()).values()) / len(G_male.nodes()),
            nx.number_connected_components(G_male)
        ],
        'Female': [
            len(G_female.nodes()),
            len(G_female.edges()),
            nx.density(G_female),
            sum(dict(G_female.degree()).values()) / len(G_female.nodes()),
            nx.number_connected_components(G_female)
        ]
    }

    # Create a nice comparison table/figure
    metrics_df = pd.DataFrame(metrics)
    metrics_df['Difference %'] = (metrics_df['Female'] - metrics_df['Male']) / metrics_df['Male'] * 100

    # Create a visual table
    plt.figure(figsize=(10, 5), dpi=300)
    plt.axis('off')
    table = plt.table(
        cellText=metrics_df.values,
        colLabels=metrics_df.columns,
        cellLoc='center',
        loc='center',
        cellColours=[['#f2f2f2'] * 4] * 5
    )
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1.2, 1.8)
    plt.title('Network Property Comparison: Male vs Female OUD', fontsize=16, pad=20)
    plt.savefig(f"{output_dir}/network_metrics_comparison.png", dpi=300, bbox_inches='tight')
    plt.close()

    # 4. STATISTICAL COMPARISON OF SCORE DISTRIBUTIONS
    print("Creating statistical comparison of score distributions...")

    # A. Kolmogorov-Smirnov test to statistically compare distributions
    ks_stat, p_value = stats.ks_2samp(male_df['prediction'], female_df['prediction'])
    print(f"KS statistic: {ks_stat}, p-value: {p_value}")

    # B. Quantile-Quantile plot
    plt.figure(figsize=(8, 8), dpi=300)
    male_scores = sorted(male_df['prediction'])
    female_scores = sorted(female_df['prediction'])

    # Make equal length if needed
    min_len = min(len(male_scores), len(female_scores))
    male_quantiles = np.linspace(0, 1, min_len)
    female_quantiles = np.linspace(0, 1, min_len)

    male_sample = np.quantile(male_scores, male_quantiles)
    female_sample = np.quantile(female_scores, female_quantiles)

    plt.scatter(male_sample, female_sample, alpha=0.5)
    plt.plot([0.4, 0.6], [0.4, 0.6], 'r--')  # diagonal line
    plt.xlabel('Male Prediction Scores')
    plt.ylabel('Female Prediction Scores')
    plt.title('Q-Q Plot: Male vs Female Prediction Scores')
    plt.grid(alpha=0.3)

    # Add KS test results to plot
    plt.figtext(0.15, 0.15,
                f"Kolmogorov-Smirnov test:\nStatistic: {ks_stat:.4f}\np-value: {p_value:.4f}",
                bbox=dict(facecolor='white', alpha=0.8),
                fontsize=10)

    plt.savefig(f"{output_dir}/qq_plot.png", dpi=300, bbox_inches='tight')
    plt.close()

    # C. Enhanced histogram/density plot
    plt.figure(figsize=(10, 6), dpi=300)

    # Create histogram/KDE of prediction scores
    sns.histplot(male_df['prediction'], kde=True, color='blue', alpha=0.5, label='Male', bins=20)
    sns.histplot(female_df['prediction'], kde=True, color='red', alpha=0.5, label='Female', bins=20)

    # Add mean lines
    plt.axvline(male_df['prediction'].mean(), color='blue', linestyle='--',
                label=f'Male mean: {male_df["prediction"].mean():.4f}')
    plt.axvline(female_df['prediction'].mean(), color='red', linestyle='--',
                label=f'Female mean: {female_df["prediction"].mean():.4f}')

    # Add KS test results to plot
    plt.text(0.42, plt.gca().get_ylim()[1] * 0.9,
             f"KS test: statistic={ks_stat:.4f}, p-value={p_value:.4f}",
             bbox=dict(facecolor='white', alpha=0.8))

    plt.xlabel('Prediction Score')
    plt.ylabel('Frequency')
    plt.title('Distribution of Prediction Scores by Sex', fontsize=16)
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/prediction_score_distribution.png", dpi=300)
    plt.close()

    print(f"All visualizations saved to {output_dir}")
    return metrics_df


def map_genes_and_download_sequences(gene_list, protein_sequences_dict=None):
    """Map gene symbols to UniProt IDs and download their sequences"""
    if protein_sequences_dict is None:
        protein_sequences_dict = {}

    gene_to_uniprot = {}

    print(f"Mapping and downloading sequences for {len(gene_list)} genes...")
    from mygene import MyGeneInfo
    import requests
    import re
    mg = MyGeneInfo()

    # Process in batches to avoid API limits
    batch_size = 100
    for i in range(0, len(gene_list), batch_size):
        batch = gene_list[i:i+batch_size]
        try:
            # Query for UniProt IDs
            results = mg.querymany(batch, scopes='symbol', fields='uniprot', species='human')

            for result in results:
                if 'uniprot' in result and not 'notfound' in result:
                    gene = result['query']
                    # Extract Swiss-Prot ID if available
                    if isinstance(result['uniprot'], dict) and 'Swiss-Prot' in result['uniprot']:
                        uniprot_id = result['uniprot']['Swiss-Prot']
                        if isinstance(uniprot_id, list) and uniprot_id:
                            uniprot_id = uniprot_id[0]  # Take the first one if multiple
                        gene_to_uniprot[gene] = uniprot_id
                    elif isinstance(result['uniprot'], list) and len(result['uniprot']) > 0:
                        uniprot_id = result['uniprot'][0]  # Default to first
                        gene_to_uniprot[gene] = uniprot_id

                    # Download the sequence if we have a UniProt ID
                    if gene in gene_to_uniprot and gene_to_uniprot[gene] not in protein_sequences_dict:
                        uniprot_id = gene_to_uniprot[gene]
                        try:
                            # Primary source: UniProt REST API
                            url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
                            response = requests.get(url)
                            if response.status_code == 200:
                                # Parse FASTA format
                                fasta_text = response.text
                                lines = fasta_text.strip().split('\n')
                                if len(lines) > 1:  # Valid FASTA has header + sequence
                                    sequence = ''.join(lines[1:])  # Skip header, join sequence lines
                                    protein_sequences_dict[uniprot_id] = sequence
                            else:
                                # Fallback source: NCBI
                                ncbi_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={uniprot_id}&rettype=fasta&retmode=text"
                                ncbi_response = requests.get(ncbi_url)
                                if ncbi_response.status_code == 200:
                                    # Parse FASTA from NCBI
                                    fasta_text = ncbi_response.text
                                    lines = fasta_text.strip().split('\n')
                                    if len(lines) > 1:
                                        sequence = ''.join(lines[1:])
                                        protein_sequences_dict[uniprot_id] = sequence
                        except Exception as e:
                            print(f"Error downloading sequence for {uniprot_id}: {e}")
        except Exception as e:
            print(f"Error in batch query: {e}")

    print(f"Successfully mapped {len(gene_to_uniprot)} genes to UniProt IDs")
    print(f"Downloaded {len(protein_sequences_dict)} protein sequences")

    return gene_to_uniprot, protein_sequences_dict

# Main function to run both sexes in a single pass
def main():
    print("Starting memory-optimized DL-PPI pipeline for both sexes...")
    run_params = {
        "model": "facebook/esm2_t33_650M_UR50D",
        "batch_size_train": 8,
        "batch_size_predict": 4,
        "epochs": 20,
        "patience": 3,
        "lr": 1e-4,
        "weight_decay": 1e-5,
        "embedding_dim": 1280,
        "timestamp": time.strftime("%Y%m%d-%H%M%S")
    }

    # Create output directories
    out_dir = "DL-PPI outputs"
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(f"{out_dir}/data", exist_ok=True)
    os.makedirs(f"{out_dir}/models", exist_ok=True)
    os.makedirs(f"{out_dir}/results", exist_ok=True)

    # Create sex-specific output directories
    male_dir = ensure_path(f"{out_dir}/results/male")
    female_dir = ensure_path(f"{out_dir}/results/female")
    comp_dir = ensure_path(f"{out_dir}/results/comparison")

    # Load inflammatory pathway data
    pathway_file = "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/KEGG outputs/kegg_inflammatory_pathways.csv"
    pathway_df = load_inflammatory_pathways(pathway_file)

    if pathway_df.empty:
        print("WARNING: Pathway data unavailable!")
        return

    # Load or create training data
    if not os.path.exists(f"{out_dir}/data/known_interactions.csv"):
        print("\nCreating training data for inflammatory pathways...")
        training_df, protein_sequences_dict = create_complement_training_data(out_dir)
    else:
        print("\nLoading existing training data...")
        training_df = pd.read_csv(f"{out_dir}/data/known_interactions.csv")

        # Load protein sequences from FASTA
        protein_sequences_dict = {}
        with open(f"{out_dir}/data/protein_sequences.fasta", 'r') as f:
            current_id = None
            current_seq = ""
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id:
                        protein_sequences_dict[current_id] = current_seq
                    current_id = line[1:]
                    current_seq = ""
                else:
                    current_seq += line
            if current_id:
                protein_sequences_dict[current_id] = current_seq

    # Create training pairs and labels
    train_pairs = list(zip(training_df['protein1'], training_df['protein2']))
    train_labels = training_df['interaction_label'].tolist()

    print(
        f"Loaded {len(train_pairs)} training pairs, {sum(train_labels)} positive, {len(train_labels) - sum(train_labels)} negative")

    # Initialize deep learning model - ONCE for both sexes
    embedding_cache_file = f"{out_dir}/data/embedding_cache.pkl"
    dl_model = DeepLearningPPI()

    # Load embedding cache if available
    dl_model.load_embedding_cache(embedding_cache_file)

    # Check if model exists or needs training
    if os.path.exists(f"{out_dir}/models/classifier_head.pt"):
        print("\nFound existing trained model.")
        # Load model for inference
        dl_model.classifier.load_state_dict(torch.load(f"{out_dir}/models/classifier_head.pt"))
    else:
        print("\nTraining new model with memory optimizations...")
        # Train model with smaller batch size for M1 Mac
        dl_model = train_ppi_model(train_pairs, train_labels, protein_sequences_dict,
                                   batch_size=8,  # Smaller batch size for M1 Mac
                                   epochs=20,
                                   patience=3,
                                   output_dir=f"{out_dir}/models")

    # Precompute embeddings for all proteins to avoid recomputing
    precompute_all_embeddings(protein_sequences_dict, dl_model, embedding_cache_file)

    # Extract all genes from inflammatory pathways
    all_pathway_genes = pathway_df.gene.unique().tolist()
    # Filter out entries that aren't gene symbols
    all_pathway_genes = [g for g in all_pathway_genes if isinstance(g, str) and not g.startswith('[')]
    print(f"Loaded {len(all_pathway_genes)} genes from inflammatory pathways")

    # Load sex-specific DEGs
    male_genes = load_sex_specific_genes("M")
    female_genes = load_sex_specific_genes("F")
    print(f"Loaded {len(male_genes)} male-specific and {len(female_genes)} female-specific genes")

    # Map genes to UniProt IDs - do this ONCE for all genes
    print("\nMapping genes to UniProt IDs...")
    from mygene import MyGeneInfo
    mg = MyGeneInfo()

    # Get mappings for all genes at once
    all_genes = set(all_pathway_genes + male_genes + female_genes)
    gene_info = mg.querymany(all_genes, scopes='symbol', fields='uniprot', species='human', returnall=True)

    # Process results
    gene_to_uniprot = {}
    for entry in gene_info['out']:
        if 'uniprot' in entry and 'Swiss-Prot' in entry['uniprot']:
            gene = entry.get('query', '')
            uniprot_id = entry['uniprot']['Swiss-Prot']
            if isinstance(uniprot_id, list):
                uniprot_id = uniprot_id[0]  # Take the first one
            gene_to_uniprot[gene] = uniprot_id

    # Create reverse mapping for visualization
    uniprot_to_gene = {v: k for k, v in gene_to_uniprot.items()}

    # Add gene mapping to model for convenient access
    dl_model.uniprot_to_gene = uniprot_to_gene

    # Process male data
    male_results = None
    print("\n========== PROCESSING MALE DATA ==========")
    if male_genes:
        # Generate pathway-based pairs for males
        male_pairs = generate_pathway_based_pairs(male_genes, pathway_df, protein_sequences_dict, gene_to_uniprot)
        male_test_pairs = [(p[0], p[1]) for p in male_pairs]

        if male_test_pairs:
            male_results = predict_ppi_deep_learning(
                male_test_pairs,
                protein_sequences_dict,
                model=dl_model,
                batch_size=run_params["batch_size_predict"],
                output_dir=male_dir
            )

            # Add metadata to results
            if not male_results.empty:
                # Add gene names
                male_results['gene1'] = male_results['acc1'].map(uniprot_to_gene)
                male_results['gene2'] = male_results['acc2'].map(uniprot_to_gene)

                # Add pathway information
                def get_pathway_info_for_row(row):
                    gene1 = row['gene1'] if pd.notna(row['gene1']) else ""
                    gene2 = row['gene2'] if pd.notna(row['gene2']) else ""
                    return get_pathway_info(gene1, gene2, pathway_df)

                male_results['pathway_info'] = male_results.apply(get_pathway_info_for_row, axis=1)
                male_results.to_csv(f"{male_dir}/dl_ppi_predictions_with_metadata.csv", index=False)

                # Save visualization-friendly format
                male_viz_df = create_viz_data(male_results, uniprot_to_gene, threshold=0.4)
                male_viz_df.to_csv(f"{out_dir}/results/dl_ppi_viz_data_male_0.4_genes.csv", index=False)

    # Force garbage collection to free memory before processing female data
    gc.collect()
    if hasattr(torch.cuda, 'empty_cache'):
        torch.cuda.empty_cache()

    # Process female data
    female_results = None
    print("\n========== PROCESSING FEMALE DATA ==========")
    if female_genes:
        # Generate pathway-based pairs for females
        female_pairs = generate_pathway_based_pairs(female_genes, pathway_df, protein_sequences_dict, gene_to_uniprot)
        female_test_pairs = [(p[0], p[1]) for p in female_pairs]

        if female_test_pairs:
            female_results = predict_ppi_deep_learning(
                female_test_pairs,
                protein_sequences_dict,
                model=dl_model,
                batch_size=run_params["batch_size_predict"],
                output_dir=female_dir
            )

            # Add metadata to results
            if not female_results.empty:
                # Add gene names
                female_results['gene1'] = female_results['acc1'].map(uniprot_to_gene)
                female_results['gene2'] = female_results['acc2'].map(uniprot_to_gene)

                # Add pathway information
                def get_pathway_info_for_row(row):
                    gene1 = row['gene1'] if pd.notna(row['gene1']) else ""
                    gene2 = row['gene2'] if pd.notna(row['gene2']) else ""
                    return get_pathway_info(gene1, gene2, pathway_df)

                female_results['pathway_info'] = female_results.apply(get_pathway_info_for_row, axis=1)
                female_results.to_csv(f"{female_dir}/dl_ppi_predictions_with_metadata.csv", index=False)

                # Save visualization-friendly format
                female_viz_df = create_viz_data(female_results, uniprot_to_gene, threshold=0.4)
                female_viz_df.to_csv(f"{out_dir}/results/dl_ppi_viz_data_female_0.4_genes.csv", index=False)

    # Compare networks if both have results
    if male_results is not None and female_results is not None and not male_results.empty and not female_results.empty:
        print("\n========== COMPARING MALE AND FEMALE NETWORKS ==========")
        # Convert UniProt IDs to gene symbols for visualization
        male_viz_df = pd.DataFrame({
            'protein1': male_results['gene1'],
            'protein2': male_results['gene2'],
            'prediction': male_results['prediction'],
            'method': male_results['method']
        })

        female_viz_df = pd.DataFrame({
            'protein1': female_results['gene1'],
            'protein2': female_results['gene2'],
            'prediction': female_results['prediction'],
            'method': female_results['method']
        })

        # Compare the networks
        network_metrics = compare_networks(male_viz_df, female_viz_df, output_dir=comp_dir)

        # Save metrics
        network_metrics.to_csv(f"{comp_dir}/network_metrics.csv", index=False)

    # Save embedding cache before finishing
    print("\nSaving embedding cache...")
    with open(embedding_cache_file, 'wb') as f:
        pickle.dump(dl_model.embedding_cache, f)

    # Save run parameters for reproducibility
    with open(f"{out_dir}/run_parameters.json", 'w') as f:
        json.dump(run_params, f, indent=2)

    print("\nDL-PPI prediction pipeline completed for both sexes.")
    print(f"Results saved to {out_dir}")
    print(f"Run parameters saved to {out_dir}/run_parameters.json")


if __name__ == "__main__":
    main()
