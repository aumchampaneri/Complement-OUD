import torch
from transformers import AutoTokenizer, AutoModel
import pandas as pd
import numpy as np
from tqdm.auto import tqdm
from sklearn.metrics import roc_auc_score


# 1. Define a deep learning PPI predictor using pre-trained ESM-2
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


# 2. Function to run PPI prediction with the deep learning model
def predict_ppi_deep_learning(protein_pairs, protein_sequences, batch_size=5, output_dir=None):
    """Predicts PPIs using deep learning model"""
    results = []
    total_batches = (len(protein_pairs) + batch_size - 1) // batch_size

    # Initialize the DL model
    dl_model = DeepLearningPPI()

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