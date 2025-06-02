# ==============================================================================
# GSE225158 Pseudobulk Data Preparation
# ==============================================================================
# Purpose: Load h5ad single-cell data and create pseudobulk counts for R analysis
# Dataset: GSE225158 - Single-nucleus RNA-seq from human postmortem striatum
# Design:  Paired samples (Caudate + Putamen) from 10 subjects (OUD vs CTL)
# Output:  CSV files compatible with GSE174409 R analysis pipeline
# ==============================================================================

import pandas as pd
import numpy as np
import scanpy as sc
import os
from pathlib import Path

# Configuration
DATA_DIR = "/Users/aumchampaneri/Complement-OUD/Translational_Study/Homo_sapiens/Data/GSE225158"
OUTPUT_DIR = "/Users/aumchampaneri/Complement-OUD/Translational_Study/Homo_sapiens/Results/GSE225158"

def load_h5ad_data():
    """
    Load GSE225158 h5ad file and extract subject IDs from sample identifiers.
    
    Returns:
        adata: AnnData object with added 'True_Subject_ID' column
    """
    print("=== Loading GSE225158 h5ad Data ===")
    
    # Find h5ad file
    h5ad_files = list(Path(DATA_DIR).glob("*.h5ad"))
    if not h5ad_files:
        raise FileNotFoundError(f"No .h5ad files found in {DATA_DIR}")
    
    # Load data
    h5ad_file = h5ad_files[0]
    print(f"Loading: {h5ad_file.name}")
    adata = sc.read_h5ad(h5ad_file)
    print(f"Loaded: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
    
    # Extract true subject IDs from sample IDs (e.g., "C-1366" -> "1366")
    print("Extracting subject IDs from sample identifiers...")
    adata.obs['True_Subject_ID'] = adata.obs['ID'].str.extract(r'[CP]-(\d+)')[0]
    
    # Verify extraction worked
    id_sample = adata.obs[['ID', 'True_Subject_ID', 'Region']].drop_duplicates().head(5)
    print("Sample ID mapping:")
    print(id_sample.to_string(index=False))
    
    return adata

def verify_pairing_structure(adata):
    """
    Verify the paired design structure of the dataset.
    
    Args:
        adata: AnnData object with subject and region information
        
    Returns:
        tuple: (paired_subjects_list, pairing_summary_dict)
    """
    print("\n=== Verifying Paired Design ===")
    
    # Get unique subject-region combinations
    subject_region_df = adata.obs[['True_Subject_ID', 'Region', 'Sex', 'level1']].drop_duplicates()
    
    # Count subjects per region
    subjects_per_region = subject_region_df.groupby('Region')['True_Subject_ID'].nunique()
    print(f"Subjects per region: {subjects_per_region.to_dict()}")
    
    # Find paired subjects (appearing in both regions)
    regions_per_subject = subject_region_df.groupby('True_Subject_ID')['Region'].nunique()
    paired_subjects = regions_per_subject[regions_per_subject == 2].index.tolist()
    
    print(f"Paired subjects found: {len(paired_subjects)}")
    
    # Show pairing structure
    if len(paired_subjects) >= 5:
        print("Sample pairings:")
        for subject in paired_subjects[:5]:
            subject_data = subject_region_df[subject_region_df['True_Subject_ID'] == subject]
            regions = subject_data['Region'].tolist()
            oud_status = subject_data['level1'].iloc[0]
            sex = subject_data['Sex'].iloc[0]
            print(f"  Subject {subject}: {regions} | {oud_status} | {sex}")
    
    # Summary statistics
    summary = {
        'total_subjects': len(paired_subjects),
        'total_samples': len(paired_subjects) * 2,
        'oud_subjects': sum(subject_region_df[subject_region_df['True_Subject_ID'].isin(paired_subjects)]['level1'] == 'OUD') // 2,
        'ctl_subjects': sum(subject_region_df[subject_region_df['True_Subject_ID'].isin(paired_subjects)]['level1'] == 'CTL') // 2
    }
    
    return paired_subjects, summary

def create_pseudobulk_data(adata, paired_subjects):
    """
    Aggregate single-cell counts to pseudobulk by subject and region.
    
    Args:
        adata: AnnData object
        paired_subjects: List of subject IDs to include
        
    Returns:
        tuple: (pseudobulk_counts_df, metadata_df)
    """
    print("\n=== Creating Pseudobulk Data ===")
    
    # Filter to paired subjects only
    paired_mask = adata.obs['True_Subject_ID'].isin(paired_subjects)
    adata_paired = adata[paired_mask, :].copy()
    print(f"Using {adata_paired.n_obs:,} cells from {len(paired_subjects)} paired subjects")
    
    # Create grouping variable (Subject_Region)
    adata_paired.obs['group_id'] = (
        adata_paired.obs['True_Subject_ID'].astype(str) + "_" + 
        adata_paired.obs['Region'].astype(str)
    )
    
    # Get count matrix (use raw counts if available)
    if hasattr(adata_paired, 'raw') and adata_paired.raw is not None:
        counts_matrix = adata_paired.raw.X
        gene_names = adata_paired.raw.var_names
        print("Using raw counts from adata.raw.X")
    else:
        counts_matrix = adata_paired.X
        gene_names = adata_paired.var_names
        print("Using counts from adata.X")
    
    # Convert sparse to dense if needed
    if hasattr(counts_matrix, 'toarray'):
        counts_matrix = counts_matrix.toarray()
    
    # Aggregate counts by group
    print("Aggregating counts by subject-region groups...")
    pseudobulk_counts = []
    sample_metadata = []
    
    for group_id in adata_paired.obs['group_id'].unique():
        # Get cells in this group
        group_mask = adata_paired.obs['group_id'] == group_id
        group_cells = counts_matrix[group_mask, :]
        
        # Sum counts across cells
        pseudobulk_sample = group_cells.sum(axis=0)
        pseudobulk_counts.append(pseudobulk_sample)
        
        # Extract metadata (from first cell in group)
        group_meta = adata_paired.obs[group_mask].iloc[0]
        sample_metadata.append({
            'Sample_ID': group_id,
            'Subject_ID': group_meta['True_Subject_ID'],
            'Region': group_meta['Region'],
            'Sex': group_meta['Sex'],
            'OUD_Status': group_meta['level1'],
            'n_cells': group_mask.sum()
        })
    
    # Create dataframes
    pseudobulk_matrix = np.array(pseudobulk_counts).T  # Genes x Samples
    pseudobulk_df = pd.DataFrame(
        pseudobulk_matrix,
        index=gene_names,
        columns=[meta['Sample_ID'] for meta in sample_metadata]
    )
    
    metadata_df = pd.DataFrame(sample_metadata)
    metadata_df.set_index('Sample_ID', inplace=True)
    
    print(f"Pseudobulk matrix: {pseudobulk_df.shape[0]:,} genes × {pseudobulk_df.shape[1]} samples")
    
    return pseudobulk_df, metadata_df

def save_data_for_r(pseudobulk_df, metadata_df, summary):
    """
    Save pseudobulk data in R-compatible format (matching GSE174409 structure).
    
    Args:
        pseudobulk_df: Pseudobulk count matrix
        metadata_df: Sample metadata
        summary: Dataset summary statistics
    """
    print("\n=== Saving Data for R Analysis ===")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Save files
    counts_file = os.path.join(OUTPUT_DIR, "GSE225158_pseudobulk_counts.csv")
    metadata_file = os.path.join(OUTPUT_DIR, "GSE225158_pseudobulk_metadata.csv")
    
    pseudobulk_df.to_csv(counts_file)
    metadata_df.to_csv(metadata_file)
    
    print(f"✓ Counts saved: {counts_file}")
    print(f"✓ Metadata saved: {metadata_file}")
    
    # Print final summary
    print(f"\n=== Dataset Summary ===")
    print(f"Paired subjects: {summary['total_subjects']}")
    print(f"Total samples: {summary['total_samples']}")
    print(f"OUD subjects: {summary['oud_subjects']}")
    print(f"CTL subjects: {summary['ctl_subjects']}")
    print(f"Genes: {pseudobulk_df.shape[0]:,}")
    print(f"Regions: {metadata_df['Region'].value_counts().to_dict()}")
    print(f"Sex distribution: {metadata_df['Sex'].value_counts().to_dict()}")
    
    print(f"\n✓ Ready for R analysis using GSE174409-compatible workflow!")

def main():
    """
    Main pipeline: Load h5ad data, create pseudobulk aggregation, save for R.
    """
    print("="*70)
    print("GSE225158 PSEUDOBULK DATA PREPARATION")
    print("="*70)
    
    try:
        # Load data and extract subject IDs
        adata = load_h5ad_data()
        
        # Verify paired design structure
        paired_subjects, summary = verify_pairing_structure(adata)
        
        if len(paired_subjects) < 3:
            raise ValueError(f"Insufficient paired subjects found: {len(paired_subjects)}")
        
        # Create pseudobulk data
        pseudobulk_df, metadata_df = create_pseudobulk_data(adata, paired_subjects)
        
        # Save for R analysis
        save_data_for_r(pseudobulk_df, metadata_df, summary)
        
        print("\n" + "="*70)
        print("SUCCESS: Pseudobulk data preparation complete!")
        print("="*70)
        
        return adata, pseudobulk_df, metadata_df
        
    except Exception as e:
        print(f"\nERROR: {e}")
        raise

if __name__ == "__main__":
    results = main()
