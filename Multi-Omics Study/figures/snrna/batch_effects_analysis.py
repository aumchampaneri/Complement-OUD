#!/usr/bin/env python3
"""
QC Plots for Nuclei Analysis

Creates:
1. Box plots showing nuclei counts by brain area and treatment condition, with gender shapes
2. UMAP colored by number of genes per cell

Dataset: N = 22 caudate and putamen samples from M = 12 individuals
"""

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from scipy import stats
import os
import warnings
warnings.filterwarnings('ignore')

# Custom color palette
COLORS = {
    'black': '#000000',
    'light_blue': '#B6DBFF',
    'medium_blue': '#7BB0DF',
    'dark_blue': '#1964B0',
    'teal': '#00C992',
    'dark_teal': '#008A69',
    'dark_green': '#386350',
    'yellow': '#E9DC6D',
    'orange': '#F4A637',
    'red_orange': '#DB5829',
    'brown': '#894B45',
    'light_purple': '#D2BBD7',
    'medium_purple': '#AE75A2',
    'dark_purple': '#882D71',
    'light_gray': '#DEDEDE'
}

def setup_publication_style():
    """Configure matplotlib for publication-ready figures"""
    plt.rcParams.update({
        'font.family': 'Helvetica',
        'font.size': 12,
        'axes.linewidth': 1.5,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.edgecolor': '#000000',
        'axes.labelcolor': '#000000',
        'text.color': '#000000',
        'xtick.color': '#000000',
        'ytick.color': '#000000',
        'xtick.major.width': 1.5,
        'ytick.major.width': 1.5,
        'legend.frameon': True,
        'legend.facecolor': 'white',
        'legend.edgecolor': '#000000',
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'pdf.fonttype': 42,
        'ps.fonttype': 42
    })

def load_data(data_path):
    """Load single-cell data from h5ad file"""
    print(f"📊 Loading data from: {data_path}")
    adata = sc.read_h5ad(data_path)
    print(f"✓ Data loaded: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")
    return adata

def create_sample_qc_dataframe(adata):
    """Create per-sample QC metrics dataframe"""
    print("📊 Creating per-sample QC metrics...")
    
    # Define the columns we know exist based on your metadata
    sample_col = 'ID'
    region_col = 'Region'  # Brain area (caudate/putamen)
    diagnosis_col = 'Dx_OUD'  # Treatment condition
    sex_col = 'Sex'  # Gender
    
    print(f"Using columns: ID={sample_col}, Region={region_col}, Diagnosis={diagnosis_col}, Sex={sex_col}")
    
    # Calculate per-sample statistics
    sample_stats = []
    
    for sample_id in adata.obs[sample_col].unique():
        sample_mask = adata.obs[sample_col] == sample_id
        sample_data = adata.obs[sample_mask]
        
        if len(sample_data) > 0:
            stats_row = {
                'sample_id': sample_id,
                'brain_area': sample_data[region_col].iloc[0],
                'treatment': sample_data[diagnosis_col].iloc[0],
                'gender': sample_data[sex_col].iloc[0],
                'n_nuclei': len(sample_data),  # Number of nuclei per sample
                'mean_genes': sample_data['n_genes_by_counts'].mean(),
                'median_genes': sample_data['n_genes_by_counts'].median(),
                'mean_total_counts': sample_data['total_counts'].mean(),
                'mean_mt_pct': sample_data['pct_counts_mt'].mean()
            }
            sample_stats.append(stats_row)
    
    sample_df = pd.DataFrame(sample_stats)
    print(f"✓ Created QC metrics for {len(sample_df)} samples")
    
    return sample_df

def create_nuclei_count_boxplots(sample_df, output_dir):
    """Create box plots showing nuclei counts by brain area and treatment, with gender shapes"""
    print("🎨 Creating nuclei count box plots...")
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    fig.suptitle('Nuclei Count Analysis\nN = 22 samples from M = 12 individuals', 
                 fontsize=16, fontweight='bold', y=0.95)
    
    # Define gender markers
    gender_markers = {'Male': 'o', 'Female': '^', 'M': 'o', 'F': '^'}
    gender_colors = {'Male': COLORS['dark_blue'], 'Female': COLORS['red_orange'], 
                    'M': COLORS['dark_blue'], 'F': COLORS['red_orange']}
    
    # Plot 1: Nuclei count by brain area
    ax1.set_title('Nuclei Count by Brain Area', fontweight='bold', fontsize=14)
    
    # Create boxplot
    brain_areas = sample_df['brain_area'].unique()
    area_data = [sample_df[sample_df['brain_area'] == area]['n_nuclei'].values for area in brain_areas]
    
    bp1 = ax1.boxplot(area_data, labels=brain_areas, patch_artist=True, 
                     notch=True, showfliers=False)
    
    # Color boxes
    colors = [COLORS['dark_blue'], COLORS['teal']]
    for patch, color in zip(bp1['boxes'], colors[:len(brain_areas)]):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
        patch.set_edgecolor('black')
    
    # Add individual points with gender shapes
    for i, area in enumerate(brain_areas):
        area_data = sample_df[sample_df['brain_area'] == area]
        
        for _, row in area_data.iterrows():
            gender = row['gender']
            marker = gender_markers.get(gender, 'o')
            color = gender_colors.get(gender, COLORS['black'])
            
            x = i + 1 + np.random.normal(0, 0.04)
            y = row['n_nuclei']
            ax1.scatter(x, y, marker=marker, s=60, color=color, alpha=0.8, 
                       edgecolors='black', linewidth=1, zorder=3)
    
    # Statistical test
    if len(area_data) == 2:
        stat, p_val = stats.mannwhitneyu(area_data[0], area_data[1])
        test_name = "Mann-Whitney U"
    else:
        stat, p_val = stats.kruskal(*area_data)
        test_name = "Kruskal-Wallis"
    
    ax1.text(0.02, 0.98, f'{test_name} p={p_val:.3e}', 
            transform=ax1.transAxes, va='top', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax1.set_ylabel('Number of Nuclei per Sample', fontweight='bold')
    ax1.set_xlabel('Brain Area', fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Nuclei count by treatment condition
    ax2.set_title('Nuclei Count by Treatment Condition', fontweight='bold', fontsize=14)
    
    # Create boxplot
    treatments = sample_df['treatment'].unique()
    treatment_data = [sample_df[sample_df['treatment'] == treat]['n_nuclei'].values for treat in treatments]
    
    bp2 = ax2.boxplot(treatment_data, labels=treatments, patch_artist=True, 
                     notch=True, showfliers=False)
    
    # Color boxes
    treat_colors = [COLORS['orange'], COLORS['medium_blue']]
    for patch, color in zip(bp2['boxes'], treat_colors[:len(treatments)]):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
        patch.set_edgecolor('black')
    
    # Add individual points with gender shapes
    for i, treat in enumerate(treatments):
        treat_data = sample_df[sample_df['treatment'] == treat]
        
        for _, row in treat_data.iterrows():
            gender = row['gender']
            marker = gender_markers.get(gender, 'o')
            color = gender_colors.get(gender, COLORS['black'])
            
            x = i + 1 + np.random.normal(0, 0.04)
            y = row['n_nuclei']
            ax2.scatter(x, y, marker=marker, s=60, color=color, alpha=0.8, 
                       edgecolors='black', linewidth=1, zorder=3)
    
    # Statistical test
    if len(treatment_data) == 2:
        stat, p_val = stats.mannwhitneyu(treatment_data[0], treatment_data[1])
        test_name = "Mann-Whitney U"
    else:
        stat, p_val = stats.kruskal(*treatment_data)
        test_name = "Kruskal-Wallis"
    
    ax2.text(0.02, 0.98, f'{test_name} p={p_val:.3e}', 
            transform=ax2.transAxes, va='top', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax2.set_ylabel('Number of Nuclei per Sample', fontweight='bold')
    ax2.set_xlabel('Treatment Condition', fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    # Add legend for gender markers
    legend_elements = []
    for gender in sample_df['gender'].unique():
        marker = gender_markers.get(gender, 'o')
        color = gender_colors.get(gender, COLORS['black'])
        legend_elements.append(plt.Line2D([0], [0], marker=marker, color='w', 
                                        markerfacecolor=color, markersize=8, 
                                        markeredgecolor='black', label=gender))
    
    ax2.legend(handles=legend_elements, title='Gender', loc='upper right', 
              title_fontsize=12, fontsize=10)
    
    # Style axes
    for ax in [ax1, ax2]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color(COLORS['black'])
        ax.spines['bottom'].set_color(COLORS['black'])
        ax.tick_params(colors=COLORS['black'], labelsize=11)
    
    plt.tight_layout()
    
    # Save figure
    nuclei_fig_path = f"{output_dir}/Figure_Nuclei_Count_Analysis"
    plt.savefig(f"{nuclei_fig_path}.png", dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(f"{nuclei_fig_path}.pdf", bbox_inches='tight', facecolor='white')
    
    print(f"💾 Nuclei count plots saved to: {nuclei_fig_path}.png/pdf")
    plt.close()
    
    return nuclei_fig_path

def create_umap_genes_plot(adata, output_dir):
    """Create UMAP plot colored by number of genes per cell"""
    print("🎨 Creating UMAP plot colored by genes per cell...")
    
    # Check if UMAP coordinates exist
    if 'X_umap' not in adata.obsm:
        print("⚠ No UMAP coordinates found. Computing UMAP...")
        sc.pp.neighbors(adata, use_rep='X_scVI')
        sc.tl.umap(adata)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Get UMAP coordinates and gene counts
    umap_coords = adata.obsm['X_umap']
    gene_counts = adata.obs['n_genes_by_counts']
    
    # Create custom colormap (yellow to red gradient)
    colors_list = [COLORS['yellow'], COLORS['orange'], COLORS['red_orange'], COLORS['brown']]
    cmap = mcolors.LinearSegmentedColormap.from_list('genes', colors_list, N=100)
    
    # Create scatter plot
    scatter = ax.scatter(umap_coords[:, 0], umap_coords[:, 1], 
                        c=gene_counts, cmap=cmap, s=1.0, alpha=0.7, 
                        rasterized=True, edgecolors='none')
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.8, aspect=20)
    cbar.set_label('Number of Genes per Cell', fontsize=12, fontweight='bold')
    cbar.ax.tick_params(labelsize=10)
    
    # Style plot
    ax.set_title('UMAP - Genes per Cell\nN = 75,942 nuclei', 
                fontsize=16, fontweight='bold', pad=20)
    ax.set_xlabel('UMAP 1', fontsize=12, fontweight='bold')
    ax.set_ylabel('UMAP 2', fontsize=12, fontweight='bold')
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color(COLORS['black'])
    ax.spines['bottom'].set_color(COLORS['black'])
    ax.tick_params(colors=COLORS['black'], labelsize=10)
    
    # Add statistics text
    mean_genes = gene_counts.mean()
    median_genes = gene_counts.median()
    std_genes = gene_counts.std()
    
    stats_text = f"Mean: {mean_genes:.0f}\nMedian: {median_genes:.0f}\nStd: {std_genes:.0f}"
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, va='top', fontsize=10,
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    
    # Save figure
    umap_fig_path = f"{output_dir}/Figure_UMAP_Genes_per_Cell"
    plt.savefig(f"{umap_fig_path}.png", dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(f"{umap_fig_path}.pdf", bbox_inches='tight', facecolor='white')
    
    print(f"💾 UMAP genes plot saved to: {umap_fig_path}.png/pdf")
    plt.close()
    
    return umap_fig_path

def print_sample_summary(sample_df):
    """Print summary of sample composition"""
    print("\n📋 SAMPLE COMPOSITION SUMMARY")
    print("="*50)
    print(f"Total samples analyzed: {len(sample_df)}")
    
    print(f"\nBrain area distribution:")
    area_counts = sample_df['brain_area'].value_counts()
    for area, count in area_counts.items():
        print(f"  • {area}: {count} samples")
    
    print(f"\nTreatment condition distribution:")
    treatment_counts = sample_df['treatment'].value_counts()
    for treatment, count in treatment_counts.items():
        print(f"  • {treatment}: {count} samples")
    
    print(f"\nGender distribution:")
    gender_counts = sample_df['gender'].value_counts()
    for gender, count in gender_counts.items():
        print(f"  • {gender}: {count} samples")
    
    print(f"\nNuclei per sample:")
    print(f"  • Mean: {sample_df['n_nuclei'].mean():.0f}")
    print(f"  • Median: {sample_df['n_nuclei'].median():.0f}")
    print(f"  • Range: {sample_df['n_nuclei'].min():.0f} - {sample_df['n_nuclei'].max():.0f}")

def main():
    """Main function to run QC analysis"""
    
    print("🎯 NUCLEI QC ANALYSIS")
    print("="*60)
    print("Creating QC plots for nuclei count analysis and UMAP visualization")
    print("="*60)
    
    # Setup
    setup_publication_style()
    print("✓ Publication-ready matplotlib settings configured")
    
    # Data paths
    scvi_path = '/Users/aumchampaneri/Complement-OUD/Multi-Omics Study/data/processed/snrna_scvi/GSE225158_annotated_scvi.h5ad'
    output_dir = '/Users/aumchampaneri/Complement-OUD/Multi-Omics Study/figures/snrna/outputs'
    os.makedirs(output_dir, exist_ok=True)
    
    # Load data
    adata = load_data(scvi_path)
    
    # Create per-sample QC dataframe
    sample_df = create_sample_qc_dataframe(adata)
    
    # Print sample summary
    print_sample_summary(sample_df)
    
    # Create nuclei count box plots
    print("\n" + "="*60)
    print("📊 CREATING NUCLEI COUNT BOX PLOTS")
    print("="*60)
    
    try:
        nuclei_fig_path = create_nuclei_count_boxplots(sample_df, output_dir)
        print(f"✅ Nuclei count analysis completed!")
        
    except Exception as e:
        print(f"❌ Error in nuclei count analysis: {e}")
        import traceback
        traceback.print_exc()
    
    # Create UMAP genes plot
    print("\n" + "="*60)
    print("📊 CREATING UMAP GENES PLOT")
    print("="*60)
    
    try:
        umap_fig_path = create_umap_genes_plot(adata, output_dir)
        print(f"✅ UMAP genes analysis completed!")
        
    except Exception as e:
        print(f"❌ Error in UMAP genes analysis: {e}")
        import traceback
        traceback.print_exc()
    
    print(f"\n✅ All analyses completed successfully!")
    print(f"📁 OUTPUT FILES:")
    print(f"  • {nuclei_fig_path}.png/pdf")
    print(f"  • {umap_fig_path}.png/pdf")

if __name__ == "__main__":
    main()