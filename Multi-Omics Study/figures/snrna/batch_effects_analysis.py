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

# Nature journal publication-ready color palette
COLORS = {
    'black': '#000000',
    'dark_gray': '#2C2C2C',
    'medium_gray': '#666666',
    'light_gray': '#CCCCCC',
    'very_light_gray': '#F5F5F5',

    # Primary colors - colorblind friendly
    'blue': '#1F77B4',        # Primary blue
    'orange': '#FF7F0E',      # Primary orange
    'green': '#2CA02C',       # Primary green
    'red': '#D62728',         # Primary red
    'purple': '#9467BD',      # Primary purple

    # Lighter variants for boxes
    'light_blue': '#AEC7E8',
    'light_orange': '#FFBB78',
    'light_green': '#98DF8A',
    'light_red': '#FF9896',
    'light_purple': '#C5B0D5',

    # Darker variants for contrast
    'dark_blue': '#0B5394',
    'dark_orange': '#CC5500',
    'dark_green': '#1B7017',
    'dark_red': '#AA1F20',
    'dark_purple': '#6B4C93',
}

def setup_publication_style():
    """Configure matplotlib for Nature journal publication standards"""
    plt.rcParams.update({
        # Font settings - Nature prefers Arial/Helvetica
        'font.family': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 8,           # Nature standard font size
        'font.weight': 'normal',

        # Figure settings
        'figure.dpi': 300,
        'figure.facecolor': 'white',
        'figure.edgecolor': 'white',
        'savefig.dpi': 300,
        'savefig.facecolor': 'white',
        'savefig.edgecolor': 'white',
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.1,

        # Axes settings
        'axes.linewidth': 0.8,    # Thinner lines for cleaner look
        'axes.edgecolor': 'black',
        'axes.facecolor': 'white',
        'axes.labelcolor': 'black',
        'axes.labelsize': 8,
        'axes.titlesize': 9,
        'axes.titleweight': 'bold',
        'axes.titlepad': 10,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.spines.left': True,
        'axes.spines.bottom': True,
        'axes.grid': False,
        'axes.axisbelow': True,

        # Tick settings
        'xtick.major.size': 3,
        'xtick.minor.size': 2,
        'xtick.major.width': 0.8,
        'xtick.minor.width': 0.6,
        'xtick.direction': 'out',
        'xtick.color': 'black',
        'xtick.labelsize': 7,
        'ytick.major.size': 3,
        'ytick.minor.size': 2,
        'ytick.major.width': 0.8,
        'ytick.minor.width': 0.6,
        'ytick.direction': 'out',
        'ytick.color': 'black',
        'ytick.labelsize': 7,

        # Legend settings
        'legend.frameon': True,
        'legend.facecolor': 'white',
        'legend.edgecolor': 'black',
        'legend.fontsize': 7,
        'legend.title_fontsize': 8,
        'legend.borderpad': 0.4,
        'legend.columnspacing': 1.0,
        'legend.handlelength': 1.5,
        'legend.handletextpad': 0.5,

        # Text settings
        'text.color': 'black',

        # PDF/vector output settings
        'pdf.fonttype': 42,       # Embed fonts as vector text
        'ps.fonttype': 42,
        'svg.fonttype': 'none',
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
    """Create Nature-style grouped box plot showing nuclei counts by brain area and treatment, with gender shapes"""
    print("🎨 Creating publication-ready grouped nuclei count box plot...")

    # Create figure with proper Nature dimensions for single grouped plot
    fig, ax = plt.subplots(1, 1, figsize=(5.0, 4.0))  # Nature single column width, slightly taller

    # Define gender markers and colors (colorblind-friendly)
    gender_markers = {'Male': 'o', 'Female': 's', 'M': 'o', 'F': 's'}  # circle and square
    gender_colors = {'Male': COLORS['blue'], 'Female': COLORS['red'],
                    'M': COLORS['blue'], 'F': COLORS['red']}

    # Create grouped box plot using seaborn for better control
    # Prepare data for grouped plotting
    brain_areas = sorted(sample_df['brain_area'].unique())
    treatments = sorted(sample_df['treatment'].unique())

    # Set up positions for grouped boxes
    positions = []
    labels = []
    group_data = []
    box_colors = []

    # Define colors for each treatment within regions
    treat_colors = {'None': COLORS['light_blue'], 'OUD': COLORS['light_orange']}
    treat_edge_colors = {'None': COLORS['dark_blue'], 'OUD': COLORS['dark_orange']}

    pos = 1
    for i, area in enumerate(brain_areas):
        for j, treat in enumerate(treatments):
            subset = sample_df[(sample_df['brain_area'] == area) & (sample_df['treatment'] == treat)]
            if len(subset) > 0:
                positions.append(pos)
                labels.append(f"{area}\n{treat}")
                group_data.append(subset['n_nuclei'].values)
                box_colors.append(treat_colors[treat])
                pos += 1
        pos += 0.5  # Add space between brain regions

    # Create grouped boxplot
    bp = ax.boxplot(group_data, positions=positions, patch_artist=True,
                   notch=False, showfliers=False, widths=0.4,
                   boxprops=dict(linewidth=0.8),
                   medianprops=dict(linewidth=1.2, color='black'),
                   whiskerprops=dict(linewidth=0.8),
                   capprops=dict(linewidth=0.8))

    # Color boxes
    for patch, color in zip(bp['boxes'], box_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.8)
        patch.set_edgecolor(COLORS['dark_gray'])

    # Add individual points with gender shapes and treatment-specific jitter
    for i, (area, treat) in enumerate([(a, t) for a in brain_areas for t in treatments]):
        subset = sample_df[(sample_df['brain_area'] == area) & (sample_df['treatment'] == treat)]
        if len(subset) > 0:
            pos_idx = [idx for idx, label in enumerate(labels) if f"{area}\n{treat}" == label][0]

            for _, row in subset.iterrows():
                gender = row['gender']
                marker = gender_markers.get(gender, 'o')
                color = gender_colors.get(gender, COLORS['dark_gray'])

                # Add slight horizontal jitter
                x = positions[pos_idx] + np.random.normal(0, 0.02)
                y = row['n_nuclei']
                ax.scatter(x, y, marker=marker, s=25, color=color, alpha=0.9,
                          edgecolors='white', linewidth=0.8, zorder=5)

    # Set labels and styling
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel('Nuclei per sample', fontweight='bold', fontsize=8)
    ax.set_title('Nuclei count by brain region and treatment', fontweight='bold', fontsize=9, pad=10)

    # Create treatment legend (boxes)
    treat_legend_elements = []
    for treat in treatments:
        treat_legend_elements.append(plt.Rectangle((0, 0), 1, 1,
                                                 facecolor=treat_colors[treat],
                                                 edgecolor=COLORS['dark_gray'],
                                                 linewidth=0.8, alpha=0.8,
                                                 label=treat))

    # Create gender legend (markers)
    gender_legend_elements = []
    unique_genders = sorted(sample_df['gender'].unique())
    for gender in unique_genders:
        marker = gender_markers.get(gender, 'o')
        color = gender_colors.get(gender, COLORS['dark_gray'])
        gender_legend_elements.append(plt.Line2D([0], [0], marker=marker, color='w',
                                               markerfacecolor=color, markersize=6,
                                               markeredgecolor='white', markeredgewidth=0.8,
                                               label=gender))

    # Position legends
    legend1 = ax.legend(handles=treat_legend_elements, title='Treatment',
                       loc='upper left', frameon=True, fancybox=False, shadow=False)
    ax.add_artist(legend1)
    ax.legend(handles=gender_legend_elements, title='Sex', loc='upper right',
             frameon=True, fancybox=False, shadow=False)

    # Statistical testing - perform 2-way ANOVA-style tests
    from scipy.stats import f_oneway

    # Test main effect of brain region
    caudate_data = sample_df[sample_df['brain_area'] == 'Caudate']['n_nuclei']
    putamen_data = sample_df[sample_df['brain_area'] == 'Putamen']['n_nuclei']
    region_stat, region_p = stats.mannwhitneyu(caudate_data, putamen_data, alternative='two-sided')

    # Test main effect of treatment
    none_data = sample_df[sample_df['treatment'] == 'None']['n_nuclei']
    oud_data = sample_df[sample_df['treatment'] == 'OUD']['n_nuclei']
    treat_stat, treat_p = stats.mannwhitneyu(none_data, oud_data, alternative='two-sided')

    # Add statistical annotations
    stats_text = []
    if region_p < 0.001:
        stats_text.append("Region: P < 0.001")
    else:
        stats_text.append(f"Region: P = {region_p:.3f}")

    if treat_p < 0.001:
        stats_text.append("Treatment: P < 0.001")
    else:
        stats_text.append(f"Treatment: P = {treat_p:.3f}")

    ax.text(0.02, 0.98, '\n'.join(stats_text), transform=ax.transAxes,
           va='top', ha='left', fontsize=7,
           bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                    edgecolor=COLORS['light_gray'], linewidth=0.5, alpha=0.95))

    # Style axes according to Nature standards
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color(COLORS['black'])
    ax.spines['bottom'].set_color(COLORS['black'])
    ax.spines['left'].set_linewidth(0.8)
    ax.spines['bottom'].set_linewidth(0.8)

    # Style ticks
    ax.tick_params(axis='both', which='major', labelsize=7,
                  colors=COLORS['black'], width=0.8, length=3)
    ax.tick_params(axis='both', which='minor', width=0.6, length=2)

    # Add subtle grid
    ax.grid(True, alpha=0.2, linewidth=0.5, color=COLORS['light_gray'])
    ax.set_axisbelow(True)

    # Adjust layout
    plt.tight_layout(pad=0.5)

    # Save figure with Nature specifications
    nuclei_fig_path = f"{output_dir}/Figure_Nuclei_Count_Combined_Analysis"
    plt.savefig(f"{nuclei_fig_path}.png", dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none', pad_inches=0.1)
    plt.savefig(f"{nuclei_fig_path}.pdf", bbox_inches='tight',
                facecolor='white', edgecolor='none', pad_inches=0.1)
    plt.savefig(f"{nuclei_fig_path}.eps", bbox_inches='tight',
                facecolor='white', edgecolor='none', pad_inches=0.1)  # EPS for vector graphics

    print(f"💾 Publication-ready combined nuclei count plot saved to: {nuclei_fig_path}.[png/pdf/eps]")
    plt.close()

    return nuclei_fig_path

def create_umap_genes_plot(adata, output_dir):
    """Create Nature-style UMAP plot colored by number of genes per cell"""
    print("🎨 Creating publication-ready UMAP plot colored by genes per cell...")

    # Check if UMAP coordinates exist
    if 'X_umap' not in adata.obsm:
        print("⚠️ No UMAP coordinates found. Computing UMAP...")
        sc.pp.neighbors(adata, use_rep='X_scVI')
        sc.tl.umap(adata)

    # Create figure with Nature single-column width
    fig, ax = plt.subplots(figsize=(3.5, 3.5))  # Square format, Nature single column

    # Get UMAP coordinates and gene counts
    umap_coords = adata.obsm['X_umap']
    gene_counts = adata.obs['n_genes_by_counts']

    # Create perceptually uniform colormap (viridis-like but Nature-friendly)
    colors_list = ['#440154', '#31688e', '#35b779', '#fde725']  # Viridis-inspired
    cmap = mcolors.LinearSegmentedColormap.from_list('genes_nature', colors_list, N=256)

    # Create scatter plot with appropriate point size for publication
    scatter = ax.scatter(umap_coords[:, 0], umap_coords[:, 1],
                        c=gene_counts, cmap=cmap, s=0.8, alpha=0.8,
                        rasterized=True, edgecolors='none')

    # Add colorbar with proper styling
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.8, aspect=15, pad=0.02)
    cbar.set_label('Genes per cell', fontsize=8, fontweight='bold', labelpad=8)
    cbar.ax.tick_params(labelsize=7, width=0.8, length=3)

    # Style colorbar
    cbar.ax.yaxis.set_ticks_position('right')
    cbar.outline.set_linewidth(0.8)
    cbar.outline.set_edgecolor('black')

    # Style plot according to Nature standards
    ax.set_title('Gene expression diversity', fontsize=9, fontweight='bold', pad=8)
    ax.set_xlabel('UMAP 1', fontsize=8, fontweight='bold')
    ax.set_ylabel('UMAP 2', fontsize=8, fontweight='bold')

    # Remove spines and ticks for cleaner look (common for UMAP plots)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])

    # Add scale bar or sample size annotation
    n_cells = adata.shape[0]
    ax.text(0.02, 0.98, f'n = {n_cells:,} cells',
           transform=ax.transAxes, va='top', ha='left', fontsize=7,
           bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                    edgecolor=COLORS['light_gray'], linewidth=0.5, alpha=0.9))

    # Add gene count statistics in a clean format
    mean_genes = gene_counts.mean()
    median_genes = gene_counts.median()

    stats_text = f'Mean: {mean_genes:.0f}\nMedian: {median_genes:.0f}'
    ax.text(0.98, 0.02, stats_text, transform=ax.transAxes,
           va='bottom', ha='right', fontsize=7,
           bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                    edgecolor=COLORS['light_gray'], linewidth=0.5, alpha=0.9))

    # Ensure equal aspect ratio for UMAP
    ax.set_aspect('equal', adjustable='box')

    plt.tight_layout(pad=0.5)

    # Save figure with multiple formats for Nature submission
    umap_fig_path = f"{output_dir}/Figure_UMAP_Genes_per_Cell"
    plt.savefig(f"{umap_fig_path}.png", dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none', pad_inches=0.05)
    plt.savefig(f"{umap_fig_path}.pdf", bbox_inches='tight',
                facecolor='white', edgecolor='none', pad_inches=0.05)
    plt.savefig(f"{umap_fig_path}.eps", bbox_inches='tight',
                facecolor='white', edgecolor='none', pad_inches=0.05)

    print(f"💾 Publication-ready UMAP plot saved to: {umap_fig_path}.[png/pdf/eps]")
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
    print(f"  • Total nuclei: {sample_df['n_nuclei'].sum():,}")

def main():
    """Main function to run publication-ready QC analysis"""

    print("🎯 PUBLICATION-READY NUCLEI QC ANALYSIS")
    print("="*60)
    print("Creating Nature journal-style QC plots")
    print("• Combined box plot: nuclei count by region and treatment with statistical testing")
    print("• UMAP: gene expression diversity visualization")
    print("="*60)

    # Setup
    setup_publication_style()
    print("✓ Nature journal matplotlib settings configured")

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
    print("📊 CREATING PUBLICATION-READY COMBINED BOX PLOT")
    print("="*60)

    nuclei_fig_path = None
    try:
        nuclei_fig_path = create_nuclei_count_boxplots(sample_df, output_dir)
        print(f"✅ Combined box plot analysis completed successfully!")

    except Exception as e:
        print(f"❌ Error in combined nuclei count analysis: {e}")
        import traceback
        traceback.print_exc()

    # Create UMAP genes plot
    print("\n" + "="*60)
    print("📊 CREATING PUBLICATION-READY UMAP")
    print("="*60)

    umap_fig_path = None
    try:
        umap_fig_path = create_umap_genes_plot(adata, output_dir)
        print(f"✅ UMAP analysis completed successfully!")

    except Exception as e:
        print(f"❌ Error in UMAP analysis: {e}")
        import traceback
        traceback.print_exc()

    # Final summary
    print(f"\n" + "="*60)
    print("✅ PUBLICATION-READY ANALYSIS COMPLETED")
    print("="*60)

    if nuclei_fig_path or umap_fig_path:
        print("📁 OUTPUT FILES (Nature journal ready):")
        if nuclei_fig_path:
            print(f"  Combined box plot:")
            print(f"    • {nuclei_fig_path}.png (raster, 300 DPI)")
            print(f"    • {nuclei_fig_path}.pdf (vector)")
            print(f"    • {nuclei_fig_path}.eps (vector, Nature preferred)")
        if umap_fig_path:
            print(f"  UMAP:")
            print(f"    • {umap_fig_path}.png (raster, 300 DPI)")
            print(f"    • {umap_fig_path}.pdf (vector)")
            print(f"    • {umap_fig_path}.eps (vector, Nature preferred)")

        print(f"\n📊 FIGURE SPECIFICATIONS:")
        print(f"  • Font: Arial/Helvetica, 8pt body, 9pt titles")
        print(f"  • Colors: Colorblind-friendly palette")
        print(f"  • Format: EPS/PDF vector + 300 DPI PNG")
        print(f"  • Style: Nature journal standards")
        print(f"  • Statistics: Non-parametric tests with p-values")

    else:
        print("❌ No figures were generated successfully")

    print("\n🎉 Ready for Nature journal submission!")

if __name__ == "__main__":
    main()
