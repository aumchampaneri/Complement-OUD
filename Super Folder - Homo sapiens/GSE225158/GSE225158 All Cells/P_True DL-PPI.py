import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap, Normalize, ListedColormap
from matplotlib.cm import ScalarMappable
import community as community_louvain
from sklearn.metrics import roc_curve, auc
from itertools import cycle

# Create output directory
output_dir = "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158 All Cells/DL-PPI outputs/visualizations"
os.makedirs(output_dir, exist_ok=True)

# Load data - skip both header rows
male_csv = "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158 All Cells/DL-PPI outputs/dl_ppi_viz_data_male_0.4_genes.csv"
female_csv = "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158 All Cells/DL-PPI outputs/dl_ppi_viz_data_female_0.4_genes.csv"

# Load data properly, skipping both problematic header rows
male_df = pd.read_csv(male_csv,
                     skiprows=2,  # Skip both header rows
                     names=["protein1", "protein2", "prediction", "method"])

female_df = pd.read_csv(female_csv,
                       skiprows=2,  # Skip both header rows
                       names=["protein1", "protein2", "prediction", "method"])

print(f"Loaded {len(male_df)} male and {len(female_df)} female interactions")

# 1. COMBINED NETWORK VISUALIZATION WITH SEX-SPECIFIC COLORING
print("Creating combined network visualization...")
G_combined = nx.Graph()

# Add male edges (blue)
for _, row in male_df.iterrows():
    G_combined.add_edge(row['protein1'], row['protein2'],
                        weight=float(row['prediction']),
                        male=1, female=0,
                        color='blue')

# Add/update female edges (red or green if shared)
for _, row in female_df.iterrows():
    if G_combined.has_edge(row['protein1'], row['protein2']):
        # Edge exists in both networks - update to shared
        G_combined[row['protein1']][row['protein2']]['female'] = 1
        G_combined[row['protein1']][row['protein2']]['color'] = 'green'
    else:
        # Female-specific edge
        G_combined.add_edge(row['protein1'], row['protein2'],
                            weight=float(row['prediction']),
                            male=0, female=1,
                            color='red')

plt.figure(figsize=(16, 14), dpi=300)

# Calculate node sizes based on degree centrality
node_sizes = [G_combined.degree(n) * 50 + 100 for n in G_combined.nodes()]

# Position using force-directed layout
try:
    pos = nx.spring_layout(G_combined, seed=42, k=0.3)
except:
    print("Spring layout failed, trying Kamada-Kawai layout...")
    pos = nx.kamada_kawai_layout(G_combined)

# Draw edges first, colored by sex-specificity
edge_colors = [G_combined[u][v]['color'] for u, v in G_combined.edges()]
nx.draw_networkx_edges(G_combined, pos, width=1.5, alpha=0.6, edge_color=edge_colors)

# Draw nodes
nx.draw_networkx_nodes(G_combined, pos,
                       node_size=node_sizes,
                       node_color='lightgray',
                       edgecolors='black',
                       linewidths=0.5,
                       alpha=0.8)

# Label only hub nodes (top 15% by degree)
degree_threshold = np.percentile([G_combined.degree(n) for n in G_combined.nodes()], 85)
hub_labels = {n: n for n in G_combined.nodes() if G_combined.degree(n) >= degree_threshold}
nx.draw_networkx_labels(G_combined, pos, labels=hub_labels, font_size=10, font_weight='bold')

# Create legend with correct colors
male_line = plt.Line2D([0], [0], color='blue', lw=2, label='Male-specific')
female_line = plt.Line2D([0], [0], color='red', lw=2, label='Female-specific')
both_line = plt.Line2D([0], [0], color='green', lw=2, label='Shared')
plt.legend(handles=[male_line, female_line, both_line], loc='lower right', fontsize=12)

# Network stats
male_only = sum(1 for _, _, d in G_combined.edges(data=True) if d['color'] == 'blue')
female_only = sum(1 for _, _, d in G_combined.edges(data=True) if d['color'] == 'red')
shared = sum(1 for _, _, d in G_combined.edges(data=True) if d['color'] == 'green')

# Add stats box
stats_text = (
    f"Network Comparison:\n"
    f"• Total proteins: {len(G_combined.nodes())}\n"
    f"• Male-specific interactions: {male_only}\n"
    f"• Female-specific interactions: {female_only}\n"
    f"• Shared interactions: {shared}\n"
)
plt.figtext(0.02, 0.02, stats_text, bbox=dict(facecolor='white', alpha=0.8), fontsize=12)

plt.title("Sex Differences in OUD Complement System Interactions", fontsize=18)
plt.axis('off')
plt.tight_layout()
plt.savefig(f"{output_dir}/sex_comparison_network.png", dpi=300, bbox_inches='tight')
plt.close()

# 2. HUB PROTEIN COMPARISON
print("Creating hub protein comparison...")
def get_top_proteins(df, n=15):
    all_proteins = pd.concat([df['protein1'], df['protein2']]).value_counts()
    return all_proteins.head(n).index.tolist()

male_top = get_top_proteins(male_df)
female_top = get_top_proteins(female_df)
all_top = list(set(male_top + female_top))

# Create comparison dataframe
comparison_data = []
for protein in all_top:
    male_connections = sum((male_df['protein1'] == protein) | (male_df['protein2'] == protein))
    female_connections = sum((female_df['protein1'] == protein) | (female_df['protein2'] == protein))
    comparison_data.append({
        'Protein': protein,
        'Male Connections': male_connections,
        'Female Connections': female_connections,
        'Difference': female_connections - male_connections
    })

comp_df = pd.DataFrame(comparison_data)
comp_df = comp_df.sort_values('Difference', ascending=False)

# Create lollipop chart for top hub proteins comparison
plt.figure(figsize=(12, 8), dpi=300)

# Calculate y positions
y_pos = np.arange(len(comp_df))

# Create horizontal lines between points
for i, row in enumerate(comp_df.iterrows()):
    plt.plot([row[1]['Male Connections'], row[1]['Female Connections']],
             [i, i], color='gray', linestyle='-', linewidth=1.5, alpha=0.6)

# Male points (blue)
plt.scatter(comp_df['Male Connections'], y_pos, s=80, color='blue', alpha=0.7, label='Male')

# Female points (red)
plt.scatter(comp_df['Female Connections'], y_pos, s=80, color='red', alpha=0.7, label='Female')

# Add protein labels
plt.yticks(y_pos, comp_df['Protein'])
plt.grid(axis='x', linestyle='--', alpha=0.7)

plt.xlabel('Number of Protein Interactions', fontsize=12)
plt.title('Hub Protein Comparison: Male vs Female OUD', fontsize=16)
plt.legend()
plt.tight_layout()
plt.savefig(f"{output_dir}/hub_protein_comparison.png", dpi=300)
plt.close()

# 3. INTERACTION STRENGTH HEATMAP
print("Creating interaction strength comparison...")
# Get shared interactions to compare strength
shared_interactions = []
for _, male_row in male_df.iterrows():
    m_p1, m_p2 = male_row['protein1'], male_row['protein2']
    pair_key = frozenset([m_p1, m_p2])

    # Look for the same interaction in female df
    for _, female_row in female_df.iterrows():
        f_p1, f_p2 = female_row['protein1'], female_row['protein2']
        if frozenset([f_p1, f_p2]) == pair_key:
            shared_interactions.append({
                'Protein1': m_p1,
                'Protein2': m_p2,
                'Male Score': male_row['prediction'],
                'Female Score': female_row['prediction'],
                'Difference': female_row['prediction'] - male_row['prediction']
            })
            break

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
    plt.yticks([i + bar_width/2 for i in index], interaction_labels)
    plt.xlabel('Interaction Confidence Score')
    plt.title('Top Differential Protein Interactions by Sex', fontsize=16)
    plt.legend()
    plt.grid(axis='x', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/differential_interactions.png", dpi=300)
    plt.close()

# 4. NETWORK STATISTICS COMPARISON
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
        sum(dict(G_male.degree()).values())/len(G_male.nodes()),
        nx.number_connected_components(G_male)
    ],
    'Female': [
        len(G_female.nodes()),
        len(G_female.edges()),
        nx.density(G_female),
        sum(dict(G_female.degree()).values())/len(G_female.nodes()),
        nx.number_connected_components(G_female)
    ]
}

# Create a nice comparison table/figure
metrics_df = pd.DataFrame(metrics)
metrics_df['Difference %'] = (metrics_df['Female'] - metrics_df['Male'])/metrics_df['Male']*100

# Create a visual table
plt.figure(figsize=(10, 5), dpi=300)
plt.axis('off')
table = plt.table(
    cellText=metrics_df.values,
    colLabels=metrics_df.columns,
    cellLoc='center',
    loc='center',
    cellColours=[['#f2f2f2']*4]*5
)
table.auto_set_font_size(False)
table.set_fontsize(12)
table.scale(1.2, 1.8)
plt.title('Network Property Comparison: Male vs Female OUD', fontsize=16, pad=20)
plt.savefig(f"{output_dir}/network_metrics_comparison.png", dpi=300, bbox_inches='tight')
plt.close()

# 5. STATISTICAL COMPARISON OF SCORE DISTRIBUTIONS
print("Creating statistical comparison of score distributions...")

# A. Kolmogorov-Smirnov test to statistically compare distributions
from scipy import stats
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

# C. Enhanced histogram/density plot (already implemented in section 6)
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
plt.text(0.42, plt.gca().get_ylim()[1]*0.9,
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