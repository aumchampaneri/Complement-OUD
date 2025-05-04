import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.patches import Patch
import plotly.graph_objects as go
import matplotlib.gridspec as gridspec
from matplotlib_venn import venn2


# Load the predictions
df = pd.read_csv("/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/DL-PPI outputs/full_ppi_predictions.csv")


# 1. NETWORK VISUALIZATION
def create_ppi_network(df, min_score=0.5):
    # Filter by prediction score if desired
    filtered_df = df[df['prediction'] >= min_score]

    # Check if we have any interactions after filtering
    if len(filtered_df) == 0:
        print(f"No interactions found with score >= {min_score}. Try lowering the threshold.")
        return

    # Create a graph
    G = nx.Graph()

    # Add edges with prediction score as weight
    for _, row in filtered_df.iterrows():
        G.add_edge(row['gene1'], row['gene2'],
                   weight=row['prediction'],
                   color=row['prediction'])  # Store the value, not the color object

    # Create figure with proper axes
    fig, ax = plt.subplots(figsize=(12, 10))

    # Position nodes using force-directed layout
    pos = nx.spring_layout(G, k=0.2, iterations=50)

    # Draw edges with colors based on weight
    edges = nx.draw_networkx_edges(G, pos, ax=ax,
                                   edge_color=[G[u][v]['color'] for u, v in G.edges()],
                                   width=[G[u][v]['weight'] * 2 for u, v in G.edges()],
                                   edge_cmap=plt.cm.plasma)

    # Only include nodes that exist in the graph
    comp_genes = [node for node in G.nodes() if node in set(df['gene1'])]
    sig_genes = [node for node in G.nodes() if node in set(df['gene2'])]

    nx.draw_networkx_nodes(G, pos, ax=ax,
                           nodelist=comp_genes,
                           node_color='lightblue',
                           node_size=300,
                           alpha=0.8)

    nx.draw_networkx_nodes(G, pos, ax=ax,
                           nodelist=sig_genes,
                           node_color='lightgreen',
                           node_size=200,
                           alpha=0.8)

    # Add labels to nodes
    nx.draw_networkx_labels(G, pos, font_size=8, ax=ax)

    ax.set_title(f"Protein-Protein Interaction Network (score ≥ {min_score})")
    ax.axis('off')

    # Create a colorbar with proper normalization
    from matplotlib.colors import Normalize
    norm = Normalize(vmin=min_score, vmax=1.0)
    sm = plt.cm.ScalarMappable(cmap=plt.cm.plasma, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label="Interaction Score")

    plt.tight_layout()
    plt.savefig("DL-PPI outputs/ppi_network.png", dpi=300)
    plt.show()


# 2. HEATMAP VISUALIZATION
def create_ppi_heatmap(df, max_proteins=50):
    # If there are too many proteins, select top interactions
    if len(set(df['gene1'])) * len(set(df['gene2'])) > max_proteins ** 2:
        # Sort by prediction score and take top interactions
        filtered_df = df.sort_values('prediction', ascending=False)
        filtered_df = filtered_df.head(max_proteins ** 2)
    else:
        filtered_df = df.copy()

    # Create a pivot table for the heatmap
    pivot = filtered_df.pivot_table(index='gene1', columns='gene2',
                                    values='prediction', aggfunc='max')

    # Plot the heatmap
    plt.figure(figsize=(12, 10))
    sns.heatmap(pivot, cmap='viridis', linewidths=0.5,
                annot=False, cbar_kws={'label': 'Interaction Score'})
    plt.title("PPI Prediction Scores")
    plt.tight_layout()
    plt.savefig("DL-PPI outputs/ppi_heatmap.png", dpi=300)
    plt.show()


# Call the visualization functions
create_ppi_network(df, min_score=0.6)  # Adjust min_score as needed
create_ppi_heatmap(df)


# 3. TOP INTERACTIONS TABLE
def show_top_interactions(df, n=20):
    # Sort by prediction score
    top_df = df.sort_values('prediction', ascending=False).head(n)

    # Select relevant columns
    top_df = top_df[['gene1', 'gene2', 'prediction', 'method']]

    # Display the top interactions
    print(f"Top {n} predicted protein interactions:")
    return top_df


# 1. CIRCULAR LAYOUT NETWORK
def create_circular_ppi_network(df, min_score=0.5):
    """Creates a circular layout network visualization."""
    filtered_df = df[df['prediction'] >= min_score]

    if len(filtered_df) == 0:
        print(f"No interactions found with score >= {min_score}.")
        return

    G = nx.Graph()
    for _, row in filtered_df.iterrows():
        G.add_edge(row['gene1'], row['gene2'], weight=row['prediction'])

    # Create circular layout
    fig, ax = plt.subplots(figsize=(12, 12))
    pos = nx.circular_layout(G)

    # Draw network
    comp_genes = [n for n in G.nodes() if n in set(filtered_df['gene1'])]
    sig_genes = [n for n in G.nodes() if n in set(filtered_df['gene2'])]

    # Edges with width proportional to weight
    nx.draw_networkx_edges(G, pos, alpha=0.3, width=[G[u][v]['weight'] * 3 for u, v in G.edges()])

    # Draw nodes
    nx.draw_networkx_nodes(G, pos, nodelist=comp_genes, node_color='dodgerblue',
                           node_size=200, label='Complement genes')
    nx.draw_networkx_nodes(G, pos, nodelist=sig_genes, node_color='tomato',
                           node_size=200, label='DEG genes')

    # Add labels for only the most important nodes (option to reduce cluttering)
    important_nodes = [n for n in G.nodes() if G.degree(n) > 2]
    nx.draw_networkx_labels(G, pos, {n: n for n in important_nodes}, font_size=8)

    plt.legend()
    plt.title(f"Circular PPI Network (score ≥ {min_score})")
    plt.axis('off')
    plt.tight_layout()
    plt.savefig("DL-PPI outputs/circular_network.png", dpi=300)
    plt.show()


# 2. CONNECTION DEGREE ANALYSIS
def analyze_node_degrees(df, min_score=0.5):
    """Shows the most connected proteins in the network."""
    filtered_df = df[df['prediction'] >= min_score]

    # Build graph
    G = nx.Graph()
    for _, row in filtered_df.iterrows():
        G.add_edge(row['gene1'], row['gene2'], weight=row['prediction'])

    # Get degree for each node
    degree_dict = dict(G.degree())

    # Sort by degree
    sorted_degrees = sorted(degree_dict.items(), key=lambda x: x[1], reverse=True)

    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # Top connected complement genes
    comp_genes_degrees = [(node, deg) for node, deg in sorted_degrees if node in set(df['gene1'])]
    comp_df = pd.DataFrame(comp_genes_degrees[:15], columns=['Gene', 'Connections'])
    comp_df.plot(kind='barh', x='Gene', y='Connections', ax=ax1, color='dodgerblue')
    ax1.set_title('Top Connected Complement Genes')

    # Top connected significant genes
    sig_genes_degrees = [(node, deg) for node, deg in sorted_degrees if node in set(df['gene2'])]
    sig_df = pd.DataFrame(sig_genes_degrees[:15], columns=['Gene', 'Connections'])
    sig_df.plot(kind='barh', x='Gene', y='Connections', ax=ax2, color='tomato')
    ax2.set_title('Top Connected DEG Genes')

    plt.tight_layout()
    plt.savefig("DL-PPI outputs/degree_analysis.png", dpi=300)
    plt.show()

    return comp_df, sig_df


# 3. INTERACTION SCORE DISTRIBUTION
def plot_score_distribution(df):
    """Shows the distribution of interaction scores."""
    plt.figure(figsize=(10, 6))

    # Plot histogram with kernel density estimate
    sns.histplot(df['prediction'], kde=True, bins=20)

    plt.axvline(x=0.5, color='red', linestyle='--', label='Threshold = 0.5')
    plt.axvline(x=0.7, color='green', linestyle='--', label='Threshold = 0.7')

    plt.title('Distribution of PPI Prediction Scores')
    plt.xlabel('Prediction Score')
    plt.ylabel('Frequency')
    plt.legend()
    plt.tight_layout()
    plt.savefig("DL-PPI outputs/score_distribution.png", dpi=300)
    plt.show()


# 4. COMMUNITY DETECTION
def detect_communities(df, min_score=0.6):
    """Detects and visualizes communities within the network."""
    filtered_df = df[df['prediction'] >= min_score]

    if len(filtered_df) == 0:
        print(f"No interactions found with score >= {min_score}.")
        return

    G = nx.Graph()
    for _, row in filtered_df.iterrows():
        G.add_edge(row['gene1'], row['gene2'], weight=row['prediction'])

    # Use Louvain method for community detection
    try:
        import community as community_louvain
        partition = community_louvain.best_partition(G)

        # Create figure
        plt.figure(figsize=(12, 10))

        # Create position layout
        pos = nx.spring_layout(G, k=0.3)

        # Color nodes according to communities
        cmap = plt.cm.get_cmap("tab20", max(partition.values()) + 1)
        nx.draw_networkx_nodes(G, pos, partition.keys(),
                               node_size=200,
                               cmap=cmap,
                               node_color=list(partition.values()))

        # Draw edges
        nx.draw_networkx_edges(G, pos, alpha=0.5)

        # Add labels to nodes
        nx.draw_networkx_labels(G, pos, font_size=8)

        plt.title(f"Community Detection in PPI Network (score ≥ {min_score})")
        plt.axis('off')
        plt.tight_layout()
        plt.savefig("DL-PPI outputs/community_detection.png", dpi=300)
        plt.show()

        return partition
    except ImportError:
        print("Please install python-louvain: pip install python-louvain")
        return None


# Call the new visualization functions
create_circular_ppi_network(df, min_score=0.5)
analyze_node_degrees(df, min_score=0.5)
plot_score_distribution(df)
# Optional: detect_communities(df, min_score=0.5) - requires python-louvain

top_interactions = show_top_interactions(df)
print(top_interactions)
top_interactions.to_csv("DL-PPI outputs/top_interactions.csv", index=False)