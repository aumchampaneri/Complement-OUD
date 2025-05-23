import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.patches import Patch
import plotly.graph_objects as go
import matplotlib.gridspec as gridspec
from matplotlib_venn import venn2
from matplotlib.colors import Normalize
from plotly.offline import plot
import plotly.express as px



# Load the predictions
df = pd.read_csv("/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158 All Cells/DL-PPI outputs/full_ppi_predictions.csv")


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
detect_communities(df, min_score=0.5)

top_interactions = show_top_interactions(df)
print(top_interactions)
top_interactions.to_csv("DL-PPI outputs/top_interactions.csv", index=False)


def interactive_ppi_network(df, min_score=0.5):
    """Creates an interactive network visualization using Plotly."""
    filtered_df = df[df['prediction'] >= min_score]
    if len(filtered_df) == 0:
        print(f"No interactions found with score >= {min_score}.")
        return

    # Create graph
    G = nx.Graph()
    for _, row in filtered_df.iterrows():
        G.add_edge(row['gene1'], row['gene2'], weight=row['prediction'])

    # Get positions for nodes
    pos = nx.spring_layout(G)

    # Create edge traces
    edge_x = []
    edge_y = []
    edge_text = []
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
        # Add edge weight to hover text (3 times because of the None)
        weight = G[edge[0]][edge[1]]['weight']
        edge_text.extend([f"{edge[0]} - {edge[1]} (Score: {weight:.3f})"] * 3)

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='text',
        mode='lines',
        text=edge_text)

    # Create node traces for complement genes
    comp_genes = [node for node in G.nodes() if node in set(df['gene1'])]

    comp_x = []
    comp_y = []
    comp_text = []
    for node in comp_genes:
        x, y = pos[node]
        comp_x.append(x)
        comp_y.append(y)
        comp_text.append(f"Gene: {node}<br>Type: Complement<br>Connections: {G.degree[node]}")

    comp_trace = go.Scatter(
        x=comp_x, y=comp_y,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            size=10,
            color='dodgerblue',
            line_width=2),
        text=comp_text,
        name='Complement genes')

    # Create node traces for significant genes
    sig_genes = [node for node in G.nodes() if node in set(df['gene2'])]

    sig_x = []
    sig_y = []
    sig_text = []
    for node in sig_genes:
        x, y = pos[node]
        sig_x.append(x)
        sig_y.append(y)
        sig_text.append(f"Gene: {node}<br>Type: DEG<br>Connections: {G.degree[node]}")

    sig_trace = go.Scatter(
        x=sig_x, y=sig_y,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            size=10,
            color='tomato',
            line_width=2),
        text=sig_text,
        name='DEG genes')

    # Create figure
    fig = go.Figure(data=[edge_trace, comp_trace, sig_trace],
                    layout=go.Layout(
                        title=f"Interactive PPI Network (score ≥ {min_score})",
                        showlegend=True,
                        hovermode='closest',
                        margin=dict(b=20, l=5, r=5, t=40),
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)))

    # Save the interactive plot as an HTML file
    plot(fig, filename="DL-PPI outputs/interactive_network.html", auto_open=False)
    fig.show()
    return fig


def create_3d_ppi_network(df, min_score=0.5):
    """Creates a 3D network visualization using Plotly."""
    filtered_df = df[df['prediction'] >= min_score]
    if len(filtered_df) == 0:
        print(f"No interactions found with score >= {min_score}.")
        return

    # Create graph
    G = nx.Graph()
    for _, row in filtered_df.iterrows():
        G.add_edge(row['gene1'], row['gene2'], weight=row['prediction'])

    # Get 3D positions
    pos = nx.spring_layout(G, dim=3)

    # Create edge traces
    edge_x = []
    edge_y = []
    edge_z = []
    edge_text = []

    for edge in G.edges():
        x0, y0, z0 = pos[edge[0]]
        x1, y1, z1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
        edge_z.extend([z0, z1, None])

        weight = G[edge[0]][edge[1]]['weight']
        edge_text.extend([f"{edge[0]} - {edge[1]} (Score: {weight:.3f})"] * 3)

    edge_trace = go.Scatter3d(
        x=edge_x, y=edge_y, z=edge_z,
        line=dict(width=2, color='#888'),
        hoverinfo='text',
        text=edge_text,
        mode='lines')

    # Create node traces for complement genes
    comp_genes = [node for node in G.nodes() if node in set(df['gene1'])]

    comp_x = []
    comp_y = []
    comp_z = []
    comp_text = []

    for node in comp_genes:
        x, y, z = pos[node]
        comp_x.append(x)
        comp_y.append(y)
        comp_z.append(z)
        comp_text.append(f"Gene: {node}<br>Type: Complement<br>Connections: {G.degree[node]}")

    comp_trace = go.Scatter3d(
        x=comp_x, y=comp_y, z=comp_z,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            size=5,
            color='dodgerblue',
            line_width=1),
        text=comp_text,
        name='Complement genes')

    # Create node traces for significant genes
    sig_genes = [node for node in G.nodes() if node in set(df['gene2'])]

    sig_x = []
    sig_y = []
    sig_z = []
    sig_text = []

    for node in sig_genes:
        x, y, z = pos[node]
        sig_x.append(x)
        sig_y.append(y)
        sig_z.append(z)
        sig_text.append(f"Gene: {node}<br>Type: DEG<br>Connections: {G.degree[node]}")

    sig_trace = go.Scatter3d(
        x=sig_x, y=sig_y, z=sig_z,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            size=5,
            color='tomato',
            line_width=1),
        text=sig_text,
        name='DEG genes')

    # Create figure
    fig = go.Figure(data=[edge_trace, comp_trace, sig_trace],
                    layout=go.Layout(
                        title=f"3D PPI Network (score ≥ {min_score})",
                        showlegend=True,
                        scene=dict(
                            xaxis=dict(showticklabels=False),
                            yaxis=dict(showticklabels=False),
                            zaxis=dict(showticklabels=False),
                        )))

    # Save the interactive plot
    plot(fig, filename="DL-PPI outputs/3d_network.html", auto_open=False)
    fig.show()
    return fig


def create_sankey_diagram(df, min_score=0.7, max_nodes=30):
    """Creates a Sankey diagram showing gene interactions."""
    # We need a higher threshold for Sankey to be readable
    filtered_df = df[df['prediction'] >= min_score]

    if len(filtered_df) == 0:
        print(f"No interactions found with score >= {min_score}.")
        return

    # Limit number of nodes if needed
    if len(set(filtered_df['gene1']).union(set(filtered_df['gene2']))) > max_nodes:
        filtered_df = filtered_df.sort_values('prediction', ascending=False)
        filtered_df = filtered_df.head(max_nodes)

    # Get all unique nodes
    all_nodes = list(set(filtered_df['gene1']).union(set(filtered_df['gene2'])))
    node_map = {node: i for i, node in enumerate(all_nodes)}

    # Create source and target indices
    sources = [node_map[src] for src in filtered_df['gene1']]
    targets = [node_map[tgt] for tgt in filtered_df['gene2']]

    # Use prediction scores as values
    values = filtered_df['prediction'].tolist()

    # Determine node colors based on type
    node_colors = []
    for node in all_nodes:
        if node in set(df['gene1']) and node in set(df['gene2']):
            node_colors.append('rgba(100, 200, 100, 0.8)')  # Green for both types
        elif node in set(df['gene1']):
            node_colors.append('rgba(65, 105, 225, 0.8)')  # Blue for complement
        else:
            node_colors.append('rgba(255, 99, 71, 0.8)')  # Red for DEG

    # Create Sankey diagram
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=all_nodes,
            color=node_colors
        ),
        link=dict(
            source=sources,
            target=targets,
            value=values,
            color=[f'rgba(0, 0, 255, {val / 2 + 0.3})' for val in values]
        ))])

    fig.update_layout(
        title_text=f"PPI Sankey Diagram (score ≥ {min_score})",
        font_size=10
    )

    # Save the interactive plot
    plot(fig, filename="DL-PPI outputs/sankey_diagram.html", auto_open=False)
    fig.show()
    return fig


def plot_network_metrics(df, min_score=0.5):
    """Analyzes and visualizes various network metrics."""
    filtered_df = df[df['prediction'] >= min_score]

    if len(filtered_df) == 0:
        print(f"No interactions found with score >= {min_score}.")
        return

    # Create graph
    G = nx.Graph()
    for _, row in filtered_df.iterrows():
        G.add_edge(row['gene1'], row['gene2'], weight=row['prediction'])

    # Calculate network metrics
    metrics = {
        'Degree Centrality': nx.degree_centrality(G),
        'Betweenness Centrality': nx.betweenness_centrality(G),
        'Closeness Centrality': nx.closeness_centrality(G),
        'Eigenvector Centrality': nx.eigenvector_centrality(G, max_iter=1000)
    }

    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    axes = axes.flatten()

    # Plot each metric
    for i, (metric_name, metric_values) in enumerate(metrics.items()):
        # Sort nodes by metric value
        sorted_nodes = sorted(metric_values.items(), key=lambda x: x[1], reverse=True)[:15]
        nodes = [n[0] for n in sorted_nodes]
        values = [n[1] for n in sorted_nodes]

        # Determine node colors
        colors = ['dodgerblue' if node in set(df['gene1']) else 'tomato' for node in nodes]

        # Create bar chart
        axes[i].bar(range(len(nodes)), values, color=colors)
        axes[i].set_xticks(range(len(nodes)))
        axes[i].set_xticklabels(nodes, rotation=45, ha='right')
        axes[i].set_title(f'Top Genes by {metric_name}')
        axes[i].set_ylabel(metric_name)

        # Add a horizontal line for the average
        avg = np.mean(list(metric_values.values()))
        axes[i].axhline(y=avg, color='green', linestyle='--', alpha=0.7,
                        label=f'Average: {avg:.3f}')
        axes[i].legend()

    plt.tight_layout()
    plt.savefig("DL-PPI outputs/network_metrics.png", dpi=300)
    plt.show()

    return metrics


# Update function calls at the end of the script
# Keep your existing function calls and add:

# Run advanced visualization functions
print("\nGenerating interactive PPI network...")
interactive_ppi_network(df, min_score=0.5)

print("\nGenerating 3D PPI network...")
create_3d_ppi_network(df, min_score=0.5)

print("\nGenerating Sankey diagram...")
create_sankey_diagram(df, min_score=0.7)  # Higher threshold for readability

print("\nAnalyzing network metrics...")
plot_network_metrics(df, min_score=0.5)