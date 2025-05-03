# GSE225158/4_PPI_Network.py
import pandas as pd
import numpy as np
import requests
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import os


# 1. Load DEG results
def load_deseq_results(file_path):
    """Load differential expression results from DESeq2 analysis"""
    deg_results = pd.read_csv(file_path)

    if deg_results.empty:
        raise ValueError(f"The DESeq2 results file at {file_path} is empty.")

    required_columns = {'log2FoldChange', 'padj', 'gene'}
    if not required_columns.issubset(deg_results.columns):
        raise ValueError(
            f"The DESeq2 results file at {file_path} is missing required columns: {required_columns - set(deg_results.columns)}")

    return deg_results


# 2. Load complement genes from CSV file
def load_complement_genes(file_path):
    """Load complement genes from CSV file"""
    complement_df = pd.read_csv(file_path)

    if complement_df.empty:
        raise ValueError(f"The complement genes file at {file_path} is empty.")

    # If the CSV has a column named 'Gene' or similar
    if 'Gene' in complement_df.columns:
        complement_genes = complement_df['Gene'].dropna().tolist()
    else:
        complement_genes = complement_df.iloc[:, 0].dropna().tolist()

    if not complement_genes:
        raise ValueError(f"No valid genes found in the complement genes file at {file_path}.")

    return complement_genes


# 3. Get protein interactions from STRING
def get_string_interactions(genes, species=9606, score_threshold=700):
    """Fetch protein-protein interactions from STRING database"""
    string_api_url = "https://string-db.org/api"
    output_format = "json"
    method = "network"

    request_url = "/".join([string_api_url, output_format, method])

    params = {
        "identifiers": "%0d".join(genes),
        "species": species,
        "caller_identity": "complement_oud_analysis",
        "required_score": score_threshold
    }

    response = requests.post(request_url, data=params)

    if response.status_code != 200:
        raise ConnectionError(f"STRING API request failed with status code {response.status_code}: {response.text}")

    try:
        return response.json()
    except ValueError:
        raise ValueError("Failed to parse STRING API response as JSON.")


# 4. Create and analyze network
def create_ppi_network(ppi_data, deg_data=None):
    """Create NetworkX graph from STRING PPI data and analyze"""
    G = nx.Graph(name="Protein Interaction Network")

    # Add edges with confidence scores
    for interaction in ppi_data:
        source = interaction['preferredName_A']
        target = interaction['preferredName_B']
        score = float(interaction['score'])
        G.add_edge(source, target, weight=score / 1000.0)

    # Calculate network metrics
    degree_dict = dict(G.degree())
    betweenness_dict = nx.betweenness_centrality(G)

    # Create node attribute dict with expression data if available
    node_attr = {}
    if deg_data is not None:
        for node in G.nodes():
            if node in deg_data['gene'].values:
                idx = deg_data[deg_data['gene'] == node].index[0]
                node_attr[node] = {
                    'log2FC': deg_data.loc[idx, 'log2FoldChange'],
                    'padj': deg_data.loc[idx, 'padj'],
                    'degree': degree_dict[node],
                    'betweenness': betweenness_dict[node]
                }
            else:
                node_attr[node] = {
                    'log2FC': 0,
                    'padj': 1,
                    'degree': degree_dict[node],
                    'betweenness': betweenness_dict[node]
                }

        # Set node attributes
        nx.set_node_attributes(G, node_attr)

    return G


# 5. Visualize the network
def visualize_network(G, output_path):
    """Create and save network visualization"""
    plt.figure(figsize=(12, 12))

    # Set up layout
    pos = nx.spring_layout(G, seed=42)

    # Get node attributes for visualization
    node_size = [300 * (0.1 + G.nodes[n].get('degree', 1)) for n in G.nodes()]

    # Create custom colormap for log2FC (blue-white-red)
    custom_cmap = LinearSegmentedColormap.from_list(
        "blue_white_red",
        [(0, "blue"), (0.5, "white"), (1, "red")],
    )

    # Node colors based on log2FC if available
    if 'log2FC' in G.nodes[list(G.nodes())[0]]:
        node_colors = [G.nodes[n].get('log2FC', 0) for n in G.nodes()]
        vmin, vmax = -2, 2  # Range for log2FC visualization
        nx.draw_networkx_nodes(G, pos, node_size=node_size, node_color=node_colors,
                               cmap=custom_cmap, vmin=vmin, vmax=vmax)

        # Add colorbar
        sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=plt.gca())
        cbar.set_label('log2 Fold Change')
    else:
        nx.draw_networkx_nodes(G, pos, node_size=node_size)

    # Draw edges with transparency based on weight
    edge_alphas = [G[u][v].get('weight', 0.1) for u, v in G.edges()]
    nx.draw_networkx_edges(G, pos, alpha=edge_alphas, edge_color='gray')

    # Draw labels for nodes with high degree
    high_degree_nodes = {n: n for n in G.nodes() if G.nodes[n].get('degree', 0) > 5}
    nx.draw_networkx_labels(G, pos, labels=high_degree_nodes, font_size=10)

    plt.title("Protein-Protein Interaction Network")
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


# 6. Identify hub genes and save results
def identify_hubs(G, output_path, top_n=20):
    """Identify hub genes and save results to file"""
    # Create DataFrame with node metrics
    nodes_df = pd.DataFrame([
        {
            'Gene': node,
            'Degree': G.nodes[node].get('degree', 0),
            'Betweenness': G.nodes[node].get('betweenness', 0),
            'log2FC': G.nodes[node].get('log2FC', None),
            'padj': G.nodes[node].get('padj', None)
        }
        for node in G.nodes()
    ])

    # Sort by degree (hub score)
    nodes_df = nodes_df.sort_values('Degree', ascending=False)

    # Save results
    nodes_df.to_csv(output_path, index=False)

    return nodes_df.head(top_n)


# Main execution
def main():
    # Create output directory if it doesn't exist
    output_dir = "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/PPI outputs"
    os.makedirs(output_dir, exist_ok=True)

    # 1. Load DESeq2 results
    deg_results = load_deseq_results(
        "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/DESeq2 outputs/deseq2_results_M_OUD_vs_M_None.csv")

    # 2. Filter for significant DEGs
    sig_genes = deg_results[(deg_results['padj'] < 0.05) &
                            (abs(deg_results['log2FoldChange']) > 1)].copy()
    print(f"Found {len(sig_genes)} significant DEGs")

    # 3. Load complement genes from CSV
    complement_genes = load_complement_genes(
        "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/KEGG outputs/kegg_complement_unique_genes.csv")
    print(f"Loaded {len(complement_genes)} complement genes from file")

    # 4. Combine gene lists (DEGs + complement genes)
    # Using set to remove duplicates
    genes_to_analyze = list(set(sig_genes['gene'].tolist() + complement_genes))

    # If the list is too large for API, take top genes by significance
    if len(genes_to_analyze) > 400:
        sig_genes = sig_genes.sort_values('padj').head(380)
        genes_to_analyze = list(set(sig_genes['gene'].tolist() + complement_genes))

    print(f"Analyzing PPI network for {len(genes_to_analyze)} genes")

    # 5. Get PPI data from STRING
    ppi_data = get_string_interactions(genes_to_analyze)
    print(f"Retrieved {len(ppi_data)} interactions")

    # 6. Create network
    G = create_ppi_network(ppi_data, deg_results)
    print(f"Created network with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")

    # 7. Visualize network
    visualize_network(G, f"{output_dir}/PPI_network.png")

    # 8. Identify and save hub genes
    top_hubs = identify_hubs(G, f"{output_dir}/network_metrics.csv")
    print("Top hub genes in the network:")
    print(top_hubs)


if __name__ == "__main__":
    main()