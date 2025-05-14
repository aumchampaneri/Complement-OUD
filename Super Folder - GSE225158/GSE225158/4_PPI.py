# GSE225158/4_PPI_Network.py
import pandas as pd
import numpy as np
import requests
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import os
import warnings
import seaborn as sns
import math
from scipy import stats
import argparse
import io
import logging
import sys
from datetime import datetime

def setup_logging(output_dir, pathway_name=None):
    """Set up logging to both console and file"""
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Create a unique log file name with timestamp, using .txt extension
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"{pathway_name}_analysis_{timestamp}.txt" if pathway_name else f"ppi_analysis_{timestamp}.txt"
    log_path = os.path.join(output_dir, log_filename)

    # Configure logging to write to both console and file
    file_handler = logging.FileHandler(log_path)
    console_handler = logging.StreamHandler(sys.stdout)

    handlers = [file_handler, console_handler]

    logging.basicConfig(
        level=logging.INFO,
        format='%(message)s',
        handlers=handlers
    )

    logging.info(f"Logging to: {log_path}")
    return log_path


def setup_output_capture(output_dir, pathway_name=None):
    """Set up output capture to file using redirection"""
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Create a unique output file name with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{pathway_name}_analysis_{timestamp}.txt" if pathway_name else f"ppi_analysis_{timestamp}.txt"
    filepath = os.path.join(output_dir, filename)

    # Create and open the file for writing
    output_file = open(filepath, 'w')

    # Save original stdout/stderr
    original_stdout = sys.stdout
    original_stderr = sys.stderr

    # Create a custom stdout/stderr that writes to both file and console
    class OutputCapture:
        def __init__(self, file, original):
            self.file = file
            self.original = original

        def write(self, text):
            self.file.write(text)
            self.original.write(text)

        def flush(self):
            self.file.flush()
            self.original.flush()

    # Replace stdout and stderr with our capturing versions
    sys.stdout = OutputCapture(output_file, original_stdout)
    sys.stderr = OutputCapture(output_file, original_stderr)

    print(f"Output capturing started. Saving to: {filepath}")
    return filepath, original_stdout, original_stderr, output_file


def end_output_capture(original_stdout, original_stderr, output_file):
    """Restore original stdout/stderr and close the output file"""
    sys.stdout = original_stdout
    sys.stderr = original_stderr
    output_file.close()
    print("Output capture ended.")

warnings.filterwarnings("ignore", category=UserWarning)


# STRING API base URLs
STRING_API_URL = "https://string-db.org/api"
OUTPUT_FORMAT = "tsv"
METHOD_GET_STRING_IDS = "get_string_ids"
METHOD_NETWORK = "network"


def get_string_ids(genes, species):
    params = {
        "identifiers": "%0d".join(genes),
        "species": species,
        "limit": 1,
        "echo_query": 1
    }
    request_url = f"{STRING_API_URL}/{OUTPUT_FORMAT}/{METHOD_GET_STRING_IDS}"
    try:
        response = requests.post(request_url, data=params)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        logging.error(f"Error fetching STRING IDs: {e}")
        return None

    df = pd.read_csv(io.StringIO(response.text), sep="\t")
    if df.empty:
        logging.warning("No STRING IDs were found for the provided genes.")
    return df


def get_interactions(string_ids, species, required_score):
    identifiers = "%0d".join(string_ids)
    params = {
        "identifiers": identifiers,
        "species": species,
        "required_score": required_score
    }
    request_url = f"{STRING_API_URL}/{OUTPUT_FORMAT}/{METHOD_NETWORK}"
    try:
        response = requests.post(request_url, data=params)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        logging.error(f"Error fetching interaction network: {e}")
        return None

    df = pd.read_csv(io.StringIO(response.text), sep="\t")
    if df.empty:
        logging.warning("No interactions were retrieved.")
    return df


def plot_score_distribution(df, output_file):
    plt.hist(df["score"], bins=50, color='skyblue', edgecolor='black')
    plt.title("STRING Interaction Score Distribution")
    plt.xlabel("Score")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(output_file)
    print(f"Histogram saved to {output_file}")


def export_graphml(df, output_file):
    G = nx.Graph()
    for _, row in df.iterrows():
        G.add_edge(row['preferredName_A'], row['preferredName_B'], weight=row['score'])
    nx.write_graphml(G, output_file)
    print(f"GraphML network exported to {output_file}")


def plot_interaction_network(df, output_file):
    G = nx.Graph()
    for _, row in df.iterrows():
        G.add_edge(row['preferredName_A'], row['preferredName_B'], weight=row['score'])

    plt.figure(figsize=(10, 10))
    pos = nx.spring_layout(G, k=0.5)
    nx.draw(G, pos, with_labels=True, node_size=50, font_size=8, edge_color='gray', node_color='skyblue')
    plt.title("STRING PPI Network")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    print(f"Network plot saved to {output_file}")

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


# 2. Load pathway genes from CSV file
def load_pathway_genes(file_path):
    """Load pathway genes from CSV file"""
    pathway_df = pd.read_csv(file_path)

    if pathway_df.empty:
        raise ValueError(f"The pathway genes file at {file_path} is empty.")

    # If the CSV has a column named 'Gene' or similar
    if 'Gene' in pathway_df.columns:
        pathway_genes = pathway_df['Gene'].dropna().tolist()
    elif 'gene' in pathway_df.columns:
        pathway_genes = pathway_df['gene'].dropna().tolist()
    else:
        pathway_genes = pathway_df.iloc[:, 0].dropna().tolist()

    if not pathway_genes:
        raise ValueError(f"No valid genes found in the pathway genes file at {file_path}.")

    return pathway_genes


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
        "caller_identity": "ppi_network_analysis",
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


# 7. Detect communities in the network
def detect_communities(G, output_path):
    """Detect communities in the network and visualize them"""
    from networkx.algorithms.community import greedy_modularity_communities

    # Find communities
    communities = list(greedy_modularity_communities(G))
    print(f"Detected {len(communities)} communities")

    # Assign community to each node
    node_community = {}
    for i, community in enumerate(communities):
        for node in community:
            node_community[node] = i

    nx.set_node_attributes(G, node_community, 'community')

    # Visualize with communities as colors
    plt.figure(figsize=(12, 12))
    pos = nx.spring_layout(G, seed=42)

    # Color nodes by community
    node_colors = [G.nodes[n]['community'] for n in G.nodes()]

    # Size nodes by degree
    node_size = [300 * (0.1 + G.nodes[n].get('degree', 1)) for n in G.nodes()]

    nx.draw_networkx_nodes(G, pos, node_size=node_size, node_color=node_colors, cmap=plt.cm.tab20)
    nx.draw_networkx_edges(G, pos, alpha=0.3, edge_color='gray')

    # Label key nodes
    key_nodes = {n: n for n in G.nodes() if G.nodes[n].get('degree', 0) > 10}
    nx.draw_networkx_labels(G, pos, labels=key_nodes, font_size=10)

    plt.title("Network Communities")
    plt.axis('off')
    plt.savefig(output_path, dpi=300)
    plt.close()

    # Save community information
    community_data = []
    for i, community in enumerate(communities):
        for gene in community:
            community_data.append({
                'Gene': gene,
                'Community': i,
                'CommunitySize': len(community)
            })

    pd.DataFrame(community_data).to_csv(output_path.replace('.png', '.csv'), index=False)

    return communities


# 8. Perform pathway enrichment analysis
def pathway_enrichment(genes, output_path):
    """Perform pathway enrichment analysis using g:Profiler"""
    try:
        from gprofiler import GProfiler

        gp = GProfiler(return_dataframe=True)
        result = gp.profile(organism='hsapiens', query=genes)

        # Filter for relevant sources
        filtered = result[result['source'].isin(['GO:BP', 'KEGG', 'REAC'])]

        # Sort by p-value
        filtered = filtered.sort_values('p_value')

        # Save results
        filtered.to_csv(output_path, index=False)

        print(f"Top 5 enriched pathways:")
        print(filtered[['name', 'p_value', 'source']].head())

        return filtered
    except ImportError:
        print("gprofiler-python package not found. Install with: pip install gprofiler-python")
        print("Continuing without pathway enrichment analysis...")
        return None


# 9. Analyze pathway-specific subnetwork
def analyze_pathway_subnetwork(G, pathway_genes, pathway_name, output_path):
    """Extract and analyze subnetwork of pathway genes"""
    # Get pathway genes that exist in the network
    pathway_in_network = [gene for gene in pathway_genes if gene in G.nodes()]

    # Create subgraph of pathway genes plus their direct neighbors
    subgraph_nodes = set(pathway_in_network)
    for gene in pathway_in_network:
        subgraph_nodes.update(G.neighbors(gene))

    subgraph = G.subgraph(subgraph_nodes)

    print(f"{pathway_name} subnetwork: {len(pathway_in_network)}/{len(pathway_genes)} pathway genes")
    print(f"Total subnetwork size: {subgraph.number_of_nodes()} nodes, {subgraph.number_of_edges()} edges")

    # Visualize the subnetwork
    plt.figure(figsize=(12, 12))
    pos = nx.spring_layout(subgraph, seed=42)

    # Color nodes (red for pathway genes, blue for others)
    node_colors = ['red' if node in pathway_genes else 'blue' for node in subgraph.nodes()]

    # Size nodes by degree
    node_size = [300 * (0.1 + subgraph.degree(n)) for n in subgraph.nodes()]

    nx.draw_networkx_nodes(subgraph, pos, node_size=node_size, node_color=node_colors)
    nx.draw_networkx_edges(subgraph, pos, alpha=0.4, edge_color='gray')
    nx.draw_networkx_labels(subgraph, pos, font_size=8)

    plt.title(f"{pathway_name} Pathway Subnetwork")
    plt.axis('off')
    plt.savefig(output_path, dpi=300)
    plt.close()

    # Save subnetwork metrics
    subgraph_metrics = pd.DataFrame([
        {
            'Gene': node,
            'InPathway': node in pathway_genes,
            'Degree': subgraph.degree(node),
            'log2FC': G.nodes[node].get('log2FC', 0),
            'padj': G.nodes[node].get('padj', 1)
        }
        for node in subgraph.nodes()
    ])

    subgraph_metrics.to_csv(output_path.replace('.png', '_metrics.csv'), index=False)

    return subgraph


# 10. Analyze hub gene expression
def analyze_hub_expression(G, output_path):
    """Analyze relationship between hub status and expression changes"""
    # Extract data
    nodes_data = [(n,
                   G.nodes[n].get('degree', 0),
                   G.nodes[n].get('log2FC', 0),
                   G.nodes[n].get('padj', 1))
                  for n in G.nodes()]

    df = pd.DataFrame(nodes_data, columns=['Gene', 'Degree', 'log2FC', 'padj'])

    # Categorize by significance
    df['Significant'] = (df['padj'] < 0.05) & (abs(df['log2FC']) > 1)

    # Create scatter plot
    plt.figure(figsize=(10, 8))
    plt.scatter(df[df['Significant']]['Degree'], df[df['Significant']]['log2FC'],
                color='red', alpha=0.7, label='Significant DEGs')
    plt.scatter(df[~df['Significant']]['Degree'], df[~df['Significant']]['log2FC'],
                color='blue', alpha=0.3, label='Non-significant')

    # Label top genes
    for i, row in df.nlargest(10, 'Degree').iterrows():
        plt.annotate(row['Gene'], (row['Degree'], row['log2FC']),
                     fontsize=9, alpha=0.8)

    plt.xlabel('Degree (Hub Score)')
    plt.ylabel('log2 Fold Change')
    plt.axhline(y=0, color='gray', linestyle='--', alpha=0.3)
    plt.legend()
    plt.title('Hub Score vs. Expression Change')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

    # Calculate correlation
    correlation = df['Degree'].corr(df['log2FC'].abs())
    print(f"Correlation between hub score and expression change magnitude: {correlation:.3f}")

    return df


# 11. Compare different centrality measures
def visualize_centrality_comparison(G, output_path, top_n=15):
    """Compare different centrality measures for top genes"""
    try:
        centrality_metrics = {
            'Degree': dict(G.degree()),
            'Betweenness': nx.betweenness_centrality(G),
            'Closeness': nx.closeness_centrality(G)
        }

        # Eigenvector centrality might not converge for some networks
        try:
            centrality_metrics['Eigenvector'] = nx.eigenvector_centrality(G, max_iter=300)
        except:
            print("Eigenvector centrality calculation failed, using only 3 centrality measures")

        # Combine into dataframe
        df = pd.DataFrame({metric: values for metric, values in centrality_metrics.items()})
        df.index.name = 'Gene'
        df = df.reset_index()

        # Select top genes by average rank across metrics
        df['avg_rank'] = df[list(centrality_metrics.keys())].rank(ascending=False).mean(axis=1)
        top_genes = df.nsmallest(top_n, 'avg_rank')

        # Plot
        plt.figure(figsize=(12, 8))
        x = np.arange(len(top_genes))
        width = 0.2
        metrics = list(centrality_metrics.keys())

        # Normalize values for comparison
        for metric in metrics:
            top_genes[f'norm_{metric}'] = top_genes[metric] / max(top_genes[metric])

        # Plot bars for each metric
        for i, metric in enumerate(metrics):
            plt.bar(x + width * (i - len(metrics) / 2 + 0.5),
                    top_genes[f'norm_{metric}'],
                    width,
                    label=metric)

        plt.xlabel('Genes')
        plt.ylabel('Normalized Centrality')
        plt.title('Comparison of Centrality Measures')
        plt.xticks(x, top_genes['Gene'], rotation=45, ha='right')
        plt.legend()
        plt.tight_layout()
        plt.savefig(output_path, dpi=300)
        plt.close()

        # Save detailed metrics
        top_genes[['Gene'] + metrics].to_csv(output_path.replace('.png', '.csv'), index=False)

        return top_genes
    except Exception as e:
        print(f"Error in centrality comparison: {e}")
        return None


# 12. Visualize pathway enrichment results
def visualize_pathway_enrichment(enrichment_results, output_path, top_n=15):
    """Create bar chart of top enriched pathways"""
    if enrichment_results is None or enrichment_results.empty:
        print("No pathway enrichment results to visualize")
        return None

    try:
        # Filter and get top pathways
        filtered = enrichment_results[enrichment_results['source'].isin(['GO:BP', 'KEGG', 'REAC'])]
        top_pathways = filtered.nsmallest(top_n, 'p_value')

        # Create bar chart
        plt.figure(figsize=(12, 8))
        bars = plt.barh(range(len(top_pathways)), -np.log10(top_pathways['p_value']))

        # Color by source
        sources = top_pathways['source'].unique()
        colors = plt.cm.tab10(range(len(sources)))
        source_color = {source: color for source, color in zip(sources, colors)}

        for i, bar in enumerate(bars):
            bar.set_color(source_color[top_pathways.iloc[i]['source']])

        # Create shorter pathway names for display
        display_names = []
        for name in top_pathways['name']:
            if len(name) > 50:
                display_names.append(name[:47] + '...')
            else:
                display_names.append(name)

        plt.yticks(range(len(top_pathways)), display_names, fontsize=9)
        plt.xlabel('-log10(p-value)')
        plt.title('Top Enriched Pathways')

        # Add legend
        handles = [plt.Rectangle((0, 0), 1, 1, color=color) for color in source_color.values()]
        plt.legend(handles, source_color.keys(), title='Source')

        plt.tight_layout()
        plt.savefig(output_path, dpi=300)
        plt.close()

        return top_pathways
    except Exception as e:
        print(f"Error visualizing pathway enrichment: {e}")
        return None


# 13. Visualize community expression patterns
def visualize_community_heatmap(G, communities, deg_data, output_path):
    """Create heatmap of gene expression by community"""
    try:
        # Extract community data
        community_gene_data = []

        for i, community in enumerate(communities):
            for gene in community:
                if gene in deg_data['gene'].values:
                    idx = deg_data[deg_data['gene'] == gene].index[0]
                    community_gene_data.append({
                        'Gene': gene,
                        'Community': f"Community {i}",
                        'log2FC': deg_data.loc[idx, 'log2FoldChange'],
                        'padj': deg_data.loc[idx, 'padj'],
                        'Significant': deg_data.loc[idx, 'padj'] < 0.05
                    })

        df = pd.DataFrame(community_gene_data)
        if df.empty:
            print("No overlapping genes between communities and DEG data")
            return None

        # Filter to reasonable size for visualization
        if len(df) > 100:
            # Keep significant genes and a few from each community
            sig_genes = df[df['Significant']].copy()
            other_genes = df[~df['Significant']].groupby('Community').head(3)
            df = pd.concat([sig_genes, other_genes])

        # Sort communities by average log2FC
        community_order = df.groupby('Community')['log2FC'].mean().sort_values().index.tolist()

        # Filter to top genes if still too many
        if len(df) > 50:
            top_genes = []
            for comm in community_order:
                comm_genes = pd.concat([
                    df[df['Community'] == comm].nlargest(5, 'log2FC'),
                    df[df['Community'] == comm].nsmallest(5, 'log2FC')
                ])
                top_genes.append(comm_genes)
            df = pd.concat(top_genes)

        # Create pivoted data for heatmap
        pivot_data = df.pivot_table(index='Gene', columns='Community', values='log2FC', aggfunc='mean')

        # Reorder columns based on community order
        pivot_data = pivot_data[community_order]

        # Plot heatmap
        plt.figure(figsize=(10, max(8, len(pivot_data) / 4)))
        sns.heatmap(pivot_data, cmap='coolwarm', center=0,
                    linewidths=0.5, linecolor='lightgray',
                    cbar_kws={'label': 'log2 Fold Change'})
        plt.title('Gene Expression by Network Community')
        plt.tight_layout()
        plt.savefig(output_path, dpi=300)
        plt.close()

        # Save data
        df.to_csv(output_path.replace('.png', '.csv'), index=False)

        return df
    except Exception as e:
        print(f"Error in community heatmap visualization: {e}")
        return None


# 14. Create volcano plot with hub genes highlighted
def visualize_volcano_with_hubs(deg_data, hub_genes, output_path):
    """Create volcano plot highlighting hub genes"""
    try:
        plt.figure(figsize=(10, 8))

        # Plot all genes
        plt.scatter(deg_data['log2FoldChange'], -np.log10(deg_data['padj']),
                    alpha=0.3, color='gray', s=20)

        # Highlight hub genes
        hub_data = deg_data[deg_data['gene'].isin(hub_genes)]
        plt.scatter(hub_data['log2FoldChange'], -np.log10(hub_data['padj']),
                    alpha=1.0, color='red', s=50)

        # Add labels for top hub genes
        for i, row in hub_data.iterrows():
            if row['gene'] in hub_genes[:10]:  # Top 10 hubs
                plt.annotate(row['gene'],
                             (row['log2FoldChange'], -np.log10(row['padj'])),
                             fontsize=9, alpha=0.8)

        # Add significance lines
        plt.axhline(-np.log10(0.05), color='blue', linestyle='--', alpha=0.3)
        plt.axvline(-1, color='blue', linestyle='--', alpha=0.3)
        plt.axvline(1, color='blue', linestyle='--', alpha=0.3)

        plt.xlabel('log2 Fold Change')
        plt.ylabel('-log10(p-value)')
        plt.title('Volcano Plot with Hub Genes Highlighted')
        plt.tight_layout()
        plt.savefig(output_path, dpi=300)
        plt.close()

        return hub_data
    except Exception as e:
        print(f"Error in volcano plot visualization: {e}")
        return None


# 15. Visualize network motifs
def visualize_network_motifs(G, pathway_genes, pathway_name, output_path):
    """Identify and visualize common network motifs"""
    try:
        # Find triangles (3-node cliques)
        triangles = list(nx.triangles(G).items())
        triangles = sorted(triangles, key=lambda x: x[1], reverse=True)

        # Find nodes participating in most triangles
        top_nodes = [node for node, count in triangles[:15] if count > 0]

        if not top_nodes:
            print("No significant triangular motifs found in network")
            return None

        # Extract subgraph containing these nodes
        motif_graph = G.subgraph(top_nodes)

        # Visualize
        plt.figure(figsize=(10, 10))
        pos = nx.spring_layout(motif_graph, seed=42)

        # Color nodes (red for pathway genes)
        node_colors = ['red' if node in pathway_genes else 'skyblue' for node in motif_graph.nodes()]

        nx.draw_networkx_nodes(motif_graph, pos, node_color=node_colors, node_size=500)
        nx.draw_networkx_edges(motif_graph, pos, width=2, alpha=0.7)
        nx.draw_networkx_labels(motif_graph, pos, font_size=12)

        plt.title(f'Key Network Motifs ({pathway_name})')
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(output_path, dpi=300)
        plt.close()

        # Save motif nodes info
        motif_data = pd.DataFrame([
            {
                'Gene': node,
                'TriangleCount': dict(triangles)[node],
                'IsPathwayGene': node in pathway_genes,
                'Degree': G.degree(node),
                'log2FC': G.nodes[node].get('log2FC', 0)
            }
            for node in motif_graph.nodes()
        ])

        motif_data.to_csv(output_path.replace('.png', '.csv'), index=False)

        return motif_data
    except Exception as e:
        print(f"Error in network motif visualization: {e}")
        return None


# 16. Shortest Path Analysis
def analyze_shortest_paths(G, pathway_genes, hub_genes, output_path, max_paths=10):
    """Analyze shortest paths between pathway genes and hub genes"""
    try:
        print("Analyzing shortest paths between pathway and hub genes...")

        # Filter for genes actually in the network
        pathway_in_network = [gene for gene in pathway_genes if gene in G.nodes()]
        hub_in_network = [gene for gene in hub_genes if gene in G.nodes()]

        if not pathway_in_network or not hub_in_network:
            print("No pathway genes or hub genes found in network")
            return None

        # Find all shortest paths between pathway and hub genes
        all_paths = []

        for path_gene in pathway_in_network:
            for hub_gene in hub_in_network:
                if path_gene != hub_gene:
                    try:
                        paths = list(nx.all_shortest_paths(G, path_gene, hub_gene))
                        for path in paths:
                            all_paths.append({
                                'source': path_gene,
                                'target': hub_gene,
                                'path': path,
                                'length': len(path) - 1  # Subtract 1 to get edge count
                            })
                    except nx.NetworkXNoPath:
                        continue

        if not all_paths:
            print("No paths found between pathway and hub genes")
            return None

        # Convert to DataFrame
        paths_df = pd.DataFrame(all_paths)

        # Calculate statistics
        avg_path_length = paths_df['length'].mean()
        print(f"Average shortest path length: {avg_path_length:.2f} edges")

        # Get most common intermediary genes
        intermediaries = []
        for path in paths_df['path']:
            if len(path) > 2:  # If there are intermediaries
                intermediaries.extend(path[1:-1])

        intermediary_counts = pd.Series(intermediaries).value_counts()

        # Visualize most frequent paths
        sorted_paths = paths_df.sort_values('length')
        top_paths = sorted_paths.head(max_paths)

        plt.figure(figsize=(12, 10))

        for i, (_, row) in enumerate(top_paths.iterrows()):
            # Create subgraph for this path
            path_graph = G.subgraph(row['path'])
            pos = nx.spring_layout(path_graph, seed=42)

            plt.subplot(min(5, max_paths), max(1, math.ceil(max_paths / 5)), i + 1)

            # Color nodes (red for pathway, blue for hub, green for intermediaries)
            node_colors = []
            for node in path_graph.nodes():
                if node == row['source']:
                    node_colors.append('red')
                elif node == row['target']:
                    node_colors.append('blue')
                else:
                    node_colors.append('green')

            nx.draw_networkx_nodes(path_graph, pos, node_color=node_colors, node_size=300)
            nx.draw_networkx_edges(path_graph, pos, width=2)
            nx.draw_networkx_labels(path_graph, pos, font_size=8)

            plt.title(f"{row['source']} → {row['target']} ({row['length']} steps)")
            plt.axis('off')

        plt.tight_layout()
        plt.savefig(output_path, dpi=300)
        plt.close()

        # Save path data
        paths_df.to_csv(output_path.replace('.png', '.csv'), index=False)

        # Save intermediary statistics
        if not intermediary_counts.empty:
            pd.DataFrame({'Gene': intermediary_counts.index, 'PathCount': intermediary_counts.values}).to_csv(
                output_path.replace('.png', '_intermediaries.csv'), index=False)

        return paths_df
    except Exception as e:
        print(f"Error in shortest path analysis: {e}")
        return None

# 17. Compare to random networks
def compare_to_random_networks(G, n_random=100):
    """Compare network properties to random networks"""
    try:
        print("Comparing network properties to random models...")

        # Network properties to compare
        props = {
            'Clustering': nx.average_clustering(G),
            'Assortativity': nx.degree_assortativity_coefficient(G),
            'Number of Communities': len(list(nx.algorithms.community.greedy_modularity_communities(G)))
        }

        # Add avg path length only if connected
        if nx.is_connected(G):
            props['Avg Shortest Path'] = nx.average_shortest_path_length(G)

        # Generate random networks with same degree sequence
        random_props = {p: [] for p in props}
        degree_sequence = [d for _, d in G.degree()]

        for i in range(n_random):
            try:
                R = nx.configuration_model(degree_sequence, seed=i)
                R = nx.Graph(R)  # Remove parallel edges
                R.remove_edges_from(nx.selfloop_edges(R))  # Remove self-loops

                random_props['Clustering'].append(nx.average_clustering(R))
                random_props['Assortativity'].append(nx.degree_assortativity_coefficient(R))
                random_props['Number of Communities'].append(
                    len(list(nx.algorithms.community.greedy_modularity_communities(R))))

                if 'Avg Shortest Path' in props and nx.is_connected(R):
                    random_props['Avg Shortest Path'].append(nx.average_shortest_path_length(R))
            except Exception as e:
                pass

        # Calculate statistics
        results = []
        for prop in props:
            random_values = random_props[prop]
            if not random_values:
                continue

            random_mean = np.mean(random_values)
            random_std = np.std(random_values)
            z_score = (props[prop] - random_mean) / random_std if random_std > 0 else 0

            # One-sample t-test
            if len(random_values) > 1:
                _, p_value = stats.ttest_1samp(random_values, props[prop])
            else:
                p_value = 1

            results.append({
                'Property': prop,
                'Observed': props[prop],
                'Random_Mean': random_mean,
                'Random_Std': random_std,
                'Z_score': z_score,
                'P_value': p_value
            })

        results_df = pd.DataFrame(results)
        print("Network comparison results:")
        print(results_df)

        return results_df
    except Exception as e:
        print(f"Error in random network comparison: {e}")
        return None

# Main function to run the analysis
def main_pathway_analysis():
    # Set up paths
    base_dir = "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158"
    data_dir = os.path.join(base_dir, "DESeq2 outputs")
    output_dir = os.path.join(base_dir, "PPI outputs")

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Load DESeq2 results
    deg_results = load_deseq_results(
        os.path.join(data_dir, "deseq2_results_F_OUD_vs_F_None.csv"))

    # Filter for significant genes (adjust thresholds as needed)
    sig_genes = deg_results[(deg_results['padj'] < 0.05) & (abs(deg_results['log2FoldChange']) > 1)]
    print(f"Loaded {len(sig_genes)} significant DEGs")

    # Load pathway genes - rename this variable to be pathway-agnostic
    pathway_name = "Complement"  # Change this for different pathways
    log_path = setup_logging(output_dir, pathway_name)
    output_path, orig_stdout, orig_stderr, output_file = setup_output_capture(output_dir, pathway_name)
    pathway_genes = load_pathway_genes(
        "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/KEGG outputs/kegg_complement_unique_genes.csv")
    try:
        print(f"Loaded {len(pathway_genes)} {pathway_name} genes from file")

        # Combine DEGs and pathway genes for PPI analysis
        genes_to_analyze = list(set(sig_genes['gene'].tolist() + pathway_genes))

        # If too many genes, prioritize most significant DEGs and pathway genes
        if len(genes_to_analyze) > 400:
            # Equal prioritization - balance between DEGs and pathway genes
            max_pathway_genes = min(100, len(pathway_genes))
            max_deg_genes = 400 - max_pathway_genes

            sig_genes_subset = sig_genes.sort_values('padj').head(max_deg_genes)
            genes_to_analyze = list(set(sig_genes_subset['gene'].tolist() + pathway_genes[:max_pathway_genes]))

        print(f"Analyzing PPI network for {len(genes_to_analyze)} genes")

        # Get protein interactions from STRING
        try:
            ppi_data = get_string_interactions(genes_to_analyze)
            print(f"Retrieved {len(ppi_data)} interactions from STRING")

            # Save raw PPI data
            pd.DataFrame(ppi_data).to_csv(os.path.join(output_dir, "string_ppi_data.csv"), index=False)
        except Exception as e:
            print(f"Error fetching STRING data: {e}")
            return

        # Create and analyze network
        G = create_ppi_network(ppi_data, deg_data=deg_results)
        print(f"Created network with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")

        # Basic network visualization
        visualize_network(G, os.path.join(output_dir, f"{pathway_name}_ppi_network.png"))

        # Identify hub genes
        top_hubs = identify_hubs(G, os.path.join(output_dir, f"{pathway_name}_network_metrics.csv"))
        print(f"Top hub genes: {', '.join(top_hubs['Gene'].head(5).tolist())}")

        # Detect communities
        communities = detect_communities(G, os.path.join(output_dir, f"{pathway_name}_network_communities.png"))

        # Pathway enrichment
        enrichment_results = pathway_enrichment(list(G.nodes()),
                                                os.path.join(output_dir, f"{pathway_name}_pathway_enrichment.csv"))

        # Analyze pathway subnetwork
        pathway_subgraph = analyze_pathway_subnetwork(G, pathway_genes, pathway_name,
                                                      os.path.join(output_dir, f"{pathway_name}_subnetwork.png"))

        # Hub gene expression analysis
        hub_expression = analyze_hub_expression(G, os.path.join(output_dir, f"{pathway_name}_hub_expression.png"))

        # Compare centrality measures
        visualize_centrality_comparison(G, os.path.join(output_dir, f"{pathway_name}_centrality_comparison.png"))

        # Visualize pathway enrichment
        if enrichment_results is not None:
            visualize_pathway_enrichment(enrichment_results,
                                         os.path.join(output_dir, f"{pathway_name}_top_pathways.png"))

        # Community expression heatmap
        visualize_community_heatmap(G, communities, deg_results,
                                    os.path.join(output_dir, f"{pathway_name}_community_heatmap.png"))

        # Volcano plot with hubs
        visualize_volcano_with_hubs(deg_results, top_hubs['Gene'].tolist(),
                                    os.path.join(output_dir, f"{pathway_name}_volcano_with_hubs.png"))

        # Network motifs
        visualize_network_motifs(G, pathway_genes, pathway_name,
                                 os.path.join(output_dir, f"{pathway_name}_network_motifs.png"))

        # Shortest path analysis
        shortest_paths = analyze_shortest_paths(G, pathway_genes, top_hubs['Gene'].tolist()[:20],
                                                os.path.join(output_dir, f"{pathway_name}_shortest_paths.png"))

        # Compare to random networks
        random_comparison = compare_to_random_networks(G)
        if random_comparison is not None:
            random_comparison.to_csv(os.path.join(output_dir, f"{pathway_name}_random_network_comparison.csv"),
                                     index=False)

        print(f"{pathway_name} network analysis complete!")

        print(f"{pathway_name} network analysis complete!")
        print(f"Log saved to: {log_path}")

    finally:
        # Always restore the original stdout/stderr
        end_output_capture(orig_stdout, orig_stderr, output_file)

def main_cli():
    parser = argparse.ArgumentParser(description="Run STRING PPI analysis from gene list.")
    parser.add_argument("--genes", required=True, help="Path to gene list file (one gene per line)")
    parser.add_argument("--species", default="9606", help="NCBI species ID (default: 9606 for human)")
    parser.add_argument("--score", type=int, default=400, help="Minimum interaction score (default: 400)")
    parser.add_argument("--outdir", default="ppi_output", help="Output directory")
    parser.add_argument("--plot", action="store_true", help="Plot interaction score distribution")
    parser.add_argument("--graphml", action="store_true", help="Export interaction network as GraphML")
    parser.add_argument("--networkplot", action="store_true", help="Plot the interaction network using NetworkX")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Set up logging
    log_path = setup_logging(args.outdir)

    with open(args.genes) as f:
        genes = [line.strip() for line in f if line.strip()]
    print(f"Loaded {len(genes)} genes from {args.genes}")

    string_ids_df = get_string_ids(genes, args.species)
    if string_ids_df is None or string_ids_df.empty:
        return

    mapped_ids = string_ids_df['stringId'].tolist()
    print(f"Mapped {len(mapped_ids)} genes to STRING IDs")

    interactions_df = get_interactions(mapped_ids, args.species, args.score)
    if interactions_df is None or interactions_df.empty:
        return

    output_file = os.path.join(args.outdir, "ppi_interactions.tsv")
    interactions_df.to_csv(output_file, sep='\t', index=False)
    print(f"Interactions saved to {output_file}")

    if args.plot:
        plot_path = os.path.join(args.outdir, "score_distribution.png")
        plot_score_distribution(interactions_df, plot_path)

    if args.graphml:
        graphml_path = os.path.join(args.outdir, "ppi_network.graphml")
        export_graphml(interactions_df, graphml_path)

    if args.networkplot:
        network_plot_path = os.path.join(args.outdir, "ppi_network_plot.png")
        plot_interaction_network(interactions_df, network_plot_path)



# Execute main function if script is run directly
if __name__ == "__main__":
    import sys
    if len(sys.argv) == 1:
        main_pathway_analysis()  # Run the pathway analysis
    else:
        main_cli()  # Run the CLI version