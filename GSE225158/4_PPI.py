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

warnings.filterwarnings("ignore", category=UserWarning)


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


# 9. Analyze complement-specific subnetwork
def analyze_complement_subnetwork(G, complement_genes, output_path):
    """Extract and analyze subnetwork of complement genes"""
    # Get complement genes that exist in the network
    complement_in_network = [gene for gene in complement_genes if gene in G.nodes()]

    # Create subgraph of complement genes plus their direct neighbors
    subgraph_nodes = set(complement_in_network)
    for gene in complement_in_network:
        subgraph_nodes.update(G.neighbors(gene))

    subgraph = G.subgraph(subgraph_nodes)

    print(f"Complement subnetwork: {len(complement_in_network)}/{len(complement_genes)} complement genes")
    print(f"Total subnetwork size: {subgraph.number_of_nodes()} nodes, {subgraph.number_of_edges()} edges")

    # Visualize the subnetwork
    plt.figure(figsize=(12, 12))
    pos = nx.spring_layout(subgraph, seed=42)

    # Color nodes (red for complement genes, blue for others)
    node_colors = ['red' if node in complement_genes else 'blue' for node in subgraph.nodes()]

    # Size nodes by degree
    node_size = [300 * (0.1 + subgraph.degree(n)) for n in subgraph.nodes()]

    nx.draw_networkx_nodes(subgraph, pos, node_size=node_size, node_color=node_colors)
    nx.draw_networkx_edges(subgraph, pos, alpha=0.4, edge_color='gray')
    nx.draw_networkx_labels(subgraph, pos, font_size=8)

    plt.title("Complement System Subnetwork")
    plt.axis('off')
    plt.savefig(output_path, dpi=300)
    plt.close()

    # Save subnetwork metrics
    subgraph_metrics = pd.DataFrame([
        {
            'Gene': node,
            'InComplement': node in complement_genes,
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
                comm_genes = df[df['Community'] == comm].nlargest(5, 'log2FC').append(
                    df[df['Community'] == comm].nsmallest(5, 'log2FC'))
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
def visualize_network_motifs(G, complement_genes, output_path):
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

        # Color nodes (red for complement genes)
        node_colors = ['red' if node in complement_genes else 'skyblue' for node in motif_graph.nodes()]

        nx.draw_networkx_nodes(motif_graph, pos, node_color=node_colors, node_size=500)
        nx.draw_networkx_edges(motif_graph, pos, width=2, alpha=0.7)
        nx.draw_networkx_labels(motif_graph, pos, font_size=12)

        plt.title('Key Network Motifs')
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(output_path, dpi=300)
        plt.close()

        # Save motif nodes info
        motif_data = pd.DataFrame([
            {
                'Gene': node,
                'TriangleCount': dict(triangles)[node],
                'IsComplement': node in complement_genes,
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


# 18. Shortest Path Analysis
def analyze_shortest_paths(G, complement_genes, hub_genes, output_path, max_paths=10):
    """Analyze shortest paths between complement genes and hub genes"""
    try:
        print("Analyzing shortest paths between complement and hub genes...")

        # Filter for genes actually in the network
        complement_in_network = [gene for gene in complement_genes if gene in G.nodes()]
        hub_in_network = [gene for gene in hub_genes if gene in G.nodes()]

        if not complement_in_network or not hub_in_network:
            print("No complement genes or hub genes found in network")
            return None

        # Find all shortest paths between complement and hub genes
        all_paths = []

        for comp_gene in complement_in_network:
            for hub_gene in hub_in_network:
                if comp_gene != hub_gene:
                    try:
                        paths = list(nx.all_shortest_paths(G, comp_gene, hub_gene))
                        for path in paths:
                            all_paths.append({
                                'source': comp_gene,
                                'target': hub_gene,
                                'path': path,
                                'length': len(path) - 1  # Subtract 1 to get edge count
                            })
                    except nx.NetworkXNoPath:
                        continue

        if not all_paths:
            print("No paths found between complement and hub genes")
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

            # Color nodes (red for complement, blue for hub, green for intermediaries)
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
        intermediary_counts.head(20).to_csv(output_path.replace('.png', '_intermediaries.csv'))

        return {
            'paths': paths_df,
            'avg_length': avg_path_length,
            'intermediaries': intermediary_counts
        }
    except Exception as e:
        print(f"Error in shortest path analysis: {e}")
        return None


# 19. Random Network Comparison
def compare_to_random_networks(G, output_path, n_random=100):
    """Compare network properties with random networks"""
    try:
        print("Comparing network to random networks...")

        # Original network metrics
        original_metrics = {
            'Clustering': nx.average_clustering(G),
            'Assortativity': nx.degree_assortativity_coefficient(G),
            'Avg Shortest Path': nx.average_shortest_path_length(G),
            'Diameter': nx.diameter(G)
        }

        print(f"Original network metrics calculated")

        # Generate random networks with same number of nodes and edges
        n_nodes = G.number_of_nodes()
        n_edges = G.number_of_edges()

        random_metrics = []
        for i in range(n_random):
            if i % 10 == 0:
                print(f"Generating random network {i}/{n_random}...")

            # Create random graph with same node/edge count
            random_G = nx.gnm_random_graph(n_nodes, n_edges, seed=i)

            # Calculate metrics
            metrics = {}
            try:
                metrics['Clustering'] = nx.average_clustering(random_G)
                metrics['Assortativity'] = nx.degree_assortativity_coefficient(random_G)

                # Only calculate these if graph is connected
                if nx.is_connected(random_G):
                    metrics['Avg Shortest Path'] = nx.average_shortest_path_length(random_G)
                    metrics['Diameter'] = nx.diameter(random_G)
                else:
                    largest_cc = max(nx.connected_components(random_G), key=len)
                    largest_cc_graph = random_G.subgraph(largest_cc)
                    metrics['Avg Shortest Path'] = nx.average_shortest_path_length(largest_cc_graph)
                    metrics['Diameter'] = nx.diameter(largest_cc_graph)
            except:
                continue

            random_metrics.append(metrics)

        # Convert to DataFrame
        random_df = pd.DataFrame(random_metrics)

        # Calculate z-scores
        z_scores = {}
        p_values = {}
        for metric in original_metrics:
            mean = random_df[metric].mean()
            std = random_df[metric].std()
            z_scores[metric] = (original_metrics[metric] - mean) / std

            # Calculate p-value (two-tailed)
            if original_metrics[metric] > mean:
                p_values[metric] = 2 * (1 - stats.norm.cdf(z_scores[metric]))
            else:
                p_values[metric] = 2 * stats.norm.cdf(z_scores[metric])

        # Visualize comparison
        plt.figure(figsize=(14, 8))

        # Create subplots for each metric
        for i, metric in enumerate(original_metrics):
            plt.subplot(2, 2, i + 1)

            # Plot histogram of random networks
            sns.histplot(random_df[metric], kde=True)

            # Add vertical line for original network
            plt.axvline(original_metrics[metric], color='red', linestyle='--')

            plt.title(f"{metric}\nZ-score: {z_scores[metric]:.2f}, p-value: {p_values[metric]:.3e}")
            plt.xlabel(metric)
            plt.ylabel('Frequency')

        plt.tight_layout()
        plt.savefig(output_path, dpi=300)
        plt.close()

        # Save comparison data
        comparison_df = pd.DataFrame({
            'Metric': list(original_metrics.keys()),
            'Original': list(original_metrics.values()),
            'Random_Mean': [random_df[m].mean() for m in original_metrics],
            'Random_Std': [random_df[m].std() for m in original_metrics],
            'Z_Score': list(z_scores.values()),
            'P_Value': list(p_values.values())
        })

        comparison_df.to_csv(output_path.replace('.png', '.csv'), index=False)

        print("Network comparison completed")
        return {
            'original': original_metrics,
            'random': random_df,
            'z_scores': z_scores,
            'p_values': p_values
        }
    except Exception as e:
        print(f"Error in random network comparison: {e}")
        return None

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

    # 9. Detect network communities
    print("\nPerforming community detection analysis...")
    communities = detect_communities(G, f"{output_dir}/network_communities.png")

    # 10. Analyze complement-specific subnetwork
    print("\nAnalyzing complement-specific subnetwork...")
    complement_subgraph = analyze_complement_subnetwork(G, complement_genes, f"{output_dir}/complement_subnetwork.png")

    # 11. Analyze hub gene expression patterns
    print("\nAnalyzing hub gene expression patterns...")
    hub_expression = analyze_hub_expression(G, f"{output_dir}/hub_expression.png")

    # 12. Compare different centrality measures
    print("\nComparing different centrality measures...")
    centrality_comparison = visualize_centrality_comparison(G, f"{output_dir}/centrality_comparison.png")

    # 13. Perform pathway enrichment on hub genes (top 50)
    print("\nPerforming pathway enrichment analysis on hub genes...")
    try:
        top_hub_genes = top_hubs['Gene'].tolist()
        pathway_results = pathway_enrichment(top_hub_genes, f"{output_dir}/pathway_enrichment.csv")

# 14. Visualize pathway enrichment results
        if pathway_results is not None:
            print("\nVisualizing pathway enrichment results...")
            visualize_pathway_enrichment(pathway_results, f"{output_dir}/pathway_enrichment_plot.png")
    except Exception as e:
        print(f"Error in pathway enrichment analysis: {e}")
        print("Continuing without pathway enrichment...")

    # 15. Create community expression heatmap
    print("\nCreating community expression heatmap...")
    community_expression = visualize_community_heatmap(G, communities, deg_results, f"{output_dir}/community_expression.png")

    # 16. Create volcano plot with hub genes
    print("\nCreating volcano plot highlighting hub genes...")
    volcano_hubs = visualize_volcano_with_hubs(deg_results, top_hubs['Gene'].tolist(), f"{output_dir}/volcano_with_hubs.png")

    # 17. Analyze network motifs
    print("\nAnalyzing network motifs...")
    network_motifs = visualize_network_motifs(G, complement_genes, f"{output_dir}/network_motifs.png")

    # 18. Analyze shortest paths between complement and hub genes
    print("\nAnalyzing shortest paths between complement and hub genes...")
    shortest_paths = analyze_shortest_paths(G, complement_genes, top_hubs['Gene'].tolist()[:20],
                                            f"{output_dir}/shortest_paths.png")

    # 19. Compare network to random networks
    print("\nComparing network properties to random networks...")
    random_comparison = compare_to_random_networks(G, f"{output_dir}/random_network_comparison.png", n_random=50)

    print("\nAll analyses completed!")
    print(f"Results saved to: {output_dir}")


if __name__ == "__main__":
    main()