import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import os
import requests
import time

# Create output directory if it doesn't exist
os.makedirs("/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/DL-PPI outputs/visualizations",
            exist_ok=True)

# Load data
male_df = pd.read_csv(
    "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/DL-PPI outputs/dl_ppi_viz_data_male_0.4.csv")
female_df = pd.read_csv(
    "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/DL-PPI outputs/dl_ppi_viz_data_female_0.4.csv")

# Print column names to debug
print("Male dataframe columns:", male_df.columns.tolist())

# Determine protein identifier columns
protein1_col = male_df.columns[0]  # First column
protein2_col = male_df.columns[1]  # Second column
print(f"Using columns {protein1_col} and {protein2_col} for protein identifiers")

# Get all unique accession numbers
all_accessions = set(pd.concat([
    male_df[protein1_col], male_df[protein2_col],
    female_df[protein1_col], female_df[protein2_col]
]).unique())


# Function to convert UniProt accession to gene name (one at a time to be more reliable)
def get_gene_name(uniprot_id):
    try:
        url = f"https://www.uniprot.org/uniprot/{uniprot_id}.xml"
        response = requests.get(url)
        if response.status_code == 200:
            # Extract gene name from XML - looking for name in first gene tag
            xml = response.text
            gene_start = xml.find("<gene")
            if gene_start > 0:
                name_start = xml.find("<name", gene_start)
                if name_start > 0:
                    name_content_start = xml.find(">", name_start) + 1
                    name_content_end = xml.find("</name>", name_content_start)
                    if name_content_start > 0 and name_content_end > 0:
                        return xml[name_content_start:name_content_end]
        return uniprot_id  # Return original ID if extraction fails
    except:
        return uniprot_id


# Create mapping from accession to gene name with progress tracking
print(f"Converting {len(all_accessions)} UniProt IDs to gene names...")
acc_to_gene = {}
for i, acc in enumerate(all_accessions):
    if i % 10 == 0:  # Progress update every 10 proteins
        print(f"Processed {i}/{len(all_accessions)} proteins")

    # Skip problematic or non-UniProt IDs
    if len(acc) < 5 or not acc[0].isalpha() or not acc[1].isdigit():
        acc_to_gene[acc] = acc
        continue

    gene_name = get_gene_name(acc)
    acc_to_gene[acc] = gene_name
    time.sleep(0.5)  # Avoid overloading the server


def create_network_visualization(df, title, output_filename, acc_to_gene):
    # Create graph with gene names
    G = nx.Graph()

    # Map of edges with weights
    edge_weights = {}

    # Add edges with weights based on prediction confidence
    print(f"Building network with {len(df)} interactions")
    for _, row in df.iterrows():
        acc1 = row[protein1_col]
        acc2 = row[protein2_col]

        # Get gene names, fallback to accession if not found
        gene1 = acc_to_gene.get(acc1, acc1)
        gene2 = acc_to_gene.get(acc2, acc2)

        # Use gene names as nodes
        G.add_edge(gene1, gene2, weight=row['prediction'])
        edge_weights[(gene1, gene2)] = row['prediction']

    # Print some gene names to verify conversion
    print(f"Sample of gene names used in network: {list(G.nodes())[:5]}")

    # Get statistics
    num_nodes = G.number_of_nodes()
    num_edges = G.number_of_edges()
    print(f"{title} Network Statistics:")
    print(f"  Nodes: {num_nodes}")
    print(f"  Edges: {num_edges}")

    # Get the top 10 nodes by degree (most connected)
    top_nodes = sorted(dict(G.degree()).items(), key=lambda x: x[1], reverse=True)[:10]
    print(f"  Top 10 connected proteins: {top_nodes}\n")

    # Create visualization
    plt.figure(figsize=(12, 10))

    # Node sizes based on degree
    node_sizes = [G.degree(node) * 50 + 100 for node in G.nodes()]

    # Edge widths based on prediction confidence
    edge_widths = [G[u][v]['weight'] * 2 for u, v in G.edges()]

    # Use spring layout for node positioning
    pos = nx.spring_layout(G, seed=42)

    # Draw the network
    nx.draw_networkx(
        G, pos,
        node_size=node_sizes,
        width=edge_widths,
        alpha=0.7,
        with_labels=False,
        node_color='skyblue',
        edge_color='gray'
    )

    # Draw labels for top connected nodes only
    top_node_names = [node for node, _ in top_nodes]
    labels = {node: node for node in G.nodes() if node in top_node_names}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=8)

    plt.title(f"{title} (Threshold ≥ 0.4)")
    plt.axis('off')
    plt.tight_layout()

    # Save figure
    plt.savefig(output_filename, dpi=300)
    print(f"Network visualization saved to {output_filename}")
    plt.close()


# Create visualizations
create_network_visualization(
    male_df,
    "Male OUD Complement Interactions",
    "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/DL-PPI outputs/visualizations/male_network_0.4.png",
    acc_to_gene
)

create_network_visualization(
    female_df,
    "Female OUD Complement Interactions",
    "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/DL-PPI outputs/visualizations/female_network_0.4.png",
    acc_to_gene
)

# Create comparison chart
plt.figure(figsize=(12, 5))

# Get sets of proteins using gene names
male_proteins = set([acc_to_gene.get(acc, acc) for acc in
                     pd.concat([male_df[protein1_col], male_df[protein2_col]]).unique()])
female_proteins = set([acc_to_gene.get(acc, acc) for acc in
                       pd.concat([female_df[protein1_col], female_df[protein2_col]]).unique()])

# Calculate values for comparison
only_male = len(male_proteins - female_proteins)
only_female = len(female_proteins - male_proteins)
common = len(male_proteins.intersection(female_proteins))

plt.subplot(1, 2, 1)
plt.title('Proteins in Networks')
venn_colors = ['lightblue', 'lightpink']
plt.bar(['Male-specific', 'Common', 'Female-specific'], [only_male, common, only_female],
        color=venn_colors + ['lightgreen'])
plt.ylabel('Number of proteins')

# Create subplot for edge count comparison
plt.subplot(1, 2, 2)
plt.title('Predicted Interactions')
plt.bar(['Male', 'Female'], [len(male_df), len(female_df)], color=venn_colors)
plt.ylabel('Number of interactions')

plt.tight_layout()
plt.savefig(
    "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/DL-PPI outputs/visualizations/comparison_0.4.png",
    dpi=300)
print("Comparison visualization saved")