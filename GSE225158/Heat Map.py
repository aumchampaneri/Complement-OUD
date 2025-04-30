import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

# Load the data
adata = sc.read_h5ad("/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE233279/GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5ad")

# Rename cell types
rename_dict = {
    "Microglia": "Microglia",
    "Oligos": "Oligodendrocytes",
    "Oligos_Pre": "OPCs",
    "Astrocytes": "Astrocytes",
    "Endothelial": "Endothelial cells",
    "Mural": "Mural cells",
    "D1-Matrix": "D1 MSNs (Matrix)",
    "D2-Matrix": "D2 MSNs (Matrix)",
    "D1-Striosome": "D1 MSNs (Striosome)",
    "D2-Striosome": "D2 MSNs (Striosome)",
    "D1/D2-Hybrid": "D1/D2 Hybrid MSNs",
    "Int-PTHLH": "PTHLH+ Interneurons",
    "Int-CCK": "CCK+ Interneurons",
    "Int-SST": "SST+ Interneurons",
    "Int-TH": "TH+ Interneurons"
}
adata.obs['celltype3_renamed'] = adata.obs['celltype3'].map(rename_dict)

# Load gene list
gene_list = pd.read_csv("/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158 Gene Library/gene_names.csv")
gene_names = gene_list['Gene Name'].tolist()
genes_present = [gene for gene in gene_names if gene in adata.var_names]

# Subset data
adata_subset = adata[:, genes_present].copy()

# Create group label
adata_subset.obs['group'] = (
    adata_subset.obs['celltype3_renamed'].astype(str) + "_" +
    adata_subset.obs['Sex'].astype(str) + "_" +
    adata_subset.obs['Dx_OUD'].astype(str)
)

# Average expression per group
mean_expr = adata_subset.to_df().groupby(adata_subset.obs['group']).mean().T

# Z-score normalize per gene
mean_expr_z = mean_expr.sub(mean_expr.mean(axis=1), axis=0).div(mean_expr.std(axis=1), axis=0)

# Extract group components
group_split = mean_expr.columns.str.split("_", expand=True)
group_split = pd.DataFrame(group_split.tolist(), index=mean_expr.columns, columns=['Cell Type', 'Sex', 'Diagnosis'])

# Define the desired order for cell types (alphabetically or as desired)
cell_type_order = [
    "Microglia", "Oligodendrocytes", "OPCs", "Astrocytes",
    "Endothelial cells", "Mural cells", "D1 MSNs (Matrix)",
    "D2 MSNs (Matrix)", "D1 MSNs (Striosome)", "D2 MSNs (Striosome)",
    "D1/D2 Hybrid MSNs", "PTHLH+ Interneurons", "CCK+ Interneurons",
    "SST+ Interneurons", "TH+ Interneurons"
]

# Reorder `group_split` based on the defined cell_type_order
group_split['Cell Type'] = pd.Categorical(group_split['Cell Type'], categories=cell_type_order, ordered=True)
sorted_columns = group_split.sort_values(['Cell Type', 'Diagnosis', 'Sex']).index
mean_expr_z = mean_expr_z[sorted_columns]
col_colors = pd.DataFrame({
    'Diagnosis': group_split['Diagnosis'].map({'None': 'indianred', 'OUD': 'skyblue'}),
    'Sex': group_split['Sex'].map({'F': 'gray', 'M': 'black'}),
    'Cell Type': group_split['Cell Type'].map(dict(zip(cell_type_order, sns.color_palette("tab20", len(cell_type_order)))))
}, index=mean_expr.columns)

# Plot the heatmap
g = sns.clustermap(
    mean_expr_z,
    col_cluster=False,  # Disable column clustering
    row_cluster=True,   # Keep row clustering
    col_colors=col_colors,
    cmap=sns.diverging_palette(145, 300, s=80, as_cmap=True),
    center=0,
    figsize=(14, 11),
    xticklabels=True,
    yticklabels=True
)

# Create legend patches
legend_patches = (
    [Patch(color=v, label=f'Cell Type: {k}') for k, v in dict(zip(cell_type_order, sns.color_palette("tab20", len(cell_type_order)))).items()] +
    [Patch(color=v, label=f'Diagnosis: {k}') for k, v in {'None': 'indianred', 'OUD': 'skyblue'}.items()] +
    [Patch(color=v, label=f'Sex: {k}') for k, v in {'F': 'gray', 'M': 'black'}.items()]
)

# Add legend to the heatmap
g.ax_heatmap.legend(
    handles=legend_patches,
    loc='upper left',
    bbox_to_anchor=(1.05, 1),
    borderaxespad=0,
    frameon=False
)

# Display plot
plt.show()

# Save the figure
g.savefig("heatmap_with_annotations.png", dpi=300, bbox_inches='tight')
