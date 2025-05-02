import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from sklearn.preprocessing import StandardScaler
'''
The heatmap shows how each group's average expression of a gene deviates from the overall mean 
'''
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

# List of cell types to plot (use the renamed labels)
celltypes_to_plot = [
    "Oligodendrocytes",
    "Astrocytes",
    "OPCs",
    "Microglia"
]

# Subset AnnData
adata_subset = adata[adata.obs['celltype3_renamed'].isin(celltypes_to_plot), :].copy()

# Load gene list
gene_list = pd.read_csv("/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158 Gene Library/gene_names.csv")
gene_names = gene_list['Gene Name'].tolist()
genes_present = [gene for gene in gene_names if gene in adata.var_names]

# Subset data
adata_subset = adata_subset[:, genes_present].copy()

# Create group label
adata_subset.obs['group'] = (
    adata_subset.obs['celltype3_renamed'].astype(str) + "_" +
    adata_subset.obs['Sex'].astype(str) + "_" +
    adata_subset.obs['Dx_OUD'].astype(str)
)

# Average expression per group
mean_expr = adata_subset.to_df().groupby(adata_subset.obs['group']).mean().T

# Z-score normalize using sklearn StandardScaler
scaler = StandardScaler()
mean_expr_z = pd.DataFrame(
    scaler.fit_transform(mean_expr.T).T,
    index=mean_expr.index,
    columns=mean_expr.columns
)

# Extract group components
group_split = mean_expr.columns.str.split("_", expand=True)
group_split = pd.DataFrame(group_split.tolist(), index=mean_expr.columns, columns=['Cell Type', 'Sex', 'Diagnosis'])

# Custom colors for Diagnosis and Sex
diagnosis_palette = {
    'None': '#d95f02',  # Burnt orange
    'OUD': '#1b9e77'    # Teal
}

sex_palette = {
    'F': '#7570b3',  # Purple
    'M': '#e7298a'   # Magenta
}

# Define cell type colors using Set3 and extend it for 15 categories
set3_base = sns.color_palette("Set3", 12)
extra_pastels = ["#e6ccff", "#ffcccc", "#ccffe6"]  # Light pastel purple, pink, green
extended_set3 = set3_base + extra_pastels

unique_celltypes = group_split['Cell Type'].unique()
celltype_palette = dict(zip(
    unique_celltypes,
    extended_set3[:len(unique_celltypes)]
))

# Create col_colors DataFrame
col_colors = pd.DataFrame({
    'Diagnosis': group_split['Diagnosis'].map(diagnosis_palette),
    'Sex': group_split['Sex'].map(sex_palette),
    'Cell Type': group_split['Cell Type'].map(celltype_palette)
}, index=mean_expr.columns)

# Sort by Sex → Diagnosis → Cell Type
sorted_columns = group_split.sort_values(['Sex', 'Diagnosis', 'Cell Type']).index
mean_expr_z = mean_expr_z[sorted_columns]
col_colors = col_colors.loc[sorted_columns]

# Add Cell Type to col_colors
col_colors['Cell Type'] = group_split['Cell Type'].map(celltype_palette)
col_colors = col_colors[['Diagnosis', 'Sex', 'Cell Type']]  # Reorder if needed

# Plot clustermap
g = sns.clustermap(
    mean_expr_z,
    col_cluster=False,
    row_cluster=True,
    col_colors=col_colors,
    cmap=sns.diverging_palette(145, 300, s=80, as_cmap=True),
    center=0,
    figsize=(14, 11),
    xticklabels=False,
    yticklabels=True,
    dendrogram_ratio=(0.2, 0.2)  # Adjust this to control the size of the dendrograms (row, column)
)

# Adjust tick labels
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=8)

# Create legend patches
legend_patches_diag = [Patch(color=v, label=f'Diagnosis: {k}') for k, v in diagnosis_palette.items()]
legend_patches_sex = [Patch(color=v, label=f'Sex: {k}') for k, v in sex_palette.items()]
legend_patches_cell = [Patch(color=v, label=f'Cell Type: {k}') for k, v in celltype_palette.items()]

# Place legends separately

# Diagnosis + Sex
leg1 = g.figure.legend(
    handles=legend_patches_diag + legend_patches_sex,  # Specifies the handles (the items to be shown in the legend)
    title="Group Annotations",  # Title for this legend
    loc='upper right',  # Specifies the position of the legend. 'upper right' places it at the top right corner of the plot
    bbox_to_anchor=(1.00, 0.92),  # Adjusts the exact position of the legend.
                               # The first value (1.02) shifts the legend slightly right from the plot area (1.0 is the edge),
                               # The second value (1) places it at the top (1 is the top edge, 0 is the bottom).
                               # You can change these values to move the legend further or closer.
    frameon=False  # Removes the frame around the legend for a cleaner look
)

# Cell Types (move to top center)
leg2 = g.figure.legend(
    handles=legend_patches_cell,  # Specifies the handles (the items to be shown in the legend)
    title="Cell Types",  # Title for this legend
    loc='center',  # Specifies the position of the legend. 'center' places it in the center of the plot area.
    bbox_to_anchor=(0.52, 0.87),  # Adjusts the exact position of the legend.
                               # The first value (0.3) positions the legend closer to the left (1 would be fully to the right).
                               # The second value (1.0) keeps it at the top of the plot.
                               # Adjust these to fine-tune the top-center positioning.
    ncol=3,  # Specifies the number of columns for the legend items. Set to 3 for 3 columns.
    frameon=False  # Removes the frame around the legend for a cleaner look
)

# Adjust colorbar position manually (place it to the right of the plot)
cbar = g.ax_heatmap.collections[0].colorbar  # Access the colorbar for the heatmap
cbar.set_label('Z-Score Expression', rotation=270, labelpad=20)  # Adds a label to the colorbar, rotates it 270° (vertical), and adds padding
cbar.ax.tick_params(labelsize=8)  # Adjust the size of the colorbar ticks for better readability
cbar.ax.set_position([0.1, 0.8, 0.03, 0.15])  # Manually adjusts the colorbar's position:
                                                # [1.03] - x position of the colorbar (shift to the right)
                                                # [0.2] - y position of the colorbar (higher is closer to the top)
                                                # [0.03] - width of the colorbar
                                                # [0.3] - height of the colorbar


# Save plot
plt.show()
g.savefig("heatmap_with_separated_legends_and_adjusted_colorbar.png", dpi=300, bbox_inches='tight')
