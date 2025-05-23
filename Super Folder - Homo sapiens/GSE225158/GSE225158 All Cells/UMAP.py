import scanpy as sc
'''
This script generates UMAP plots for the GSE225158 All Cells dataset.
'''

# Load the data
adata = sc.read_h5ad("/Users/aumchampaneri/PycharmProjects/Complement-OUD/Super Series - GSE233279/GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5ad")

# # Rename the cell types to more descriptive names
# rename_dict = {
#     "Microglia": "Microglia",
#     "Oligos": "Oligodendrocytes",
#     "Oligos_Pre": "OPCs",
#     "Astrocytes": "Astrocytes",
#     "Endothelial": "Endothelial cells",
#     "Mural": "Mural cells",
#     "D1-Matrix": "D1 MSNs (Matrix)",
#     "D2-Matrix": "D2 MSNs (Matrix)",
#     "D1-Striosome": "D1 MSNs (Striosome)",
#     "D2-Striosome": "D2 MSNs (Striosome)",
#     "D1/D2-Hybrid": "D1/D2 Hybrid MSNs",
#     "Int-PTHLH": "PTHLH+ Interneurons",
#     "Int-CCK": "CCK+ Interneurons",
#     "Int-SST": "SST+ Interneurons",
#     "Int-TH": "TH+ Interneurons"
# }
# adata.obs['celltype3_renamed'] = adata.obs['celltype3'].map(rename_dict)

# # Visualize broad cell type structure
# sc.pl.umap(adata, color='celltype3_renamed', legend_loc='on data', legend_fontsize='xx-small', legend_fontweight='light', title='Cell Types', save='broad_cell_types.pdf')
# sc.pl.umap(adata, color='celltype3_renamed', title='Cell Types', save='broad_cell_types_side_legend.pdf')
#
# # Check for sex-related structure or batch effects
# sc.pl.umap(adata, color='Sex', title='Sex', save='sex.png')

# Assess disease state distribution across UMAP
sc.pl.umap(adata, color='Dx_OUD', title='OUD Diagnosis', save='OUD_diagnosis.pdf')

# # Assess Tissue Types distribution across UMAP.
# sc.pl.umap(adata, color='Region', title='Brain Region', save='brain_region.png')