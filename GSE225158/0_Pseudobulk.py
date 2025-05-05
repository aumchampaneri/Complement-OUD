import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE233279/GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5ad')

if adata.raw is not None:
    counts = adata.raw.X.T  # Transpose: now (n_genes, n_cells)
    var_names = adata.raw.var['features']
else:
    counts = adata.X.T
    var_names = adata.var['features']

group_cols = ['ID', 'celltype3']
adata.obs['group'] = adata.obs[group_cols].agg('_'.join, axis=1)
groups = adata.obs['group'].astype('category')
group_names = groups.cat.categories

pseudobulk = []
for group in group_names:
    idx = np.where(groups == group)[0]
    summed = counts[:, idx].sum(axis=1)
    pseudobulk.append(np.asarray(summed).flatten())
pseudobulk = np.stack(pseudobulk, axis=1)

pseudobulk_df = pd.DataFrame(
    pseudobulk, index=var_names, columns=group_names
)
pseudobulk_df.index.name = 'gene'

meta = adata.obs.drop_duplicates('group').set_index('group')
meta = meta.loc[group_names, ['ID', 'celltype3', 'Dx_OUD', 'Sex', 'Age', 'PMI']]

pseudobulk_df.to_csv('counts.csv')
meta.to_csv('meta.csv')