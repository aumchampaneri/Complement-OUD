# run_ssgsea_pipeline.py

import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse
import gseapy as gp
import seaborn as sns
import matplotlib.pyplot as plt

# === USER CONFIGURATION ===
ANADTA_PATH   = "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE233279/GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5ad"
GENE_LIST_CSV = "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158 Gene Library/gene_names.csv"
# ==========================

# 1. Load AnnData
adata = sc.read_h5ad(ANADTA_PATH)

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
adata.obs['celltype'] = adata.obs['celltype3'].map(rename_dict)

# 2. Build grouping key
adata.obs['group'] = adata.obs[['Dx_OUD','Sex','celltype']].agg('_'.join, axis=1)

# 3. Subset genes
genes = pd.read_csv(GENE_LIST_CSV)['Gene Name'].tolist()
adata = adata[:, adata.var_names.isin(genes)].copy()

# 4. Extract expression matrix
expr = adata.X.toarray() if scipy.sparse.issparse(adata.X) else adata.X
expr_df = pd.DataFrame(expr, index=adata.obs_names, columns=adata.var_names)

# 5. Pseudobulk & CPM
pb = expr_df.groupby(adata.obs['group']).sum().T  # genes × samples
cpm = pb.div(pb.sum(axis=0), axis=1) * 1e6

# 6. Define pathways
pathway_dict = {
    "classical_pathway": ["C1QA","C1QB","C1QC","C1R","C1S","C2","C4A","C4B","C3"],
    "alternative_pathway": ["CFB","CFD","CFP","CFH","CFI","CR1","CR2"],
    "complosome": ["C3","C5","CTSL","C3AR1","C5AR1","C5AR2","CD46"],
    "nfkb_pathway": ["NFKB1","RELA","NFKB2","RELB","NFKBIA"],
    "inflammasome": ["NLRP3","PYCARD","CASP1","IL1B","IL18","IL1R1","IL1R2","IL1RN"],
    "tnf_pathway": ["TNF","TNFRSF1A","TNFRSF1B","TNFAIP3","MAPK8","MAPK14"],
    "interferon_pathway": ["IFNG","IFNGR1","IFNGR2","STAT1","IRF1","IRF7",
                           "IFNA2","IFNA5","IFNA6","IFNA13","IFNA21","IFNW1",
                           "IL6","IL10"]
}

# 7. Run ssGSEA with relaxed size filters
ss = gp.ssgsea(
    data=cpm,
    gene_sets=pathway_dict,
    sample_norm_method="rank",
    outdir=None,
    permutation_num=0,
    no_plot=True,
    min_size=1,
    max_size=1000
)

# 8. Extract results (handle different versions)
if hasattr(ss, 'resultsOnSamples'):
    df = ss.resultsOnSamples
elif hasattr(ss, 'res2d'):
    df = ss.res2d
elif hasattr(ss, 'results'):
    df = ss.results
else:
    raise AttributeError("Can't find ssGSEA results on the SingleSampleGSEA object")

# 9. Pivot into a numeric matrix
matrix = df.pivot(index='Name', columns='Term', values='NES')

# 9b. Coerce everything to float
matrix = matrix.apply(pd.to_numeric, errors='coerce')
matrix = matrix.dropna(axis=0, how='all').dropna(axis=1, how='all')

print("Final numeric matrix shape:", matrix.shape)
print(matrix.dtypes.unique())  # should show only float64

# --- AFTER pivoting & coercion, before sorting & plotting ---

# Build metadata DataFrame from the matrix index
meta = pd.DataFrame(index=matrix.index)

# Split only on the first two underscores → exactly 3 columns
parts = meta.index.to_series().str.split('_', n=2, expand=True)
parts.columns = ['Dx_OUD','Sex','celltype']

# Attach to meta
meta = meta.join(parts)

# Now sort however you like
meta_sorted = meta.sort_values(['Sex','Dx_OUD','celltype'])

# Reindex the matrix
matrix_sorted = matrix.loc[meta_sorted.index]

# Plot without reclustering rows
sns.clustermap(
    matrix_sorted,
    cmap="vlag",
    z_score=0,
    figsize=(10, 6),
    row_cluster=False,
    col_cluster=True
)
plt.title("ssGSEA Heatmap (sorted by Dx → CellType → Sex)")
plt.show()
