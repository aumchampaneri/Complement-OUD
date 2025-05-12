import os
import scanpy as sc
import liana as li
import scipy.sparse as sp

# ─── Load & subset ───────────────────────────────────────────
adata = sc.read_h5ad("/Users/aumchampaneri/PycharmProjects/Complement-OUD/Super Series - GSE233279/GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5ad")

adata_subset = adata[
    adata.obs['Dx_OUD'].isin(['OUD','None']) &
    adata.obs['Sex'].isin(['M','F']),
    :
].copy()

glia = ['Microglia','Oligos','Oligos_Pre','Astrocytes']
adata_glial = adata_subset[adata_subset.obs['celltype3'].isin(glia), :].copy()

adata_glial.obs['celltype'] = adata_glial.obs['celltype3'].map({
    'Microglia':    'Microglia',
    'Oligos':       'Oligodendrocytes',
    'Oligos_Pre':   'OPCs',
    'Astrocytes':   'Astrocytes',
})

# Replace 'None' with 'CTL'
adata_glial.obs['Dx_OUD'] = adata_glial.obs['Dx_OUD'].replace('None','CTL')
adata_glial.obs['group'] = (
    adata_glial.obs['Dx_OUD'].astype(str) + "_" +
    adata_glial.obs['Sex'].astype(str)
).str.lower()

# ─── Materialize sparse view → CSR matrix ─────────────────────
adata_glial.X = sp.csr_matrix(adata_glial.X)

# ─── Run LIANA by_sample ───────────────────────────────────────
li.mt.rank_aggregate.by_sample(
    adata=adata_glial,
    sample_key='group',
    groupby='celltype',
    expr_prop=0.1,
    min_cells=5,
    resource_name='consensus',
    n_perms=100,
    use_raw=False,
    verbose=True
)

# ─── Save results ──────────────────────────────────────────────
out_dir = "LIANA+ output"
os.makedirs(out_dir, exist_ok=True)
adata_glial.write(os.path.join(out_dir, "Glial_LIANA_output.h5ad"))
print("Done — results saved to", out_dir)

