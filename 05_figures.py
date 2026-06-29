"""
05_figures.py  --  paper figures: volcano + CD8 UMAP.

  volcano.png         -- responder vs non-responder DE; signature genes in red
  umap_cd8_states.png -- CD8 cells colored by response + key marker genes

Run:  .venv\\Scripts\\python.exe 05_figures.py
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import h5py
import matplotlib.pyplot as plt

MEL_H5   = "data/SKCM_GSE120575_aPD1aCTLA4_expression.h5"
MEL_META = "data/SKCM_GSE120575_aPD1aCTLA4_CellMetainfo_table.tsv"
CELLTYPE_COL = "Celltype (major-lineage)"
FOCUS = "CD8"
sc.settings.figdir = "figures"


def _decode(a):
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in a]


def load(h5, meta):
    with h5py.File(h5, "r") as f:
        g = f["matrix"]
        X = sp.csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]),
                          shape=tuple(int(x) for x in g["shape"][:])).T.tocsr()
        bc, gn = _decode(g["barcodes"][:]), _decode(g["features/name"][:])
    a = ad.AnnData(X=X)
    a.obs_names = bc
    a.var_names = gn
    a.var_names_make_unique()
    m = pd.read_csv(meta, sep="\t").set_index("Cell")
    common = a.obs_names.intersection(m.index)
    a = a[common].copy()
    a.obs = m.loc[a.obs_names].copy()
    if float(a.X.max()) > 50:
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
    return a


# baseline CD8 T cells
mel = load(MEL_H5, MEL_META)
mel = mel[mel.obs["TimePoint"].astype(str) == "Pre"].copy()
cd8 = mel[mel.obs[CELLTYPE_COL].astype(str).str.contains(FOCUS, case=False, na=False)].copy()
print("CD8 cells:", cd8.n_obs)

# ---- volcano: responder vs non-responder ----
sc.tl.rank_genes_groups(cd8, "Response", groups=["Responder"],
                        reference="Non-responder", method="wilcoxon")
df = sc.get.rank_genes_groups_df(cd8, group="Responder")
df["nlp"] = -np.log10(df["pvals_adj"].clip(lower=1e-300))
sig = {g.strip() for g in open("response_signature.txt") if g.strip()}

plt.figure(figsize=(5, 4))
plt.scatter(df["logfoldchanges"], df["nlp"], s=4, c="lightgrey")
hit = df[df["names"].isin(sig)]
plt.scatter(hit["logfoldchanges"], hit["nlp"], s=12, c="crimson")
for _, r in hit.nlargest(8, "nlp").iterrows():
    plt.annotate(r["names"], (r["logfoldchanges"], r["nlp"]), fontsize=7)
plt.axvline(0, c="k", lw=0.5)
plt.xlabel("log2 fold-change (responder vs non-responder)")
plt.ylabel("-log10 adjusted p")
plt.title("Baseline CD8 T: responder vs non-responder")
plt.tight_layout()
plt.savefig("figures/volcano.png", dpi=150)
plt.close()
print("Saved -> figures/volcano.png")

# ---- UMAP of CD8 states ----
sc.tl.pca(cd8, n_comps=30)
sc.pp.neighbors(cd8, n_neighbors=15, n_pcs=30)
sc.tl.umap(cd8)
sc.pl.umap(cd8, color=["Response", "TCF7", "CD8A"], save="_cd8_states.png", show=False)
print("Saved -> figures/umap_cd8_states.png")
