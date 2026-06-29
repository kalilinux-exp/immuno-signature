"""
06_rigor.py  --  statistical honesty for the cross-cancer AUC.

For each signature, at the PATIENT level (the honest unit, n=11):
  - point AUC
  - bootstrap 95% confidence interval (resample patients) -- shows how shaky
  - permutation p-value (shuffle labels) -- is it better than chance?

Small n => expect a WIDE CI. Reporting that yourself = the maturity signal.

Run:  .venv\\Scripts\\python.exe 06_rigor.py
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import h5py
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score

BCC_H5   = "data/BCC_GSE123813_aPD1_expression.h5"
BCC_META = "data/BCC_GSE123813_aPD1_CellMetainfo_table.tsv"
SIGS = {"cell-level": "response_signature.txt", "pseudobulk": "pseudobulk_signature.txt"}
CELLTYPE_COL = "Celltype (major-lineage)"
PATIENT_COL = "Patient"
FOCUS = "CD8"
N_BOOT, N_PERM = 2000, 10000
rng = np.random.default_rng(0)
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


def patient_scores(adata, signature):
    present = [g for g in signature if g in adata.var_names]
    sc.tl.score_genes(adata, present, score_name="sig")
    pt = adata.obs.groupby(PATIENT_COL).agg(score=("sig", "mean"), resp=("Response", "first"))
    y = (pt["resp"] == "R").astype(int).values
    return y, pt["score"].values, len(present), len(signature)


def bootstrap_ci(y, s):
    idx = np.arange(len(y))
    aucs = []
    for _ in range(N_BOOT):
        b = rng.choice(idx, size=len(idx), replace=True)
        if len(set(y[b])) == 2:           # need both classes to define AUC
            aucs.append(roc_auc_score(y[b], s[b]))
    return np.percentile(aucs, [2.5, 97.5]), np.array(aucs)


def perm_p(y, s, obs):
    perm = np.array([roc_auc_score(rng.permutation(y), s) for _ in range(N_PERM)])
    return (np.sum(perm >= obs) + 1) / (N_PERM + 1), perm


# ---- load BCC once, build baseline pre CD8 with labels ----
bcc = load(BCC_H5, BCC_META)
lab = bcc.obs.dropna(subset=["Response"]).groupby(PATIENT_COL)["Response"].agg(lambda x: x.mode().iloc[0])
bcc.obs["Response"] = bcc.obs[PATIENT_COL].map(lab.to_dict())
bcc = bcc[bcc.obs["TimePoint"].astype(str) == "pre"].copy()
cd8 = bcc[bcc.obs[CELLTYPE_COL].astype(str).str.contains(FOCUS, case=False, na=False)].copy()
cd8 = cd8[cd8.obs["Response"].astype(str).isin(["R", "NR"])].copy()
print(f"BCC pre CD8 cells: {cd8.n_obs}, patients: {cd8.obs[PATIENT_COL].nunique()}")

for name, fn in SIGS.items():
    sig = [g.strip() for g in open(fn) if g.strip()]
    y, s, present, total = patient_scores(cd8.copy(), sig)
    obs = roc_auc_score(y, s)
    (lo, hi), boot = bootstrap_ci(y, s)
    p, perm = perm_p(y, s, obs)
    print(f"\n=== {name} signature ({present}/{total} genes present) ===")
    print(f"patients: n={len(y)}  (R={int(y.sum())}, NR={int((1 - y).sum())})")
    print(f"AUC = {obs:.3f}   95% CI [{lo:.3f}, {hi:.3f}]   permutation p = {p:.3f}")
    if name == "pseudobulk":
        plt.figure(figsize=(5, 3))
        plt.hist(perm, bins=30, color="lightgrey", label="null (shuffled labels)")
        plt.axvline(obs, c="crimson", label=f"observed {obs:.2f}")
        plt.xlabel("AUC"); plt.ylabel("count"); plt.legend()
        plt.tight_layout(); plt.savefig("figures/permutation_null.png", dpi=150); plt.close()
        print("Saved -> figures/permutation_null.png")
