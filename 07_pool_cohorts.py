"""
07_pool_cohorts.py  --  validate the signature across MANY anti-PD-1 cohorts and
pool them for a larger n, a tighter confidence interval, and a stronger p-value.

The signature is NOT rebuilt here. It is the melanoma-derived signature, scored on
each validation cohort. Each TISCH dataset labels things slightly differently, so each
cohort gets one small CONFIG entry below. Add a new dataset by downloading its two files
into data/ and appending a dict. Run it; if a column or label is wrong it prints what is
actually in the file so you can fix the entry.

Scores are standardized within each cohort before pooling, so a dataset's overall
brightness does not bias the pooled AUC (a basic guard against batch effects).

Run:  .venv\\Scripts\\python.exe 07_pool_cohorts.py
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import h5py
from sklearn.metrics import roc_auc_score

SIGNATURE_FILE = "pseudobulk_signature.txt"
FOCUS = "CD8"
PATIENT_COL = "Patient"
N_BOOT, N_PERM = 2000, 10000
rng = np.random.default_rng(0)

# ---------------------------------------------------------------------------
# COHORTS  --  one entry per validation dataset. The first is your existing BCC.
# To add a cohort: download its 2 TISCH files into data/, copy an entry, and fill
# the column/label values. Run once; mismatches print the real columns to copy from.
# ---------------------------------------------------------------------------
COHORTS = [
    dict(name="BCC_GSE123813",
         h5="data/BCC_GSE123813_aPD1_expression.h5",
         meta="data/BCC_GSE123813_aPD1_CellMetainfo_table.tsv",
         celltype_col="Celltype (major-lineage)",
         response_col="Response", responder="R", nonresponder="NR",
         timepoint_col="TimePoint", timepoint_keep="pre"),

    dict(name="SCC_GSE123813",
         h5="data/SCC_GSE123813_aPD1_expression.h5",
         meta="data/SCC_GSE123813_aPD1_CellMetainfo_table.tsv",
         celltype_col="Celltype (major-lineage)",
         response_col="Response", responder="R", nonresponder="NR",
         timepoint_col="TimePoint", timepoint_keep="pre"),

    # Add more cohorts here (download 2 files into data/, copy the block above, adjust).
]


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


def show_columns(name, a):
    print(f"  {name}: config mismatch. Columns actually present:")
    for c in a.obs.columns:
        u = a.obs[c].unique()
        if len(u) <= 20:
            print(f"    {c}: {list(u)}")
        else:
            print(f"    {c}: ({len(u)} unique)")


def process_cohort(c, signature):
    """Return (within-cohort standardized patient scores, patient labels 0/1, AUC) or None."""
    a = load(c["h5"], c["meta"])
    for col in (c["celltype_col"], c["response_col"], PATIENT_COL):
        if col not in a.obs.columns:
            show_columns(c["name"], a)
            return None

    # response is a per-patient outcome; some datasets only label some cells (e.g. post),
    # so propagate each patient's label to all of their cells
    lab = a.obs.dropna(subset=[c["response_col"]]).groupby(PATIENT_COL)[c["response_col"]].agg(
        lambda s: s.mode().iloc[0])
    a.obs[c["response_col"]] = a.obs[PATIENT_COL].map(lab.to_dict())

    if c.get("timepoint_keep") and c.get("timepoint_col") in a.obs.columns:
        a = a[a.obs[c["timepoint_col"]].astype(str) == c["timepoint_keep"]].copy()

    a = a[a.obs[c["celltype_col"]].astype(str).str.contains(FOCUS, case=False, na=False)].copy()
    a = a[a.obs[c["response_col"]].astype(str).isin([c["responder"], c["nonresponder"]])].copy()
    if a.n_obs == 0:
        print(f"  {c['name']}: 0 usable cells after filtering. Check labels/timepoint in config.")
        return None

    present = [g for g in signature if g in a.var_names]
    sc.tl.score_genes(a, present, score_name="sig")
    pt = a.obs.groupby(PATIENT_COL).agg(score=("sig", "mean"), resp=(c["response_col"], "first"))
    y = (pt["resp"].astype(str) == c["responder"]).astype(int).values
    s = pt["score"].values
    if len(set(y)) < 2:
        print(f"  {c['name']}: only one response class among patients. Skipping.")
        return None

    auc = roc_auc_score(y, s)
    s_std = (s - s.mean()) / (s.std() + 1e-9)            # standardize within cohort
    print(f"  {c['name']}: n={len(y)} patients (R={int(y.sum())}, NR={int((1 - y).sum())}), "
          f"AUC={auc:.3f}, genes {len(present)}/{len(signature)}")
    return s_std, y, auc


def main():
    signature = [g.strip() for g in open(SIGNATURE_FILE) if g.strip()]
    print(f"Signature: {len(signature)} genes\nPer-cohort results:")

    z_all, y_all, per_cohort = [], [], []
    for c in COHORTS:
        out = process_cohort(c, signature)
        if out is not None:
            s_std, y, auc = out
            z_all.append(s_std)
            y_all.append(y)
            per_cohort.append((c["name"], auc, len(y)))

    if not per_cohort:
        print("\nNo usable cohorts. Fix the configs above and re-run.")
        return

    z = np.concatenate(z_all)
    y = np.concatenate(y_all)
    pooled = roc_auc_score(y, z)

    idx = np.arange(len(y))
    boot = []
    for _ in range(N_BOOT):
        b = rng.choice(idx, len(idx), replace=True)
        if len(set(y[b])) == 2:
            boot.append(roc_auc_score(y[b], z[b]))
    lo, hi = np.percentile(boot, [2.5, 97.5])
    perm = np.array([roc_auc_score(rng.permutation(y), z) for _ in range(N_PERM)])
    p = (np.sum(perm >= pooled) + 1) / (N_PERM + 1)

    print(f"\nPOOLED across {len(per_cohort)} cohort(s): n={len(y)} patients "
          f"(R={int(y.sum())}, NR={int((1 - y).sum())})")
    print(f"pooled AUC = {pooled:.3f}   95% CI [{lo:.3f}, {hi:.3f}]   permutation p = {p:.4f}")
    print("\nThe more cohorts you add, the larger n gets and the tighter the interval. "
          "Watch whether the signal holds or fades. Both are honest outcomes worth reporting.")


if __name__ == "__main__":
    main()
