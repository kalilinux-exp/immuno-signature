"""
04_pseudobulk.py  --  rigorous (pseudobulk) signature + cross-cancer validation.

Upgrade over 02/03: instead of treating each CELL as an independent sample
(pseudoreplication -- inflates significance), aggregate each PATIENT's baseline
CD8 T cells into ONE profile, then do differential expression at the PATIENT
level. Then validate the resulting signature on BCC, same metric as 03.

Run:  .venv\\Scripts\\python.exe 04_pseudobulk.py
"""

import re
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import h5py
from scipy.stats import ttest_ind
from sklearn.metrics import roc_auc_score

# ---- files ----
MEL_H5   = "data/SKCM_GSE120575_aPD1aCTLA4_expression.h5"
MEL_META = "data/SKCM_GSE120575_aPD1aCTLA4_CellMetainfo_table.tsv"
BCC_H5   = "data/BCC_GSE123813_aPD1_expression.h5"
BCC_META = "data/BCC_GSE123813_aPD1_CellMetainfo_table.tsv"

CELLTYPE_COL = "Celltype (major-lineage)"
PATIENT_COL  = "Patient"
FOCUS = "CD8"
N_SIG = 50

# ---- noise filter (same as 02/03) ----
_NOISE_PREFIXES = ("RPS", "RPL", "MRPS", "MRPL", "MT-")
_NOISE_EXACT = {"MALAT1", "NEAT1", "XIST", "TMSB4X", "TMSB10", "FOS", "FOSB",
                "JUN", "JUNB", "JUND", "EGR1", "HSPA1A", "HSPA1B", "HSPA8", "DNAJB1"}
_CLONE = re.compile(r"^(AC|AL|AP|LINC)\d{3,}")


def is_noise_gene(g):
    g = g.upper()
    return (g.startswith(_NOISE_PREFIXES) or "." in g or bool(_CLONE.match(g))
            or g in _NOISE_EXACT)


def _decode(arr):
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in arr]


def load(h5, meta):
    """Load 10x HDF5 + merge metadata + log-normalize."""
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


def cd8_at(a, timepoint):
    """Subset to baseline (pre) CD8 T cells."""
    a = a[a.obs["TimePoint"].astype(str) == timepoint].copy()
    mask = a.obs[CELLTYPE_COL].astype(str).str.contains(FOCUS, case=False, na=False)
    return a[mask].copy()


def pseudobulk(a, resp_map):
    """Mean log-expression per patient -> patients x genes, labelled by response."""
    pats = a.obs[PATIENT_COL].astype(str)
    cats = sorted(pats.unique())
    M = np.vstack([np.asarray(a.X[(pats == c).values].mean(axis=0)).ravel() for c in cats])
    pb = ad.AnnData(X=M)
    pb.obs_names = cats
    pb.var_names = a.var_names
    pb.obs["resp"] = [resp_map.get(c, None) for c in cats]
    return pb[pb.obs["resp"].notna()].copy()


# ==== 1. melanoma -> pseudobulk -> patient-level DE -> signature ====
print("=== melanoma (discovery) ===")
mel = load(MEL_H5, MEL_META)
mel_lab = mel.obs.groupby(PATIENT_COL)["Response"].agg(lambda s: s.dropna().mode().iloc[0])
mel_resp = mel_lab.to_dict()
mel_cd8 = cd8_at(mel, "Pre")
pb = pseudobulk(mel_cd8, mel_resp)
print("pseudobulk patients:", pb.obs["resp"].value_counts().to_dict())

resp = pb.X[(pb.obs["resp"] == "Responder").values]
nonr = pb.X[(pb.obs["resp"] == "Non-responder").values]
t, _ = ttest_ind(resp, nonr, axis=0, equal_var=False)
t = np.nan_to_num(t, nan=-np.inf)
# floor on MEAN expression magnitude -- real genes (TCF7~1-3) sit well above the
# low-expression flukes (~0.01-0.1) that get huge t-stats at small n (n=19)
expressed = np.asarray(pb.X.mean(axis=0)).ravel() >= 0.5
t = np.where(expressed, t, -np.inf)
print(f"expressed genes kept: {int(expressed.sum())}")
ranked = [str(pb.var_names[i]) for i in np.argsort(-t)]
signature = [g for g in ranked if not is_noise_gene(g)][:N_SIG]
print(f"\npseudobulk signature ({len(signature)} genes):\n{signature}")
with open("pseudobulk_signature.txt", "w") as fh:
    fh.write("\n".join(signature))
print("Saved -> pseudobulk_signature.txt")

# ==== 2. validate on BCC (patient-level AUC, same as 03) ====
print("\n=== BCC (validation) ===")
bcc = load(BCC_H5, BCC_META)
bcc_lab = bcc.obs.dropna(subset=["Response"]).groupby(PATIENT_COL)["Response"].agg(lambda s: s.mode().iloc[0])
bcc.obs["Response"] = bcc.obs[PATIENT_COL].map(bcc_lab.to_dict())
bcc_cd8 = cd8_at(bcc, "pre")
bcc_cd8 = bcc_cd8[bcc_cd8.obs["Response"].astype(str).isin(["R", "NR"])].copy()
present = [g for g in signature if g in bcc_cd8.var_names]
print(f"signature genes present in BCC: {len(present)}/{len(signature)}")
sc.tl.score_genes(bcc_cd8, present, score_name="sig")
pt = bcc_cd8.obs.groupby(PATIENT_COL).agg(score=("sig", "mean"), resp=("Response", "first"))
y = (pt["resp"] == "R").astype(int).values
auc = roc_auc_score(y, pt["score"].values)
print(f"\nBCC PATIENT-level AUC (pseudobulk signature) = {auc:.3f}  (n={len(pt)} patients)")
print("(compare to the cell-level signature's 0.70)")
