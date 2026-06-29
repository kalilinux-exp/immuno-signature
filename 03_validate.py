"""
03_validate.py  --  cross-cancer validation of the melanoma responder signature.

Question: does the signature built on MELANOMA (GSE120575) ALSO separate
responders from non-responders in a DIFFERENT cancer -- basal/squamous cell
carcinoma (Yost et al. 2019, GSE123813)? If yes => generalizable biomarker.

Download from TISCH2 (Skin -> BCC), dataset GSE123813, the two files into data/:
    <name>_expression.h5
    <name>_CellMetainfo_table.tsv
Then set the two *_PATH values + check the column config below.
If a column/label is wrong, the script PRINTS the available ones and stops.

Run:  .venv\\Scripts\\python.exe 03_validate.py
"""

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# CONFIG  -- set the two paths to YOUR downloaded BCC files, then run.
# ---------------------------------------------------------------------------
H5_PATH   = "data/BCC_GSE123813_aPD1_expression.h5"            # <-- adjust to real filename
META_PATH = "data/BCC_GSE123813_aPD1_CellMetainfo_table.tsv"   # <-- adjust to real filename
SIGNATURE_FILE = "response_signature.txt"

RESPONSE_COL       = "Response"
CELLTYPE_COL       = "Celltype (major-lineage)"
PATIENT_COL        = "Patient"
RESPONDER_LABEL    = "R"            # BCC labels responders "R"
NONRESPONDER_LABEL = "NR"           # ...and non-responders "NR"
FOCUS_CELLTYPE     = "CD8"
TIMEPOINT_COL      = "TimePoint"    # used only if present
TIMEPOINT_KEEP     = "pre"          # BCC uses lowercase "pre"

sc.settings.verbosity = 1
sc.settings.figdir = "figures"


def _decode(arr):
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in arr]


def load_matrix(h5_path):
    with h5py.File(h5_path, "r") as f:
        g = f["matrix"]
        data, indices, indptr = g["data"][:], g["indices"][:], g["indptr"][:]
        shape = tuple(int(x) for x in g["shape"][:])
        barcodes = _decode(g["barcodes"][:])
        genes = _decode(g["features/name"][:])
    mat = sp.csc_matrix((data, indices, indptr), shape=shape).T.tocsr()
    adata = ad.AnnData(X=mat)
    adata.obs_names = barcodes
    adata.var_names = genes
    adata.var_names_make_unique()
    print(f"Matrix: {adata.n_obs} cells x {adata.n_vars} genes")
    return adata


def attach_metadata(adata, meta_path):
    meta = pd.read_csv(meta_path, sep="\t").set_index("Cell")
    common = adata.obs_names.intersection(meta.index)
    print(f"Barcode overlap: {len(common)} / {adata.n_obs}")
    adata = adata[common].copy()
    adata.obs = meta.loc[adata.obs_names].copy()
    return adata


def ensure_lognorm(adata):
    if float(adata.X.max()) > 50:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    return adata


def check_config(adata):
    """Stop early with a helpful printout if the column config doesn't match."""
    ok = True
    if RESPONSE_COL not in adata.obs.columns or CELLTYPE_COL not in adata.obs.columns:
        ok = False
    elif RESPONDER_LABEL not in set(adata.obs[RESPONSE_COL].astype(str)):
        ok = False
    if not ok:
        print("\n!!! CONFIG MISMATCH -- here is what's actually in this file:")
        for c in adata.obs.columns:
            u = adata.obs[c].unique()
            if len(u) <= 20:
                print(f"  {c}: {list(u)}")
            else:
                print(f"  {c}: ({len(u)} unique)")
        print("\nFix RESPONSE_COL / CELLTYPE_COL / RESPONDER_LABEL above, then re-run.")
    return ok


def main():
    signature = [g.strip() for g in open(SIGNATURE_FILE) if g.strip()]
    print(f"Loaded {len(signature)}-gene melanoma signature")

    adata = load_matrix(H5_PATH)
    adata = attach_metadata(adata, META_PATH)
    if not check_config(adata):
        return
    adata = ensure_lognorm(adata)

    # BCC tags Response only on post-treatment cells, but response is a per-PATIENT
    # outcome -> propagate each patient's label to ALL their cells (incl. baseline/pre)
    labeled = adata.obs.dropna(subset=[RESPONSE_COL])
    pat_resp = labeled.groupby(PATIENT_COL)[RESPONSE_COL].agg(lambda s: s.mode().iloc[0])
    adata.obs[RESPONSE_COL] = adata.obs[PATIENT_COL].map(pat_resp).astype(object)
    print("patient -> response:", pat_resp.to_dict())

    if TIMEPOINT_KEEP is not None and TIMEPOINT_COL in adata.obs.columns:
        adata = adata[adata.obs[TIMEPOINT_COL] == TIMEPOINT_KEEP].copy()
        print(f"after TimePoint=={TIMEPOINT_KEEP}: {adata.n_obs} cells")

    # focus on CD8 T cells, score them with the melanoma signature
    mask = adata.obs[CELLTYPE_COL].astype(str).str.contains(FOCUS_CELLTYPE, case=False, na=False)
    cd8 = adata[mask].copy()
    # keep only cells with a known response label (BCC has NaN-labeled cells)
    valid = cd8.obs[RESPONSE_COL].astype(str).isin([RESPONDER_LABEL, NONRESPONDER_LABEL])
    cd8 = cd8[valid].copy()
    print(f"CD8 cells with response label: {cd8.n_obs}")
    present = [g for g in signature if g in cd8.var_names]
    print(f"signature genes present in BCC data: {len(present)}/{len(signature)}")
    sc.tl.score_genes(cd8, present, score_name="resp_sig")

    y_cell = (cd8.obs[RESPONSE_COL].astype(str) == RESPONDER_LABEL).astype(int).values
    auc_cell = roc_auc_score(y_cell, cd8.obs["resp_sig"].values)
    print(f"\nCELL-level AUC  = {auc_cell:.3f}")

    # patient-level (the rigorous metric -- one value per patient)
    pt = cd8.obs.groupby(PATIENT_COL).agg(
        score=("resp_sig", "mean"),
        resp=(RESPONSE_COL, "first")).dropna()
    y_pt = (pt["resp"].astype(str) == RESPONDER_LABEL).astype(int).values
    if len(set(y_pt)) == 2:
        auc_pt = roc_auc_score(y_pt, pt["score"].values)
        print(f"PATIENT-level AUC = {auc_pt:.3f}  (n={len(pt)} patients)")
    else:
        auc_pt = None
        print("Only one response class among patients -- can't compute patient AUC")

    # figures: violin + ROC
    sc.pl.violin(cd8, "resp_sig", groupby=RESPONSE_COL, save="_bcc_score.png", show=False)
    fpr, tpr, _ = roc_curve(y_cell, cd8.obs["resp_sig"].values)
    plt.figure(figsize=(4, 4))
    plt.plot(fpr, tpr, label=f"cell AUC={auc_cell:.2f}")
    plt.plot([0, 1], [0, 1], "k--", lw=0.8)
    plt.xlabel("False positive rate"); plt.ylabel("True positive rate")
    plt.title("Melanoma signature -> BCC response"); plt.legend()
    plt.tight_layout(); plt.savefig("figures/roc_validation.png", dpi=150)
    print("Saved -> figures/violin_bcc_score.png, figures/roc_validation.png")

    print("\nVERDICT: AUC > 0.5 means the melanoma signature carries cross-cancer "
          "signal. The closer to 1.0, the stronger your headline finding.")


if __name__ == "__main__":
    main()
