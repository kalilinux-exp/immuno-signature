"""
02_real_data.py  --  responder vs non-responder signature on REAL data.

Dataset: Sade-Feldman et al. 2018 melanoma anti-PD-1/CTLA-4 (GSE120575),
downloaded from TISCH2 as two files in data/:
    SKCM_GSE120575_aPD1aCTLA4_expression.h5           (10x-format gene matrix)
    SKCM_GSE120575_aPD1aCTLA4_CellMetainfo_table.tsv  (per-cell labels)

Pipeline:
  load matrix + merge metadata -> log-normalize -> keep BASELINE (pre-treatment)
  CD8 T cells -> differential expression responder vs non-responder
  -> filter technical noise genes -> top genes = response signature
  -> score every cell -> plot.

Run:  .venv\\Scripts\\python.exe 02_real_data.py
"""

import re
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc

# ---------------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------------
H5_PATH   = "data/SKCM_GSE120575_aPD1aCTLA4_expression.h5"
META_PATH = "data/SKCM_GSE120575_aPD1aCTLA4_CellMetainfo_table.tsv"

RESPONSE_COL       = "Response"
CELLTYPE_COL       = "Celltype (major-lineage)"
TIMEPOINT_COL      = "TimePoint"
PATIENT_COL        = "Patient"
RESPONDER_LABEL    = "Responder"
NONRESPONDER_LABEL = "Non-responder"

FOCUS_CELLTYPE    = "CD8"    # substring -> matches CD8T and CD8Tex (exhausted)
TIMEPOINT_KEEP    = "Pre"    # baseline only = predictive-biomarker story; None = all
N_SIGNATURE_GENES = 50

sc.settings.verbosity = 1
sc.settings.figdir = "figures"

# --- noise filter (a real research judgment call; edit to defend your choice) ---
# Genes that are technical artifacts in scRNA-seq, not biology:
_NOISE_PREFIXES = ("RPS", "RPL", "MRPS", "MRPL", "MT-")          # ribosomal + mitochondrial
_NOISE_EXACT = {"MALAT1", "NEAT1", "XIST", "TMSB4X", "TMSB10",   # lncRNA / sex / housekeeping
                "FOS", "FOSB", "JUN", "JUNB", "JUND", "EGR1",    # immediate-early / stress
                "HSPA1A", "HSPA1B", "HSPA8", "DNAJB1"}
_CLONE = re.compile(r"^(AC|AL|AP|LINC)\d{3,}")                   # uncharacterized clones/lncRNAs


def is_noise_gene(gene):
    """True = drop from signature. Ribosomal/mito/stress/clone artifacts."""
    g = gene.upper()
    if g.startswith(_NOISE_PREFIXES):
        return True
    if "." in g:                 # versioned clone names like AC016596.1
        return True
    if _CLONE.match(g):
        return True
    return g in _NOISE_EXACT


def _decode(arr):
    """h5py returns byte strings; decode to normal python str."""
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in arr]


def load_matrix(h5_path):
    """Read the TISCH 10x-format HDF5 into a cells x genes AnnData."""
    with h5py.File(h5_path, "r") as f:
        g = f["matrix"]
        data    = g["data"][:]
        indices = g["indices"][:]
        indptr  = g["indptr"][:]
        shape   = tuple(int(x) for x in g["shape"][:])   # (genes, cells)
        barcodes = _decode(g["barcodes"][:])
        genes    = _decode(g["features/name"][:])
    # 10x stores CSC as genes x cells; transpose to cells x genes
    mat = sp.csc_matrix((data, indices, indptr), shape=shape).T.tocsr()
    adata = ad.AnnData(X=mat)
    adata.obs_names = barcodes
    adata.var_names = genes
    adata.var_names_make_unique()
    print(f"Matrix: {adata.n_obs} cells x {adata.n_vars} genes")
    return adata


def attach_metadata(adata, meta_path):
    """Merge per-cell labels (response, cell type, ...) onto the matrix."""
    meta = pd.read_csv(meta_path, sep="\t").set_index("Cell")
    common = adata.obs_names.intersection(meta.index)
    print(f"Barcode overlap: {len(common)} / {adata.n_obs}")
    adata = adata[common].copy()
    adata.obs = meta.loc[adata.obs_names].copy()
    return adata


def ensure_lognorm(adata):
    mx = float(adata.X.max())
    print(f"max expression value = {mx:.1f}")
    if mx > 50:
        print("-> looks linear; normalize_total + log1p")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        print("-> already log-scaled; leaving as is")
    return adata


def build_signature(adata):
    # baseline (pre-treatment) only -> predict response BEFORE treatment starts
    if TIMEPOINT_KEEP is not None:
        adata = adata[adata.obs[TIMEPOINT_COL] == TIMEPOINT_KEEP].copy()
        print(f"after TimePoint=={TIMEPOINT_KEEP}: {adata.n_obs} cells")

    # focus on CD8 T cells (the effectors of anti-PD-1 response)
    mask = adata.obs[CELLTYPE_COL].astype(str).str.contains(FOCUS_CELLTYPE, case=False, na=False)
    sub = adata[mask].copy()
    print(f"\n{FOCUS_CELLTYPE} T cells: {sub.n_obs}")
    print("response breakdown:")
    print(sub.obs[RESPONSE_COL].value_counts())
    print(f"patients: {sub.obs[PATIENT_COL].nunique()}")

    # differential expression: genes UP in responders' CD8 T cells vs non-responders
    sc.tl.rank_genes_groups(sub, RESPONSE_COL, groups=[RESPONDER_LABEL],
                            reference=NONRESPONDER_LABEL, method="wilcoxon")
    names = sub.uns["rank_genes_groups"]["names"][RESPONDER_LABEL]

    # rank deep, then drop technical noise, then keep the top N real genes
    ranked = [str(n) for n in names[:300]]
    dropped = [g for g in ranked if is_noise_gene(g)]
    signature = [g for g in ranked if not is_noise_gene(g)][:N_SIGNATURE_GENES]
    print(f"\ndropped {len(dropped)} noise genes from top 300, e.g. {dropped[:12]}")
    print(f"\nTop {N_SIGNATURE_GENES}-gene CLEAN RESPONDER signature (CD8 T, baseline):")
    print(signature)

    with open("response_signature.txt", "w") as fh:
        fh.write("\n".join(signature))
    print("Saved -> response_signature.txt")

    # score every cell, then check responders score higher (a figure for your paper)
    sc.tl.score_genes(adata, signature, score_name="resp_sig")
    sc.pl.violin(adata, "resp_sig", groupby=RESPONSE_COL,
                 save="_signature_score.png", show=False)
    print("Saved -> figures/violin_signature_score.png")
    return signature


def main():
    adata = load_matrix(H5_PATH)
    adata = attach_metadata(adata, META_PATH)
    adata = ensure_lognorm(adata)
    build_signature(adata)
    print("\nNOTE (rigor): cells from the same patient aren't independent "
          "(pseudoreplication). A stronger version pseudobulks per patient before "
          "testing. This cell-level pass is the first look.")


if __name__ == "__main__":
    main()
