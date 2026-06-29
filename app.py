"""
app.py  --  ImmunoScore: score single-cell tumor data against a CD8 T-cell
signature linked to anti-PD-1 (checkpoint immunotherapy) response.

Run:  .venv\\Scripts\\streamlit run app.py
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import h5py
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score
import streamlit as st

st.set_page_config(page_title="ImmunoScore", layout="centered")

# hide Streamlit's default chrome for a cleaner, less-templated look
st.markdown(
    "<style>#MainMenu{visibility:hidden} footer{visibility:hidden} "
    "header{visibility:hidden} .block-container{padding-top:3rem;max-width:760px}</style>",
    unsafe_allow_html=True,
)

SIG_FILE = "pseudobulk_signature.txt"
BCC_H5 = "data/BCC_GSE123813_aPD1_expression.h5"
BCC_META = "data/BCC_GSE123813_aPD1_CellMetainfo_table.tsv"


@st.cache_data
def load_signature(path=SIG_FILE):
    return [g.strip() for g in open(path) if g.strip()]


def _decode(a):
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in a]


@st.cache_data
def load_demo():
    with h5py.File(BCC_H5, "r") as f:
        g = f["matrix"]
        X = sp.csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]),
                          shape=tuple(int(x) for x in g["shape"][:])).T.tocsr()
        bc, gn = _decode(g["barcodes"][:]), _decode(g["features/name"][:])
    a = ad.AnnData(X=X)
    a.obs_names = bc
    a.var_names = gn
    a.var_names_make_unique()
    m = pd.read_csv(BCC_META, sep="\t").set_index("Cell")
    common = a.obs_names.intersection(m.index)
    a = a[common].copy()
    a.obs = m.loc[a.obs_names].copy()
    lab = a.obs.dropna(subset=["Response"]).groupby("Patient")["Response"].agg(lambda x: x.mode().iloc[0])
    a.obs["Response"] = a.obs["Patient"].map(lab.to_dict())
    a = a[a.obs["TimePoint"].astype(str) == "pre"].copy()
    a = a[a.obs["Celltype (major-lineage)"].astype(str).str.contains("CD8", case=False, na=False)].copy()
    a = a[a.obs["Response"].astype(str).isin(["R", "NR"])].copy()
    return a


def score(adata, sig):
    present = [g for g in sig if g in adata.var_names]
    sc.tl.score_genes(adata, present, score_name="immuno_score")
    return present


st.title("ImmunoScore")
st.write(
    "Score single-cell tumor data against a 50-gene CD8 T-cell signature linked to "
    "anti-PD-1 response. The signature was built on melanoma patients (Sade-Feldman 2018) "
    "and tested on basal cell carcinoma (Yost 2019)."
)

sig = load_signature()
with st.expander(f"The {len(sig)} signature genes"):
    st.write(
        "Ranked by differential expression between responders and non-responders in "
        "baseline CD8 T cells, after removing ribosomal, mitochondrial, and stress-gene "
        "artifacts."
    )
    st.write(", ".join(sig))

st.divider()
tab1, tab2 = st.tabs(["Demo", "Upload data"])

with tab1:
    st.write("Score the basal cell carcinoma cohort the signature was never trained on.")
    if st.button("Run demo"):
        with st.spinner("Loading and scoring..."):
            a = load_demo()
            present = score(a, sig)
            pt = a.obs.groupby("Patient").agg(score=("immuno_score", "mean"), resp=("Response", "first"))
            y = (pt["resp"] == "R").astype(int).values
            auc = roc_auc_score(y, pt["score"].values)
        st.metric("Patient-level AUC", f"{auc:.3f}",
                  help="11 patients. Suggestive, not statistically significant.")
        st.write(f"{len(present)} of {len(sig)} genes found. {a.n_obs:,} baseline CD8 T cells across "
                 f"{a.obs['Patient'].nunique()} patients.")
        fig, ax = plt.subplots(figsize=(5, 3))
        for label, color in [("R", "#2e7d6f"), ("NR", "#b5524a")]:
            ax.hist(a.obs.loc[a.obs["Response"] == label, "immuno_score"],
                    bins=30, alpha=0.65, label=label, color=color)
        ax.set_xlabel("ImmunoScore"); ax.set_ylabel("cells")
        ax.legend(title="Response"); ax.spines[["top", "right"]].set_visible(False)
        st.pyplot(fig)
        st.caption("Responders sit slightly higher. With 11 patients this is a hypothesis, not a "
                   "result. The repository reports the confidence interval and a permutation test.")

with tab2:
    up = st.file_uploader("AnnData .h5ad, log-normalized, gene symbols in var_names", type=["h5ad"])
    if up:
        with st.spinner("Scoring..."):
            a = sc.read_h5ad(up)
            present = score(a, sig)
        st.write(f"{a.n_obs:,} cells scored. {len(present)} of {len(sig)} signature genes found.")
        st.dataframe(a.obs[["immuno_score"]].describe())
        fig, ax = plt.subplots(figsize=(5, 3))
        ax.hist(a.obs["immuno_score"], bins=40, color="#3b6ea5")
        ax.set_xlabel("ImmunoScore"); ax.set_ylabel("cells")
        ax.spines[["top", "right"]].set_visible(False)
        st.pyplot(fig)

st.divider()
st.caption("Data: Sade-Feldman et al. 2018; Yost et al. 2019. For research and education, "
           "not clinical use.")
