"""
01_learn_pipeline.py
--------------------
Goal: learn + verify the single-cell RNA-seq workflow on a tiny built-in
dataset (pbmc3k, ~2,700 immune cells) BEFORE touching the real immunotherapy
data. Every step here is identical to what you'll run on the real cohort --
only the data-loading line at the top changes later.

Pipeline: load -> QC -> filter -> normalize -> find variable genes ->
PCA -> neighbors -> UMAP -> cluster (Leiden) -> find marker genes.

Run:  .venv\\Scripts\\python.exe 01_learn_pipeline.py
"""

import scanpy as sc

# ----------------------------------------------------------------------------
# CONFIG  --  your QC tuning knobs (this is a real research judgment call)
# ----------------------------------------------------------------------------
# These thresholds decide which cells are "real" vs junk/dying/doublets.
# Too strict = you throw away good data. Too loose = noise pollutes results.
# Defaults below are the field-standard starting point for pbmc3k.
MIN_GENES_PER_CELL = 200    # drop cells expressing fewer genes (likely empty droplets)
MIN_CELLS_PER_GENE = 3      # drop genes seen in too few cells (noise)
MAX_PCT_MITO       = 5.0    # drop cells with >5% mitochondrial reads (dying cells)
LEIDEN_RESOLUTION  = 1.0    # higher = more, smaller clusters; lower = fewer, broader

sc.settings.verbosity = 1
sc.settings.figdir = "figures"


def load_data():
    """Load the learning dataset. Swap this function out for the real data later."""
    adata = sc.datasets.pbmc3k()          # downloads ~5 MB the first time
    adata.var_names_make_unique()
    print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")
    return adata


def quality_control(adata):
    """Flag and remove low-quality cells/genes using the CONFIG thresholds."""
    # mitochondrial genes start with "MT-" in human data; high % = dying cell
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)

    n_before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=MIN_GENES_PER_CELL)
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)
    adata = adata[adata.obs["pct_counts_mt"] < MAX_PCT_MITO].copy()
    print(f"QC: kept {adata.n_obs}/{n_before} cells")
    return adata


def normalize(adata):
    """Make cells comparable (depth) and stabilize variance, then pick variable genes."""
    sc.pp.normalize_total(adata, target_sum=1e4)   # each cell -> 10k counts
    sc.pp.log1p(adata)                              # log transform
    adata.raw = adata                              # stash full log-data for marker plots
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    return adata


def cluster(adata):
    """Reduce dimensions, build the cell graph, embed (UMAP), and cluster."""
    sc.tl.pca(adata, n_comps=50)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=LEIDEN_RESOLUTION,
                 flavor="igraph", n_iterations=2, directed=False)
    print(f"Found {adata.obs['leiden'].nunique()} clusters")
    return adata


def find_markers(adata):
    """Rank the genes that define each cluster -- this is how you ID cell types."""
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
    # print top 5 marker genes per cluster
    for cl in adata.obs["leiden"].cat.categories:
        top = [adata.uns["rank_genes_groups"]["names"][cl][i] for i in range(5)]
        print(f"  cluster {cl}: {', '.join(top)}")


def main():
    adata = load_data()
    adata = quality_control(adata)
    adata = normalize(adata)
    adata = cluster(adata)
    find_markers(adata)

    # save the two headline figures
    sc.pl.umap(adata, color="leiden", save="_clusters.png", show=False)
    sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False,
                            save="_markers.png", show=False)
    print("\nDone. Figures saved in figures/  ->  umap_clusters.png, ...markers.png")


if __name__ == "__main__":
    main()
