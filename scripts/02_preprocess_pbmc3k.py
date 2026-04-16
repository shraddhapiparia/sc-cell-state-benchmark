#!/usr/bin/env python3
"""Preprocess PBMC3k dataset.

QC thresholds below were chosen by inspecting figures/pbmc3k_qc_violin_prefilter.png.
Adjust MAX_GENES and MAX_PCT_MT if the pre-filter violin shows a different distribution
in your run.

  MAX_GENES  = 2500   -- PBMC3k has only 5 cells above this (likely doublets or empty
                          drops); 99th percentile is ~1750.
  MAX_PCT_MT = 5      -- median is 2%, 95th percentile is 4%; 57 cells exceed 5%.
                          The single outlier at 22.6% is visible as a clear tail.
"""

import scanpy as sc

from sc_cell_state_benchmark.config import (
    DATA_PROCESSED,
    FIGURES,
    PBMC3K_PREPROCESSED,
    PBMC3K_RAW,
    RANDOM_SEED,
)
from sc_cell_state_benchmark.data import ensure_dirs, load_anndata, save_anndata
from sc_cell_state_benchmark.plotting import plot_qc_violin

# QC filter thresholds -- inspect pbmc3k_qc_violin_prefilter.png before changing
MAX_GENES: int = 2500
MAX_PCT_MT: float = 5.0

if __name__ == "__main__":
    ensure_dirs([DATA_PROCESSED, FIGURES])

    adata = load_anndata(PBMC3K_RAW)
    print(f"[preprocess] loaded {PBMC3K_RAW.name} | cells={adata.n_obs} genes={adata.n_vars}")

    # --- mitochondrial gene annotation ----------------------------------
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    n_mt = int(adata.var["mt"].sum())
    print(f"[preprocess] {n_mt} mitochondrial genes annotated")

    # --- basic cell and gene filters ------------------------------------
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"[preprocess] after min-gene/min-cell filter: cells={adata.n_obs}")

    # --- QC metrics -----------------------------------------------------
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    # --- pre-filter QC violin -------------------------------------------
    # Inspect this figure to verify the threshold choices above are appropriate.
    qc_fig_path = FIGURES / "pbmc3k_qc_violin_prefilter.png"
    plot_qc_violin(
        adata,
        qc_fig_path,
        thresholds={"n_genes_by_counts": MAX_GENES, "pct_counts_mt": MAX_PCT_MT},
    )
    print(f"[preprocess] saved pre-filter QC figure to {qc_fig_path}")

    # --- upper-bound and MT filters -------------------------------------
    n_before = adata.n_obs
    sc.pp.filter_cells(adata, max_genes=MAX_GENES)
    n_after_genes = adata.n_obs
    adata = adata[adata.obs["pct_counts_mt"] < MAX_PCT_MT].copy()
    n_after_mt = adata.n_obs
    print(
        f"[preprocess] upper-bound filter (>{MAX_GENES} genes): removed {n_before - n_after_genes} cells"
    )
    print(
        f"[preprocess] MT filter (>={MAX_PCT_MT}% MT): removed {n_after_genes - n_after_mt} cells"
    )
    print(f"[preprocess] cells retained: {adata.n_obs} / {n_before}")

    # --- normalisation and log transform --------------------------------
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Store the full log-normalised matrix before HVG subsetting.
    # Scoring functions (score_cell_states.py) use adata.raw to access
    # the complete gene vocabulary.
    adata.raw = adata.copy()

    # --- highly variable genes ------------------------------------------
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
    n_hvg = int(adata.var["highly_variable"].sum())
    adata = adata[:, adata.var.highly_variable].copy()
    print(f"[preprocess] {n_hvg} HVGs selected; matrix reduced to {adata.n_vars} genes")

    # --- dimensionality reduction and clustering ------------------------
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, random_state=RANDOM_SEED)
    sc.pp.neighbors(adata, random_state=RANDOM_SEED)
    sc.tl.umap(adata, random_state=RANDOM_SEED)
    sc.tl.leiden(
        adata, flavor="igraph", directed=False, n_iterations=2, random_state=RANDOM_SEED
    )

    save_anndata(adata, PBMC3K_PREPROCESSED)

    n_clusters = adata.obs["leiden"].nunique()
    print(
        f"[preprocess] {PBMC3K_RAW.name} -> {PBMC3K_PREPROCESSED.name} | "
        f"cells={adata.n_obs} genes={adata.n_vars} hvg={n_hvg} clusters={n_clusters}"
    )
