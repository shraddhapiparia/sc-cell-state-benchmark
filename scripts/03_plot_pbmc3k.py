#!/usr/bin/env python3
"""Plot PBMC3k UMAP visualizations."""

from pathlib import Path
from sc_cell_state_benchmark.config import FIGURES, PBMC3K_PREPROCESSED
from sc_cell_state_benchmark.data import load_anndata
from sc_cell_state_benchmark.plotting import plot_umap_clusters, plot_umap_qc

if __name__ == "__main__":
    # Load processed data
    adata = load_anndata(PBMC3K_PREPROCESSED)

    # Create and save cluster UMAP
    cluster_path = FIGURES / "pbmc3k_umap_clusters.png"
    plot_umap_clusters(adata, cluster_path)

    # Create and save QC UMAP
    qc_path = FIGURES / "pbmc3k_umap_qc.png"
    plot_umap_qc(adata, qc_path)

    # Print summary
    print(f"[plot] saved UMAP figures to {FIGURES}/")