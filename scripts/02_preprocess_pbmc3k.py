#!/usr/bin/env python3
"""Preprocess PBMC3k dataset."""

import scanpy as sc
from sc_cell_state_benchmark.config import DATA_PROCESSED, PBMC3K_PREPROCESSED, PBMC3K_RAW
from sc_cell_state_benchmark.data import ensure_dirs, load_anndata, save_anndata

if __name__ == "__main__":
    # Ensure processed directory exists
    ensure_dirs([DATA_PROCESSED])

    # Load raw data
    adata = load_anndata(PBMC3K_RAW)

    # Filter low-quality cells and genes
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, inplace=True)

    # Normalize total counts
    sc.pp.normalize_total(adata, target_sum=1e4)

    # Log1p transform
    sc.pp.log1p(adata)

    # Identify highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')

    # Subset to HVGs and copy to avoid view warnings
    adata = adata[:, adata.var.highly_variable].copy()

    # Scale data
    sc.pp.scale(adata, max_value=10)

    # PCA
    sc.tl.pca(adata, n_comps=50)

    # Compute neighbors
    sc.pp.neighbors(adata)

    # UMAP
    sc.tl.umap(adata)

    # Leiden clustering
    sc.tl.leiden(adata, flavor="igraph", directed=False, n_iterations=2)

    # Save processed data
    save_anndata(adata, PBMC3K_PREPROCESSED)

    # Print summary
    n_clusters = len(adata.obs['leiden'].unique())
    n_hvg = adata.var.highly_variable.sum()
    print(f"[preprocess] {PBMC3K_RAW.name} -> {PBMC3K_PREPROCESSED.name} | cells={adata.shape[0]} genes={adata.shape[1]} hvg={n_hvg} clusters={n_clusters}")
