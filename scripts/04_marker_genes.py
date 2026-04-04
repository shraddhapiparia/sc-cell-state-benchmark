#!/usr/bin/env python3
"""Extract marker genes for PBMC3k clusters."""

import pandas as pd
import scanpy as sc
from pathlib import Path
from sc_cell_state_benchmark.config import PBMC3K_PREPROCESSED, RESULTS_TABLES
from sc_cell_state_benchmark.data import load_anndata

if __name__ == "__main__":
    # Load processed data
    adata = load_anndata(PBMC3K_PREPROCESSED)

    # Run differential expression
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

    # Extract top 10 markers per cluster
    markers = []
    for cluster in adata.obs['leiden'].cat.categories:
        cluster_markers = sc.get.rank_genes_groups_df(adata, group=cluster).head(10)
        cluster_markers['cluster'] = cluster
        markers.append(cluster_markers)

    # Combine and save
    markers_df = pd.concat(markers)
    markers_df = markers_df[['cluster', 'names', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj']]
    output_path = RESULTS_TABLES / "pbmc3k_top_markers.csv"
    markers_df.to_csv(output_path, index=False)

    # Print summary
    n_clusters = len(adata.obs['leiden'].cat.categories)
    print(f"[markers] saved top marker genes for {n_clusters} clusters to {output_path}")