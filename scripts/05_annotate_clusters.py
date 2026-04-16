#!/usr/bin/env python3
"""Annotate PBMC3k Leiden clusters with likely cell types."""

import pandas as pd
from pathlib import Path
from sc_cell_state_benchmark.config import RESULTS_TABLES, FIGURES, PBMC3K_PREPROCESSED, DATA_PROCESSED, PBMC3K_ANNOTATED
from sc_cell_state_benchmark.data import load_anndata, ensure_dirs, save_anndata
from sc_cell_state_benchmark.plotting import plot_marker_heatmap, plot_umap_annotation

# T cells are split into CD4 and CD8 subsets.
# CD3D is the shared T-cell anchor; IL7R enriches for CD4, CD8A/NKG7 for CD8.
# MALAT1 removed -- it is a broadly expressed lncRNA with no T-cell specificity.
CELL_TYPE_RULES = {
    'CD4 T cell': {'CD3D', 'IL7R', 'CD4', 'LTB'},
    'CD8 T cell': {'CD3D', 'CD8A', 'NKG7', 'CST7'},
    'B cell':     {'MS4A1', 'CD79A', 'CD74', 'HLA-DRA'},
    'NK cell':    {'NKG7', 'GNLY', 'KLRD1'},
    'Monocyte':   {'LYZ', 'CST3', 'FCN1', 'S100A8', 'FCGR3A'},
    'Dendritic cell': {'FCER1A', 'CST3', 'HLA-DRA'},
    'Platelet':   {'PPBP', 'PF4'},
}

MARKER_GENES = [
    'CD3D', 'IL7R', 'CD4', 'CD8A',
    'NKG7', 'MS4A1', 'LST1', 'FCGR3A', 'HLA-DRA',
]


def predict_cell_type(marker_genes):
    """Predict the cell type for a cluster using simple marker rules."""
    scores = {}
    for cell_type, markers in CELL_TYPE_RULES.items():
        scores[cell_type] = len(set(marker_genes) & markers)

    best_type = max(scores, key=scores.get)
    best_score = scores[best_type]
    counts = list(scores.values())

    if best_score == 0 or counts.count(best_score) > 1:
        return 'Unknown'
    return best_type


def annotate_clusters(marker_table: pd.DataFrame) -> pd.DataFrame:
    """Create cluster annotations from top marker genes."""
    records = []
    for cluster, group in marker_table.groupby('cluster'):
        top_genes = group['names'].head(10).tolist()
        predicted = predict_cell_type(top_genes)
        records.append(
            {
                'cluster': cluster,
                'top_marker_genes': ','.join(top_genes),
                'predicted_cell_type': predicted,
            }
        )
    return pd.DataFrame(records)


if __name__ == '__main__':
    ensure_dirs([RESULTS_TABLES, FIGURES, DATA_PROCESSED])

    marker_path = RESULTS_TABLES / 'pbmc3k_top_markers.csv'
    markers = pd.read_csv(marker_path)
    annotations = annotate_clusters(markers)

    output_table = RESULTS_TABLES / 'pbmc3k_cluster_annotations.csv'
    annotations.to_csv(output_table, index=False)

    adata = load_anndata(PBMC3K_PREPROCESSED)
    cluster_map = annotations.assign(cluster=annotations['cluster'].astype(str))
    cluster_to_type = cluster_map.set_index('cluster')['predicted_cell_type'].to_dict()

    target_clusters = adata.obs['leiden'].astype(str)
    predicted = target_clusters.map(cluster_to_type).fillna('Unknown')
    predicted = pd.Categorical(predicted, categories=sorted(predicted.unique()))

    adata.obs['predicted_cell_type'] = predicted

    save_anndata(adata, PBMC3K_ANNOTATED)
    print(f"[annotate] saved annotated AnnData to {PBMC3K_ANNOTATED}")

    print(f"[annotate-debug] applied mapping: {cluster_to_type}")
    print(
        f"[annotate-debug] unique predicted_cell_type values: {', '.join(map(str, adata.obs['predicted_cell_type'].cat.categories))}"
    )

    annotated_path = FIGURES / 'pbmc3k_umap_annotated.png'
    plot_umap_annotation(adata, annotated_path)

    marker_path = FIGURES / 'pbmc3k_marker_heatmap.png'
    plot_marker_heatmap(
        adata,
        MARKER_GENES,
        groupby='predicted_cell_type',
        save_path=marker_path,
        title='Canonical marker expression by predicted cell type (PBMC3k)',
    )
    print(f"[annotate] saved marker heatmap to {marker_path}")

    n_clusters = len(annotations)
    print(f"[annotate] assigned likely cell types to {n_clusters} clusters")