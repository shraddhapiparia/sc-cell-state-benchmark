#!/usr/bin/env python3
"""Compute PBMC3k interferon cell-state scores and random control comparisons."""

import numpy as np
import pandas as pd
from sc_cell_state_benchmark.config import PBMC3K_ANNOTATED, RESULTS_TABLES, FIGURES, RANDOM_SEED
from sc_cell_state_benchmark.data import ensure_dirs, load_anndata
from sc_cell_state_benchmark.plotting import (
    plot_random_control_comparison,
    plot_score_violin,
    plot_umap_score,
)
from sc_cell_state_benchmark.scoring import (
    average_expression_score,
    matched_random_gene_sets,
    rank_based_score,
    scanpy_score_genes,
)

INTERFERON_GENES = [
    'IFIT1',
    'IFIT2',
    'IFIT3',
    'ISG15',
    'MX1',
    'OAS1',
    'OASL',
    'RSAD2',
    'IRF7',
    'STAT1',
]

RANDOM_CONTROL_SETS = 25


def flatten_summary(summary_df: pd.DataFrame) -> pd.DataFrame:
    """Flatten a multi-index aggregation summary frame."""
    summary_df.columns = ['_'.join(col).strip() if col[1] else col[0] for col in summary_df.columns.values]
    summary_df = summary_df.reset_index()
    return summary_df


if __name__ == '__main__':
    ensure_dirs([RESULTS_TABLES, FIGURES])

    adata = load_anndata(PBMC3K_ANNOTATED)
    print('[score] loaded annotated AnnData with predicted_cell_type')

    if 'predicted_cell_type' not in adata.obs.columns:
        raise ValueError(
            'predicted_cell_type not found in AnnData. Please run scripts/05_annotate_clusters.py first.'
        )

    available_genes = [g for g in INTERFERON_GENES if g in adata.var_names]
    if len(available_genes) == 0:
        raise RuntimeError('No interferon genes were found in the annotated dataset.')

    # Compute real scores
    scanpy_score_genes(adata, available_genes, score_name='ifn_scanpy')
    adata.obs['ifn_meanexpr'] = average_expression_score(adata, available_genes)
    adata.obs['ifn_rank'] = rank_based_score(adata, available_genes)
    print('[score] computed scanpy score')

    # Generate matched random control sets
    control_sets = matched_random_gene_sets(adata, available_genes, n_sets=RANDOM_CONTROL_SETS, seed=RANDOM_SEED)
    print(f'[score] generated {len(control_sets)} random control gene sets')

    # Compute random control distributions for each method
    random_control_means = {'scanpy': [], 'meanexpr': [], 'rank': []}
    for control_genes in control_sets:
        adata_ctrl = adata.copy()
        scanpy_score_genes(adata_ctrl, control_genes, score_name='ifn_scanpy_ctrl')
        random_control_means['scanpy'].append(adata_ctrl.obs['ifn_scanpy_ctrl'].to_numpy().mean())
        random_control_means['meanexpr'].append(average_expression_score(adata, control_genes).mean())
        random_control_means['rank'].append(rank_based_score(adata, control_genes).mean())

    # Save score tables
    score_table = adata.obs[['leiden', 'predicted_cell_type', 'ifn_scanpy', 'ifn_meanexpr', 'ifn_rank']].copy()
    score_table.insert(0, 'cell_id', adata.obs_names)
    score_path = RESULTS_TABLES / 'interferon_scores.csv'
    score_table.to_csv(score_path, index=False)

    grouping = 'predicted_cell_type' if 'predicted_cell_type' in adata.obs.columns else 'leiden'
    summary = adata.obs.groupby(grouping)[['ifn_scanpy', 'ifn_meanexpr', 'ifn_rank']].agg(['mean', 'median', 'std'])
    summary_df = flatten_summary(summary)
    summary_path = RESULTS_TABLES / 'interferon_summary_by_cell_type.csv'
    summary_df.to_csv(summary_path, index=False)

    # Create figures
    plot_umap_score(adata, 'ifn_scanpy', FIGURES / 'interferon_umap_scanpy.png', 'Interferon score (Scanpy)')
    plot_umap_score(adata, 'ifn_meanexpr', FIGURES / 'interferon_umap_meanexpr.png', 'Interferon score (mean expression)')
    plot_score_violin(
        score_table,
        ['ifn_scanpy', 'ifn_meanexpr', 'ifn_rank'],
        grouping,
        FIGURES / 'interferon_violin_by_cell_type.png',
    )
    plot_random_control_comparison(
        adata.obs['ifn_scanpy'].to_numpy(),
        np.array(random_control_means['scanpy']),
        FIGURES / 'interferon_vs_random_controls.png',
    )

    print(f'[score] saved figures to {FIGURES}/')
    print(f'[score] saved score tables to {RESULTS_TABLES}/')
    print('[score] completed interferon scoring benchmark')
