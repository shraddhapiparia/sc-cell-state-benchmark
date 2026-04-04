#!/usr/bin/env python3
"""Score a perturbation single-cell dataset and compare stimulated vs control."""

import argparse
import numpy as np
import pandas as pd
from sc_cell_state_benchmark.config import (
    KANG_PBMC_PREPROCESSED,
    RESULTS_TABLES,
    FIGURES,
    RANDOM_SEED,
)
from sc_cell_state_benchmark.data import ensure_dirs, load_anndata
from sc_cell_state_benchmark.evaluation import compute_binary_comparison
from sc_cell_state_benchmark.plotting import (
    plot_random_control_comparison,
    plot_score_by_cell_type_and_condition,
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
CONDITION_KEYWORDS = ['stim', 'condition', 'treat', 'group', 'status']
CELL_TYPE_KEYWORDS = ['cell_type', 'celltype', 'celllabels', 'cell_labels', 'annotation', 'subclass']


def parse_args():
    parser = argparse.ArgumentParser(description='Score a perturbation single-cell dataset.')
    parser.add_argument(
        '--input',
        type=str,
        default=None,
        help='Path to a preprocessed AnnData file to score. Defaults to the standard preprocessed path.',
    )
    parser.add_argument(
        '--condition',
        type=str,
        default=None,
        help='Name of the condition column in adata.obs for stimulated/control comparison.',
    )
    parser.add_argument(
        '--cell-type',
        type=str,
        default=None,
        help='Optional cell type column name in adata.obs for additional summaries.',
    )
    return parser.parse_args()


def find_metadata_column(obs, keywords):
    for col in obs.columns:
        low = col.lower()
        if any(keyword in low for keyword in keywords):
            return col
    return None


def infer_positive_label(values):
    labels = [str(v).lower() for v in values]
    if 'stim' in ' '.join(labels):
        for v in values:
            if 'stim' in str(v).lower():
                return v
    if 'control' in ' '.join(labels) or 'ctrl' in ' '.join(labels):
        for v in values:
            if 'control' not in str(v).lower() and 'ctrl' not in str(v).lower():
                return v
    if len(values) == 2:
        return values[1]
    raise ValueError('Could not infer stimulated label from condition values: {}'.format(values))


def flatten_summary(summary_df: pd.DataFrame) -> pd.DataFrame:
    summary_df.columns = ['_'.join(col).strip() if col[1] else col[0] for col in summary_df.columns.values]
    return summary_df.reset_index()


if __name__ == '__main__':
    args = parse_args()
    ensure_dirs([RESULTS_TABLES, FIGURES])

    input_path = args.input if args.input is not None else KANG_PBMC_PREPROCESSED
    adata = load_anndata(input_path)
    print(f'[score] loaded preprocessed dataset from {input_path}')
    print(f'[score] available obs columns: {list(adata.obs.columns)}')

    condition_col = args.condition
    if condition_col is None:
        condition_col = find_metadata_column(adata.obs, CONDITION_KEYWORDS)
        if condition_col is not None:
            print(f'[score] auto-detected condition column: {condition_col}')
    else:
        if condition_col not in adata.obs.columns:
            raise ValueError(f'Condition column {condition_col} not found in adata.obs')

    if condition_col is None:
        raise ValueError(
            'Could not find a stimulation condition metadata column. '
            'Provide --condition <column_name> for the dataset.'
        )

    if args.cell_type is not None:
        cell_type_col = args.cell_type
        if cell_type_col not in adata.obs.columns:
            raise ValueError(f'Cell type column {cell_type_col} not found in adata.obs')
    else:
        cell_type_col = find_metadata_column(adata.obs, CELL_TYPE_KEYWORDS)
        if cell_type_col is not None:
            print(f'[score] auto-detected cell type column: {cell_type_col}')

    condition_values = adata.obs[condition_col].astype(str)
    unique_conditions = condition_values.unique()
    print(f'[score] using condition column: {condition_col} with values {list(unique_conditions)}')

    if len(unique_conditions) != 2:
        raise ValueError(
            f'Expected exactly 2 stimulation groups for comparison, found {len(unique_conditions)}: {list(unique_conditions)}'
        )

    positive_label = infer_positive_label(unique_conditions)
    print(f'[score] inferred stimulated label: {positive_label}')

    available_genes = [g for g in INTERFERON_GENES if g in adata.var_names]
    if len(available_genes) == 0:
        raise RuntimeError('No interferon genes were found in the provided dataset.')

    scanpy_score_genes(adata, available_genes, score_name='ifn_scanpy')
    adata.obs['ifn_meanexpr'] = average_expression_score(adata, available_genes)
    adata.obs['ifn_rank'] = rank_based_score(adata, available_genes)
    print('[score] computed interferon scores')

    control_sets = matched_random_gene_sets(adata, available_genes, n_sets=RANDOM_CONTROL_SETS, seed=RANDOM_SEED)
    print(f'[score] generated {len(control_sets)} random control gene sets')

    random_control_means = {'scanpy': [], 'meanexpr': [], 'rank': []}
    for control_genes in control_sets:
        adata_ctrl = adata.copy()
        scanpy_score_genes(adata_ctrl, control_genes, score_name='ifn_scanpy_ctrl')
        random_control_means['scanpy'].append(adata_ctrl.obs['ifn_scanpy_ctrl'].to_numpy().mean())
        random_control_means['meanexpr'].append(average_expression_score(adata, control_genes).mean())
        random_control_means['rank'].append(rank_based_score(adata, control_genes).mean())

    score_table = adata.obs[[condition_col, 'ifn_scanpy', 'ifn_meanexpr', 'ifn_rank']].copy()
    score_table.insert(0, 'cell_id', adata.obs_names)
    score_table.rename(columns={condition_col: 'condition'}, inplace=True)
    score_path = RESULTS_TABLES / 'kang_interferon_scores.csv'
    score_table.to_csv(score_path, index=False)

    summary_cond = adata.obs.groupby(condition_col, observed=True)[['ifn_scanpy', 'ifn_meanexpr', 'ifn_rank']].agg(['mean', 'median', 'std'])
    summary_cond_df = flatten_summary(summary_cond)
    summary_cond_path = RESULTS_TABLES / 'kang_interferon_summary_by_condition.csv'
    summary_cond_df.to_csv(summary_cond_path, index=False)

    if cell_type_col is not None:
        summary_cell_type = adata.obs.groupby(cell_type_col, observed=True)[['ifn_scanpy', 'ifn_meanexpr', 'ifn_rank']].agg(['mean', 'median', 'std'])
        summary_cell_type_df = flatten_summary(summary_cell_type)
        summary_cell_type_path = RESULTS_TABLES / 'kang_interferon_summary_by_cell_type.csv'
        summary_cell_type_df.to_csv(summary_cell_type_path, index=False)

        # Summary table stratified by cell type AND condition
        summary_by_ct_cond = (
            adata.obs.groupby([cell_type_col, condition_col], observed=True)
            .agg(n_cells=('ifn_scanpy', 'count'), mean=('ifn_scanpy', 'mean'), median=('ifn_scanpy', 'median'), std=('ifn_scanpy', 'std'))
            .reset_index()
        )
        # Also add meanexpr stats
        meanexpr_stats = (
            adata.obs.groupby([cell_type_col, condition_col], observed=True)['ifn_meanexpr']
            .agg(meanexpr_mean='mean', meanexpr_median='median', meanexpr_std='std')
            .reset_index()
        )
        summary_by_ct_cond = summary_by_ct_cond.merge(meanexpr_stats, on=[cell_type_col, condition_col])
        summary_ct_cond_path = RESULTS_TABLES / 'kang_interferon_summary_by_cell_type_and_condition.csv'
        summary_by_ct_cond.to_csv(summary_ct_cond_path, index=False)
        print(f'[score] saved per-cell-type × condition summary to {summary_ct_cond_path}')

        # Console summary: top 3 cell types by stimulated score and by ctrl-vs-stim delta
        stim_means = (
            summary_by_ct_cond[summary_by_ct_cond[condition_col] == positive_label]
            .set_index(cell_type_col)['mean']
            .sort_values(ascending=False)
        )
        ctrl_means = (
            summary_by_ct_cond[summary_by_ct_cond[condition_col] != positive_label]
            .set_index(cell_type_col)['mean']
        )
        delta = (stim_means - ctrl_means).dropna().sort_values(ascending=False)

        print('\n[score] Top 3 cell types by stimulated interferon score (scanpy):')
        for ct, val in stim_means.head(3).items():
            print(f'  {ct}: {val:.4f}')
        print('[score] Top 3 cell types by ctrl-vs-stim delta (scanpy):')
        for ct, val in delta.head(3).items():
            print(f'  {ct}: {val:+.4f}')
        print()
    else:
        summary_cell_type_path = None

    method_comparisons = []
    for method in ['ifn_scanpy', 'ifn_meanexpr', 'ifn_rank']:
        metrics = compute_binary_comparison(adata.obs[method], condition_values, positive_label)
        metrics.update({'method': method, 'condition_column': condition_col})
        method_comparisons.append(metrics)
    comparison_df = pd.DataFrame(method_comparisons)
    comparison_path = RESULTS_TABLES / 'kang_method_comparison.csv'
    comparison_df.to_csv(comparison_path, index=False)

    plot_umap_score(adata, 'ifn_scanpy', FIGURES / 'kang_interferon_umap_scanpy.png', 'Interferon score (Scanpy)')
    plot_umap_score(adata, 'ifn_meanexpr', FIGURES / 'kang_interferon_umap_meanexpr.png', 'Interferon score (mean expression)')
    plot_score_violin(
        score_table,
        ['ifn_scanpy', 'ifn_meanexpr', 'ifn_rank'],
        'condition',
        FIGURES / 'kang_interferon_violin_by_condition.png',
    )
    plot_random_control_comparison(
        adata.obs['ifn_scanpy'].to_numpy(),
        np.array(random_control_means['scanpy']),
        FIGURES / 'kang_interferon_vs_random_controls.png',
    )

    if cell_type_col is not None:
        plot_score_by_cell_type_and_condition(
            adata.obs,
            score_field='ifn_scanpy',
            cell_type_col=cell_type_col,
            condition_col=condition_col,
            save_path=FIGURES / 'kang_interferon_by_cell_type_and_condition.png',
            title='Interferon score (Scanpy) by cell type and condition',
        )

    print(f'[score] saved figures to {FIGURES}/')
    print(f'[score] saved score tables to {RESULTS_TABLES}/')
    print('[score] completed perturbation interferon benchmark')
