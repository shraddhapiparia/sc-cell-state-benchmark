#!/usr/bin/env python3
"""Within-cell-type differential expression: stim vs ctrl (Kang PBMC).

For each cell type the dataset is subsetted and Wilcoxon rank-sum DE is run
comparing stimulated cells against control cells.  Results are combined into a
single table and a summary heatmap is produced.

Prerequisite: scripts/08_preprocess_kang_pbmc.py must have been run.

Outputs
-------
  results/tables/kang_de_stim_vs_ctrl_per_cell_type.csv
      Columns: cell_type, names, scores, logfoldchanges, pvals, pvals_adj,
               n_stim, n_ctrl

  figures/kang_de_top5_per_cell_type.png
      Heatmap of log-fold-changes for the top-5 DE genes per cell type.

  figures/kang_de_volcano_cd14_monocytes.png
  figures/kang_de_volcano_cd4_t_cells.png
      Volcano plots for the two most interpretable cell types.
"""

import argparse
import re

import pandas as pd
import scanpy as sc

from sc_cell_state_benchmark.config import (
    FIGURES,
    KANG_PBMC_PREPROCESSED,
    RESULTS_TABLES,
)
from sc_cell_state_benchmark.data import ensure_dirs, load_anndata
from sc_cell_state_benchmark.plotting import plot_de_summary, plot_volcano

# Minimum cells required per condition group to attempt DE for a cell type.
MIN_CELLS_PER_GROUP = 30

# Cell types for which individual volcano plots are produced.
VOLCANO_CELL_TYPES = ['CD14+ Monocytes', 'CD4 T cells']

CONDITION_KEYWORDS = ['stim', 'condition', 'treat', 'group', 'status', 'label']
CELL_TYPE_KEYWORDS = ['cell_type', 'celltype', 'celllabels', 'cell_labels', 'annotation', 'subclass']


def parse_args():
    parser = argparse.ArgumentParser(
        description='Within-cell-type DE: stim vs ctrl for the Kang PBMC dataset.'
    )
    parser.add_argument('--input', type=str, default=None,
                        help='Path to preprocessed AnnData. Defaults to standard path.')
    parser.add_argument('--condition', type=str, default=None,
                        help='Obs column for condition labels.')
    parser.add_argument('--cell-type', type=str, default=None,
                        help='Obs column for cell type labels.')
    return parser.parse_args()


def find_metadata_column(obs, keywords):
    for col in obs.columns:
        if any(kw in col.lower() for kw in keywords):
            return col
    return None


def safe_filename(name: str) -> str:
    """Convert a cell-type label to a safe filename fragment."""
    return re.sub(r'[^a-zA-Z0-9]+', '_', name).strip('_').lower()


if __name__ == '__main__':
    args = parse_args()
    ensure_dirs([RESULTS_TABLES, FIGURES])

    input_path = args.input or KANG_PBMC_PREPROCESSED
    adata = load_anndata(input_path)
    print(f'[de] loaded dataset from {input_path}')
    print(f'[de] {adata.n_obs} cells | {adata.n_vars} genes (HVG subset)')
    if adata.raw is not None:
        print(f'[de] {adata.raw.n_vars} genes in adata.raw (used for scoring context)')

    condition_col = args.condition or find_metadata_column(adata.obs, CONDITION_KEYWORDS)
    if condition_col is None or condition_col not in adata.obs.columns:
        raise ValueError('Could not find condition column. Provide --condition <col>.')

    cell_type_col = args.cell_type or find_metadata_column(adata.obs, CELL_TYPE_KEYWORDS)
    if cell_type_col is None or cell_type_col not in adata.obs.columns:
        raise ValueError('Could not find cell type column. Provide --cell-type <col>.')

    print(f'[de] condition column: {condition_col!r}')
    print(f'[de] cell type column: {cell_type_col!r}')

    unique_conditions = list(adata.obs[condition_col].astype(str).unique())
    if 'stim' not in unique_conditions or 'ctrl' not in unique_conditions:
        raise ValueError(
            f'Expected condition values "ctrl" and "stim". Found: {unique_conditions}. '
            'Pass --condition with the correct column name.'
        )

    cell_types = sorted(adata.obs[cell_type_col].astype(str).unique())
    print(f'[de] {len(cell_types)} cell types: {cell_types}\n')

    all_results = []

    for ct in cell_types:
        mask = adata.obs[cell_type_col].astype(str) == ct
        subset = adata[mask].copy()

        n_stim = int((subset.obs[condition_col] == 'stim').sum())
        n_ctrl = int((subset.obs[condition_col] == 'ctrl').sum())

        if n_stim < MIN_CELLS_PER_GROUP or n_ctrl < MIN_CELLS_PER_GROUP:
            print(
                f'[de] SKIP {ct!r}: n_stim={n_stim}, n_ctrl={n_ctrl} '
                f'(minimum {MIN_CELLS_PER_GROUP} per group required)'
            )
            continue

        sc.tl.rank_genes_groups(
            subset,
            groupby=condition_col,
            groups=['stim'],
            reference='ctrl',
            method='wilcoxon',
            use_raw=True,
        )

        result = sc.get.rank_genes_groups_df(subset, group='stim')
        result['cell_type'] = ct
        result['n_stim'] = n_stim
        result['n_ctrl'] = n_ctrl
        all_results.append(result)

        top3 = result.nlargest(3, 'scores')['names'].tolist()
        print(f'[de] {ct!r}: n_stim={n_stim} n_ctrl={n_ctrl} | top genes: {top3}')

    if not all_results:
        raise RuntimeError('No cell types passed the minimum cell count filter.')

    de_df = pd.concat(all_results, ignore_index=True)
    de_df = de_df[['cell_type', 'names', 'scores', 'logfoldchanges',
                   'pvals', 'pvals_adj', 'n_stim', 'n_ctrl']]

    out_table = RESULTS_TABLES / 'kang_de_stim_vs_ctrl_per_cell_type.csv'
    de_df.to_csv(out_table, index=False)
    print(f'\n[de] saved combined DE table: {out_table}')
    print(f'[de] {len(de_df)} rows across {de_df["cell_type"].nunique()} cell types')

    summary_path = FIGURES / 'kang_de_top5_per_cell_type.png'
    plot_de_summary(de_df, summary_path, top_n=5)
    print(f'[de] saved summary heatmap: {summary_path}')

    for ct in VOLCANO_CELL_TYPES:
        if ct not in de_df['cell_type'].values:
            print(f'[de] volcano skipped for {ct!r}: not in results')
            continue
        fname = f'kang_de_volcano_{safe_filename(ct)}.png'
        vpath = FIGURES / fname
        plot_volcano(de_df, cell_type=ct, save_path=vpath)
        print(f'[de] saved volcano: {vpath}')

    print('\n[de] completed within-cell-type DE analysis')
