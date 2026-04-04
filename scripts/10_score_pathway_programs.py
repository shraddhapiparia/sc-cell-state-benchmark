#!/usr/bin/env python3
"""Score curated immune biological programs and summarise by cell type and condition.

Reads the preprocessed perturbation dataset produced by script 08, scores five
curated immune programs using Scanpy score_genes, and produces:

  Per-cell output
    results/tables/kang_program_scores.csv

  Per-cell-type x condition summary
    results/tables/kang_program_summary_by_cell_type_and_condition.csv

  Program AUROC ranking
    results/tables/kang_program_auc.csv

  Option B (primary figure): delta heatmap (stim - ctrl mean per program x cell type)
    figures/kang_program_heatmap_delta.png

  Option A (supporting figure): side-by-side heatmaps (ctrl | stim raw scores)
    figures/kang_program_heatmap_sidebyside.png

  AUROC ranking bar chart
    figures/kang_program_auc_barplot.png

Prerequisite: scripts/08_preprocess_kang_pbmc.py must have been run first.
"""

import argparse
import numpy as np
import pandas as pd

from sc_cell_state_benchmark.config import (
    FIGURES,
    KANG_PBMC_PREPROCESSED,
    RESULTS_TABLES,
)
from sc_cell_state_benchmark.data import ensure_dirs, load_anndata
from sc_cell_state_benchmark.evaluation import compute_binary_comparison
from sc_cell_state_benchmark.gene_sets import (
    IMMUNE_PROGRAMS,
    MIN_GENES_REQUIRED,
    PROGRAM_DISPLAY_NAMES,
)
from sc_cell_state_benchmark.plotting import (
    plot_program_auc_barplot,
    plot_program_heatmap_delta,
    plot_program_heatmap_sidebyside,
)
from sc_cell_state_benchmark.scoring import scanpy_score_genes

# Keyword lists match script 09 for consistent auto-detection behaviour.
CONDITION_KEYWORDS = ['stim', 'condition', 'treat', 'group', 'status']
CELL_TYPE_KEYWORDS = ['cell_type', 'celltype', 'celllabels', 'cell_labels', 'annotation', 'subclass']


def parse_args():
    parser = argparse.ArgumentParser(
        description='Score curated immune programs and summarise by cell type and condition.'
    )
    parser.add_argument(
        '--input',
        type=str,
        default=None,
        help='Path to preprocessed AnnData (.h5ad). Defaults to the standard preprocessed path.',
    )
    parser.add_argument(
        '--condition',
        type=str,
        default=None,
        help='Obs column for condition (stimulated/control). Auto-detected if not provided.',
    )
    parser.add_argument(
        '--cell-type',
        type=str,
        default=None,
        help='Obs column for cell type labels. Auto-detected if not provided.',
    )
    return parser.parse_args()


def find_metadata_column(obs, keywords):
    """Return the first obs column name that contains any of the given keywords.

    Matches the logic in script 09 for consistent auto-detection across scripts.
    """
    for col in obs.columns:
        if any(keyword in col.lower() for keyword in keywords):
            return col
    return None


def infer_positive_label(values):
    """Infer the stimulated condition label from the set of unique condition values.

    Matches the logic in script 09.
    """
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
    raise ValueError('Could not infer stimulated label from condition values: {}'.format(list(values)))


if __name__ == '__main__':
    args = parse_args()
    ensure_dirs([RESULTS_TABLES, FIGURES])

    # ------------------------------------------------------------------
    # Load data
    # ------------------------------------------------------------------
    input_path = args.input if args.input is not None else KANG_PBMC_PREPROCESSED
    adata = load_anndata(input_path)
    print(f'[programs] loaded dataset from {input_path}')
    print(f'[programs] {adata.n_obs} cells | {adata.n_vars} genes (HVG subset)')
    if adata.raw is not None:
        print(f'[programs] {adata.raw.n_vars} genes in adata.raw (full matrix used for scoring)')

    # ------------------------------------------------------------------
    # Resolve metadata columns
    # ------------------------------------------------------------------
    condition_col = args.condition
    if condition_col is None:
        condition_col = find_metadata_column(adata.obs, CONDITION_KEYWORDS)
        if condition_col is not None:
            print(f'[programs] auto-detected condition column: {condition_col!r}')
    if condition_col is None or condition_col not in adata.obs.columns:
        raise ValueError(
            'Could not find a condition column. Provide --condition <column_name>.'
        )

    cell_type_col = args.cell_type
    if cell_type_col is None:
        cell_type_col = find_metadata_column(adata.obs, CELL_TYPE_KEYWORDS)
        if cell_type_col is not None:
            print(f'[programs] auto-detected cell type column: {cell_type_col!r}')
    if cell_type_col is not None and cell_type_col not in adata.obs.columns:
        raise ValueError(f'Cell type column {cell_type_col!r} not found in adata.obs.')

    condition_values = adata.obs[condition_col].astype(str)
    unique_conditions = condition_values.unique()
    if len(unique_conditions) != 2:
        raise ValueError(
            f'Expected exactly 2 condition groups, found {len(unique_conditions)}: {list(unique_conditions)}'
        )
    positive_label = infer_positive_label(unique_conditions)
    print(f'[programs] condition: {condition_col!r} | stimulated label: {positive_label!r}')

    if cell_type_col is not None:
        n_types = adata.obs[cell_type_col].nunique()
        print(f'[programs] cell type: {cell_type_col!r} | {n_types} types: '
              f'{sorted(adata.obs[cell_type_col].astype(str).unique())}')

    # ------------------------------------------------------------------
    # Gene availability report
    # ------------------------------------------------------------------
    gene_vocab = set(adata.raw.var_names if adata.raw is not None else adata.var_names)
    print('\n[programs] gene availability per program:')
    scored_programs = []
    gene_availability = {}
    for prog_name, gene_list in IMMUNE_PROGRAMS.items():
        found = [g for g in gene_list if g in gene_vocab]
        missing = [g for g in gene_list if g not in gene_vocab]
        gene_availability[prog_name] = len(found)
        status = f'{len(found):2d}/{len(gene_list)} found'
        if missing:
            status += f'  [not in dataset: {", ".join(missing)}]'
        print(f'  {prog_name:<35} {status}')
        if len(found) >= MIN_GENES_REQUIRED:
            scored_programs.append(prog_name)
        else:
            print(f'  *** skipping {prog_name}: fewer than {MIN_GENES_REQUIRED} genes found ***')

    if not scored_programs:
        raise RuntimeError(
            'No programs had enough genes to score. '
            'Verify the dataset uses human HGNC gene symbols.'
        )
    print()

    # ------------------------------------------------------------------
    # Score all programs using scanpy score_genes
    # ------------------------------------------------------------------
    score_cols = []
    for prog_name in scored_programs:
        gene_list = [g for g in IMMUNE_PROGRAMS[prog_name] if g in gene_vocab]
        score_col = f'prog_{prog_name}'
        scanpy_score_genes(adata, gene_list, score_name=score_col)
        score_cols.append(score_col)
    print(f'[programs] scored {len(scored_programs)} programs: {scored_programs}')

    # ------------------------------------------------------------------
    # Per-cell score table
    # ------------------------------------------------------------------
    obs_cols = [condition_col] + ([] if cell_type_col is None else [cell_type_col]) + score_cols
    cell_table = adata.obs[obs_cols].copy()
    cell_table.insert(0, 'cell_id', adata.obs_names)
    cell_table.rename(columns={condition_col: 'condition'}, inplace=True)
    if cell_type_col is not None:
        cell_table.rename(columns={cell_type_col: 'cell_type'}, inplace=True)
    cell_table_path = RESULTS_TABLES / 'kang_program_scores.csv'
    cell_table.to_csv(cell_table_path, index=False)
    print(f'[programs] saved per-cell score table: {cell_table_path}')

    # ------------------------------------------------------------------
    # Per-cell-type x condition summary table
    # ------------------------------------------------------------------
    if cell_type_col is not None:
        summary_rows = []
        for prog_name in scored_programs:
            score_col = f'prog_{prog_name}'
            grp = (
                adata.obs
                .groupby([cell_type_col, condition_col], observed=True)[score_col]
                .agg(n_cells='count', mean='mean', median='median', std='std')
                .reset_index()
                .rename(columns={cell_type_col: 'cell_type', condition_col: 'condition'})
            )
            grp['program'] = prog_name
            grp['display_name'] = PROGRAM_DISPLAY_NAMES[prog_name]
            summary_rows.append(grp)
        summary_df = pd.concat(summary_rows, ignore_index=True)
        summary_df = summary_df[['program', 'display_name', 'cell_type', 'condition',
                                  'n_cells', 'mean', 'median', 'std']]
        summary_path = RESULTS_TABLES / 'kang_program_summary_by_cell_type_and_condition.csv'
        summary_df.to_csv(summary_path, index=False)
        print(f'[programs] saved per-cell-type x condition summary: {summary_path}')
    else:
        summary_df = None
        print('[programs] no cell type column — skipping per-cell-type summary and heatmaps')

    # ------------------------------------------------------------------
    # AUROC per program (global ctrl vs stim, all cell types pooled)
    # ------------------------------------------------------------------
    auc_rows = []
    for prog_name in scored_programs:
        score_col = f'prog_{prog_name}'
        try:
            metrics = compute_binary_comparison(
                adata.obs[score_col], condition_values, positive_label
            )
            auc_rows.append({
                'program': prog_name,
                'display_name': PROGRAM_DISPLAY_NAMES[prog_name],
                'auc': metrics['auc'],
                'effect_size': metrics['effect_size'],
                'pvalue': metrics['pvalue'],
                'n_genes_found': gene_availability[prog_name],
                'n_genes_total': len(IMMUNE_PROGRAMS[prog_name]),
            })
        except ValueError as exc:
            print(f'[programs] could not compute AUC for {prog_name}: {exc}')

    auc_df = pd.DataFrame(auc_rows).sort_values('auc', ascending=False).reset_index(drop=True)
    auc_path = RESULTS_TABLES / 'kang_program_auc.csv'
    auc_df.to_csv(auc_path, index=False)
    print(f'[programs] saved program AUC table: {auc_path}')

    # Console AUROC summary
    print('\n[programs] AUROC ranking (stimulated vs control, all cell types pooled):')
    for _, row in auc_df.iterrows():
        bar = '\u2588' * int(round(row['auc'] * 20))
        print(f'  {row["display_name"]:<38}  AUC={row["auc"]:.3f}  {bar}')
    print()

    # ------------------------------------------------------------------
    # Build program x cell_type matrices for heatmaps
    # ------------------------------------------------------------------
    if cell_type_col is not None:
        # Sort cell types alphabetically for consistent figure layout.
        cell_types = sorted(adata.obs[cell_type_col].astype(str).unique())
        program_labels = [PROGRAM_DISPLAY_NAMES[p] for p in scored_programs]

        ctrl_matrix = pd.DataFrame(np.nan, index=program_labels, columns=cell_types)
        stim_matrix = pd.DataFrame(np.nan, index=program_labels, columns=cell_types)

        ct_values = adata.obs[cell_type_col].astype(str)
        for prog_name in scored_programs:
            score_col = f'prog_{prog_name}'
            disp = PROGRAM_DISPLAY_NAMES[prog_name]
            for ct in cell_types:
                mask_ctrl = (ct_values == ct) & (condition_values != positive_label)
                mask_stim = (ct_values == ct) & (condition_values == positive_label)
                n_ctrl = mask_ctrl.sum()
                n_stim = mask_stim.sum()
                if n_ctrl > 0:
                    ctrl_matrix.loc[disp, ct] = float(adata.obs.loc[mask_ctrl, score_col].mean())
                if n_stim > 0:
                    stim_matrix.loc[disp, ct] = float(adata.obs.loc[mask_stim, score_col].mean())

        delta_matrix = stim_matrix - ctrl_matrix

    # ------------------------------------------------------------------
    # Figures
    # ------------------------------------------------------------------

    # AUC barplot (no cell type required)
    plot_program_auc_barplot(
        auc_df,
        FIGURES / 'kang_program_auc_barplot.png',
        title='Immune program separation: stimulated vs control (AUROC)',
    )
    print('[programs] saved AUROC barplot: kang_program_auc_barplot.png')

    if cell_type_col is not None:
        # Option B — primary summary figure
        plot_program_heatmap_delta(
            delta_matrix,
            FIGURES / 'kang_program_heatmap_delta.png',
            title='Immune program score change: stimulated \u2212 control (Kang PBMC)',
        )
        print('[programs] saved delta heatmap (Option B): kang_program_heatmap_delta.png')

        # Option A — supporting figure showing absolute scores per condition
        plot_program_heatmap_sidebyside(
            ctrl_matrix,
            stim_matrix,
            FIGURES / 'kang_program_heatmap_sidebyside.png',
            ctrl_label='Control',
            stim_label='Stimulated',
            title='Immune program scores by condition (Kang PBMC)',
        )
        print('[programs] saved side-by-side heatmap (Option A): kang_program_heatmap_sidebyside.png')

    print(f'\n[programs] figures saved to {FIGURES}/')
    print(f'[programs] tables saved to {RESULTS_TABLES}/')
    print('[programs] completed pathway/program scoring')
