#!/usr/bin/env python3
"""Exploratory cell-cell communication analysis for the Kang PBMC perturbation dataset.

Scores a curated panel of 16 ligand-receptor pairs using mean-product scoring
(mean log-expr of ligand in sender × mean log-expr of receptor in receiver)
separately for ctrl and stim conditions. Produces:

  Per-interaction score table
    results/tables/kang_communication_scores.csv

  Top interactions in the stimulated condition
    results/tables/kang_communication_top_stim.csv

  Sender × receiver summed-score matrix, both conditions
  (Option A — supporting figure)
    figures/kang_communication_heatmap_sidebyside.png

  Sender × receiver delta matrix (stim − ctrl)
  (Option B — primary summary figure)
    figures/kang_communication_heatmap_delta.png

  Arc-network diagram of increased communication under stimulation
    figures/kang_communication_network.png

Prerequisite: scripts/08_preprocess_kang_pbmc.py must have been run.

Important caveat
----------------
This is an exploratory analysis. Scores reflect gene expression levels, not
validated ligand-receptor binding or downstream signalling activity. Several
pairs (IFNG, IL10) are near-absent in this dataset and score near zero; they
are retained as transparent internal negative controls.
"""

import argparse

from sc_cell_state_benchmark.communication import (
    LR_PAIRS,
    build_sender_receiver_matrix,
    score_lr_pairs,
)
from sc_cell_state_benchmark.config import (
    FIGURES,
    KANG_PBMC_PREPROCESSED,
    RESULTS_TABLES,
)
from sc_cell_state_benchmark.data import ensure_dirs, load_anndata
from sc_cell_state_benchmark.plotting import (
    plot_communication_heatmap_delta,
    plot_communication_heatmap_sidebyside,
    plot_communication_network,
)

CONDITION_KEYWORDS = ['stim', 'condition', 'treat', 'group', 'status', 'label']
CELL_TYPE_KEYWORDS = ['cell_type', 'celltype', 'celllabels', 'cell_labels', 'annotation', 'subclass']


def parse_args():
    parser = argparse.ArgumentParser(
        description='Exploratory cell-cell communication analysis for a perturbation dataset.'
    )
    parser.add_argument('--input', type=str, default=None,
                        help='Path to preprocessed AnnData (.h5ad). Defaults to the Kang preprocessed path.')
    parser.add_argument('--condition', type=str, default=None,
                        help='Obs column for condition labels. Auto-detected if not provided.')
    parser.add_argument('--cell-type', type=str, default=None,
                        help='Obs column for cell type labels. Auto-detected if not provided.')
    parser.add_argument('--top-n', type=int, default=20,
                        help='Number of top interactions to include in the summary table.')
    return parser.parse_args()


def find_metadata_column(obs, keywords):
    for col in obs.columns:
        if any(kw in col.lower() for kw in keywords):
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
    raise ValueError(f'Could not infer stimulated label from: {list(values)}')


if __name__ == '__main__':
    args = parse_args()
    ensure_dirs([RESULTS_TABLES, FIGURES])

    # ------------------------------------------------------------------
    # Load data
    # ------------------------------------------------------------------
    input_path = args.input if args.input is not None else KANG_PBMC_PREPROCESSED
    adata = load_anndata(input_path)
    n_raw = adata.raw.n_vars if adata.raw is not None else adata.n_vars
    print(f'[comm] loaded dataset from {input_path}')
    print(f'[comm] {adata.n_obs} cells | {n_raw} genes available for scoring (raw matrix)')

    # ------------------------------------------------------------------
    # Resolve metadata columns
    # ------------------------------------------------------------------
    condition_col = args.condition or find_metadata_column(adata.obs, CONDITION_KEYWORDS)
    if condition_col is None or condition_col not in adata.obs.columns:
        raise ValueError('Could not find condition column. Provide --condition <col>.')
    print(f'[comm] condition column: {condition_col!r}')

    cell_type_col = args.cell_type or find_metadata_column(adata.obs, CELL_TYPE_KEYWORDS)
    if cell_type_col is None or cell_type_col not in adata.obs.columns:
        raise ValueError('Could not find cell type column. Provide --cell-type <col>.')
    print(f'[comm] cell type column:  {cell_type_col!r}')

    condition_values = adata.obs[condition_col].astype(str)
    unique_conditions = list(condition_values.unique())
    if len(unique_conditions) != 2:
        raise ValueError(f'Expected exactly 2 conditions, found {len(unique_conditions)}: {unique_conditions}')
    positive_label = infer_positive_label(unique_conditions)
    negative_label = [c for c in unique_conditions if c != positive_label][0]
    print(f'[comm] stimulated label: {positive_label!r} | control label: {negative_label!r}')

    cell_types = sorted(adata.obs[cell_type_col].astype(str).unique())
    ct_counts = adata.obs[cell_type_col].value_counts()
    print(f'[comm] {len(cell_types)} cell types: {cell_types}')
    if 'Megakaryocytes' in cell_types:
        n_mega = ct_counts.get('Megakaryocytes', 0)
        print(f'[comm] NOTE: Megakaryocytes n={n_mega} — scores for this type should not be interpreted confidently.')

    # ------------------------------------------------------------------
    # Gene-availability and mean-expression report
    # ------------------------------------------------------------------
    raw = adata.raw if adata.raw is not None else adata
    raw_genes = set(raw.var_names)
    all_lr_genes = sorted(set(g for pair in LR_PAIRS for g in (pair[0], pair[1])))

    print(f'\n[comm] LR gene availability ({len(LR_PAIRS)} pairs, {len(all_lr_genes)} unique genes):')
    absent = []
    for gene in all_lr_genes:
        tag = 'FOUND' if gene in raw_genes else 'ABSENT'
        if tag == 'ABSENT':
            absent.append(gene)
        print(f'  {tag}  {gene}')
    if absent:
        print(f'\n[comm] Genes absent from dataset (will score 0): {absent}')
    else:
        print('[comm] All LR genes present in dataset.')

    print('\n[comm] Per-pair gene status:')
    for ligand, receptor, name in LR_PAIRS:
        l_ok = ligand in raw_genes
        r_ok = receptor in raw_genes
        flag = '' if (l_ok and r_ok) else '  *** ligand absent ***' if not l_ok else '  *** receptor absent ***'
        print(f'  {name:<28} L:{ligand} {"OK" if l_ok else "ABSENT"}  |  R:{receptor} {"OK" if r_ok else "ABSENT"}{flag}')
    print()

    # ------------------------------------------------------------------
    # Score all LR pairs
    # ------------------------------------------------------------------
    scores_df = score_lr_pairs(adata, ct_col=cell_type_col, cond_col=condition_col)
    print(f'[comm] scored {len(LR_PAIRS)} pairs × {len(cell_types)} senders × {len(cell_types)} receivers × 2 conditions'
          f' = {len(scores_df)} rows')

    scores_path = RESULTS_TABLES / 'kang_communication_scores.csv'
    scores_df.to_csv(scores_path, index=False)
    print(f'[comm] saved full score table: {scores_path}')

    # Top interactions in stimulated condition
    top_stim = (
        scores_df[scores_df['condition'] == positive_label]
        .sort_values('lr_score', ascending=False)
        .head(args.top_n)
        .reset_index(drop=True)
    )
    top_stim_path = RESULTS_TABLES / 'kang_communication_top_stim.csv'
    top_stim.to_csv(top_stim_path, index=False)
    print(f'[comm] saved top-{args.top_n} stim interactions: {top_stim_path}')

    # ------------------------------------------------------------------
    # Build sender × receiver matrices
    # ------------------------------------------------------------------
    ctrl_matrix = build_sender_receiver_matrix(scores_df, negative_label, cell_types)
    stim_matrix = build_sender_receiver_matrix(scores_df, positive_label, cell_types)
    delta_matrix = stim_matrix - ctrl_matrix

    # ------------------------------------------------------------------
    # Console summaries
    # ------------------------------------------------------------------
    print(f'\n[comm] Top 5 interactions in stimulated condition:')
    for _, row in top_stim.head(5).iterrows():
        print(f'  {row["interaction_name"]:<28}  {row["sender"]:<22} → {row["receiver"]:<22}  score={row["lr_score"]:.4f}')

    # Top delta sender→receiver pairs
    delta_flat = (
        delta_matrix.stack()
        .reset_index()
        .rename(columns={'level_0': 'sender', 'level_1': 'receiver', 0: 'delta'})
    )
    delta_flat = delta_flat[delta_flat['sender'] != delta_flat['receiver']].sort_values('delta', ascending=False)

    print(f'\n[comm] Top 5 sender→receiver pairs by stim−ctrl delta (summed across all pairs):')
    for _, row in delta_flat.head(5).iterrows():
        print(f'  {row["sender"]:<22} → {row["receiver"]:<22}  delta={row["delta"]:.4f}')

    # Top per-pair deltas
    scores_pivot = scores_df.pivot_table(
        index=['interaction_name', 'sender', 'receiver'],
        columns='condition',
        values='lr_score',
        aggfunc='first',
    ).reset_index()
    scores_pivot.columns.name = None
    if positive_label in scores_pivot.columns and negative_label in scores_pivot.columns:
        scores_pivot['delta'] = scores_pivot[positive_label] - scores_pivot[negative_label]
        top_pair_delta = (
            scores_pivot[scores_pivot['sender'] != scores_pivot['receiver']]
            .sort_values('delta', ascending=False)
            .head(5)
        )
        print(f'\n[comm] Top 5 individual pair deltas (stim−ctrl):')
        for _, row in top_pair_delta.iterrows():
            print(f'  {row["interaction_name"]:<28}  {row["sender"]:<22} → {row["receiver"]:<22}  delta={row["delta"]:.4f}')
    print()

    # ------------------------------------------------------------------
    # Figures
    # ------------------------------------------------------------------
    plot_communication_heatmap_sidebyside(
        ctrl_matrix, stim_matrix,
        FIGURES / 'kang_communication_heatmap_sidebyside.png',
        ctrl_label='Control', stim_label='Stimulated',
        title='Cell-cell communication scores by condition (Kang PBMC)',
    )
    print('[comm] saved side-by-side heatmap: kang_communication_heatmap_sidebyside.png')

    plot_communication_heatmap_delta(
        delta_matrix,
        FIGURES / 'kang_communication_heatmap_delta.png',
        title='Cell-cell communication change: stimulated \u2212 control (Kang PBMC)',
    )
    print('[comm] saved delta heatmap: kang_communication_heatmap_delta.png')

    plot_communication_network(
        delta_matrix,
        FIGURES / 'kang_communication_network.png',
        title='Increased cell-cell communication: stimulated vs control (Kang PBMC)',
    )
    print('[comm] saved network figure: kang_communication_network.png')

    print(f'\n[comm] figures saved to {FIGURES}/')
    print(f'[comm] tables saved to {RESULTS_TABLES}/')
    print('[comm] completed exploratory cell-cell communication analysis')
