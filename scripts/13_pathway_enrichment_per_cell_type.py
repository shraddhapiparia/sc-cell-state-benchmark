#!/usr/bin/env python3
"""Pathway enrichment for within-cell-type DE results (Kang stim vs ctrl).

Three complementary analyses are run for each cell type:

ORA (upregulated)
    Over-representation analysis on upregulated DEGs (padj < 0.05, logFC > 0.5).
    Uses the full set of genes tested for that cell type as the background
    universe -- the experiment-specific background required for a valid
    hypergeometric test (Wijesooriya et al. 2024, PMC11557902).

ORA (downregulated)
    Same approach for downregulated DEGs (padj < 0.05, logFC < -0.5).
    Running ORA separately for each direction has substantially higher power
    than a combined or one-directional analysis (Zyla et al. 2014, PMC3899863).

Preranked GSEA
    Enrichment analysis on the complete ranked gene list.
    Ranking metric: sign(logFC) * -log10(padj + epsilon).
    Distinct from ORA: uses the full ranked list and permutation statistics,
    not a binary foreground/background split.

Requires internet access (Enrichr API + MSigDB Hallmark download).
Requires gseapy >= 1.0:  pip install "gseapy>=1.0"

Input
-----
  results/tables/kang_de_stim_vs_ctrl_per_cell_type.csv
      Produced by scripts/12_de_per_cell_type.py.
      Columns: cell_type, names, scores, logfoldchanges, pvals, pvals_adj

Outputs
-------
  results/tables/kang_ora_up_results.csv
      ORA results for upregulated genes, adj p < 0.25.

  results/tables/kang_ora_down_results.csv
      ORA results for downregulated genes, adj p < 0.25.

  results/tables/kang_gsea_preranked_results.csv
      Preranked GSEA results, FDR q-value < 0.25.

  figures/kang_ora_up_top_pathways.png
  figures/kang_ora_down_top_pathways.png
  figures/kang_gsea_preranked_top_pathways.png
"""

import sys

import numpy as np
import pandas as pd

try:
    import gseapy as gp
except ImportError:
    sys.exit(
        '[enrichment] gseapy is not installed.  Install it with:\n'
        '    pip install "gseapy>=1.0"\n'
        'then re-run this script.'
    )

from sc_cell_state_benchmark.config import FIGURES, RESULTS_TABLES
from sc_cell_state_benchmark.data import ensure_dirs
from sc_cell_state_benchmark.plotting import plot_enrichment_dotplot, plot_gsea_prerank_dotplot

# ── ORA thresholds ────────────────────────────────────────────────────────────
PADJ_THRESH   = 0.05   # foreground: adjusted p-value cutoff
LFC_THRESH    = 0.5    # foreground: absolute log-fold-change cutoff
ORA_SAVE_PADJ = 0.25   # save rows with ORA adj p below this
MIN_GENES     = 5      # skip ORA if fewer foreground genes than this

# ── Preranked GSEA parameters ─────────────────────────────────────────────────
GSEA_FDR_SAVE  = 0.25   # save rows with FDR q-value below this (Broad recommendation)
PRERANK_PERMS  = 1000
PRERANK_SEED   = 42
MIN_SET_SIZE   = 15     # Broad Institute default
MAX_SET_SIZE   = 500    # Broad Institute default

# ── Gene set library ──────────────────────────────────────────────────────────
GENE_SET_LIBRARY = 'MSigDB_Hallmark_2020'


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _ranking_metric(df: pd.DataFrame) -> pd.Series:
    """Compute sign(logFC) * (-log10(padj + epsilon) + 1e-6 * abs(logFC)).

    Positive values rank at the top (upregulated, small p-value).
    Negative values rank at the bottom (downregulated, small p-value).
    The tiny abs(logFC) term is a deterministic tie-breaker that reduces
    rank collisions among genes with equal adjusted p-values.
    """
    lfc = df['logfoldchanges'].values
    sign = np.sign(lfc)
    neg_log_p = -np.log10(df['pvals_adj'].values + 1e-300)
    return pd.Series(sign * (neg_log_p + 1e-6 * np.abs(lfc)),
                     index=df['names'].values, name='rank')


def _normalise_ora_columns(df: pd.DataFrame) -> pd.DataFrame:
    df.columns = [c.strip() for c in df.columns]
    if 'Overlap' in df.columns and 'n_overlap' not in df.columns:
        df['n_overlap'] = df['Overlap'].str.split('/').str[0].astype(int)
    return df


# ─────────────────────────────────────────────────────────────────────────────
# ORA
# ─────────────────────────────────────────────────────────────────────────────

def run_ora(
    gene_list: list,
    background: list,
    gene_sets: dict,
    direction: str,
    cell_type: str,
) -> 'pd.DataFrame | None':
    """Run ORA with an explicit experiment-specific background.

    Uses gseapy.enrichr with background=<list> which triggers a local
    Fisher exact test rather than the Enrichr API, ensuring the universe
    is restricted to the genes actually tested in this experiment.
    """
    if len(gene_list) < MIN_GENES:
        print(f'[ora] SKIP {cell_type!r} ({direction}): {len(gene_list)} genes '
              f'< minimum {MIN_GENES}')
        return None

    print(f'[ora] {cell_type!r} ({direction}): {len(gene_list)} foreground genes, '
          f'background={len(background)} ...')

    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=gene_sets,
            background=background,
            outdir=None,
        )
        res = enr.results.copy()
    except Exception as exc:
        print(f'[ora] WARNING: ORA failed for {cell_type!r} ({direction}): {exc}')
        return None

    if res.empty:
        print(f'[ora] no results for {cell_type!r} ({direction})')
        return None

    res = _normalise_ora_columns(res)
    res['cell_type'] = cell_type
    res['direction'] = direction
    res['n_input'] = len(gene_list)

    top3 = res.nsmallest(3, 'Adjusted P-value')[['Term', 'Adjusted P-value']]
    for _, row in top3.iterrows():
        print(f'  {row["Term"]:<50}  adj_p={row["Adjusted P-value"]:.3e}')

    return res


# ─────────────────────────────────────────────────────────────────────────────
# Preranked GSEA
# ─────────────────────────────────────────────────────────────────────────────

def run_prerank(
    de_grp: pd.DataFrame,
    gene_sets: dict,
    cell_type: str,
) -> 'pd.DataFrame | None':
    """Run preranked GSEA on the full tested gene list for one cell type."""
    rnk = _ranking_metric(de_grp)
    # Deduplicate gene names; keep entry with highest absolute rank value.
    rnk = rnk.groupby(level=0).apply(lambda s: s.iloc[s.abs().argmax()])
    rnk = rnk.sort_values(ascending=False)

    if len(rnk) < MIN_SET_SIZE * 2:
        print(f'[gsea] SKIP {cell_type!r}: only {len(rnk)} genes in ranked list')
        return None

    print(f'[gsea] {cell_type!r}: {len(rnk)} ranked genes -> prerank '
          f'({PRERANK_PERMS} permutations) ...')

    try:
        pre = gp.prerank(
            rnk=rnk,
            gene_sets=gene_sets,
            min_size=MIN_SET_SIZE,
            max_size=MAX_SET_SIZE,
            permutation_num=PRERANK_PERMS,
            outdir=None,
            seed=PRERANK_SEED,
            no_plot=True,
            threads=1,
        )
        res = pre.res2d.copy() if hasattr(pre, 'res2d') else pd.DataFrame(pre.results).T.copy()
    except Exception as exc:
        print(f'[gsea] WARNING: prerank failed for {cell_type!r}: {exc}')
        return None

    if res.empty:
        print(f'[gsea] no results for {cell_type!r}')
        return None

    res.columns = [c.strip() for c in res.columns]
    # Term may be the index rather than a column depending on gseapy version.
    if 'Term' not in res.columns:
        res = res.reset_index().rename(columns={'index': 'Term'})
    res['cell_type'] = cell_type

    fdr_col = next((c for c in res.columns if 'FDR' in c.upper()), None)
    if fdr_col:
        res[fdr_col] = pd.to_numeric(res[fdr_col], errors='coerce')
        res['NES'] = pd.to_numeric(res['NES'], errors='coerce')
        preview = res.dropna(subset=[fdr_col]).nsmallest(3, fdr_col)
        for _, row in preview.iterrows():
            try:
                print(f'  {str(row["Term"]):<50}  NES={float(row["NES"]):+.2f}  '
                      f'fdr={float(row[fdr_col]):.3e}')
            except Exception:
                pass

    return res


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    ensure_dirs([RESULTS_TABLES, FIGURES])

    de_path = RESULTS_TABLES / 'kang_de_stim_vs_ctrl_per_cell_type.csv'
    if not de_path.exists():
        sys.exit(
            f'[enrichment] DE table not found at {de_path}.\n'
            'Run scripts/12_de_per_cell_type.py first.'
        )

    de = pd.read_csv(de_path)
    print(f'[enrichment] loaded DE table: {len(de):,} rows, '
          f'{de["cell_type"].nunique()} cell types')

    print(f'[enrichment] fetching gene set library: {GENE_SET_LIBRARY} ...')
    try:
        gene_sets = gp.get_library(name=GENE_SET_LIBRARY, organism='Human')
    except Exception as exc:
        sys.exit(f'[enrichment] failed to fetch gene set library: {exc}')
    print(f'[enrichment] {len(gene_sets)} gene sets loaded')

    ora_up_rows:   list = []
    ora_down_rows: list = []
    gsea_rows:     list = []

    for ct, grp in de.groupby('cell_type'):
        background = list(dict.fromkeys(grp['names'].tolist()))  # unique, order-preserving

        up = grp[
            (grp['pvals_adj'] < PADJ_THRESH) & (grp['logfoldchanges'] > LFC_THRESH)
        ]['names'].tolist()

        down = grp[
            (grp['pvals_adj'] < PADJ_THRESH) & (grp['logfoldchanges'] < -LFC_THRESH)
        ]['names'].tolist()

        res_up = run_ora(up, background, gene_sets, 'up', ct)
        if res_up is not None:
            ora_up_rows.append(res_up)

        res_dn = run_ora(down, background, gene_sets, 'down', ct)
        if res_dn is not None:
            ora_down_rows.append(res_dn)

        res_gsea = run_prerank(grp, gene_sets, ct)
        if res_gsea is not None:
            gsea_rows.append(res_gsea)

    def _flag_low_confidence(df: pd.DataFrame, overlap_threshold: int = 3) -> pd.DataFrame:
        """Add low_confidence_overlap and confidence_note columns."""
        flag = df['n_overlap'] <= overlap_threshold if 'n_overlap' in df.columns else False
        df['low_confidence_overlap'] = flag
        df['confidence_note'] = ''
        if isinstance(flag, pd.Series):
            df.loc[flag, 'confidence_note'] = f'overlap <= {overlap_threshold} genes'
        return df

    # ── ORA up ────────────────────────────────────────────────────────────────
    _ORA_KEEP = ['cell_type', 'direction', 'Term', 'Overlap', 'P-value',
                 'Adjusted P-value', 'Combined Score', 'Genes', 'n_input', 'n_overlap']

    if ora_up_rows:
        ora_up = pd.concat(ora_up_rows, ignore_index=True)
        ora_up = ora_up[[c for c in _ORA_KEEP if c in ora_up.columns]]
        ora_up = _flag_low_confidence(ora_up)
        sig_up = ora_up[ora_up['Adjusted P-value'] < ORA_SAVE_PADJ]

        out = RESULTS_TABLES / 'kang_ora_up_results.csv'
        sig_up.to_csv(out, index=False)
        print(f'\n[ora] saved ORA (up) table: {out}')
        print(f'[ora] {len(sig_up)} rows across '
              f'{sig_up["cell_type"].nunique()} cell types (adj p < {ORA_SAVE_PADJ})')

        out_fig = FIGURES / 'kang_ora_up_top_pathways.png'
        plot_enrichment_dotplot(
            sig_up,
            out_fig,
            top_n=5,
            max_padj=0.05,
            title='ORA -- upregulated genes: top Hallmark pathways per cell type (stim vs ctrl)',
        )
        print(f'[ora] saved dot plot: {out_fig}')
    else:
        print('[ora] no significant ORA (up) results')

    # ── ORA down ──────────────────────────────────────────────────────────────
    if ora_down_rows:
        ora_dn = pd.concat(ora_down_rows, ignore_index=True)
        ora_dn = ora_dn[[c for c in _ORA_KEEP if c in ora_dn.columns]]
        ora_dn = _flag_low_confidence(ora_dn)
        sig_dn = ora_dn[ora_dn['Adjusted P-value'] < ORA_SAVE_PADJ]

        out = RESULTS_TABLES / 'kang_ora_down_results.csv'
        sig_dn.to_csv(out, index=False)
        print(f'\n[ora] saved ORA (down) table: {out}')
        print(f'[ora] {len(sig_dn)} rows across '
              f'{sig_dn["cell_type"].nunique()} cell types (adj p < {ORA_SAVE_PADJ})')

        out_fig = FIGURES / 'kang_ora_down_top_pathways.png'
        plot_enrichment_dotplot(
            sig_dn,
            out_fig,
            top_n=5,
            max_padj=0.05,
            title='ORA -- downregulated genes: top Hallmark pathways per cell type (stim vs ctrl)',
        )
        print(f'[ora] saved dot plot: {out_fig}')
    else:
        print('[ora] no significant ORA (down) results')

    # ── Preranked GSEA ────────────────────────────────────────────────────────
    if gsea_rows:
        gsea = pd.concat(gsea_rows, ignore_index=True)
        gsea.columns = [c.strip() for c in gsea.columns]

        fdr_col = next((c for c in gsea.columns if 'FDR' in c.upper()), 'FDR q-val')
        gsea[fdr_col] = pd.to_numeric(gsea[fdr_col], errors='coerce')
        sig_gsea = gsea[gsea[fdr_col].fillna(1.0) < GSEA_FDR_SAVE]

        out = RESULTS_TABLES / 'kang_gsea_preranked_results.csv'
        sig_gsea.to_csv(out, index=False)
        print(f'\n[gsea] saved preranked GSEA table: {out}')
        print(f'[gsea] {len(sig_gsea)} rows across '
              f'{sig_gsea["cell_type"].nunique()} cell types (FDR < {GSEA_FDR_SAVE})')

        out_fig = FIGURES / 'kang_gsea_preranked_top_pathways.png'
        plot_gsea_prerank_dotplot(
            sig_gsea,
            out_fig,
            fdr_col=fdr_col,
            top_n=5,
            fdr_threshold=GSEA_FDR_SAVE,
            title='Preranked GSEA: top Hallmark pathways per cell type (stim vs ctrl)',
        )
        print(f'[gsea] saved dot plot: {out_fig}')
    else:
        print('[gsea] no significant preranked GSEA results')

    print('\n[enrichment] pathway enrichment analysis complete')
