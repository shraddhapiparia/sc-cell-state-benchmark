#!/usr/bin/env python3
"""Pathway enrichment for within-cell-type DE results (Kang stim vs ctrl).

Uses Enrichr (via gseapy) to test the upregulated genes from each cell type
against the MSigDB Hallmark 2020 gene set collection.  Hallmark is a good
choice here: 50 well-curated gene sets covering the major immune and stress
response programs, with minimal redundancy.

Requires internet access (Enrichr API) and gseapy:
    pip install gseapy

Input
-----
  results/tables/kang_de_stim_vs_ctrl_per_cell_type.csv
      Produced by scripts/12_de_per_cell_type.py.

Upregulation filter
-------------------
  pvals_adj < 0.05  AND  logfoldchanges > 0.5
  This is more conservative than the broad filter used for the DE table so
  that ORA input genes are confidently upregulated rather than borderline.

Outputs
-------
  results/tables/kang_gsea_results.csv
      All enrichment results across cell types with adj p-value < 0.25.
      Columns: cell_type, Term, Overlap, P-value, Adjusted P-value,
               Combined Score, Genes, n_input, n_overlap

  figures/kang_gsea_top_pathways.png
      Dot plot: rows = top-5 pathways per cell type, columns = cell types.
      Dot colour = -log10(adj p-value), dot size = overlapping gene count.
"""

import sys
import pandas as pd

try:
    import gseapy as gp
except ImportError:
    sys.exit(
        '[gsea] gseapy is not installed.  Install it with:\n'
        '    pip install gseapy\n'
        'then re-run this script.'
    )

from sc_cell_state_benchmark.config import FIGURES, RESULTS_TABLES
from sc_cell_state_benchmark.data import ensure_dirs
from sc_cell_state_benchmark.plotting import plot_gsea_dotplot

# Filter applied to DE results to define the "upregulated" gene list per cell type
PADJ_THRESH = 0.05
LFC_THRESH  = 0.5

# Enrichr gene set library
GENE_SET_LIBRARY = 'MSigDB_Hallmark_2020'

# Enrichment significance threshold for saving results
ENRICH_PADJ_SAVE = 0.25   # broad save threshold; dotplot uses 0.05

MIN_GENES = 5  # minimum input genes to attempt enrichment for a cell type


if __name__ == '__main__':
    ensure_dirs([RESULTS_TABLES, FIGURES])

    de_path = RESULTS_TABLES / 'kang_de_stim_vs_ctrl_per_cell_type.csv'
    if not de_path.exists():
        sys.exit(
            f'[gsea] DE table not found at {de_path}.\n'
            'Run scripts/12_de_per_cell_type.py first.'
        )

    de = pd.read_csv(de_path)
    print(f'[gsea] loaded DE table: {len(de)} rows, '
          f'{de["cell_type"].nunique()} cell types')

    all_results = []

    for ct, grp in de.groupby('cell_type'):
        up = grp[
            (grp['pvals_adj'] < PADJ_THRESH) & (grp['logfoldchanges'] > LFC_THRESH)
        ]['names'].tolist()

        if len(up) < MIN_GENES:
            print(f'[gsea] SKIP {ct!r}: only {len(up)} upregulated genes '
                  f'(minimum {MIN_GENES})')
            continue

        print(f'[gsea] {ct!r}: {len(up)} upregulated genes -> Enrichr ...')

        try:
            enr = gp.enrichr(
                gene_list=up,
                gene_sets=GENE_SET_LIBRARY,
                organism='human',
                outdir=None,
            )
            res = enr.results.copy()
        except ValueError as exc:
            # Raised by gseapy for invalid parameters (e.g. bad organism value,
            # unrecognised gene set library name).
            print(f'[gsea] ERROR: invalid parameter for {ct!r}: {exc}')
            print(f'[gsea]   organism must be one of: human, mouse, yeast, fly, fish, worm')
            print(f'[gsea]   gene_sets must be a valid Enrichr library name')
            continue
        except Exception as exc:
            # Network errors, Enrichr API timeouts, or other connectivity issues.
            print(f'[gsea] WARNING: Enrichr API call failed for {ct!r}: {exc}')
            print(f'[gsea]   Check internet connection and Enrichr availability '
                  f'at https://maayanlab.cloud/Enrichr/')
            continue

        if res.empty:
            print(f'[gsea] no results returned for {ct!r}')
            continue

        # Normalise column names across gseapy versions
        res.columns = [c.strip() for c in res.columns]

        res['cell_type'] = ct
        res['n_input']   = len(up)

        # Parse overlap string "3/200" -> integer numerator
        if 'Overlap' in res.columns:
            res['n_overlap'] = res['Overlap'].str.split('/').str[0].astype(int)
        else:
            res['n_overlap'] = 0

        all_results.append(res)

        top3 = res.nsmallest(3, 'Adjusted P-value')[['Term', 'Adjusted P-value']]
        for _, row in top3.iterrows():
            print(f'  {row["Term"]:<50}  adj_p={row["Adjusted P-value"]:.3e}')

    if not all_results:
        sys.exit('[gsea] No enrichment results obtained. '
                 'Check internet access and Enrichr availability.')

    combined = pd.concat(all_results, ignore_index=True)

    # Keep standard columns; drop anything library-version-specific
    keep_cols = ['cell_type', 'Term', 'Overlap', 'P-value', 'Adjusted P-value',
                 'Combined Score', 'Genes', 'n_input', 'n_overlap']
    keep_cols = [c for c in keep_cols if c in combined.columns]
    combined = combined[keep_cols]

    sig_combined = combined[combined['Adjusted P-value'] < ENRICH_PADJ_SAVE]

    out_table = RESULTS_TABLES / 'kang_gsea_results.csv'
    sig_combined.to_csv(out_table, index=False)
    print(f'\n[gsea] saved enrichment table: {out_table}')
    print(f'[gsea] {len(sig_combined)} rows (adj p < {ENRICH_PADJ_SAVE}) '
          f'across {sig_combined["cell_type"].nunique()} cell types')

    out_fig = FIGURES / 'kang_gsea_top_pathways.png'
    plot_gsea_dotplot(
        sig_combined,
        out_fig,
        top_n=5,
        max_padj=0.05,
        title='Top enriched Hallmark pathways per cell type (stim upregulated genes)',
    )
    print(f'[gsea] saved dot plot: {out_fig}')
    print('[gsea] completed pathway enrichment analysis')
