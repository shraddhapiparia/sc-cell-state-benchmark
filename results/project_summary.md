# Project Summary: sc-cell-state-benchmark

## Objective

Benchmark three gene-set scoring methods on single-cell RNA-seq data and evaluate
whether each method can distinguish stimulated cells from controls using a biologically
meaningful signature gene set.

---

## Datasets

| Dataset | Role | Cells | Notes |
|---|---|---|---|
| PBMC3k (10x) | Warm-up / dry run | ~2,700 | Standard clustering pipeline; no perturbation signal |
| Kang-style interferon PBMC | Main benchmark | ~29,000 | Paired ctrl/stim; interferon-beta stimulation |

---

## Methods

- **Scanpy `score_genes`:** Background-corrected score comparing target genes to a
  random reference set of equal size. Centered near zero.
- **Mean expression score:** Mean normalized log expression of signature genes per cell.
  Simple and interpretable.
- **Rank-based score (AUCell-inspired):** Scores each cell by the expression rank of
  signature genes within that cell's overall gene ranking.
- **Matched random controls:** 25 randomly drawn gene sets of the same size as the
  interferon signature. Used to confirm that separation is gene-set-specific, not an
  artifact of set size or data structure.

Interferon signature used (10 genes):
`IFIT1, IFIT2, IFIT3, ISG15, MX1, OAS1, OASL, RSAD2, IRF7, STAT1`

---

## Key Findings

1. All three scoring methods cleanly separate stimulated from control cells in the
   Kang interferon dataset.
2. The real interferon signature achieves higher mean scores than any of the 25
   matched random gene sets, confirming the signal is biologically specific.
3. Score separation is consistent across cell types, with CD14 monocytes and NK cells
   showing among the largest ctrl-vs-stim deltas.
4. Score scales differ across methods (scanpy: centered ~0; mean-expression and rank:
   positive, larger magnitude), but group separation direction is consistent.

---

## Next Steps

- Extend to additional perturbation datasets (e.g., different cytokines, drug treatments).
- Add more scoring methods: GSVA, AUCell (full implementation), ssGSEA.
- Benchmark on datasets with weaker or noisier perturbation signals.
- Add statistical testing (e.g., Mann-Whitney U) to the condition comparison output.
- Package gene-set loading to support MSigDB or other curated collections.
