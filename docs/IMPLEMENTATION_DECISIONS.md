# Implementation Decisions

This document records non-obvious design choices and the reasoning behind them.
It is intended to help future contributors understand why the code is structured
the way it is, and to make it easier to change decisions that no longer apply.

---

## Scoring methods

### Why three scoring methods instead of one

The benchmark currently implements Scanpy `score_genes`, mean expression, and a
rank-based score. The goal is not to pick a winner — at this dataset size and with
this gene set, all three methods show similar AUROC (~0.99). The value of having
multiple methods is to check whether a finding is method-specific or consistent
across implementations. Conclusions that only hold for one method should be treated
with more caution than those that replicate across all three.

### Why `adata.raw` is used when available

Scanpy's `score_genes` and the other scoring functions access `adata.raw` when it
is set. This is because the `.raw` attribute stores the full gene matrix before
HVG filtering, giving scoring methods access to the complete gene vocabulary. If
`adata.raw` is not set, the functions fall back to `adata.X`. This is important for
the rank-based score in particular, where ranking genes across the full matrix is
different from ranking across a reduced HVG subset.

### Why the rank-based score is an approximation, not AUCell

The current rank-based score is a simplified approximation of the AUCell logic: it
ranks genes per cell and averages the ranks of target genes. It is not the full
AUCell implementation (which uses a recovery curve and an AUC per cell). A proper
AUCell implementation is noted in the roadmap as a future addition.

---

## Negative controls

### Why matched random gene sets instead of permuted labels

Two types of negative controls are used in single-cell scoring benchmarks:
1. **Random gene sets:** shuffles the gene identity while keeping cell structure intact.
   Answers: is this score explained by the number of genes alone?
2. **Permuted labels:** shuffles the condition labels while keeping scores intact.
   Answers: is the ctrl/stim separation better than chance?

The current implementation uses random gene sets (in `scoring.py:matched_random_gene_sets`).
Permuted-label controls are stubbed in `controls.py:permute_labels` but not yet wired
into the benchmark pipeline. Both should eventually be run.

### `controls.py` vs `scoring.py`

The `controls.py` module contains two early stubs (`random_gene_set`, `permute_labels`)
from the initial scaffold. `random_gene_set` has been superseded by the more complete
`matched_random_gene_sets` in `scoring.py`, which uses a seeded RNG, respects `adata.raw`,
and returns multiple sets. `controls.py` has not been removed to preserve the original
scaffold, but its contents should be either promoted to the main pipeline or removed
before a v1.1 release.

---

## Evaluation metrics

### Why AUROC as the primary metric

AUROC (area under the ROC curve) was chosen because:
- It is threshold-free and does not depend on choosing a cutoff.
- It measures rank-based separation, which is what gene-set scores are designed to produce.
- It is easy to interpret: 0.5 = random, 1.0 = perfect separation.

Caveats:
- AUROC is sensitive to class imbalance when classes are very unequal in size.
  The Kang dataset is approximately balanced (ctrl/stim), so this is not a concern here.
- AUROC alone does not capture effect size. The `effect_size` field in the comparison
  table (Cohen's d-like standardized mean difference) is reported alongside AUROC.

### Why Mann-Whitney U alongside AUROC

AUROC and the Mann-Whitney U statistic are mathematically equivalent for two-group
comparisons. Both are reported for completeness and to make the output compatible
with tools that expect one or the other.

---

## Data structure

### Why the pipeline uses numbered scripts instead of a workflow manager

The numbered script approach (01 through 09) keeps the pipeline explicit and readable
without requiring Snakemake, Nextflow, or similar infrastructure. Each script is
self-contained and can be re-run independently. The tradeoff is that dependency
tracking is manual. If the pipeline grows substantially, migrating to a DAG-based
workflow manager would be appropriate.

### Why `observed=True` in groupby

Pandas will raise a `FutureWarning` when grouping by a categorical column without
specifying `observed`. The `observed=True` argument requests only groups that actually
appear in the data (as opposed to all levels of the categorical). This is the correct
behavior here and will become the default in a future pandas version.

Note: script `06_score_cell_states.py` still uses `groupby` without `observed=True`
and should be updated.

---

## Gene set

### Why this specific 10-gene interferon signature

The 10 genes (IFIT1, IFIT2, IFIT3, ISG15, MX1, OAS1, OASL, RSAD2, IRF7, STAT1)
are a commonly cited subset of interferon-stimulated genes (ISGs) that are strongly
induced by type I interferon. They were chosen as a minimal, well-validated set for
benchmarking purposes — not as a comprehensive interferon signature. The benchmark
framework is gene-set agnostic; any gene list can be passed in.

The intersection of these 10 genes with the Kang dataset is reported at runtime.
Not all genes may be present in every dataset.

---

## Megakaryocytes caveat

The Kang dataset contains very few Megakaryocytes (n=10 ctrl, n=12 stim in the
preprocessed data). Summary statistics for this cell type should not be interpreted
with the same confidence as for cell types with hundreds or thousands of cells.
This caveat should be noted in any downstream analysis that references per-cell-type results.
