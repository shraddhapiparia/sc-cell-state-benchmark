# Next Steps

This document describes what to do before extending the benchmark to new biology.
Items are grouped by urgency and scope.

---

## Before any new biology is added

These are existing issues in the current codebase that should be resolved first
so the v1 state is clean and reproducible.

### Code cleanup
- **Fix `groupby` in `scripts/06_score_cell_states.py`:** add `observed=True` to match
  the fix already applied in `scripts/09_score_kang_interferon.py`.
- **Resolve `controls.py`:** the `random_gene_set` function is superseded by
  `matched_random_gene_sets` in `scoring.py`. Decide: remove `controls.py`, or promote
  `permute_labels` into the benchmark pipeline and keep the module with a clear purpose.
- **Remove TODO stubs in `evaluation.py`:** `compare_scores` and `robustness_test` are
  empty functions. Remove them or implement them — empty stubs with TODO comments
  suggest incomplete work to new contributors.
- **Update `copilot-instructions.md`:** it says "Install dependencies with
  `pip install -r requirements.txt` (to be created)" — `requirements.txt` now exists.

### Data hygiene
- **`data/kang.h5ad`:** there is a `kang.h5ad` file at `data/kang.h5ad` (outside the
  `raw/` and `processed/` subdirectory structure). Verify whether this is the same
  file as `data/raw/kang_pbmc_raw.h5ad`. If so, one copy can be removed.
- **`.gitignore` and figures:** the current `.gitignore` does not explicitly exclude
  figure PNGs. If figures should not be committed (they are generated outputs), add
  `figures/*.png` with `!figures/.gitkeep` to the gitignore.

### Documentation
- **`pyproject.toml` author field:** currently reads `"Your Name"` and
  `"your.email@example.com"` — update before making the repo public.

---

## Robustness tests (v1.1)

Before claiming that a method is reliable, these tests should be run:

1. **Permuted-label null:** shuffle ctrl/stim labels, re-score, and report AUROC.
   Expected: AUROC ≈ 0.5. If not, the scoring is leaking label information.

2. **Subsampling:** randomly subsample cells to 75%, 50%, 25% of the dataset and
   re-run scoring. Check whether AUROC degrades gracefully or drops sharply.

3. **Gene set size sensitivity:** run scoring with 5, 10, 15, and 20 genes (subsets
   of the full interferon signature). Report how AUROC changes with set size.

4. **Normalization sensitivity:** re-run scoring on unnormalized counts vs
   log-normalized counts to check whether preprocessing choices affect separation.

---

## Additional scoring methods (v1.1)

The benchmark currently has three methods. Obvious additions to consider:

- **ssGSEA:** widely used in bulk RNA-seq; single-cell adaptation is straightforward.
- **Full AUCell:** the current rank-based score is an approximation. The actual AUCell
  uses the area under the recovery curve, which is more principled.
- **Classifier-based upper bound:** train a logistic regression on ctrl/stim labels and
  report cross-validated AUROC. This sets a ceiling for what any scoring method can achieve
  on this dataset.

---

## Biological extensions (v2.0)

These are appropriate after v1.1 robustness work is complete and the codebase is stable.

### Additional datasets
- A second perturbation dataset using a different stimulus type (e.g., LPS, TNF-alpha,
  or a drug compound) to test whether the benchmark generalizes beyond interferon.
- A dataset where perturbation effects are weaker or cell-type-specific, to test
  method sensitivity in more challenging conditions.
- A disease-state atlas (e.g., published COVID-19 or autoimmune scRNA-seq) where
  "condition" is patient group rather than in-vitro stimulation.

### Gene set support
- Add a loader for MSigDB Hallmark gene sets so the benchmark can be run on any
  pathway, not just the hardcoded interferon signature.
- Evaluate whether AUROC varies substantially across different gene sets of similar size.

### Multi-gene-set comparison
- Run all methods on multiple gene sets in a single script, and produce a
  method × gene-set AUROC heatmap. This would show whether method rankings are
  consistent or gene-set dependent.

---

## What this project is not trying to do

To avoid scope creep, the following are explicitly out of scope for this benchmark:

- **Biological discovery:** the interferon gene set is used as a known positive control,
  not to make new claims about interferon biology.
- **Cell type deconvolution:** the benchmark assumes cell type annotations are provided
  as metadata; it does not attempt to infer them.
- **Trajectory or pseudotime analysis:** the focus is binary condition comparison
  (ctrl vs stim), not continuous state transitions.
- **Production-grade pipeline:** the numbered script approach is intentional for clarity.
  A workflow manager (Snakemake, Nextflow) is not a near-term goal unless the pipeline
  grows substantially in complexity.
