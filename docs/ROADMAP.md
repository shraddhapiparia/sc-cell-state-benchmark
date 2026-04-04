# Roadmap

This document tracks what has been completed, what is in progress, and what is
planned for future development. It is updated manually to reflect the actual state
of the repository.

---

## v1.0 — Interferon perturbation benchmark (current)

**Status: complete**

### Warm-up pipeline (PBMC3k)
- [x] Download and preprocess PBMC3k (scripts 01–02)
- [x] UMAP, Leiden clustering, QC plots (script 03)
- [x] Marker gene identification (script 04)
- [x] Cluster annotation (script 05)
- [x] Interferon gene scoring on PBMC3k as a sanity check (script 06)

### Main benchmark (Kang interferon PBMC)
- [x] Download and preprocess Kang-style interferon PBMC dataset (scripts 07–08)
- [x] Score cells with three methods: Scanpy `score_genes`, mean expression, rank-based (script 09)
- [x] Generate 25 matched random control gene sets and compare against real signature
- [x] Compute AUROC and Mann-Whitney U for each method
- [x] Per-cell-type and per-condition summary tables
- [x] Violin figures by condition and by cell type × condition
- [x] UMAP figures colored by score

---

## v1.1 — Robustness and cleanup (next)

**Status: not started**

### Immediate technical debt
- [ ] Fix `groupby` in `scripts/06_score_cell_states.py` to use `observed=True`
- [ ] Resolve `controls.py` — either remove the superseded stubs or promote them
- [ ] Remove the unimplemented TODO stubs from `evaluation.py` or implement them
- [ ] Clarify the `data/kang.h5ad` file at the root data level (outside `raw/processed/` structure)
- [ ] Update `.gitignore` to explicitly exclude figure PNGs if they should not be committed

### Robustness tests (planned, not implemented)
- [ ] Permuted-label negative control (label permutation as a second null)
- [ ] Subsampling robustness: does AUROC hold at 50%, 25% of cells?
- [ ] Sensitivity to gene set size: score 5 genes vs 10 vs 20
- [ ] Sensitivity to preprocessing: raw counts vs normalized vs log-normalized

### Additional scoring methods (planned)
- [ ] ssGSEA (single-sample GSEA)
- [ ] AUCell (full implementation, not the rank approximation currently used)
- [ ] Classifier-based score (logistic regression as an upper bound)

---

## v2.0 — Multi-dataset benchmark (future)

**Status: not started**

- [ ] Add a second perturbation dataset (different stimulus, different cell type composition)
- [ ] Add a disease-state dataset (e.g., a published COVID-19 or autoimmune scRNA-seq atlas)
- [ ] Evaluate scoring on cell types with weaker or more heterogeneous responses
- [ ] Add support for MSigDB or other curated gene set collections as inputs
- [ ] Structured comparison table across datasets and methods

---

## Notes on scope

This project benchmarks **scoring methods**, not biological pathways. The interferon
gene set is used as a well-validated positive control because the ctrl/stim labels are
known ground truth. The goal is to evaluate whether a scoring method can recover a
known signal — not to make new biological discoveries about interferon biology.
