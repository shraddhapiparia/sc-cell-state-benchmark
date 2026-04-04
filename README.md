# sc-cell-state-benchmark

Benchmarking cell-state scoring methods in single-cell RNA-seq, with a focus on
perturbation responses (e.g. interferon stimulation).

---

## Overview

This project has three layers, built on top of each other:

**Layer 1 — Preprocessing and annotation (PBMC3k warm-up, scripts 01–06):**
Standard 10x PBMC dataset used to validate the preprocessing and clustering pipeline.
Normalisation, HVG selection, UMAP, Leiden clustering, marker gene identification, and
cluster annotation. No perturbation signal — used as a dry run.

**Layer 2 — Scoring method benchmark (Kang interferon PBMC, scripts 07–09):**
A paired ctrl/stim dataset where PBMCs were stimulated with interferon-beta. Three
gene-set scoring methods are compared head-to-head: can each one separate stimulated
from control cells using a 10-gene interferon signature? Matched random gene sets
are used as a null control.

**Layer 3 — Biological program interpretation (scripts 10):**
Five curated immune programs are scored across all cells, producing a cell-type ×
condition score matrix. This layer shifts from asking "which method works?" to asking
"which biological programs change and in which cell types?" using the best-validated
scoring method from Layer 2.

---

## Methods benchmarked

| Method | Description |
|---|---|
| **Scanpy `score_genes`** | Background-corrected gene-set score (Seurat-style) |
| **Mean expression score** | Simple mean of normalized log expression across signature genes |
| **Rank-based score** | AUCell-inspired score using per-cell gene expression ranks |
| **Matched random controls** | 25 randomly sampled gene sets of equal size; used to confirm the real signal is not noise |

---

## Key findings

- Stimulated cells show consistently stronger interferon scores than control cells
  across all three scoring methods.
- The real 10-gene interferon signature (IFIT1, IFIT2, ISG15, MX1, OAS1, etc.)
  separates ctrl vs stim groups substantially better than matched random gene sets.
- The separation is visible both at the global level and within individual cell types
  (CD14 monocytes and NK cells show the largest ctrl-vs-stim delta).
- Score scaling differs across methods (scanpy scores are centered near zero;
  mean-expression and rank scores are on their own scales), but the relative ordering
  of stimulated vs control is consistent.

---

## Repo structure

```
sc-cell-state-benchmark/
├── data/
│   ├── raw/                  # Raw .h5ad files (not committed)
│   └── processed/            # Preprocessed AnnData files (not committed)
├── figures/                  # All output plots
├── results/
│   ├── tables/               # CSV score tables and summaries
│   ├── metrics/              # Evaluation metrics
│   ├── project_summary.md    # High-level project narrative
│   └── program_scoring_summary.md  # Biological interpretation of program scoring results
├── scripts/                  # Numbered pipeline scripts (run in order)
├── src/sc_cell_state_benchmark/   # Python package
│   ├── config.py             # Paths and constants
│   ├── data.py               # Data loading utilities
│   ├── gene_sets.py          # Curated immune program gene sets
│   ├── scoring.py            # Scoring method implementations
│   ├── evaluation.py         # Binary comparison metrics (AUROC, etc.)
│   └── plotting.py           # All figure-generation functions
├── tests/                    # Unit tests
├── environment.yml
└── requirements.txt
```

---

## Usage

Install:

```bash
pip install -e .
# or
conda env create -f environment.yml && conda activate sc-cell-state-benchmark
```

Run the full pipeline in script order:

```bash
# Warm-up: PBMC3k preprocessing (dry run)
python scripts/01_download_pbmc3k.py
python scripts/02_preprocess_pbmc3k.py
python scripts/03_plot_pbmc3k.py
python scripts/04_marker_genes.py
python scripts/05_annotate_clusters.py
python scripts/06_score_cell_states.py

# Main benchmark: perturbation dataset
python scripts/07_download_kang_pbmc.py
python scripts/08_preprocess_kang_pbmc.py
python scripts/09_score_kang_interferon.py

# Pathway/program scoring (interpretation layer)
python scripts/10_score_pathway_programs.py
```

Scripts 09 and 10 auto-detect the condition and cell-type columns by keyword matching,
or you can specify them explicitly. For the Kang dataset the condition column is named
`label` (values: `ctrl`, `stim`) and is not matched by the auto-detector — pass it
explicitly:

```bash
python scripts/09_score_kang_interferon.py \
  --input data/processed/kang_pbmc_preprocessed.h5ad \
  --condition label \
  --cell-type cell_type

python scripts/10_score_pathway_programs.py \
  --input data/processed/kang_pbmc_preprocessed.h5ad \
  --condition label \
  --cell-type cell_type
```

---

## Output files

After running `09_score_kang_interferon.py`:

| File | Description |
|---|---|
| `figures/kang_interferon_violin_by_condition.png` | Score distributions: ctrl vs stim |
| `figures/kang_interferon_by_cell_type_and_condition.png` | Per-cell-type violin, hue = condition |
| `figures/kang_interferon_umap_scanpy.png` | UMAP colored by Scanpy interferon score |
| `figures/kang_interferon_vs_random_controls.png` | Real gene set vs 25 random controls |
| `results/tables/kang_interferon_scores.csv` | Per-cell scores for all three methods |
| `results/tables/kang_interferon_summary_by_condition.csv` | Mean/median/std by condition |
| `results/tables/kang_interferon_summary_by_cell_type_and_condition.csv` | Per cell type × condition summary |
| `results/tables/kang_method_comparison.csv` | AUROC and other metrics per method |

After running `10_score_pathway_programs.py`:

| File | Description |
|---|---|
| `figures/kang_program_heatmap_delta.png` | **Primary figure:** stim − ctrl score per program × cell type (diverging) |
| `figures/kang_program_heatmap_sidebyside.png` | Supporting: raw scores for ctrl and stim side by side |
| `figures/kang_program_auc_barplot.png` | AUROC ranking: which programs best separate ctrl vs stim |
| `results/tables/kang_program_scores.csv` | Per-cell scores for all five programs |
| `results/tables/kang_program_summary_by_cell_type_and_condition.csv` | Mean/median/std per program × cell type × condition |
| `results/tables/kang_program_auc.csv` | AUROC, effect size, gene count per program |
| `results/program_scoring_summary.md` | Biological interpretation of program scoring results |

---

## Limitations

- **Dataset metadata assumptions:** The pipeline auto-detects condition and cell-type
  columns by keyword matching. If your dataset uses non-standard column names, pass
  them explicitly via `--condition` and `--cell-type`.
- **Score scaling differences:** The three methods use different scales and cannot be
  directly compared numerically — only their relative separation of groups is comparable.
- **No built-in data loader:** The Kang-style interferon dataset must be downloaded
  separately. `scripts/07_download_kang_pbmc.py` handles this, but requires an internet
  connection and access to the source repository.
- **Interferon gene set:** The 10-gene signature used here is a minimal representative
  set. Scores will vary with gene set choice; the benchmark framework is gene-set agnostic.
- **Curated program gene sets:** The five biological programs in `gene_sets.py` are
  manually curated minimal sets (7–10 genes each) intended to capture broad cell-state
  variation in human PBMCs. They are not derived from any licensed database and are not
  exhaustive. Cross-program gene membership has been minimised in this version to reduce
  score correlation and aid interpretability, not because biological overlap is incorrect
  — some immune genes (e.g., CXCL10 as both an ISG and an inflammatory chemokine)
  legitimately belong to more than one context. The module docstring documents specific
  cases where assignment was a pragmatic choice rather than a biological claim.
  Three programs (`antigen_presentation_mhc`, `cytotoxic_nk_like`, `monocyte_activation`)
  are expected to behave primarily as cell-type identity markers in this dataset rather
  than as perturbation-responsive programs — this is a feature of including them for
  biological context, not a limitation. Program scores should be read as approximate,
  qualitative indicators of biological state, not as precise pathway activity measurements.

---

## Planned Extensions

- Permuted-label negative control to complement the current random gene set null
- Subsampling and gene-set-size sensitivity tests
- Additional scoring methods: ssGSEA, full AUCell, classifier-based upper bound
- Second perturbation dataset with a different stimulus type
- MSigDB gene set loader so the benchmark is not tied to the hardcoded interferon signature

See [docs/ROADMAP.md](docs/ROADMAP.md) and [docs/NEXT_STEPS.md](docs/NEXT_STEPS.md) for details.

---

## Development

```bash
pytest tests/
```
