# sc-cell-state-benchmark

Benchmarking cell-state scoring methods in single-cell RNA-seq, with a focus on
perturbation responses (e.g. interferon stimulation).

---

## Overview

This project benchmarks three gene-set scoring methods on single-cell RNA-seq data
and evaluates whether they can separate stimulated from control cells using a
biologically meaningful gene set.

**Warm-up dataset (PBMC3k):** Standard 10x PBMC dataset used to validate the
preprocessing and clustering pipeline. No perturbation signal — used as a dry run only.

**Main benchmark dataset (Kang-style interferon PBMC):** A paired ctrl/stim dataset
where PBMCs were stimulated with interferon-beta. The benchmark asks: can each scoring
method separate stimulated from control cells using a 10-gene interferon signature?

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
│   └── metrics/              # Evaluation metrics
├── scripts/                  # Numbered pipeline scripts (run in order)
├── src/sc_cell_state_benchmark/   # Python package
│   ├── config.py             # Paths and constants
│   ├── data.py               # Data loading utilities
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
```

The scoring script auto-detects the condition and cell-type columns, or you can
specify them explicitly:

```bash
python scripts/09_score_kang_interferon.py \
  --input data/processed/kang_pbmc_preprocessed.h5ad \
  --condition condition \
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
| `results/tables/kang_interferon_scores.csv` | Per-cell scores for all methods |
| `results/tables/kang_interferon_summary_by_condition.csv` | Mean/median/std by condition |
| `results/tables/kang_interferon_summary_by_cell_type_and_condition.csv` | Per cell type × condition summary |
| `results/tables/kang_method_comparison.csv` | AUROC and other metrics per method |

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

---

## Development

```bash
pytest tests/
```
