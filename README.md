# sc-cell-state-benchmark

Benchmarking RNA-based cell-state scoring methods in single-cell RNA-seq, with layered
biological interpretation using a paired interferon-stimulation PBMC dataset.

---

## What this repo demonstrates

This repository focuses on RNA-only benchmarking and interpretation; paired RNA+ATAC integration, TF/regulon inference, and dynamic regulatory analysis are reserved for a follow-up multiome repository.

- End-to-end single-cell RNA-seq preprocessing, clustering, and annotation
- Head-to-head comparison of three gene-set scoring methods on a perturbation dataset
- Null controls using matched random gene sets to confirm signal specificity
- Biological program scoring across cell types and conditions
- Exploratory cell-cell communication analysis using curated ligand-receptor pairs
- Reproducible, numbered pipeline with modular Python package

---

## Overview

The project is structured as four layers, each building on the previous:

**Layer 1 — RNA-only preprocessing and annotation (PBMC3k warm-up, scripts 01–06)**
Standard 10x PBMC dataset used to validate the preprocessing pipeline. Normalisation,
HVG selection, UMAP, Leiden clustering, marker gene identification, and cluster
annotation. No perturbation signal — used as a dry run.

**Layer 2 — RNA-only scoring method benchmark (Kang interferon PBMC, scripts 07–09)**
A paired ctrl/stim dataset where PBMCs were stimulated with interferon-beta. Three
gene-set scoring methods are compared: can each one separate stimulated from control
cells using a 10-gene interferon signature? 25 matched random gene sets are used as
null controls.

**Layer 3 — RNA-only biological program interpretation (script 10)**
Five curated immune programs scored across all cells, producing a cell-type × condition
score matrix. Shifts from "which method works?" to "which biological programs change
and in which cell types?"

**Layer 4 — Optional exploratory cell-cell communication (script 11)**
16 curated ligand-receptor pairs scored using mean-product expression (mean ligand
expression in sender cells × mean receptor expression in receiver cells) per condition.
**This is an exploratory layer only** — scores reflect co-expression, not validated
intercellular signalling.

---

## Datasets

| Dataset | Role | Cells | Notes |
|---|---|---|---|
| PBMC3k (10x Genomics) | Warm-up / dry run | ~2,700 | No perturbation |
| Kang PBMC (IFN-beta stimulation) | Main benchmark | 24,562 | Paired ctrl/stim |

---

## Methods benchmarked

| Method | Description |
|---|---|
| **Scanpy `score_genes`** | Background-corrected gene-set score (Seurat-style) |
| **Mean expression score** | Mean normalised log expression across signature genes |
| **Rank-based score** | AUCell-inspired per-cell expression rank score |
| **Matched random controls** | 25 randomly sampled gene sets of equal size as null |

---

## Key findings

- All three scoring methods separate stimulated from control cells with AUROC ≥ 0.985.
- The 10-gene interferon signature achieves ~1.7 SD effect size; matched random gene
  sets produce near-zero effect sizes, confirming the signal is gene-set-specific.
- Program scoring shows the IFN-α/β response (AUC 0.991) and inflammatory response
  (AUC 0.655) as perturbation-driven; antigen presentation, cytotoxic, and monocyte
  activation programs primarily reflect cell-type identity in this dataset.
- Exploratory communication scoring suggests monocytes and DCs as dominant senders of increased
  signal under stimulation, driven by ISG-induced CXCL10 expression.

---

## Repo structure

```
sc-cell-state-benchmark/
├── data/
│   ├── raw/                      # Raw .h5ad files (not committed)
│   └── processed/                # Preprocessed AnnData files (not committed)
├── docs/
│   ├── ROADMAP.md                # Versioned roadmap and status
│   ├── IMPLEMENTATION_DECISIONS.md  # Design rationale
│   └── NEXT_STEPS.md             # Ordered pre-extension task list
├── figures/                      # All generated plots (not committed)
├── results/
│   ├── tables/                   # CSV score tables (not committed)
│   ├── project_summary.md        # High-level project summary
│   ├── program_scoring_summary.md  # Program scoring interpretation
│   └── communication_summary.md  # Cell-cell communication interpretation
├── scripts/                      # Numbered pipeline scripts (run in order)
│   ├── 01–06                     # PBMC3k warm-up
│   ├── 07–09                     # Kang benchmark
│   ├── 10                        # Program scoring
│   └── 11                        # Cell-cell communication
├── src/sc_cell_state_benchmark/  # Python package
│   ├── config.py                 # Paths and constants
│   ├── data.py                   # Data loading utilities
│   ├── gene_sets.py              # Curated immune program gene sets
│   ├── communication.py          # Curated LR pairs and mean-product scoring
│   ├── scoring.py                # Scoring method implementations
│   ├── evaluation.py             # Binary comparison metrics (AUROC, etc.)
│   └── plotting.py               # All figure-generation functions
├── tests/                        # Import smoke tests
├── environment.yml
├── requirements.txt
└── pyproject.toml
```

---

## Setup and reproducibility

```bash
# Option 1: conda (recommended)
conda env create -f environment.yml
conda activate sc-cell-state-benchmark

# Option 2: pip
pip install -e .
```

All scripts use `RANDOM_SEED = 42` (defined in `config.py`). Results are fully
reproducible given the same input data files and package versions pinned in
`environment.yml`.

---

## Usage

Run scripts in numbered order. Each script is self-contained and can be re-run
independently once its prerequisites have been run.

```bash
# Layer 1 — PBMC3k warm-up
python scripts/01_download_pbmc3k.py
python scripts/02_preprocess_pbmc3k.py
python scripts/03_plot_pbmc3k.py
python scripts/04_marker_genes.py
python scripts/05_annotate_clusters.py
python scripts/06_score_cell_states.py

# Layer 2 — Kang benchmark
python scripts/07_download_kang_pbmc.py
python scripts/08_preprocess_kang_pbmc.py
python scripts/09_score_kang_interferon.py --condition label --cell-type cell_type

# Layer 3 — Program scoring
python scripts/10_score_pathway_programs.py --condition label --cell-type cell_type

# Layer 4 — Cell-cell communication
python scripts/11_cell_communication.py --condition label --cell-type cell_type
```

Scripts 09–11 auto-detect the condition and cell-type columns by keyword matching.
For the Kang dataset the condition column is `label` (values: `ctrl`, `stim`) — pass
it explicitly as shown above.

Full CLI options:

```bash
python scripts/09_score_kang_interferon.py \
  --input data/processed/kang_pbmc_preprocessed.h5ad \
  --condition label \
  --cell-type cell_type

python scripts/10_score_pathway_programs.py \
  --input data/processed/kang_pbmc_preprocessed.h5ad \
  --condition label \
  --cell-type cell_type

python scripts/11_cell_communication.py \
  --input data/processed/kang_pbmc_preprocessed.h5ad \
  --condition label \
  --cell-type cell_type
```

---

## Output files

**After `09_score_kang_interferon.py`:**

| File | Description |
|---|---|
| `figures/kang_interferon_violin_by_condition.png` | Score distributions: ctrl vs stim |
| `figures/kang_interferon_by_cell_type_and_condition.png` | Per-cell-type violin, hue = condition |
| `figures/kang_interferon_umap_scanpy.png` | UMAP coloured by Scanpy interferon score |
| `figures/kang_interferon_vs_random_controls.png` | Real gene set vs 25 random controls |
| `results/tables/kang_interferon_scores.csv` | Per-cell scores for all three methods |
| `results/tables/kang_interferon_summary_by_condition.csv` | Mean/median/std by condition |
| `results/tables/kang_interferon_summary_by_cell_type_and_condition.csv` | Per cell type × condition |
| `results/tables/kang_method_comparison.csv` | AUROC and effect size per method |

**After `10_score_pathway_programs.py`:**

| File | Description |
|---|---|
| `figures/kang_program_heatmap_delta.png` | **Primary:** stim − ctrl score per program × cell type |
| `figures/kang_program_heatmap_sidebyside.png` | Supporting: raw scores ctrl and stim side by side |
| `figures/kang_program_auc_barplot.png` | AUROC ranking by program |
| `results/tables/kang_program_scores.csv` | Per-cell scores for all five programs |
| `results/tables/kang_program_summary_by_cell_type_and_condition.csv` | Per program × cell type × condition |
| `results/tables/kang_program_auc.csv` | AUROC, effect size, gene count per program |
| `results/program_scoring_summary.md` | Biological interpretation of program scoring |

**After `11_cell_communication.py`:**

| File | Description |
|---|---|
| `figures/kang_communication_heatmap_delta.png` | **Primary:** stim − ctrl summed LR score per sender × receiver |
| `figures/kang_communication_heatmap_sidebyside.png` | Supporting: absolute LR scores ctrl and stim |
| `figures/kang_communication_network.png` | Arc-network of top increased sender→receiver pairs |
| `results/tables/kang_communication_scores.csv` | Per-pair × sender × receiver × condition (2,048 rows) |
| `results/tables/kang_communication_top_stim.csv` | Top 20 interactions in stimulated condition |
| `results/communication_summary.md` | Interpretation, assumptions, and caveats |

---

## Limitations

- **Dataset metadata assumptions:** Auto-detection of condition and cell-type columns
  uses keyword matching. For the Kang dataset pass `--condition label --cell-type cell_type`
  explicitly.
- **Score scaling:** The three scoring methods use different scales and cannot be
  compared numerically — only their relative separation of groups is meaningful.
- **Data access:** The Kang PBMC dataset must be provided locally or downloaded via
  `scripts/07_download_kang_pbmc.py`, which requires an internet connection.
- **Gene set size:** All gene signatures are minimal curated sets (7–10 genes). Scores
  will vary with gene set choice; the framework is gene-set agnostic.
- **Program gene sets:** Cross-program gene membership is minimised for interpretability,
  not because biological overlap is incorrect. Some genes (e.g., CXCL10 as both an ISG
  and an inflammatory chemokine) legitimately belong to multiple programs. The
  `gene_sets.py` module docstring documents specific placement decisions.
- **Cell-type identity programs:** `antigen_presentation_mhc`, `cytotoxic_nk_like`,
  and `monocyte_activation` are expected to behave primarily as cell-type identity
  markers in this dataset, not as perturbation-responsive programs. This is informative
  biological context, not a limitation.
- **Cell-cell communication is exploratory:** Mean-product LR scores reflect
  co-expression of ligands and receptors, not validated intercellular signalling.
  Results should not be cited as evidence of specific interactions without independent
  experimental validation.
- **Megakaryocytes:** n=22 total cells (10 ctrl, 12 stim). Scores for this cell type
  should not be interpreted with confidence.
- **Single dataset:** All results are from one in-vitro IFN-beta stimulation experiment.
  Generalisability requires additional datasets.

---

## Future extensions

- Permuted-label null controls for program and communication scores
- Additional scoring methods: ssGSEA, full AUCell, classifier-based upper bound
- Subsampling and gene-set-size sensitivity tests
- Second perturbation dataset (different stimulus or disease context)
- MSigDB or GO gene set loader
- Richer communication statistics with permutation-based significance testing
- Trajectory/pseudotime analysis for continuous cell-state transitions

## Future extensions
Within this repo:
- additional RNA scoring methods
- permutation nulls
- subsampling robustness
- extra perturbation datasets
- gene-set loaders

Planned follow-up repository:
- paired RNA+ATAC integration
- WNN-based multimodal analysis
- TF / regulon inference
- dynamic regulatory modeling and pseudotime

See [docs/ROADMAP.md](docs/ROADMAP.md) and [docs/NEXT_STEPS.md](docs/NEXT_STEPS.md) for details.

---

## Development

```bash
pytest tests/
```
