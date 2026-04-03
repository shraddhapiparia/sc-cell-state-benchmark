# Project Guidelines

## Architecture
This is a single-cell RNA-seq benchmarking project evaluating cell-state scoring methods. The project uses Python with Scanpy and AnnData for analysis.

Key components:
- `data/`: Datasets and metadata
- `notebooks/`: Analysis notebooks (01_pbmc3k_basics through planned benchmarking)
- `src/`: Reusable scoring and evaluation functions
- `results/`: Summary tables and metrics
- `figures/`: Publication-quality plots

Focus on distinguishing true biological signal from noise using known perturbations, negative controls, and robustness tests.

## Build and Test
Install dependencies with `pip install -r requirements.txt` (to be created).

Run Jupyter notebooks in sequence for analysis.

## Conventions
- Notebook naming: `{NN}_{descriptive_name}` (e.g., 01_pbmc3k_scanpy_basics)
- Implement scoring methods: average expression, Scanpy module score, rank-based, classifier-based
- Always include negative controls (random gene sets, permuted labels) and robustness tests (normalization variations, downsampling)
- Pin dependencies and seed randomness for reproducibility
- Document preprocessing parameters and gene set sources

See [README.md](README.md) for detailed project goals, biological rationale, dataset strategy, and workflow phases.