"""Configuration file for paths and constants."""

from pathlib import Path

# Project root
PROJECT_ROOT = Path(__file__).parent.parent.parent

# Data paths
DATA_RAW = PROJECT_ROOT / "data" / "raw"
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"

# Results paths
RESULTS_TABLES = PROJECT_ROOT / "results" / "tables"
RESULTS_METRICS = PROJECT_ROOT / "results" / "metrics"

# Figures path
FIGURES = PROJECT_ROOT / "figures"

# File path constants
PBMC3K_RAW = DATA_RAW / "pbmc3k_raw.h5ad"
PBMC3K_PREPROCESSED = DATA_PROCESSED / "pbmc3k_preprocessed.h5ad"
PBMC3K_ANNOTATED = DATA_PROCESSED / "pbmc3k_annotated.h5ad"
KANG_PBMC_RAW = DATA_RAW / "kang_pbmc_raw.h5ad"
KANG_PBMC_PREPROCESSED = DATA_PROCESSED / "kang_pbmc_preprocessed.h5ad"

# Constants
RANDOM_SEED = 42
