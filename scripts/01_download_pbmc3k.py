#!/usr/bin/env python3
"""Download PBMC3k dataset."""

from sc_cell_state_benchmark.config import DATA_RAW, PBMC3K_RAW
from sc_cell_state_benchmark.data import ensure_dirs, load_pbmc3k, save_anndata

if __name__ == "__main__":
    # Ensure data directory exists
    ensure_dirs([DATA_RAW])

    # Load dataset
    adata = load_pbmc3k()

    # Save to disk
    save_anndata(adata, PBMC3K_RAW)

    # Print summary
    print(f"[download] saved PBMC3k to {PBMC3K_RAW} | cells={adata.shape[0]} genes={adata.shape[1]}")
