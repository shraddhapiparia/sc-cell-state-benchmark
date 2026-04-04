#!/usr/bin/env python3
"""Download or load a perturbation single-cell dataset for benchmarking."""

import argparse
from pathlib import Path

from sc_cell_state_benchmark.config import DATA_RAW, KANG_PBMC_RAW
from sc_cell_state_benchmark.data import ensure_dirs, load_kang_dataset, save_anndata


def parse_args():
    parser = argparse.ArgumentParser(description='Download or load a perturbation single-cell dataset.')
    parser.add_argument(
        '--input',
        type=Path,
        default=None,
        help='Optional local perturbation dataset h5ad file to load instead of Scanpy loader.',
    )
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    ensure_dirs([DATA_RAW])
    adata = load_kang_dataset(path=args.input)
    save_anndata(adata, KANG_PBMC_RAW)
    print(f"[download] saved perturbation dataset to {KANG_PBMC_RAW} | cells={adata.shape[0]} genes={adata.shape[1]}")