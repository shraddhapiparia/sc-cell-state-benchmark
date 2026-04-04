#!/usr/bin/env python3
"""Preprocess a perturbation single-cell dataset for benchmarking."""

import argparse
import scanpy as sc
from pathlib import Path
from sc_cell_state_benchmark.config import DATA_PROCESSED, KANG_PBMC_PREPROCESSED, KANG_PBMC_RAW
from sc_cell_state_benchmark.data import ensure_dirs, load_anndata, save_anndata

METADATA_KEYWORDS = {
    'condition': ['stim', 'condition', 'treatment', 'group', 'status'],
    'cell_type': ['cell_type', 'celltype', 'celllabels', 'cell_labels', 'annotation'],
    'donor': ['donor', 'batch', 'patient', 'sample'],
}


def parse_args():
    parser = argparse.ArgumentParser(description='Preprocess a perturbation single-cell dataset.')
    parser.add_argument(
        '--input',
        type=Path,
        default=None,
        help='Local perturbation dataset h5ad file to preprocess. Defaults to the downloaded raw dataset.',
    )
    parser.add_argument(
        '--output',
        type=Path,
        default=KANG_PBMC_PREPROCESSED,
        help='Path to save the preprocessed AnnData file.',
    )
    return parser.parse_args()


def find_metadata_columns(obs):
    matches = {key: [] for key in METADATA_KEYWORDS}
    for col in obs.columns:
        low = col.lower()
        for key, keywords in METADATA_KEYWORDS.items():
            if any(keyword in low for keyword in keywords):
                matches[key].append(col)
    return matches


if __name__ == '__main__':
    args = parse_args()
    ensure_dirs([DATA_PROCESSED])

    raw_path = args.input if args.input is not None else KANG_PBMC_RAW
    adata = load_anndata(raw_path)
    metadata_matches = find_metadata_columns(adata.obs)

    print(f'[preprocess] loaded raw dataset from {raw_path}')
    print(f'[preprocess] available obs columns: {list(adata.obs.columns)}')
    print(
        f'[preprocess] likely metadata: condition={metadata_matches["condition"]}, '
        f'cell_type={metadata_matches["cell_type"]}, donor/batch={metadata_matches["donor"]}'
    )

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')
    
    # Store full normalized matrix for scoring
    adata.raw = adata.copy()
    
    # Subset to HVGs for dimensionality reduction
    adata = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, flavor='igraph', directed=False, n_iterations=2)

    save_anndata(adata, args.output)

    print(
        f"[preprocess] {raw_path.name} -> {args.output.name} | cells={adata.shape[0]} genes={adata.shape[1]} hvg={adata.var.highly_variable.sum()} clusters={adata.obs['leiden'].nunique()}"
    )