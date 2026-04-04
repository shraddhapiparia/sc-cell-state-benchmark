"""Cell-state scoring functions."""

from typing import List

import numpy as np
import scanpy as sc


def intersect_genes(adata, gene_list: List[str]) -> List[str]:
    """Return genes present in the AnnData object."""
    return [gene for gene in gene_list if gene in adata.var_names]


def average_expression_score(adata, gene_list: List[str]) -> np.ndarray:
    """Compute average expression of a gene list using normalized log data.
    
    Uses adata.raw if available (full gene matrix), otherwise falls back to adata.X.
    """
    present_genes = intersect_genes(adata, gene_list)
    if len(present_genes) == 0:
        return np.zeros(adata.n_obs)
    
    # Use full matrix from adata.raw if available, else use adata
    matrix = adata.raw if adata.raw is not None else adata
    expr = matrix[:, present_genes].X
    if hasattr(expr, 'toarray'):
        expr = expr.toarray()
    return np.asarray(expr).mean(axis=1).ravel()


def scanpy_score_genes(adata, gene_list: List[str], score_name: str = 'interferon_score') -> np.ndarray:
    """Compute Scanpy score_genes and return the resulting score array.
    
    Uses adata.raw if available (full gene matrix), otherwise uses main matrix.
    """
    present_genes = intersect_genes(adata, gene_list)
    if len(present_genes) == 0:
        adata.obs[score_name] = 0.0
        return adata.obs[score_name].to_numpy()
    # Use raw matrix if available, else use main matrix
    use_raw = adata.raw is not None
    sc.tl.score_genes(adata, present_genes, score_name=score_name, use_raw=use_raw)
    return adata.obs[score_name].to_numpy()


def rank_based_score(adata, gene_list: List[str]) -> np.ndarray:
    """Compute a simple rank-based score inspired by AUCell.
    
    Uses adata.raw if available (full gene matrix), otherwise falls back to adata.
    """
    present_genes = intersect_genes(adata, gene_list)
    if len(present_genes) == 0:
        return np.zeros(adata.n_obs)

    # Use full matrix from adata.raw if available, else use adata
    matrix_obj = adata.raw if adata.raw is not None else adata
    
    if matrix_obj.X.shape[1] == 0:
        return np.zeros(matrix_obj.n_obs)
    
    matrix = matrix_obj.X.toarray() if hasattr(matrix_obj.X, 'toarray') else matrix_obj.X
    ranks = np.argsort(np.argsort(-matrix, axis=1), axis=1)
    gene_idx = [matrix_obj.var_names.get_loc(g) for g in present_genes]
    selected_ranks = ranks[:, gene_idx]
    return np.mean(matrix.shape[1] - selected_ranks, axis=1)


def matched_random_gene_sets(adata, gene_list: List[str], n_sets: int = 25, seed: int = 42) -> List[List[str]]:
    """Generate matched random control gene sets with the same size.
    
    Uses adata.raw if available (full gene matrix), otherwise falls back to adata.
    """
    rng = np.random.default_rng(seed)
    # Use full matrix from adata.raw if available, else use adata
    matrix_obj = adata.raw if adata.raw is not None else adata
    gene_pool = [g for g in matrix_obj.var_names if g not in gene_list]
    size = len(gene_list)
    random_sets = []
    for _ in range(n_sets):
        random_sets.append(list(rng.choice(gene_pool, size, replace=False)))
    return random_sets

