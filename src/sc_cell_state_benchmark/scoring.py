"""Cell-state scoring functions."""

from typing import List

import numpy as np
import scanpy as sc


def ucell_score(
    adata,
    gene_list: List[str],
    max_rank: int = 1500,
) -> np.ndarray:
    """Compute a UCell-style rank-statistic gene signature score per cell.

    This is inspired by UCell logic:
    - Rank genes within each cell by expression, highest expression = best rank.
    - Use the ranks of the signature genes.
    - Convert the rank-sum statistic into a normalized score.

    Higher score means the signature genes are ranked closer to the top of the
    cell's expression profile.

    Uses adata.raw if available, otherwise falls back to adata.
    """

    matrix_obj = adata.raw if adata.raw is not None else adata

    present_genes = [g for g in gene_list if g in matrix_obj.var_names]
    if len(present_genes) == 0:
        return np.zeros(adata.n_obs)

    if matrix_obj.X.shape[1] == 0:
        return np.zeros(matrix_obj.n_obs)

    matrix = matrix_obj.X.toarray() if hasattr(matrix_obj.X, "toarray") else matrix_obj.X

    n_cells, n_genes = matrix.shape
    max_rank = min(max_rank, n_genes)

    gene_idx = np.array([matrix_obj.var_names.get_loc(g) for g in present_genes])
    n_signature_genes = len(gene_idx)

    best_rank_sum = n_signature_genes * (n_signature_genes + 1) / 2
    worst_rank_sum = n_signature_genes * max_rank

    if worst_rank_sum == best_rank_sum:
        return np.zeros(n_cells)

    scores = np.zeros(n_cells)

    for i in range(n_cells):
        # Highest expression gets rank 1.
        order = np.argsort(-matrix[i, :])
        ranks = np.empty(n_genes, dtype=float)
        ranks[order] = np.arange(1, n_genes + 1)

        # UCell-style truncation: genes below max_rank are treated as max_rank.
        ranks = np.minimum(ranks, max_rank)

        signature_rank_sum = ranks[gene_idx].sum()

        score = 1 - (
            (signature_rank_sum - best_rank_sum)
            / (worst_rank_sum - best_rank_sum)
        )

        scores[i] = max(score, 0.0)

    return scores

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

# New AUC score computation
def aucell_score(
    adata,
    gene_list: List[str],
    auc_max_rank: int | None = None,
) -> np.ndarray:
    """Compute AUCell-style gene-set enrichment score per cell.

    For each cell:
    1. Rank genes from highest to lowest expression.
    2. Check where the signature genes appear in the ranked list.
    3. Compute AUC of the recovery curve within the top auc_max_rank genes.
    """

    matrix_obj = adata.raw if adata.raw is not None else adata

    present_genes = [g for g in gene_list if g in matrix_obj.var_names]
    if len(present_genes) == 0:
        return np.zeros(matrix_obj.n_obs)

    matrix = matrix_obj.X.toarray() if hasattr(matrix_obj.X, "toarray") else matrix_obj.X
    n_cells, n_genes = matrix.shape

    if auc_max_rank is None:
        auc_max_rank = int(0.05 * n_genes)  # top 5% genes, AUCell-like default idea

    auc_max_rank = min(auc_max_rank, n_genes)

    gene_idx = np.array([matrix_obj.var_names.get_loc(g) for g in present_genes])
    gene_idx_set = set(gene_idx)

    scores = np.zeros(n_cells)

    for i in range(n_cells):
        ranked_genes = np.argsort(-matrix[i, :])[:auc_max_rank]

        hits = np.array([1 if g in gene_idx_set else 0 for g in ranked_genes])

        if hits.sum() == 0:
            scores[i] = 0.0
            continue

        recovery_curve = np.cumsum(hits)

        auc = np.trapz(recovery_curve, dx=1)

        max_auc = len(present_genes) * auc_max_rank
        scores[i] = auc / max_auc

    return scores

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

