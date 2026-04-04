"""Negative control functions."""

import numpy as np


def random_gene_set(adata, size):
    """Generate random gene set of given size."""
    genes = adata.var_names
    return np.random.choice(genes, size, replace=False)


def permute_labels(labels):
    """Permute condition labels."""
    return np.random.permutation(labels)
