"""Benchmarking and evaluation functions."""

import numpy as np
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score


def compare_scores(scores1, scores2):
    """Compare two scoring methods."""
    # TODO: Implement comparison metrics
    pass


def robustness_test(adata, perturbations):
    """Test score robustness under perturbations."""
    # TODO: Implement
    pass


def compute_binary_comparison(scores, labels, positive_label):
    """Compute simple benchmark metrics for a binary comparison."""
    labels = np.asarray(labels)
    positive = labels == positive_label
    if positive.sum() == 0 or positive.sum() == len(labels):
        raise ValueError('Binary comparison requires both positive and negative samples.')

    real_scores = np.asarray(scores)
    auc = roc_auc_score(positive.astype(int), real_scores)
    group1 = real_scores[positive]
    group0 = real_scores[~positive]
    stat, pvalue = mannwhitneyu(group1, group0, alternative='two-sided')
    effect_size = (group1.mean() - group0.mean()) / real_scores.std(ddof=1)
    return {
        'auc': float(auc),
        'mannwhitney_u': float(stat),
        'pvalue': float(pvalue),
        'effect_size': float(effect_size),
    }
