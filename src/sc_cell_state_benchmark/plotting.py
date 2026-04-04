"""Plotting functions."""

from pathlib import Path

import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns


def plot_umap_clusters(adata, save_path: Path) -> None:
    """Plot UMAP colored by Leiden clusters."""
    fig, ax = plt.subplots(figsize=(6, 5))
    sc.pl.umap(adata, color='leiden', ax=ax, show=False, legend_loc='on data')
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_umap_qc(adata, save_path: Path) -> None:
    """Plot UMAP colored by QC metrics (total counts and n_genes)."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    sc.pl.umap(adata, color='total_counts', ax=axes[0], show=False)
    sc.pl.umap(adata, color='n_genes_by_counts', ax=axes[1], show=False)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_umap_annotation(adata, save_path: Path) -> None:
    """Plot UMAP colored by predicted cell-type labels."""
    fig, ax = plt.subplots(figsize=(6, 5))
    sc.pl.umap(
        adata,
        color='predicted_cell_type',
        ax=ax,
        show=False,
        legend_loc='on data',
        title='Annotated cell types',
    )
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_umap_score(adata, score_field: str, save_path: Path, title: str) -> None:
    """Plot UMAP colored by a continuous score field."""
    fig, ax = plt.subplots(figsize=(6, 5))
    sc.pl.umap(adata, color=score_field, ax=ax, show=False, cmap='viridis', title=title)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_score_violin(scores_df, score_fields, group_field, save_path: Path, standardize: bool = True) -> None:
    """Plot violin distributions of scores grouped by cell type."""
    plot_df = scores_df.copy()
    if standardize:
        for field in score_fields:
            plot_df[field] = (plot_df[field] - plot_df[field].mean()) / plot_df[field].std(ddof=1)
    melted = plot_df.melt(id_vars=[group_field], value_vars=score_fields, var_name='method', value_name='score')
    plt.figure(figsize=(10, 6))
    sns.violinplot(data=melted, x=group_field, y='score', hue='method', split=False, inner='quartile')
    plt.xticks(rotation=45, ha='right')
    plt.title('Interferon scores by group')
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_score_by_cell_type_and_condition(
    scores_df,
    score_field: str,
    cell_type_col: str,
    condition_col: str,
    save_path: Path,
    title: str = None,
) -> None:
    """Violin plot of a score stratified by cell type and condition (ctrl vs stim)."""
    plot_df = scores_df[[cell_type_col, condition_col, score_field]].dropna()
    # Order cell types by median stimulated score descending
    stim_label = plot_df[condition_col].unique()
    order = (
        plot_df.groupby(cell_type_col, observed=True)[score_field]
        .median()
        .sort_values(ascending=False)
        .index.tolist()
    )
    n_types = len(order)
    fig_width = max(10, n_types * 1.4)
    fig, ax = plt.subplots(figsize=(fig_width, 5))
    sns.violinplot(
        data=plot_df,
        x=cell_type_col,
        y=score_field,
        hue=condition_col,
        order=order,
        split=False,
        inner='quartile',
        ax=ax,
    )
    ax.set_xlabel('Cell type')
    ax.set_ylabel(score_field)
    ax.set_title(title or f'{score_field} by cell type and condition')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_random_control_comparison(real_score, random_means, save_path: Path) -> None:
    """Plot the real score distribution versus random control distributions."""
    plt.figure(figsize=(8, 5))
    sns.histplot(random_means, bins=20, color='gray', kde=False, label='random control mean', stat='density')
    plt.axvline(real_score.mean(), color='red', linestyle='--', label='real gene set mean')
    plt.xlabel('Score')
    plt.ylabel('Density')
    plt.title('Real interferon score versus random control distributions')
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
