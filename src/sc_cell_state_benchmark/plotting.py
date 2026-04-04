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


def plot_program_heatmap_delta(delta_df, save_path: Path, title: str = 'Immune program score change: stimulated \u2212 control') -> None:
    """Option B (primary): diverging heatmap of per-program mean score delta (stim − ctrl) by cell type.

    Parameters
    ----------
    delta_df : DataFrame, shape (programs, cell_types)
        Pre-computed stim_mean − ctrl_mean per program × cell type.
        Index = program display names; columns = cell type labels.
    """
    fig_w = max(8, delta_df.shape[1] * 1.1 + 2.5)
    fig_h = max(3.5, delta_df.shape[0] * 0.85 + 2.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    abs_max = delta_df.abs().max().max()
    if abs_max == 0 or (hasattr(abs_max, '__float__') and abs_max != abs_max):
        abs_max = 1.0  # fallback for all-zero or all-NaN

    sns.heatmap(
        delta_df,
        ax=ax,
        cmap='RdBu_r',
        center=0,
        vmin=-abs_max,
        vmax=abs_max,
        annot=True,
        fmt='.2f',
        linewidths=0.5,
        linecolor='#cccccc',
        cbar_kws={'label': 'stim \u2212 ctrl (Scanpy score)', 'shrink': 0.8},
    )
    ax.set_title(title)
    ax.set_xlabel('Cell type')
    ax.set_ylabel('')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_program_heatmap_sidebyside(
    ctrl_df,
    stim_df,
    save_path: Path,
    ctrl_label: str = 'Control',
    stim_label: str = 'Stimulated',
    title: str = 'Immune program scores by condition (Scanpy score)',
) -> None:
    """Option A (supporting): side-by-side heatmaps of raw program scores for each condition.

    Each panel has its own colour scale so that the within-condition pattern of which
    programs are high in which cell types is clearly visible. Use the delta heatmap
    (plot_program_heatmap_delta) to compare conditions directly.

    Parameters
    ----------
    ctrl_df, stim_df : DataFrame, shape (programs, cell_types)
        Mean program score per cell type for the control and stimulated conditions.
        Index = program display names; columns = cell type labels.
    """
    n_cts = ctrl_df.shape[1]
    n_progs = ctrl_df.shape[0]
    fig_w = max(16, n_cts * 1.1 * 2 + 3.0)
    fig_h = max(3.5, n_progs * 0.85 + 2.0)
    fig, axes = plt.subplots(1, 2, figsize=(fig_w, fig_h))

    for ax, df, label in zip(axes, [ctrl_df, stim_df], [ctrl_label, stim_label]):
        vmin = df.min().min()
        vmax = df.max().max()
        if vmin == vmax:
            vmax = vmin + 0.01
        sns.heatmap(
            df,
            ax=ax,
            cmap='YlOrRd',
            vmin=vmin,
            vmax=vmax,
            annot=True,
            fmt='.2f',
            linewidths=0.5,
            linecolor='#cccccc',
            cbar_kws={'label': 'Scanpy score', 'shrink': 0.8},
        )
        ax.set_title(label)
        ax.set_xlabel('Cell type')
        ax.set_ylabel('')
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right')

    axes[0].set_ylabel('Program')
    fig.suptitle(title, y=1.02)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_program_auc_barplot(auc_df, save_path: Path, title: str = 'Program separation of stimulated vs control (AUROC)') -> None:
    """Horizontal bar chart ranking programs by AUROC for ctrl vs stim separation.

    Parameters
    ----------
    auc_df : DataFrame
        Must contain columns: display_name, auc, n_genes_found, n_genes_total.
        Rows sorted descending by auc before plotting.
    """
    df = auc_df.sort_values('auc', ascending=True).reset_index(drop=True)
    fig_h = max(3.0, len(df) * 0.75 + 1.5)
    fig, ax = plt.subplots(figsize=(7.5, fig_h))

    bar_colors = ['#d73027' if v >= 0.7 else '#74add1' if v >= 0.5 else '#cccccc'
                  for v in df['auc']]
    bars = ax.barh(df['display_name'], df['auc'], color=bar_colors, edgecolor='white', height=0.6)

    ax.axvline(0.5, color='#333333', linestyle='--', linewidth=1.2, label='random (AUC = 0.5)')
    ax.set_xlim(0, 1.12)
    ax.set_xlabel('AUROC (stimulated vs control, all cell types)')
    ax.set_title(title)

    for bar, row in zip(bars, df.itertuples()):
        # AUC value to the right of the bar
        ax.text(
            row.auc + 0.015, bar.get_y() + bar.get_height() / 2,
            f'{row.auc:.3f}', va='center', ha='left', fontsize=9,
        )
        # Gene count inside the bar
        if row.auc > 0.15:
            ax.text(
                0.01, bar.get_y() + bar.get_height() / 2,
                f'n={int(row.n_genes_found)}/{int(row.n_genes_total)}',
                va='center', ha='left', fontsize=7.5, color='white',
            )

    ax.legend(loc='lower right', fontsize=9)
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
