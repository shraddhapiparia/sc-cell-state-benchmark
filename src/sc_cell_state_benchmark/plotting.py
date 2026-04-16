"""Plotting functions."""

from pathlib import Path
from typing import List

import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
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


def plot_qc_violin(adata, save_path: Path, thresholds: dict = None) -> None:
    """Plot per-cell QC metric distributions before filtering.

    Produces one violin per metric: n_genes_by_counts, total_counts, and
    pct_counts_mt (only when MT genes were found and the column is non-zero).
    Optional threshold lines are drawn as red dashed horizontal rules so the
    figure documents the chosen filter values.

    Parameters
    ----------
    adata : AnnData
        Object with QC metrics already computed via sc.pp.calculate_qc_metrics.
    save_path : Path
        Destination file for the saved figure.
    thresholds : dict, optional
        Mapping of obs column name to threshold value, e.g.
        {'n_genes_by_counts': 2500, 'pct_counts_mt': 5}.
        Each key must match a column in adata.obs.
    """
    thresholds = thresholds or {}

    # Decide which panels to show
    metrics = ['n_genes_by_counts', 'total_counts']
    labels = ['Genes per cell', 'Total counts per cell']
    if 'pct_counts_mt' in adata.obs.columns and adata.obs['pct_counts_mt'].max() > 0:
        metrics.append('pct_counts_mt')
        labels.append('% mitochondrial counts')

    n_panels = len(metrics)
    fig, axes = plt.subplots(1, n_panels, figsize=(4 * n_panels, 5))
    if n_panels == 1:
        axes = [axes]

    for ax, metric, label in zip(axes, metrics, labels):
        data = adata.obs[metric].values
        ax.violinplot(data, positions=[0], showmedians=True, showextrema=True)
        ax.set_xticks([])
        ax.set_ylabel(label)
        ax.set_title(label)

        if metric in thresholds:
            ax.axhline(thresholds[metric], color='red', linestyle='--', linewidth=1.2,
                       label=f'threshold = {thresholds[metric]}')
            ax.legend(fontsize=8)

    fig.suptitle('Pre-filter QC metrics (one dot = one cell)', y=1.01, fontsize=10)
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


def plot_umap_categorical(adata, color: str, save_path: Path, title: str = None) -> None:
    """Plot UMAP colored by any categorical obs column."""
    fig, ax = plt.subplots(figsize=(7, 5))
    sc.pl.umap(
        adata,
        color=color,
        ax=ax,
        show=False,
        title=title or color,
    )
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_marker_heatmap(
    adata,
    genes: List[str],
    groupby: str,
    save_path: Path,
    title: str = 'Marker gene expression by cell type',
) -> None:
    """Heatmap of mean log-normalised expression for canonical marker genes.

    Uses adata.raw when available so genes not in the HVG subset are still
    accessible.  Values are z-scored per gene (column-wise) so that differences
    across cell types are visible regardless of absolute expression level.

    Parameters
    ----------
    adata : AnnData
    genes : list of str
        Genes to show as columns.  Genes absent from the dataset are silently
        dropped with a printed warning.
    groupby : str
        obs column to group rows by (e.g. 'predicted_cell_type').
    save_path : Path
    title : str
    """
    import pandas as pd

    matrix = adata.raw if adata.raw is not None else adata
    var_names = list(matrix.var_names)

    present = [g for g in genes if g in var_names]
    missing = [g for g in genes if g not in var_names]
    if missing:
        print(f'[plot_marker_heatmap] genes not found, skipped: {missing}')
    if not present:
        print('[plot_marker_heatmap] no genes found -- skipping figure')
        return

    import scipy.sparse as sp

    X = matrix[:, present].X
    if sp.issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=float)

    groups = adata.obs[groupby].astype(str)
    unique_groups = sorted(groups.unique())

    mean_expr = np.zeros((len(unique_groups), len(present)))
    for i, grp in enumerate(unique_groups):
        mask = groups == grp
        mean_expr[i] = X[mask].mean(axis=0)

    df = pd.DataFrame(mean_expr, index=unique_groups, columns=present)

    # z-score per gene so scale differences do not dominate
    col_std = df.std(axis=0).replace(0, 1)
    df_z = (df - df.mean(axis=0)) / col_std

    fig_w = max(6, len(present) * 0.75 + 2.0)
    fig_h = max(3, len(unique_groups) * 0.55 + 1.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(
        df_z,
        ax=ax,
        cmap='RdBu_r',
        center=0,
        annot=True,
        fmt='.1f',
        linewidths=0.4,
        linecolor='#dddddd',
        cbar_kws={'label': 'z-scored mean log-expr', 'shrink': 0.7},
    )
    ax.set_title(title)
    ax.set_xlabel('Gene')
    ax.set_ylabel('Cell type')
    plt.xticks(rotation=45, ha='right')
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
    # Order cell types by median score descending
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


def plot_communication_heatmap_sidebyside(
    ctrl_matrix,
    stim_matrix,
    save_path: Path,
    ctrl_label: str = 'Control',
    stim_label: str = 'Stimulated',
    title: str = 'Cell-cell communication scores by condition (summed LR mean-product)',
) -> None:
    """Option A (supporting): side-by-side heatmaps of summed LR scores per condition.

    Rows = sender cell types; columns = receiver cell types. Each cell is the sum
    of mean-product LR scores across all curated pairs for that sender→receiver
    combination. Each panel uses its own colour scale so within-condition structure
    is visible independently.

    Parameters
    ----------
    ctrl_matrix, stim_matrix : DataFrame, shape (n_cell_types, n_cell_types)
        Sender × receiver summed LR score matrices, one per condition.
    """
    n_ct = ctrl_matrix.shape[0]
    fig_panel_w = max(6.0, n_ct * 0.85 + 2.0)
    fig_h = max(5.0, n_ct * 0.75 + 2.5)
    fig, axes = plt.subplots(1, 2, figsize=(fig_panel_w * 2 + 1.0, fig_h))

    for ax, df, label in zip(axes, [ctrl_matrix, stim_matrix], [ctrl_label, stim_label]):
        vmin = float(df.min().min())
        vmax = float(df.max().max())
        if vmax <= vmin:
            vmax = vmin + 0.01
        sns.heatmap(
            df,
            ax=ax,
            cmap='YlOrRd',
            vmin=vmin,
            vmax=vmax,
            annot=True,
            fmt='.2f',
            linewidths=0.4,
            linecolor='#dddddd',
            cbar_kws={'label': 'Summed LR score', 'shrink': 0.8},
        )
        ax.set_title(label, fontsize=12)
        ax.set_xlabel('Receiver cell type')
        ax.set_ylabel('Sender cell type')
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right', fontsize=8)
        plt.setp(ax.get_yticklabels(), rotation=0, fontsize=8)

    fig.suptitle(title, y=1.01, fontsize=11)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_communication_heatmap_delta(
    delta_matrix,
    save_path: Path,
    title: str = 'Cell-cell communication change: stimulated \u2212 control (summed LR score)',
) -> None:
    """Option B (primary): diverging heatmap of summed LR score change (stim \u2212 ctrl).

    Rows = sender cell types; columns = receiver cell types. Positive values (red)
    indicate increased communication potential in the stimulated condition;
    negative values (blue) indicate decreased potential.

    Parameters
    ----------
    delta_matrix : DataFrame, shape (n_cell_types, n_cell_types)
        stim_matrix \u2212 ctrl_matrix.
    """
    n_ct = delta_matrix.shape[0]
    fig_w = max(7.0, n_ct * 0.9 + 2.5)
    fig_h = max(5.0, n_ct * 0.8 + 2.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    abs_max = float(delta_matrix.abs().max().max())
    if abs_max == 0:
        abs_max = 0.01

    sns.heatmap(
        delta_matrix,
        ax=ax,
        cmap='RdBu_r',
        center=0,
        vmin=-abs_max,
        vmax=abs_max,
        annot=True,
        fmt='.2f',
        linewidths=0.4,
        linecolor='#dddddd',
        cbar_kws={'label': 'stim \u2212 ctrl (summed LR score)', 'shrink': 0.8},
    )
    ax.set_title(title, fontsize=11)
    ax.set_xlabel('Receiver cell type')
    ax.set_ylabel('Sender cell type')
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right', fontsize=8)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=8)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_communication_network(
    delta_matrix,
    save_path: Path,
    top_n: int = 14,
    title: str = 'Increased cell-cell communication: stimulated vs control',
) -> None:
    """Circular arc-network diagram of communication changes (stim \u2212 ctrl).

    Nodes represent cell types, positioned on a circle. Directed arcs go from
    sender to receiver. Only arcs with a positive delta (increased communication
    under stimulation) are drawn. Arc width and opacity are proportional to the
    delta magnitude. The top_n arcs by delta are shown.

    No external dependencies beyond matplotlib and numpy.

    Parameters
    ----------
    delta_matrix : DataFrame, shape (n_cell_types, n_cell_types)
        stim \u2212 ctrl summed LR score per sender × receiver.
    top_n : int
        Maximum number of arcs to draw.
    """
    cell_types: List[str] = list(delta_matrix.index)
    n = len(cell_types)

    # Short labels for readability inside the figure
    short: dict = {}
    for ct in cell_types:
        parts = ct.split()
        if len(parts) == 1:
            short[ct] = ct
        elif parts[0].startswith('CD') or parts[0].startswith('FCGR') or parts[0].startswith('NK'):
            short[ct] = parts[0]
        else:
            short[ct] = ' '.join(p[0] for p in parts)  # initials fallback
    # Override ambiguous initials
    short['CD14+ Monocytes'] = 'CD14+ Mono'
    short['FCGR3A+ Monocytes'] = 'FCGR3A+ Mono'
    short['Dendritic cells'] = 'DC'
    short['Megakaryocytes'] = 'Mega'

    # Node positions on unit circle (start from top, clockwise)
    angles = np.linspace(np.pi / 2, np.pi / 2 - 2 * np.pi, n, endpoint=False)
    pos = {ct: np.array([np.cos(a), np.sin(a)]) for ct, a in zip(cell_types, angles)}

    # Palette for nodes — consistent with seaborn tab10
    palette = plt.cm.tab10(np.linspace(0, 1, n))
    node_color = {ct: palette[i] for i, ct in enumerate(cell_types)}

    # Collect positive-delta sender→receiver pairs, sorted descending
    edges = []
    for sender in cell_types:
        for receiver in cell_types:
            if sender == receiver:
                continue
            val = float(delta_matrix.loc[sender, receiver])
            if val > 0:
                edges.append((sender, receiver, val))
    edges.sort(key=lambda x: x[2], reverse=True)
    edges = edges[:top_n]

    if not edges:
        # Nothing to draw — save blank figure with message
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.text(0.5, 0.5, 'No positive-delta interactions to display.',
                ha='center', va='center', transform=ax.transAxes)
        ax.axis('off')
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close()
        return

    max_delta = edges[0][2]

    fig, ax = plt.subplots(figsize=(9, 9))
    ax.set_aspect('equal')
    ax.set_xlim(-1.75, 1.75)
    ax.set_ylim(-1.75, 1.75)
    ax.axis('off')

    # Draw arcs using FancyArrowPatch with arc3 connectionstyle
    # Alternate arc curvature sign so overlapping pairs are distinguishable
    for i, (sender, receiver, delta_val) in enumerate(edges):
        norm = delta_val / max_delta           # 0–1
        lw = 0.8 + norm * 5.5
        alpha = 0.35 + norm * 0.55
        rad = 0.28 if i % 2 == 0 else -0.28   # alternate curvature
        color = plt.cm.YlOrRd(0.3 + norm * 0.65)

        arrow = mpatches.FancyArrowPatch(
            posA=tuple(pos[sender] * 0.88),
            posB=tuple(pos[receiver] * 0.88),
            arrowstyle=mpatches.ArrowStyle('->', head_length=8, head_width=5),
            connectionstyle=f'arc3,rad={rad}',
            linewidth=lw,
            color=color,
            alpha=alpha,
            zorder=2,
        )
        ax.add_patch(arrow)

    # Draw nodes on top of arcs
    node_r = 0.13
    for ct in cell_types:
        x, y = pos[ct]
        circle = plt.Circle((x, y), node_r, color=node_color[ct],
                             zorder=4, linewidth=1.5, ec='white')
        ax.add_patch(circle)

    # Labels just outside the node ring
    label_r = 1.32
    for ct in cell_types:
        x, y = pos[ct]
        ang = np.arctan2(y, x)
        lx, ly = label_r * np.cos(ang), label_r * np.sin(ang)
        ha = 'left' if np.cos(ang) >= 0 else 'right'
        ax.text(lx, ly, short[ct], ha=ha, va='center', fontsize=9,
                path_effects=[pe.withStroke(linewidth=2, foreground='white')])

    # Colour-scale legend (proxy patches)
    legend_handles = [
        mpatches.Patch(color=plt.cm.YlOrRd(0.95), label=f'High delta (~{max_delta:.1f})'),
        mpatches.Patch(color=plt.cm.YlOrRd(0.50), label='Medium delta'),
        mpatches.Patch(color=plt.cm.YlOrRd(0.30), label='Low delta'),
    ]
    ax.legend(handles=legend_handles, loc='lower right', fontsize=8,
              title='Summed LR\nscore increase', title_fontsize=8, framealpha=0.85)

    ax.set_title(title, fontsize=11, pad=14)
    fig.text(0.5, 0.01,
             f'Top {len(edges)} sender\u2192receiver pairs by stim\u2212ctrl delta. '
             'Arc width \u221d magnitude.',
             ha='center', fontsize=7.5, color='#555555')
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
