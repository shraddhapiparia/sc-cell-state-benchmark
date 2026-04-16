#!/usr/bin/env python3
"""Preprocess a perturbation single-cell dataset for benchmarking.

HVG selection
-------------
Uses flavor='seurat_v3' with batch_key set to the donor/replicate column.
seurat_v3 requires raw integer counts in adata.X, so HVG selection is performed
BEFORE normalize_total.  The Kang dataset has a paired design (each of 8 donors
has both ctrl and stim cells), so using the donor column as batch_key selects
genes that are variable within donors -- reducing the chance that ISGs dominate
the HVG list purely because of the ctrl/stim split.

If no donor column is detected, the script falls back to batch_key set to the
condition column.  If neither is available, HVG selection proceeds without a
batch_key and a warning is printed.

QC thresholds
-------------
The Kang raw dataset has been pre-filtered by the original Seurat workflow
(nCount_RNA / nFeature_RNA are already in obs).  The n_genes distribution is
tight (median ~519, 99th percentile ~1267, max 2757; only 2 cells exceed 2500).
MT genes are absent from this dataset's gene vocabulary, so MT filtering is
skipped automatically.

  MAX_GENES = 2500  -- upper bound on n_genes_by_counts (doublet proxy)

Inspect figures/kang_qc_violin_prefilter.png before changing thresholds.
"""

import argparse
import scanpy as sc
from pathlib import Path

from sc_cell_state_benchmark.config import (
    DATA_PROCESSED,
    FIGURES,
    KANG_PBMC_PREPROCESSED,
    KANG_PBMC_RAW,
    RANDOM_SEED,
)
from sc_cell_state_benchmark.data import ensure_dirs, load_anndata, save_anndata
from sc_cell_state_benchmark.plotting import plot_qc_violin, plot_umap_categorical

# QC filter threshold -- inspect kang_qc_violin_prefilter.png before changing
MAX_GENES: int = 2500

METADATA_KEYWORDS = {
    "condition": ["stim", "condition", "treatment", "group", "status"],
    "cell_type": ["cell_type", "celltype", "celllabels", "cell_labels", "annotation"],
    # 'replicate' added: Kang dataset uses that column name for donor/patient ID
    "donor": ["donor", "batch", "patient", "sample", "replicate"],
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Preprocess a perturbation single-cell dataset."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help="Local h5ad file to preprocess. Defaults to the downloaded raw dataset.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=KANG_PBMC_PREPROCESSED,
        help="Path to save the preprocessed AnnData file.",
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


if __name__ == "__main__":
    args = parse_args()
    ensure_dirs([DATA_PROCESSED, FIGURES])

    raw_path = args.input if args.input is not None else KANG_PBMC_RAW
    adata = load_anndata(raw_path)
    metadata_matches = find_metadata_columns(adata.obs)

    print(f"[preprocess] loaded raw dataset from {raw_path}")
    print(f"[preprocess] cells={adata.n_obs} genes={adata.n_vars}")
    print(f"[preprocess] available obs columns: {list(adata.obs.columns)}")
    print(
        f"[preprocess] detected metadata -- condition={metadata_matches['condition']}, "
        f"cell_type={metadata_matches['cell_type']}, donor={metadata_matches['donor']}"
    )

    # Identify batch key for HVG selection (donor preferred over condition)
    if metadata_matches["donor"]:
        hvg_batch_key = metadata_matches["donor"][0]
        print(f"[preprocess] HVG batch_key: {hvg_batch_key!r} (donor/replicate)")
    elif metadata_matches["condition"]:
        hvg_batch_key = metadata_matches["condition"][0]
        print(
            f"[preprocess] WARNING: no donor column found. "
            f"Using condition column {hvg_batch_key!r} as HVG batch_key."
        )
    else:
        hvg_batch_key = None
        print(
            "[preprocess] WARNING: no donor or condition column found. "
            "HVG selection will proceed without a batch_key."
        )

    # --- mitochondrial gene annotation ----------------------------------
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    n_mt = int(adata.var["mt"].sum())
    print(f"[preprocess] {n_mt} mitochondrial genes annotated")
    if n_mt == 0:
        print(
            "[preprocess] NOTE: no MT- genes found in this dataset. "
            "MT filtering will be skipped."
        )

    # --- basic cell and gene filters ------------------------------------
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"[preprocess] after min-gene/min-cell filter: cells={adata.n_obs}")

    # --- QC metrics -----------------------------------------------------
    qc_vars = ["mt"] if n_mt > 0 else []
    sc.pp.calculate_qc_metrics(adata, qc_vars=qc_vars, inplace=True)

    # --- pre-filter QC violin -------------------------------------------
    # Inspect this figure to verify the threshold choices above are appropriate.
    qc_fig_path = FIGURES / "kang_qc_violin_prefilter.png"
    qc_thresholds = {"n_genes_by_counts": MAX_GENES}
    if n_mt > 0:
        qc_thresholds["pct_counts_mt"] = 5.0  # placeholder; update from the figure
    plot_qc_violin(adata, qc_fig_path, thresholds=qc_thresholds)
    print(f"[preprocess] saved pre-filter QC figure to {qc_fig_path}")

    # --- upper-bound filter (MT filter skipped if no MT genes) ----------
    n_before = adata.n_obs
    sc.pp.filter_cells(adata, max_genes=MAX_GENES)
    n_after_genes = adata.n_obs
    print(
        f"[preprocess] upper-bound filter (>{MAX_GENES} genes): "
        f"removed {n_before - n_after_genes} cells | retained {adata.n_obs}"
    )
    if n_mt > 0:
        adata = adata[adata.obs["pct_counts_mt"] < 5.0].copy()
        print(f"[preprocess] MT filter: retained {adata.n_obs} cells")

    # --- HVG selection on raw counts ------------------------------------
    # seurat_v3 requires raw integer counts in adata.X.  This must run
    # BEFORE normalize_total.
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=2000,
        flavor="seurat_v3",
        batch_key=hvg_batch_key,
    )
    n_hvg = int(adata.var["highly_variable"].sum())
    print(
        f"[preprocess] seurat_v3 HVG selection (batch_key={hvg_batch_key!r}): "
        f"{n_hvg} HVGs selected"
    )

    # --- normalisation and log transform --------------------------------
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Store the full log-normalised matrix before HVG subsetting.
    # Scoring functions (scripts 09-11) use adata.raw to access the
    # complete gene vocabulary; this must be set after normalisation.
    adata.raw = adata.copy()
    print(f"[preprocess] adata.raw set: {adata.raw.n_vars} genes stored")

    # --- dimensionality reduction and clustering ------------------------
    adata = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, random_state=RANDOM_SEED)
    sc.pp.neighbors(adata, random_state=RANDOM_SEED)
    sc.tl.umap(adata, random_state=RANDOM_SEED)
    sc.tl.leiden(
        adata, flavor="igraph", directed=False, n_iterations=2, random_state=RANDOM_SEED
    )

    save_anndata(adata, args.output)

    n_clusters = adata.obs["leiden"].nunique()
    print(
        f"[preprocess] {raw_path.name} -> {args.output.name} | "
        f"cells={adata.n_obs} genes={adata.n_vars} hvg={n_hvg} clusters={n_clusters}"
    )

    # UMAP coloured by condition (label) and by donor (replicate).
    # These figures are useful for visually checking whether the embedding
    # separates cell types (expected) or is dominated by condition/donor
    # effects (would suggest batch correction is needed).
    if "label" in adata.obs.columns:
        plot_umap_categorical(
            adata,
            color="label",
            save_path=FIGURES / "kang_umap_by_condition.png",
            title="Kang PBMC -- UMAP by condition (ctrl / stim)",
        )
        print("[preprocess] saved kang_umap_by_condition.png")

    if "replicate" in adata.obs.columns:
        plot_umap_categorical(
            adata,
            color="replicate",
            save_path=FIGURES / "kang_umap_by_replicate.png",
            title="Kang PBMC -- UMAP by donor (replicate)",
        )
        print("[preprocess] saved kang_umap_by_replicate.png")
