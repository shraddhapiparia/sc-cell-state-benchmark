"""Data loading functions."""

from pathlib import Path
from typing import List, Optional

import anndata as ad
import scanpy as sc


def ensure_dirs(paths: List[Path]) -> None:
    """Ensure that the given directories exist, creating them if necessary."""
    for path in paths:
        path.mkdir(parents=True, exist_ok=True)


def load_pbmc3k() -> ad.AnnData:
    """Load the PBMC3k dataset from Scanpy."""
    return sc.datasets.pbmc3k()


def save_anndata(adata: ad.AnnData, path: Path) -> None:
    """Save an AnnData object to disk in HDF5 format."""
    adata.write(path)


def load_anndata(path: Path) -> ad.AnnData:
    """Load an AnnData object from disk."""
    return ad.read_h5ad(path)


def load_kang_dataset(path: Optional[Path] = None) -> ad.AnnData:
    """Load the Kang PBMC interferon-beta stimulation dataset.

    If a local h5ad path is provided, it is loaded from disk. Otherwise, this
    function attempts to use Scanpy's built-in kang2018 loader.
    """
    if path is not None:
        return ad.read_h5ad(path)

    try:
        return sc.datasets.kang2018()
    except AttributeError as exc:
        raise RuntimeError(
            'Scanpy does not provide kang2018 dataset loader in this version. '
            'Please provide a local Kang h5ad file via --input /path/to/kang_pbmc_raw.h5ad.'
        ) from exc
