#!/usr/bin/env bash
# Run the full sc-cell-state-benchmark pipeline in order.
#
# Prerequisites:
#   conda env create -f environment.yml
#   conda activate sc-benchmark
#
# Steps 01 and 07 require internet access (dataset downloads).
#   Step 01 downloads PBMC3k via scanpy.
#   Step 07 downloads the Kang 2018 IFN-beta PBMC dataset (~37 MB) from the
#   scverse example data CDN (exampledata.scverse.org).  If the download fails,
#   fetch manually and pass via --input:
#     wget -O data/raw/kang_pbmc_raw.h5ad https://exampledata.scverse.org/pertpy/kang_2018.h5ad
#     python scripts/07_download_kang_pbmc.py --input data/raw/kang_pbmc_raw.h5ad
# Step 13 requires internet access (MSigDB Hallmark gene set download via gseapy).
#
# Usage:
#   bash run_pipeline.sh

set -euo pipefail

# ── PBMC3k warm-up (scripts 01–06) ───────────────────────────────────────────
python scripts/01_download_pbmc3k.py
python scripts/02_preprocess_pbmc3k.py
python scripts/03_plot_pbmc3k.py
python scripts/04_marker_genes.py
python scripts/05_annotate_clusters.py
python scripts/06_score_cell_states.py

# ── Kang IFN-β benchmark (scripts 07–13) ─────────────────────────────────────
KANG_INPUT="${KANG_INPUT:-data/raw/kang_pbmc_raw.h5ad}"

if [[ -f "$KANG_INPUT" ]]; then
  python scripts/07_download_kang_pbmc.py --input "$KANG_INPUT"
else
  echo "Kang dataset not found at $KANG_INPUT"
  echo "Please place the file at data/raw/kang_pbmc_raw.h5ad or set KANG_INPUT=/path/to/file"
  exit 1
fi

python scripts/08_preprocess_kang_pbmc.py

# For the Kang dataset: condition column = "label" (ctrl/stim),
#                       cell-type column = "cell_type"
python scripts/09_score_kang_interferon.py --condition label --cell-type cell_type
python scripts/10_score_pathway_programs.py --condition label --cell-type cell_type
python scripts/11_cell_communication.py --condition label --cell-type cell_type
python scripts/12_de_per_cell_type.py --condition label --cell-type cell_type

# Step 13 runs preranked GSEA with 1000 permutations per cell type (~5–10 min).
python scripts/13_pathway_enrichment_per_cell_type.py
