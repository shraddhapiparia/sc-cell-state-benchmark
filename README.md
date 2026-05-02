# sc-cell-state-benchmark

Benchmarking three RNA-based cell-state scoring methods (scanpy `score_genes`,
AUCell-style ranking, and a custom implementation) on the Kang PBMC IFN-β
dataset, with downstream biological interpretation via within-cell-type
differential expression, Hallmark pathway enrichment, and exploratory
ligand–receptor analysis.

RNA-only by design. Paired RNA+ATAC integration, TF/regulon inference, and
dynamic regulatory modeling live in a separate follow-up multiome repository.

---

## Main findings

- All three scoring methods separate stimulated from control cells with
  AUROC ≥ 0.985 — method choice matters less than is sometimes assumed on a
  strong perturbation signal.
- Interferon-stimulated genes (IFIT1, IFIT3, ISG15, MX1, IFI6) are induced
  across every immune cell type, with the largest response in CD14+ monocytes
  and dendritic cells.
- Hallmark Interferon Alpha and Interferon Gamma Response dominate enrichment
  in every cell type. Monocytes and DCs additionally show TNF-α / NF-κB and
  inflammatory-response enrichment, consistent with broader innate activation.

---

## Figures

<img src="figures/kang_umap_by_condition.png" width="500">

*Kang PBMC UMAP: strong ctrl vs IFN-β separation while preserving expected
immune cell types. Donor mixing is acceptable; condition is the dominant
source of variation.*

<img src="figures/kang_de_top5_per_cell_type.png" width="500">

*Top stim-vs-ctrl genes within each cell type. Classical ISGs are
consistently induced; monocytes and DCs show the strongest response.*

<img src="figures/kang_de_volcano_cd14_monocytes.png" width="500">

*CD14+ monocytes volcano. Broad transcriptional response dominated by
interferon and inflammatory genes.*

<img src="figures/kang_ora_up_top_pathways.png" width="500">

*Hallmark ORA on upregulated DEGs per cell type. IFN-α/γ response is
universal; TNF/NF-κB and inflammatory signaling are monocyte- and
DC-specific.*

---

## Pipeline

| Layer                   | Scripts | Purpose                                                  |
| ----------------------- | ------- | -------------------------------------------------------- |
| PBMC3k warm-up          | 01–06   | QC, clustering, annotation on a standard PBMC dataset    |
| Kang benchmark          | 07–09   | IFN scoring; comparison of three methods                 |
| Program scoring         | 10      | Curated immune programs across cell types and conditions |
| Cell-cell communication | 11      | Exploratory ligand–receptor co-expression                |
| Differential expression | 12      | Stim vs ctrl DE within each cell type                    |
| Pathway enrichment      | 13      | ORA (up/down) + preranked GSEA from DE genes             |

Numbered scripts, fixed seeds, modular `src/` package, single
`bash run_pipeline.sh` entrypoint.

---

## Datasets

| Dataset             | Role                  | Cells  | Source                                |
| ------------------- | --------------------- | ------ | ------------------------------------- |
| PBMC3k (10x)        | Pipeline validation   | ~2,700 | `scanpy.datasets.pbmc3k()`            |
| Kang PBMC IFN-β     | Main benchmark        | 24,562 | GEO GSE96583, expected at `data/raw/kang_pbmc_raw.h5ad` |

---

## Quick start

### Environment setup

```bash
# 1. Recreate the conda environment (Python + compiled dependencies)
conda env create -f environment.yml
conda activate sc-benchmark

# 2. Install the src/ package in editable mode so scripts can import sc_cell_state_benchmark modules
pip install -e .

# 3. Run the full pipeline
bash run_pipeline.sh
```

`environment.yml` pins the full conda environment for reproducibility.
`pyproject.toml` installs the `src/sc_cell_state_benchmark/` package so scripts can do
`from sc_cell_state_benchmark import scoring` without path hacks. Both files are intentional.

The pipeline expects the Kang PBMC dataset at `data/raw/kang_pbmc_raw.h5ad`
(GEO accession GSE96583). To use a file at a different path:

```bash
KANG_INPUT=/path/to/kang_pbmc_raw.h5ad bash run_pipeline.sh
```

Or run steps individually:

```bash
conda env create -f environment.yml
conda activate sc-benchmark
pip install -e .

python scripts/01_download_pbmc3k.py
python scripts/02_preprocess_pbmc3k.py
python scripts/03_plot_pbmc3k.py
python scripts/04_marker_genes.py
python scripts/05_annotate_clusters.py
python scripts/06_score_cell_states.py

python scripts/07_download_kang_pbmc.py
python scripts/08_preprocess_kang_pbmc.py
python scripts/09_score_kang_interferon.py --condition label --cell-type cell_type
python scripts/10_score_pathway_programs.py --condition label --cell-type cell_type
python scripts/11_cell_communication.py --condition label --cell-type cell_type
python scripts/12_de_per_cell_type.py --condition label --cell-type cell_type
python scripts/13_pathway_enrichment_per_cell_type.py
# canonical pathway script; 13_gsea_per_cell_type.py is deprecated
```

For the Kang dataset, the condition column is `label` (`ctrl`, `stim`) and the cell-type column is `cell_type`.

---

## Limitations

- Single in-vitro IFN-β stimulation dataset; AUROC ceiling effects likely
  mean the benchmark does not discriminate between methods on harder signals.
- Cell-cell communication is exploratory ligand–receptor co-expression, not
  validated signaling.
- Megakaryocytes are present at very low numbers and should not be
  interpreted strongly.
- Further caveats in `docs/IMPLEMENTATION_DECISIONS.md`.

---

## Roadmap

In this repository: additional scoring methods, permutation-based null
controls, subsampling robustness, more perturbation datasets, richer
gene-set loaders.

Follow-up multiome repository: paired RNA+ATAC integration, WNN-based
multimodal analysis, TF/regulon inference, dynamic regulatory modeling
and pseudotime.
