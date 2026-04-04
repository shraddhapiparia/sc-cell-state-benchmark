# Project Summary: sc-cell-state-benchmark

## Objective

Benchmark gene-set scoring methods for single-cell RNA-seq perturbation data and
build a layered analysis from raw preprocessing through biological interpretation.

---

## Datasets

| Dataset | Role | Cells | Notes |
|---|---|---|---|
| PBMC3k (10x Genomics) | Warm-up / dry run | ~2,700 | Standard clustering pipeline; no perturbation signal |
| Kang PBMC (IFN-beta) | Main benchmark | 24,562 | Paired ctrl/stim; interferon-beta stimulation |

---

## Pipeline layers

**Layer 1 — Preprocessing and annotation (scripts 01–06)**
PBMC3k: QC filtering, normalisation, HVG selection, PCA, UMAP, Leiden clustering,
marker gene identification, cluster annotation. Used to validate the preprocessing
pipeline before applying it to the perturbation dataset.

**Layer 2 — Scoring method benchmark (scripts 07–09)**
Kang PBMC: three gene-set scoring methods (Scanpy `score_genes`, mean expression,
rank-based) compared head-to-head on a 10-gene interferon signature. 25 matched
random gene sets used as null controls. AUROC and Mann-Whitney U computed per method.

**Layer 3 — Biological program interpretation (script 10)**
Five curated immune programs scored across all cells using Scanpy `score_genes`.
Programs: IFN-α/β response, inflammatory response, antigen presentation (MHC-II),
cytotoxic/NK-like, monocyte activation. Per-cell-type × condition summaries and
AUROC ranking produced.

**Layer 4 — Exploratory cell-cell communication (script 11)**
16 curated ligand-receptor pairs scored using mean-product expression (mean ligand
expression in sender × mean receptor expression in receiver) per condition. Sender ×
receiver heatmaps (side-by-side and delta) and arc-network figure produced.

---

## Key findings

1. All three scoring methods separate stimulated from control cells with AUROC ≥ 0.985.
   The interferon signature achieves effect sizes around 1.7 SD; matched random gene
   sets produce near-zero effect sizes, confirming signal specificity.

2. Score scales differ across methods but group separation is consistent. Scanpy
   `score_genes` is the most background-corrected and was used for all downstream layers.

3. Program scoring shows IFN-α/β response (AUC 0.991) and inflammatory response
   (AUC 0.655) as the perturbation-responsive programs. Antigen presentation, cytotoxic,
   and monocyte activation programs behave primarily as cell-type identity signals —
   this is expected and informative for biological context.

4. Communication scoring identifies monocytes and dendritic cells as the dominant
   senders of increased signalling potential under IFN-beta stimulation, driven mainly
   by CXCL10 expression. IFNG and IL10 pairs score near zero, consistent with the
   absence of T-cell IFN-gamma secretion in this short in-vitro protocol.

---

## Limitations and caveats

- All results are from a single in-vitro IFN-beta stimulation dataset. Generalisation
  to other stimuli, time points, or disease contexts requires additional validation.
- Gene sets are manually curated minimal sets, not comprehensive pathway signatures.
- Cell-cell communication scores reflect ligand/receptor co-expression, not validated
  intercellular signalling.
- Megakaryocytes (n=22 total) are included in summaries but have very low cell counts
  and should not be interpreted with confidence.
- No permuted-label statistical null testing has been applied to either program or
  communication scores.

---

## Future extensions

- Additional scoring methods: ssGSEA, full AUCell, classifier-based upper bound
- Second perturbation dataset with a different stimulus for cross-dataset comparison
- Permuted-label null controls for program and communication scores
- MSigDB gene set loader
- Richer cell-cell communication statistics (permutation-based significance testing)
- Trajectory or pseudotime analysis for continuous state transitions
