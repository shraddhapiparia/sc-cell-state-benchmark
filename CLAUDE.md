# CLAUDE.md

This file provides guidance to Claude Code when working in this repository.

## Project Overview

This repository benchmarks RNA-based cell-state scoring methods on an IFN-β stimulated PBMC dataset (Kang et al., GSE96583). It is RNA-only. Paired RNA+ATAC multiome analysis is developed in a separate follow-up repository.

The project is organized as numbered scripts (01-13). Each script should:

* read outputs from earlier scripts
* write clearly named files to `results/` or `figures/`
* avoid modifying previous outputs in-place
* save both tables and figures when possible

## Environment and Package Setup

Canonical reproduction path:

```bash
conda env create -f environment.yml   # recreate pinned conda environment
conda activate sc-benchmark
pip install -e .                      # install src/ package in editable mode
bash run_pipeline.sh
```

`environment.yml` -- full pinned conda environment for exact reproduction. Generated with `conda env export`. To regenerate portably (no build hashes, no prefix): `conda env export --no-builds | grep -v "^prefix:" > environment.yml`.

`pyproject.toml` -- declares the `sc_cell_state_benchmark` package and its direct dependencies. Needed so scripts can `from sc_cell_state_benchmark import scoring` without `sys.path` hacks. Both files are intentional.

## Pipeline Scripts

| Scripts | Layer | Purpose |
|---------|-------|---------|
| 01-06 | PBMC3k warm-up | QC, clustering, annotation, cell-state scoring |
| 07-09 | Kang benchmark | IFN-β scoring, method comparison |
| 10 | Program scoring | Curated immune programs across conditions |
| 11 | Cell communication | Exploratory ligand-receptor co-expression |
| 12 | Differential expression | Stim vs ctrl DE per cell type |
| 13 | Pathway enrichment | ORA and preranked GSEA on DE genes |

Key outputs:

* `figures/kang_umap_by_condition.png`
* `results/kang_method_comparison.csv`
* `results/kang_gsea_preranked_results.csv`

## Do Not Commit

* `data/raw/` -- large input files (Kang h5ad is ~several hundred MB)
* `data/processed/` -- generated intermediate AnnData objects
* `results/` and `figures/` -- generated outputs (checked in selectively for portfolio)
* `.env` or any credentials

## Coding Style

* Prefer short, modular scripts
* Keep reusable functions inside `src/`
* Avoid hidden side effects
* Every script should print:

  * what file was loaded
  * number of cells/features retained
  * what output file was written
* Use deterministic seeds when possible
* Save intermediate CSVs for debugging

## Figure Expectations

Every major phase should produce at least:

* one summary table
* one publication-style figure
* one markdown summary file

Preferred figure types:

* annotated UMAPs
* violin plots by cell type
* heatmaps
* TF activity barplots
* RNA vs ATAC scatter plots
* pseudotime trend curves

## Common Pitfalls

* Seurat assay accidentally switched from RNA to ATAC
* Peak names not matching genome build
* Slingshot failing because only one cluster is present
* ChromVAR motif names not matching JASPAR motif IDs
* Broad labels missing after subsetting
* Row mismatch when replacing ATAC count matrices

When debugging, first verify:

1. assay names
2. number of cells
3. number of peaks/features
4. metadata columns
5. rownames consistency between matrices and Seurat object

## Pathway analysis best practices

When adding or modifying pathway-enrichment code:

1. Distinguish ORA from GSEA clearly in script names, docstrings, outputs, and plots.
2. For GSEA/prerank:
   - use the full tested gene list
   - rank genes by a signed statistic when available
   - do not prefilter to significant genes before GSEA
3. For ORA:
   - analyze upregulated and downregulated genes separately
   - use an explicit experiment-specific background/universe of tested or expressed genes
4. Prefer low-redundancy libraries for compact summaries, but document that this is not equivalent to formal redundancy collapse.
5. Before implementing pathway analysis changes, check current best-practice documentation from authoritative sources and record the links and key decisions in `docs/ANALYSIS_DECISIONS.md`.
6. In `docs/ANALYSIS_DECISIONS.md`, include:
   - method used (ORA vs GSEA)
   - ranking metric
   - thresholds used
   - background definition
   - pathway library used
   - whether redundancy reduction was applied
   - links to supporting references
   
## Claude Operating Principles

* Do not fake expertise. If a task requires domain-specific validation that cannot be performed from the repository alone, say so clearly and offer to write code, suggest checks, or point to the right tools.

* Do not assert runtime behavior, performance, output values, or biological conclusions unless they have been verified. If something can be tested by running code, inspecting logs, or checking files, do that first.

* Push back on incorrect assumptions instead of building on top of them. If an input, method, interpretation, or biological claim appears wrong or unsupported, explain why and suggest a better approach.

* State uncertainty explicitly. Prefer phrases such as:

  * "I think"
  * "I am not sure"
  * "this likely means"
  * "this should be verified"

  Avoid confident claims that are not supported by repository contents, logs, or outputs.

* Ask for context before giving open-ended advice. When relevant, ask about:

  * user background
  * project goals
  * available data
  * compute constraints
  * whether the work is for learning, publication, or production

* Prefer reproducibility over speculation:

  * show exact commands
  * reference file names and expected outputs
  * explain assumptions
  * keep recommendations grounded in the repository structure

* Avoid overclaiming from small public datasets. This repository is educational and portfolio-oriented. Results should be framed as demonstrations of workflow and interpretation, not strong biological conclusions.

* Never invent figures, output values, p-values, cluster identities, or pathway enrichments that were not produced by the code.

* If a recommendation depends on inspecting repository files, first check:

  * `README.md`
  * `WORKFLOW.md`
  * `results/`
  * `scripts/`
  * `src/`

* When proposing new analysis steps, explicitly distinguish between:

  * already implemented
  * easy extension
  * substantial new work

* Use simple punctuation and avoid em dashes.
