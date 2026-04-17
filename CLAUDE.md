# CLAUDE.md

This file provides guidance to Claude Code when working in this repository.

## Project Overview

This repository demonstrates progressively more advanced single-cell analysis workflows:

1. RNA-only PBMC cell-state benchmarking
2. Paired RNA+ATAC multiome integration with WNN
3. TF / regulon activity inference
4. Cross-modal RNA vs ATAC validation
5. Optional pseudotime and dynamic regulatory analysis

The project is organized as a sequence of numbered scripts. Each script should:

* read outputs from earlier phases
* write clearly named files
* avoid modifying previous outputs in-place
* save both tables and figures when possible

## Main Analysis Phases

### Phase 1: RNA-only Benchmarking

Goal:
Establish baseline single-cell analysis and cell-state scoring using PBMC RNA datasets.

Typical workflow:

```text
01_import_rna_baseline.py
→ QC, normalization, HVGs, PCA, UMAP, Leiden
→ cluster annotation
→ cell-state scoring
```

Important outputs:

* `results/figures/annotated_umap.png`
* `results/cell_states/*.csv`
* `results/project_summary.md`

### Phase 2: RNA + ATAC Integration

Goal:
Integrate paired multiome data using Seurat WNN.

Workflow:

```text
02_download_or_prepare_multiome.R
03_preprocess_multiome.R
04_integrate_rna_atac_wnn.R
05_compare_rna_vs_wnn.py
```

Expected outputs:

* `data/processed/multiome_wnn_object.rds`
* `results/comparison/rna_vs_wnn_broad_summary.md`
* WNN UMAP figures and broad-label comparisons

Important conventions:

* RNA assay should remain default unless temporarily switching to ATAC
* Broad cell-type labels should be saved in metadata
* Keep mappings between original WNN labels and broad labels in TSV files

## Phase 3: TF / Regulon Activity

Scripts:

```text
06_tf_activity.py
07_plot_tf_activity.py
```

Goal:
Estimate TF activity from RNA programs for a small set of biologically interpretable axes:

* interferon_alpha_beta_response
* inflammatory_response
* antigen_presentation_mhc
* cytotoxic_nk_like
* monocyte_activation

Outputs:

* `results/regulons/tf_activity_rna_by_celltype.csv`
* `results/figures/tf_activity/*.png`

Guidelines:

* Prefer small interpretable TF sets rather than large black-box regulons
* Save both per-cell and per-cell-type summaries if possible
* Always print top TFs and strongest cell-type enrichments

## Phase 4: RNA vs ATAC Cross-Validation

Scripts:

```text
08_motif_enrichment_atac.R
09_compare_rna_atac_tf.py
```

Goal:
Validate RNA-derived TF activity using ATAC motif accessibility.

Expected outputs:

* `results/regulons/atac_motif_deviations.csv`
* `results/comparison/tf_rna_atac_correlation.csv`
* scatter plots and correlation summary tables

Interpretation:

* Strong positive correlation supports the TF signal biologically
* Correlations around 0.8–0.9 are strong
* Interferon signatures may be weaker due to limited matching motifs

## Optional Phase 5: Pseudotime

Goal:
Study dynamic transitions, especially in monocyte lineage cells.

Expected workflow:

```text
subset monocyte-related cells
→ run Slingshot or Monocle
→ project RNA programs and ATAC motif activity along pseudotime
```

Expected outputs:

* `results/dynamics/pseudotime_monocyte_cells.csv`
* `results/figures/dynamics/pseudotime_umap_monocyte.png`
* RNA and ATAC trend plots along pseudotime

Important rule:
Slingshot requires at least two cluster labels within the subset. If only one label is present, stop and print a useful error.

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
