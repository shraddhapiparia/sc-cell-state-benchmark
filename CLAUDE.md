# CLAUDE.md — sc-cell-state-benchmark

## Goal
Support reliable development of a reproducible single-cell RNA-seq benchmarking workflow for perturbation-aware cell-state scoring in the Kang IFN-β PBMC dataset.

Correctness matters more than speed.  
Verification matters more than confident explanation.

## Core Rules

- Do not fake expertise. If a task needs biological validation beyond repository contents, say so clearly and suggest the right checks.
- Do not hallucinate metadata, file names, assay names, cell-type labels, or output values.
- Do not claim runtime behavior, biological findings, or figure contents unless they were verified from code, logs, or files.
- If something can be checked by running code, inspecting AnnData objects, or reading outputs, do that before explaining.
- Push back on incorrect assumptions instead of building on top of them.
- State uncertainty explicitly when unsure.
- Verify current library syntax and API behavior against documentation instead of relying on memory.
- Update `.memory/` after every successful pipeline run or major debugging milestone.
- Avoid em dashes. Use simple sentence structure instead.

## Repo Scope

This repository is RNA-only benchmarking on the Kang PBMC IFN-β dataset.

Do not introduce ATAC-specific guidance or multiome assumptions unless explicitly requested.

## Active Workflow

Canonical reproduction:

    conda env create -f environment.yml
    conda activate sc-benchmark
    pip install -e .
    bash run_pipeline.sh

Expected workflow structure:

- 01–06: PBMC3k warm-up, QC, clustering, annotation, scoring
- 07–09: Kang IFN-β benchmarking
- 10: curated immune programs
- 11: cell communication (exploratory)
- 12: differential expression
- 13: pathway enrichment

Each script should:

- read outputs from prior steps  
- write outputs to `results/` or `figures/`  
- not overwrite previous outputs  
- print input path, cell counts, and output path  

## Project Structure

- scripts/ → workflow steps  
- src/sc_cell_state_benchmark/ → reusable functions  
- results/ → tables  
- figures/ → plots  
- docs/ → decisions and notes  
- data/ → local data (usually not committed)  

## Single-Cell Constraints

- Always verify cell counts after filtering  
- Always check `adata.obs.columns` before grouping  
- Distinguish:
  - cluster labels
  - cell types
  - conditions
  - scores  
- Do not over-interpret UMAP structure  
- Use negative controls where relevant  
- Keep ORA and GSEA separate  

## Validation Checklist

- Cell counts consistent?  
- Metadata columns correct?  
- Gene IDs mapped correctly?  
- Score distributions reasonable?  
- Known biology (IFN response) visible?  
- Outputs saved correctly?  

## Pathway Analysis Rules

- GSEA → use full ranked list  
- ORA → use defined gene sets + background  
- Document decisions in docs/  

## Memory

Maintain `.memory/`:

Failure:
- metadata column mismatch  

Fix:
- standardized column names  

Lesson:
- inspect metadata before use  

## Example Prompts

Design:  
"List failure modes before writing code."

Implementation:  
"Write minimal correct code given verified schema."

Debugging:  
"This result looks wrong. Give likely causes and checks."

## Final Principle

Reproducibility and biological plausibility matter more than clean-looking plots.
