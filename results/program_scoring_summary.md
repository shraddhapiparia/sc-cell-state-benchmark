# Program Scoring Results Summary

**Script:** `scripts/10_score_pathway_programs.py`
**Dataset:** Kang PBMC IFN-beta stimulation (24,562 cells, 8 annotated cell types)
**Condition column:** `label` (values: `ctrl`, `stim`)
**Scoring method:** Scanpy `score_genes` (background-corrected)

---

## AUROC ranking: program separation of stimulated vs control

| Program | AUROC | Genes found | Interpretation |
|---|---|---|---|
| IFN-α/β response | **0.991** | 10/10 | Strong perturbation signal — expected positive control |
| Inflammatory response | **0.655** | 7/8 (CXCL8 absent) | Moderate perturbation signal — IFN-beta induces a subset of inflammatory genes |
| Antigen presentation (MHC-II) | 0.487 | 7/7 | Near-random — cell-type identity program, not perturbation-responsive |
| Cytotoxic / NK-like | 0.487 | 7/7 | Near-random globally — cell-type identity program |
| Monocyte activation | 0.452 | 7/7 | Below 0.5 — monocyte identity markers slightly decrease under IFN-beta stimulation |

---

## Key observations by program

### IFN-α/β response (AUC 0.991, effect size 1.72)

The strongest and cleanest perturbation signal in the dataset, consistent with this
being a matched IFN-beta stimulation experiment. The response is present across all
cell types, with the largest absolute score difference in CD14+ monocytes and FCGR3A+
monocytes (~4 score units delta), followed by dendritic cells. Lymphocytes show a
smaller but clear response. This is consistent with monocytes being potent innate
immune responders to type I interferon.

Cell type means (ctrl → stim):
- CD14+ Monocytes: -1.07 → +2.90 (delta +3.97)
- FCGR3A+ Monocytes: -0.74 → +2.61 (delta +3.35)
- Dendritic cells: -1.29 → +2.11 (delta +3.40)
- NK cells: -0.72 → +1.68 (delta +2.40)
- B cells: -1.01 → +1.43 (delta +2.44)
- CD4 T cells: -1.11 → +1.07 (delta +2.18)
- CD8 T cells: -0.88 → +1.25 (delta +2.13)
- Megakaryocytes: -0.78 → +0.98 (delta +1.76, n=10/12 — interpret with caution)

### Inflammatory response (AUC 0.655, effect size 0.68)

Moderate perturbation signal, driven mainly by monocyte subtypes and dendritic cells.
The global signal is weaker than the IFN program because IFN-beta is not the primary
driver of classical NF-κB inflammatory genes (IL6, IL1B, TNF). The signal captured here
likely reflects: (1) CXCL10, which is also an ISG and dominates the inflammatory score
in stimulated cells; (2) CCL2 and CCL4 chemokines, which are modestly induced by IFN-beta
in monocytes. Note CXCL8 (IL-8) was absent from this dataset's gene vocabulary, reducing
the set to 7/8 genes.

CD14+ Monocytes show the largest inflammatory delta (+1.29). NK cells, CD8 T cells, and
CD4 T cells show minimal change, consistent with these cell types being less responsive
to IFN-beta-driven inflammatory signalling at 6 hours.

### Antigen presentation / MHC-II (AUC 0.487 — cell-type context program)

As expected, this program primarily tracks cell-type identity rather than perturbation
state. Dendritic cells have the highest scores (~3.0), followed by B cells (~2.2) and
monocyte subtypes (~1.2). T cells and NK cells are near or below zero. The slightly lower
scores in stimulated cells for some types (B cells: 2.19 → 1.83; DCs: 3.03 → 2.67) are
small in magnitude and may reflect the background correction in Scanpy's score_genes
shifting when the transcriptional landscape changes under stimulation. No strong
directional perturbation effect is expected or observed. This program is included to
provide cell-type context in the heatmap, not to measure a perturbation response.

### Cytotoxic / NK-like (AUC 0.487 — cell-type context program)

Clearly delineates NK cells and CD8 T cells from all other types (NK ctrl: 1.85, NK
stim: 2.28; CD8 ctrl: 1.32, CD8 stim: 1.57). All other cell types cluster near -0.4
to -0.5 regardless of condition. The slight upward shift in NK and CD8 scores under
stimulation (approximately +0.4) is consistent with IFN-beta having documented
activity-enhancing effects on NK cells, but this is a small signal relative to between
cell-type differences and should not be over-interpreted from this data alone.

### Monocyte activation (AUC 0.452 — cell-type context program, inverse perturbation effect)

This program marks monocytes (CD14+ ctrl: 1.54; FCGR3A+ ctrl: 0.61) and DCs (ctrl: 0.65)
clearly above lymphocytes, confirming its utility as a cell-type identity signal. The
AUROC below 0.5 reflects a small but systematic *decrease* in monocyte activation scores
under IFN-beta stimulation (CD14+ monocytes: 1.54 → 1.08; FCGR3A+ monocytes: 0.61 → 0.36;
DCs: 0.65 → 0.38). A plausible interpretation is that IFN-beta stimulation shifts monocyte
gene expression toward ISG programs and away from constitutive identity/alarmin expression
(S100A8, S100A9, LYZ), producing a relative downward shift in monocyte activation scores.
This interpretation is speculative and not confirmed by this dataset alone.

---

## Internal consistency checks

- IFN-α/β response scores are consistent with those in `kang_interferon_scores.csv`
  produced by script 09 (same gene set, same scoring method, same dataset).
- The three cell-type identity programs produce score patterns that match the known
  biology of the annotated cell types, providing an independent cross-check that the
  cell type annotations in `adata.obs.cell_type` are plausible.
- CXCL8 is absent from this dataset's gene vocabulary. All other 38/39 genes were found.
- Megakaryocytes have very low cell counts (n=10 ctrl, n=12 stim). Summary statistics
  for this cell type are included in the tables but should not be interpreted with
  confidence.

---

## Assumptions and caveats

- Scores are produced by Scanpy `score_genes`, which is background-corrected against
  randomly sampled reference genes. This reduces library-size confounding but means
  scores are not on an absolute expression scale.
- Gene sets are minimal curated proxies (7–10 genes), not comprehensive pathway
  signatures. Scores are approximate qualitative indicators.
- The three cell-type identity programs (antigen presentation, cytotoxic, monocyte
  activation) have near-random AUROC for ctrl-vs-stim comparison. This is the expected
  result and does not indicate a problem with the scoring method or the data.
- Cross-program gene membership was minimised to reduce score correlation. CXCL10 is
  biologically an ISG but was placed in the inflammatory_response set; this contributes
  to that program's AUC of 0.655 rather than 0.5.
- The dataset uses a 6-hour IFN-beta stimulation protocol. Longer time points or
  different IFN subtypes may produce different program patterns.
