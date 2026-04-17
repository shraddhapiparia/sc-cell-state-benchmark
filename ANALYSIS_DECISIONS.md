## Pathway analysis decision log

### Script
`scripts/13_pathway_enrichment_per_cell_type.py`

### References

- [Subramanian et al. 2005 PNAS -- original GSEA paper](https://pubmed.ncbi.nlm.nih.gov/16199517/)
- [GSEAPreranked v6 module documentation (Broad Institute)](https://gsea-msigdb.github.io/gseapreranked-gpmodule/v6/index.html)
- [Pathway Commons GSEA primer](https://www.pathwaycommons.org/guide/primers/data_analysis/gsea/)
- [Irizarry et al. 2009 -- Gene Set Enrichment Analysis Made Simple](https://pmc.ncbi.nlm.nih.gov/articles/PMC3134237/)
- [Wijesooriya et al. 2024 -- Two subtle problems with overrepresentation analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC11557902/)
- [Zyla et al. 2014 -- Separate enrichment analysis for up- and downregulated genes](https://pmc.ncbi.nlm.nih.gov/articles/PMC3899863/)
- [Peset et al. 2022 -- Nine quick tips for pathway enrichment analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC9371296/)

---

### Best-practice rules

**ORA vs GSEA are different methods and must not be conflated.**

| Dimension | ORA | Preranked GSEA |
|---|---|---|
| Input | Binary foreground list vs background | Continuous ranked list of all tested genes |
| Null model | Hypergeometric/Fisher (random sampling) | Gene-set permutation of the ranked list |
| Question | Is this pathway over-represented among my DEGs? | Is this pathway coordinately shifted across the full ranking? |
| Direction | Requires separate up / down runs | Natively captures both ends in one run |
| Key risk | Wrong background inflates false positives | Gene-set permutation is anti-conservative vs phenotype permutation |

**ORA rules**

1. The foreground must be split: run ORA separately on upregulated genes and on
   downregulated genes. Combining them dilutes directional signal (Zyla et al. 2014
   showed combined analysis found only 1 pathway where directional analysis found 20).
2. The background universe must be the set of genes actually tested in the DE
   analysis, not the whole genome. Using the whole genome artificially inflates the
   denominator and distorts p-values (Wijesooriya et al. 2024).
3. Report Benjamini-Hochberg adjusted p-values. Save results with adj p < 0.25 as
   the broad archive; use adj p < 0.05 for high-confidence interpretation.
4. Flag any result where the observed overlap is 3 or fewer genes as low-confidence
   regardless of p-value.

**Preranked GSEA rules**

1. Use the full, unfiltered gene list as input. Filtering to DEGs before ranking
   defeats the purpose of the method.
2. Recommended ranking metric for RNA-seq: `sign(logFC) * -log10(padj + epsilon)`.
   This captures both direction and statistical confidence.
3. Remove duplicate gene names before ranking; keep the entry with the highest
   absolute rank value.
4. Use a minimum gene set size of 15 and maximum of 500 (Broad Institute defaults).
5. Use gene-set permutation (the only option available for preranked mode). Report
   FDR q-value, not nominal p-value. Save results with FDR q < 0.25 (Broad
   recommendation for hypothesis generation).
6. Preranked GSEA on single-cell pseudo-bulk data lacks biological replicates.
   Results should be framed as methodology demonstration, not strong biological
   conclusions.

**Naming rules**

- Files, plot titles, and markdown must distinguish ORA from GSEA.
- ORA outputs use the prefix `ora_up` or `ora_down`.
- GSEA outputs use the prefix `gsea_preranked`.
- Never name an ORA output file or plot title with the word "GSEA".

---

### What this repo had before this change

- Only ORA on upregulated genes was implemented (no downregulated ORA, no GSEA).
- The Enrichr API was called without a custom background, so the background
  defaulted to Enrichr's internal gene universe rather than the experiment-specific
  tested-gene list.
- All outputs were mislabeled as GSEA (`kang_gsea_results.csv`,
  `kang_gsea_top_pathways.png`, `[gsea]` log prefix).
- The plotting helper was named `plot_gsea_dotplot` even though it visualised ORA results.

---

### What was changed and why

| Change | Reason |
|---|---|
| Added ORA for downregulated genes | Separate directional ORA has substantially higher power than combined or upregulated-only analysis (Zyla et al. 2014) |
| Added preranked GSEA on full ranked list | GSEA is complementary to ORA: it detects coordinated weak signals invisible to ORA and does not require a p-value cutoff |
| Switched to explicit background for ORA | Using Enrichr's internal universe instead of the tested-gene list is a documented source of false positives and false negatives (Wijesooriya et al. 2024) |
| Renamed outputs from `gsea_*` to `ora_up_*`, `ora_down_*`, `gsea_preranked_*` | Prevents confusion between two fundamentally different statistical methods |
| Renamed `plot_gsea_dotplot` to `plot_enrichment_dotplot`; added `plot_gsea_prerank_dotplot` | Correct method labels in function names |
| Updated log prefixes to `[ora]` vs `[gsea]` per analysis type | Makes console output readable and methodologically correct |

---

### Methods (current implementation)

- ORA on upregulated genes (padj < 0.05, logFC > 0.5) with explicit tested-gene background
- ORA on downregulated genes (padj < 0.05, logFC < -0.5) with explicit tested-gene background
- Preranked GSEA on full tested gene list (ranking metric: sign(logFC) * -log10(padj + 1e-300))

### Ranking metric

`sign(logfoldchanges) * -log10(pvals_adj + 1e-300)`

Positive values = upregulated with small p-value (top of ranked list).
Negative values = downregulated with small p-value (bottom of ranked list).

### Background universe

All genes tested in DE for that cell type. Passed explicitly to `gp.enrichr(background=...)`.

### Gene set libraries

- MSigDB Hallmark 2020 for compact, low-redundancy interpretation.
- Optional: GO BP / Reactome for broader biological follow-up (not yet implemented).

### Redundancy handling

Hallmark reduces but does not eliminate redundancy. No additional filtering is applied.
