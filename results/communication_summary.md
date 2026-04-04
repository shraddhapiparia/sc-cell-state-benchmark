# Cell-Cell Communication Results Summary

**Script:** `scripts/11_cell_communication.py`
**Dataset:** Kang PBMC IFN-beta stimulation (24,562 cells, 8 annotated cell types)
**Method:** Mean-product LR scoring —
`lr_score = mean_log_expr(ligand, sender) × mean_log_expr(receptor, receiver)`
**Panel:** 16 curated ligand-receptor pairs covering interferon, chemokine,
inflammatory, and regulatory axes.

> **Important caveat:** This is an exploratory analysis. Scores reflect expression
> levels in the log-normalised matrix; they are not measurements of actual
> ligand-receptor binding, paracrine signalling, or functional outcomes.
> High scores indicate co-expression of a ligand in sender cells and its receptor
> in receiver cells under a given condition — consistent with a potential interaction,
> but not proof of one.

---

## Gene availability

All 24 unique genes across the 16 pairs were found in the dataset. No pairs were
dropped. Several ligands are near-absent by expression (see internal controls below).

---

## Top interactions in the stimulated condition

The five highest-scoring individual interactions in the stimulated condition are all
driven by the **CXCL10 → CXCR4** pair:

| Rank | Interaction | Sender | Receiver | Score |
|---|---|---|---|---|
| 1 | CXCL10 → CXCR4 | FCGR3A+ Monocytes | B cells | 13.19 |
| 2 | CXCL10 → CXCR4 | CD14+ Monocytes | B cells | 13.08 |
| 3 | CXCL10 → CXCR4 | Dendritic cells | B cells | 10.66 |
| 4 | CXCL10 → CXCR4 | FCGR3A+ Monocytes | CD8 T cells | 10.59 |
| 5 | CXCL10 → CXCR4 | CD14+ Monocytes | CD8 T cells | 10.51 |

CXCL10 is one of the most strongly induced ISGs in this dataset (mean expression in
stimulated CD14+ monocytes: ~4.98, FCGR3A+ monocytes: ~5.02). CXCR4 is broadly
expressed across most cell types (B cells, CD4 T cells, CD8 T cells, NK cells all
above 60% expressing in ctrl), making monocytes and DCs — the dominant CXCL10
producers — high-scoring senders to essentially all lymphocyte populations via CXCR4.

The canonical CXCL10 → CXCR3 pair (CXCR3 is considered the primary CXCL10 receptor
for T/NK chemotaxis) also scores strongly in stim, but ranks behind CXCL10→CXCR4
because CXCR4 is more broadly expressed than CXCR3 in this dataset.

---

## Top sender → receiver pairs by stim − ctrl delta

Aggregating across all 16 LR pairs, the largest increases in communication potential
under stimulation are:

| Sender | Receiver | Delta |
|---|---|---|
| CD14+ Monocytes | CD8 T cells | +13.16 |
| CD14+ Monocytes | B cells | +12.98 |
| FCGR3A+ Monocytes | B cells | +10.997 |
| Dendritic cells | CD8 T cells | +10.77 |
| FCGR3A+ Monocytes | CD8 T cells | +10.74 |

The monocyte subtypes (CD14+ and FCGR3A+) and dendritic cells emerge as the
dominant senders of increased communication signal under IFN-beta stimulation. This
is consistent with their role as innate immune responders: they are the primary
producers of ISG-driven chemokines (CXCL10, CCL2, CCL3) in this dataset. Lymphocytes
(B cells, CD8 T cells, NK cells) appear predominantly as receivers, reflecting their
expression of the cognate receptors (CXCR4, CXCR3, CCR5).

---

## Program-level interpretation

### Strongly increased under stimulation (expected)

- **CXCL10 → CXCR4 and CXCL10 → CXCR3:** Dominant signal. CXCL10 is one of the
  most strongly ISG-induced genes; its potential interactions with CXCR4 (broad) and
  CXCR3 (NK/CD8-enriched) increase sharply in stim. Biologically consistent with
  IFN-beta-driven lymphocyte chemotaxis.
- **CCL2 → CCR2:** CCL2 increases strongly in stimulated monocytes (mean 1.76 ctrl
  → 4.57 stim in CD14+). CCR2 is very lowly expressed across all types in this dataset
  (consistently < 5%), so the score is dominated by the ligand. The CCL2 signal is
  real; its interpretation as a CCR2-mediated interaction is limited by the absence of
  detectable CCR2.
- **CCL3 → CCR5 and CCL4 → CCR5:** Moderate increases. CCR5 is upregulated in
  stimulated monocytes (CD14+: 13% ctrl → 33% stim), contributing to both the sender
  side (monocytes expressing CCL3/CCL4) and receiver side.

### Decreased or unchanged under stimulation

- **IL1B → IL1R1/IL1R2:** IL1B expression is slightly lower in stimulated monocytes
  (CD14+ ctrl ~1.03 mean, stim lower). This produces a slight negative delta for the
  IL1B pairs — consistent with IFN-beta suppressing some constitutive monocyte
  inflammatory tone.
- **CCL5 → CCR5 and CCL5 → CCR1:** CCL5 is constitutively high in NK and CD8 T cells
  and does not change substantially under IFN-beta stimulation. NK cells' CCR5
  expression is slightly lower in stim. This pair tracks constitutive NK/CD8 effector
  identity, not perturbation response.

### Internal negative controls (near-zero as expected)

- **IFNG → IFNGR1/IFNGR2:** IFNG is expressed in < 1% of cells in both conditions.
  Scores are near-zero, confirming no detectable T-cell IFN-gamma secretion in this
  short in-vitro IFN-beta stimulation protocol. These pairs behave as expected negative
  controls.
- **IL10 → IL10RA:** IL10 is similarly absent from the expression matrix. Score ≈ 0.
- **TGFB1 → TGFBR1:** TGFB1 is sparsely expressed (< 20% of cells in any group);
  scores are low but non-zero in monocyte populations.

---

## Assumptions and caveats

1. **Mean-product is a proxy.** Scores are the product of mean log-normalised expression
   values. They do not account for receptor saturation, co-receptor requirements,
   subcellular localisation, or secretion vs surface presentation of the ligand.

2. **CXCL10 → CXCR4 dominates by score magnitude.** CXCL10 is among the most
   strongly induced genes in the dataset; because CXCR4 is broadly expressed, the
   CXCL10→CXCR4 pair produces very high absolute scores in stim. This makes the
   heatmaps and rankings heavily CXCL10-driven. This is a real biological signal
   (ISG-driven chemokine) but means the analysis is less informative for subtler pairs.
   Normalising scores per pair before aggregating would change the ranking and is a
   straightforward future extension.

3. **No statistical testing.** Scores are not compared to a permuted-label null.
   High scores do not imply statistical significance.

4. **Megakaryocytes (n=10 ctrl, n=12 stim):** Cell count is very low. Scores involving
   Megakaryocytes as sender or receiver should not be interpreted confidently.

5. **CCR2 is present in the gene vocabulary but near-absent by expression.** The
   CCL2 → CCR2 score is dominated by CCL2 ligand expression. The receptor side provides
   a near-zero multiplier; the pair is included to show the strong CCL2 induction but
   does not represent a credible CCR2-mediated interaction in this dataset.

6. **In-vitro, short-term protocol.** Results may differ at longer time points,
   different IFN subtypes, or in in-vivo contexts.
