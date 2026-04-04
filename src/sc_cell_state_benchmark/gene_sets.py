"""Curated immune gene sets for pathway/program scoring in human PBMCs.

These are manually assembled minimal signatures intended to capture broad
cell-state variation in human PBMC scRNA-seq data. They are not derived from
any licensed database and no external database is required to run this code.
Gene symbols are human HGNC, consistent with 10x PBMC datasets.

Assumptions and limitations
----------------------------
On cross-program gene membership:

  Some immune genes legitimately participate in more than one biological program.
  For example, CXCL10 is both a canonical interferon-stimulated gene and a
  pro-inflammatory chemokine; HLA-DRA is both a component of the antigen
  presentation machinery and a constitutive marker of monocytes. Placing a gene
  in multiple programs is biologically valid when its expression genuinely spans
  those contexts.

  In this v1 implementation gene overlap across programs has been minimised as a
  pragmatic choice for interpretability: fewer shared genes reduces correlation
  between program scores, making the heatmap and AUROC comparisons easier to read.
  This is a design decision, not a biological claim that these programs are
  mutually exclusive. Future versions may relax this constraint where biologically
  justified (see docs/NEXT_STEPS.md).

  Specific cases where the biological overlap is real but was resolved by placing
  the gene in only one program:

  CXCL10 -- assigned to inflammatory_response. It is equally a canonical ISG
             induced by type I IFN. The interferon_alpha_beta_response set does
             not include it in order to keep the two program scores more
             independent, not because it does not belong there biologically.
  IRF7   -- assigned to interferon_alpha_beta_response. It is a master regulator
             of type I IFN and relevant to broader innate immune activation.
  HLA-DRA, CD74 -- assigned to antigen_presentation_mhc. Both are also
             constitutive monocyte/DC markers, so antigen_presentation_mhc and
             monocyte_activation scores will correlate in those cell types even
             though no gene is shared.

Program-level caveats:

  antigen_presentation_mhc, cytotoxic_nk_like, monocyte_activation -- these
      three programs are expected to behave primarily as cell-type identity or
      context markers in the Kang IFN-beta dataset, not as perturbation-responsive
      programs. That is a feature of the design, not a flaw: they provide
      interpretable biological context that complements the two programs that
      do show perturbation signal (interferon_alpha_beta_response and, to a
      lesser extent, inflammatory_response). Specifically:

      antigen_presentation_mhc: MHC class II genes are constitutively expressed
          on monocytes, DCs, and B cells. IFN-gamma (not IFN-beta) is the
          canonical inducer of MHC-II upregulation, so the perturbation effect
          here is expected to be modest.
      cytotoxic_nk_like: cytotoxic granule genes mark NK and CD8 T cells. They
          are not substantially perturbed by short-term IFN-beta stimulation.
          Expected pattern: high in NK/CD8, low elsewhere, regardless of condition.
      monocyte_activation: constitutive monocyte identity genes plus alarmins
          (S100A8/S100A9). The alarmins can increase under inflammatory activation
          but the overall set reads primarily as a cell-type signature.

  inflammatory_response -- some genes (IL6, IL1B, TNF, CXCL8) are weakly
      expressed at rest in PBMCs and may not be detectable in all cells. Gene
      availability is reported at runtime.

Set sizes are kept small (7–10 genes) to stay within the scoring regime
validated in scripts 06 and 09. Larger sets reduce per-cell score variance
but may include genes with lower specificity.

Extending this module
----------------------
To add a program, add an entry to IMMUNE_PROGRAMS and PROGRAM_DISPLAY_NAMES.
The scoring pipeline in script 10 will pick it up automatically.
"""

from typing import Dict, List

# Minimum genes from a program that must be present in the dataset for scoring
# to proceed. Below this threshold the score is unreliable and the program is
# skipped with a warning.
MIN_GENES_REQUIRED: int = 3

# Human-readable labels used in all figures and tables.
PROGRAM_DISPLAY_NAMES: Dict[str, str] = {
    'interferon_alpha_beta_response': 'IFN-\u03b1/\u03b2 response',
    'inflammatory_response':          'Inflammatory response',
    'antigen_presentation_mhc':       'Antigen presentation (MHC-II)',
    'cytotoxic_nk_like':              'Cytotoxic / NK-like',
    'monocyte_activation':            'Monocyte activation',
}

# Curated gene sets. Key = program identifier used in column names and tables.
# Value = list of HGNC gene symbols.
IMMUNE_PROGRAMS: Dict[str, List[str]] = {

    # Interferon alpha/beta response
    # Classic interferon-stimulated genes (ISGs) induced via JAK-STAT1/2
    # downstream of type I interferon receptors. This is the primary perturbation
    # signal in the Kang IFN-beta dataset and is used as a positive control.
    # Note: CXCL10 is biologically an ISG but is assigned to inflammatory_response
    # to reduce score correlation between the two programs (see module docstring).
    # Informed by: Schoggins & Rice (2011) Nat Rev Immunol; MSigDB Hallmark
    # INTERFERON_ALPHA_RESPONSE (used as a reference, not directly copied).
    'interferon_alpha_beta_response': [
        'IFIT1',  # IFN-induced protein with tetratricopeptide repeats 1, antiviral
        'IFIT2',  # IFN-induced protein with tetratricopeptide repeats 2
        'IFIT3',  # IFN-induced protein with tetratricopeptide repeats 3
        'ISG15',  # ISG15 ubiquitin-like modifier, one of the most strongly induced ISGs
        'MX1',    # Dynamin-like GTPase, antiviral restriction factor
        'OAS1',   # 2'-5' oligoadenylate synthetase 1, activates RNase L
        'OASL',   # OAS-like, antiviral
        'RSAD2',  # Radical SAM domain 2 (viperin), broad-spectrum antiviral effector
        'IRF7',   # Interferon regulatory factor 7, master regulator of type I IFN response
        'STAT1',  # Signal transducer and activator of transcription 1, JAK-STAT pathway
    ],

    # Inflammatory response
    # Pro-inflammatory cytokines and chemokines regulated by NF-κB and related
    # pathways. Some genes (IL6, IL1B, TNF, CXCL8) are weakly expressed at
    # baseline in resting PBMCs; gene availability is reported at runtime.
    # Note: CXCL10 overlaps with interferon_alpha_beta_response.
    # Informed by: MSigDB Hallmark INFLAMMATORY_RESPONSE (reference only).
    'inflammatory_response': [
        'IL6',    # Interleukin-6, pleiotropic pro-inflammatory cytokine
        'TNF',    # Tumor necrosis factor, NF-κB target
        'IL1B',   # Interleukin-1 beta, inflammasome product
        'CXCL8',  # IL-8, neutrophil chemoattractant
        'CCL2',   # MCP-1, monocyte recruitment chemokine
        'CCL4',   # MIP-1β, inflammatory chemokine released by monocytes and NK cells
        'NFKBIA', # IκBα, NF-κB pathway target and negative regulator
        'CXCL10', # IP-10, IFN-inducible chemokine and ISG; assigned here rather than
                  # interferon_alpha_beta_response to reduce cross-program correlation
                  # (see module docstring -- biologically it belongs to both contexts)
    ],

    # Antigen presentation / MHC class II
    # MHC class II genes constitutively expressed on monocytes, dendritic cells,
    # and B cells. IFN-gamma is the canonical inducer of MHC-II upregulation;
    # changes under IFN-beta may be modest. Primarily a cell-type identity marker
    # in this dataset.
    # Note: HLA-DRA and CD74 also overlap with monocyte identity (see module docstring).
    # Informed by: Gene Ontology GO:0019882 (antigen processing and presentation).
    'antigen_presentation_mhc': [
        'HLA-DRA',  # MHC class II alpha chain, constitutive on antigen-presenting cells
        'HLA-DRB1', # MHC class II DR beta 1 chain
        'HLA-DQA1', # MHC class II DQ alpha 1 chain
        'HLA-DPB1', # MHC class II DP beta 1 chain
        'CD74',     # Invariant chain, required for MHC-II assembly and peptide loading
        'HLA-DMA',  # HLA-DM alpha, catalyses peptide exchange on MHC-II
        'HLA-DMB',  # HLA-DM beta
    ],

    # Cytotoxic / NK-like effector program
    # Cytotoxic granule proteins expressed in NK cells and CD8 T cells.
    # Expected pattern: high in NK and CD8 clusters, low elsewhere.
    # Not strongly perturbed by short-term IFN-beta stimulation; this program
    # functions primarily as a cell-type identity signature in this dataset.
    # Informed by: Crinier et al. (2018) Immunity; canonical NK/CTL literature.
    'cytotoxic_nk_like': [
        'GNLY',   # Granulysin, cytotoxic granule protein with antimicrobial activity
        'PRF1',   # Perforin-1, pore-forming cytotoxic effector molecule
        'GZMA',   # Granzyme A, serine protease in cytotoxic granules
        'GZMB',   # Granzyme B, primary CTL/NK effector protease (induces apoptosis)
        'GZMH',   # Granzyme H
        'NKG7',   # NK cell granule protein 7, cytotoxic marker
        'KLRD1',  # CD94, C-type lectin NK receptor
    ],

    # Monocyte activation
    # Constitutive monocyte identity genes and alarmins. Low in lymphocytes.
    # S100A8/S100A9 are pro-inflammatory alarmins that can increase under
    # activation, but the overall set is a cell-type identity marker.
    # Informed by: Villani et al. (2017) Science; Ziegler-Heitbrock et al. (2010)
    # monocyte nomenclature.
    'monocyte_activation': [
        'LYZ',    # Lysozyme, canonical monocyte/neutrophil marker
        'S100A8', # S100 calcium-binding protein A8, alarmin (calprotectin subunit)
        'S100A9', # S100 calcium-binding protein A9, alarmin (heterodimerises with S100A8)
        'FCN1',   # Ficolin-1, monocyte pattern recognition receptor
        'VCAN',   # Versican, extracellular matrix proteoglycan, classical monocyte marker
        'CST3',   # Cystatin C, monocyte/DC marker
        'CD14',   # LPS co-receptor, defines CD14+ classical monocytes
    ],
}
