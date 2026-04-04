"""Curated ligand-receptor pairs and mean-product scoring for cell-cell communication.

Method
------
Each potential sender–receiver interaction is scored as:

    lr_score = mean_log_expr(ligand, sender_cells) × mean_log_expr(receptor, receiver_cells)

where the mean is taken across all cells of a given type and condition, using
the log-normalised expression matrix stored in adata.raw.

This mean-product approach is the scoring strategy used in CellPhoneDB v1 and
early NicheNet analyses. It is a proxy for interaction potential, not a direct
measurement of binding, signalling activity, or biological outcome.

Limitations
-----------
- No statistical null testing is performed. Scores reflect expression levels only.
- Several ligands (IFNG, IL10) are near-absent in this dataset; their pairs score
  near zero and serve as transparent internal negative controls.
- Scores are not corrected for cell-type abundance differences across conditions.
- mean-product can be dominated by the higher-expressing partner (usually the ligand
  for ISG-driven chemokines such as CXCL10 in the stimulated condition).
- Correlation between expression and communication potential ≠ causation.

Extending this module
---------------------
Add entries to LR_PAIRS as (ligand, receptor, interaction_name) tuples.
Script 11 will pick them up automatically.
"""

from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp


# ---------------------------------------------------------------------------
# Curated ligand-receptor panel
# ---------------------------------------------------------------------------
# 16 pairs covering the main interferon/inflammatory PBMC axes.
# Sources: literature-supported interactions for human PBMCs, consistent with
# CellPhoneDB v2 and published interferon/chemokine biology.
#
# Near-zero ligand pairs (IFNG, IL10) are intentionally retained so that the
# output honestly shows which pairs lack detectable ligand expression in this
# dataset. They function as internal negative controls.
#
# Format: (ligand, receptor, display_name)
LR_PAIRS: List[Tuple[str, str, str]] = [
    # IFN-driven chemokine axes — expected to increase strongly under IFN-beta stim
    ('CXCL10', 'CXCR3',    'CXCL10 \u2192 CXCR3'),    # canonical IP-10→T/NK chemotaxis
    ('CXCL10', 'CXCR4',    'CXCL10 \u2192 CXCR4'),    # secondary CXCL10 receptor
    ('CXCL9',  'CXCR3',    'CXCL9 \u2192 CXCR3'),     # MIG→T/NK, ISG-driven

    # CC-chemokine axes — monocyte/NK/T communication
    ('CCL2',   'CCR2',     'CCL2 \u2192 CCR2'),        # MCP-1 monocyte recruitment; CCR2 lowly expressed
    ('CCL3',   'CCR5',     'CCL3 \u2192 CCR5'),        # MIP-1α→monocyte/NK
    ('CCL4',   'CCR5',     'CCL4 \u2192 CCR5'),        # MIP-1β→monocyte/NK
    ('CCL5',   'CCR5',     'CCL5 \u2192 CCR5'),        # RANTES→monocyte/NK; constitutive NK/CD8 source
    ('CCL5',   'CCR1',     'CCL5 \u2192 CCR1'),        # secondary CCL5 receptor

    # TNF axis — inflammatory; TNF expression is near-zero in most cell types
    ('TNF',    'TNFRSF1A', 'TNF \u2192 TNFRSF1A'),    # TNF-R1, ubiquitous; TNF ligand sparse
    ('TNF',    'TNFRSF1B', 'TNF \u2192 TNFRSF1B'),    # TNF-R2

    # IL-1 axis — monocyte-driven; mainly ctrl signal (IL1B decreases in stim monocytes)
    ('IL1B',   'IL1R1',    'IL1B \u2192 IL1R1'),
    ('IL1B',   'IL1R2',    'IL1B \u2192 IL1R2'),       # decoy receptor

    # IFN-gamma axis — IFNG is near-absent in this IFN-beta stimulation dataset;
    # retained as negative control to confirm lack of T-cell IFN-gamma secretion
    ('IFNG',   'IFNGR1',   'IFNG \u2192 IFNGR1'),     # near-zero ligand — internal control
    ('IFNG',   'IFNGR2',   'IFNG \u2192 IFNGR2'),     # near-zero ligand — internal control

    # Anti-inflammatory / regulatory axes — both ligands very lowly expressed
    ('IL10',   'IL10RA',   'IL10 \u2192 IL10RA'),     # near-zero ligand — internal control
    ('TGFB1',  'TGFBR1',   'TGFB1 \u2192 TGFBR1'),   # low ligand expression in this dataset
]


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def _mean_expr_by_group(
    adata,
    gene: str,
    cell_types: List[str],
    conditions: List[str],
    ct_col: str,
    cond_col: str,
) -> Dict[Tuple[str, str], float]:
    """Return mean log-normalised expression of `gene` for every (cell_type, condition) pair.

    Uses adata.raw when available (full gene matrix). Returns 0.0 for absent genes
    or empty groups rather than raising an error, so near-zero pairs score honestly.
    """
    raw = adata.raw if adata.raw is not None else adata
    var_names = list(raw.var_names)

    result: Dict[Tuple[str, str], float] = {
        (ct, c): 0.0 for ct in cell_types for c in conditions
    }
    if gene not in var_names:
        return result

    gi = var_names.index(gene)
    X = raw.X
    col = X[:, gi]
    if sp.issparse(col):
        col = np.asarray(col.todense()).ravel()
    else:
        col = np.asarray(col).ravel()

    ct_arr = adata.obs[ct_col].astype(str).values
    cond_arr = adata.obs[cond_col].astype(str).values

    for ct in cell_types:
        for c in conditions:
            mask = (ct_arr == ct) & (cond_arr == c)
            if mask.sum() > 0:
                result[(ct, c)] = float(col[mask].mean())
    return result


def score_lr_pairs(
    adata,
    ct_col: str,
    cond_col: str,
) -> pd.DataFrame:
    """Compute mean-product LR scores for all pairs, senders, receivers, and conditions.

    Parameters
    ----------
    adata : AnnData
        Preprocessed dataset. Uses adata.raw for expression values.
    ct_col : str
        Column in adata.obs containing cell type labels.
    cond_col : str
        Column in adata.obs containing condition labels (e.g. 'ctrl'/'stim').

    Returns
    -------
    DataFrame with columns:
        ligand, receptor, interaction_name, sender, receiver, condition,
        mean_ligand_expr, mean_receptor_expr, lr_score
    """
    cell_types = sorted(adata.obs[ct_col].astype(str).unique())
    conditions = sorted(adata.obs[cond_col].astype(str).unique())

    rows = []
    for ligand, receptor, interaction_name in LR_PAIRS:
        lig_expr = _mean_expr_by_group(adata, ligand, cell_types, conditions, ct_col, cond_col)
        rec_expr = _mean_expr_by_group(adata, receptor, cell_types, conditions, ct_col, cond_col)

        for cond in conditions:
            for sender in cell_types:
                for receiver in cell_types:
                    l_val = lig_expr[(sender, cond)]
                    r_val = rec_expr[(receiver, cond)]
                    rows.append({
                        'ligand': ligand,
                        'receptor': receptor,
                        'interaction_name': interaction_name,
                        'sender': sender,
                        'receiver': receiver,
                        'condition': cond,
                        'mean_ligand_expr': round(l_val, 5),
                        'mean_receptor_expr': round(r_val, 5),
                        'lr_score': round(l_val * r_val, 6),
                    })

    df = pd.DataFrame(rows)
    # Rank within each condition (1 = highest score)
    df['rank_within_condition'] = (
        df.groupby('condition')['lr_score']
        .rank(ascending=False, method='min')
        .astype(int)
    )
    return df


def build_sender_receiver_matrix(
    scores_df: pd.DataFrame,
    condition: str,
    cell_types: List[str],
) -> pd.DataFrame:
    """Sum LR scores across all pairs for a given condition into a sender × receiver matrix."""
    sub = scores_df[scores_df['condition'] == condition]
    matrix = (
        sub.groupby(['sender', 'receiver'])['lr_score']
        .sum()
        .unstack(fill_value=0.0)
    )
    # Reindex to ensure consistent ordering even if some groups are missing
    matrix = matrix.reindex(index=cell_types, columns=cell_types, fill_value=0.0)
    return matrix
