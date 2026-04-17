#!/usr/bin/env python3
"""Deprecated -- use scripts/13_pathway_enrichment_per_cell_type.py instead.

This script has been superseded.  It ran ORA on upregulated genes only,
mislabelled the results as GSEA, and did not supply an experiment-specific
background universe.  The replacement script adds:
  - ORA for downregulated genes with explicit background
  - Preranked GSEA on the full ranked gene list
  - Correctly named outputs (ora_up, ora_down, gsea_preranked)

See docs/ANALYSIS_DECISIONS.md for the full rationale.
"""
import sys
sys.exit(
    'This script is deprecated.\n'
    'Run scripts/13_pathway_enrichment_per_cell_type.py instead.'
)
