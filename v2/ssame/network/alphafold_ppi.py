"""Loader for AlphaFold-predicted PPI (Baker Lab, RoseTTAFold2-PPI).

Data source: Zhang et al., Science 2024 (doi:10.1126/science.adt1630)
Download: https://datadryad.org/dataset/doi:10.5061/dryad.15dv41p84
File: final_predictions_90.tsv (17,849 interactions at 90% precision)
"""

from __future__ import annotations

from pathlib import Path

import networkx as nx
import pandas as pd

_CONFIDENCE_CANDIDATES = ("max_prob", "confidence", "score", "probability", "combined_score")
_SOURCE_CANDIDATES = ("source", "prediction_source", "type")


def load_alphafold_ppi(
    path: Path,
    *,
    min_confidence: float = 0.9,
    gene_a_col: str | None = None,
    gene_b_col: str | None = None,
    confidence_col: str | None = None,
) -> nx.Graph:
    """Load Baker Lab AlphaFold-predicted PPI network.

    Args:
        path: Path to the interaction file (TSV).
        min_confidence: Minimum confidence threshold (default: 0.9 = 90% precision).
        gene_a_col: Column name for gene A (auto-detected if None).
        gene_b_col: Column name for gene B (auto-detected if None).
        confidence_col: Column name for confidence score (auto-detected if None).

    Returns:
        Undirected NetworkX graph with Hugo gene symbol nodes.
    """
    if not path.exists():
        import warnings

        warnings.warn(
            f"AlphaFold PPI file not found at {path}. Returning empty graph. "
            "Download from https://datadryad.org/dataset/doi:10.5061/dryad.15dv41p84",
            UserWarning,
            stacklevel=2,
        )
        return nx.Graph()

    df = pd.read_csv(path, sep="\t")

    # Auto-detect confidence column
    conf_col = confidence_col
    if conf_col is None:
        for candidate in _CONFIDENCE_CANDIDATES:
            if candidate in df.columns:
                conf_col = candidate
                break

    # Auto-detect source column (DNS/DNU/PPI)
    src_col = None
    for candidate in _SOURCE_CANDIDATES:
        if candidate in df.columns:
            src_col = candidate
            break

    # Auto-detect gene columns
    a_col = gene_a_col or df.columns[0]
    b_col = gene_b_col or df.columns[1]

    # Filter by confidence
    if conf_col and conf_col in df.columns:
        df = df[df[conf_col] >= min_confidence]

    G = nx.Graph()
    for _, row in df.iterrows():
        a = str(row[a_col]).strip()
        b = str(row[b_col]).strip()
        if a == b:
            continue

        conf = float(row[conf_col]) if conf_col and conf_col in df.columns else 1.0
        attrs: dict = {"source": "AlphaFold_Baker2025", "confidence": conf}

        if src_col and src_col in df.columns:
            attrs["prediction_source"] = str(row[src_col])

        G.add_edge(a, b, **attrs)

    return G
