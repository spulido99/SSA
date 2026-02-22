"""Loader for PIONEER oncoPPI structural interaction predictions.

Data source: Xiong et al., Nature Biotechnology 2024 (doi:10.1038/s41587-024-02428-4)
Download: https://pioneer.yulab.org/downloads
File: predictions/very_high/human.txt (166 MB)
"""

from __future__ import annotations

from pathlib import Path

import networkx as nx
import pandas as pd

_CONFIDENCE_CANDIDATES = (
    "confidence", "pioneer_score", "score", "probability",
)


def load_pioneer(
    path: Path,
    *,
    min_confidence: float | None = None,
    gene_a_col: str | None = None,
    gene_b_col: str | None = None,
    confidence_col: str | None = None,
) -> nx.Graph:
    """Load PIONEER structural interaction predictions.

    Args:
        path: Path to PIONEER prediction file (TSV).
        min_confidence: Optional confidence threshold for filtering.
        gene_a_col: Column name for gene/protein A (auto-detected if None).
        gene_b_col: Column name for gene/protein B (auto-detected if None).
        confidence_col: Column name for confidence score (auto-detected if None).

    Returns:
        Undirected NetworkX graph with gene symbol nodes.
        Edge attributes include interface residue data when available.
    """
    if not path.exists():
        import warnings

        warnings.warn(
            f"PIONEER file not found at {path}. Returning empty graph. "
            "Download from https://pioneer.yulab.org/downloads",
            UserWarning,
            stacklevel=2,
        )
        return nx.Graph()

    df = pd.read_csv(path, sep="\t")

    # Auto-detect gene columns
    a_col = gene_a_col or df.columns[0]
    b_col = gene_b_col or df.columns[1]

    # Auto-detect confidence column
    conf_col = confidence_col
    if conf_col is None:
        for candidate in _CONFIDENCE_CANDIDATES:
            if candidate in df.columns:
                conf_col = candidate
                break

    # Filter by confidence if threshold given
    if min_confidence is not None and conf_col and conf_col in df.columns:
        df = df[df[conf_col] >= min_confidence]

    G = nx.Graph()
    for _, row in df.iterrows():
        a = str(row[a_col]).strip()
        b = str(row[b_col]).strip()
        if a == b:
            continue

        conf = float(row[conf_col]) if conf_col and conf_col in df.columns else 1.0
        attrs: dict = {"source": "PIONEER", "confidence": conf}

        # Preserve all additional columns as edge attributes
        for col in df.columns:
            if col not in (a_col, b_col, conf_col):
                attrs[col] = row[col]

        G.add_edge(a, b, **attrs)

    return G
