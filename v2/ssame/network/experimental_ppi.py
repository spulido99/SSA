"""Loader for experimental (unbiased) PPI networks like HI-II-14."""

from __future__ import annotations

from pathlib import Path

import networkx as nx
import pandas as pd


def load_hi_union(path: Path) -> nx.Graph:
    """Load HI-II-14 / HI-Union interactome.

    Expected format: TSV with columns Symbol_A, Symbol_B.

    Args:
        path: Path to HI-II-14 TSV file.

    Returns:
        Undirected NetworkX graph with Hugo gene symbol nodes.
    """
    df = pd.read_csv(path, sep="\t")

    G = nx.Graph()
    for _, row in df.iterrows():
        gene_a = str(row.iloc[0]).strip()
        gene_b = str(row.iloc[1]).strip()
        if gene_a != gene_b:
            G.add_edge(gene_a, gene_b, source="HI-II-14", confidence=1.0)

    return G
