"""Loader for literature-biased PPI networks (STRING, BioGRID) — for comparison only."""

from __future__ import annotations

from pathlib import Path

import networkx as nx
import pandas as pd


def _load_string_aliases(aliases_path: Path) -> dict[str, str]:
    """Load STRING protein aliases to map ENSP IDs to Hugo symbols.

    Returns dict mapping ENSP ID -> Hugo symbol.
    """
    ensp_to_hugo: dict[str, str] = {}
    df = pd.read_csv(aliases_path, sep="\t", comment="#", header=None,
                      names=["string_protein_id", "alias", "source"])
    # Prefer HGNC aliases
    hgnc = df[df["source"].str.contains("HGNC", case=False, na=False)]
    for _, row in hgnc.iterrows():
        ensp_to_hugo[row["string_protein_id"]] = row["alias"]
    return ensp_to_hugo


def load_string(
    links_path: Path,
    aliases_path: Path,
    min_score: int = 0,
) -> nx.Graph:
    """Load STRING PPI network with Hugo symbol mapping.

    Args:
        links_path: Path to STRING protein links file.
        aliases_path: Path to STRING protein aliases file.
        min_score: Minimum combined score (0-1000). Default: include all.

    Returns:
        Undirected NetworkX graph with Hugo gene symbol nodes.
    """
    ensp_to_hugo = _load_string_aliases(aliases_path)

    df = pd.read_csv(links_path, sep="\t")

    if min_score > 0:
        df = df[df["combined_score"] >= min_score]

    G = nx.Graph()
    for _, row in df.iterrows():
        protein1 = row["protein1"]
        protein2 = row["protein2"]
        score = row["combined_score"]

        hugo_a = ensp_to_hugo.get(protein1)
        hugo_b = ensp_to_hugo.get(protein2)

        if hugo_a and hugo_b and hugo_a != hugo_b:
            G.add_edge(
                hugo_a,
                hugo_b,
                source="STRING",
                confidence=score / 1000.0,
                raw_score=int(score),
            )

    return G
