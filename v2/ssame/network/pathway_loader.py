"""Loader for pathway-based networks (Reactome, KEGG)."""

from __future__ import annotations

from itertools import combinations
from pathlib import Path

import networkx as nx


def load_reactome(path: Path) -> nx.Graph:
    """Load Reactome pathways and build a gene-gene network.

    Genes in the same pathway are connected by edges.
    Expected format: GMT (Gene Matrix Transposed) —
    each line: pathway_id<TAB>pathway_name<TAB>gene1<TAB>gene2<TAB>...

    Args:
        path: Path to Reactome GMT file.

    Returns:
        Undirected NetworkX graph where edges represent shared pathway membership.
    """
    G = nx.Graph()

    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            pathway_id = parts[0]
            pathway_name = parts[1]
            genes = parts[2:]

            for gene_a, gene_b in combinations(genes, 2):
                if G.has_edge(gene_a, gene_b):
                    existing = G[gene_a][gene_b]
                    existing.setdefault("pathways", []).append(pathway_id)
                else:
                    G.add_edge(
                        gene_a,
                        gene_b,
                        source="Reactome",
                        pathway=pathway_id,
                        pathway_name=pathway_name,
                        pathways=[pathway_id],
                        confidence=1.0,
                    )

    return G
