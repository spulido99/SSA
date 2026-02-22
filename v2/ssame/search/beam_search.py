"""Beam search for subnetwork discovery."""

from __future__ import annotations

import networkx as nx

from ssame.types import Gene, Interaction, Subnetwork


class BeamSearch:
    """Greedy beam search for high-scoring connected subnetworks.

    Starts from each node, extends by adding highest-scoring neighbors,
    keeping top-k candidates at each step (beam width).
    """

    def __init__(
        self,
        graph: nx.Graph,
        beam_width: int = 10,
        subnetwork_size: int = 3,
        top_k: int = 50,
    ) -> None:
        self.graph = graph
        self.beam_width = beam_width
        self.subnetwork_size = subnetwork_size
        self.top_k = top_k

    def search(self, gene_scores: dict[str, float]) -> list[Subnetwork]:
        """Run beam search from all seed nodes.

        Args:
            gene_scores: Per-gene driverness scores.

        Returns:
            List of subnetworks ranked by total score, deduplicated.
        """
        all_discovered: dict[frozenset[str], float] = {}

        for seed in self.graph.nodes():
            candidates = self._beam_from_seed(seed, gene_scores)
            for gene_set, score in candidates:
                key = frozenset(gene_set)
                if key not in all_discovered or all_discovered[key] < score:
                    all_discovered[key] = score

        # Rank and convert
        ranked = sorted(all_discovered.items(), key=lambda x: -x[1])
        results = []
        for gene_set, _score in ranked[: self.top_k]:
            results.append(self._make_subnetwork(gene_set))
        return results

    def _beam_from_seed(
        self,
        seed: str,
        gene_scores: dict[str, float],
    ) -> list[tuple[set[str], float]]:
        """Beam search expanding from a single seed node."""
        initial_score = gene_scores.get(seed, 0.0)
        beam: list[tuple[set[str], float]] = [({seed}, initial_score)]

        for _ in range(self.subnetwork_size - 1):
            candidates: list[tuple[set[str], float]] = []
            for current_set, current_score in beam:
                # Find frontier: unvisited neighbors
                for node in current_set:
                    for neighbor in self.graph.neighbors(node):
                        if neighbor not in current_set:
                            new_set = current_set | {neighbor}
                            new_score = current_score + gene_scores.get(neighbor, 0.0)
                            candidates.append((new_set, new_score))

            if not candidates:
                break

            # Deduplicate and keep top beam_width
            deduped: dict[frozenset[str], tuple[set[str], float]] = {}
            for gene_set, score in candidates:
                key = frozenset(gene_set)
                if key not in deduped or deduped[key][1] < score:
                    deduped[key] = (gene_set, score)

            sorted_candidates = sorted(deduped.values(), key=lambda x: -x[1])
            beam = sorted_candidates[: self.beam_width]

        return beam

    def _make_subnetwork(self, gene_set: frozenset[str]) -> Subnetwork:
        """Convert a gene set to a Subnetwork object."""
        gene_list = sorted(gene_set)
        genes = [Gene(symbol=g) for g in gene_list]
        interactions = []
        for i, g1 in enumerate(gene_list):
            for g2 in gene_list[i + 1 :]:
                if self.graph.has_edge(g1, g2):
                    interactions.append(
                        Interaction(
                            gene_a=Gene(symbol=g1),
                            gene_b=Gene(symbol=g2),
                            source="beam_search",
                        )
                    )
        return Subnetwork(genes=genes, interactions=interactions)
