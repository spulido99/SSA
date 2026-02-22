"""Hybrid search: ant colony with neural pheromone init + beam search."""

from __future__ import annotations

import networkx as nx

from ssame.search.ant_colony import AntColonySearch
from ssame.search.beam_search import BeamSearch
from ssame.types import Subnetwork


class HybridSearch:
    """Combines ant colony search (neural-initialized pheromone) with beam search.

    The neural encoder scores initialize the pheromone levels, biasing
    ant colony exploration toward high-scoring genes. Beam search results
    are merged for additional diversity.
    """

    def __init__(
        self,
        graph: nx.Graph,
        n_ants: int = 50,
        n_iterations: int = 1000,
        subnetwork_size: int = 3,
        beam_width: int = 10,
        seed: int = 42,
    ) -> None:
        self.graph = graph
        self.n_ants = n_ants
        self.n_iterations = n_iterations
        self.subnetwork_size = subnetwork_size
        self.beam_width = beam_width
        self.seed = seed

    def search(self, gene_scores: dict[str, float]) -> list[Subnetwork]:
        """Run hybrid search combining ant colony and beam search.

        Args:
            gene_scores: Per-gene driverness scores from neural encoder.

        Returns:
            Merged list of subnetworks from both search strategies.
        """
        # Run ant colony with neural pheromone initialization
        aco = AntColonySearch(
            graph=self.graph,
            n_ants=self.n_ants,
            n_iterations=self.n_iterations,
            subnetwork_size=self.subnetwork_size,
            seed=self.seed,
        )
        aco_results = aco.search(gene_scores)

        # Run beam search
        beam = BeamSearch(
            graph=self.graph,
            beam_width=self.beam_width,
            subnetwork_size=self.subnetwork_size,
        )
        beam_results = beam.search(gene_scores)

        # Merge and deduplicate
        return self._merge_results(aco_results, beam_results, gene_scores)

    @staticmethod
    def _merge_results(
        aco_results: list[Subnetwork],
        beam_results: list[Subnetwork],
        gene_scores: dict[str, float],
    ) -> list[Subnetwork]:
        """Merge results from both strategies, deduplicate, rank by score."""
        seen: dict[frozenset[str], Subnetwork] = {}

        for subnet in aco_results + beam_results:
            key = frozenset(g.symbol for g in subnet.genes)
            if key not in seen:
                seen[key] = subnet

        # Rank by total gene score
        def subnet_score(s: Subnetwork) -> float:
            return sum(gene_scores.get(g.symbol, 0.0) for g in s.genes)

        return sorted(seen.values(), key=subnet_score, reverse=True)
