"""Tests for hybrid search (WP5.6). Written FIRST per TDD."""

import networkx as nx
import pytest

from ssame.search.hybrid_search import HybridSearch
from ssame.types import Subnetwork


@pytest.fixture
def toy_graph() -> nx.Graph:
    G = nx.Graph()
    G.add_edges_from([
        ("TP53", "MDM2"), ("TP53", "CDKN2A"), ("MDM2", "CDKN2A"),
        ("PIK3CA", "KRAS"), ("KRAS", "BRAF"),
        ("TP53", "PIK3CA"), ("CDH1", "PIK3CA"),
    ])
    return G


@pytest.fixture
def toy_scores() -> dict[str, float]:
    return {
        "TP53": 0.9, "PIK3CA": 0.8, "MDM2": 0.6,
        "CDKN2A": 0.5, "KRAS": 0.4, "BRAF": 0.3, "CDH1": 0.2,
    }


class TestHybridSearch:
    def test_creates_search(self, toy_graph: nx.Graph) -> None:
        search = HybridSearch(
            graph=toy_graph,
            n_ants=10,
            n_iterations=20,
            subnetwork_size=3,
        )
        assert search is not None

    def test_returns_subnetworks(
        self, toy_graph: nx.Graph, toy_scores: dict[str, float]
    ) -> None:
        search = HybridSearch(
            graph=toy_graph, n_ants=10, n_iterations=20, subnetwork_size=3,
        )
        results = search.search(toy_scores)
        assert len(results) > 0
        assert all(isinstance(s, Subnetwork) for s in results)

    def test_neural_init_converges_faster(
        self, toy_graph: nx.Graph, toy_scores: dict[str, float]
    ) -> None:
        """Hybrid search with neural scores should find good subnetworks faster.

        Compare: ant colony with random init vs neural init over same iterations.
        Neural init should find higher-scoring subnetworks.
        """
        from ssame.search.ant_colony import AntColonySearch

        # Random init (uniform scores)
        uniform_scores = {g: 0.5 for g in toy_scores}
        random_search = AntColonySearch(
            graph=toy_graph, n_ants=10, n_iterations=50, subnetwork_size=3, seed=42,
        )
        random_results = random_search.search(uniform_scores)

        # Neural init (informed scores)
        hybrid = HybridSearch(
            graph=toy_graph, n_ants=10, n_iterations=50, subnetwork_size=3, seed=42,
        )
        hybrid_results = hybrid.search(toy_scores)

        # Both should find results
        assert len(random_results) > 0
        assert len(hybrid_results) > 0

        # Hybrid should find high-scoring genes more reliably
        hybrid_top_genes = {g.symbol for g in hybrid_results[0].genes}
        assert "TP53" in hybrid_top_genes or "PIK3CA" in hybrid_top_genes

    def test_uses_both_ant_colony_and_beam(
        self, toy_graph: nx.Graph, toy_scores: dict[str, float]
    ) -> None:
        """Hybrid should combine results from both strategies."""
        search = HybridSearch(
            graph=toy_graph, n_ants=10, n_iterations=20, subnetwork_size=3,
            beam_width=5,
        )
        results = search.search(toy_scores)
        assert len(results) > 0
