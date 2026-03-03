"""Tests for ant colony search (WP5.3). Written FIRST per TDD."""

import networkx as nx
import numpy as np
import pytest

from ssame.search.ant_colony import AntColonySearch
from ssame.types import Subnetwork


@pytest.fixture
def toy_graph() -> nx.Graph:
    """Small network with known structure."""
    G = nx.Graph()
    # A cluster of known "driver" genes
    G.add_edges_from([
        ("TP53", "MDM2"), ("TP53", "CDKN2A"), ("MDM2", "CDKN2A"),
        ("PIK3CA", "PIK3R1"), ("PIK3CA", "AKT1"), ("PIK3R1", "AKT1"),
        ("TP53", "PIK3CA"),  # bridge
    ])
    # Some peripheral genes
    G.add_edges_from([
        ("NOISE1", "NOISE2"), ("NOISE2", "NOISE3"),
        ("TP53", "NOISE1"),
    ])
    return G


@pytest.fixture
def gene_scores() -> dict[str, float]:
    """High scores for known drivers, low for noise."""
    return {
        "TP53": 0.95,
        "MDM2": 0.85,
        "CDKN2A": 0.80,
        "PIK3CA": 0.90,
        "PIK3R1": 0.75,
        "AKT1": 0.70,
        "NOISE1": 0.10,
        "NOISE2": 0.05,
        "NOISE3": 0.02,
    }


class TestAntColonySearch:
    def test_creates_search(self, toy_graph: nx.Graph) -> None:
        search = AntColonySearch(graph=toy_graph, n_ants=10, n_iterations=50)
        assert search is not None

    def test_returns_subnetworks(
        self, toy_graph: nx.Graph, gene_scores: dict[str, float]
    ) -> None:
        search = AntColonySearch(
            graph=toy_graph, n_ants=10, n_iterations=50,
            subnetwork_size=3, seed=42,
        )
        results = search.search(gene_scores)
        assert isinstance(results, list)
        assert len(results) > 0
        for subnet in results:
            assert isinstance(subnet, Subnetwork)

    def test_subnetworks_are_connected(
        self, toy_graph: nx.Graph, gene_scores: dict[str, float]
    ) -> None:
        search = AntColonySearch(
            graph=toy_graph, n_ants=10, n_iterations=50,
            subnetwork_size=3, seed=42,
        )
        results = search.search(gene_scores)
        for subnet in results:
            # All genes in a subnetwork should be connected in the original graph
            sub_g = toy_graph.subgraph([g.symbol for g in subnet.genes])
            assert nx.is_connected(sub_g)

    def test_prefers_high_score_genes(
        self, toy_graph: nx.Graph, gene_scores: dict[str, float]
    ) -> None:
        search = AntColonySearch(
            graph=toy_graph, n_ants=20, n_iterations=100,
            subnetwork_size=3, seed=42,
        )
        results = search.search(gene_scores)
        # Count how often each gene appears across all discovered subnetworks
        gene_counts: dict[str, int] = {}
        for subnet in results:
            for g in subnet.genes:
                gene_counts[g.symbol] = gene_counts.get(g.symbol, 0) + 1

        # TP53 (score=0.95) should appear more than NOISE3 (score=0.02)
        if "TP53" in gene_counts and "NOISE3" in gene_counts:
            assert gene_counts["TP53"] > gene_counts["NOISE3"]

    def test_pheromone_initialization_from_scores(
        self, toy_graph: nx.Graph, gene_scores: dict[str, float]
    ) -> None:
        """Neural pheromone initialization: scores initialize pheromone levels."""
        search = AntColonySearch(
            graph=toy_graph, n_ants=10, n_iterations=10,
            subnetwork_size=3, seed=42,
        )
        search._init_pheromone(gene_scores)
        tp53_idx = search._node_to_idx["TP53"]
        noise3_idx = search._node_to_idx["NOISE3"]
        assert search._pheromone[tp53_idx] > search._pheromone[noise3_idx]

    def test_evaporation_reduces_pheromone(
        self, toy_graph: nx.Graph
    ) -> None:
        search = AntColonySearch(
            graph=toy_graph, n_ants=5, n_iterations=5,
            subnetwork_size=3, seed=42,
        )
        search._init_pheromone({"TP53": 0.9, "MDM2": 0.8})
        tp53_idx = search._node_to_idx["TP53"]
        before = search._pheromone[tp53_idx]
        search._evaporate()
        after = search._pheromone[tp53_idx]
        assert after < before

    def test_reproducible_with_seed(
        self, toy_graph: nx.Graph, gene_scores: dict[str, float]
    ) -> None:
        search1 = AntColonySearch(
            graph=toy_graph, n_ants=10, n_iterations=50,
            subnetwork_size=3, seed=42,
        )
        results1 = search1.search(gene_scores)

        search2 = AntColonySearch(
            graph=toy_graph, n_ants=10, n_iterations=50,
            subnetwork_size=3, seed=42,
        )
        results2 = search2.search(gene_scores)

        symbols1 = {frozenset(g.symbol for g in s.genes) for s in results1}
        symbols2 = {frozenset(g.symbol for g in s.genes) for s in results2}
        assert symbols1 == symbols2


class TestVectorizedPheromone:
    def test_pheromone_is_numpy_array(
        self, toy_graph: nx.Graph, gene_scores: dict[str, float]
    ) -> None:
        """After init, internal pheromone should be numpy array."""
        search = AntColonySearch(
            graph=toy_graph, n_ants=10, n_iterations=10,
            subnetwork_size=3, seed=42,
        )
        search._node_list = list(toy_graph.nodes())
        search._init_pheromone(gene_scores)
        assert isinstance(search._pheromone, np.ndarray)

    def test_vectorized_evaporation_matches_original(
        self, toy_graph: nx.Graph, gene_scores: dict[str, float]
    ) -> None:
        """Vectorized evaporation should match dict-based."""
        search = AntColonySearch(
            graph=toy_graph, n_ants=10, n_iterations=10,
            subnetwork_size=3, seed=42,
        )
        search._node_list = list(toy_graph.nodes())
        search._init_pheromone(gene_scores)

        # Capture pre-evaporation values
        before = search._pheromone.copy()

        search._evaporate()

        # Each value should be: val * evap_rate * 0.99 + initial * 0.01, clamped
        for i, node in enumerate(search._node_list):
            expected = before[i] * search.evaporation_rate
            expected = expected * 0.99 + search.initial_prob * 0.01
            expected = max(search.min_prob, expected)
            assert search._pheromone[i] == pytest.approx(expected, abs=1e-10)
