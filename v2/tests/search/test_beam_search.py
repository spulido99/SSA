"""Tests for beam search (WP5.5). Written FIRST per TDD."""

import networkx as nx
import pytest

from ssame.search.beam_search import BeamSearch
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
        "TP53": 0.9,
        "PIK3CA": 0.8,
        "MDM2": 0.6,
        "CDKN2A": 0.5,
        "KRAS": 0.4,
        "BRAF": 0.3,
        "CDH1": 0.2,
    }


class TestBeamSearch:
    def test_creates_search(self, toy_graph: nx.Graph) -> None:
        search = BeamSearch(graph=toy_graph, beam_width=5, subnetwork_size=3)
        assert search is not None

    def test_returns_subnetworks(
        self, toy_graph: nx.Graph, toy_scores: dict[str, float]
    ) -> None:
        search = BeamSearch(graph=toy_graph, beam_width=5, subnetwork_size=3)
        results = search.search(toy_scores)
        assert len(results) > 0
        assert all(isinstance(s, Subnetwork) for s in results)

    def test_subnetwork_size(
        self, toy_graph: nx.Graph, toy_scores: dict[str, float]
    ) -> None:
        search = BeamSearch(graph=toy_graph, beam_width=5, subnetwork_size=3)
        results = search.search(toy_scores)
        for subnet in results:
            assert len(subnet.genes) >= 2
            assert len(subnet.genes) <= 3

    def test_top_results_contain_high_scorers(
        self, toy_graph: nx.Graph, toy_scores: dict[str, float]
    ) -> None:
        """Top subnetworks should contain high-scoring genes."""
        search = BeamSearch(graph=toy_graph, beam_width=10, subnetwork_size=3)
        results = search.search(toy_scores)
        top_genes = {g.symbol for g in results[0].genes}
        # TP53 should be in the top result (highest scoring gene)
        assert "TP53" in top_genes

    def test_beam_width_affects_diversity(
        self, toy_graph: nx.Graph, toy_scores: dict[str, float]
    ) -> None:
        """Wider beam should produce more diverse results."""
        narrow = BeamSearch(graph=toy_graph, beam_width=2, subnetwork_size=3)
        wide = BeamSearch(graph=toy_graph, beam_width=10, subnetwork_size=3)
        results_narrow = narrow.search(toy_scores)
        results_wide = wide.search(toy_scores)
        assert len(results_wide) >= len(results_narrow)

    def test_top_k_limit(
        self, toy_graph: nx.Graph, toy_scores: dict[str, float]
    ) -> None:
        search = BeamSearch(
            graph=toy_graph, beam_width=10, subnetwork_size=3, top_k=3
        )
        results = search.search(toy_scores)
        assert len(results) <= 3

    def test_results_are_connected(
        self, toy_graph: nx.Graph, toy_scores: dict[str, float]
    ) -> None:
        """Each subnetwork should be a connected subgraph."""
        search = BeamSearch(graph=toy_graph, beam_width=5, subnetwork_size=3)
        results = search.search(toy_scores)
        for subnet in results:
            gene_set = {g.symbol for g in subnet.genes}
            subgraph = toy_graph.subgraph(gene_set)
            assert nx.is_connected(subgraph)
