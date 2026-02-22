"""Tests for multi-view graph builder (WP2.8). Written FIRST per TDD."""

import networkx as nx
import numpy as np
import pytest

from ssame.network.multi_view_graph import MultiViewGraph, build_multi_view


def _make_graph(edges: list[tuple[str, str]], source: str) -> nx.Graph:
    G = nx.Graph()
    for a, b in edges:
        G.add_edge(a, b, source=source, confidence=1.0)
    return G


@pytest.fixture
def two_views() -> dict[str, nx.Graph]:
    view1 = _make_graph(
        [("TP53", "MDM2"), ("BRCA1", "BARD1"), ("TP53", "CDKN2A")],
        source="HI-II-14",
    )
    view2 = _make_graph(
        [("TP53", "MDM2"), ("KRAS", "BRAF"), ("PIK3CA", "PIK3R1")],
        source="Reactome",
    )
    return {"hi_union": view1, "reactome": view2}


class TestMultiViewGraph:
    def test_build_from_dict(self, two_views: dict[str, nx.Graph]) -> None:
        mv = build_multi_view(two_views)
        assert isinstance(mv, MultiViewGraph)

    def test_shared_node_set(self, two_views: dict[str, nx.Graph]) -> None:
        mv = build_multi_view(two_views)
        # Union of all genes across views
        expected_nodes = {
            "TP53", "MDM2", "BRCA1", "BARD1", "CDKN2A",
            "KRAS", "BRAF", "PIK3CA", "PIK3R1",
        }
        assert mv.node_set == expected_nodes

    def test_view_names(self, two_views: dict[str, nx.Graph]) -> None:
        mv = build_multi_view(two_views)
        assert set(mv.view_names) == {"hi_union", "reactome"}

    def test_per_view_adjacency(self, two_views: dict[str, nx.Graph]) -> None:
        mv = build_multi_view(two_views)
        adj = mv.get_adjacency("hi_union")
        assert isinstance(adj, np.ndarray)
        assert adj.shape[0] == adj.shape[1] == len(mv.node_set)

    def test_adjacency_is_symmetric(self, two_views: dict[str, nx.Graph]) -> None:
        mv = build_multi_view(two_views)
        for view in mv.view_names:
            adj = mv.get_adjacency(view)
            np.testing.assert_array_equal(adj, adj.T)

    def test_adjacency_values_are_binary(self, two_views: dict[str, nx.Graph]) -> None:
        mv = build_multi_view(two_views)
        adj = mv.get_adjacency("hi_union")
        assert set(np.unique(adj)).issubset({0.0, 1.0})

    def test_node_index_mapping(self, two_views: dict[str, nx.Graph]) -> None:
        mv = build_multi_view(two_views)
        idx = mv.node_to_index
        assert "TP53" in idx
        assert len(idx) == len(mv.node_set)
        # All indices are unique
        assert len(set(idx.values())) == len(idx)

    def test_edge_in_view(self, two_views: dict[str, nx.Graph]) -> None:
        mv = build_multi_view(two_views)
        adj = mv.get_adjacency("hi_union")
        idx = mv.node_to_index
        # TP53-MDM2 should be connected in hi_union
        assert adj[idx["TP53"], idx["MDM2"]] == 1.0
        # KRAS-BRAF should NOT be in hi_union
        assert adj[idx["KRAS"], idx["BRAF"]] == 0.0

    def test_n_nodes_property(self, two_views: dict[str, nx.Graph]) -> None:
        mv = build_multi_view(two_views)
        assert mv.n_nodes == 9

    def test_n_views_property(self, two_views: dict[str, nx.Graph]) -> None:
        mv = build_multi_view(two_views)
        assert mv.n_views == 2

    def test_get_edge_index(self, two_views: dict[str, nx.Graph]) -> None:
        mv = build_multi_view(two_views)
        edge_index = mv.get_edge_index("hi_union")
        assert edge_index.shape[0] == 2
        # 3 undirected edges -> 6 directed edges
        assert edge_index.shape[1] == 6
