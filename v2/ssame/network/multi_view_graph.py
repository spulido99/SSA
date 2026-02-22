"""Multi-view graph construction from multiple biological networks."""

from __future__ import annotations

from dataclasses import dataclass, field

import networkx as nx
import numpy as np


@dataclass
class MultiViewGraph:
    """Multi-view graph combining multiple biological networks.

    All views share the same node set (union of genes across all views).
    Adjacency matrices are computed lazily on first access.
    """

    node_list: list[str]
    node_to_index: dict[str, int]
    view_names: list[str]
    _adjacencies: dict[str, np.ndarray] = field(repr=False)
    _graphs: dict[str, nx.Graph] = field(default_factory=dict, repr=False)

    @property
    def node_set(self) -> set[str]:
        return set(self.node_list)

    @property
    def n_nodes(self) -> int:
        return len(self.node_list)

    @property
    def n_views(self) -> int:
        return len(self.view_names)

    def get_adjacency(self, view: str) -> np.ndarray:
        """Get adjacency matrix for a specific view (built lazily)."""
        if view not in self._adjacencies:
            G = self._graphs[view]
            n = len(self.node_list)
            adj = np.zeros((n, n), dtype=np.float64)
            for u, v in G.edges():
                if u in self.node_to_index and v in self.node_to_index:
                    i = self.node_to_index[u]
                    j = self.node_to_index[v]
                    adj[i, j] = 1.0
                    adj[j, i] = 1.0
            self._adjacencies[view] = adj
        return self._adjacencies[view]

    def get_graph(self, view: str) -> nx.Graph:
        """Get NetworkX graph for a specific view.

        Uses stored original graph if available, otherwise reconstructs
        from adjacency matrix (for views created without build_multi_view).
        """
        if view in self._graphs:
            return self._graphs[view]
        # Reconstruct from adjacency
        adj = self._adjacencies[view]
        G = nx.from_numpy_array(adj)
        node_map = dict(enumerate(self.node_list))
        G = nx.relabel_nodes(G, node_map)
        self._graphs[view] = G
        return G

    def get_edge_index(self, view: str) -> np.ndarray:
        """Get edge index in COO format (2 x n_edges) for a specific view.

        Returns directed edges (both directions for each undirected edge).
        """
        adj = self.get_adjacency(view)
        rows, cols = np.nonzero(adj)
        return np.stack([rows, cols], axis=0)


def build_multi_view(views: dict[str, nx.Graph]) -> MultiViewGraph:
    """Build a multi-view graph from a dictionary of NetworkX graphs.

    Adjacency matrices are NOT pre-computed; they are built lazily
    when get_adjacency() or get_edge_index() is called (saves memory
    for large graphs like PIONEER with ~16K nodes).

    Args:
        views: Dict mapping view name -> NetworkX graph.
            All graphs should use Hugo gene symbols as node IDs.

    Returns:
        MultiViewGraph with shared node set and per-view graphs.
    """
    # Compute union of all nodes
    all_nodes: set[str] = set()
    for G in views.values():
        all_nodes.update(G.nodes())

    node_list = sorted(all_nodes)
    node_to_index = {n: i for i, n in enumerate(node_list)}

    return MultiViewGraph(
        node_list=node_list,
        node_to_index=node_to_index,
        view_names=list(views.keys()),
        _adjacencies={},
        _graphs=dict(views),
    )
