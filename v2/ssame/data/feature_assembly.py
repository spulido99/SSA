"""Feature assembly: builds per-gene feature matrix from multiple sources."""

from __future__ import annotations

import networkx as nx
import numpy as np

from ssame.network.multi_view_graph import MultiViewGraph
from ssame.types import AlterationMatrix


class FeatureAssembler:
    """Assembles all features into per-gene feature matrix X in R^{G x d}.

    Features are split into two modalities (SGCD-inspired bimodal encoding):
    - Omics features: mutation frequency, etc.
    - Topology features: degree, betweenness per network view.
    """

    def __init__(
        self,
        alteration_matrix: AlterationMatrix,
        multi_view: MultiViewGraph,
        additional_features: dict[str, dict[str, float]] | None = None,
    ) -> None:
        self.alteration_matrix = alteration_matrix
        self.multi_view = multi_view
        self.additional_features = additional_features or {}
        self._feature_names: list[str] = []
        self._omics_indices: list[int] = []
        self._topo_indices: list[int] = []

    @property
    def feature_names(self) -> list[str]:
        return self._feature_names

    def get_bimodal_indices(self) -> tuple[list[int], list[int]]:
        """Return (omics_feature_indices, topology_feature_indices)."""
        return self._omics_indices, self._topo_indices

    def assemble(
        self,
        *,
        standardize: bool = False,
    ) -> tuple[np.ndarray, list[str]]:
        """Assemble features for all genes.

        Args:
            standardize: If True, z-score standardize each feature column.

        Returns:
            (X, gene_list) where X is shape (n_genes, n_features).
        """
        # Use genes present in both alteration matrix and network
        matrix_genes = set(self.alteration_matrix.genes)
        network_genes = self.multi_view.node_set
        common_genes = sorted(matrix_genes & network_genes)

        if not common_genes:
            # Fallback: use matrix genes if no overlap
            common_genes = sorted(matrix_genes)

        features_list: list[np.ndarray] = []
        self._feature_names = []
        self._omics_indices = []
        self._topo_indices = []
        col_idx = 0

        # --- Omics features ---
        omics_start = col_idx

        # 1. Mutation frequency
        freq = np.array([
            self.alteration_matrix.gene_frequency(g)
            if g in matrix_genes else 0.0
            for g in common_genes
        ])
        features_list.append(freq.reshape(-1, 1))
        self._feature_names.append("mutation_frequency")
        col_idx += 1

        # 2. Log mutation frequency (useful for skewed distributions)
        log_freq = np.log1p(freq * self.alteration_matrix.n_samples)
        features_list.append(log_freq.reshape(-1, 1))
        self._feature_names.append("log_mutation_count")
        col_idx += 1

        # 3. Additional features (e.g., DepMap, pLI, etc.)
        for feat_name, feat_dict in sorted(self.additional_features.items()):
            vals = np.array([feat_dict.get(g, 0.0) for g in common_genes])
            features_list.append(vals.reshape(-1, 1))
            self._feature_names.append(feat_name)
            col_idx += 1

        self._omics_indices = list(range(omics_start, col_idx))

        # --- Topology features ---
        topo_start = col_idx

        for view_name in self.multi_view.view_names:
            G = self.multi_view.get_graph(view_name)

            # Degree
            degree_dict = dict(G.degree())
            degrees = np.array([
                degree_dict.get(g, 0) for g in common_genes
            ], dtype=np.float64)
            features_list.append(degrees.reshape(-1, 1))
            self._feature_names.append(f"{view_name}_degree")
            col_idx += 1

            # Betweenness centrality (approximate for large graphs, skip for very large)
            n_nodes = G.number_of_nodes()
            n_edges = G.number_of_edges()
            if n_edges > 500_000:
                # Too many edges for betweenness — use zeros
                betweenness = {}
            else:
                k = min(500, n_nodes) if n_nodes > 500 else None
                betweenness = nx.betweenness_centrality(G, k=k)
            betw_vals = np.array([
                betweenness.get(g, 0.0) for g in common_genes
            ])
            features_list.append(betw_vals.reshape(-1, 1))
            self._feature_names.append(f"{view_name}_betweenness")
            col_idx += 1

        self._topo_indices = list(range(topo_start, col_idx))

        X = np.hstack(features_list)

        # Handle any NaN
        X = np.nan_to_num(X, nan=0.0)

        if standardize:
            mean = X.mean(axis=0)
            std = X.std(axis=0)
            std[std == 0] = 1.0
            X = (X - mean) / std

        return X, common_genes
