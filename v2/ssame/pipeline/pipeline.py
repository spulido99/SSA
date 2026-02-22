"""Pipeline orchestration: Encoder -> Search -> Score -> Output."""

from __future__ import annotations

from dataclasses import dataclass, field

import networkx as nx
import numpy as np
import torch

from ssame.data.feature_assembly import FeatureAssembler
from ssame.encoder.encoder import SSAEncoder
from ssame.encoder.tabular_baseline import TabularBaseline
from ssame.network.multi_view_graph import MultiViewGraph
from ssame.scorer.me_scorer import MEScorer
from ssame.search.ant_colony import AntColonySearch
from ssame.types import AlterationMatrix, Subnetwork


@dataclass
class PipelineResult:
    """Output of a pipeline run."""

    ranked_genes: list[tuple[str, float]]
    modules: list[Subnetwork]
    gene_scores: dict[str, float] = field(default_factory=dict)
    module_scores: list[float] = field(default_factory=list)


class Pipeline:
    """Full SSA-ME v2 pipeline: features -> encode -> search -> score -> rank.

    Supports two encoder levels:
      Level 0: XGBoost (tabular baseline, no graph)
      Level 1: Graph Transformer (bimodal encoder + graph attention)
    """

    def __init__(
        self,
        alteration_matrix: AlterationMatrix,
        multi_view: MultiViewGraph,
        labels: dict[str, int],
        *,
        encoder_level: int = 0,
        encoder_n_estimators: int = 100,
        encoder_hidden_dim: int = 64,
        encoder_n_epochs: int = 100,
        search_n_ants: int = 50,
        search_n_iterations: int = 1000,
        search_subnetwork_size: int = 3,
        seed: int = 42,
    ) -> None:
        self.alteration_matrix = alteration_matrix
        self.multi_view = multi_view
        self.labels = labels
        self.encoder_level = encoder_level
        self.encoder_n_estimators = encoder_n_estimators
        self.encoder_hidden_dim = encoder_hidden_dim
        self.encoder_n_epochs = encoder_n_epochs
        self.search_n_ants = search_n_ants
        self.search_n_iterations = search_n_iterations
        self.search_subnetwork_size = search_subnetwork_size
        self.seed = seed

    def run(self) -> PipelineResult:
        """Execute the full pipeline."""

        # Step 1: Feature assembly
        assembler = FeatureAssembler(
            alteration_matrix=self.alteration_matrix,
            multi_view=self.multi_view,
        )
        X, gene_list = assembler.assemble(standardize=True)

        # Step 2: Prepare labels (PU: 1=known driver, 0=unlabeled)
        y = np.array([self.labels.get(g, 0) for g in gene_list])

        # Step 3: Train encoder
        if self.encoder_level == 0:
            scores = self._run_l0(X, y)
        elif self.encoder_level == 1:
            omics_idx, topo_idx = assembler.get_bimodal_indices()
            scores = self._run_l1(X, y, gene_list, omics_idx, topo_idx)
        else:
            msg = f"Unknown encoder level: {self.encoder_level}"
            raise ValueError(msg)

        gene_scores = {g: float(s) for g, s in zip(gene_list, scores, strict=True)}

        # Step 4: Build search graph from first network view,
        # restricted to genes with encoder scores (huge speedup).
        first_view = self.multi_view.view_names[0]
        search_graph = self.multi_view.get_graph(first_view)

        scored_nodes = set(gene_scores.keys()) & set(search_graph.nodes())
        if scored_nodes:
            search_graph = search_graph.subgraph(scored_nodes).copy()

        # Step 5: Ant colony search
        search = AntColonySearch(
            graph=search_graph,
            n_ants=self.search_n_ants,
            n_iterations=self.search_n_iterations,
            subnetwork_size=self.search_subnetwork_size,
            seed=self.seed,
        )
        modules = search.search(gene_scores)

        # Step 6: Score modules with ME scorer
        scorer = MEScorer(alteration_matrix=self.alteration_matrix)
        module_scores = [scorer.score(m) for m in modules]

        # Step 7: Rank genes by encoder scores
        ranked_genes = sorted(gene_scores.items(), key=lambda x: -x[1])

        return PipelineResult(
            ranked_genes=ranked_genes,
            modules=modules,
            gene_scores=gene_scores,
            module_scores=module_scores,
        )

    def _run_l0(self, X: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Level 0: XGBoost tabular baseline."""
        encoder = TabularBaseline(
            n_estimators=self.encoder_n_estimators,
            seed=self.seed,
        )
        encoder.fit(X, y)
        return encoder.predict(X)

    def _run_l1(
        self,
        X: np.ndarray,
        y: np.ndarray,
        gene_list: list[str],
        omics_idx: list[int],
        topo_idx: list[int],
    ) -> np.ndarray:
        """Level 1: Graph Transformer encoder."""
        omics_dim = len(omics_idx)
        topo_dim = len(topo_idx)

        omics_features = torch.from_numpy(X[:, omics_idx]).float()
        topo_features = torch.from_numpy(X[:, topo_idx]).float()
        labels_tensor = torch.from_numpy(y).float()

        # Build edge indices for each view
        edge_indices = []
        for view_name in self.multi_view.view_names:
            edge_index = self.multi_view.get_edge_index(view_name)
            edge_indices.append(torch.from_numpy(edge_index).long())

        # Build alteration tensor for biological priors
        alt_data = np.zeros(
            (len(gene_list), len(self.alteration_matrix.samples)),
            dtype=np.float32,
        )
        for g_idx, gene in enumerate(gene_list):
            if gene in self.alteration_matrix.genes:
                src_idx = self.alteration_matrix.genes.index(gene)
                alt_data[g_idx, :] = self.alteration_matrix.data[src_idx, :]
        alt_tensor = torch.from_numpy(alt_data)

        # Create and train encoder
        encoder = SSAEncoder(
            omics_dim=omics_dim,
            topo_dim=topo_dim,
            n_views=len(self.multi_view.view_names),
            hidden_dim=self.encoder_hidden_dim,
            alteration_matrix=alt_tensor,
        )

        optimizer = torch.optim.Adam(encoder.parameters(), lr=0.01)

        encoder.train()
        for _epoch in range(self.encoder_n_epochs):
            scores = encoder(omics_features, topo_features, edge_indices)
            loss = encoder.compute_loss(scores, labels_tensor)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

        # Predict
        encoder.eval()
        with torch.no_grad():
            scores = encoder(omics_features, topo_features, edge_indices)

        return scores.numpy()
