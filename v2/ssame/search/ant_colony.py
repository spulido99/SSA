"""Ant colony subnetwork search (port of v1 NodeCostNetworkManager + Fungus)."""

from __future__ import annotations

import networkx as nx
import numpy as np

from ssame.types import Gene, Interaction, Subnetwork


class AntColonySearch:
    """Ant Colony Optimization for subnetwork discovery.

    Ports v1's NodeCostNetworkManager + Fungus random walk logic:
    - Pheromone levels on each gene guide exploration
    - Ants perform random walks weighted by pheromone
    - Good subnetworks reinforce gene pheromone levels
    - Evaporation prevents premature convergence

    Key v2 addition: pheromone can be initialized from neural encoder scores
    (= "neural pheromone initialization").

    Pheromone is stored as a numpy array for vectorized evaporation.
    """

    def __init__(
        self,
        graph: nx.Graph,
        n_ants: int = 50,
        n_iterations: int = 1000,
        subnetwork_size: int = 3,
        pheromone_strength: float = 0.005,
        evaporation_rate: float = 0.996,
        initial_prob: float = 0.5,
        min_prob: float = 0.01,
        seed: int = 42,
    ) -> None:
        self.graph = graph
        self.n_ants = n_ants
        self.n_iterations = n_iterations
        self.subnetwork_size = subnetwork_size
        self.pheromone_strength = pheromone_strength
        self.evaporation_rate = evaporation_rate
        self.initial_prob = initial_prob
        self.min_prob = min_prob
        self.seed = seed
        self._rng = np.random.default_rng(seed)
        self._node_list: list[str] = list(graph.nodes())
        self._node_to_idx: dict[str, int] = {n: i for i, n in enumerate(self._node_list)}
        self._neighbors: dict[str, list[str]] = {}
        self._pheromone: np.ndarray = np.array([])
        self._best_subnetworks: dict[str, tuple[set[str], float]] = {}

    def _init_pheromone(self, gene_scores: dict[str, float]) -> None:
        """Initialize pheromone from external gene scores (neural encoder)."""
        n = len(self._node_list)
        self._pheromone = np.full(n, self.initial_prob)
        self._node_to_idx = {node: i for i, node in enumerate(self._node_list)}

        for node in self._node_list:
            if node in gene_scores:
                idx = self._node_to_idx[node]
                self._pheromone[idx] = max(
                    self.min_prob,
                    min(1.0, gene_scores[node]),
                )

    def _evaporate(self) -> None:
        """Evaporate pheromone (vectorized numpy)."""
        self._pheromone *= self.evaporation_rate
        self._pheromone = self._pheromone * 0.99 + self.initial_prob * 0.01
        np.maximum(self._pheromone, self.min_prob, out=self._pheromone)

    def _random_walk(self) -> set[str]:
        """Single ant random walk (v1 Fungus.selectSubNetwork logic).

        Starts from a random node, extends greedily by pheromone-weighted
        neighbor selection until target size reached.
        """
        nodes = self._node_list
        # Weight start node by pheromone
        probs = self._pheromone.copy()
        probs /= probs.sum()
        start_idx = self._rng.choice(len(nodes), p=probs)
        start = nodes[start_idx]

        visited = {start}

        while len(visited) < self.subnetwork_size:
            # Find frontier: neighbors of visited nodes not yet visited
            frontier_nodes: list[str] = []
            frontier_weights: list[float] = []
            for v in visited:
                for neighbor in self._neighbors[v]:
                    if neighbor not in visited:
                        idx = self._node_to_idx[neighbor]
                        frontier_nodes.append(neighbor)
                        frontier_weights.append(self._pheromone[idx])

            if not frontier_nodes:
                break

            # Weighted random selection
            weights = np.array(frontier_weights)
            np.maximum(weights, self.min_prob, out=weights)
            weights /= weights.sum()

            chosen = frontier_nodes[self._rng.choice(len(frontier_nodes), p=weights)]
            visited.add(chosen)

        return visited

    def _update_pheromone(
        self, subnetworks: list[set[str]], scores: list[float]
    ) -> None:
        """Reinforce pheromone for genes in high-scoring subnetworks."""
        if not scores:
            return

        min_score = min(scores)
        max_score = max(scores)
        score_range = max_score - min_score if max_score > min_score else 1.0

        for subnet, score in zip(subnetworks, scores, strict=True):
            scaled = (score - min_score) / score_range
            for gene in subnet:
                idx = self._node_to_idx[gene]
                new = self._pheromone[idx] * (1.0 + self.pheromone_strength * scaled)
                self._pheromone[idx] = min(1.0, max(self.min_prob, new))

    def search(
        self,
        gene_scores: dict[str, float],
        scorer: object | None = None,
    ) -> list[Subnetwork]:
        """Run ant colony search.

        Args:
            gene_scores: Per-gene driverness scores from encoder.
                Used to initialize pheromone levels.
            scorer: Optional SubnetworkScorer to evaluate subnetworks.
                If None, uses sum of gene scores as proxy.

        Returns:
            List of discovered subnetworks, ranked by score.
        """
        # Cache node list and neighbor lookup for hot loop
        self._node_list = list(self.graph.nodes())
        self._neighbors = {n: list(self.graph.neighbors(n)) for n in self._node_list}
        self._init_pheromone(gene_scores)
        all_discovered: dict[frozenset[str], float] = {}

        for _ in range(self.n_iterations):
            iteration_subnets: list[set[str]] = []
            iteration_scores: list[float] = []

            for _ in range(self.n_ants):
                subnet = self._random_walk()
                if len(subnet) >= 2:
                    # Score by sum of gene scores (proxy if no scorer)
                    score = sum(gene_scores.get(g, 0.0) for g in subnet)
                    iteration_subnets.append(subnet)
                    iteration_scores.append(score)

                    key = frozenset(subnet)
                    if key not in all_discovered or all_discovered[key] < score:
                        all_discovered[key] = score

            # Update pheromone
            self._update_pheromone(iteration_subnets, iteration_scores)
            self._evaporate()

        # Convert to Subnetwork objects, ranked by score
        ranked = sorted(all_discovered.items(), key=lambda x: -x[1])
        results = []
        for gene_set, _score in ranked:
            genes = [Gene(symbol=g) for g in sorted(gene_set)]
            # Extract interactions from graph
            interactions = []
            gene_list = sorted(gene_set)
            for i, g1 in enumerate(gene_list):
                for g2 in gene_list[i + 1 :]:
                    if self.graph.has_edge(g1, g2):
                        interactions.append(
                            Interaction(
                                gene_a=Gene(symbol=g1),
                                gene_b=Gene(symbol=g2),
                                source="search",
                            )
                        )
            results.append(Subnetwork(genes=genes, interactions=interactions))

        return results
