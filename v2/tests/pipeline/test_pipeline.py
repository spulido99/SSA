"""Tests for pipeline orchestration (WP7.1). Written FIRST per TDD."""

import networkx as nx
import numpy as np
import pytest

from ssame.network.multi_view_graph import build_multi_view
from ssame.pipeline.pipeline import Pipeline, PipelineResult
from ssame.types import AlterationMatrix


@pytest.fixture
def toy_alteration_matrix() -> AlterationMatrix:
    genes = ["TP53", "PIK3CA", "CDH1", "KRAS", "BRAF", "MDM2", "CDKN2A"]
    samples = [f"S{i}" for i in range(30)]
    data = np.zeros((7, 30), dtype=np.int8)
    # ME pattern: TP53 and PIK3CA
    data[0, :15] = 1   # TP53
    data[1, 15:25] = 1  # PIK3CA
    data[2, :5] = 1     # CDH1
    data[3, 10:16] = 1  # KRAS
    data[4, 20:23] = 1  # BRAF
    data[5, 2:8] = 1    # MDM2
    data[6, 8:12] = 1   # CDKN2A
    return AlterationMatrix(data=data, genes=genes, samples=samples)


@pytest.fixture
def toy_network():
    G = nx.Graph()
    G.add_edges_from([
        ("TP53", "MDM2"), ("TP53", "CDKN2A"), ("MDM2", "CDKN2A"),
        ("PIK3CA", "KRAS"), ("KRAS", "BRAF"),
        ("TP53", "PIK3CA"), ("CDH1", "PIK3CA"),
    ])
    return build_multi_view({"ppi": G})


@pytest.fixture
def toy_labels() -> dict[str, int]:
    return {
        "TP53": 1,
        "PIK3CA": 1,
        "CDH1": 0,
        "KRAS": 0,
        "BRAF": 0,
        "MDM2": 0,
        "CDKN2A": 0,
    }


class TestPipeline:
    def test_creates_pipeline(
        self, toy_alteration_matrix, toy_network, toy_labels
    ) -> None:
        pipeline = Pipeline(
            alteration_matrix=toy_alteration_matrix,
            multi_view=toy_network,
            labels=toy_labels,
        )
        assert pipeline is not None

    def test_run_returns_result(
        self, toy_alteration_matrix, toy_network, toy_labels
    ) -> None:
        pipeline = Pipeline(
            alteration_matrix=toy_alteration_matrix,
            multi_view=toy_network,
            labels=toy_labels,
            encoder_n_estimators=10,
            search_n_ants=5,
            search_n_iterations=20,
        )
        result = pipeline.run()
        assert isinstance(result, PipelineResult)

    def test_result_has_ranked_genes(
        self, toy_alteration_matrix, toy_network, toy_labels
    ) -> None:
        pipeline = Pipeline(
            alteration_matrix=toy_alteration_matrix,
            multi_view=toy_network,
            labels=toy_labels,
            encoder_n_estimators=10,
            search_n_ants=5,
            search_n_iterations=20,
        )
        result = pipeline.run()
        assert len(result.ranked_genes) > 0
        # Should be a list of (gene, score) tuples
        for gene, score in result.ranked_genes:
            assert isinstance(gene, str)
            assert isinstance(score, float)

    def test_result_has_modules(
        self, toy_alteration_matrix, toy_network, toy_labels
    ) -> None:
        pipeline = Pipeline(
            alteration_matrix=toy_alteration_matrix,
            multi_view=toy_network,
            labels=toy_labels,
            encoder_n_estimators=10,
            search_n_ants=5,
            search_n_iterations=20,
        )
        result = pipeline.run()
        assert len(result.modules) > 0

    def test_known_drivers_rank_high(
        self, toy_alteration_matrix, toy_network, toy_labels
    ) -> None:
        pipeline = Pipeline(
            alteration_matrix=toy_alteration_matrix,
            multi_view=toy_network,
            labels=toy_labels,
            encoder_n_estimators=50,
            search_n_ants=10,
            search_n_iterations=50,
        )
        result = pipeline.run()
        top_genes = [g for g, _ in result.ranked_genes[:5]]
        # At least one known driver should appear in top 5
        assert "TP53" in top_genes or "PIK3CA" in top_genes

    def test_search_graph_restricted_to_scored_genes(
        self, toy_alteration_matrix, toy_labels
    ) -> None:
        """Graph with extra nodes not in alteration matrix should be pruned."""
        G = nx.Graph()
        # Genes in alteration matrix
        G.add_edges_from([
            ("TP53", "MDM2"), ("TP53", "CDKN2A"), ("MDM2", "CDKN2A"),
            ("PIK3CA", "KRAS"), ("KRAS", "BRAF"),
            ("TP53", "PIK3CA"), ("CDH1", "PIK3CA"),
        ])
        # 50 extra nodes NOT in alteration matrix (simulates large PPI)
        for i in range(50):
            G.add_edge(f"EXTRA_{i}", f"EXTRA_{i + 1}")
            G.add_edge(f"EXTRA_{i}", "TP53")
        multi_view = build_multi_view({"ppi": G})

        pipeline = Pipeline(
            alteration_matrix=toy_alteration_matrix,
            multi_view=multi_view,
            labels=toy_labels,
            encoder_n_estimators=10,
            search_n_ants=5,
            search_n_iterations=20,
        )
        result = pipeline.run()

        # Modules should only contain genes from alteration matrix
        for module in result.modules:
            for gene in module.genes:
                assert gene.symbol in toy_alteration_matrix.genes, (
                    f"Module contains gene {gene.symbol} not in alteration matrix"
                )
