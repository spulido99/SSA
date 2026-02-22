"""Tests for Level 1 (Graph Transformer) pipeline integration. Written FIRST per TDD."""

import networkx as nx
import numpy as np
import pytest

from ssame.network.multi_view_graph import build_multi_view
from ssame.pipeline.pipeline import Pipeline, PipelineResult
from ssame.types import AlterationMatrix


@pytest.fixture
def toy_data():
    genes = ["TP53", "PIK3CA", "CDH1", "KRAS", "BRAF", "MDM2", "CDKN2A"]
    samples = [f"S{i}" for i in range(30)]
    data = np.zeros((7, 30), dtype=np.int8)
    data[0, :15] = 1
    data[1, 15:25] = 1
    data[2, :8] = 1
    data[3, 5:12] = 1
    data[4, 20:27] = 1
    data[5, :3] = 1
    data[6, 25:30] = 1

    alt_matrix = AlterationMatrix(data=data, genes=genes, samples=samples)

    graph = nx.Graph()
    graph.add_edges_from([
        ("TP53", "MDM2"), ("TP53", "CDKN2A"), ("MDM2", "CDKN2A"),
        ("PIK3CA", "KRAS"), ("KRAS", "BRAF"),
        ("TP53", "PIK3CA"), ("CDH1", "PIK3CA"),
    ])
    multi_view = build_multi_view({"ppi": graph})

    labels = {"TP53": 1, "PIK3CA": 1}

    return alt_matrix, multi_view, labels


class TestPipelineL1:
    def test_pipeline_with_encoder_level_1(self, toy_data) -> None:
        """Pipeline should work with encoder_level=1 (Graph Transformer)."""
        alt, mv, labels = toy_data
        pipe = Pipeline(
            alteration_matrix=alt,
            multi_view=mv,
            labels=labels,
            encoder_level=1,
            encoder_hidden_dim=16,
            encoder_n_epochs=10,
            search_n_iterations=20,
            search_n_ants=10,
        )
        result = pipe.run()
        assert isinstance(result, PipelineResult)
        assert len(result.ranked_genes) > 0

    def test_l1_produces_valid_scores(self, toy_data) -> None:
        """L1 encoder should produce scores in [0, 1]."""
        alt, mv, labels = toy_data
        pipe = Pipeline(
            alteration_matrix=alt,
            multi_view=mv,
            labels=labels,
            encoder_level=1,
            encoder_hidden_dim=16,
            encoder_n_epochs=20,
            search_n_iterations=10,
            search_n_ants=5,
        )
        result = pipe.run()
        for _, score in result.ranked_genes:
            assert 0 <= score <= 1

    def test_l0_still_works(self, toy_data) -> None:
        """L0 (XGBoost) should still work with encoder_level parameter."""
        alt, mv, labels = toy_data
        pipe = Pipeline(
            alteration_matrix=alt,
            multi_view=mv,
            labels=labels,
            encoder_level=0,
            search_n_iterations=10,
            search_n_ants=5,
        )
        result = pipe.run()
        assert isinstance(result, PipelineResult)

    def test_default_is_l0(self, toy_data) -> None:
        """Default encoder_level should be 0."""
        alt, mv, labels = toy_data
        pipe = Pipeline(
            alteration_matrix=alt,
            multi_view=mv,
            labels=labels,
            search_n_iterations=10,
            search_n_ants=5,
        )
        assert pipe.encoder_level == 0
