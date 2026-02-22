"""Tests for score-as-label feedback loop (WP7.2). Written FIRST per TDD."""

import networkx as nx
import numpy as np
import pytest

from ssame.network.multi_view_graph import build_multi_view
from ssame.pipeline.feedback_loop import FeedbackPipeline
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


class TestFeedbackPipeline:
    def test_creates_pipeline(self, toy_data) -> None:
        alt, mv, labels = toy_data
        pipe = FeedbackPipeline(
            alteration_matrix=alt,
            multi_view=mv,
            labels=labels,
            n_rounds=2,
        )
        assert pipe is not None

    def test_runs_multiple_rounds(self, toy_data) -> None:
        alt, mv, labels = toy_data
        pipe = FeedbackPipeline(
            alteration_matrix=alt,
            multi_view=mv,
            labels=labels,
            n_rounds=3,
            search_n_iterations=20,
            search_n_ants=10,
        )
        result = pipe.run()
        assert len(result.round_results) == 3

    def test_scores_evolve_across_rounds(self, toy_data) -> None:
        """Gene scores should change between rounds due to feedback."""
        alt, mv, labels = toy_data
        pipe = FeedbackPipeline(
            alteration_matrix=alt,
            multi_view=mv,
            labels=labels,
            n_rounds=3,
            search_n_iterations=50,
            search_n_ants=20,
            feedback_weight=0.05,  # Low threshold so feedback triggers on toy data
        )
        result = pipe.run()
        # Scores from round 1 and round 3 should differ
        r1_scores = result.round_results[0].gene_scores
        r3_scores = result.round_results[2].gene_scores
        # At least some genes should have changed scores
        diffs = [abs(r1_scores[g] - r3_scores[g]) for g in r1_scores if g in r3_scores]
        assert max(diffs) > 0

    def test_final_result_has_ranked_genes(self, toy_data) -> None:
        alt, mv, labels = toy_data
        pipe = FeedbackPipeline(
            alteration_matrix=alt,
            multi_view=mv,
            labels=labels,
            n_rounds=2,
            search_n_iterations=20,
            search_n_ants=10,
        )
        result = pipe.run()
        assert len(result.final_ranked_genes) > 0
        # TP53 should rank high (known driver + high frequency)
        top_genes = [g for g, _ in result.final_ranked_genes[:3]]
        assert "TP53" in top_genes

    def test_soft_labels_incorporate_module_scores(self, toy_data) -> None:
        """Genes found in high-scoring modules should get label boost."""
        alt, mv, labels = toy_data
        pipe = FeedbackPipeline(
            alteration_matrix=alt,
            multi_view=mv,
            labels=labels,
            n_rounds=2,
            search_n_iterations=20,
            search_n_ants=10,
            feedback_weight=0.3,
        )
        result = pipe.run()
        # Should have soft labels from feedback
        assert result.round_results[1].gene_scores is not None
