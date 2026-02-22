"""Tests for Pan-Cancer framework (WP10). Written FIRST per TDD.

Runs the pipeline independently per cancer type, then aggregates
results using rank aggregation (Borda count or similar).
"""

import numpy as np
import pytest

from ssame.pipeline.pan_cancer import PanCancerResult, PanCancerRunner
from ssame.types import AlterationMatrix


def _make_alteration_matrix(n_genes: int = 50, n_samples: int = 30, seed: int = 42):
    rng = np.random.default_rng(seed)
    genes = [f"GENE{i}" for i in range(n_genes)]
    samples = [f"S{i}" for i in range(n_samples)]
    data = (rng.random((n_genes, n_samples)) < 0.15).astype(np.int8)
    return AlterationMatrix(data=data, genes=genes, samples=samples)


@pytest.fixture
def cancer_types():
    return {
        "BRCA": _make_alteration_matrix(50, 30, seed=1),
        "LUAD": _make_alteration_matrix(50, 25, seed=2),
        "COAD": _make_alteration_matrix(50, 20, seed=3),
    }


@pytest.fixture
def shared_labels():
    """Shared known driver labels across cancer types."""
    return {"GENE0": 1, "GENE1": 1, "GENE2": 1}


class TestPanCancerRunner:
    def test_creates_runner(self, cancer_types, shared_labels) -> None:
        runner = PanCancerRunner(
            alteration_matrices=cancer_types,
            labels=shared_labels,
            encoder_level=0,
            encoder_n_estimators=10,
            search_n_iterations=10,
            seed=42,
        )
        assert runner is not None

    def test_run_returns_result(self, cancer_types, shared_labels) -> None:
        runner = PanCancerRunner(
            alteration_matrices=cancer_types,
            labels=shared_labels,
            encoder_level=0,
            encoder_n_estimators=10,
            search_n_iterations=10,
            seed=42,
        )
        result = runner.run()
        assert isinstance(result, PanCancerResult)

    def test_result_has_per_cancer_results(self, cancer_types, shared_labels) -> None:
        runner = PanCancerRunner(
            alteration_matrices=cancer_types,
            labels=shared_labels,
            encoder_level=0,
            encoder_n_estimators=10,
            search_n_iterations=10,
            seed=42,
        )
        result = runner.run()
        assert len(result.per_cancer) == 3
        assert "BRCA" in result.per_cancer
        assert "LUAD" in result.per_cancer

    def test_result_has_aggregated_ranking(self, cancer_types, shared_labels) -> None:
        runner = PanCancerRunner(
            alteration_matrices=cancer_types,
            labels=shared_labels,
            encoder_level=0,
            encoder_n_estimators=10,
            search_n_iterations=10,
            seed=42,
        )
        result = runner.run()
        assert len(result.aggregated_ranking) > 0
        # Should be sorted: (gene, score) tuples, descending score
        scores = [s for _, s in result.aggregated_ranking]
        assert scores == sorted(scores, reverse=True)

    def test_aggregated_ranking_covers_all_genes(self, cancer_types, shared_labels) -> None:
        runner = PanCancerRunner(
            alteration_matrices=cancer_types,
            labels=shared_labels,
            encoder_level=0,
            encoder_n_estimators=10,
            search_n_iterations=10,
            seed=42,
        )
        result = runner.run()
        ranked_genes = {g for g, _ in result.aggregated_ranking}
        # All genes from any cancer type should appear
        all_genes = set()
        for am in cancer_types.values():
            all_genes.update(am.genes)
        assert ranked_genes == all_genes

    def test_recurrence_count(self, cancer_types, shared_labels) -> None:
        runner = PanCancerRunner(
            alteration_matrices=cancer_types,
            labels=shared_labels,
            encoder_level=0,
            encoder_n_estimators=10,
            search_n_iterations=10,
            seed=42,
        )
        result = runner.run()
        # recurrence_count: gene → how many cancer types it was flagged in
        assert isinstance(result.recurrence_count, dict)
        for _gene, count in result.recurrence_count.items():
            assert 0 <= count <= 3
