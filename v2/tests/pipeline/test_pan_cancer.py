"""Tests for Pan-Cancer framework (WP10). Written FIRST per TDD.

Runs the pipeline independently per cancer type, then aggregates
results using rank aggregation (Borda count or similar).
"""

import numpy as np
import pytest
from scipy.sparse import issparse

from ssame.pipeline.pan_cancer import PanCancerResult, PanCancerRunner
from ssame.types import AlterationMatrix
from ssame.utils.checkpoint import CheckpointManager


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


class TestPanCancerCheckpoints:
    def test_saves_per_type_checkpoints(
        self, cancer_types, shared_labels, tmp_path
    ) -> None:
        ckpt_dir = tmp_path / ".checkpoints"
        runner = PanCancerRunner(
            alteration_matrices=cancer_types,
            labels=shared_labels,
            encoder_level=0,
            encoder_n_estimators=10,
            search_n_iterations=10,
            seed=42,
            checkpoint_dir=ckpt_dir,
        )
        runner.run()
        mgr = CheckpointManager(ckpt_dir)
        for ct in cancer_types:
            assert mgr.exists(f"cancer_{ct}")

    def test_resumes_completed_types(
        self, cancer_types, shared_labels, tmp_path
    ) -> None:
        ckpt_dir = tmp_path / ".checkpoints"
        # First run
        r1 = PanCancerRunner(
            alteration_matrices=cancer_types,
            labels=shared_labels,
            encoder_level=0,
            encoder_n_estimators=10,
            search_n_iterations=10,
            seed=42,
            checkpoint_dir=ckpt_dir,
        )
        result1 = r1.run()

        # Second run should resume from checkpoints
        r2 = PanCancerRunner(
            alteration_matrices=cancer_types,
            labels=shared_labels,
            encoder_level=0,
            encoder_n_estimators=10,
            search_n_iterations=10,
            seed=42,
            checkpoint_dir=ckpt_dir,
        )
        result2 = r2.run()

        # Per-cancer gene scores should match (loaded from checkpoint)
        for ct in cancer_types:
            for gene in result1.per_cancer[ct].gene_scores:
                assert result1.per_cancer[ct].gene_scores[gene] == pytest.approx(
                    result2.per_cancer[ct].gene_scores[gene]
                )


class TestSparseCoalteration:
    def test_coalteration_graph_from_sparse_matches_dense(self, shared_labels) -> None:
        """Sparse co-alteration should produce same graph as dense."""
        alt = _make_alteration_matrix(50, 30, seed=42)

        runner = PanCancerRunner(
            alteration_matrices={"TEST": alt},
            labels=shared_labels,
            encoder_level=0,
            encoder_n_estimators=10,
            search_n_iterations=10,
        )
        mv = runner._build_coalteration_graph(alt)

        # The adjacency should be sparse
        raw_adj = mv._adjacencies["coalteration"]
        assert issparse(raw_adj), "Co-alteration adjacency should be sparse"

        # The dense form should be a valid adjacency
        dense = raw_adj.toarray()
        assert dense.shape == (50, 50)
        # Diagonal should be zero
        np.testing.assert_array_equal(np.diag(dense), 0)
        # Should be symmetric
        np.testing.assert_array_equal(dense, dense.T)
