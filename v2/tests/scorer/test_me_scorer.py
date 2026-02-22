"""Tests for MEScorer (WP6.2). Written FIRST per TDD."""

import numpy as np
import pytest

from ssame.scorer.me_scorer import MEScorer
from ssame.types import AlterationMatrix, Gene, Subnetwork


@pytest.fixture
def me_matrix() -> AlterationMatrix:
    """Alteration matrix with clear ME pattern for TP53/PIK3CA."""
    genes = ["TP53", "PIK3CA", "CDH1", "RANDOM"]
    samples = [f"S{i}" for i in range(20)]
    data = np.zeros((4, 20), dtype=np.int8)
    # ME pattern: TP53 and PIK3CA rarely co-occur
    data[0, :10] = 1   # TP53 in samples 0-9
    data[1, 10:18] = 1  # PIK3CA in samples 10-17
    # CDH1 overlaps with both (not ME)
    data[2, 5:15] = 1
    # RANDOM sparse
    data[3, [0, 5, 15]] = 1
    return AlterationMatrix(data=data, genes=genes, samples=samples)


@pytest.fixture
def me_module() -> Subnetwork:
    return Subnetwork(
        genes=[Gene(symbol="TP53"), Gene(symbol="PIK3CA")],
        interactions=[],
    )


@pytest.fixture
def non_me_module() -> Subnetwork:
    return Subnetwork(
        genes=[Gene(symbol="TP53"), Gene(symbol="CDH1")],
        interactions=[],
    )


class TestMEScorer:
    def test_creates_scorer(self, me_matrix: AlterationMatrix) -> None:
        scorer = MEScorer(alteration_matrix=me_matrix)
        assert scorer is not None

    def test_me_module_scores_higher(
        self,
        me_matrix: AlterationMatrix,
        me_module: Subnetwork,
        non_me_module: Subnetwork,
    ) -> None:
        scorer = MEScorer(alteration_matrix=me_matrix)
        me_score = scorer.score(me_module)
        non_me_score = scorer.score(non_me_module)
        assert me_score > non_me_score

    def test_score_is_nonnegative(
        self, me_matrix: AlterationMatrix, me_module: Subnetwork
    ) -> None:
        scorer = MEScorer(alteration_matrix=me_matrix)
        assert scorer.score(me_module) >= 0

    def test_explain_returns_explanation(
        self, me_matrix: AlterationMatrix, me_module: Subnetwork
    ) -> None:
        scorer = MEScorer(alteration_matrix=me_matrix)
        explanation = scorer.explain(me_module)
        assert explanation.score > 0
        assert "coverage" in explanation.evidence

    def test_coverage_in_explanation(
        self, me_matrix: AlterationMatrix, me_module: Subnetwork
    ) -> None:
        scorer = MEScorer(alteration_matrix=me_matrix)
        explanation = scorer.explain(me_module)
        coverage = explanation.evidence["coverage"]
        # TP53 covers 0-9, PIK3CA covers 10-17, union = 18/20
        assert coverage == pytest.approx(18 / 20)

    def test_perfect_me_module(self) -> None:
        """Perfectly ME module: no co-occurrence at all."""
        genes = ["A", "B"]
        samples = [f"S{i}" for i in range(10)]
        data = np.zeros((2, 10), dtype=np.int8)
        data[0, :5] = 1   # A in first 5
        data[1, 5:] = 1   # B in last 5
        matrix = AlterationMatrix(data=data, genes=genes, samples=samples)
        scorer = MEScorer(alteration_matrix=matrix)
        module = Subnetwork(
            genes=[Gene(symbol="A"), Gene(symbol="B")],
            interactions=[],
        )
        score = scorer.score(module)
        assert score > 0

    def test_v1_dendrix_style_scoring(
        self, me_matrix: AlterationMatrix
    ) -> None:
        """Dendrix-style: score = coverage - overlap (normalized)."""
        scorer = MEScorer(alteration_matrix=me_matrix, method="dendrix")
        module = Subnetwork(
            genes=[Gene(symbol="TP53"), Gene(symbol="PIK3CA")],
            interactions=[],
        )
        score = scorer.score(module)
        # coverage(TP53 | PIK3CA) = 18, overlap(TP53 & PIK3CA) = 0
        # Dendrix-like: (18 - 0) / 20 = 0.9 (normalized)
        assert score > 0.5
