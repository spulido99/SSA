"""Tests for StructuralDisruptionScorer (WP6.4). Written FIRST per TDD."""

import pytest

from ssame.scorer.protocol import SubnetworkScorer
from ssame.scorer.structural_scorer import StructuralDisruptionScorer
from ssame.types import Gene, Interaction, ScorerExplanation, Subnetwork


@pytest.fixture
def interface_data() -> dict[tuple[str, str], list[int]]:
    """Per-edge interface residue positions."""
    return {
        ("TP53", "MDM2"): [17, 18, 19, 20, 21, 22, 23],
        ("KRAS", "BRAF"): [12, 13, 14, 60, 61],
        ("PIK3CA", "PTEN"): [345, 346, 347],
    }


@pytest.fixture
def mutation_positions() -> dict[str, list[int]]:
    """Per-gene mutation positions from TCGA."""
    return {
        "TP53": [18, 19, 175, 248, 273],  # 2 of 5 at interface
        "MDM2": [50, 51, 100],  # 2 of 3 at interface with TP53
        "KRAS": [12, 13, 61],  # 3 of 3 at interface (hotspot!)
        "BRAF": [600],  # V600E, NOT at interface
        "PIK3CA": [100, 200, 345],  # 1 of 3 at interface
        "PTEN": [130, 136],  # 0 at interface
    }


@pytest.fixture
def scorer(interface_data, mutation_positions) -> StructuralDisruptionScorer:
    return StructuralDisruptionScorer(
        interface_residues=interface_data,
        mutation_positions=mutation_positions,
    )


@pytest.fixture
def kras_braf_module() -> Subnetwork:
    return Subnetwork(
        genes=[Gene("KRAS"), Gene("BRAF")],
        interactions=[Interaction(Gene("KRAS"), Gene("BRAF"), source="test")],
    )


@pytest.fixture
def tp53_mdm2_module() -> Subnetwork:
    return Subnetwork(
        genes=[Gene("TP53"), Gene("MDM2")],
        interactions=[Interaction(Gene("TP53"), Gene("MDM2"), source="test")],
    )


class TestStructuralDisruptionScorer:
    def test_implements_protocol(self, scorer) -> None:
        assert isinstance(scorer, SubnetworkScorer)

    def test_scores_module(self, scorer, kras_braf_module) -> None:
        score = scorer.score(kras_braf_module)
        assert isinstance(score, float)
        assert score >= 0.0

    def test_high_interface_mutation_high_score(self, scorer, kras_braf_module) -> None:
        """KRAS has 3/3 mutations at KRAS-BRAF interface → high disruption score."""
        score = scorer.score(kras_braf_module)
        assert score > 0.5

    def test_no_interface_data_zero_score(self, scorer) -> None:
        """Module with genes not in interface data should score ~0."""
        module = Subnetwork(
            genes=[Gene("NOVEL1"), Gene("NOVEL2")],
            interactions=[Interaction(Gene("NOVEL1"), Gene("NOVEL2"), source="test")],
        )
        score = scorer.score(module)
        assert score == pytest.approx(0.0)

    def test_explain_returns_explanation(self, scorer, tp53_mdm2_module) -> None:
        explanation = scorer.explain(tp53_mdm2_module)
        assert isinstance(explanation, ScorerExplanation)
        assert explanation.narrative != ""

    def test_explain_has_per_gene_data(self, scorer, tp53_mdm2_module) -> None:
        explanation = scorer.explain(tp53_mdm2_module)
        assert "TP53" in explanation.components or "TP53" in str(explanation.evidence)

    def test_empty_module(self, scorer) -> None:
        module = Subnetwork(genes=[], interactions=[])
        score = scorer.score(module)
        assert score == 0.0
