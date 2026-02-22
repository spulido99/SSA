"""Tests for core domain types (WP0.5). Written FIRST per TDD."""

import numpy as np
import pytest

from ssame.types import (
    AlterationMatrix,
    Gene,
    Interaction,
    ScorerExplanation,
    Subnetwork,
)


class TestGene:
    def test_create_gene(self) -> None:
        gene = Gene(symbol="TP53")
        assert gene.symbol == "TP53"

    def test_gene_equality(self) -> None:
        assert Gene(symbol="TP53") == Gene(symbol="TP53")
        assert Gene(symbol="TP53") != Gene(symbol="BRCA1")

    def test_gene_hashable(self) -> None:
        genes = {Gene(symbol="TP53"), Gene(symbol="BRCA1"), Gene(symbol="TP53")}
        assert len(genes) == 2

    def test_gene_str(self) -> None:
        assert str(Gene(symbol="TP53")) == "TP53"

    def test_gene_with_entrez_id(self) -> None:
        gene = Gene(symbol="TP53", entrez_id=7157)
        assert gene.entrez_id == 7157

    def test_gene_rejects_empty_symbol(self) -> None:
        with pytest.raises(ValueError, match="symbol"):
            Gene(symbol="")


class TestInteraction:
    def test_create_interaction(self) -> None:
        interaction = Interaction(
            gene_a=Gene(symbol="TP53"),
            gene_b=Gene(symbol="MDM2"),
            source="HI-II-14",
        )
        assert interaction.gene_a.symbol == "TP53"
        assert interaction.gene_b.symbol == "MDM2"
        assert interaction.source == "HI-II-14"

    def test_interaction_confidence_default(self) -> None:
        interaction = Interaction(
            gene_a=Gene(symbol="TP53"),
            gene_b=Gene(symbol="MDM2"),
            source="STRING",
        )
        assert interaction.confidence == 1.0

    def test_interaction_with_confidence(self) -> None:
        interaction = Interaction(
            gene_a=Gene(symbol="TP53"),
            gene_b=Gene(symbol="MDM2"),
            source="STRING",
            confidence=0.9,
        )
        assert interaction.confidence == 0.9

    def test_interaction_as_tuple(self) -> None:
        interaction = Interaction(
            gene_a=Gene(symbol="TP53"),
            gene_b=Gene(symbol="MDM2"),
            source="HI-II-14",
        )
        assert interaction.as_tuple() == ("TP53", "MDM2")

    def test_interaction_rejects_self_loop(self) -> None:
        with pytest.raises(ValueError, match="self-loop"):
            Interaction(
                gene_a=Gene(symbol="TP53"),
                gene_b=Gene(symbol="TP53"),
                source="test",
            )


class TestSubnetwork:
    def test_create_subnetwork(self) -> None:
        genes = [Gene(symbol="TP53"), Gene(symbol="MDM2"), Gene(symbol="CDKN2A")]
        interactions = [
            Interaction(gene_a=genes[0], gene_b=genes[1], source="test"),
        ]
        subnet = Subnetwork(genes=genes, interactions=interactions)
        assert len(subnet.genes) == 3
        assert len(subnet.interactions) == 1

    def test_subnetwork_gene_symbols(self) -> None:
        genes = [Gene(symbol="TP53"), Gene(symbol="MDM2")]
        subnet = Subnetwork(genes=genes, interactions=[])
        assert subnet.gene_symbols == {"TP53", "MDM2"}

    def test_subnetwork_size(self) -> None:
        genes = [Gene(symbol="TP53"), Gene(symbol="MDM2"), Gene(symbol="RB1")]
        subnet = Subnetwork(genes=genes, interactions=[])
        assert subnet.size == 3

    def test_subnetwork_contains(self) -> None:
        genes = [Gene(symbol="TP53"), Gene(symbol="MDM2")]
        subnet = Subnetwork(genes=genes, interactions=[])
        assert subnet.contains("TP53")
        assert not subnet.contains("BRCA1")

    def test_subnetwork_merge(self) -> None:
        g1 = [Gene(symbol="TP53"), Gene(symbol="MDM2")]
        g2 = [Gene(symbol="MDM2"), Gene(symbol="RB1")]
        s1 = Subnetwork(genes=g1, interactions=[])
        s2 = Subnetwork(genes=g2, interactions=[])
        merged = s1.merge(s2)
        assert merged.gene_symbols == {"TP53", "MDM2", "RB1"}

    def test_subnetwork_intersection(self) -> None:
        g1 = [Gene(symbol="TP53"), Gene(symbol="MDM2")]
        g2 = [Gene(symbol="MDM2"), Gene(symbol="RB1")]
        s1 = Subnetwork(genes=g1, interactions=[])
        s2 = Subnetwork(genes=g2, interactions=[])
        common = s1.intersection(s2)
        assert common.gene_symbols == {"MDM2"}


class TestAlterationMatrix:
    def test_create_from_dense(self) -> None:
        genes = ["TP53", "PIK3CA", "CDH1"]
        samples = ["S1", "S2", "S3"]
        data = np.array([
            [1, 0, 1],
            [0, 1, 0],
            [1, 1, 0],
        ], dtype=np.int8)
        matrix = AlterationMatrix(data=data, genes=genes, samples=samples)
        assert matrix.n_genes == 3
        assert matrix.n_samples == 3

    def test_gene_frequency(self) -> None:
        genes = ["TP53", "PIK3CA"]
        samples = ["S1", "S2", "S3", "S4"]
        data = np.array([
            [1, 1, 0, 1],
            [0, 1, 0, 0],
        ], dtype=np.int8)
        matrix = AlterationMatrix(data=data, genes=genes, samples=samples)
        assert matrix.gene_frequency("TP53") == 0.75
        assert matrix.gene_frequency("PIK3CA") == 0.25

    def test_sample_coverage(self) -> None:
        genes = ["TP53", "PIK3CA"]
        samples = ["S1", "S2", "S3"]
        data = np.array([
            [1, 0, 0],
            [0, 1, 0],
        ], dtype=np.int8)
        matrix = AlterationMatrix(data=data, genes=genes, samples=samples)
        coverage = matrix.coverage(["TP53", "PIK3CA"])
        assert coverage == pytest.approx(2 / 3)

    def test_co_occurrence_count(self) -> None:
        genes = ["TP53", "PIK3CA"]
        samples = ["S1", "S2", "S3"]
        data = np.array([
            [1, 1, 0],
            [1, 0, 1],
        ], dtype=np.int8)
        matrix = AlterationMatrix(data=data, genes=genes, samples=samples)
        assert matrix.co_occurrence_count("TP53", "PIK3CA") == 1

    def test_filter_by_min_frequency(self) -> None:
        genes = ["TP53", "RARE_GENE"]
        samples = ["S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10"]
        data = np.zeros((2, 10), dtype=np.int8)
        data[0, :8] = 1  # TP53: 80% frequency
        data[1, 0] = 1   # RARE_GENE: 10% frequency
        matrix = AlterationMatrix(data=data, genes=genes, samples=samples)
        filtered = matrix.filter_by_min_frequency(0.2)
        assert filtered.n_genes == 1
        assert "TP53" in filtered.genes

    def test_shape_validation(self) -> None:
        with pytest.raises(ValueError, match="shape"):
            AlterationMatrix(
                data=np.zeros((2, 3), dtype=np.int8),
                genes=["A", "B", "C"],  # 3 genes but data has 2 rows
                samples=["S1", "S2", "S3"],
            )

    def test_submatrix(self) -> None:
        genes = ["TP53", "PIK3CA", "CDH1"]
        samples = ["S1", "S2"]
        data = np.array([[1, 0], [0, 1], [1, 1]], dtype=np.int8)
        matrix = AlterationMatrix(data=data, genes=genes, samples=samples)
        sub = matrix.submatrix(gene_subset=["TP53", "CDH1"])
        assert sub.n_genes == 2
        assert sub.genes == ["TP53", "CDH1"]


class TestScorerExplanation:
    def test_create_explanation(self) -> None:
        explanation = ScorerExplanation(
            score=0.85,
            components={"me_score": 0.7, "coverage": 0.95},
            narrative="Module of 4 genes with strong ME pattern",
            evidence={"p_value": 0.001},
        )
        assert explanation.score == 0.85
        assert explanation.components["me_score"] == 0.7
        assert "ME pattern" in explanation.narrative

    def test_explanation_defaults(self) -> None:
        explanation = ScorerExplanation(score=0.5)
        assert explanation.components == {}
        assert explanation.narrative == ""
        assert explanation.evidence == {}
