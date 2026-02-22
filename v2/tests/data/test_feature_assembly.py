"""Tests for feature assembly (WP3.1, 3.5, 3.6). Written FIRST per TDD."""

import networkx as nx
import numpy as np
import pytest

from ssame.data.feature_assembly import FeatureAssembler
from ssame.network.multi_view_graph import MultiViewGraph, build_multi_view
from ssame.types import AlterationMatrix


@pytest.fixture
def sample_matrix() -> AlterationMatrix:
    genes = ["TP53", "PIK3CA", "CDH1", "KRAS", "BRAF"]
    samples = [f"S{i}" for i in range(20)]
    data = np.zeros((5, 20), dtype=np.int8)
    # TP53: high frequency
    data[0, :15] = 1
    # PIK3CA: moderate frequency
    data[1, :8] = 1
    # CDH1: low frequency
    data[2, :3] = 1
    # KRAS: moderate
    data[3, :6] = 1
    # BRAF: very low
    data[4, :1] = 1
    return AlterationMatrix(data=data, genes=genes, samples=samples)


@pytest.fixture
def sample_multi_view() -> MultiViewGraph:
    g1 = nx.Graph()
    g1.add_edges_from([
        ("TP53", "MDM2"), ("TP53", "CDKN2A"), ("PIK3CA", "PIK3R1"),
        ("KRAS", "BRAF"), ("CDH1", "CTNNB1"),
    ])
    g2 = nx.Graph()
    g2.add_edges_from([
        ("TP53", "BRCA1"), ("PIK3CA", "AKT1"), ("KRAS", "MEK1"),
    ])
    return build_multi_view({"ppi": g1, "pathway": g2})


class TestFeatureAssembler:
    def test_creates_feature_matrix(
        self, sample_matrix: AlterationMatrix, sample_multi_view: MultiViewGraph
    ) -> None:
        assembler = FeatureAssembler(
            alteration_matrix=sample_matrix,
            multi_view=sample_multi_view,
        )
        X, gene_list = assembler.assemble()
        assert isinstance(X, np.ndarray)
        assert X.ndim == 2
        assert len(gene_list) > 0

    def test_feature_matrix_shape(
        self, sample_matrix: AlterationMatrix, sample_multi_view: MultiViewGraph
    ) -> None:
        assembler = FeatureAssembler(
            alteration_matrix=sample_matrix,
            multi_view=sample_multi_view,
        )
        X, gene_list = assembler.assemble()
        # Rows = genes in the intersection of alteration matrix and network
        assert X.shape[0] == len(gene_list)
        # Columns = number of features
        assert X.shape[1] > 0

    def test_no_nan_in_output(
        self, sample_matrix: AlterationMatrix, sample_multi_view: MultiViewGraph
    ) -> None:
        assembler = FeatureAssembler(
            alteration_matrix=sample_matrix,
            multi_view=sample_multi_view,
        )
        X, _ = assembler.assemble()
        assert not np.isnan(X).any()

    def test_known_drivers_have_higher_frequency_feature(
        self, sample_matrix: AlterationMatrix, sample_multi_view: MultiViewGraph
    ) -> None:
        assembler = FeatureAssembler(
            alteration_matrix=sample_matrix,
            multi_view=sample_multi_view,
        )
        X, gene_list = assembler.assemble()
        feature_names = assembler.feature_names
        if "mutation_frequency" in feature_names:
            freq_idx = feature_names.index("mutation_frequency")
            # TP53 should have highest frequency if present
            if "TP53" in gene_list and "BRAF" in gene_list:
                tp53_idx = gene_list.index("TP53")
                braf_idx = gene_list.index("BRAF")
                assert X[tp53_idx, freq_idx] > X[braf_idx, freq_idx]

    def test_feature_names_match_columns(
        self, sample_matrix: AlterationMatrix, sample_multi_view: MultiViewGraph
    ) -> None:
        assembler = FeatureAssembler(
            alteration_matrix=sample_matrix,
            multi_view=sample_multi_view,
        )
        X, _ = assembler.assemble()
        assert len(assembler.feature_names) == X.shape[1]

    def test_bimodal_split(
        self, sample_matrix: AlterationMatrix, sample_multi_view: MultiViewGraph
    ) -> None:
        """Feature assembly should separate omics vs topology features."""
        assembler = FeatureAssembler(
            alteration_matrix=sample_matrix,
            multi_view=sample_multi_view,
        )
        X, _ = assembler.assemble()
        omics_idx, topo_idx = assembler.get_bimodal_indices()
        # Both should be non-empty
        assert len(omics_idx) > 0
        assert len(topo_idx) > 0
        # Together they cover all features
        assert len(omics_idx) + len(topo_idx) == X.shape[1]
        # No overlap
        assert len(set(omics_idx) & set(topo_idx)) == 0

    def test_standardize(
        self, sample_matrix: AlterationMatrix, sample_multi_view: MultiViewGraph
    ) -> None:
        assembler = FeatureAssembler(
            alteration_matrix=sample_matrix,
            multi_view=sample_multi_view,
        )
        X, _ = assembler.assemble(standardize=True)
        # After standardization, mean should be ~0 and std ~1
        # (with tolerance for small datasets)
        col_means = X.mean(axis=0)
        assert np.allclose(col_means, 0, atol=0.1)
