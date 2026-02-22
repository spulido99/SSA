"""Tests for PIONEER oncoPPI loader (WP2.7). Written FIRST per TDD."""

from pathlib import Path

import networkx as nx
import pytest

from ssame.network.pioneer_oncoppi import load_pioneer


@pytest.fixture
def pioneer_tsv(tmp_path: Path) -> Path:
    """Create a fixture mimicking PIONEER very_high confidence predictions.

    PIONEER predicts partner-specific protein interfaces.
    Format likely: protein_a, protein_b, interface_residues_a, interface_residues_b,
    confidence, etc.
    """
    content = (
        "protein_a\tprotein_b\tinterface_residues_a\tinterface_residues_b\tconfidence\n"
        "TP53\tMDM2\t17,18,19,20,21,22,23\t50,51,52,53,54\t0.95\n"
        "KRAS\tBRAF\t12,13,14,60,61\t100,101,102\t0.88\n"
        "PIK3CA\tPTEN\t345,346,347\t123,124,125,126\t0.92\n"
        "EGFR\tGRB2\t1068,1069\t49,50,51\t0.75\n"
        "BRCA1\tBRCA2\t1650,1651,1652\t3280,3281\t0.60\n"
    )
    f = tmp_path / "pioneer_very_high_human.txt"
    f.write_text(content)
    return f


@pytest.fixture
def pioneer_uniprot_tsv(tmp_path: Path) -> Path:
    """Alternative PIONEER format with UniProt IDs."""
    content = (
        "uniprot_a\tuniprot_b\tgene_a\tgene_b\tpioneer_score\tinterface_a\tinterface_b\n"
        "P04637\tQ00987\tTP53\tMDM2\t0.95\t17;18;19;20\t50;51;52\n"
        "P01116\tP15056\tKRAS\tBRAF\t0.88\t12;13;14\t100;101\n"
    )
    f = tmp_path / "pioneer_uniprot.tsv"
    f.write_text(content)
    return f


class TestLoadPioneer:
    def test_loads_into_graph(self, pioneer_tsv: Path) -> None:
        G = load_pioneer(pioneer_tsv)
        assert isinstance(G, nx.Graph)
        assert G.number_of_nodes() > 0

    def test_all_edges_loaded(self, pioneer_tsv: Path) -> None:
        G = load_pioneer(pioneer_tsv)
        assert G.number_of_edges() == 5

    def test_confidence_filter(self, pioneer_tsv: Path) -> None:
        """Should filter by confidence when threshold given."""
        G = load_pioneer(pioneer_tsv, min_confidence=0.85)
        assert G.has_edge("TP53", "MDM2")  # 0.95
        assert G.has_edge("KRAS", "BRAF")  # 0.88
        assert G.has_edge("PIK3CA", "PTEN")  # 0.92
        assert not G.has_edge("EGFR", "GRB2")  # 0.75
        assert not G.has_edge("BRCA1", "BRCA2")  # 0.60

    def test_edge_has_source_attribute(self, pioneer_tsv: Path) -> None:
        G = load_pioneer(pioneer_tsv)
        edge_data = G.edges["TP53", "MDM2"]
        assert edge_data["source"] == "PIONEER"

    def test_edge_has_confidence(self, pioneer_tsv: Path) -> None:
        G = load_pioneer(pioneer_tsv)
        assert G.edges["TP53", "MDM2"]["confidence"] == pytest.approx(0.95)

    def test_edge_has_interface_residues(self, pioneer_tsv: Path) -> None:
        """Interface residue data should be preserved on edges."""
        G = load_pioneer(pioneer_tsv)
        edge_data = G.edges["TP53", "MDM2"]
        assert "interface_residues_a" in edge_data
        assert "interface_residues_b" in edge_data

    def test_missing_file_returns_empty_graph(self, tmp_path: Path) -> None:
        with pytest.warns(UserWarning, match="PIONEER file not found"):
            G = load_pioneer(tmp_path / "nonexistent.tsv")
        assert G.number_of_nodes() == 0

    def test_undirected(self, pioneer_tsv: Path) -> None:
        G = load_pioneer(pioneer_tsv)
        assert not G.is_directed()
        assert G.has_edge("MDM2", "TP53")

    def test_get_interface_genes(self, pioneer_tsv: Path) -> None:
        """Should be able to extract genes with interface residue annotations."""
        G = load_pioneer(pioneer_tsv)
        interface_genes = {n for n in G.nodes() if any(
            "interface_residues" in str(k)
            for edge_data in G[n].values()
            for k in edge_data
        )}
        assert "TP53" in interface_genes
        assert "MDM2" in interface_genes

    def test_custom_columns(self, pioneer_uniprot_tsv: Path) -> None:
        """Should handle UniProt-format files with explicit column names."""
        G = load_pioneer(
            pioneer_uniprot_tsv,
            gene_a_col="gene_a",
            gene_b_col="gene_b",
            confidence_col="pioneer_score",
        )
        assert G.has_edge("TP53", "MDM2")
        assert G.number_of_edges() == 2
