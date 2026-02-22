"""Tests for AlphaFold PPI loader (WP2.6). Written FIRST per TDD."""

from pathlib import Path

import networkx as nx
import pytest

from ssame.network.alphafold_ppi import load_alphafold_ppi


@pytest.fixture
def rf2_ppi_tsv(tmp_path: Path) -> Path:
    """Create a fixture mimicking RoseTTAFold2-PPI final_predictions_90.tsv.

    Based on the Dryad dataset (doi:10.5061/dryad.15dv41p84).
    Format: protein_a, protein_b, max_prob, source
    (source: DNS=de novo shared localization, DNU=de novo unknown, PPI=known database)
    """
    content = (
        "protein_a\tprotein_b\tmax_prob\tsource\n"
        "TP53\tMDM2\t0.98\tPPI\n"
        "KRAS\tBRAF\t0.95\tDNS\n"
        "PIK3CA\tPTEN\t0.92\tDNU\n"
        "CDH1\tCTNNB1\t0.88\tDNS\n"  # Below 0.9 threshold
        "EGFR\tGRB2\t0.91\tPPI\n"
        "TP53\tTP53\t0.99\tPPI\n"  # Self-interaction, should be excluded
        "BRCA1\tBRCA2\t0.85\tDNU\n"  # Below 0.9 threshold
    )
    f = tmp_path / "final_predictions_90.tsv"
    f.write_text(content)
    return f


@pytest.fixture
def uniprot_format_tsv(tmp_path: Path) -> Path:
    """Alternative format with UniProt IDs and gene name columns."""
    content = (
        "uniprot_a\tuniprot_b\tgene_a\tgene_b\tscore\n"
        "P04637\tQ00987\tTP53\tMDM2\t0.98\n"
        "P01116\tP15056\tKRAS\tBRAF\t0.95\n"
        "P42336\tP60484\tPIK3CA\tPTEN\t0.80\n"
    )
    f = tmp_path / "predictions_uniprot.tsv"
    f.write_text(content)
    return f


class TestLoadAlphaFoldPPI:
    def test_loads_tsv_into_graph(self, rf2_ppi_tsv: Path) -> None:
        G = load_alphafold_ppi(rf2_ppi_tsv, min_confidence=0.0)
        assert isinstance(G, nx.Graph)
        assert G.number_of_nodes() > 0

    def test_default_confidence_filter_90(self, rf2_ppi_tsv: Path) -> None:
        """Default min_confidence=0.9 should filter out low-confidence pairs."""
        G = load_alphafold_ppi(rf2_ppi_tsv)
        # TP53-MDM2 (0.98), KRAS-BRAF (0.95), PIK3CA-PTEN (0.92), EGFR-GRB2 (0.91) pass
        # CDH1-CTNNB1 (0.88), BRCA1-BRCA2 (0.85) filtered out
        assert G.number_of_edges() == 4
        assert G.has_edge("TP53", "MDM2")
        assert G.has_edge("KRAS", "BRAF")
        assert not G.has_edge("CDH1", "CTNNB1")
        assert not G.has_edge("BRCA1", "BRCA2")

    def test_excludes_self_interactions(self, rf2_ppi_tsv: Path) -> None:
        G = load_alphafold_ppi(rf2_ppi_tsv, min_confidence=0.0)
        assert not G.has_edge("TP53", "TP53")

    def test_edge_attributes(self, rf2_ppi_tsv: Path) -> None:
        G = load_alphafold_ppi(rf2_ppi_tsv)
        edge_data = G.edges["TP53", "MDM2"]
        assert edge_data["source"] == "AlphaFold_Baker2025"
        assert edge_data["confidence"] == pytest.approx(0.98)

    def test_missing_file_returns_empty_graph(self, tmp_path: Path) -> None:
        with pytest.warns(UserWarning, match="AlphaFold PPI file not found"):
            G = load_alphafold_ppi(tmp_path / "nonexistent.tsv")
        assert G.number_of_nodes() == 0

    def test_custom_confidence_threshold(self, rf2_ppi_tsv: Path) -> None:
        G = load_alphafold_ppi(rf2_ppi_tsv, min_confidence=0.95)
        # Only TP53-MDM2 (0.98) and KRAS-BRAF (0.95) pass
        assert G.number_of_edges() == 2

    def test_undirected_graph(self, rf2_ppi_tsv: Path) -> None:
        G = load_alphafold_ppi(rf2_ppi_tsv)
        assert not G.is_directed()
        # Both directions work
        assert G.has_edge("MDM2", "TP53")

    def test_preserves_source_column(self, rf2_ppi_tsv: Path) -> None:
        """The 'source' column from the TSV (DNS/DNU/PPI) should be preserved."""
        G = load_alphafold_ppi(rf2_ppi_tsv)
        edge_data = G.edges["KRAS", "BRAF"]
        assert edge_data.get("prediction_source") == "DNS"

    def test_no_confidence_column(self, tmp_path: Path) -> None:
        """Files without a confidence column should load all edges."""
        content = "gene_a\tgene_b\n" "TP53\tMDM2\n" "KRAS\tBRAF\n"
        f = tmp_path / "simple.tsv"
        f.write_text(content)
        G = load_alphafold_ppi(f)
        assert G.number_of_edges() == 2

    def test_gene_name_columns_auto_detected(self, uniprot_format_tsv: Path) -> None:
        """Should handle files where gene names are in gene_a/gene_b columns."""
        G = load_alphafold_ppi(
            uniprot_format_tsv,
            min_confidence=0.0,
            gene_a_col="gene_a",
            gene_b_col="gene_b",
            confidence_col="score",
        )
        assert G.has_edge("TP53", "MDM2")
        assert G.number_of_edges() == 3
