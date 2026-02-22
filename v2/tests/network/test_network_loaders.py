"""Tests for network loaders (WP2.5). Written FIRST per TDD."""

from pathlib import Path

import networkx as nx
import pytest

from ssame.network.biased_ppi import load_string
from ssame.network.experimental_ppi import load_hi_union
from ssame.network.pathway_loader import load_reactome


@pytest.fixture
def hi_union_tsv(tmp_path: Path) -> Path:
    """Create synthetic HI-II-14 format file."""
    content = (
        "Symbol_A\tSymbol_B\n"
        "TP53\tMDM2\n"
        "BRCA1\tBARD1\n"
        "PIK3CA\tPIK3R1\n"
        "KRAS\tBRAF\n"
        "TP53\tCDKN2A\n"
    )
    f = tmp_path / "HI-II-14.tsv"
    f.write_text(content)
    return f


@pytest.fixture
def reactome_gmt(tmp_path: Path) -> Path:
    """Create synthetic Reactome GMT format file."""
    content = (
        "R-HSA-109582\tHemostasis\tTP53\tBRCA1\tPIK3CA\n"
        "R-HSA-162582\tSignal Transduction\tKRAS\tBRAF\tMEK1\tERK1\n"
        "R-HSA-73857\tRNA Polymerase\tPOLR2A\tPOLR2B\n"
    )
    f = tmp_path / "ReactomePathways.gmt"
    f.write_text(content)
    return f


@pytest.fixture
def string_tsv(tmp_path: Path) -> Path:
    """Create synthetic STRING format file."""
    content = (
        "protein1\tprotein2\tcombined_score\n"
        "9606.ENSP00000269305\t9606.ENSP00000270202\t999\n"
        "9606.ENSP00000275493\t9606.ENSP00000418960\t800\n"
        "9606.ENSP00000256078\t9606.ENSP00000263967\t700\n"
        "9606.ENSP00000256078\t9606.ENSP00000318219\t400\n"
    )
    f = tmp_path / "9606.protein.links.v12.0.txt"
    f.write_text(content)
    return f


@pytest.fixture
def string_aliases(tmp_path: Path) -> Path:
    """Create STRING protein aliases for ID mapping."""
    content = (
        "#string_protein_id\talias\tsource\n"
        "9606.ENSP00000269305\tTP53\tEnsembl_HGNC\n"
        "9606.ENSP00000270202\tMDM2\tEnsembl_HGNC\n"
        "9606.ENSP00000275493\tBRCA1\tEnsembl_HGNC\n"
        "9606.ENSP00000418960\tBARD1\tEnsembl_HGNC\n"
        "9606.ENSP00000256078\tPIK3CA\tEnsembl_HGNC\n"
        "9606.ENSP00000263967\tPIK3R1\tEnsembl_HGNC\n"
        "9606.ENSP00000318219\tAKT1\tEnsembl_HGNC\n"
    )
    f = tmp_path / "9606.protein.aliases.v12.0.txt"
    f.write_text(content)
    return f


class TestHIUnionLoader:
    def test_returns_networkx_graph(self, hi_union_tsv: Path) -> None:
        G = load_hi_union(hi_union_tsv)
        assert isinstance(G, nx.Graph)

    def test_correct_node_count(self, hi_union_tsv: Path) -> None:
        G = load_hi_union(hi_union_tsv)
        # 6 unique genes: TP53, MDM2, BRCA1, BARD1, PIK3CA, PIK3R1, KRAS, BRAF, CDKN2A
        assert G.number_of_nodes() == 9

    def test_correct_edge_count(self, hi_union_tsv: Path) -> None:
        G = load_hi_union(hi_union_tsv)
        assert G.number_of_edges() == 5

    def test_edge_has_source_attribute(self, hi_union_tsv: Path) -> None:
        G = load_hi_union(hi_union_tsv)
        for _, _, data in G.edges(data=True):
            assert data["source"] == "HI-II-14"

    def test_known_edge_exists(self, hi_union_tsv: Path) -> None:
        G = load_hi_union(hi_union_tsv)
        assert G.has_edge("TP53", "MDM2")

    def test_no_self_loops(self, hi_union_tsv: Path) -> None:
        G = load_hi_union(hi_union_tsv)
        assert nx.number_of_selfloops(G) == 0


class TestReactomeLoader:
    def test_returns_networkx_graph(self, reactome_gmt: Path) -> None:
        G = load_reactome(reactome_gmt)
        assert isinstance(G, nx.Graph)

    def test_pathway_genes_connected(self, reactome_gmt: Path) -> None:
        G = load_reactome(reactome_gmt)
        # Genes in same pathway should be connected
        assert G.has_edge("KRAS", "BRAF")

    def test_cross_pathway_not_connected(self, reactome_gmt: Path) -> None:
        G = load_reactome(reactome_gmt)
        # POLR2A is only in RNA Polymerase pathway, not connected to KRAS
        assert not G.has_edge("POLR2A", "KRAS")

    def test_edge_has_source(self, reactome_gmt: Path) -> None:
        G = load_reactome(reactome_gmt)
        for _, _, data in G.edges(data=True):
            assert data["source"] == "Reactome"

    def test_edge_has_pathway(self, reactome_gmt: Path) -> None:
        G = load_reactome(reactome_gmt)
        # At least some edges should have pathway annotation
        edge_data = G.get_edge_data("KRAS", "BRAF")
        assert "pathway" in edge_data


class TestSTRINGLoader:
    def test_returns_networkx_graph(
        self, string_tsv: Path, string_aliases: Path
    ) -> None:
        G = load_string(string_tsv, string_aliases)
        assert isinstance(G, nx.Graph)

    def test_uses_hugo_symbols(
        self, string_tsv: Path, string_aliases: Path
    ) -> None:
        G = load_string(string_tsv, string_aliases)
        # Nodes should be Hugo symbols, not ENSP IDs
        assert "TP53" in G.nodes()
        assert "9606.ENSP00000269305" not in G.nodes()

    def test_edge_has_confidence(
        self, string_tsv: Path, string_aliases: Path
    ) -> None:
        G = load_string(string_tsv, string_aliases)
        for _, _, data in G.edges(data=True):
            assert "confidence" in data
            assert 0 <= data["confidence"] <= 1

    def test_confidence_filter(
        self, string_tsv: Path, string_aliases: Path
    ) -> None:
        G = load_string(string_tsv, string_aliases, min_score=700)
        # Edge with score 400 should be filtered out
        assert not G.has_edge("PIK3CA", "AKT1")
        # Edge with score 999 should remain
        assert G.has_edge("TP53", "MDM2")

    def test_edge_has_source(
        self, string_tsv: Path, string_aliases: Path
    ) -> None:
        G = load_string(string_tsv, string_aliases)
        for _, _, data in G.edges(data=True):
            assert data["source"] == "STRING"
