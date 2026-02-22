"""Tests for CGC ground truth loader (WP2.1). Written FIRST per TDD."""

from pathlib import Path

import pytest

from ssame.data.ground_truth import load_cgc


@pytest.fixture
def cgc_csv(tmp_path: Path) -> Path:
    """Create a synthetic CGC-format CSV for testing."""
    content = (
        "Gene Symbol,Name,Entrez GeneId,Genome Location,"
        "Tier,Hallmark,Chr Band,Somatic,Germline,"
        "Tumour Types(Somatic),Tumour Types(Germline),"
        "Cancer Syndrome,Tissue Type,Molecular Genetics,"
        "Role in Cancer,Mutation Types,Translocation Partner,"
        "Other Germline Mut,Other Syndrome,Synonyms,"
        "First Curated Date\n"
        'TP53,Tumor protein p53,7157,"17p13.1",'
        "1,Yes,17p13.1,Yes,Yes,"
        '"breast, lung, ovarian","Li-Fraumeni",'
        "Li-Fraumeni,S+E,Dom,"
        '"TSG, oncogene","Mis, N, F",,,,"TRP53",'
        "2004-01-01\n"
        'BRCA1,BRCA1 DNA repair,672,"17q21.31",'
        "1,Yes,17q21.31,Yes,Yes,"
        '"breast, ovarian","BRCA1",'
        "BRCA1,E,Dom,"
        "TSG,Mis,,,,"
        ",2005-06-15\n"
        'PIK3CA,PI3K catalytic alpha,5290,"3q26.32",'
        "1,Yes,3q26.32,Yes,No,"
        '"breast, colon, glioblastoma",,'
        ",S,Dom,"
        "oncogene,Mis,,,,"
        ",2010-03-01\n"
        'CDH1,Cadherin 1,999,"16q22.1",'
        "2,Yes,16q22.1,Yes,Yes,"
        '"breast, gastric","HDGC",'
        "HDGC,S+E,Dom,"
        "TSG,Mis,,,,"
        ",2018-07-20\n"
        'NEW_GENE,New discovery,99999,"1p36",'
        "1,Yes,1p36,Yes,No,"
        '"lung",,'
        ",S,Dom,"
        "oncogene,Mis,,,,"
        ",2022-01-15\n"
        'RECENT_GENE,Recent addition,99998,"2q33",'
        "1,Yes,2q33,Yes,No,"
        '"breast",,'
        ",S,Dom,"
        "TSG,Mis,,,,"
        ",2024-06-01\n"
    )
    csv_file = tmp_path / "cancer_gene_census.csv"
    csv_file.write_text(content)
    return csv_file


class TestLoadCGC:
    def test_load_returns_gene_set(self, cgc_csv: Path) -> None:
        result = load_cgc(cgc_csv)
        assert isinstance(result.positives, set)
        assert len(result.positives) > 0

    def test_known_genes_present(self, cgc_csv: Path) -> None:
        result = load_cgc(cgc_csv)
        assert "TP53" in result.positives
        assert "BRCA1" in result.positives
        assert "PIK3CA" in result.positives

    def test_tier1_filter(self, cgc_csv: Path) -> None:
        result = load_cgc(cgc_csv, tier=1)
        # CDH1 is Tier 2 in our test data
        assert "CDH1" not in result.positives
        assert "TP53" in result.positives

    def test_all_tiers(self, cgc_csv: Path) -> None:
        result = load_cgc(cgc_csv)
        # By default, all tiers included
        assert "CDH1" in result.positives

    def test_temporal_split_before(self, cgc_csv: Path) -> None:
        result = load_cgc(cgc_csv, before_year=2020)
        assert "TP53" in result.positives       # 2004
        assert "PIK3CA" in result.positives     # 2010
        assert "CDH1" in result.positives       # 2018
        assert "NEW_GENE" not in result.positives    # 2022
        assert "RECENT_GENE" not in result.positives  # 2024

    def test_temporal_holdout(self, cgc_csv: Path) -> None:
        result = load_cgc(cgc_csv, before_year=2020)
        assert "NEW_GENE" in result.temporal_holdout      # 2022
        assert "RECENT_GENE" in result.temporal_holdout   # 2024
        assert "TP53" not in result.temporal_holdout      # 2004

    def test_pu_framing_no_negatives(self, cgc_csv: Path) -> None:
        """PU framing: unlabeled genes are NOT labeled as negative."""
        result = load_cgc(cgc_csv)
        # There should be no negatives attribute
        assert not hasattr(result, "negatives")

    def test_role_filter(self, cgc_csv: Path) -> None:
        result = load_cgc(cgc_csv, role="TSG")
        assert "TP53" in result.positives
        assert "PIK3CA" not in result.positives  # oncogene only

    def test_result_has_metadata(self, cgc_csv: Path) -> None:
        result = load_cgc(cgc_csv)
        # Should be able to get gene metadata
        meta = result.get_metadata("TP53")
        assert meta is not None
        assert meta["entrez_id"] == 7157
