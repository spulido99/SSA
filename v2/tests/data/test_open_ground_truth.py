"""Tests for open-access ground truth loaders (IntOGen, NCG, CancerMine) + consensus.

Replaces CGC/OncoKB dependency with freely available databases.
Written FIRST per TDD.
"""

from pathlib import Path

import pytest

from ssame.data.ground_truth import GroundTruthResult
from ssame.data.open_ground_truth import (
    ConsensusGroundTruth,
    build_consensus_ground_truth,
    load_cancermine,
    load_intogen,
    load_ncg,
)

# ---------------------------------------------------------------------------
# IntOGen fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def intogen_tsv(tmp_path: Path) -> Path:
    """Mimics IntOGen drivers.tsv from IntOGen-Drivers-20240920.zip."""
    content = (
        "SYMBOL\tTUMOR_TYPE\tCOHORT\tSAMPLES\tQVALUE_COMBINATION\tROLE\t"
        "CGC_GENE\tMUTATIONS\tSIG_METHODS\tALL_METHODS\tRANKING\n"
        "TP53\tBRCA\tTCGA_BRCA\t1000\t1.2e-50\tLoF\tTrue\t350\tdndscv;hotmaps\tdndscv;hotmaps;cbase\t1\n"
        "KRAS\tLUAD\tTCGA_LUAD\t500\t3.5e-40\tAct\tTrue\t200\tdndscv;hotmaps\tdndscv;hotmaps;cbase\t1\n"
        "PIK3CA\tBRCA\tTCGA_BRCA\t1000\t1.1e-30\tAct\tTrue\t280\tdndscv\tdndscv;cbase\t2\n"
        "BRAF\tSKCM\tTCGA_SKCM\t300\t5.0e-25\tAct\tTrue\t180\thotmaps\tdndscv;hotmaps\t1\n"
        "PTEN\tUCEC\tTCGA_UCEC\t400\t2.0e-20\tLoF\tTrue\t150\tdndscv\tdndscv;cbase\t2\n"
        "NOVELGENE1\tBRCA\tTCGA_BRCA\t1000\t0.01\tLoF\tFalse\t5\tdndscv\tdndscv\t50\n"
    )
    f = tmp_path / "drivers.tsv"
    f.write_text(content)
    return f


# ---------------------------------------------------------------------------
# NCG fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def ncg_tsv(tmp_path: Path) -> Path:
    """Mimics NCG7 cancer driver annotation file."""
    content = (
        "ncg_symbol\tentrez\tcancer_type\tNCG_oncogene\tNCG_tsg\tmethod\tpmid\n"
        "TP53\t7157\tpan-cancer\t0\t1\tvogelstein\t8139\n"
        "KRAS\t3845\tpan-cancer\t1\t0\tvogelstein\t8139\n"
        "PIK3CA\t5290\tbreast\t1\t0\tscreen\t25079552\n"
        "BRCA1\t672\tbreast\t0\t1\tcgc\t7545954\n"
        "RB1\t5925\tretinoblastoma\t0\t1\tvogelstein\t8139\n"
        "CDK4\t1019\tmelanoma\t1\t0\tscreen\t25079552\n"
    )
    f = tmp_path / "NCG_cancerdrivers_annotation.tsv"
    f.write_text(content)
    return f


# ---------------------------------------------------------------------------
# CancerMine fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def cancermine_tsv(tmp_path: Path) -> Path:
    """Mimics cancermine_collated.tsv from Zenodo."""
    content = (
        "matching_id\tgene_hugo_id\tgene_entrez_id\tcancer_id\t"
        "cancer_normalized\trole\tcitation_count\n"
        "1\tTP53\t7157\tDOID:1612\tbreast cancer\tTumor_Suppressor\t500\n"
        "2\tTP53\t7157\tDOID:3571\tliver cancer\tTumor_Suppressor\t300\n"
        "3\tKRAS\t3845\tDOID:1324\tlung cancer\tOncogene\t450\n"
        "4\tPIK3CA\t5290\tDOID:1612\tbreast cancer\tOncogene\t200\n"
        "5\tBRAF\t673\tDOID:8923\tmelanoma\tOncogene\t350\n"
        "6\tPTEN\t5728\tDOID:1612\tbreast cancer\tTumor_Suppressor\t250\n"
        "7\tMYC\t4609\tDOID:9256\tcolorectal cancer\tOncogene\t180\n"
        "8\tBRCA1\t672\tDOID:1612\tbreast cancer\tTumor_Suppressor\t400\n"
        "9\tNOVEL_CM\t99999\tDOID:1324\tlung cancer\tDriver\t3\n"
    )
    f = tmp_path / "cancermine_collated.tsv"
    f.write_text(content)
    return f


# ===================================================================
# IntOGen loader tests
# ===================================================================

class TestLoadIntOGen:
    def test_returns_ground_truth_result(self, intogen_tsv: Path) -> None:
        result = load_intogen(intogen_tsv)
        assert isinstance(result, GroundTruthResult)

    def test_extracts_driver_genes(self, intogen_tsv: Path) -> None:
        result = load_intogen(intogen_tsv)
        assert "TP53" in result.positives
        assert "KRAS" in result.positives
        assert "PIK3CA" in result.positives

    def test_qvalue_filter(self, intogen_tsv: Path) -> None:
        """Should filter by q-value threshold."""
        result = load_intogen(intogen_tsv, max_qvalue=0.001)
        assert "NOVELGENE1" not in result.positives  # q=0.01
        assert "TP53" in result.positives  # q=1.2e-50

    def test_role_metadata(self, intogen_tsv: Path) -> None:
        result = load_intogen(intogen_tsv)
        meta = result.get_metadata("TP53")
        assert meta is not None
        assert meta["role"] == "TSG"  # LoF is normalized to TSG

    def test_missing_file_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            load_intogen(tmp_path / "nonexistent.tsv")

    def test_unique_genes(self, intogen_tsv: Path) -> None:
        """Same gene in multiple cohorts should appear once."""
        result = load_intogen(intogen_tsv)
        # TP53 appears once even if in multiple cohorts
        assert isinstance(result.positives, set)


# ===================================================================
# NCG loader tests
# ===================================================================

class TestLoadNCG:
    def test_returns_ground_truth_result(self, ncg_tsv: Path) -> None:
        result = load_ncg(ncg_tsv)
        assert isinstance(result, GroundTruthResult)

    def test_extracts_driver_genes(self, ncg_tsv: Path) -> None:
        result = load_ncg(ncg_tsv)
        assert "TP53" in result.positives
        assert "KRAS" in result.positives
        assert "BRCA1" in result.positives

    def test_role_metadata_oncogene(self, ncg_tsv: Path) -> None:
        result = load_ncg(ncg_tsv)
        meta = result.get_metadata("KRAS")
        assert meta is not None
        assert meta["role"] == "oncogene"

    def test_role_metadata_tsg(self, ncg_tsv: Path) -> None:
        result = load_ncg(ncg_tsv)
        meta = result.get_metadata("TP53")
        assert meta is not None
        assert meta["role"] == "TSG"

    def test_missing_file_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            load_ncg(tmp_path / "nonexistent.tsv")


# ===================================================================
# CancerMine loader tests
# ===================================================================

class TestLoadCancerMine:
    def test_returns_ground_truth_result(self, cancermine_tsv: Path) -> None:
        result = load_cancermine(cancermine_tsv)
        assert isinstance(result, GroundTruthResult)

    def test_extracts_driver_genes(self, cancermine_tsv: Path) -> None:
        result = load_cancermine(cancermine_tsv)
        assert "TP53" in result.positives
        assert "KRAS" in result.positives

    def test_citation_threshold(self, cancermine_tsv: Path) -> None:
        """Low-citation genes should be filtered out."""
        result = load_cancermine(cancermine_tsv, min_citations=10)
        assert "NOVEL_CM" not in result.positives  # Only 3 citations
        assert "TP53" in result.positives  # 500+300=800 citations

    def test_role_metadata(self, cancermine_tsv: Path) -> None:
        result = load_cancermine(cancermine_tsv)
        meta = result.get_metadata("TP53")
        assert meta is not None
        assert meta["role"] == "TSG"  # Tumor_Suppressor normalized to TSG

    def test_aggregates_across_cancers(self, cancermine_tsv: Path) -> None:
        """TP53 in breast + liver cancer → total citations = 800."""
        result = load_cancermine(cancermine_tsv)
        meta = result.get_metadata("TP53")
        assert meta["total_citations"] == 800

    def test_missing_file_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            load_cancermine(tmp_path / "nonexistent.tsv")


# ===================================================================
# Consensus ground truth tests
# ===================================================================

class TestConsensusGroundTruth:
    def test_returns_consensus(
        self, intogen_tsv: Path, ncg_tsv: Path, cancermine_tsv: Path
    ) -> None:
        result = build_consensus_ground_truth(
            intogen_path=intogen_tsv,
            ncg_path=ncg_tsv,
            cancermine_path=cancermine_tsv,
        )
        assert isinstance(result, ConsensusGroundTruth)

    def test_high_confidence_in_multiple_sources(
        self, intogen_tsv: Path, ncg_tsv: Path, cancermine_tsv: Path
    ) -> None:
        """TP53, KRAS, PIK3CA appear in all 3 sources → high confidence."""
        result = build_consensus_ground_truth(
            intogen_path=intogen_tsv,
            ncg_path=ncg_tsv,
            cancermine_path=cancermine_tsv,
        )
        assert "TP53" in result.high_confidence  # in all 3
        assert "KRAS" in result.high_confidence  # in all 3
        assert "PIK3CA" in result.high_confidence  # in all 3

    def test_min_sources_filter(
        self, intogen_tsv: Path, ncg_tsv: Path, cancermine_tsv: Path
    ) -> None:
        result = build_consensus_ground_truth(
            intogen_path=intogen_tsv,
            ncg_path=ncg_tsv,
            cancermine_path=cancermine_tsv,
            min_sources=2,
        )
        # Genes in >=2 sources
        for gene in result.positives:
            assert result.source_count[gene] >= 2

    def test_source_count_tracking(
        self, intogen_tsv: Path, ncg_tsv: Path, cancermine_tsv: Path
    ) -> None:
        result = build_consensus_ground_truth(
            intogen_path=intogen_tsv,
            ncg_path=ncg_tsv,
            cancermine_path=cancermine_tsv,
            min_sources=1,
        )
        # NOVELGENE1 only in IntOGen
        assert result.source_count.get("NOVELGENE1", 0) == 1
        # TP53 in all 3
        assert result.source_count["TP53"] == 3

    def test_as_labels_dict(
        self, intogen_tsv: Path, ncg_tsv: Path, cancermine_tsv: Path
    ) -> None:
        """Should provide a dict[str, int] for pipeline consumption."""
        result = build_consensus_ground_truth(
            intogen_path=intogen_tsv,
            ncg_path=ncg_tsv,
            cancermine_path=cancermine_tsv,
        )
        labels = result.as_labels()
        assert all(v == 1 for v in labels.values())
        assert "TP53" in labels

    def test_role_consensus(
        self, intogen_tsv: Path, ncg_tsv: Path, cancermine_tsv: Path
    ) -> None:
        """Should determine consensus role (oncogene/TSG)."""
        result = build_consensus_ground_truth(
            intogen_path=intogen_tsv,
            ncg_path=ncg_tsv,
            cancermine_path=cancermine_tsv,
        )
        assert result.consensus_role["TP53"] == "TSG"
        assert result.consensus_role["KRAS"] == "oncogene"


# ===================================================================
# Real data format tests (covering formats from actual downloads)
# ===================================================================

class TestRealDataFormats:
    """Test loaders with formats matching actual downloaded data."""

    def test_cancermine_gene_normalized_column(self, tmp_path: Path) -> None:
        """CancerMine real data uses gene_normalized for gene symbols."""
        content = (
            "matching_id\trole\tcancer_id\tcancer_normalized\t"
            "gene_hugo_id\tgene_entrez_id\tgene_normalized\tcitation_count\n"
            "abc\tOncogene\tDOID:1612\tbreast cancer\tHGNC:3430\t2064\tERBB2\t1280\n"
            "def\tTumor_Suppressor\tDOID:1612\tbreast cancer\tHGNC:11998\t7157\tTP53\t185\n"
            "ghi\tTumor_Suppressor\tDOID:684\thcc\tHGNC:11998\t7157\tTP53\t135\n"
            "jkl\tOncogene\tDOID:8923\tmelanoma\tHGNC:1097\t673\tBRAF\t686\n"
        )
        f = tmp_path / "cancermine_real.tsv"
        f.write_text(content)

        result = load_cancermine(f, min_citations=5)
        # Should use gene_normalized (TP53), not gene_hugo_id (HGNC:11998)
        assert "TP53" in result.positives
        assert "ERBB2" in result.positives
        assert "HGNC:11998" not in result.positives
        # TP53 citations = 185 + 135 = 320
        meta = result.get_metadata("TP53")
        assert meta is not None
        assert meta["total_citations"] == 320

    def test_ncg_springer_supplementary_format(self, tmp_path: Path) -> None:
        """NCG data from Springer supplementary has role column."""
        content = (
            "ncg_symbol\tdriver_annotation\trole\tis_canonical\tn_screens\n"
            "TP53\tcanonical cancer\tTumour Suppressor\tTrue\t150\n"
            "KRAS\tcanonical cancer\tOncogene\tTrue\t100\n"
            "ABL1\tcanonical cancer\tOncogene\tTrue\t10\n"
            "NOVEL1\tcandidate cancer\tunknown\tFalse\t1\n"
            "DUAL1\tcanonical cancer\tDual Role\tTrue\t5\n"
        )
        f = tmp_path / "ncg_springer.tsv"
        f.write_text(content)

        result = load_ncg(f)
        assert "TP53" in result.positives
        assert "KRAS" in result.positives
        assert "NOVEL1" in result.positives

        # Role mapping
        assert result.get_metadata("TP53")["role"] == "TSG"
        assert result.get_metadata("KRAS")["role"] == "oncogene"
        assert result.get_metadata("NOVEL1")["role"] == "unknown"
        assert result.get_metadata("DUAL1")["role"] == "both"

    def test_intogen_compendium_is_driver_filter(self, tmp_path: Path) -> None:
        """IntOGen Compendium has IS_DRIVER column to filter entries."""
        content = (
            "SYMBOL\tTRANSCRIPT\tCOHORT\tCANCER_TYPE\tMETHODS\tMUTATIONS\t"
            "SAMPLES\tQVALUE_COMBINATION\tROLE\tIS_DRIVER\n"
            "TP53\tENST001\tTCGA_BRCA\tBRCA\tdndscv\t100\t1000\t1e-50\tLoF\tTrue\n"
            "KRAS\tENST002\tTCGA_LUAD\tLUAD\thotmaps\t80\t500\t3e-30\tAct\tTrue\n"
            "FAKEGENE\tENST003\tTCGA_BRCA\tBRCA\tcbase\t2\t100\t0.04\tAct\tFalse\n"
        )
        f = tmp_path / "compendium.tsv"
        f.write_text(content)

        result = load_intogen(f)
        assert "TP53" in result.positives
        assert "KRAS" in result.positives
        assert "FAKEGENE" not in result.positives  # IS_DRIVER=False

    def test_intogen_unfiltered_driver_column(self, tmp_path: Path) -> None:
        """IntOGen Unfiltered_drivers has DRIVER column (True/False)."""
        content = (
            "SYMBOL\tTRANSCRIPT\tCOHORT\tCANCER_TYPE\tMUTATIONS\t"
            "QVALUE_COMBINATION\tROLE\tDRIVER\tFILTER\n"
            "TP53\tENST001\tTCGA_BRCA\tBRCA\t100\t1e-50\tLoF\tTrue\tPASS\n"
            "BRAF\tENST002\tTCGA_SKCM\tSKCM\t80\t2e-20\tAct\tTrue\tPASS\n"
            "NONGENE\tENST003\tTCGA_BRCA\tBRCA\t5\t0.04\tAct\tFalse\tNo driver\n"
        )
        f = tmp_path / "unfiltered.tsv"
        f.write_text(content)

        result = load_intogen(f)
        assert "TP53" in result.positives
        assert "BRAF" in result.positives
        assert "NONGENE" not in result.positives  # DRIVER=False
