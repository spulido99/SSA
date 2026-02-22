"""Tests for DepMap CRISPR essentiality loader (WP2.4). Written FIRST per TDD."""

from pathlib import Path

import pytest

from ssame.data.depmap_loader import DepMapResult, load_depmap_essentiality


@pytest.fixture
def depmap_csv(tmp_path: Path) -> Path:
    """Synthetic DepMap CRISPRGeneEffect CSV."""
    content = (
        ",TP53 (7157),PIK3CA (5290),KRAS (3845),CDH1 (999),BRAF (673)\n"
        "ACH-000001,-1.5,-0.8,-2.0,-0.3,0.1\n"
        "ACH-000002,-1.2,-1.0,-1.8,-0.5,-0.2\n"
        "ACH-000003,-1.8,-0.5,-1.5,-0.1,0.3\n"
        "ACH-000004,-0.9,-1.2,-2.2,-0.4,0.0\n"
    )
    csv_path = tmp_path / "CRISPRGeneEffect.csv"
    csv_path.write_text(content)
    return csv_path


class TestDepMapLoader:
    def test_loads_csv(self, depmap_csv: Path) -> None:
        result = load_depmap_essentiality(depmap_csv)
        assert isinstance(result, DepMapResult)

    def test_returns_per_gene_scores(self, depmap_csv: Path) -> None:
        result = load_depmap_essentiality(depmap_csv)
        assert "TP53" in result.mean_essentiality
        assert "PIK3CA" in result.mean_essentiality
        assert "KRAS" in result.mean_essentiality

    def test_essential_genes_have_negative_scores(self, depmap_csv: Path) -> None:
        """More essential genes should have more negative mean scores."""
        result = load_depmap_essentiality(depmap_csv)
        # KRAS is most essential in our synthetic data
        assert result.mean_essentiality["KRAS"] < result.mean_essentiality["BRAF"]

    def test_gene_names_cleaned(self, depmap_csv: Path) -> None:
        """Column format 'GENE (ENTREZ)' should be parsed to just 'GENE'."""
        result = load_depmap_essentiality(depmap_csv)
        assert "TP53" in result.mean_essentiality
        # Should NOT have the raw column name
        assert "TP53 (7157)" not in result.mean_essentiality

    def test_entrez_ids_extracted(self, depmap_csv: Path) -> None:
        """Entrez IDs should be extracted from column names."""
        result = load_depmap_essentiality(depmap_csv)
        assert result.entrez_ids["TP53"] == 7157
        assert result.entrez_ids["KRAS"] == 3845

    def test_per_cell_line_matrix(self, depmap_csv: Path) -> None:
        """Should provide raw matrix for per-cell-line analysis."""
        result = load_depmap_essentiality(depmap_csv)
        assert result.matrix.shape == (4, 5)  # 4 cell lines, 5 genes

    def test_selective_dependency(self, depmap_csv: Path) -> None:
        """Should compute selective dependency (variance across cell lines)."""
        result = load_depmap_essentiality(depmap_csv)
        # Genes with variable essentiality across cell lines = selective
        assert "TP53" in result.selective_dependency
        assert result.selective_dependency["TP53"] >= 0

    def test_threshold_filter(self, depmap_csv: Path) -> None:
        """Should filter essential genes by threshold."""
        result = load_depmap_essentiality(depmap_csv, threshold=-0.5)
        essential = result.get_essential_genes(threshold=-0.5)
        # TP53 and KRAS have mean < -0.5
        assert "TP53" in essential
        assert "KRAS" in essential
        # BRAF has positive mean
        assert "BRAF" not in essential
