"""Tests for gnomAD constraint loader (WP2.9). Written FIRST per TDD."""

from pathlib import Path

import pytest

from ssame.data.gnomad_loader import GnomADResult, load_gnomad_constraints


@pytest.fixture
def gnomad_tsv(tmp_path: Path) -> Path:
    """Synthetic gnomAD constraint TSV."""
    header = "gene\ttranscript\tobs_lof\texp_lof\toe_lof\toe_lof_upper\tpLI\tloeuf\n"
    rows = [
        "TP53\tENST00000269305\t2\t50.3\t0.04\t0.08\t1.00\t0.08\n",
        "PIK3CA\tENST00000263967\t3\t42.1\t0.07\t0.12\t0.99\t0.12\n",
        "TTN\tENST00000589042\t150\t180.0\t0.83\t0.95\t0.00\t0.95\n",
        "BRAF\tENST00000288602\t5\t38.5\t0.13\t0.22\t0.95\t0.22\n",
        "RANDOM1\tENST00000999999\t80\t85.0\t0.94\t1.10\t0.00\t1.10\n",
    ]
    tsv_path = tmp_path / "gnomad_constraints.tsv"
    tsv_path.write_text(header + "".join(rows))
    return tsv_path


class TestGnomADLoader:
    def test_loads_tsv(self, gnomad_tsv: Path) -> None:
        result = load_gnomad_constraints(gnomad_tsv)
        assert isinstance(result, GnomADResult)

    def test_pli_scores(self, gnomad_tsv: Path) -> None:
        result = load_gnomad_constraints(gnomad_tsv)
        assert result.pli["TP53"] == pytest.approx(1.0)
        assert result.pli["TTN"] == pytest.approx(0.0)

    def test_loeuf_scores(self, gnomad_tsv: Path) -> None:
        result = load_gnomad_constraints(gnomad_tsv)
        # Lower LOEUF = more constrained = more likely driver
        assert result.loeuf["TP53"] < result.loeuf["TTN"]

    def test_constrained_genes(self, gnomad_tsv: Path) -> None:
        """Genes with high pLI should be identified as constrained."""
        result = load_gnomad_constraints(gnomad_tsv)
        constrained = result.get_constrained_genes(pli_threshold=0.9)
        assert "TP53" in constrained
        assert "PIK3CA" in constrained
        assert "TTN" not in constrained

    def test_loeuf_constrained(self, gnomad_tsv: Path) -> None:
        """Genes with low LOEUF should be constrained."""
        result = load_gnomad_constraints(gnomad_tsv)
        constrained = result.get_constrained_genes(loeuf_threshold=0.35)
        assert "TP53" in constrained
        assert "BRAF" in constrained
        assert "TTN" not in constrained

    def test_all_genes_present(self, gnomad_tsv: Path) -> None:
        result = load_gnomad_constraints(gnomad_tsv)
        assert len(result.pli) == 5
        assert len(result.loeuf) == 5

    def test_feature_vector(self, gnomad_tsv: Path) -> None:
        """Should produce feature vectors for use in the encoder."""
        result = load_gnomad_constraints(gnomad_tsv)
        features = result.as_features(["TP53", "PIK3CA", "TTN"])
        # (3 genes, 2 features: pLI, LOEUF)
        assert features.shape == (3, 2)
        # TP53 pLI should be 1.0
        assert features[0, 0] == pytest.approx(1.0)
