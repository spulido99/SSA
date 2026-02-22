"""Tests for alteration matrix parsing and construction (WP2.2). Written FIRST per TDD."""

from pathlib import Path

import numpy as np
import pytest

from ssame.data.alteration_matrix import (
    export_m2,
    load_m2,
    parse_maf_to_matrix,
)
from ssame.types import AlterationMatrix


class TestLoadM2:
    """Test loading v1-format .m2 files."""

    def test_load_brca_m2(self, v1_test_data: Path) -> None:
        brca_file = v1_test_data / "BRCA.m2"
        matrix = load_m2(brca_file)
        assert isinstance(matrix, AlterationMatrix)
        assert matrix.n_samples > 0
        assert matrix.n_genes > 0

    def test_m2_contains_known_drivers(self, v1_test_data: Path) -> None:
        brca_file = v1_test_data / "BRCA.m2"
        matrix = load_m2(brca_file)
        # TP53 and PIK3CA are the most frequent BRCA drivers
        assert "TP53" in matrix.genes
        assert "PIK3CA" in matrix.genes

    def test_m2_sample_ids_are_tcga_format(self, v1_test_data: Path) -> None:
        brca_file = v1_test_data / "BRCA.m2"
        matrix = load_m2(brca_file)
        # TCGA sample IDs start with "TCGA-"
        assert all(s.startswith("TCGA-") for s in matrix.samples)

    def test_m2_matrix_is_binary(self, v1_test_data: Path) -> None:
        brca_file = v1_test_data / "BRCA.m2"
        matrix = load_m2(brca_file)
        unique_vals = np.unique(matrix.data)
        assert set(unique_vals).issubset({0, 1})

    def test_m2_strips_alteration_type_annotations(self, v1_test_data: Path) -> None:
        """Gene names like PAK1(A), PTEN(D) should have annotation stripped."""
        brca_file = v1_test_data / "BRCA.m2"
        matrix = load_m2(brca_file)
        # No gene name should end with (A) or (D)
        for gene in matrix.genes:
            assert not gene.endswith("(A)"), f"Gene {gene} still has (A) annotation"
            assert not gene.endswith("(D)"), f"Gene {gene} still has (D) annotation"

    def test_m2_skips_comment_lines(self, v1_test_data: Path) -> None:
        brca_file = v1_test_data / "BRCA.m2"
        matrix = load_m2(brca_file)
        # Comment line should not be a sample
        assert not any(s.startswith("#") for s in matrix.samples)

    def test_load_gbm_m2(self, v1_test_data: Path) -> None:
        gbm_file = v1_test_data / "GBM.m2"
        matrix = load_m2(gbm_file)
        assert matrix.n_samples > 0
        assert matrix.n_genes > 0

    def test_load_lung_m2(self, v1_test_data: Path) -> None:
        lung_file = v1_test_data / "Lung.m2"
        matrix = load_m2(lung_file)
        assert matrix.n_samples > 0


class TestLoadM2WithAnnotations:
    """Test that alteration type annotations (A/D) are tracked."""

    def test_load_with_alteration_types(self, v1_test_data: Path) -> None:
        brca_file = v1_test_data / "BRCA.m2"
        _matrix, annotations = load_m2(brca_file, return_annotations=True)
        assert isinstance(annotations, dict)
        # Should have annotation entries for genes that had (A) or (D)
        assert len(annotations) > 0

    def test_annotations_contain_amp_and_del(self, v1_test_data: Path) -> None:
        brca_file = v1_test_data / "BRCA.m2"
        _, annotations = load_m2(brca_file, return_annotations=True)
        all_types = set()
        for types_per_sample in annotations.values():
            for type_val in types_per_sample.values():
                all_types.add(type_val)
        # Should have A (amplification) and D (deletion) types
        assert "A" in all_types or "D" in all_types


class TestExportM2:
    """Test exporting back to .m2 format."""

    def test_roundtrip(self, v1_test_data: Path, tmp_path: Path) -> None:
        brca_file = v1_test_data / "BRCA.m2"
        original = load_m2(brca_file)
        output_file = tmp_path / "roundtrip.m2"
        export_m2(original, output_file)
        reloaded = load_m2(output_file)
        assert reloaded.n_samples == original.n_samples
        assert reloaded.n_genes == original.n_genes


class TestParseMafToMatrix:
    """Test MAF file parsing to alteration matrix."""

    def test_parse_minimal_maf(self, tmp_path: Path) -> None:
        """Test with a minimal synthetic MAF file."""
        maf_content = (
            "Hugo_Symbol\tTumor_Sample_Barcode\tVariant_Classification\t"
            "Chromosome\tStart_Position\tEnd_Position\tReference_Allele\t"
            "Tumor_Seq_Allele2\n"
            "TP53\tTCGA-A1-0001\tMissense_Mutation\tchr17\t7577120\t7577120\tC\tT\n"
            "PIK3CA\tTCGA-A1-0001\tMissense_Mutation\tchr3\t178936091\t178936091\tA\tG\n"
            "TP53\tTCGA-A1-0002\tNonsense_Mutation\tchr17\t7578190\t7578190\tG\tA\n"
            "CDH1\tTCGA-A1-0002\tFrame_Shift_Del\tchr16\t68811195\t68811195\tA\t-\n"
            "TP53\tTCGA-A1-0003\tSplice_Site\tchr17\t7577019\t7577019\tG\tT\n"
        )
        maf_file = tmp_path / "test.maf"
        maf_file.write_text(maf_content)
        matrix = parse_maf_to_matrix(maf_file)
        assert isinstance(matrix, AlterationMatrix)
        assert "TP53" in matrix.genes
        assert "PIK3CA" in matrix.genes
        assert "CDH1" in matrix.genes
        assert matrix.n_samples == 3
        # TP53 mutated in all 3 samples
        assert matrix.gene_frequency("TP53") == pytest.approx(1.0)
        # PIK3CA only in sample 1
        assert matrix.gene_frequency("PIK3CA") == pytest.approx(1 / 3)

    def test_parse_maf_filters_silent_mutations(self, tmp_path: Path) -> None:
        maf_content = (
            "Hugo_Symbol\tTumor_Sample_Barcode\tVariant_Classification\t"
            "Chromosome\tStart_Position\tEnd_Position\tReference_Allele\t"
            "Tumor_Seq_Allele2\n"
            "TP53\tTCGA-A1-0001\tMissense_Mutation\tchr17\t7577120\t7577120\tC\tT\n"
            "SILENT_GENE\tTCGA-A1-0001\tSilent\tchr1\t100\t100\tA\tG\n"
        )
        maf_file = tmp_path / "test.maf"
        maf_file.write_text(maf_content)
        matrix = parse_maf_to_matrix(maf_file)
        assert "SILENT_GENE" not in matrix.genes

    def test_parse_maf_min_samples(self, tmp_path: Path) -> None:
        maf_content = (
            "Hugo_Symbol\tTumor_Sample_Barcode\tVariant_Classification\t"
            "Chromosome\tStart_Position\tEnd_Position\tReference_Allele\t"
            "Tumor_Seq_Allele2\n"
            "TP53\tTCGA-A1-0001\tMissense_Mutation\tchr17\t1\t1\tC\tT\n"
            "TP53\tTCGA-A1-0002\tMissense_Mutation\tchr17\t1\t1\tC\tT\n"
            "TP53\tTCGA-A1-0003\tMissense_Mutation\tchr17\t1\t1\tC\tT\n"
            "RARE\tTCGA-A1-0001\tMissense_Mutation\tchr1\t1\t1\tA\tG\n"
        )
        maf_file = tmp_path / "test.maf"
        maf_file.write_text(maf_content)
        matrix = parse_maf_to_matrix(maf_file, min_samples=2)
        assert "TP53" in matrix.genes
        assert "RARE" not in matrix.genes

    def test_parse_maf_sample_barcodes_filter(self, tmp_path: Path) -> None:
        """Only mutations in specified sample_barcodes are included."""
        maf_content = (
            "Hugo_Symbol\tTumor_Sample_Barcode\tVariant_Classification\t"
            "Chromosome\tStart_Position\tEnd_Position\tReference_Allele\t"
            "Tumor_Seq_Allele2\n"
            "TP53\tTCGA-A1-0001\tMissense_Mutation\tchr17\t1\t1\tC\tT\n"
            "KRAS\tTCGA-A1-0001\tMissense_Mutation\tchr12\t1\t1\tG\tA\n"
            "TP53\tTCGA-A1-0002\tNonsense_Mutation\tchr17\t1\t1\tG\tA\n"
            "BRAF\tTCGA-A1-0002\tMissense_Mutation\tchr7\t1\t1\tT\tA\n"
            "EGFR\tTCGA-A1-0003\tMissense_Mutation\tchr7\t1\t1\tC\tT\n"
        )
        maf_file = tmp_path / "test.maf"
        maf_file.write_text(maf_content)

        # Only keep samples 0001 and 0003
        subset = {"TCGA-A1-0001", "TCGA-A1-0003"}
        matrix = parse_maf_to_matrix(maf_file, sample_barcodes=subset)

        assert matrix.n_samples == 2
        assert "TCGA-A1-0001" in matrix.samples
        assert "TCGA-A1-0003" in matrix.samples
        assert "TCGA-A1-0002" not in matrix.samples

        # BRAF only in excluded sample 0002
        assert "BRAF" not in matrix.genes
        # TP53 still present (in sample 0001)
        assert "TP53" in matrix.genes
        # EGFR still present (in sample 0003)
        assert "EGFR" in matrix.genes

    def test_parse_maf_sample_barcodes_empty_set(self, tmp_path: Path) -> None:
        """Empty sample_barcodes set returns empty matrix."""
        maf_content = (
            "Hugo_Symbol\tTumor_Sample_Barcode\tVariant_Classification\t"
            "Chromosome\tStart_Position\tEnd_Position\tReference_Allele\t"
            "Tumor_Seq_Allele2\n"
            "TP53\tTCGA-A1-0001\tMissense_Mutation\tchr17\t1\t1\tC\tT\n"
        )
        maf_file = tmp_path / "test.maf"
        maf_file.write_text(maf_content)

        matrix = parse_maf_to_matrix(maf_file, sample_barcodes=set())
        assert matrix.n_genes == 0
        assert matrix.n_samples == 0
