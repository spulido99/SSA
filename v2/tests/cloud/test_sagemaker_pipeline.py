"""Tests for SageMaker pipeline definition (no AWS dependencies needed)."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from scripts.cloud.sagemaker_pipeline import build_pipeline_definition, select_instance


class TestSelectInstance:
    def test_gpu_for_l1(self) -> None:
        assert select_instance(100, 1000, encoder_level=1) == "ml.g4dn.xlarge"

    def test_large_for_small_l0(self) -> None:
        assert select_instance(200, 1000, encoder_level=0) == "ml.m5.large"

    def test_xlarge_for_medium_l0(self) -> None:
        # 5000 genes * 500 samples * 4 bytes * 3 = ~30MB — should be large
        assert select_instance(500, 5000, encoder_level=0) == "ml.m5.large"

    def test_2xlarge_for_huge_l0(self) -> None:
        # Very large: 20000 * 2000 * 4 * 3 = ~480MB * 3 = ~1.4GB -> still large
        # Need extremely large: gene * sample to push above 8GB
        inst = select_instance(50000, 50000, encoder_level=0)
        assert inst == "ml.m5.2xlarge"


class TestBuildPipelineDefinition:
    def test_returns_dict_with_steps(self) -> None:
        cancer_types = [
            {"name": "BRCA", "n_samples": 1000, "n_genes": 5000},
            {"name": "LUAD", "n_samples": 600, "n_genes": 4000},
        ]
        defn = build_pipeline_definition(
            role="arn:aws:iam::123:role/test",
            bucket="test-bucket",
            cancer_types=cancer_types,
        )
        assert defn["pipeline_name"] == "ssa-me-pan-cancer"
        # Should have: ParseMAF + 2 cancer types + Aggregate = 4 steps
        assert len(defn["steps"]) == 4

    def test_first_step_is_parse_maf(self) -> None:
        defn = build_pipeline_definition(
            role="arn:aws:iam::123:role/test",
            bucket="test-bucket",
            cancer_types=[{"name": "BRCA", "n_samples": 100, "n_genes": 1000}],
        )
        assert defn["steps"][0]["name"] == "ParseMAF"
        assert defn["steps"][0]["depends_on"] == []

    def test_last_step_is_aggregate(self) -> None:
        cancer_types = [
            {"name": "BRCA", "n_samples": 100, "n_genes": 1000},
            {"name": "LUAD", "n_samples": 100, "n_genes": 1000},
        ]
        defn = build_pipeline_definition(
            role="arn:aws:iam::123:role/test",
            bucket="test-bucket",
            cancer_types=cancer_types,
        )
        agg = defn["steps"][-1]
        assert agg["name"] == "Aggregate"
        assert "Run_BRCA" in agg["depends_on"]
        assert "Run_LUAD" in agg["depends_on"]

    def test_cancer_steps_depend_on_parse_maf(self) -> None:
        defn = build_pipeline_definition(
            role="arn:aws:iam::123:role/test",
            bucket="test-bucket",
            cancer_types=[
                {"name": "BRCA", "n_samples": 100, "n_genes": 1000},
                {"name": "COAD", "n_samples": 100, "n_genes": 1000},
            ],
        )
        for step in defn["steps"][1:-1]:
            assert step["depends_on"] == ["ParseMAF"]

    def test_parameters_stored(self) -> None:
        defn = build_pipeline_definition(
            role="arn:aws:iam::123:role/test",
            bucket="test-bucket",
            cancer_types=[{"name": "BRCA", "n_samples": 100, "n_genes": 1000}],
            encoder_level=1,
            n_iter=500,
            n_ants=30,
        )
        assert defn["parameters"]["EncoderLevel"] == 1
        assert defn["parameters"]["NumACOIterations"] == 500
        assert defn["parameters"]["NumAnts"] == 30

    def test_spot_enabled_by_default(self) -> None:
        defn = build_pipeline_definition(
            role="arn:aws:iam::123:role/test",
            bucket="test-bucket",
            cancer_types=[{"name": "BRCA", "n_samples": 100, "n_genes": 1000}],
        )
        cancer_step = defn["steps"][1]
        assert cancer_step.get("use_spot") is True

    def test_gpu_instance_for_l1(self) -> None:
        defn = build_pipeline_definition(
            role="arn:aws:iam::123:role/test",
            bucket="test-bucket",
            cancer_types=[{"name": "BRCA", "n_samples": 100, "n_genes": 1000}],
            encoder_level=1,
        )
        cancer_step = defn["steps"][1]
        assert "g4dn" in cancer_step["instance_type"]
