"""Tests for configuration system (WP0.6). Written FIRST per TDD."""

from pathlib import Path
from textwrap import dedent

import pytest

from ssame.config import (
    DataConfig,
    EncoderConfig,
    PipelineConfig,
    ScorerConfig,
    SearchConfig,
    load_config,
)


class TestDataConfig:
    def test_defaults(self) -> None:
        cfg = DataConfig()
        assert cfg.raw_dir == Path("data/raw")
        assert cfg.processed_dir == Path("data/processed")
        assert cfg.cache_dir == Path("data/cache")

    def test_custom_paths(self) -> None:
        cfg = DataConfig(raw_dir=Path("/custom/raw"))
        assert cfg.raw_dir == Path("/custom/raw")


class TestEncoderConfig:
    def test_defaults(self) -> None:
        cfg = EncoderConfig()
        assert cfg.level == 0
        assert cfg.pu_prior == 0.05
        assert cfg.hidden_dim == 128
        assert cfg.n_heads == 8
        assert cfg.n_layers == 3
        assert cfg.learning_rate == 1e-3
        assert cfg.epochs == 100

    def test_level_validation(self) -> None:
        with pytest.raises(ValueError, match="level"):
            EncoderConfig(level=5)


class TestSearchConfig:
    def test_defaults(self) -> None:
        cfg = SearchConfig()
        assert cfg.strategy == "ant_colony"
        assert cfg.n_iterations == 1000
        assert cfg.n_ants == 50
        assert cfg.subnetwork_size == 3

    def test_strategy_validation(self) -> None:
        with pytest.raises(ValueError, match="strategy"):
            SearchConfig(strategy="invalid_strategy")


class TestScorerConfig:
    def test_defaults(self) -> None:
        cfg = ScorerConfig()
        assert "me" in cfg.scorers
        assert cfg.aggregation == "weighted_sum"

    def test_custom_scorers(self) -> None:
        cfg = ScorerConfig(scorers=["me", "functional"], aggregation="min")
        assert cfg.aggregation == "min"


class TestPipelineConfig:
    def test_defaults(self) -> None:
        cfg = PipelineConfig()
        assert isinstance(cfg.data, DataConfig)
        assert isinstance(cfg.encoder, EncoderConfig)
        assert isinstance(cfg.search, SearchConfig)
        assert isinstance(cfg.scorer, ScorerConfig)

    def test_cancer_type(self) -> None:
        cfg = PipelineConfig(cancer_type="BRCA")
        assert cfg.cancer_type == "BRCA"


class TestLoadConfig:
    def test_load_from_yaml(self, tmp_path: Path) -> None:
        yaml_content = dedent("""\
            cancer_type: BRCA
            encoder:
              level: 1
              hidden_dim: 256
              epochs: 200
            search:
              strategy: beam_search
              n_iterations: 500
            scorer:
              scorers:
                - me
                - functional
              aggregation: min
        """)
        config_file = tmp_path / "config.yaml"
        config_file.write_text(yaml_content)
        cfg = load_config(config_file)
        assert cfg.cancer_type == "BRCA"
        assert cfg.encoder.level == 1
        assert cfg.encoder.hidden_dim == 256
        assert cfg.encoder.epochs == 200
        assert cfg.search.strategy == "beam_search"
        assert cfg.scorer.aggregation == "min"

    def test_load_preserves_defaults(self, tmp_path: Path) -> None:
        yaml_content = dedent("""\
            cancer_type: GBM
        """)
        config_file = tmp_path / "config.yaml"
        config_file.write_text(yaml_content)
        cfg = load_config(config_file)
        assert cfg.cancer_type == "GBM"
        assert cfg.encoder.level == 0  # default preserved
        assert cfg.search.strategy == "ant_colony"  # default preserved

    def test_load_nonexistent_file(self) -> None:
        with pytest.raises(FileNotFoundError):
            load_config(Path("/nonexistent/config.yaml"))

    def test_load_invalid_yaml_values(self, tmp_path: Path) -> None:
        yaml_content = dedent("""\
            encoder:
              level: 99
        """)
        config_file = tmp_path / "config.yaml"
        config_file.write_text(yaml_content)
        with pytest.raises(ValueError):
            load_config(config_file)
