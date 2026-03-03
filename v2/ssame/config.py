"""Configuration system for SSA-ME v2 pipeline."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import yaml

VALID_ENCODER_LEVELS = {0, 1, 2}
VALID_SEARCH_STRATEGIES = {"ant_colony", "beam_search", "hybrid"}
VALID_SCORERS = {"me", "functional", "structural", "metabolic"}
VALID_AGGREGATIONS = {"weighted_sum", "min", "vote"}

# Parameter presets for memory/speed trade-off
PRESETS: dict[str, dict[str, int]] = {
    "lite": {"n_iter": 50, "n_ants": 10, "n_estimators": 50},
    "standard": {"n_iter": 200, "n_ants": 20, "n_estimators": 100},
    "full": {"n_iter": 1000, "n_ants": 50, "n_estimators": 200},
}


@dataclass
class DataConfig:
    raw_dir: Path = field(default_factory=lambda: Path("data/raw"))
    processed_dir: Path = field(default_factory=lambda: Path("data/processed"))
    cache_dir: Path = field(default_factory=lambda: Path("data/cache"))
    networks: list[str] = field(
        default_factory=lambda: ["hi_union", "reactome", "alphafold_baker"]
    )


@dataclass
class EncoderConfig:
    level: int = 0
    pu_prior: float = 0.05
    hidden_dim: int = 128
    n_heads: int = 8
    n_layers: int = 3
    learning_rate: float = 1e-3
    epochs: int = 100
    dropout: float = 0.1

    def __post_init__(self) -> None:
        if self.level not in VALID_ENCODER_LEVELS:
            raise ValueError(
                f"Encoder level must be one of {VALID_ENCODER_LEVELS}, got {self.level}"
            )


@dataclass
class SearchConfig:
    strategy: str = "ant_colony"
    n_iterations: int = 1000
    n_ants: int = 50
    subnetwork_size: int = 3
    beam_width: int = 10
    evaporation_rate: float = 0.1

    def __post_init__(self) -> None:
        if self.strategy not in VALID_SEARCH_STRATEGIES:
            raise ValueError(
                f"Search strategy must be one of {VALID_SEARCH_STRATEGIES}, got {self.strategy}"
            )


@dataclass
class ScorerConfig:
    scorers: list[str] = field(default_factory=lambda: ["me"])
    aggregation: str = "weighted_sum"
    weights: dict[str, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.aggregation not in VALID_AGGREGATIONS:
            raise ValueError(
                f"Aggregation must be one of {VALID_AGGREGATIONS}, got {self.aggregation}"
            )


@dataclass
class PipelineConfig:
    cancer_type: str = "BRCA"
    data: DataConfig = field(default_factory=DataConfig)
    encoder: EncoderConfig = field(default_factory=EncoderConfig)
    search: SearchConfig = field(default_factory=SearchConfig)
    scorer: ScorerConfig = field(default_factory=ScorerConfig)
    seed: int = 42


def _build_sub_config(cls: type, raw: dict) -> object:
    """Build a dataclass config from a dict, converting Path fields."""
    import dataclasses

    field_types = {f.name: f.type for f in dataclasses.fields(cls)}
    kwargs = {}
    for k, v in raw.items():
        if k in field_types and "Path" in str(field_types[k]):
            v = Path(v)
        kwargs[k] = v
    return cls(**kwargs)


def load_config(path: Path) -> PipelineConfig:
    """Load pipeline config from a YAML file."""
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")

    with open(path) as f:
        raw = yaml.safe_load(f) or {}

    data_cfg = _build_sub_config(DataConfig, raw.get("data", {}))
    encoder_cfg = _build_sub_config(EncoderConfig, raw.get("encoder", {}))
    search_cfg = _build_sub_config(SearchConfig, raw.get("search", {}))
    scorer_cfg = _build_sub_config(ScorerConfig, raw.get("scorer", {}))

    return PipelineConfig(
        cancer_type=raw.get("cancer_type", "BRCA"),
        data=data_cfg,  # type: ignore[arg-type]
        encoder=encoder_cfg,  # type: ignore[arg-type]
        search=search_cfg,  # type: ignore[arg-type]
        scorer=scorer_cfg,  # type: ignore[arg-type]
        seed=raw.get("seed", 42),
    )
