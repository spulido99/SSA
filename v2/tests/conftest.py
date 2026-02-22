"""Shared test fixtures for SSA-ME v2."""

from pathlib import Path

import pytest


@pytest.fixture
def project_root() -> Path:
    """Return the v2 project root directory."""
    return Path(__file__).parent.parent


@pytest.fixture
def data_dir(project_root: Path) -> Path:
    """Return the data directory."""
    return project_root / "data"


@pytest.fixture
def v1_test_data() -> Path:
    """Return path to v1 test data (dendrix .m2 files)."""
    return Path(__file__).parent.parent.parent / "src" / "test" / "resources" / "dendrix"


@pytest.fixture
def fixtures_dir() -> Path:
    """Return path to test fixtures."""
    return Path(__file__).parent / "fixtures"
