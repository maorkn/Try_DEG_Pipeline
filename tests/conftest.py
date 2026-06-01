"""Pytest configuration and shared fixtures for the DEG pipeline tests."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
import pytest
import yaml


# Paths to test fixtures
FIXTURES_DIR = Path(__file__).parent.parent / "data" / "fixtures"


@pytest.fixture
def fixtures_dir() -> Path:
    """Return the path to the test fixtures directory."""
    return FIXTURES_DIR


@pytest.fixture
def toy_counts_path(fixtures_dir: Path) -> Path:
    """Return the path to the toy counts fixture."""
    return fixtures_dir / "toy_counts.tsv"


@pytest.fixture
def toy_metadata_path(fixtures_dir: Path) -> Path:
    """Return the path to the toy metadata fixture."""
    return fixtures_dir / "toy_metadata.tsv"


@pytest.fixture
def toy_annotation_path(fixtures_dir: Path) -> Path:
    """Return the path to the toy annotation fixture."""
    return fixtures_dir / "toy_annotation.tsv"


@pytest.fixture
def toy_go_map_path(fixtures_dir: Path) -> Path:
    """Return the path to the toy GO map fixture."""
    return fixtures_dir / "toy_go_map.tsv"


@pytest.fixture
def toy_counts_df(toy_counts_path: Path) -> pd.DataFrame:
    """Load and return the toy counts DataFrame."""
    return pd.read_csv(toy_counts_path, sep="\t", index_col=0)


@pytest.fixture
def toy_metadata_df(toy_metadata_path: Path) -> pd.DataFrame:
    """Load and return the toy metadata DataFrame."""
    df = pd.read_csv(toy_metadata_path, sep="\t")
    return df.set_index("sample_id")


@pytest.fixture
def toy_annotation_df(toy_annotation_path: Path) -> pd.DataFrame:
    """Load and return the toy annotation DataFrame."""
    df = pd.read_csv(toy_annotation_path, sep="\t")
    return df.set_index("gene_id")


@pytest.fixture
def toy_config() -> dict[str, Any]:
    """Return a minimal toy configuration for testing."""
    return {
        "project_name": "test_project",
        "species": "human",
        "gene_id_type": "ensembl",
        "condition_column": "condition",
        "reference_level": "control",
        "results_dir": "results/test_project",
        "design": {
            "formula": "~ condition",
            "variables": ["condition"],
        },
        "contrasts": [
            {
                "name": "treated_vs_control",
                "variable": "condition",
                "numerator": "treated",
                "denominator": "control",
            }
        ],
        "filtering": {
            "min_count": 10,
            "min_samples": 2,
        },
        "pydeseq2": {
            "alpha": 0.05,
            "cooks_filter": True,
            "independent_filter": True,
        },
        "deg": {
            "padj_threshold": 0.05,
            "log2fc_threshold": 1.0,
        },
        "go": {
            "ontology": ["BP"],
            "method": "offline_ora",
            "padj_method": "fdr_bh",
            "min_genes_per_term": 2,
            "max_genes_per_term": 100,
        },
    }


@pytest.fixture
def toy_config_file(tmp_path: Path, toy_config: dict[str, Any]) -> Path:
    """Write a toy config to a temp file and return the path."""
    config_path = tmp_path / "test_config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(toy_config, f)
    return config_path


@pytest.fixture
def minimal_counts_df() -> pd.DataFrame:
    """Create a minimal counts DataFrame for unit tests."""
    return pd.DataFrame(
        {
            "sample_1": [100, 50, 0, 10],
            "sample_2": [120, 55, 1, 12],
            "sample_3": [10, 60, 200, 8],
            "sample_4": [12, 52, 210, 15],
        },
        index=["gene_high", "gene_med", "gene_diff", "gene_low"],
    )


@pytest.fixture
def minimal_metadata_df() -> pd.DataFrame:
    """Create a minimal metadata DataFrame for unit tests."""
    return pd.DataFrame(
        {
            "condition": ["control", "control", "treated", "treated"],
        },
        index=pd.Index(["sample_1", "sample_2", "sample_3", "sample_4"], name="sample_id"),
    )


# Markers for test selection
def pytest_configure(config: pytest.Config) -> None:
    """Register custom pytest markers."""
    config.addinivalue_line("markers", "requires_r: marks tests that require R")
    config.addinivalue_line("markers", "airway: marks tests using airway dataset")
    config.addinivalue_line("markers", "slow: marks slow tests")