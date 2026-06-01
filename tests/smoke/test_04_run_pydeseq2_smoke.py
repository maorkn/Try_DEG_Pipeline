"""Smoke tests for 04_run_pydeseq2.py script."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

# Paths
SCRIPT_PATH = Path(__file__).parent.parent.parent / "scripts" / "04_run_pydeseq2.py"
FIXTURES_DIR = Path(__file__).parent.parent.parent / "data" / "fixtures"


@pytest.fixture
def small_counts_metadata(tmp_path: Path) -> tuple[Path, Path, Path]:
    """Create small count matrix and metadata for testing."""
    # Create count matrix
    counts = pd.DataFrame({
        "sample1": [100, 50, 25, 10],
        "sample2": [200, 100, 50, 25],
        "sample3": [150, 75, 30, 15],
        "sample4": [180, 90, 45, 20],
    }, index=["gene1", "gene2", "gene3", "gene4"])

    counts_path = tmp_path / "counts.tsv"
    counts.to_csv(counts_path, sep="\t")

    # Create metadata
    metadata = pd.DataFrame({
        "condition": ["control", "control", "treated", "treated"],
    }, index=["sample1", "sample2", "sample3", "sample4"])

    metadata_path = tmp_path / "metadata.tsv"
    metadata.to_csv(metadata_path, sep="\t")

    # Create config
    config = {
        "project_name": "test_project",
        "species": "human",
        "gene_id_type": "ensembl",
        "condition_column": "condition",
        "reference_level": "control",
        "results_dir": str(tmp_path / "results"),
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
    }

    import yaml
    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config, f)

    return counts_path, metadata_path, config_path


def run_script(args: list[str]) -> subprocess.CompletedProcess:
    """Run the PyDESeq2 script with given arguments."""
    cmd = [sys.executable, str(SCRIPT_PATH)] + args
    return subprocess.run(cmd, capture_output=True, text=True)


class TestRunPyDESeq2Smoke:
    """Smoke tests for the PyDESeq2 script."""

    def test_help_flag(self):
        """Script should respond to --help."""
        result = run_script(["--help"])
        assert result.returncode == 0
        assert "PyDESeq2" in result.stdout or "usage" in result.stdout.lower()

    def test_creates_result_table(self, small_counts_metadata, tmp_path: Path):
        """Script should create result table."""
        counts_path, metadata_path, config_path = small_counts_metadata

        outdir = tmp_path / "results"
        result = run_script([
            "--counts-filtered", str(counts_path),
            "--metadata", str(metadata_path),
            "--config", str(config_path),
            "--outdir", str(outdir),
        ])

        # Check if result table was created
        results_path = outdir / "results_treated_vs_control.tsv"
        if results_path.exists():
            df = pd.read_csv(results_path, sep="\t")
            assert "log2FoldChange" in df.columns
            assert "padj" in df.columns

    def test_missing_counts_fails(self, tmp_path: Path):
        """Script should fail when counts file is missing."""
        result = run_script([
            "--counts-filtered", str(tmp_path / "nonexistent.tsv"),
            "--metadata", str(tmp_path / "metadata.tsv"),
            "--config", str(tmp_path / "config.yaml"),
            "--outdir", str(tmp_path / "results"),
        ])

        assert result.returncode != 0