"""Smoke tests for 01_validate_inputs.py script."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

# Paths
SCRIPT_PATH = Path(__file__).parent.parent.parent / "scripts" / "01_validate_inputs.py"
FIXTURES_DIR = Path(__file__).parent.parent.parent / "data" / "fixtures"


@pytest.fixture
def toy_config_file(tmp_path: Path) -> Path:
    """Create a toy config file for smoke tests."""
    config = {
        "project_name": "smoke_test",
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
        "filtering": {"min_count": 10, "min_samples": 2},
        "pydeseq2": {"alpha": 0.05, "cooks_filter": True, "independent_filter": True},
        "deg": {"padj_threshold": 0.05, "log2fc_threshold": 1.0},
        "go": {
            "ontology": ["BP"],
            "method": "offline_ora",
            "padj_method": "fdr_bh",
            "min_genes_per_term": 2,
            "max_genes_per_term": 100,
        },
    }
    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config, f)
    return config_path


def run_script(args: list[str]) -> subprocess.CompletedProcess:
    """Run the validation script with given arguments."""
    cmd = [sys.executable, str(SCRIPT_PATH)] + args
    return subprocess.run(cmd, capture_output=True, text=True)


class TestValidateInputsSmoke:
    """Smoke tests for the validate inputs script."""

    def test_help_flag(self):
        """Script should respond to --help."""
        result = run_script(["--help"])
        assert result.returncode == 0
        assert "Validate input files" in result.stdout or "usage" in result.stdout.lower()

    def test_valid_inputs_creates_outputs(self, tmp_path: Path, toy_config_file: Path):
        """Script should create expected output files with valid inputs."""
        outdir = tmp_path / "results" / "smoke_test"
        result = run_script([
            "--counts", str(FIXTURES_DIR / "toy_counts.tsv"),
            "--metadata", str(FIXTURES_DIR / "toy_metadata.tsv"),
            "--config", str(toy_config_file),
            "--annotation", str(FIXTURES_DIR / "toy_annotation.tsv"),
            "--outdir", str(outdir),
        ])

        assert result.returncode == 0, f"Script failed: {result.stderr}"

        # Check output files exist
        assert (outdir / "validation" / "input_validation.json").exists()
        assert (outdir / "validation" / "input_validation.txt").exists()
        assert (outdir / "intermediate" / "counts_ordered.tsv").exists()
        assert (outdir / "intermediate" / "metadata_ordered.tsv").exists()
        assert (outdir / "validation" / "manifest.json").exists()

    def test_validation_json_has_correct_structure(
        self, tmp_path: Path, toy_config_file: Path
    ):
        """Validation JSON should have expected structure."""
        outdir = tmp_path / "results" / "smoke_test"
        result = run_script([
            "--counts", str(FIXTURES_DIR / "toy_counts.tsv"),
            "--metadata", str(FIXTURES_DIR / "toy_metadata.tsv"),
            "--config", str(toy_config_file),
            "--outdir", str(outdir),
        ])

        assert result.returncode == 0, f"Script failed: {result.stderr}"

        json_path = outdir / "validation" / "input_validation.json"
        with open(json_path) as f:
            data = json.load(f)

        assert "is_valid" in data
        assert "errors" in data
        assert "warnings" in data
        assert "metrics" in data
        assert data["is_valid"] is True

    def test_ordered_metadata_matches_counts(
        self, tmp_path: Path, toy_config_file: Path
    ):
        """Ordered metadata sample order should match counts."""
        import pandas as pd

        outdir = tmp_path / "results" / "smoke_test"
        result = run_script([
            "--counts", str(FIXTURES_DIR / "toy_counts.tsv"),
            "--metadata", str(FIXTURES_DIR / "toy_metadata.tsv"),
            "--config", str(toy_config_file),
            "--outdir", str(outdir),
        ])

        assert result.returncode == 0

        counts = pd.read_csv(
            outdir / "intermediate" / "counts_ordered.tsv", sep="\t", index_col=0
        )
        metadata = pd.read_csv(
            outdir / "intermediate" / "metadata_ordered.tsv", sep="\t", index_col=0
        )

        assert list(counts.columns) == list(metadata.index)

    def test_missing_counts_fails(self, tmp_path: Path, toy_config_file: Path):
        """Script should fail when counts file is missing."""
        outdir = tmp_path / "results" / "smoke_test"
        result = run_script([
            "--counts", str(tmp_path / "nonexistent_counts.tsv"),
            "--metadata", str(FIXTURES_DIR / "toy_metadata.tsv"),
            "--config", str(toy_config_file),
            "--outdir", str(outdir),
        ])

        assert result.returncode != 0

    def test_invalid_counts_fails(self, tmp_path: Path, toy_config_file: Path):
        """Script should fail when counts have invalid values."""
        # Create counts file with negative values
        bad_counts = tmp_path / "bad_counts.tsv"
        bad_counts.write_text(
            "gene_id\tsample_1\tsample_2\n"
            "gene_1\t-100\t50\n"
            "gene_2\t120\t60\n"
        )

        outdir = tmp_path / "results" / "smoke_test"
        result = run_script([
            "--counts", str(bad_counts),
            "--metadata", str(FIXTURES_DIR / "toy_metadata.tsv"),
            "--config", str(toy_config_file),
            "--outdir", str(outdir),
        ])

        assert result.returncode != 0

        # Check that validation report indicates failure
        json_path = outdir / "validation" / "input_validation.json"
        if json_path.exists():
            with open(json_path) as f:
                data = json.load(f)
            assert data["is_valid"] is False