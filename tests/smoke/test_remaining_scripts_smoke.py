"""Smoke tests for pipeline scripts not covered elsewhere."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).parent.parent.parent
SCRIPTS = ROOT / "scripts"


def run_script(script_name: str, args: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(SCRIPTS / script_name), *args],
        capture_output=True,
        text=True,
    )


def test_00_fetch_airway_help():
    result = run_script("00_fetch_airway.py", ["--help"])
    assert result.returncode == 0
    assert "airway" in result.stdout.lower()


def test_02_qc_counts_cli_writes_outputs(tmp_path, toy_counts_path, toy_metadata_path, toy_config_file):
    outdir = tmp_path / "results"
    result = run_script(
        "02_qc_counts.py",
        [
            "--counts",
            str(toy_counts_path),
            "--metadata",
            str(toy_metadata_path),
            "--config",
            str(toy_config_file),
            "--outdir",
            str(outdir),
        ],
    )

    assert result.returncode == 0, result.stderr
    assert (outdir / "qc" / "library_sizes.tsv").exists()
    assert (outdir / "qc" / "pca_coordinates.tsv").exists()


def test_03_filter_normalize_cli_writes_outputs(
    tmp_path, toy_counts_path, toy_metadata_path, toy_config_file
):
    outdir = tmp_path / "results"
    result = run_script(
        "03_filter_normalize.py",
        [
            "--counts",
            str(toy_counts_path),
            "--metadata",
            str(toy_metadata_path),
            "--config",
            str(toy_config_file),
            "--outdir",
            str(outdir),
        ],
    )

    assert result.returncode == 0, result.stderr
    assert (outdir / "intermediate" / "counts_filtered.tsv").exists()
    assert (outdir / "filtering" / "gene_filtering_manifest.json").exists()


def test_05_validate_airway_pydeseq2_cli_passes_on_tiny_reference(tmp_path, toy_metadata_path):
    pipeline_dir = tmp_path / "pipeline"
    reference_dir = tmp_path / "reference"
    pipeline_dir.mkdir()
    reference_dir.mkdir()
    table = (
        "gene_id\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\n"
        "g1\t10\t1.0\t0.1\t10\t0.001\t0.01\n"
        "g2\t20\t-1.0\t0.1\t-10\t0.001\t0.01\n"
        "g3\t30\t0.1\t0.1\t1\t0.3\t0.5\n"
    )
    (pipeline_dir / "results_toy.tsv").write_text(table)
    (reference_dir / "reference_deseq2_toy.tsv").write_text(table)
    outdir = tmp_path / "out"

    result = run_script(
        "05_validate_airway_pydeseq2.py",
        [
            "--pipeline-results",
            str(pipeline_dir),
            "--reference-results",
            str(reference_dir),
            "--metadata",
            str(toy_metadata_path),
            "--outdir",
            str(outdir),
            "--contrast-name",
            "toy",
        ],
    )

    assert result.returncode == 0, result.stderr
    report = json.loads((outdir / "validation" / "airway_pydeseq2_validation.json").read_text())
    assert report["passed"] is True


def test_06_make_deg_outputs_cli_writes_tables(tmp_path, toy_annotation_path, toy_config_file):
    pydeseq2_dir = tmp_path / "pydeseq2"
    pydeseq2_dir.mkdir()
    (pydeseq2_dir / "results_treated_vs_control.tsv").write_text(
        "gene_id\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\n"
        "ENSG000001\t100\t2\t0.2\t10\t0.001\t0.01\n"
        "ENSG000002\t50\t-2\t0.2\t-10\t0.001\t0.01\n"
    )
    outdir = tmp_path / "results"

    result = run_script(
        "06_make_deg_outputs.py",
        [
            "--pydeseq2-results-dir",
            str(pydeseq2_dir),
            "--annotation",
            str(toy_annotation_path),
            "--config",
            str(toy_config_file),
            "--outdir",
            str(outdir),
        ],
    )

    assert result.returncode == 0, result.stderr
    assert (outdir / "deg" / "treated_vs_control_up.tsv").exists()
    assert (outdir / "deg" / "treated_vs_control_down.tsv").exists()


def test_07_run_go_enrichment_cli_writes_summary(tmp_path, toy_annotation_path, toy_go_map_path, toy_config_file):
    deg_dir = tmp_path / "deg"
    deg_dir.mkdir()
    pd.DataFrame({"gene_id": ["ENSG000001", "ENSG000002"]}).to_csv(
        deg_dir / "treated_vs_control_all_genes.tsv", sep="\t", index=False
    )
    pd.DataFrame({"gene_id": ["ENSG000001"]}).to_csv(
        deg_dir / "treated_vs_control_up.tsv", sep="\t", index=False
    )
    pd.DataFrame({"gene_id": ["ENSG000002"]}).to_csv(
        deg_dir / "treated_vs_control_down.tsv", sep="\t", index=False
    )
    outdir = tmp_path / "results"

    result = run_script(
        "07_run_go_enrichment.py",
        [
            "--deg-dir",
            str(deg_dir),
            "--annotation",
            str(toy_annotation_path),
            "--gene2go",
            str(toy_go_map_path),
            "--config",
            str(toy_config_file),
            "--outdir",
            str(outdir),
        ],
    )

    assert result.returncode == 0, result.stderr
    assert (outdir / "go" / "go_enrichment_summary.json").exists()


def test_08_make_report_cli_writes_html(tmp_path, toy_config_file):
    project_dir = tmp_path / "project"
    deg_dir = project_dir / "deg"
    deg_dir.mkdir(parents=True)
    (deg_dir / "treated_vs_control_summary.json").write_text(
        json.dumps({"genes_up": 1, "genes_down": 1, "total_significant": 2})
    )
    outdir = tmp_path / "results"

    result = run_script(
        "08_make_report.py",
        ["--project-dir", str(project_dir), "--config", str(toy_config_file), "--outdir", str(outdir)],
    )

    assert result.returncode == 0, result.stderr
    assert (outdir / "report" / "report.html").exists()


def test_run_pipeline_dry_run_lists_commands(
    tmp_path, toy_counts_path, toy_metadata_path, toy_annotation_path, toy_go_map_path, toy_config_file
):
    result = run_script(
        "run_pipeline.py",
        [
            "--config",
            str(toy_config_file),
            "--counts",
            str(toy_counts_path),
            "--metadata",
            str(toy_metadata_path),
            "--annotation",
            str(toy_annotation_path),
            "--gene2go",
            str(toy_go_map_path),
            "--outdir",
            str(tmp_path / "results"),
            "--dry-run",
        ],
    )

    assert result.returncode == 0, result.stderr
    assert "04_run_pydeseq2.py" in result.stdout
    assert "08_make_report.py" in result.stdout
