"""Unit tests for report generation."""

from __future__ import annotations

import json

from deg_pipeline.config import load_and_validate_config
from deg_pipeline.report import generate_report


def test_generate_report_includes_project_name(tmp_path, toy_config_file):
    project_dir = tmp_path / "project"
    deg_dir = project_dir / "deg"
    deg_dir.mkdir(parents=True)
    (deg_dir / "treated_vs_control_summary.json").write_text(
        json.dumps(
            {
                "genes_up": 1,
                "genes_down": 2,
                "genes_not_significant": 3,
                "total_significant": 3,
            }
        )
    )
    config = load_and_validate_config(toy_config_file)

    html, assets = generate_report(project_dir, config)

    assert config.project_name in html
    assert "Differential Expression Results" in html
    assert assets == {}
