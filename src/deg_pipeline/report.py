"""HTML report generation utilities for the DEG pipeline."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

from jinja2 import Template

logger = logging.getLogger(__name__)

# HTML report template
REPORT_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{{ project_name }} - DEG Pipeline Report</title>
    <style>
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif;
            line-height: 1.6;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f5f5f5;
        }
        .container {
            background-color: white;
            padding: 30px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        h1, h2, h3 {
            color: #333;
        }
        h1 {
            border-bottom: 2px solid #007bff;
            padding-bottom: 10px;
        }
        h2 {
            color: #007bff;
            margin-top: 30px;
        }
        table {
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }
        th, td {
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }
        th {
            background-color: #007bff;
            color: white;
        }
        tr:hover {
            background-color: #f5f5f5;
        }
        .metric-card {
            display: inline-block;
            background-color: #f8f9fa;
            padding: 15px 25px;
            margin: 10px;
            border-radius: 5px;
            border-left: 4px solid #007bff;
        }
        .metric-value {
            font-size: 24px;
            font-weight: bold;
            color: #007bff;
        }
        .metric-label {
            font-size: 14px;
            color: #666;
        }
        .status-pass {
            color: #28a745;
            font-weight: bold;
        }
        .status-fail {
            color: #dc3545;
            font-weight: bold;
        }
        .plot-container {
            margin: 20px 0;
            text-align: center;
        }
        .plot-container img {
            max-width: 100%;
            height: auto;
            border: 1px solid #ddd;
            border-radius: 5px;
        }
        .section {
            margin-bottom: 40px;
        }
        code {
            background-color: #f4f4f4;
            padding: 2px 6px;
            border-radius: 3px;
            font-family: 'Courier New', monospace;
        }
    </style>
</head>
<body>
    <div class="container">
        <h1>{{ project_name }} - DEG Pipeline Report</h1>
        
        <div class="section">
            <h2>Project Summary</h2>
            <div class="metric-card">
                <div class="metric-value">{{ n_samples }}</div>
                <div class="metric-label">Samples</div>
            </div>
            <div class="metric-card">
                <div class="metric-value">{{ n_genes_filtered }}</div>
                <div class="metric-label">Genes After Filtering</div>
            </div>
            <div class="metric-card">
                <div class="metric-value">{{ n_contrasts }}</div>
                <div class="metric-label">Contrasts</div>
            </div>
        </div>

        {% if validation %}
        <div class="section">
            <h2>Input Validation</h2>
            <p>Status: <span class="{% if validation.is_valid %}status-pass{% else %}status-fail{% endif %}">
                {% if validation.is_valid %}PASSED{% else %}FAILED{% endif %}
            </span></p>
            {% if validation.warnings %}
            <h3>Warnings</h3>
            <ul>
            {% for warning in validation.warnings %}
                <li>{{ warning }}</li>
            {% endfor %}
            </ul>
            {% endif %}
        </div>
        {% endif %}

        {% if qc_summary %}
        <div class="section">
            <h2>Quality Control</h2>
            <table>
                <tr><th>Metric</th><th>Value</th></tr>
                <tr><td>Total Counts</td><td>{{ qc_summary.total_counts | default('N/A') }}</td></tr>
                <tr><td>Median Library Size</td><td>{{ qc_summary.median_library_size | default('N/A') }}</td></tr>
                <tr><td>Median Detected Genes</td><td>{{ qc_summary.median_detected_genes | default('N/A') }}</td></tr>
            </table>
        </div>
        {% endif %}

        {% if filtering_summary %}
        <div class="section">
            <h2>Gene Filtering</h2>
            <table>
                <tr><th>Metric</th><th>Value</th></tr>
                <tr><td>Genes Before Filtering</td><td>{{ filtering_summary.n_genes_before_filtering | default('N/A') }}</td></tr>
                <tr><td>Genes After Filtering</td><td>{{ filtering_summary.n_genes_after_filtering | default('N/A') }}</td></tr>
                <tr><td>Genes Removed</td><td>{{ filtering_summary.n_genes_removed | default('N/A') }}</td></tr>
                <tr><td>Fraction Removed</td><td>{{ "%.1f%%" | format(filtering_summary.fraction_removed * 100 if filtering_summary.fraction_removed else 0) }}</td></tr>
            </table>
        </div>
        {% endif %}

        {% if deg_summaries %}
        <div class="section">
            <h2>Differential Expression Results</h2>
            {% for contrast_name, summary in deg_summaries.items() %}
            <h3>{{ contrast_name }}</h3>
            <table>
                <tr><th>Category</th><th>Count</th></tr>
                <tr><td>Up-regulated</td><td>{{ summary.genes_up | default(0) }}</td></tr>
                <tr><td>Down-regulated</td><td>{{ summary.genes_down | default(0) }}</td></tr>
                <tr><td>Not Significant</td><td>{{ summary.genes_not_significant | default(0) }}</td></tr>
                <tr><td>Total Significant</td><td>{{ summary.total_significant | default(0) }}</td></tr>
            </table>
            {% endfor %}
        </div>
        {% endif %}

        {% if go_summaries %}
        <div class="section">
            <h2>GO Enrichment</h2>
            {% for contrast_name, summary in go_summaries.items() %}
            <h3>{{ contrast_name }}</h3>
            <p>Gene Universe Size: {{ summary.universe_size | default('N/A') }}</p>
            {% endfor %}
        </div>
        {% endif %}

        <div class="section">
            <h2>Configuration</h2>
            <table>
                <tr><th>Parameter</th><th>Value</th></tr>
                <tr><td>Species</td><td>{{ config.species | default('N/A') }}</td></tr>
                <tr><td>Gene ID Type</td><td>{{ config.gene_id_type | default('N/A') }}</td></tr>
                <tr><td>Condition Column</td><td>{{ config.condition_column | default('N/A') }}</td></tr>
                <tr><td>Reference Level</td><td>{{ config.reference_level | default('N/A') }}</td></tr>
                <tr><td>padj Threshold</td><td>{{ config.deg.padj_threshold if config.deg else 'N/A' }}</td></tr>
                <tr><td>log2FC Threshold</td><td>{{ config.deg.log2fc_threshold if config.deg else 'N/A' }}</td></tr>
            </table>
        </div>

        <footer style="margin-top: 40px; padding-top: 20px; border-top: 1px solid #ddd; color: #666; font-size: 12px;">
            Generated by DEG Pipeline v0.1.0
        </footer>
    </div>
</body>
</html>
"""


def generate_report(
    project_dir: Path,
    config: Any,
) -> tuple[str, dict[str, Path]]:
    """Generate an HTML report for the pipeline results.

    Args:
        project_dir: Directory containing all pipeline results.
        config: Pipeline configuration object.

    Returns:
        Tuple of (HTML content, dict of asset paths).
    """
    logger.info(f"Generating report for project: {project_dir}")

    # Load validation results
    validation = None
    validation_path = project_dir / "validation" / "input_validation.json"
    if validation_path.exists():
        with open(validation_path) as f:
            validation = json.load(f)

    # Load QC summary
    qc_summary = None
    qc_summary_path = project_dir / "qc" / "qc_summary.json"
    if qc_summary_path.exists():
        with open(qc_summary_path) as f:
            qc_summary = json.load(f)

    # Load filtering summary
    filtering_summary = None
    filtering_path = project_dir / "filtering" / "gene_filtering_manifest.json"
    if filtering_path.exists():
        with open(filtering_path) as f:
            filtering_summary = json.load(f)

    # Load DEG summaries
    deg_summaries = {}
    deg_dir = project_dir / "deg"
    if deg_dir.exists():
        for summary_file in deg_dir.glob("*_summary.json"):
            contrast_name = summary_file.stem.replace("_summary", "")
            with open(summary_file) as f:
                deg_summaries[contrast_name] = json.load(f)

    # Load GO summaries
    go_summaries = {}
    go_dir = project_dir / "go"
    if go_dir.exists():
        for summary_file in go_dir.glob("*_summary.json"):
            contrast_name = summary_file.stem.replace("_summary", "")
            with open(summary_file) as f:
                go_summaries[contrast_name] = json.load(f)

    # Collect assets (plots)
    assets = {}
    plots_dir = project_dir / "qc" / "plots"
    if plots_dir.exists():
        for plot_file in plots_dir.glob("*.png"):
            assets[plot_file.stem] = plot_file

    deg_plots_dir = project_dir / "deg" / "plots"
    if deg_plots_dir.exists():
        for plot_file in deg_plots_dir.glob("*.png"):
            assets[plot_file.stem] = plot_file

    # Render template
    template = Template(REPORT_TEMPLATE)
    html_content = template.render(
        project_name=config.project_name,
        n_samples=qc_summary.get("n_samples", "N/A") if qc_summary else "N/A",
        n_genes_filtered=filtering_summary.get("n_genes_after_filtering", "N/A") if filtering_summary else "N/A",
        n_contrasts=len(deg_summaries),
        validation=validation,
        qc_summary=qc_summary,
        filtering_summary=filtering_summary,
        deg_summaries=deg_summaries,
        go_summaries=go_summaries,
        config=config.model_dump() if hasattr(config, "model_dump") else config,
    )

    logger.info(f"Generated report with {len(assets)} assets")
    return html_content, assets