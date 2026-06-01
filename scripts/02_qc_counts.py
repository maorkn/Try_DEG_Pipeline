#!/usr/bin/env python3
"""Generate QC metrics and plots from count matrices.

This script computes sample-level and gene-level QC metrics,
generates diagnostic plots, and identifies potential outliers.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# Add src to path for development
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from deg_pipeline.config import load_and_validate_config
from deg_pipeline.io import load_counts, load_metadata, save_dataframe, save_json
from deg_pipeline.logging_utils import configure_logging, write_manifest
from deg_pipeline.plotting import (
    plot_detected_genes,
    plot_library_sizes,
    plot_pca,
    plot_sample_correlation_heatmap,
)
from deg_pipeline.qc import (
    compute_gene_metrics,
    compute_library_metrics,
    compute_pca,
    compute_qc_summary,
    compute_sample_correlation,
)


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--counts",
        required=True,
        type=Path,
        help="Path to ordered count matrix TSV file (from validation stage)",
    )
    parser.add_argument(
        "--metadata",
        required=True,
        type=Path,
        help="Path to ordered metadata TSV file (from validation stage)",
    )
    parser.add_argument(
        "--config",
        required=True,
        type=Path,
        help="Path to YAML configuration file",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Output directory for QC results",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output files",
    )
    return parser


def run(args: argparse.Namespace) -> int:
    """Run the QC pipeline stage.

    Args:
        args: Parsed command line arguments.

    Returns:
        Exit code (0 for success).
    """
    logger = configure_logging(args.outdir, "02_qc_counts")
    logger.info("Starting count QC")

    # Create output directories
    qc_dir = args.outdir / "qc"
    plots_dir = qc_dir / "plots"
    qc_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    # Load inputs
    logger.info(f"Loading counts from {args.counts}")
    counts_df = load_counts(args.counts)

    logger.info(f"Loading metadata from {args.metadata}")
    metadata_df = load_metadata(args.metadata)

    logger.info(f"Loading config from {args.config}")
    config = load_and_validate_config(args.config)
    config_dict = config.model_dump()

    # Compute metrics
    logger.info("Computing library metrics")
    library_metrics = compute_library_metrics(counts_df)
    library_metrics_path = qc_dir / "library_sizes.tsv"
    save_dataframe(library_metrics, library_metrics_path, index_label="sample_id")

    logger.info("Computing gene metrics")
    gene_metrics = compute_gene_metrics(counts_df)
    gene_detection_path = qc_dir / "gene_detection.tsv"
    save_dataframe(gene_metrics, gene_detection_path, index_label="gene_id")

    logger.info("Computing sample correlation")
    correlation = compute_sample_correlation(counts_df)
    correlation_path = qc_dir / "sample_correlation.tsv"
    save_dataframe(correlation, correlation_path)

    logger.info("Computing PCA")
    condition_col = config.condition_column
    pca_df, variance_explained = compute_pca(
        counts_df,
        metadata_df=metadata_df,
    )
    pca_path = qc_dir / "pca_coordinates.tsv"
    save_dataframe(pca_df, pca_path, index_label="sample_id")

    # Generate plots
    logger.info("Generating plots")

    library_plot_path = plots_dir / "library_sizes.png"
    plot_library_sizes(library_metrics, library_plot_path)

    detected_genes_plot_path = plots_dir / "detected_genes.png"
    plot_detected_genes(gene_metrics, detected_genes_plot_path)

    correlation_heatmap_path = plots_dir / "sample_correlation_heatmap.png"
    plot_sample_correlation_heatmap(correlation, correlation_heatmap_path)

    pca_plot_path = plots_dir / "pca_raw_or_vst.png"
    plot_pca(pca_df, pca_plot_path, color_by=condition_col)

    # Compute and write QC summary
    logger.info("Computing QC summary")
    qc_summary = compute_qc_summary(counts_df, library_metrics, gene_metrics)
    qc_summary["variance_explained"] = variance_explained.tolist()
    qc_summary["condition_column"] = condition_col

    qc_summary_path = qc_dir / "qc_summary.json"
    save_json(qc_summary, qc_summary_path)

    # Write manifest
    output_files = [
        library_metrics_path,
        gene_detection_path,
        correlation_path,
        pca_path,
        library_plot_path,
        detected_genes_plot_path,
        correlation_heatmap_path,
        pca_plot_path,
        qc_summary_path,
    ]
    manifest_path = write_manifest(
        qc_dir,
        "02_qc_counts",
        config=config_dict,
        input_files=[args.counts, args.metadata, args.config],
        output_files=output_files,
        extra_info=qc_summary,
    )
    logger.info(f"Written manifest to {manifest_path}")

    logger.info("Count QC completed successfully")
    return 0


def main(argv: list[str] | None = None) -> int:
    """Main entry point.

    Args:
        argv: Command line arguments (uses sys.argv if None).

    Returns:
        Exit code.
    """
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())