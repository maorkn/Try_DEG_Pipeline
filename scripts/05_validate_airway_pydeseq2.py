#!/usr/bin/env python3
"""Validate PyDESeq2 results against R/Bioconductor DESeq2 reference.

This script compares PyDESeq2 outputs to pinned R DESeq2 reference outputs
for the airway dataset.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Add src to path for development
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from deg_pipeline.logging_utils import configure_logging, write_manifest
from deg_pipeline.validation import (
    align_results,
    compute_airway_validation_metrics,
    evaluate_validation,
    load_deseq2_results,
    load_pydeseq2_results,
    write_validation_report,
)


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--pipeline-results",
        required=True,
        type=Path,
        help="Directory containing PyDESeq2 pipeline results",
    )
    parser.add_argument(
        "--reference-results",
        required=True,
        type=Path,
        help="Directory containing R/Bioconductor DESeq2 reference results",
    )
    parser.add_argument(
        "--metadata",
        required=True,
        type=Path,
        help="Path to sample metadata TSV file",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Output directory for validation results",
    )
    parser.add_argument(
        "--contrast-name",
        required=True,
        help="Name of contrast to validate",
    )
    parser.add_argument(
        "--lfc-correlation-threshold",
        type=float,
        default=0.95,
        help="Minimum correlation threshold for log2 fold changes",
    )
    parser.add_argument(
        "--stat-correlation-threshold",
        type=float,
        default=0.95,
        help="Minimum correlation threshold for test statistics",
    )
    parser.add_argument(
        "--significant-count-tolerance",
        type=float,
        default=0.15,
        help="Tolerance for significant gene count comparison",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output files",
    )
    return parser


def run(args: argparse.Namespace) -> int:
    """Run the airway validation pipeline stage."""
    logger = configure_logging(args.outdir, "05_validate_airway_pydeseq2")
    logger.info("Starting airway PyDESeq2 validation")

    # Create output directory
    validation_dir = args.outdir / "validation"
    validation_dir.mkdir(parents=True, exist_ok=True)

    # Load pipeline results
    logger.info(f"Loading pipeline results from {args.pipeline_results}")
    pipeline_results = load_pydeseq2_results(args.pipeline_results, args.contrast_name)

    # Load reference results
    logger.info(f"Loading reference results from {args.reference_results}")
    reference_results = load_deseq2_results(args.reference_results, args.contrast_name)

    # Align results
    logger.info("Aligning results by gene ID")
    aligned_df = align_results(pipeline_results, reference_results)

    # Compute validation metrics
    logger.info("Computing validation metrics")
    metrics = compute_airway_validation_metrics(
        aligned_df,
        lfc_threshold=args.lfc_correlation_threshold,
        stat_threshold=args.stat_correlation_threshold,
    )

    # Evaluate validation
    logger.info("Evaluating validation results")
    validation_report = evaluate_validation(metrics, args)

    # Write validation report
    json_report_path = validation_dir / "airway_pydeseq2_validation.json"
    txt_report_path = validation_dir / "airway_pydeseq2_validation.txt"

    write_validation_report(validation_report, json_report_path, txt_report_path)

    # Write manifest
    output_files = [json_report_path, txt_report_path]
    manifest_path = write_manifest(
        validation_dir,
        "05_validate_airway_pydeseq2",
        input_files=[args.pipeline_results, args.reference_results, args.metadata],
        output_files=output_files,
        extra_info={
            "contrast_name": args.contrast_name,
            "lfc_correlation_threshold": args.lfc_correlation_threshold,
            "stat_correlation_threshold": args.stat_correlation_threshold,
            "significant_count_tolerance": args.significant_count_tolerance,
            "validation_passed": validation_report["passed"],
        },
    )
    logger.info(f"Written manifest to {manifest_path}")

    # Return exit code based on validation result
    return 0 if validation_report["passed"] else 1


def main(argv: list[str] | None = None) -> int:
    """Main entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
