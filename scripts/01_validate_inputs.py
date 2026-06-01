#!/usr/bin/env python3
"""Validate input files for the DEG pipeline.

This script validates count matrices, metadata, annotation, and configuration
before expensive analysis. It produces JSON and text validation reports.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# Add src to path for development
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from deg_pipeline.config import load_and_validate_config
from deg_pipeline.io import load_annotation, load_counts, load_metadata, save_dataframe
from deg_pipeline.logging_utils import configure_logging, write_manifest
from deg_pipeline.validation import reorder_samples, validate_inputs


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
        help="Path to gene count matrix TSV file",
    )
    parser.add_argument(
        "--metadata",
        required=True,
        type=Path,
        help="Path to sample metadata TSV file",
    )
    parser.add_argument(
        "--config",
        required=True,
        type=Path,
        help="Path to YAML configuration file",
    )
    parser.add_argument(
        "--annotation",
        type=Path,
        help="Path to gene annotation TSV file (optional)",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Output directory for validation reports",
    )
    parser.add_argument(
        "--allow-extra-metadata",
        action="store_true",
        help="Allow metadata samples not present in counts",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output files",
    )
    return parser


def run(args: argparse.Namespace) -> int:
    """Run the validation pipeline stage.

    Args:
        args: Parsed command line arguments.

    Returns:
        Exit code (0 for success, 1 for validation errors).
    """
    logger = configure_logging(args.outdir, "01_validate_inputs")
    logger.info("Starting input validation")

    # Create output directories
    validation_dir = args.outdir / "validation"
    intermediate_dir = args.outdir / "intermediate"
    validation_dir.mkdir(parents=True, exist_ok=True)
    intermediate_dir.mkdir(parents=True, exist_ok=True)

    # Load inputs
    logger.info(f"Loading counts from {args.counts}")
    try:
        counts_df = load_counts(args.counts)
    except (FileNotFoundError, ValueError) as e:
        logger.error(f"Failed to load counts: {e}")
        return 1

    logger.info(f"Loading metadata from {args.metadata}")
    try:
        metadata_df = load_metadata(args.metadata)
    except (FileNotFoundError, ValueError) as e:
        logger.error(f"Failed to load metadata: {e}")
        return 1

    logger.info(f"Loading config from {args.config}")
    try:
        config = load_and_validate_config(args.config)
        config_dict = config.model_dump()
    except Exception as e:
        logger.error(f"Failed to load config: {e}")
        return 1

    # Load optional annotation
    annotation_df = None
    if args.annotation:
        logger.info(f"Loading annotation from {args.annotation}")
        try:
            annotation_df = load_annotation(args.annotation)
        except (FileNotFoundError, ValueError) as e:
            logger.warning(f"Failed to load annotation: {e}")

    # Run validation
    logger.info("Running input validation")
    result = validate_inputs(
        counts_df,
        metadata_df,
        config_dict,
        annotation_df,
        allow_extra_metadata=args.allow_extra_metadata,
    )

    # Write validation reports
    json_report_path = validation_dir / "input_validation.json"
    txt_report_path = validation_dir / "input_validation.txt"

    report_data = result.to_dict()
    report_data["input_files"] = {
        "counts": str(args.counts),
        "metadata": str(args.metadata),
        "config": str(args.config),
        "annotation": str(args.annotation) if args.annotation else None,
    }

    with open(json_report_path, "w") as f:
        json.dump(report_data, f, indent=2, default=str)
    logger.info(f"Written JSON validation report to {json_report_path}")

    # Write human-readable text report
    with open(txt_report_path, "w") as f:
        f.write("=" * 60 + "\n")
        f.write("INPUT VALIDATION REPORT\n")
        f.write("=" * 60 + "\n\n")

        f.write(f"Status: {'PASSED' if result.is_valid else 'FAILED'}\n\n")

        if result.errors:
            f.write("ERRORS:\n")
            for error in result.errors:
                f.write(f"  - {error}\n")
            f.write("\n")

        if result.warnings:
            f.write("WARNINGS:\n")
            for warning in result.warnings:
                f.write(f"  - {warning}\n")
            f.write("\n")

        f.write("METRICS:\n")
        for key, value in result.metrics.items():
            f.write(f"  {key}: {value}\n")

    logger.info(f"Written text validation report to {txt_report_path}")

    # If validation failed, return error
    if not result.is_valid:
        logger.error(f"Validation failed with {len(result.errors)} errors")
        return 1

    # Reorder samples and write intermediate files
    logger.info("Reordering samples and writing intermediate files")
    counts_ordered, metadata_ordered = reorder_samples(counts_df, metadata_df)

    counts_out_path = intermediate_dir / "counts_ordered.tsv"
    metadata_out_path = intermediate_dir / "metadata_ordered.tsv"

    save_dataframe(counts_ordered, counts_out_path, index_label="gene_id")
    save_dataframe(metadata_ordered, metadata_out_path, index_label="sample_id")

    logger.info(f"Written ordered counts to {counts_out_path}")
    logger.info(f"Written ordered metadata to {metadata_out_path}")

    # Write manifest
    output_files = [
        json_report_path,
        txt_report_path,
        counts_out_path,
        metadata_out_path,
    ]
    manifest_path = write_manifest(
        validation_dir,
        "01_validate_inputs",
        config=config_dict,
        input_files=[args.counts, args.metadata, args.config] + ([args.annotation] if args.annotation else []),
        output_files=output_files,
        extra_info={"validation_passed": result.is_valid},
    )
    logger.info(f"Written manifest to {manifest_path}")

    logger.info("Input validation completed successfully")
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