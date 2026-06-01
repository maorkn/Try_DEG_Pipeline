#!/usr/bin/env python3
"""Filter low-count genes and prepare data for PyDESeq2.

This script applies reproducible low-count gene filtering before PyDESeq2.
It preserves raw counts for PyDESeq2 size factor estimation.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Add src to path for development
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from deg_pipeline.config import load_and_validate_config
from deg_pipeline.filtering import (
    compute_filtering_summary,
    determine_min_samples,
    filter_low_count_genes,
)
from deg_pipeline.io import load_counts, load_metadata, save_dataframe, save_json
from deg_pipeline.logging_utils import configure_logging, write_manifest


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
        help="Output directory for filtering results",
    )
    parser.add_argument(
        "--min-count",
        type=int,
        help="Override minimum count threshold from config",
    )
    parser.add_argument(
        "--min-samples",
        type=int,
        help="Override minimum samples threshold from config",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output files",
    )
    return parser


def run(args: argparse.Namespace) -> int:
    """Run the filtering pipeline stage.

    Args:
        args: Parsed command line arguments.

    Returns:
        Exit code (0 for success, 1 for failure).
    """
    logger = configure_logging(args.outdir, "03_filter_normalize")
    logger.info("Starting gene filtering")

    # Create output directories
    filtering_dir = args.outdir / "filtering"
    intermediate_dir = args.outdir / "intermediate"
    filtering_dir.mkdir(parents=True, exist_ok=True)
    intermediate_dir.mkdir(parents=True, exist_ok=True)

    # Load inputs
    logger.info(f"Loading counts from {args.counts}")
    counts_df = load_counts(args.counts)

    logger.info(f"Loading metadata from {args.metadata}")
    metadata_df = load_metadata(args.metadata)

    logger.info(f"Loading config from {args.config}")
    config = load_and_validate_config(args.config)
    config_dict = config.model_dump()

    # Determine filtering thresholds
    min_count = args.min_count if args.min_count is not None else config.filtering.min_count
    min_samples = determine_min_samples(
        metadata_df,
        config.condition_column,
        args.min_samples if args.min_samples is not None else config.filtering.min_samples,
    )

    logger.info(f"Filtering with min_count={min_count}, min_samples={min_samples}")

    # Apply filtering
    filtered_counts, gene_summary = filter_low_count_genes(
        counts_df, min_count=min_count, min_samples=min_samples
    )

    # Check if no genes remain; small fixtures may retain fewer than 10 genes.
    if filtered_counts.shape[0] == 0:
        logger.error(
            "No genes remain after filtering. "
            "Consider relaxing filtering thresholds."
        )
        return 1
    if filtered_counts.shape[0] < 10:
        logger.warning(
            "Only %s genes remain after filtering; this is acceptable for smoke tests "
            "but too small for real differential expression.",
            filtered_counts.shape[0],
        )

    # Write filtered counts
    filtered_counts_path = intermediate_dir / "counts_filtered.tsv"
    save_dataframe(filtered_counts, filtered_counts_path, index_label="gene_id")

    # Write gene filtering summary
    gene_summary_path = filtering_dir / "gene_filtering_summary.tsv"
    save_dataframe(gene_summary, gene_summary_path, index_label="gene_id")

    # Compute and write filtering manifest
    filtering_summary = compute_filtering_summary(
        counts_df, filtered_counts, min_count, min_samples
    )
    filtering_manifest_path = filtering_dir / "gene_filtering_manifest.json"
    save_json(filtering_summary, filtering_manifest_path)

    # Write manifest
    output_files = [
        filtered_counts_path,
        gene_summary_path,
        filtering_manifest_path,
    ]
    manifest_path = write_manifest(
        filtering_dir,
        "03_filter_normalize",
        config=config_dict,
        input_files=[args.counts, args.metadata, args.config],
        output_files=output_files,
        extra_info=filtering_summary,
    )
    logger.info(f"Written manifest to {manifest_path}")

    logger.info("Gene filtering completed successfully")
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
