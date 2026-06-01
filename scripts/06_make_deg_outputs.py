#!/usr/bin/env python3
"""Generate DEG output tables and plots from PyDESeq2 results.

This script converts raw PyDESeq2 results into analysis deliverables,
joins annotation, classifies DEGs, and generates plots.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Add src to path for development
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from deg_pipeline.config import load_and_validate_config
from deg_pipeline.deg import (
    classify_degs,
    compute_deg_summary,
    join_annotation,
    make_ranked_file,
)
from deg_pipeline.io import load_annotation, load_counts, load_metadata, save_dataframe, save_json
from deg_pipeline.logging_utils import configure_logging, write_manifest
from deg_pipeline.plotting import plot_ma, plot_volcano


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--pydeseq2-results-dir",
        required=True,
        type=Path,
        help="Directory containing PyDESeq2 results",
    )
    parser.add_argument(
        "--normalized-counts",
        type=Path,
        help="Path to normalized counts TSV file",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        help="Path to metadata TSV file",
    )
    parser.add_argument(
        "--annotation",
        type=Path,
        help="Path to gene annotation TSV file (optional)",
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
        help="Output directory for DEG results",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output files",
    )
    return parser


def run(args: argparse.Namespace) -> int:
    """Run the DEG outputs pipeline stage."""
    logger = configure_logging(args.outdir, "06_make_deg_outputs")
    logger.info("Starting DEG outputs generation")

    # Create output directories
    deg_dir = args.outdir / "deg"
    plots_dir = deg_dir / "plots"
    deg_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    # Load config
    config = load_and_validate_config(args.config)
    config_dict = config.model_dump()

    # Load annotation if provided
    annotation_df = None
    if args.annotation:
        logger.info(f"Loading annotation from {args.annotation}")
        annotation_df = load_annotation(args.annotation)

    # Process each contrast
    output_files = []
    all_summaries = {}

    for contrast in config.contrasts:
        logger.info(f"Processing contrast: {contrast.name}")

        # Load results
        results_path = args.pydeseq2_results_dir / f"results_{contrast.name}.tsv"
        if not results_path.exists():
            logger.warning(f"Results file not found: {results_path}, skipping")
            continue

        import pandas as pd
        results_df = pd.read_csv(results_path, sep="\t")

        # Join annotation
        results_df = join_annotation(results_df, annotation_df)

        # Classify DEGs
        results_df = classify_degs(
            results_df,
            padj_threshold=config.deg.padj_threshold,
            log2fc_threshold=config.deg.log2fc_threshold,
        )

        # Write output tables
        all_genes_path = deg_dir / f"{contrast.name}_all_genes.tsv"
        save_dataframe(results_df, all_genes_path, include_index=False)
        output_files.append(all_genes_path)

        sig_path = deg_dir / f"{contrast.name}_significant.tsv"
        sig_df = results_df[results_df["deg_class"].isin(["up", "down", "significant_low_lfc"])]
        save_dataframe(sig_df, sig_path, include_index=False)
        output_files.append(sig_path)

        up_path = deg_dir / f"{contrast.name}_up.tsv"
        up_df = results_df[results_df["deg_class"] == "up"]
        save_dataframe(up_df, up_path, include_index=False)
        output_files.append(up_path)

        down_path = deg_dir / f"{contrast.name}_down.tsv"
        down_df = results_df[results_df["deg_class"] == "down"]
        save_dataframe(down_df, down_path, include_index=False)
        output_files.append(down_path)

        # Write ranked file
        ranked_df = make_ranked_file(results_df)
        ranked_path = deg_dir / f"{contrast.name}_ranked_for_gsea.rnk"
        ranked_df.to_csv(ranked_path, sep="\t", index=False, header=False)
        output_files.append(ranked_path)

        # Compute and write summary
        summary = compute_deg_summary(results_df, contrast.name)
        all_summaries[contrast.name] = summary

        summary_path = deg_dir / f"{contrast.name}_summary.json"
        save_json(summary, summary_path)
        output_files.append(summary_path)

        # Generate plots
        volcano_path = plots_dir / f"{contrast.name}_volcano.png"
        plot_volcano(
            results_df,
            volcano_path,
            padj_threshold=config.deg.padj_threshold,
            log2fc_threshold=config.deg.log2fc_threshold,
            title=f"Volcano: {contrast.name}",
        )
        output_files.append(volcano_path)

        ma_path = plots_dir / f"{contrast.name}_ma_plot.png"
        plot_ma(results_df, ma_path, padj_threshold=config.deg.padj_threshold, title=f"MA: {contrast.name}")
        output_files.append(ma_path)

    # Write manifest
    manifest_path = write_manifest(
        deg_dir,
        "06_make_deg_outputs",
        config=config_dict,
        input_files=[args.pydeseq2_results_dir],
        output_files=output_files,
        extra_info={"summaries": all_summaries},
    )
    logger.info(f"Written manifest to {manifest_path}")

    logger.info("DEG outputs generation completed successfully")
    return 0


def main(argv: list[str] | None = None) -> int:
    """Main entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())