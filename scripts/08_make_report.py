#!/usr/bin/env python3
"""Generate HTML report bundle for the DEG pipeline.

This script generates a comprehensive HTML report summarizing all
pipeline outputs including validation, QC, DEG analysis, and GO enrichment.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Add src to path for development
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from deg_pipeline.config import load_and_validate_config
from deg_pipeline.io import save_json
from deg_pipeline.logging_utils import configure_logging, write_manifest
from deg_pipeline.report import generate_report


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--project-dir",
        required=True,
        type=Path,
        help="Project directory containing all results",
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
        help="Output directory for report",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output files",
    )
    return parser


def run(args: argparse.Namespace) -> int:
    """Run the report generation pipeline stage."""
    logger = configure_logging(args.outdir, "08_make_report")
    logger.info("Starting report generation")

    # Create output directory
    report_dir = args.outdir / "report"
    report_dir.mkdir(parents=True, exist_ok=True)

    # Load config
    config = load_and_validate_config(args.config)
    config_dict = config.model_dump()

    # Generate report
    logger.info("Generating HTML report")
    report_html, report_assets = generate_report(args.project_dir, config)

    # Write report
    report_path = report_dir / "report.html"
    report_path.write_text(report_html)

    # Copy assets
    assets_dir = report_dir / "assets"
    assets_dir.mkdir(exist_ok=True)
    for asset_name, asset_path in report_assets.items():
        asset_dest = assets_dir / Path(asset_path).name
        if asset_path.exists():
            import shutil
            shutil.copy2(asset_path, asset_dest)

    # Write manifest
    output_files = [report_path] + list(report_assets.values())
    manifest_path = write_manifest(
        report_dir,
        "08_make_report",
        config=config_dict,
        input_files=[args.project_dir, args.config],
        output_files=output_files,
        extra_info={"report_generated": True},
    )
    logger.info(f"Written manifest to {manifest_path}")

    logger.info("Report generation completed successfully")
    return 0


def main(argv: list[str] | None = None) -> int:
    """Main entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())