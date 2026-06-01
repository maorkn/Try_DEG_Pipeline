#!/usr/bin/env python3
"""Fetch and prepare the airway dataset for validation.

This script downloads or extracts the Bioconductor airway dataset
and generates R/Bioconductor DESeq2 reference outputs for validation.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

# Add src to path for development
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from deg_pipeline.logging_utils import configure_logging, write_manifest


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Output directory for external airway files",
    )
    parser.add_argument(
        "--reduced-fixture-outdir",
        required=True,
        type=Path,
        help="Output directory for reduced airway fixture files",
    )
    parser.add_argument(
        "--n-top-variable-genes",
        type=int,
        default=200,
        help="Number of top variable genes to include in the reduced fixture",
    )
    parser.add_argument(
        "--make-reference",
        action="store_true",
        help="Generate R/Bioconductor DESeq2 reference output",
    )
    parser.add_argument("--force", action="store_true", help="Overwrite existing outputs")
    return parser


def check_rscript_available() -> bool:
    """Return whether Rscript is available on PATH."""
    return shutil.which("Rscript") is not None


def render_airway_export_r_script(outdir: Path, make_reference: bool = True) -> str:
    """Render R script for airway dataset export.

    Args:
        outdir: Output directory for exported files.
        make_reference: Whether to generate DESeq2 reference outputs.

    Returns:
        R script as string.
    """
    script = f'''
# Airway dataset export script
# This script exports the airway dataset for pipeline validation

library(airway)
library(SummarizedExperiment)
library(DESeq2)

# Load airway dataset
data(airway)

# Create output directory
dir.create("{outdir}", recursive = TRUE, showWarnings = FALSE)

# Extract counts and metadata
counts <- assay(airway)
metadata <- colData(airway)
metadata <- as.data.frame(metadata)
metadata$sample_id <- rownames(metadata)
metadata <- metadata[, c("sample_id", setdiff(colnames(metadata), "sample_id"))]

# Write counts
write.table(
    counts,
    file = file.path("{outdir}", "counts.tsv"),
    sep = "\\t",
    quote = FALSE,
    row.names = TRUE,
    col.names = TRUE
)

# Write metadata
write.table(
    metadata,
    file = file.path("{outdir}", "metadata.tsv"),
    sep = "\\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)

# Generate DESeq2 reference if requested
if ({str(make_reference).upper()}) {{
    library(DESeq2)

    dds <- DESeqDataSetFromMatrix(
        countData = counts,
        colData = metadata,
        design = ~ dex
    )

    dds <- DESeq(dds)
    res <- results(dds, contrast = c("dex", "trt", "untrt"))

    write.table(
        res,
        file = file.path("{outdir}", "reference_deseq2_dex_trt_vs_untrt.tsv"),
        sep = "\\t",
        quote = FALSE,
        row.names = TRUE,
        col.names = TRUE
    )
}}

# Write session info
writeLines(
    capture.output(sessionInfo()),
    file = file.path("{outdir}", "session_info.txt")
)
'''
    return script


def run(args: argparse.Namespace) -> int:
    """Run the airway fetch pipeline stage."""
    logger = configure_logging(args.outdir, "00_fetch_airway")
    logger.info("Starting airway dataset fetch")

    # Create output directories
    outdir = args.outdir / "airway"
    outdir.mkdir(parents=True, exist_ok=True)

    fixture_outdir = args.reduced_fixture_outdir
    fixture_outdir.mkdir(parents=True, exist_ok=True)

    # Check Rscript availability
    rscript_available = check_rscript_available()
    if not rscript_available:
        logger.error("Rscript not available; airway export requires R/Bioconductor")
        return 1

    # Render R script
    r_script = render_airway_export_r_script(outdir, make_reference=args.make_reference)

    # Write R script
    r_script_path = outdir / "fetch_airway.R"
    r_script_path.write_text(r_script)

    # Execute R script if Rscript is available
    if rscript_available:
        result = subprocess.run(
            ["Rscript", str(r_script_path)],
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            logger.error(f"R script failed: {result.stderr}")
            return 1

    # Read exported files
    import pandas as pd

    counts_path = outdir / "counts.tsv"
    metadata_path = outdir / "metadata.tsv"

    if counts_path.exists():
        counts_df = pd.read_csv(counts_path, sep="\t", index_col=0)
    else:
        logger.warning(f"Counts file not found: {counts_path}")
        counts_df = None

    if metadata_path.exists():
        metadata_df = pd.read_csv(metadata_path, sep="\t").set_index("sample_id")
    else:
        logger.warning(f"Metadata file not found: {metadata_path}")
        metadata_df = None

    # Create reduced fixture
    logger.info("Creating reduced fixture")
    if counts_df is not None and metadata_df is not None:
        # Select top variable genes plus some low-count genes
        gene_var = counts_df.var(axis=1)
        top_genes = gene_var.nlargest(args.n_top_variable_genes).index
        low_genes = counts_df.index[counts_df.max(axis=1) < 10][:5]
        selected_genes = top_genes.append(low_genes).unique()

        fixture_counts = counts_df.loc[selected_genes]

        # Write fixture
        fixture_counts_path = fixture_outdir / "airway_small_counts.tsv"
        fixture_counts.to_csv(fixture_counts_path, sep="\t")

        fixture_metadata = metadata_df.loc[fixture_counts.columns]
        fixture_metadata_path = fixture_outdir / "airway_small_metadata.tsv"
        fixture_metadata.to_csv(fixture_metadata_path, sep="\t")

    # Write fetch manifest
    manifest = {
        "outdir": str(outdir),
        "make_reference": args.make_reference,
        "rscript_available": rscript_available,
        "counts_file": str(counts_path) if counts_path.exists() else None,
        "metadata_file": str(metadata_path) if metadata_path.exists() else None,
        "fixture_counts_file": str(fixture_counts_path) if counts_df is not None else None,
        "fixture_metadata_file": str(fixture_metadata_path) if counts_df is not None else None,
    }

    manifest_path = write_manifest(
        outdir,
        "00_fetch_airway",
        extra_info=manifest,
    )

    logger.info(f"Written fetch manifest to {manifest_path}")
    logger.info("Airway dataset fetch completed successfully")
    return 0


def main(argv: list[str] | None = None) -> int:
    """Main entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
