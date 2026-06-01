#!/usr/bin/env python3
"""Run GO enrichment analysis on DEG gene sets.

This script performs over-representation analysis (ORA) using goatools
for up-regulated and down-regulated gene sets.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Add src to path for development
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from deg_pipeline.config import load_and_validate_config
from deg_pipeline.enrichment import (
    build_gene_universe,
    load_deg_gene_sets,
    load_gene2go_table,
    run_go_ora,
)
from deg_pipeline.io import load_annotation, save_dataframe, save_json
from deg_pipeline.logging_utils import configure_logging, write_manifest


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--deg-dir",
        required=True,
        type=Path,
        help="Directory containing DEG gene sets",
    )
    parser.add_argument(
        "--annotation",
        required=True,
        type=Path,
        help="Path to gene annotation TSV file",
    )
    parser.add_argument(
        "--go-obo",
        required=False,
        type=Path,
        help="Path to go-basic.obo file (recorded for provenance; TSV mapping drives v1 ORA)",
    )
    parser.add_argument(
        "--gene2go",
        required=True,
        type=Path,
        help="Path to gene2go mapping file",
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
        help="Output directory for GO enrichment results",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output files",
    )
    return parser


def run(args: argparse.Namespace) -> int:
    """Run the GO enrichment pipeline stage."""
    logger = configure_logging(args.outdir, "07_run_go_enrichment")
    logger.info("Starting GO enrichment analysis")

    # Create output directory
    go_dir = args.outdir / "go"
    go_dir.mkdir(parents=True, exist_ok=True)

    # Load config
    config = load_and_validate_config(args.config)
    config_dict = config.model_dump()

    # Load annotation
    logger.info(f"Loading annotation from {args.annotation}")
    annotation_df = load_annotation(args.annotation)

    gene2go_df = load_gene2go_table(args.gene2go)

    output_files = []
    summary = {"contrasts": {}, "gene2go": str(args.gene2go)}

    for contrast in config.contrasts:
        logger.info(f"Loading DEG gene sets for {contrast.name}")
        gene_sets = load_deg_gene_sets(args.deg_dir, contrast.name)
        all_genes_path = args.deg_dir / f"{contrast.name}_all_genes.tsv"
        if all_genes_path.exists():
            import pandas as pd

            filtered_gene_ids = pd.read_csv(all_genes_path, sep="\t")["gene_id"].astype(str)
        else:
            filtered_gene_ids = annotation_df.index.astype(str)

        universe = build_gene_universe(set(filtered_gene_ids), gene2go_df)
        contrast_summary = {
            "gene_sets_analyzed": list(gene_sets.keys()),
            "universe_size": len(universe),
            "outputs": [],
        }

        for direction in ["up", "down"]:
            direction_genes = gene_sets.get(direction, set())
            for ontology in config.go.ontology:
                ontology_value = ontology.value if hasattr(ontology, "value") else str(ontology)
                logger.info(
                    "Running GO ORA for %s %s genes in %s",
                    contrast.name,
                    direction,
                    ontology_value,
                )
                results = run_go_ora(
                    direction_genes,
                    universe,
                    ontology_value,
                    args.go_obo,
                    gene2go_df,
                    padj_method=config.go.padj_method,
                    min_genes_per_term=config.go.min_genes_per_term,
                    max_genes_per_term=config.go.max_genes_per_term,
                )
                results_path = (
                    go_dir / f"{contrast.name}_{direction}_go_{ontology_value.lower()}.tsv"
                )
                save_dataframe(results, results_path, include_index=False)
                output_files.append(results_path)
                contrast_summary["outputs"].append(str(results_path))

        summary["contrasts"][contrast.name] = contrast_summary

    summary_path = go_dir / "go_enrichment_summary.json"
    save_json(summary, summary_path)
    output_files.append(summary_path)

    # Write manifest
    manifest_path = write_manifest(
        go_dir,
        "07_run_go_enrichment",
        config=config_dict,
        input_files=[
            path
            for path in [args.deg_dir, args.annotation, args.go_obo, args.gene2go]
            if path is not None
        ],
        output_files=output_files,
        extra_info=summary,
    )
    logger.info(f"Written manifest to {manifest_path}")

    logger.info("GO enrichment analysis completed successfully")
    return 0


def main(argv: list[str] | None = None) -> int:
    """Main entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
