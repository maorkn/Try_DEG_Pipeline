#!/usr/bin/env python3
"""Run PyDESeq2 differential expression analysis for simple two-group contrasts."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from deg_pipeline.config import load_and_validate_config
from deg_pipeline.io import load_counts, load_metadata, save_dataframe
from deg_pipeline.logging_utils import configure_logging, write_manifest
from deg_pipeline.pydeseq2_adapter import (
    extract_pydeseq2_outputs,
    prepare_pydeseq2_inputs,
    run_pydeseq2_model,
    run_pydeseq2_stats,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--counts-filtered",
        required=True,
        type=Path,
        help="Path to filtered raw count matrix TSV",
    )
    parser.add_argument(
        "--metadata",
        required=True,
        type=Path,
        help="Path to ordered sample metadata TSV",
    )
    parser.add_argument("--config", required=True, type=Path, help="Path to YAML config")
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Project output directory; outputs are written under pydeseq2/",
    )
    parser.add_argument(
        "--contrast-name",
        help="Optional contrast name to run. Defaults to the first configured contrast.",
    )
    parser.add_argument("--force", action="store_true", help="Overwrite existing outputs")
    return parser


def _select_contrast(config, contrast_name: str | None):
    if contrast_name:
        contrast = config.get_contrast(contrast_name)
        if contrast is None:
            raise ValueError(f"Contrast not found in config: {contrast_name}")
        return contrast
    if len(config.contrasts) != 1:
        raise ValueError(
            "v1 supports one simple contrast per PyDESeq2 run; pass --contrast-name"
        )
    return config.contrasts[0]


def run(args: argparse.Namespace) -> int:
    logger = configure_logging(args.outdir, "04_run_pydeseq2")
    logger.info("Starting PyDESeq2 stage")

    pydeseq2_dir = args.outdir / "pydeseq2"
    pydeseq2_dir.mkdir(parents=True, exist_ok=True)

    try:
        config = load_and_validate_config(args.config)
        contrast = _select_contrast(config, args.contrast_name)

        if contrast.variable != config.condition_column:
            raise ValueError(
                f"v1 only supports contrasts on condition_column '{config.condition_column}', "
                f"got '{contrast.variable}'"
            )

        counts_df = load_counts(args.counts_filtered)
        metadata_df = load_metadata(args.metadata)

        counts_for_pydeseq2, metadata_for_pydeseq2 = prepare_pydeseq2_inputs(
            counts_df,
            metadata_df,
            condition_column=config.condition_column,
            reference_level=config.reference_level,
        )

        dds = run_pydeseq2_model(
            counts_for_pydeseq2,
            metadata_for_pydeseq2,
            condition_column=config.condition_column,
            reference_level=config.reference_level,
            cooks_filter=config.pydeseq2.cooks_filter,
        )

        results_df = run_pydeseq2_stats(
            dds,
            condition_column=config.condition_column,
            numerator=contrast.numerator,
            denominator=contrast.denominator,
            alpha=config.pydeseq2.alpha,
            cooks_filter=config.pydeseq2.cooks_filter,
        )
    except Exception as exc:
        logger.error("PyDESeq2 stage failed: %s", exc)
        return 1

    output_files: list[Path] = []
    results_path = pydeseq2_dir / f"results_{contrast.name}.tsv"
    save_dataframe(results_df, results_path, include_index=False)
    output_files.append(results_path)

    extracted = extract_pydeseq2_outputs(dds)
    output_names = {
        "size_factors": "size_factors.tsv",
        "normalized_counts": "normalized_counts.tsv",
        "dispersion_estimates": "dispersion_estimates.tsv",
    }
    for key, filename in output_names.items():
        value = extracted.get(key)
        if value is None:
            continue
        path = pydeseq2_dir / filename
        save_dataframe(value, path, index_label=value.index.name)
        output_files.append(path)

    manifest_path = write_manifest(
        pydeseq2_dir,
        "04_run_pydeseq2",
        config=config.model_dump(),
        input_files=[args.counts_filtered, args.metadata, args.config],
        output_files=output_files,
        extra_info={
            "contrast": contrast.model_dump(),
            "condition_column": config.condition_column,
            "reference_level": config.reference_level,
        },
        filename="pydeseq2_run_manifest.json",
    )
    logger.info("Written manifest to %s", manifest_path)
    logger.info("PyDESeq2 stage completed successfully")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
