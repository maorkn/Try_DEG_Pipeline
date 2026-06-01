#!/usr/bin/env python3
"""Run the complete DEG pipeline.

This script orchestrates all pipeline stages from a config file.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

# Add src to path for development
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from deg_pipeline.config import load_and_validate_config
from deg_pipeline.logging_utils import configure_logging, write_manifest


STAGES = [
    ("01_validate_inputs", "scripts/01_validate_inputs.py"),
    ("02_qc_counts", "scripts/02_qc_counts.py"),
    ("03_filter_normalize", "scripts/03_filter_normalize.py"),
    ("04_run_pydeseq2", "scripts/04_run_pydeseq2.py"),
    ("06_make_deg_outputs", "scripts/06_make_deg_outputs.py"),
    ("07_run_go_enrichment", "scripts/07_run_go_enrichment.py"),
    ("08_make_report", "scripts/08_make_report.py"),
]


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--config",
        required=True,
        type=Path,
        help="Path to YAML configuration file",
    )
    parser.add_argument(
        "--counts",
        required=True,
        type=Path,
        help="Path to gene count matrix",
    )
    parser.add_argument(
        "--metadata",
        required=True,
        type=Path,
        help="Path to sample metadata",
    )
    parser.add_argument(
        "--annotation",
        type=Path,
        help="Path to gene annotation (optional)",
    )
    parser.add_argument(
        "--gene2go",
        type=Path,
        help="Path to local offline gene-to-GO mapping TSV for GO enrichment",
    )
    parser.add_argument(
        "--go-obo",
        type=Path,
        help="Optional path to go-basic.obo for GO provenance",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Output directory for results",
    )
    parser.add_argument(
        "--start-at",
        help="Start at a specific stage (stage name or number)",
    )
    parser.add_argument(
        "--stop-after",
        help="Stop after a specific stage (stage name or number)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output files",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing",
    )
    return parser


def run_command(cmd: list[str], log_file: Path, dry_run: bool = False) -> int:
    """Run a command and capture output.

    Args:
        cmd: Command to run as list of strings.
        log_file: Path to log file.
        dry_run: If True, only print command.

    Returns:
        Exit code.
    """
    if dry_run:
        print(f"Would run: {' '.join(cmd)}")
        return 0

    log_file.parent.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Write output to log
    with open(log_file, "w") as f:
        f.write(f"Command: {' '.join(cmd)}\n")
        f.write(f"Exit code: {result.returncode}\n")
        f.write(f"Stdout:\n{result.stdout}\n")
        f.write(f"Stderr:\n{result.stderr}\n")

    return result.returncode


def build_stage_commands(args: argparse.Namespace, config) -> list[tuple[str, list[str]]]:
    """Build explicit commands for each pipeline stage."""

    python = sys.executable
    commands: list[tuple[str, list[str]]] = []

    validate_cmd = [
        python,
        "scripts/01_validate_inputs.py",
        "--counts",
        str(args.counts),
        "--metadata",
        str(args.metadata),
        "--config",
        str(args.config),
        "--outdir",
        str(args.outdir),
    ]
    if args.annotation:
        validate_cmd.extend(["--annotation", str(args.annotation)])
    if args.force:
        validate_cmd.append("--force")
    commands.append(("01_validate_inputs", validate_cmd))

    ordered_counts = args.outdir / "intermediate" / "counts_ordered.tsv"
    ordered_metadata = args.outdir / "intermediate" / "metadata_ordered.tsv"
    filtered_counts = args.outdir / "intermediate" / "counts_filtered.tsv"

    for stage_name, script, counts_arg, counts_path in [
        ("02_qc_counts", "scripts/02_qc_counts.py", "--counts", ordered_counts),
        ("03_filter_normalize", "scripts/03_filter_normalize.py", "--counts", ordered_counts),
    ]:
        cmd = [
            python,
            script,
            counts_arg,
            str(counts_path),
            "--metadata",
            str(ordered_metadata),
            "--config",
            str(args.config),
            "--outdir",
            str(args.outdir),
        ]
        if args.force:
            cmd.append("--force")
        commands.append((stage_name, cmd))

    pydeseq2_cmd = [
        python,
        "scripts/04_run_pydeseq2.py",
        "--counts-filtered",
        str(filtered_counts),
        "--metadata",
        str(ordered_metadata),
        "--config",
        str(args.config),
        "--outdir",
        str(args.outdir),
    ]
    if len(config.contrasts) == 1:
        pydeseq2_cmd.extend(["--contrast-name", config.contrasts[0].name])
    if args.force:
        pydeseq2_cmd.append("--force")
    commands.append(("04_run_pydeseq2", pydeseq2_cmd))

    deg_cmd = [
        python,
        "scripts/06_make_deg_outputs.py",
        "--pydeseq2-results-dir",
        str(args.outdir / "pydeseq2"),
        "--config",
        str(args.config),
        "--outdir",
        str(args.outdir),
    ]
    if args.annotation:
        deg_cmd.extend(["--annotation", str(args.annotation)])
    if args.force:
        deg_cmd.append("--force")
    commands.append(("06_make_deg_outputs", deg_cmd))

    if args.annotation and args.gene2go:
        go_cmd = [
            python,
            "scripts/07_run_go_enrichment.py",
            "--deg-dir",
            str(args.outdir / "deg"),
            "--annotation",
            str(args.annotation),
            "--gene2go",
            str(args.gene2go),
            "--config",
            str(args.config),
            "--outdir",
            str(args.outdir),
        ]
        if args.go_obo:
            go_cmd.extend(["--go-obo", str(args.go_obo)])
        if args.force:
            go_cmd.append("--force")
        commands.append(("07_run_go_enrichment", go_cmd))

    report_cmd = [
        python,
        "scripts/08_make_report.py",
        "--project-dir",
        str(args.outdir),
        "--config",
        str(args.config),
        "--outdir",
        str(args.outdir),
    ]
    if args.force:
        report_cmd.append("--force")
    commands.append(("08_make_report", report_cmd))

    return commands


def _stage_index(stage_commands: list[tuple[str, list[str]]], selector: str) -> int:
    if selector.isdigit():
        idx = int(selector)
        if 0 <= idx < len(stage_commands):
            return idx
        if 1 <= idx <= len(stage_commands):
            return idx - 1
    for idx, (stage_name, _) in enumerate(stage_commands):
        if selector == stage_name or selector in stage_name:
            return idx
    raise ValueError(f"Unknown stage selector: {selector}")


def slice_stages(
    stage_commands: list[tuple[str, list[str]]],
    start_at: str | None,
    stop_after: str | None,
) -> list[tuple[str, list[str]]]:
    """Slice stage commands by optional start/stop selectors."""

    start = _stage_index(stage_commands, start_at) if start_at else 0
    stop = _stage_index(stage_commands, stop_after) if stop_after else len(stage_commands) - 1
    if start > stop:
        raise ValueError("--start-at must not come after --stop-after")
    return stage_commands[start : stop + 1]


def get_stage_outputs(stage_name: str, outdir: Path) -> list[Path]:
    """Return existing outputs likely created by a stage."""

    stage_dirs = {
        "01_validate_inputs": [outdir / "validation", outdir / "intermediate"],
        "02_qc_counts": [outdir / "qc"],
        "03_filter_normalize": [outdir / "filtering", outdir / "intermediate"],
        "04_run_pydeseq2": [outdir / "pydeseq2"],
        "06_make_deg_outputs": [outdir / "deg"],
        "07_run_go_enrichment": [outdir / "go"],
        "08_make_report": [outdir / "report"],
    }
    outputs: list[Path] = []
    for directory in stage_dirs.get(stage_name, []):
        if directory.exists():
            outputs.extend(path for path in directory.rglob("*") if path.is_file())
    return outputs


def run(args: argparse.Namespace) -> int:
    """Run the orchestrator pipeline."""
    logger = configure_logging(args.outdir, "run_pipeline")
    logger.info("Starting DEG pipeline orchestrator")

    # Load config
    config = load_and_validate_config(args.config)
    config_dict = config.model_dump()

    # Build stage commands
    stage_commands = build_stage_commands(args, config)

    # Apply start/stop slicing
    stage_commands = slice_stages(stage_commands, args.start_at, args.stop_after)

    # Run commands
    all_outputs = []
    for stage_name, cmd in stage_commands:
        logger.info(f"Running stage: {stage_name}")

        if args.dry_run:
            print(f"Would run: {' '.join(cmd)}")
            continue

        exit_code = run_command(cmd, args.outdir / f"{stage_name}.log")

        if exit_code != 0:
            logger.error(f"Stage {stage_name} failed with exit code {exit_code}")
            return 1

        all_outputs.extend(get_stage_outputs(stage_name, args.outdir))

    # Write master manifest
    manifest_path = write_manifest(
        args.outdir,
        "run_pipeline",
        config=config_dict,
        input_files=[args.counts, args.metadata] + ([args.annotation] if args.annotation else []),
        output_files=all_outputs,
        extra_info={
            "start_at": args.start_at,
            "stop_after": args.stop_after,
            "dry_run": args.dry_run,
        },
    )
    logger.info(f"Written master manifest to {manifest_path}")

    logger.info("Pipeline orchestration completed successfully")
    return 0


def main(argv: list[str] | None = None) -> int:
    """Main entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
