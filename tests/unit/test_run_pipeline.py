"""Unit tests for the pipeline orchestrator."""

from __future__ import annotations

from argparse import Namespace

from deg_pipeline.config import load_and_validate_config
from scripts.run_pipeline import build_stage_commands, slice_stages


def test_build_stage_commands_preserves_expected_order(
    tmp_path, toy_config_file, toy_counts_path, toy_metadata_path, toy_annotation_path, toy_go_map_path
):
    config = load_and_validate_config(toy_config_file)
    args = Namespace(
        config=toy_config_file,
        counts=toy_counts_path,
        metadata=toy_metadata_path,
        annotation=toy_annotation_path,
        gene2go=toy_go_map_path,
        go_obo=None,
        outdir=tmp_path / "results",
        force=False,
    )

    commands = build_stage_commands(args, config)
    stage_names = [stage for stage, _ in commands]

    assert stage_names[0] == "01_validate_inputs"
    assert "04_run_pydeseq2" in stage_names
    assert stage_names[-1] == "08_make_report"


def test_slice_stages_accepts_name_substrings():
    commands = [
        ("01_validate_inputs", ["one"]),
        ("02_qc_counts", ["two"]),
        ("03_filter_normalize", ["three"]),
    ]

    sliced = slice_stages(commands, "qc", "filter")

    assert [stage for stage, _ in sliced] == ["02_qc_counts", "03_filter_normalize"]
