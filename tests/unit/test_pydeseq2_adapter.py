"""Unit tests for PyDESeq2 adapter functions."""

from __future__ import annotations

import pandas as pd
import pytest

from deg_pipeline.pydeseq2_adapter import prepare_pydeseq2_inputs


class TestPreparePyDESeq2Inputs:
    """Tests for prepare_pydeseq2_inputs function."""

    def test_prepare_with_valid_inputs(self):
        """Test preparation with valid inputs."""
        # Create test data
        counts = pd.DataFrame({
            "sample1": [100, 50, 25],
            "sample2": [200, 100, 50],
        }, index=["gene1", "gene2", "gene3"])

        metadata = pd.DataFrame({
            "condition": ["control", "treated"],
        }, index=["sample1", "sample2"])

        # Prepare inputs
        counts_pydeseq2, metadata_pydeseq2 = prepare_pydeseq2_inputs(
            counts, metadata, condition_column="condition", reference_level="control"
        )

        # Check output structure
        assert counts_pydeseq2.shape[0] == 2  # 2 samples
        assert counts_pydeseq2.shape[1] == 3  # 3 genes
        assert metadata_pydeseq2.shape[0] == 2  # 2 samples
        assert "condition" in metadata_pydeseq2.columns

    def test_prepare_rejects_three_levels(self):
        """Test that three condition levels are rejected."""
        counts = pd.DataFrame({
            "sample1": [100, 50],
            "sample2": [200, 100],
            "sample3": [150, 75],
        }, index=["gene1", "gene2"])

        metadata = pd.DataFrame({
            "condition": ["control", "treated", "other"],
        }, index=["sample1", "sample2", "sample3"])

        with pytest.raises(ValueError, match="exactly two"):
            prepare_pydeseq2_inputs(
                counts, metadata, condition_column="condition", reference_level="control"
            )