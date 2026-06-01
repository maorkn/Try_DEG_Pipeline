"""Unit tests for filtering functions."""

from __future__ import annotations

import pandas as pd
import pytest

from deg_pipeline.filtering import filter_low_count_genes


class TestFilterLowCountGenes:
    """Tests for filter_low_count_genes function."""

    def test_filter_removes_low_count_genes(self):
        """Test that low count genes are removed."""
        # Create test data
        counts = pd.DataFrame({
            "sample1": [100, 5, 50, 1],
            "sample2": [200, 10, 100, 2],
            "sample3": [150, 8, 75, 1],
        }, index=["gene1", "gene2", "gene3", "gene4"])

        # Apply filtering with min_count=10, min_samples=2
        filtered, mask = filter_low_count_genes(counts, min_count=10, min_samples=2)

        # Check that genes with sufficient counts are retained
        assert "gene1" in filtered.index  # 100, 200, 150 - all >= 10
        assert "gene3" in filtered.index  # 50, 100, 75 - all >= 10
        assert "gene2" not in filtered.index  # 5, 10, 8 - only 1 sample >= 10
        assert "gene4" not in filtered.index  # 1, 2, 1 - no samples >= 10

    def test_filter_handles_edge_cases(self):
        """Test filtering with edge cases."""
        # Empty DataFrame
        empty_df = pd.DataFrame()
        filtered, mask = filter_low_count_genes(empty_df, min_count=10, min_samples=2)
        assert len(filtered) == 0

        # All zeros
        zeros_df = pd.DataFrame({
            "sample1": [0, 0, 0],
            "sample2": [0, 0, 0],
        }, index=["gene1", "gene2", "gene3"])
        filtered, mask = filter_low_count_genes(zeros_df, min_count=10, min_samples=2)
        assert len(filtered) == 0