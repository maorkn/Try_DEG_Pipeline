"""Unit tests for filtering functions."""

from __future__ import annotations

import pandas as pd
import pytest

from deg_pipeline.filtering import (
    compute_filtering_summary,
    determine_min_samples,
    filter_low_count_genes,
)


class TestDetermineMinSamples:
    """Tests for determine_min_samples function."""

    def test_configured_value_used(self, minimal_metadata_df):
        """Configured min_samples should be used when provided."""
        result = determine_min_samples(minimal_metadata_df, "condition", 5)
        assert result == 5

    def test_smallest_group_used_when_none(self, minimal_metadata_df):
        """Smallest group size should be used when config is None."""
        result = determine_min_samples(minimal_metadata_df, "condition", None)
        # Both groups have 2 samples
        assert result == 2

    def test_unbalanced_groups(self):
        """Should use smallest group size for unbalanced data."""
        df = pd.DataFrame(
            {"condition": ["a", "a", "a", "b"]},
            index=pd.Index(["s1", "s2", "s3", "s4"], name="sample_id"),
        )
        result = determine_min_samples(df, "condition", None)
        assert result == 1


class TestFilterLowCountGenes:
    """Tests for filter_low_count_genes function."""

    def test_filter_removes_low_genes(self, minimal_counts_df):
        """Filter should remove genes with low counts."""
        filtered, summary = filter_low_count_genes(
            minimal_counts_df, min_count=50, min_samples=2
        )
        # gene_high and gene_med have counts >= 50 in 2+ samples
        assert "gene_high" in filtered.index
        assert "gene_med" in filtered.index
        # gene_low has max count of 15, should be removed
        assert "gene_low" not in filtered.index

    def test_filter_keeps_all_with_low_threshold(self, minimal_counts_df):
        """Filter should keep all genes with very low threshold."""
        filtered, summary = filter_low_count_genes(
            minimal_counts_df, min_count=0, min_samples=1
        )
        # All genes with at least one non-zero count should pass
        assert len(filtered) >= 3

    def test_summary_has_correct_columns(self, minimal_counts_df):
        """Filtering summary should have expected columns."""
        _, summary = filter_low_count_genes(minimal_counts_df, min_count=10, min_samples=2)
        assert "samples_above_threshold" in summary.columns
        assert "passes_filter" in summary.columns
        assert "mean_count" in summary.columns

    def test_summary_passes_filter_correct(self, minimal_counts_df):
        """Summary passes_filter should match filtered output."""
        filtered, summary = filter_low_count_genes(
            minimal_counts_df, min_count=10, min_samples=2
        )
        passed_genes = summary[summary["passes_filter"]].index
        assert set(passed_genes) == set(filtered.index)


class TestComputeFilteringSummary:
    """Tests for compute_filtering_summary function."""

    def test_summary_structure(self, minimal_counts_df):
        """Filtering summary should have expected keys."""
        filtered, _ = filter_low_count_genes(minimal_counts_df, min_count=10, min_samples=2)
        summary = compute_filtering_summary(
            minimal_counts_df, filtered, min_count=10, min_samples=2
        )

        assert "n_genes_before_filtering" in summary
        assert "n_genes_after_filtering" in summary
        assert "n_genes_removed" in summary
        assert "fraction_removed" in summary
        assert "min_count_threshold" in summary
        assert "min_samples_threshold" in summary

    def test_summary_values_correct(self, minimal_counts_df):
        """Filtering summary values should be correct."""
        filtered, _ = filter_low_count_genes(minimal_counts_df, min_count=10, min_samples=2)
        summary = compute_filtering_summary(
            minimal_counts_df, filtered, min_count=10, min_samples=2
        )

        assert summary["n_genes_before_filtering"] == 4
        assert summary["n_genes_after_filtering"] == len(filtered)
        assert summary["n_genes_removed"] == 4 - len(filtered)
        assert summary["min_count_threshold"] == 10
        assert summary["min_samples_threshold"] == 2