"""Unit tests for QC functions."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from deg_pipeline.qc import (
    compute_gene_metrics,
    compute_library_metrics,
    compute_pca,
    compute_qc_summary,
    compute_sample_correlation,
)


class TestComputeLibraryMetrics:
    """Tests for compute_library_metrics function."""

    def test_basic_metrics(self, minimal_counts_df):
        """Library metrics should have correct structure."""
        metrics = compute_library_metrics(minimal_counts_df)
        assert "total_counts" in metrics.columns
        assert "detected_genes" in metrics.columns
        assert "zero_fraction" in metrics.columns
        assert "upper_quartile" in metrics.columns
        assert len(metrics) == 4  # 4 samples

    def test_total_counts_correct(self, minimal_counts_df):
        """Total counts should sum correctly."""
        metrics = compute_library_metrics(minimal_counts_df)
        expected_total = minimal_counts_df["sample_1"].sum()
        assert metrics.loc["sample_1", "total_counts"] == expected_total

    def test_detected_genes_counts_nonzero(self, minimal_counts_df):
        """Detected genes should count non-zero values."""
        metrics = compute_library_metrics(minimal_counts_df)
        # gene_diff has 0 in sample_1, so detected should be 3 for sample_1
        assert metrics.loc["sample_1", "detected_genes"] == 3

    def test_zero_fraction(self, minimal_counts_df):
        """Zero fraction should be calculated correctly."""
        metrics = compute_library_metrics(minimal_counts_df)
        # sample_1 has 1 zero out of 4 genes = 0.25
        assert metrics.loc["sample_1", "zero_fraction"] == 0.25


class TestComputeGeneMetrics:
    """Tests for compute_gene_metrics function."""

    def test_basic_metrics(self, minimal_counts_df):
        """Gene metrics should have correct structure."""
        metrics = compute_gene_metrics(minimal_counts_df)
        assert "total_counts" in metrics.columns
        assert "mean_counts" in metrics.columns
        assert "detected_samples" in metrics.columns
        assert len(metrics) == 4  # 4 genes

    def test_total_counts_correct(self, minimal_counts_df):
        """Total counts should sum correctly across samples."""
        metrics = compute_gene_metrics(minimal_counts_df)
        expected = minimal_counts_df.loc["gene_high"].sum()
        assert metrics.loc["gene_high", "total_counts"] == expected

    def test_detected_samples(self, minimal_counts_df):
        """Detected samples should count correctly."""
        metrics = compute_gene_metrics(minimal_counts_df)
        # gene_diff has 0 in sample_1, so detected in 3 samples
        assert metrics.loc["gene_diff", "detected_samples"] == 3


class TestComputeSampleCorrelation:
    """Tests for compute_sample_correlation function."""

    def test_correlation_shape(self, minimal_counts_df):
        """Correlation matrix should be square."""
        corr = compute_sample_correlation(minimal_counts_df)
        assert corr.shape == (4, 4)

    def test_diagonal_is_one(self, minimal_counts_df):
        """Diagonal of correlation matrix should be 1."""
        corr = compute_sample_correlation(minimal_counts_df)
        for sample in corr.index:
            assert corr.loc[sample, sample] == 1.0

    def test_correlation_symmetric(self, minimal_counts_df):
        """Correlation matrix should be symmetric."""
        corr = compute_sample_correlation(minimal_counts_df)
        assert np.allclose(corr.values, corr.T.values)

    def test_correlation_range(self, minimal_counts_df):
        """Correlation values should be between -1 and 1."""
        corr = compute_sample_correlation(minimal_counts_df)
        assert corr.min().min() >= -1
        assert corr.max().max() <= 1


class TestComputePCA:
    """Tests for compute_pca function."""

    def test_pca_shape(self, minimal_counts_df):
        """PCA output should have correct shape."""
        pca_df, variance = compute_pca(minimal_counts_df)
        assert "PC1" in pca_df.columns
        assert "PC2" in pca_df.columns
        assert len(pca_df) == 4  # 4 samples
        assert len(variance) == 2

    def test_variance_explained_valid(self, minimal_counts_df):
        """Variance explained should be valid probabilities."""
        _, variance = compute_pca(minimal_counts_df)
        assert all(0 <= v <= 1 for v in variance)
        assert sum(variance) <= 1

    def test_pca_with_metadata(self, minimal_counts_df, minimal_metadata_df):
        """PCA should include metadata when provided."""
        pca_df, _ = compute_pca(minimal_counts_df, metadata_df=minimal_metadata_df)
        assert "condition" in pca_df.columns


class TestComputeQCSummary:
    """Tests for compute_qc_summary function."""

    def test_summary_structure(self, minimal_counts_df):
        """QC summary should have expected keys."""
        library_metrics = compute_library_metrics(minimal_counts_df)
        gene_metrics = compute_gene_metrics(minimal_counts_df)
        summary = compute_qc_summary(minimal_counts_df, library_metrics, gene_metrics)

        assert "n_samples" in summary
        assert "n_genes" in summary
        assert "total_counts" in summary
        assert "median_library_size" in summary
        assert "genes_not_detected" in summary

    def test_summary_values(self, minimal_counts_df):
        """QC summary values should be correct."""
        library_metrics = compute_library_metrics(minimal_counts_df)
        gene_metrics = compute_gene_metrics(minimal_counts_df)
        summary = compute_qc_summary(minimal_counts_df, library_metrics, gene_metrics)

        assert summary["n_samples"] == 4
        assert summary["n_genes"] == 4