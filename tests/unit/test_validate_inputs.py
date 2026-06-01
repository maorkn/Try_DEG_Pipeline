"""Unit tests for input validation functions."""

from __future__ import annotations

import pandas as pd
import pytest

from deg_pipeline.validation import (
    ValidationResult,
    reorder_samples,
    validate_annotation,
    validate_count_matrix,
    validate_inputs,
    validate_metadata,
)


class TestValidationResult:
    """Tests for ValidationResult class."""

    def test_empty_result_is_valid(self):
        """Empty validation result should be valid."""
        result = ValidationResult()
        assert result.is_valid
        assert not result.has_errors
        assert not result.has_warnings

    def test_result_with_error_is_invalid(self):
        """Result with errors should be invalid."""
        result = ValidationResult()
        result.add_error("test error")
        assert not result.is_valid
        assert result.has_errors

    def test_result_with_warning_is_still_valid(self):
        """Result with only warnings should still be valid."""
        result = ValidationResult()
        result.add_warning("test warning")
        assert result.is_valid
        assert result.has_warnings

    def test_merge_combines_results(self):
        """Merging results should combine errors and warnings."""
        result1 = ValidationResult()
        result1.add_error("error 1")
        result1.metrics["a"] = 1

        result2 = ValidationResult()
        result2.add_warning("warning 1")
        result2.metrics["b"] = 2

        result1.merge(result2)
        assert len(result1.errors) == 1
        assert len(result1.warnings) == 1
        assert result1.metrics == {"a": 1, "b": 2}


class TestValidateCountMatrix:
    """Tests for validate_count_matrix function."""

    def test_valid_counts_passes(self, toy_counts_df):
        """Valid count matrix should pass validation."""
        result = validate_count_matrix(toy_counts_df)
        assert result.is_valid
        assert result.metrics["n_genes"] == 10
        assert result.metrics["n_samples"] == 4

    def test_empty_counts_fails(self):
        """Empty count matrix should fail."""
        df = pd.DataFrame()
        result = validate_count_matrix(df)
        assert not result.is_valid
        assert any("empty" in e.lower() for e in result.errors)

    def test_negative_counts_fails(self):
        """Count matrix with negative values should fail."""
        df = pd.DataFrame(
            {"sample_1": [100, -5], "sample_2": [120, 50]},
            index=["gene_1", "gene_2"],
        )
        result = validate_count_matrix(df)
        assert not result.is_valid
        assert any("negative" in e.lower() for e in result.errors)

    def test_duplicate_gene_ids_fails(self):
        """Count matrix with duplicate gene IDs should fail."""
        df = pd.DataFrame(
            {"sample_1": [100, 50], "sample_2": [120, 60]},
            index=["gene_1", "gene_1"],
        )
        result = validate_count_matrix(df)
        assert not result.is_valid
        assert any("duplicate" in e.lower() for e in result.errors)

    def test_duplicate_sample_names_fails(self):
        """Count matrix with duplicate sample names should fail."""
        df = pd.DataFrame(
            {"sample_1": [100, 50], "sample_1": [120, 60]},
            index=["gene_1", "gene_2"],
        )
        # Note: pandas will actually only keep one column with duplicate names
        # So we need to test differently
        df = pd.DataFrame([[100, 120], [50, 60]], index=["gene_1", "gene_2"])
        df.columns = ["sample_1", "sample_1"]
        result = validate_count_matrix(df)
        assert not result.is_valid
        assert any("duplicate" in e.lower() for e in result.errors)

    def test_nan_values_fails(self):
        """Count matrix with NaN values should fail."""
        df = pd.DataFrame(
            {"sample_1": [100, float("nan")], "sample_2": [120, 50]},
            index=["gene_1", "gene_2"],
        )
        result = validate_count_matrix(df)
        assert not result.is_valid
        assert any("nan" in e.lower() for e in result.errors)

    def test_non_integer_warns(self):
        """Count matrix with many non-integer values should warn."""
        df = pd.DataFrame(
            {"sample_1": [100.5, 50.3, 30.7], "sample_2": [120.2, 60.8, 40.1]},
            index=["gene_1", "gene_2", "gene_3"],
        )
        result = validate_count_matrix(df)
        assert result.is_valid  # Still passes, but with warning
        assert any("non-integer" in w.lower() for w in result.warnings)

    def test_metrics_recorded(self, minimal_counts_df):
        """Validation should record count metrics."""
        result = validate_count_matrix(minimal_counts_df)
        assert "n_genes" in result.metrics
        assert "n_samples" in result.metrics
        assert "total_counts" in result.metrics
        assert "zero_fraction" in result.metrics


class TestValidateMetadata:
    """Tests for validate_metadata function."""

    def test_valid_metadata_passes(self, toy_metadata_df, toy_counts_df, toy_config):
        """Valid metadata should pass validation."""
        result = validate_metadata(toy_metadata_df, toy_counts_df.columns, toy_config)
        assert result.is_valid
        assert result.metrics["n_samples"] == 4

    def test_empty_metadata_fails(self, toy_counts_df, toy_config):
        """Empty metadata should fail."""
        df = pd.DataFrame(columns=["condition"]).set_index(pd.Index([], name="sample_id"))
        result = validate_metadata(df, toy_counts_df.columns, toy_config)
        assert not result.is_valid

    def test_duplicate_sample_ids_fails(self, toy_counts_df, toy_config):
        """Metadata with duplicate sample IDs should fail."""
        df = pd.DataFrame(
            {"condition": ["control", "control"]},
            index=pd.Index(["sample_1", "sample_1"], name="sample_id"),
        )
        result = validate_metadata(df, toy_counts_df.columns, toy_config)
        assert not result.is_valid
        assert any("duplicate" in e.lower() for e in result.errors)

    def test_missing_condition_column_fails(self, toy_counts_df):
        """Metadata missing condition column should fail."""
        df = pd.DataFrame(
            {"other_col": ["a", "b", "c", "d"]},
            index=pd.Index(["sample_1", "sample_2", "sample_3", "sample_4"], name="sample_id"),
        )
        config = {"condition_column": "condition", "contrasts": []}
        result = validate_metadata(df, toy_counts_df.columns, config)
        assert not result.is_valid
        assert any("condition" in e.lower() for e in result.errors)

    def test_missing_samples_fails(self, toy_config):
        """Metadata missing samples from counts should fail."""
        df = pd.DataFrame(
            {"condition": ["control", "treated"]},
            index=pd.Index(["sample_1", "sample_3"], name="sample_id"),
        )
        sample_names = ["sample_1", "sample_2", "sample_3", "sample_4"]
        result = validate_metadata(df, sample_names, toy_config)
        assert not result.is_valid
        assert any("not in metadata" in e.lower() for e in result.errors)

    def test_invalid_contrast_level_fails(self, toy_metadata_df, toy_counts_df):
        """Contrast with invalid levels should fail."""
        config = {
            "condition_column": "condition",
            "reference_level": "control",
            "contrasts": [
                {
                    "name": "test",
                    "variable": "condition",
                    "numerator": "nonexistent",
                    "denominator": "control",
                }
            ],
        }
        result = validate_metadata(toy_metadata_df, toy_counts_df.columns, config)
        assert not result.is_valid
        assert any("numerator" in e.lower() for e in result.errors)

    def test_single_condition_level_fails(self, toy_counts_df):
        """Metadata with only one condition level should fail."""
        df = pd.DataFrame(
            {"condition": ["control", "control", "control", "control"]},
            index=pd.Index(["sample_1", "sample_2", "sample_3", "sample_4"], name="sample_id"),
        )
        config = {"condition_column": "condition", "contrasts": []}
        result = validate_metadata(df, toy_counts_df.columns, config)
        assert not result.is_valid
        assert any("fewer than 2" in e.lower() for e in result.errors)


class TestValidateAnnotation:
    """Tests for validate_annotation function."""

    def test_valid_annotation_passes(self, toy_annotation_df, toy_counts_df):
        """Valid annotation should pass validation."""
        result = validate_annotation(toy_annotation_df, toy_counts_df.index)
        assert result.is_valid
        assert result.metrics["overlap_rate"] == 1.0

    def test_none_annotation_warns(self, toy_counts_df):
        """None annotation should warn."""
        result = validate_annotation(None, toy_counts_df.index)
        assert result.is_valid  # Still valid, just a warning
        assert result.has_warnings

    def test_low_overlap_warns(self, toy_counts_df):
        """Annotation with low gene overlap should warn."""
        df = pd.DataFrame(
            {"gene_symbol": ["OTHER1", "OTHER2"]},
            index=pd.Index(["OTHER_GENE_1", "OTHER_GENE_2"], name="gene_id"),
        )
        result = validate_annotation(df, toy_counts_df.index)
        assert result.is_valid  # Still valid
        assert result.has_warnings
        assert any("overlap" in w.lower() for w in result.warnings)

    def test_missing_useful_columns_warns(self, toy_counts_df):
        """Annotation without useful columns should warn."""
        df = pd.DataFrame(
            {"other_column": [1, 2, 3]},
            index=pd.Index(["ENSG000001", "ENSG000002", "ENSG000003"], name="gene_id"),
        )
        result = validate_annotation(df, toy_counts_df.index)
        assert result.has_warnings


class TestValidateInputs:
    """Tests for validate_inputs function."""

    def test_valid_inputs_pass(self, toy_counts_df, toy_metadata_df, toy_config):
        """Valid inputs should pass validation."""
        result = validate_inputs(toy_counts_df, toy_metadata_df, toy_config)
        assert result.is_valid

    def test_invalid_counts_stops_early(self, toy_metadata_df, toy_config):
        """Invalid counts should stop validation early."""
        bad_counts = pd.DataFrame(
            {"sample_1": [-100, 50], "sample_2": [120, 60]},
            index=["gene_1", "gene_2"],
        )
        result = validate_inputs(bad_counts, toy_metadata_df, toy_config)
        assert not result.is_valid


class TestReorderSamples:
    """Tests for reorder_samples function."""

    def test_reorder_matches_counts(self, toy_counts_df, toy_metadata_df):
        """Reordered metadata should match count sample order."""
        counts, metadata = reorder_samples(toy_counts_df, toy_metadata_df)
        assert list(metadata.index) == list(counts.columns)

    def test_reorder_preserves_data(self, toy_counts_df, toy_metadata_df):
        """Reordering should preserve metadata values."""
        original_conditions = set(toy_metadata_df["condition"])
        _, reordered = reorder_samples(toy_counts_df, toy_metadata_df)
        reordered_conditions = set(reordered["condition"])
        assert original_conditions == reordered_conditions