"""Validation utilities for pipeline inputs and airway reference checks."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd


@dataclass
class ValidationResult:
    """Container for validation errors, warnings, and metrics."""

    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    metrics: dict[str, Any] = field(default_factory=dict)

    @property
    def is_valid(self) -> bool:
        return not self.errors

    @property
    def has_errors(self) -> bool:
        return bool(self.errors)

    @property
    def has_warnings(self) -> bool:
        return bool(self.warnings)

    def add_error(self, message: str) -> None:
        self.errors.append(message)

    def add_warning(self, message: str) -> None:
        self.warnings.append(message)

    def merge(self, other: "ValidationResult") -> None:
        self.errors.extend(other.errors)
        self.warnings.extend(other.warnings)
        self.metrics.update(other.metrics)

    def to_dict(self) -> dict[str, Any]:
        return {
            "is_valid": self.is_valid,
            "errors": self.errors,
            "warnings": self.warnings,
            "metrics": self.metrics,
        }


def _as_config_dict(config: Any) -> dict[str, Any]:
    if hasattr(config, "model_dump"):
        return config.model_dump()
    return dict(config)


def _finite_numeric_frame(df: pd.DataFrame) -> tuple[pd.DataFrame | None, list[str]]:
    errors: list[str] = []
    numeric = df.apply(pd.to_numeric, errors="coerce")

    if numeric.isna().any().any():
        errors.append("Count matrix contains NaN or non-numeric values")

    if np.isinf(numeric.to_numpy(dtype=float, na_value=np.nan)).any():
        errors.append("Count matrix contains infinite values")

    return (None if errors else numeric), errors


def validate_count_matrix(counts_df: pd.DataFrame) -> ValidationResult:
    """Validate raw count matrix integrity."""

    result = ValidationResult()

    if counts_df.empty:
        result.add_error("Count matrix is empty")
        return result

    result.metrics["n_genes"] = int(counts_df.shape[0])
    result.metrics["n_samples"] = int(counts_df.shape[1])

    if counts_df.index.has_duplicates:
        result.add_error("Count matrix contains duplicate gene IDs")

    if counts_df.columns.has_duplicates:
        result.add_error("Count matrix contains duplicate sample names")

    if any(str(idx).strip() == "" for idx in counts_df.index):
        result.add_error("Count matrix contains empty gene IDs")

    if any(str(col).strip() == "" for col in counts_df.columns):
        result.add_error("Count matrix contains empty sample names")

    numeric_counts, numeric_errors = _finite_numeric_frame(counts_df)
    for error in numeric_errors:
        result.add_error(error)

    if numeric_counts is None:
        return result

    values = numeric_counts.to_numpy()
    if (values < 0).any():
        result.add_error("Count matrix contains negative values")

    non_integer_mask = ~np.isclose(values, np.round(values), rtol=0, atol=1e-8)
    if non_integer_mask.any():
        result.add_warning(
            "Count matrix contains non-integer values; DESeq-style methods expect raw counts"
        )

    result.metrics["total_counts"] = int(np.nansum(values))
    result.metrics["zero_fraction"] = float((values == 0).sum() / values.size)
    result.metrics["min_count"] = float(np.nanmin(values))
    result.metrics["max_count"] = float(np.nanmax(values))

    return result


def validate_metadata(
    metadata_df: pd.DataFrame,
    sample_names: Iterable[str],
    config: Any,
    allow_extra_metadata: bool = False,
) -> ValidationResult:
    """Validate metadata/sample matching and simple two-group contrast scope."""

    cfg = _as_config_dict(config)
    result = ValidationResult()
    sample_names = [str(sample) for sample in sample_names]

    if metadata_df.empty:
        result.add_error("Metadata table is empty")
        return result

    if metadata_df.index.has_duplicates:
        result.add_error("Metadata contains duplicate sample IDs")

    condition_column = cfg.get("condition_column")
    if not condition_column or condition_column not in metadata_df.columns:
        result.add_error(f"Metadata missing condition column: {condition_column}")
        return result

    missing_metadata = [sample for sample in sample_names if sample not in metadata_df.index]
    if missing_metadata:
        result.add_error(f"Count samples not in metadata: {missing_metadata}")

    extra_metadata = [sample for sample in metadata_df.index if sample not in sample_names]
    if extra_metadata and not allow_extra_metadata:
        result.add_error(f"Metadata samples not in counts: {extra_metadata}")
    elif extra_metadata:
        result.add_warning(f"Ignoring metadata samples not in counts: {extra_metadata}")

    conditions = metadata_df.loc[
        [sample for sample in sample_names if sample in metadata_df.index],
        condition_column,
    ].dropna()
    levels = list(pd.Series(conditions).astype(str).unique())

    if len(levels) < 2:
        result.add_error(
            f"Condition '{condition_column}' has fewer than 2 levels; v1 requires exactly two"
        )
    elif len(levels) > 2:
        result.add_error(
            f"Condition '{condition_column}' has more than 2 levels; v1 supports simple two-group contrasts only"
        )

    group_sizes = pd.Series(conditions).astype(str).value_counts().to_dict()
    result.metrics["n_samples"] = int(len(metadata_df))
    result.metrics["condition_levels"] = levels
    result.metrics["condition_group_sizes"] = {str(k): int(v) for k, v in group_sizes.items()}

    for level, size in group_sizes.items():
        if size < 2:
            result.add_warning(
                f"Condition level '{level}' has {size} replicate(s); DE testing is unstable with fewer than 2"
            )

    reference_level = str(cfg.get("reference_level", ""))
    if reference_level and reference_level not in levels:
        result.add_error(f"Reference level '{reference_level}' not found in condition levels: {levels}")

    contrasts = cfg.get("contrasts") or []
    if not contrasts:
        result.add_error("Config must define at least one contrast")

    for contrast in contrasts:
        variable = contrast.get("variable")
        if variable != condition_column:
            result.add_error(
                f"Contrast '{contrast.get('name', '<unnamed>')}' uses unsupported variable '{variable}'; "
                f"v1 only models '{condition_column}'"
            )
        numerator = str(contrast.get("numerator", ""))
        denominator = str(contrast.get("denominator", ""))
        if numerator not in levels:
            result.add_error(
                f"Contrast '{contrast.get('name', '<unnamed>')}' numerator '{numerator}' "
                f"not found in condition levels: {levels}"
            )
        if denominator not in levels:
            result.add_error(
                f"Contrast '{contrast.get('name', '<unnamed>')}' denominator '{denominator}' "
                f"not found in condition levels: {levels}"
            )
        if numerator == denominator and numerator:
            result.add_error(f"Contrast '{contrast.get('name', '<unnamed>')}' has identical levels")

    return result


def validate_annotation(
    annotation_df: pd.DataFrame | None,
    gene_ids: Iterable[str],
) -> ValidationResult:
    """Validate optional gene annotation table."""

    result = ValidationResult()
    gene_ids = set(map(str, gene_ids))

    if annotation_df is None:
        result.add_warning("No annotation table provided")
        return result

    if annotation_df.empty:
        result.add_warning("Annotation table is empty")
        result.metrics["overlap_rate"] = 0.0
        return result

    if annotation_df.index.has_duplicates:
        result.add_error("Annotation contains duplicate gene IDs")

    annotation_gene_ids = set(map(str, annotation_df.index))
    overlap = gene_ids & annotation_gene_ids
    overlap_rate = len(overlap) / len(gene_ids) if gene_ids else 0.0
    result.metrics["annotation_genes"] = int(len(annotation_gene_ids))
    result.metrics["overlap_genes"] = int(len(overlap))
    result.metrics["overlap_rate"] = float(overlap_rate)

    if overlap_rate < 0.5:
        result.add_warning(
            f"Annotation overlap is low ({overlap_rate:.1%}); check gene ID versions/species"
        )

    useful_columns = {"gene_symbol", "entrez_id", "gene_name", "biotype", "go_id", "go_ids"}
    if not (useful_columns & set(annotation_df.columns)):
        result.add_warning("Annotation table lacks common useful columns")

    return result


def validate_inputs(
    counts_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    config: Any,
    annotation_df: pd.DataFrame | None = None,
    allow_extra_metadata: bool = False,
) -> ValidationResult:
    """Run all input validators and merge reports."""

    result = ValidationResult()
    result.merge(validate_count_matrix(counts_df))
    result.merge(
        validate_metadata(
            metadata_df,
            counts_df.columns,
            config,
            allow_extra_metadata=allow_extra_metadata,
        )
    )
    result.merge(validate_annotation(annotation_df, counts_df.index))
    return result


def reorder_samples(
    counts_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return counts and metadata ordered so metadata index matches count columns."""

    missing = [sample for sample in counts_df.columns if sample not in metadata_df.index]
    if missing:
        raise ValueError(f"Cannot reorder metadata; samples missing: {missing}")

    metadata_ordered = metadata_df.loc[list(counts_df.columns)].copy()
    return counts_df.copy(), metadata_ordered


def _read_result_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if "gene_id" not in df.columns:
        first_col = df.columns[0]
        df = df.rename(columns={first_col: "gene_id"})
    return df


def load_pydeseq2_results(results_dir: Path, contrast_name: str) -> pd.DataFrame:
    """Load PyDESeq2 result table for a contrast."""

    path = results_dir / f"results_{contrast_name}.tsv"
    if not path.exists():
        nested = results_dir / "pydeseq2" / f"results_{contrast_name}.tsv"
        path = nested if nested.exists() else path
    if not path.exists():
        raise FileNotFoundError(f"PyDESeq2 results not found: {path}")
    return _read_result_table(path)


def load_deseq2_results(reference_dir: Path, contrast_name: str) -> pd.DataFrame:
    """Load R DESeq2 reference table for a contrast."""

    candidates = [
        reference_dir / f"reference_deseq2_{contrast_name}.tsv",
        reference_dir / f"results_{contrast_name}.tsv",
        reference_dir / "reference_deseq2_dex_trt_vs_untrt.tsv",
    ]
    for path in candidates:
        if path.exists():
            return _read_result_table(path)
    raise FileNotFoundError(f"No R DESeq2 reference results found in {reference_dir}")


def align_results(current: pd.DataFrame, reference: pd.DataFrame) -> pd.DataFrame:
    """Align current and reference result tables by gene ID."""

    required = ["gene_id", "log2FoldChange", "stat", "padj"]
    for label, df in [("current", current), ("reference", reference)]:
        missing = [col for col in required if col not in df.columns]
        if missing:
            raise ValueError(f"{label} results missing columns: {missing}")

    return current.merge(reference, on="gene_id", suffixes=("_current", "_reference"))


def _safe_corr(left: pd.Series, right: pd.Series) -> float:
    data = pd.concat([left, right], axis=1).replace([np.inf, -np.inf], np.nan).dropna()
    if len(data) < 3:
        return float("nan")
    return float(data.iloc[:, 0].corr(data.iloc[:, 1]))


def compute_airway_validation_metrics(
    aligned_df: pd.DataFrame,
    alpha: float = 0.05,
    top_n: int = 100,
    **_: Any,
) -> dict[str, Any]:
    """Compute robust validation metrics for PyDESeq2 versus R DESeq2."""

    lfc_corr = _safe_corr(
        aligned_df["log2FoldChange_current"], aligned_df["log2FoldChange_reference"]
    )
    stat_corr = _safe_corr(aligned_df["stat_current"], aligned_df["stat_reference"])

    current_sig = set(aligned_df.loc[aligned_df["padj_current"] < alpha, "gene_id"])
    reference_sig = set(aligned_df.loc[aligned_df["padj_reference"] < alpha, "gene_id"])

    top_current = set(
        aligned_df.sort_values("padj_current", na_position="last").head(top_n)["gene_id"]
    )
    top_reference = set(
        aligned_df.sort_values("padj_reference", na_position="last").head(top_n)["gene_id"]
    )
    top_union = top_current | top_reference

    return {
        "n_aligned_genes": int(len(aligned_df)),
        "lfc_correlation": lfc_corr,
        "stat_correlation": stat_corr,
        "current_significant": int(len(current_sig)),
        "reference_significant": int(len(reference_sig)),
        "significant_overlap": int(len(current_sig & reference_sig)),
        "top_gene_overlap_rate": float(len(top_current & top_reference) / len(top_union))
        if top_union
        else 0.0,
    }


def evaluate_validation(metrics: dict[str, Any], args: Any) -> dict[str, Any]:
    """Evaluate airway validation metrics against configured thresholds."""

    reference_sig = metrics["reference_significant"]
    current_sig = metrics["current_significant"]
    if reference_sig == 0:
        sig_count_pass = current_sig == 0
    else:
        sig_count_pass = (
            abs(current_sig - reference_sig) / reference_sig
            <= args.significant_count_tolerance
        )

    checks = {
        "lfc_correlation": bool(
            not math.isnan(metrics["lfc_correlation"])
            and metrics["lfc_correlation"] >= args.lfc_correlation_threshold
        ),
        "stat_correlation": bool(
            not math.isnan(metrics["stat_correlation"])
            and metrics["stat_correlation"] >= args.stat_correlation_threshold
        ),
        "significant_count": bool(sig_count_pass),
        "aligned_gene_count": metrics["n_aligned_genes"] > 0,
    }

    return {"passed": all(checks.values()), "checks": checks, "metrics": metrics}


def write_validation_report(
    report: dict[str, Any],
    json_path: Path,
    text_path: Path,
) -> None:
    """Write machine-readable and text airway validation reports."""

    json_path.parent.mkdir(parents=True, exist_ok=True)
    with open(json_path, "w") as f:
        json.dump(report, f, indent=2, default=str)

    with open(text_path, "w") as f:
        f.write("AIRWAY PYDESEQ2 VALIDATION\n")
        f.write(f"Status: {'PASSED' if report['passed'] else 'FAILED'}\n\n")
        f.write("Checks:\n")
        for name, passed in report["checks"].items():
            f.write(f"  - {name}: {'PASS' if passed else 'FAIL'}\n")
        f.write("\nMetrics:\n")
        for name, value in report["metrics"].items():
            f.write(f"  {name}: {value}\n")


def validate_airway_pydeseq2(
    current: pd.DataFrame,
    reference: pd.DataFrame,
    args: Any,
) -> dict[str, Any]:
    """Convenience wrapper for airway validation."""

    aligned = align_results(current, reference)
    metrics = compute_airway_validation_metrics(aligned)
    return evaluate_validation(metrics, args)
