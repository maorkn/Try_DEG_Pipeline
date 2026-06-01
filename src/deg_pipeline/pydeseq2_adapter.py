"""PyDESeq2 adapter for the DEG pipeline."""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def prepare_pydeseq2_inputs(
    counts_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    condition_column: str,
    reference_level: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Prepare inputs for PyDESeq2.

    PyDESeq2 expects:
    - counts: DataFrame with samples as rows and genes as columns
    - metadata: DataFrame with samples as index

    Args:
        counts_df: Count matrix with genes as rows and samples as columns.
        metadata_df: Metadata with samples as index.
        condition_column: Name of the condition column.
        reference_level: Reference level for the condition.

    Returns:
        Tuple of (counts_for_pydeseq2, metadata_for_pydeseq2).

    Raises:
        ValueError: If condition doesn't have exactly two levels.
    """
    if condition_column not in metadata_df.columns:
        raise ValueError(f"Metadata missing condition column: {condition_column}")

    missing_metadata = [sample for sample in counts_df.columns if sample not in metadata_df.index]
    if missing_metadata:
        raise ValueError(f"Samples in counts missing from metadata: {missing_metadata}")

    extra_metadata = [sample for sample in metadata_df.index if sample not in counts_df.columns]
    if extra_metadata:
        logger.warning("Ignoring metadata samples not present in counts: %s", extra_metadata)

    metadata_df = metadata_df.loc[list(counts_df.columns)].copy()

    # Check condition levels
    condition_levels = metadata_df[condition_column].unique()
    if len(condition_levels) != 2:
        raise ValueError(
            f"Condition '{condition_column}' must have exactly two levels for v1, "
            f"found {len(condition_levels)}: {list(condition_levels)}"
        )

    # Check reference level exists
    if reference_level not in condition_levels:
        raise ValueError(
            f"Reference level '{reference_level}' not found in condition levels: "
            f"{list(condition_levels)}"
        )

    # Transpose counts: samples as rows, genes as columns
    counts_for_pydeseq2 = counts_df.T
    metadata_for_pydeseq2 = metadata_df.loc[counts_for_pydeseq2.index]

    # Ensure counts are integers
    counts_for_pydeseq2 = counts_for_pydeseq2.astype(int)

    logger.info(
        f"Prepared PyDESeq2 inputs: {counts_for_pydeseq2.shape[0]} samples, "
        f"{counts_for_pydeseq2.shape[1]} genes"
    )

    return counts_for_pydeseq2, metadata_for_pydeseq2


def run_pydeseq2_model(
    counts_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    condition_column: str,
    reference_level: str,
    cooks_filter: bool = True,
) -> Any:
    """Run PyDESeq2 differential expression model.

    Args:
        counts_df: Counts with samples as rows and genes as columns.
        metadata_df: Metadata with samples as index.
        condition_column: Name of condition column.
        reference_level: Reference level for condition.
        cooks_filter: Whether to apply Cook's distance filter.

    Returns:
        DeseqDataSet object after running deseq2().
    """
    try:
        from pydeseq2.dds import DeseqDataSet
    except ImportError as exc:
        raise RuntimeError(
            "PyDESeq2 is not installed. Install the project environment before running "
            "scripts/04_run_pydeseq2.py."
        ) from exc

    logger.info(
        f"Running PyDESeq2 with condition={condition_column}, "
        f"reference={reference_level}, cooks_filter={cooks_filter}"
    )

    try:
        dds = DeseqDataSet(
            counts=counts_df,
            metadata=metadata_df,
            design=f"~{condition_column}",
            ref_level=[condition_column, reference_level],
            refit_cooks=cooks_filter,
        )
    except TypeError:
        # Backward-compatible fallback for older PyDESeq2 releases.
        dds = DeseqDataSet(
            counts=counts_df,
            metadata=metadata_df,
            design_factors=condition_column,
            ref_level=[condition_column, reference_level],
            refit=cooks_filter,
        )

    dds.deseq2()

    logger.info("PyDESeq2 model completed")
    return dds


def run_pydeseq2_stats(
    dds: Any,
    condition_column: str,
    numerator: str,
    denominator: str,
    alpha: float = 0.05,
    cooks_filter: bool = True,
) -> pd.DataFrame:
    """Run PyDESeq2 statistical tests for a contrast.

    Args:
        dds: DeseqDataSet after running deseq2().
        condition_column: Name of condition column.
        numerator: Numerator level (test/treatment).
        denominator: Denominator level (reference/control).
        alpha: Significance level.
        cooks_filter: Whether to apply Cook's distance filter.

    Returns:
        Results DataFrame with DESeq2-compatible columns.
    """
    try:
        from pydeseq2.ds import DeseqStats
    except ImportError as exc:
        raise RuntimeError(
            "PyDESeq2 is not installed. Install the project environment before running "
            "scripts/04_run_pydeseq2.py."
        ) from exc

    logger.info(
        f"Running PyDESeq2 stats: {condition_column} {numerator} vs {denominator}"
    )

    try:
        stat_res = DeseqStats(
            dds,
            contrast=[condition_column, numerator, denominator],
            alpha=alpha,
            cooks_filter=cooks_filter,
        )
    except TypeError:
        stat_res = DeseqStats(
            dds,
            contrast=[condition_column, numerator, denominator],
            alpha=alpha,
        )

    stat_res.summary()

    # Get results DataFrame
    results_df = stat_res.results_df.copy()

    # Ensure DESeq2-compatible column names
    column_mapping = {
        "baseMean": "baseMean",
        "log2FoldChange": "log2FoldChange",
        "lfcSE": "lfcSE",
        "stat": "stat",
        "pvalue": "pvalue",
        "padj": "padj",
    }

    # Rename columns if needed (PyDESeq2 might use different names)
    for old_name, new_name in column_mapping.items():
        if old_name in results_df.columns and old_name != new_name:
            results_df = results_df.rename(columns={old_name: new_name})

    # Add gene_id as a column
    results_df = results_df.reset_index()
    if "index" in results_df.columns:
        results_df = results_df.rename(columns={"index": "gene_id"})

    logger.info(
        f"PyDESeq2 stats: {len(results_df)} genes tested, "
        f"{(results_df['padj'] < alpha).sum()} significant at padj < {alpha}"
    )

    return results_df


def extract_pydeseq2_outputs(dds: Any) -> dict[str, pd.DataFrame | np.ndarray]:
    """Extract various outputs from PyDESeq2 model.

    Args:
        dds: DeseqDataSet after running deseq2().

    Returns:
        Dictionary of outputs (size_factors, normalized_counts, etc.).
    """
    outputs = {}

    # Size factors
    if hasattr(dds, "size_factors") and dds.size_factors is not None:
        outputs["size_factors"] = pd.DataFrame(
            {"sample_id": dds.obs_names, "size_factor": dds.size_factors}
        ).set_index("sample_id")

    # Normalized counts (if available)
    if hasattr(dds, "normed_counts") and dds.normed_counts is not None:
        normed = pd.DataFrame(
            dds.normed_counts,
            index=dds.obs_names,
            columns=dds.var_names,
        )
        outputs["normalized_counts"] = normed.T  # Genes as rows

    # Dispersion estimates (if available)
    if hasattr(dds, "dispersions") and dds.dispersions is not None:
        outputs["dispersion_estimates"] = pd.DataFrame(
            {"gene_id": dds.var_names, "dispersion": dds.dispersions}
        ).set_index("gene_id")

    logger.info(f"Extracted {len(outputs)} PyDESeq2 outputs")
    return outputs
