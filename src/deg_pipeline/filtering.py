"""Filtering utilities for count matrices."""

from __future__ import annotations

import logging
from typing import Any

import pandas as pd

logger = logging.getLogger(__name__)


def determine_min_samples(
    metadata_df: pd.DataFrame,
    condition_column: str,
    configured_min_samples: int | None,
) -> int:
    """Determine the minimum samples threshold for filtering.

    If configured_min_samples is provided, use it.
    Otherwise, use the size of the smallest condition group.

    Args:
        metadata_df: Metadata DataFrame with condition column.
        condition_column: Name of the condition column.
        configured_min_samples: User-configured minimum, or None.

    Returns:
        Minimum samples threshold.
    """
    if configured_min_samples is not None:
        logger.info(f"Using configured min_samples: {configured_min_samples}")
        return configured_min_samples

    # Calculate smallest group size
    group_sizes = metadata_df[condition_column].value_counts()
    min_group_size = int(group_sizes.min())
    logger.info(
        f"Smallest condition group has {min_group_size} samples, "
        f"using as min_samples threshold"
    )
    return min_group_size


def filter_low_count_genes(
    counts_df: pd.DataFrame,
    min_count: int = 10,
    min_samples: int = 3,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Filter out genes with low counts.

    Keep genes that have at least min_count reads in at least min_samples samples.

    Args:
        counts_df: Count matrix with genes as rows and samples as columns.
        min_count: Minimum count threshold per sample.
        min_samples: Minimum number of samples meeting the count threshold.

    Returns:
        Tuple of:
        - Filtered count matrix (genes that passed filter)
        - Summary DataFrame with filtering info for each gene
    """
    # Count samples meeting threshold for each gene
    samples_above_threshold = (counts_df >= min_count).sum(axis=1)

    # Create filtering mask
    keep_mask = samples_above_threshold >= min_samples

    # Create summary DataFrame
    summary = pd.DataFrame({
        "gene_id": counts_df.index,
        "samples_above_threshold": samples_above_threshold,
        "passes_filter": keep_mask,
        "mean_count": counts_df.mean(axis=1),
        "max_count": counts_df.max(axis=1),
    }).set_index("gene_id")

    # Filter counts
    filtered_counts = counts_df.loc[keep_mask]

    logger.info(
        f"Filtered {counts_df.shape[0]} genes to {filtered_counts.shape[0]} genes "
        f"(min_count={min_count}, min_samples={min_samples})"
    )

    return filtered_counts, summary


def compute_filtering_summary(
    original_counts: pd.DataFrame,
    filtered_counts: pd.DataFrame,
    min_count: int,
    min_samples: int,
) -> dict[str, Any]:
    """Compute filtering summary statistics.

    Args:
        original_counts: Original count matrix before filtering.
        filtered_counts: Count matrix after filtering.
        min_count: Count threshold used.
        min_samples: Sample threshold used.

    Returns:
        Dictionary of filtering summary statistics.
    """
    n_original = original_counts.shape[0]
    n_filtered = filtered_counts.shape[0]
    n_removed = n_original - n_filtered

    summary = {
        "n_genes_before_filtering": n_original,
        "n_genes_after_filtering": n_filtered,
        "n_genes_removed": n_removed,
        "fraction_removed": n_removed / n_original if n_original > 0 else 0,
        "min_count_threshold": min_count,
        "min_samples_threshold": min_samples,
        "total_counts_before": int(original_counts.values.sum()),
        "total_counts_after": int(filtered_counts.values.sum()),
    }

    logger.info(
        f"Filtering summary: kept {n_filtered}/{n_original} genes "
        f"({1 - summary['fraction_removed']:.1%})"
    )

    return summary