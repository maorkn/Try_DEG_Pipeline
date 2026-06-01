"""DEG (Differentially Expressed Genes) classification and output utilities."""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def classify_degs(
    results_df: pd.DataFrame,
    padj_threshold: float = 0.05,
    log2fc_threshold: float = 1.0,
) -> pd.DataFrame:
    """Classify genes as up/down/not significant DEGs.

    Args:
        results_df: PyDESeq2 results with log2FoldChange and padj columns.
        padj_threshold: Adjusted p-value threshold for significance.
        log2fc_threshold: Log2 fold change threshold.

    Returns:
        Results DataFrame with added deg_class column.
    """
    df = results_df.copy()

    # Initialize deg_class column
    df["deg_class"] = "not_significant"

    # Handle NA values in padj
    na_mask = df["padj"].isna()
    df.loc[na_mask, "deg_class"] = "not_tested"

    # Classify significant genes
    sig_mask = df["padj"] <= padj_threshold
    up_mask = sig_mask & (df["log2FoldChange"] >= log2fc_threshold)
    down_mask = sig_mask & (df["log2FoldChange"] <= -log2fc_threshold)
    low_lfc_mask = sig_mask & ~up_mask & ~down_mask & ~na_mask

    df.loc[up_mask, "deg_class"] = "up"
    df.loc[down_mask, "deg_class"] = "down"
    df.loc[low_lfc_mask, "deg_class"] = "significant_low_lfc"

    # Log summary
    class_counts = df["deg_class"].value_counts()
    logger.info(f"DEG classification summary:\n{class_counts}")

    return df


def join_annotation(
    results_df: pd.DataFrame,
    annotation_df: pd.DataFrame | None,
    gene_id_column: str = "gene_id",
) -> pd.DataFrame:
    """Join annotation to results.

    Args:
        results_df: DEG results DataFrame.
        annotation_df: Annotation DataFrame with gene_id as index.
        gene_id_column: Name of gene ID column in results.

    Returns:
        Results DataFrame with annotation columns added.
    """
    if annotation_df is None:
        logger.warning("No annotation provided, skipping join")
        return results_df

    # Reset annotation index if needed
    if annotation_df.index.name == "gene_id":
        annotation = annotation_df.reset_index()
    else:
        annotation = annotation_df.copy()

    # Left join
    merged = results_df.merge(annotation, on=gene_id_column, how="left")

    n_matched = merged[gene_id_column].isin(annotation[gene_id_column]).sum()
    logger.info(f"Joined annotation for {n_matched}/{len(results_df)} genes")

    return merged


def make_ranked_file(
    results_df: pd.DataFrame,
    gene_id_column: str = "gene_id",
    use_symbol: bool = True,
) -> pd.DataFrame:
    """Create ranked file for GSEA-style analysis.

    Ranking metric: sign(log2FC) * -log10(pvalue)
    Handle pvalue=0 by clipping to smallest positive float.

    Args:
        results_df: DEG results DataFrame.
        gene_id_column: Name of gene ID column.
        use_symbol: Use gene_symbol if available, else gene_id.

    Returns:
        Two-column DataFrame for RNK format.
    """
    df = results_df.copy()

    # Choose gene identifier
    if use_symbol and "gene_symbol" in df.columns:
        gene_col = "gene_symbol"
        # Fall back to gene_id for missing symbols
        df["rank_id"] = df[gene_col].fillna(df[gene_id_column])
    else:
        df["rank_id"] = df[gene_id_column]

    # Compute ranking metric
    pvalue = df["pvalue"].clip(lower=np.finfo(float).tiny)
    sign = np.sign(df["log2FoldChange"])
    neg_log_pvalue = -np.log10(pvalue)
    ranking_metric = sign * neg_log_pvalue

    # Handle NA values
    ranking_metric = ranking_metric.fillna(0)

    # Create output DataFrame
    ranked = pd.DataFrame({
        "gene_id": df["rank_id"],
        "ranking_metric": ranking_metric,
    })

    # Sort by ranking metric descending
    ranked = ranked.sort_values("ranking_metric", ascending=False)

    logger.info(f"Created ranked file with {len(ranked)} genes")
    return ranked


def compute_deg_summary(results_df: pd.DataFrame, contrast_name: str) -> dict[str, Any]:
    """Compute DEG summary statistics.

    Args:
        results_df: Classified DEG results.
        contrast_name: Name of the contrast.

    Returns:
        Dictionary of DEG summary statistics.
    """
    summary = {
        "contrast": contrast_name,
        "total_genes_tested": len(results_df),
        "genes_not_tested": int((results_df["deg_class"] == "not_tested").sum()),
        "genes_up": int((results_df["deg_class"] == "up").sum()),
        "genes_down": int((results_df["deg_class"] == "down").sum()),
        "genes_significant_low_lfc": int(
            (results_df["deg_class"] == "significant_low_lfc").sum()
        ),
        "genes_not_significant": int(
            (results_df["deg_class"] == "not_significant").sum()
        ),
    }

    # Total significant
    summary["total_significant"] = (
        summary["genes_up"] + summary["genes_down"] + summary["genes_significant_low_lfc"]
    )

    logger.info(
        f"DEG summary for {contrast_name}: "
        f"{summary['genes_up']} up, {summary['genes_down']} down, "
        f"{summary['total_significant']} total significant"
    )

    return summary