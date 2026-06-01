"""Quality-control metrics for count matrices."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd


def compute_library_metrics(counts_df: pd.DataFrame) -> pd.DataFrame:
    """Compute per-sample library metrics from raw counts."""

    total_counts = counts_df.sum(axis=0)
    detected_genes = (counts_df > 0).sum(axis=0)
    zero_fraction = (counts_df == 0).sum(axis=0) / counts_df.shape[0]
    upper_quartile = counts_df.quantile(0.75, axis=0)

    metrics = pd.DataFrame(
        {
            "total_counts": total_counts.astype(int),
            "detected_genes": detected_genes.astype(int),
            "zero_fraction": zero_fraction.astype(float),
            "upper_quartile": upper_quartile.astype(float),
        }
    )
    metrics.index.name = "sample_id"
    return metrics


def compute_gene_metrics(counts_df: pd.DataFrame) -> pd.DataFrame:
    """Compute per-gene detection metrics from raw counts."""

    metrics = pd.DataFrame(
        {
            "total_counts": counts_df.sum(axis=1).astype(int),
            "mean_counts": counts_df.mean(axis=1).astype(float),
            "detected_samples": (counts_df > 0).sum(axis=1).astype(int),
        }
    )
    metrics.index.name = "gene_id"
    return metrics


def compute_sample_correlation(counts_df: pd.DataFrame) -> pd.DataFrame:
    """Compute Spearman sample correlation on log2(count + 1) values."""

    transformed = np.log2(counts_df.astype(float) + 1.0)
    correlation = transformed.corr(method="spearman")
    correlation.index.name = "sample_id"
    return correlation


def compute_pca(
    counts_df: pd.DataFrame,
    metadata_df: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, np.ndarray]:
    """Compute two principal components from log2(count + 1) values.

    Uses NumPy SVD so the QC stage does not fail when scikit-learn is not installed.
    """

    transformed = np.log2(counts_df.astype(float).T + 1.0)
    values = transformed.to_numpy()
    values = values - values.mean(axis=0, keepdims=True)
    std = values.std(axis=0, ddof=1)
    std[std == 0] = 1.0
    values = values / std

    if values.shape[0] < 2 or values.shape[1] < 2:
        coords = np.zeros((values.shape[0], 2))
        variance = np.array([0.0, 0.0])
    else:
        u, singular_values, _ = np.linalg.svd(values, full_matrices=False)
        n_components = min(2, u.shape[1])
        coords = np.zeros((values.shape[0], 2))
        coords[:, :n_components] = u[:, :n_components] * singular_values[:n_components]

        explained = singular_values**2
        total = explained.sum()
        variance = np.zeros(2)
        if total > 0:
            variance[:n_components] = explained[:n_components] / total

    pca_df = pd.DataFrame(coords, index=counts_df.columns, columns=["PC1", "PC2"])
    pca_df.index.name = "sample_id"
    pca_df.attrs["variance_explained"] = variance

    if metadata_df is not None:
        common_cols = [col for col in metadata_df.columns if col not in pca_df.columns]
        pca_df = pca_df.join(metadata_df[common_cols], how="left")

    return pca_df, variance


def compute_qc_summary(
    counts_df: pd.DataFrame,
    library_metrics: pd.DataFrame,
    gene_metrics: pd.DataFrame,
) -> dict[str, Any]:
    """Summarize QC metrics in a JSON-friendly dictionary."""

    return {
        "n_samples": int(counts_df.shape[1]),
        "n_genes": int(counts_df.shape[0]),
        "total_counts": int(counts_df.to_numpy().sum()),
        "median_library_size": float(library_metrics["total_counts"].median()),
        "median_detected_genes": float(library_metrics["detected_genes"].median()),
        "min_library_size": int(library_metrics["total_counts"].min()),
        "max_library_size": int(library_metrics["total_counts"].max()),
        "genes_not_detected": int((gene_metrics["detected_samples"] == 0).sum()),
        "median_gene_count": float(gene_metrics["mean_counts"].median()),
    }
