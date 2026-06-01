"""Plotting utilities for the DEG pipeline."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

logger = logging.getLogger(__name__)

# Set default style
sns.set_theme(style="whitegrid")
plt.rcParams["figure.dpi"] = 150


def plot_library_sizes(
    library_metrics: pd.DataFrame,
    output_path: Path,
    title: str = "Library Sizes per Sample",
) -> Path:
    """Plot library sizes as a bar chart.

    Args:
        library_metrics: DataFrame with library metrics.
        output_path: Path to save the plot.
        title: Plot title.

    Returns:
        Path to saved plot.
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    samples = library_metrics.index
    total_counts = library_metrics["total_counts"]
    detected_genes = library_metrics["detected_genes"]

    # Create twin axes
    ax2 = ax.twinx()

    # Bar plot for total counts
    bars = ax.bar(samples, total_counts, alpha=0.7, color="steelblue", label="Total counts")

    # Line plot for detected genes
    line = ax2.plot(
        samples, detected_genes, "ro-", label="Detected genes", linewidth=2, markersize=8
    )

    ax.set_xlabel("Sample")
    ax.set_ylabel("Total Counts", color="steelblue")
    ax2.set_ylabel("Detected Genes", color="red")

    ax.tick_params(axis="x", rotation=45)
    ax.tick_params(axis="y", labelcolor="steelblue")
    ax2.tick_params(axis="y", labelcolor="red")

    # Combined legend
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc="upper left")

    ax.set_title(title)
    plt.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, bbox_inches="tight")
    plt.close()

    logger.info(f"Saved library sizes plot to {output_path}")
    return output_path


def plot_detected_genes(
    gene_metrics: pd.DataFrame,
    output_path: Path,
    title: str = "Gene Detection Distribution",
) -> Path:
    """Plot distribution of genes by detection frequency.

    Args:
        gene_metrics: DataFrame with gene metrics.
        output_path: Path to save the plot.
        title: Plot title.

    Returns:
        Path to saved plot.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Histogram of detection frequency
    n_samples = gene_metrics["detected_samples"].max()
    axes[0].hist(
        gene_metrics["detected_samples"],
        bins=range(int(n_samples) + 2),
        edgecolor="black",
        alpha=0.7,
    )
    axes[0].set_xlabel("Number of Samples with Count > 0")
    axes[0].set_ylabel("Number of Genes")
    axes[0].set_title("Gene Detection Frequency")

    # Histogram of mean counts (log scale)
    nonzero_means = gene_metrics.loc[gene_metrics["mean_counts"] > 0, "mean_counts"]
    if len(nonzero_means) > 0:
        axes[1].hist(
            np.log10(nonzero_means),
            bins=50,
            edgecolor="black",
            alpha=0.7,
        )
        axes[1].set_xlabel("log10(Mean Count)")
        axes[1].set_ylabel("Number of Genes")
        axes[1].set_title("Mean Count Distribution")
    else:
        axes[1].text(
            0.5, 0.5, "No genes with non-zero counts",
            ha="center", va="center", transform=axes[1].transAxes
        )

    plt.suptitle(title)
    plt.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, bbox_inches="tight")
    plt.close()

    logger.info(f"Saved detected genes plot to {output_path}")
    return output_path


def plot_sample_correlation_heatmap(
    correlation_df: pd.DataFrame,
    output_path: Path,
    title: str = "Sample Correlation Heatmap",
) -> Path:
    """Plot sample correlation as a heatmap.

    Args:
        correlation_df: Correlation matrix.
        output_path: Path to save the plot.
        title: Plot title.

    Returns:
        Path to saved plot.
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    sns.heatmap(
        correlation_df,
        annot=True,
        fmt=".2f",
        cmap="RdYlBu_r",
        vmin=0.8,
        vmax=1.0,
        square=True,
        ax=ax,
    )

    ax.set_title(title)
    plt.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, bbox_inches="tight")
    plt.close()

    logger.info(f"Saved correlation heatmap to {output_path}")
    return output_path


def plot_pca(
    pca_df: pd.DataFrame,
    output_path: Path,
    color_by: str | None = None,
    title: str = "PCA Plot",
) -> Path:
    """Plot PCA coordinates.

    Args:
        pca_df: DataFrame with PC1, PC2 columns and optional color variable.
        output_path: Path to save the plot.
        color_by: Column name to use for coloring points.
        title: Plot title.

    Returns:
        Path to saved plot.
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    variance = pca_df.attrs.get("variance_explained", [0, 0])
    x_label = f"PC1 ({variance[0]*100:.1f}% variance)"
    y_label = f"PC2 ({variance[1]*100:.1f}% variance)"

    if color_by and color_by in pca_df.columns:
        # Color by categorical variable
        categories = pca_df[color_by].unique()
        colors = sns.color_palette("husl", len(categories))

        for cat, color in zip(categories, colors):
            mask = pca_df[color_by] == cat
            ax.scatter(
                pca_df.loc[mask, "PC1"],
                pca_df.loc[mask, "PC2"],
                c=[color],
                label=cat,
                s=100,
                alpha=0.8,
            )

            # Add sample labels
            for idx in pca_df.loc[mask].index:
                ax.annotate(
                    idx,
                    (pca_df.loc[idx, "PC1"], pca_df.loc[idx, "PC2"]),
                    fontsize=8,
                    alpha=0.7,
                )

        ax.legend(title=color_by)
    else:
        # No coloring
        ax.scatter(pca_df["PC1"], pca_df["PC2"], s=100, alpha=0.8)
        for idx in pca_df.index:
            ax.annotate(
                idx,
                (pca_df.loc[idx, "PC1"], pca_df.loc[idx, "PC2"]),
                fontsize=8,
            )

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.axhline(y=0, color="gray", linestyle="--", alpha=0.5)
    ax.axvline(x=0, color="gray", linestyle="--", alpha=0.5)

    plt.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, bbox_inches="tight")
    plt.close()

    logger.info(f"Saved PCA plot to {output_path}")
    return output_path


def plot_volcano(
    results_df: pd.DataFrame,
    output_path: Path,
    padj_threshold: float = 0.05,
    log2fc_threshold: float = 1.0,
    title: str = "Volcano Plot",
) -> Path:
    """Plot volcano plot for DEG results.

    Args:
        results_df: DEG results with log2FoldChange and padj columns.
        output_path: Path to save the plot.
        padj_threshold: Adjusted p-value threshold for significance.
        log2fc_threshold: Log2 fold change threshold.
        title: Plot title.

    Returns:
        Path to saved plot.
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    # Filter out NA values
    df = results_df.dropna(subset=["log2FoldChange", "padj"])

    # Compute -log10(padj), handling zeros
    neg_log_padj = -np.log10(df["padj"].clip(lower=1e-300))
    log2fc = df["log2FoldChange"]

    # Classify genes
    sig_up = (df["padj"] <= padj_threshold) & (log2fc >= log2fc_threshold)
    sig_down = (df["padj"] <= padj_threshold) & (log2fc <= -log2fc_threshold)
    not_sig = ~(sig_up | sig_down)

    # Plot
    ax.scatter(log2fc[not_sig], neg_log_padj[not_sig], c="gray", alpha=0.5, s=10, label="Not significant")
    ax.scatter(log2fc[sig_up], neg_log_padj[sig_up], c="red", alpha=0.7, s=20, label=f"Up (n={sig_up.sum()})")
    ax.scatter(log2fc[sig_down], neg_log_padj[sig_down], c="blue", alpha=0.7, s=20, label=f"Down (n={sig_down.sum()})")

    # Add threshold lines
    ax.axhline(-np.log10(padj_threshold), color="black", linestyle="--", alpha=0.5)
    ax.axvline(log2fc_threshold, color="black", linestyle="--", alpha=0.5)
    ax.axvline(-log2fc_threshold, color="black", linestyle="--", alpha=0.5)

    ax.set_xlabel("log2 Fold Change")
    ax.set_ylabel("-log10(adjusted p-value)")
    ax.set_title(title)
    ax.legend()

    plt.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, bbox_inches="tight")
    plt.close()

    logger.info(f"Saved volcano plot to {output_path}")
    return output_path


def plot_ma(
    results_df: pd.DataFrame,
    output_path: Path,
    padj_threshold: float = 0.05,
    title: str = "MA Plot",
) -> Path:
    """Plot MA plot for DEG results.

    Args:
        results_df: DEG results with baseMean and log2FoldChange columns.
        output_path: Path to save the plot.
        padj_threshold: Adjusted p-value threshold for significance.
        title: Plot title.

    Returns:
        Path to saved plot.
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    # Filter out NA values
    df = results_df.dropna(subset=["baseMean", "log2FoldChange"])

    # Log transform baseMean
    log_mean = np.log10(df["baseMean"].clip(lower=1))
    log2fc = df["log2FoldChange"]

    # Classify genes
    sig = df["padj"] <= padj_threshold
    not_sig = ~sig

    # Plot
    ax.scatter(log_mean[not_sig], log2fc[not_sig], c="gray", alpha=0.3, s=10, label="Not significant")
    ax.scatter(log_mean[sig], log2fc[sig], c="red", alpha=0.5, s=15, label=f"Significant (n={sig.sum()})")

    ax.axhline(0, color="black", linestyle="--", alpha=0.5)

    ax.set_xlabel("log10(baseMean)")
    ax.set_ylabel("log2 Fold Change")
    ax.set_title(title)
    ax.legend()

    plt.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, bbox_inches="tight")
    plt.close()

    logger.info(f"Saved MA plot to {output_path}")
    return output_path