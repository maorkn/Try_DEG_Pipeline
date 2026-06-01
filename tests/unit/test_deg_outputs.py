"""Unit tests for DEG output helpers."""

from __future__ import annotations

import pandas as pd

from deg_pipeline.deg import classify_degs, compute_deg_summary, join_annotation, make_ranked_file


def test_classify_degs_handles_up_down_low_lfc_and_na():
    results = pd.DataFrame(
        {
            "gene_id": ["up", "down", "low_lfc", "na", "ns"],
            "log2FoldChange": [2.0, -2.0, 0.2, 3.0, 0.1],
            "pvalue": [0.001, 0.001, 0.01, 0.01, 0.5],
            "padj": [0.01, 0.02, 0.03, None, 0.8],
        }
    )

    classified = classify_degs(results, padj_threshold=0.05, log2fc_threshold=1.0)

    assert classified.set_index("gene_id").loc["up", "deg_class"] == "up"
    assert classified.set_index("gene_id").loc["down", "deg_class"] == "down"
    assert (
        classified.set_index("gene_id").loc["low_lfc", "deg_class"]
        == "significant_low_lfc"
    )
    assert classified.set_index("gene_id").loc["na", "deg_class"] == "not_tested"
    assert classified.set_index("gene_id").loc["ns", "deg_class"] == "not_significant"


def test_join_annotation_preserves_all_result_rows(toy_annotation_df):
    results = pd.DataFrame({"gene_id": ["ENSG000001", "missing"], "padj": [0.1, 0.2]})

    merged = join_annotation(results, toy_annotation_df)

    assert len(merged) == 2
    assert "gene_symbol" in merged.columns
    assert merged.loc[merged["gene_id"] == "missing", "gene_symbol"].isna().all()


def test_make_ranked_file_clips_zero_pvalues():
    results = pd.DataFrame(
        {
            "gene_id": ["g1", "g2"],
            "log2FoldChange": [1.0, -1.0],
            "pvalue": [0.0, 0.01],
        }
    )

    ranked = make_ranked_file(results, use_symbol=False)

    assert list(ranked.columns) == ["gene_id", "ranking_metric"]
    assert ranked["ranking_metric"].notna().all()
    assert ranked.iloc[0]["gene_id"] == "g1"


def test_compute_deg_summary_counts_classes():
    classified = pd.DataFrame(
        {"deg_class": ["up", "down", "not_significant", "significant_low_lfc", "not_tested"]}
    )

    summary = compute_deg_summary(classified, "treated_vs_control")

    assert summary["genes_up"] == 1
    assert summary["genes_down"] == 1
    assert summary["total_significant"] == 3
