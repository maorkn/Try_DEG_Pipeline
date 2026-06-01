"""Unit tests for offline GO enrichment helpers."""

from __future__ import annotations

import pandas as pd

from deg_pipeline.enrichment import build_gene_universe, load_gene2go_table, run_go_ora


def test_build_gene_universe_excludes_unmapped_genes(toy_go_map_path):
    gene2go = load_gene2go_table(toy_go_map_path)

    universe = build_gene_universe({"ENSG000001", "ENSG000002", "missing"}, gene2go)

    assert universe == {"ENSG000001", "ENSG000002"}


def test_load_gene2go_table_expands_semicolon_terms(toy_go_map_path):
    gene2go = load_gene2go_table(toy_go_map_path)

    rows = gene2go[gene2go["gene_id"] == "ENSG000001"]
    assert set(rows["go_id"]) == {"GO:0000001", "GO:0000002"}


def test_run_go_ora_returns_expected_columns():
    gene2go = pd.DataFrame(
        {
            "gene_id": ["g1", "g2", "g3", "g4"],
            "go_id": ["GO:1", "GO:1", "GO:2", "GO:2"],
            "ontology": ["BP", "BP", "BP", "BP"],
        }
    )

    results = run_go_ora(
        {"g1", "g2"},
        {"g1", "g2", "g3", "g4"},
        "BP",
        None,
        gene2go,
        min_genes_per_term=1,
        max_genes_per_term=10,
    )

    assert {"go_id", "pvalue", "padj", "study_genes"}.issubset(results.columns)
    assert "GO:1" in set(results["go_id"])
