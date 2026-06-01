"""GO enrichment analysis utilities for the DEG pipeline."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

logger = logging.getLogger(__name__)


def build_gene_universe(
    filtered_gene_ids: list[str] | set[str],
    gene2go: pd.DataFrame | dict[str, set[str]],
) -> set[str]:
    """Build the gene universe for GO enrichment analysis.

    The gene universe is the set of genes that have valid GO annotations
    and can be used as background for enrichment analysis.

    Args:
        filtered_gene_ids: Gene IDs that passed filtering.
        gene2go: Mapping from gene IDs to GO term IDs.

    Returns:
        Set of gene IDs that form the universe.
    """
    if isinstance(gene2go, dict):
        # Convert dict format to set
        filtered = set(filtered_gene_ids)
        universe = set()
        for gene_id, go_terms in gene2go.items():
            if gene_id in filtered and go_terms:  # Only include filtered genes with GO terms
                universe.add(gene_id)
        return universe

    # DataFrame format
    universe = set(filtered_gene_ids)

    # Intersect with GO-annotated genes
    if isinstance(gene2go, pd.DataFrame):
        go_annotated = set(gene2go["gene_id"].unique())
        universe = universe & go_annotated

    logger.info(f"Built gene universe: {len(universe)} genes")
    return universe


def load_gene2go_table(path: Path) -> pd.DataFrame:
    """Load a local offline gene-to-GO mapping table.

    Supported schemas:
    - `gene_id`, `go_id`, optional `ontology`
    - `gene_id`, `go_ids` where GO IDs are semicolon-separated
    """

    if not path.exists():
        raise FileNotFoundError(f"gene2go mapping not found: {path}")

    df = pd.read_csv(path, sep="\t")
    if "gene_id" not in df.columns:
        raise ValueError("gene2go mapping must contain a gene_id column")

    if "go_id" in df.columns:
        mapping = df.copy()
    elif "go_ids" in df.columns:
        mapping = (
            df.assign(go_id=df["go_ids"].fillna("").astype(str).str.split(";"))
            .explode("go_id")
            .drop(columns=["go_ids"])
        )
    else:
        raise ValueError("gene2go mapping must contain either go_id or go_ids")

    mapping["go_id"] = mapping["go_id"].astype(str).str.strip()
    mapping = mapping[mapping["go_id"] != ""].copy()
    if "ontology" not in mapping.columns:
        mapping["ontology"] = "NA"
    return mapping


def load_deg_gene_sets(
    deg_dir: Path,
    contrast: str,
) -> dict[str, set[str]]:
    """Load DEG gene sets from output directory.

    Args:
        deg_dir: Directory containing DEG output files.
        contrast: Contrast name to load.

    Returns:
        Dictionary mapping direction ('up'/'down') to gene sets.
    """
    gene_sets = {}

    # Load up-regulated genes
    up_path = deg_dir / f"{contrast}_up.tsv"
    if up_path.exists():
        up_genes = set(pd.read_csv(up_path, sep="\t")["gene_id"].tolist())
        gene_sets["up"] = up_genes

    # Load down-regulated genes
    down_path = deg_dir / f"{contrast}_down.tsv"
    if down_path.exists():
        down_genes = set(pd.read_csv(down_path, sep="\t")["gene_id"].tolist())
        gene_sets["down"] = down_genes

    return gene_sets


def run_go_ora(
    gene_set: set[str],
    universe: set[str],
    ontology: str,
    go_obo: Path | None,
    gene2go: Path | pd.DataFrame,
    padj_method: str = "fdr_bh",
    min_genes_per_term: int = 5,
    max_genes_per_term: int = 500,
) -> pd.DataFrame:
    """Run GO over-representation analysis.

    Args:
        gene_set: Set of genes to test for enrichment.
        universe: Background gene universe.
        ontology: GO ontology to use ('BP', 'MF', or 'CC').
        go_obo: Path to go-basic.obo file.
        gene2go: Path to gene2go mapping file.

    Returns:
        DataFrame with GO enrichment results.
    """
    logger.info(f"Running GO ORA for {ontology} with {len(gene_set)} genes")

    mapping = load_gene2go_table(gene2go) if isinstance(gene2go, Path) else gene2go.copy()
    if "ontology" in mapping.columns and ontology != "NA":
        ont_mapping = mapping[
            (mapping["ontology"].astype(str).str.upper() == ontology.upper())
            | (mapping["ontology"].astype(str) == "NA")
        ].copy()
    else:
        ont_mapping = mapping

    universe = set(universe)
    study = set(gene_set) & universe
    ont_mapping = ont_mapping[ont_mapping["gene_id"].isin(universe)]

    population_n = len(universe)
    study_n = len(study)
    if population_n == 0 or study_n == 0 or ont_mapping.empty:
        return pd.DataFrame(
            columns=[
                "go_id",
                "ontology",
                "study_count",
                "study_n",
                "population_count",
                "population_n",
                "pvalue",
                "padj",
                "study_genes",
            ]
        )

    rows: list[dict[str, Any]] = []
    for go_id, group in ont_mapping.groupby("go_id"):
        term_genes = set(group["gene_id"])
        population_count = len(term_genes)
        if population_count < min_genes_per_term or population_count > max_genes_per_term:
            continue
        overlap = study & term_genes
        if not overlap:
            continue
        pvalue = hypergeom.sf(len(overlap) - 1, population_n, population_count, study_n)
        rows.append(
            {
                "go_id": go_id,
                "ontology": ontology,
                "study_count": len(overlap),
                "study_n": study_n,
                "population_count": population_count,
                "population_n": population_n,
                "pvalue": float(pvalue),
                "study_genes": ";".join(sorted(overlap)),
            }
        )

    results = pd.DataFrame(rows)
    if results.empty:
        results["padj"] = []
        return results

    results = results.sort_values("pvalue").reset_index(drop=True)
    if padj_method.lower() != "fdr_bh":
        raise ValueError(f"Unsupported GO p-value adjustment method: {padj_method}")
    ranks = np.arange(1, len(results) + 1)
    adjusted = results["pvalue"].to_numpy() * len(results) / ranks
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    results["padj"] = np.clip(adjusted, 0, 1)
    return results.sort_values(["padj", "pvalue"]).reset_index(drop=True)
