"""Input/output utilities for the DEG pipeline."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

logger = logging.getLogger(__name__)


def load_counts(counts_path: Path) -> pd.DataFrame:
    """Load a gene count matrix from a TSV file.
    
    Expected format:
    - Tab-separated file
    - First column is gene_id
    - Remaining columns are sample counts
    - Row index will be set to gene_id
    
    Args:
        counts_path: Path to the counts TSV file.
        
    Returns:
        DataFrame with gene_id as index and sample names as columns.
        
    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file format is invalid.
    """
    if not counts_path.exists():
        raise FileNotFoundError(f"Counts file not found: {counts_path}")
    
    df = pd.read_csv(counts_path, sep="\t", index_col=0)
    
    # Validate basic structure
    if df.empty:
        raise ValueError(f"Counts file is empty: {counts_path}")
    
    if df.index.has_duplicates:
        raise ValueError(f"Duplicate gene IDs found in counts file: {counts_path}")
    
    logger.info(
        f"Loaded counts: {df.shape[0]} genes x {df.shape[1]} samples from {counts_path}"
    )
    return df


def load_metadata(metadata_path: Path) -> pd.DataFrame:
    """Load sample metadata from a TSV file.
    
    Expected format:
    - Tab-separated file
    - Must contain 'sample_id' column
    - Row index will be set to sample_id
    
    Args:
        metadata_path: Path to the metadata TSV file.
        
    Returns:
        DataFrame with sample_id as index.
        
    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If required columns are missing.
    """
    if not metadata_path.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_path}")
    
    df = pd.read_csv(metadata_path, sep="\t")
    
    if "sample_id" not in df.columns:
        raise ValueError(
            f"Metadata file missing required 'sample_id' column: {metadata_path}"
        )
    
    df = df.set_index("sample_id")
    
    if df.index.has_duplicates:
        raise ValueError(f"Duplicate sample IDs found in metadata file: {metadata_path}")
    
    logger.info(f"Loaded metadata: {df.shape[0]} samples from {metadata_path}")
    return df


def load_annotation(annotation_path: Path) -> pd.DataFrame:
    """Load gene annotation from a TSV file.
    
    Expected format:
    - Tab-separated file
    - Must contain 'gene_id' column
    - Optional columns: gene_symbol, entrez_id, gene_name, biotype
    
    Args:
        annotation_path: Path to the annotation TSV file.
        
    Returns:
        DataFrame with gene_id as index.
        
    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If required columns are missing.
    """
    if not annotation_path.exists():
        raise FileNotFoundError(f"Annotation file not found: {annotation_path}")
    
    df = pd.read_csv(annotation_path, sep="\t")
    
    if "gene_id" not in df.columns:
        raise ValueError(
            f"Annotation file missing required 'gene_id' column: {annotation_path}"
        )
    
    df = df.set_index("gene_id")
    
    logger.info(f"Loaded annotation: {df.shape[0]} genes from {annotation_path}")
    return df


def load_config(config_path: Path) -> dict[str, Any]:
    """Load analysis configuration from a YAML file.
    
    Args:
        config_path: Path to the YAML configuration file.
        
    Returns:
        Configuration dictionary.
        
    Raises:
        FileNotFoundError: If the file does not exist.
        yaml.YAMLError: If the YAML is invalid.
    """
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    
    with open(config_path) as f:
        config = yaml.safe_load(f)
    
    if config is None:
        raise ValueError(f"Config file is empty: {config_path}")
    
    logger.info(f"Loaded config from {config_path}")
    return config


def save_dataframe(
    df: pd.DataFrame,
    output_path: Path,
    index_label: str | None = None,
    include_index: bool = True,
) -> Path:
    """Save a DataFrame to a TSV file.
    
    Args:
        df: DataFrame to save.
        output_path: Output file path.
        index_label: Label for the index column (uses index name if None).
        include_index: Whether to include the index as a column.
        
    Returns:
        Path to the saved file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=include_index, index_label=index_label)
    logger.info(f"Saved {df.shape[0]} rows to {output_path}")
    return output_path


def save_json(data: dict[str, Any], output_path: Path) -> Path:
    """Save a dictionary to a JSON file.
    
    Args:
        data: Dictionary to save.
        output_path: Output file path.
        
    Returns:
        Path to the saved file.
    """
    import json
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2, default=str)
    logger.info(f"Saved JSON to {output_path}")
    return output_path