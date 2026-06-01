"""Configuration models for the DEG pipeline using Pydantic."""

from __future__ import annotations

from enum import Enum
from pathlib import Path

from pydantic import BaseModel, Field, field_validator, model_validator


class Species(str, Enum):
    """Supported species for annotation."""

    HUMAN = "human"
    MOUSE = "mouse"


class GeneIdType(str, Enum):
    """Supported gene ID types."""

    ENSEMBL = "ensembl"


class GOntology(str, Enum):
    """Gene Ontology categories."""

    BP = "BP"  # Biological Process
    MF = "MF"  # Molecular Function
    CC = "CC"  # Cellular Component


class DesignConfig(BaseModel):
    """Design formula configuration."""

    formula: str = Field(..., description="Design formula, e.g., '~ condition'")
    variables: list[str] = Field(..., description="List of design variables")

    @field_validator("formula")
    @classmethod
    def validate_formula(cls, v: str) -> str:
        """Validate that formula starts with ~."""
        if not v.strip().startswith("~"):
            raise ValueError(f"Formula must start with '~': {v}")
        return v.strip()


class ContrastConfig(BaseModel):
    """Contrast configuration for differential expression."""

    name: str = Field(..., description="Contrast name for output files")
    variable: str = Field(..., description="Variable name from design")
    numerator: str = Field(..., description="Numerator level (test/treatment)")
    denominator: str = Field(..., description="Denominator level (reference/control)")


class FilteringConfig(BaseModel):
    """Low-count filtering configuration."""

    min_count: int = Field(default=10, ge=0, description="Minimum count threshold")
    min_samples: int | None = Field(
        default=None,
        description="Minimum samples with min_count. If None, uses smallest group size.",
    )


class PyDESeq2Config(BaseModel):
    """PyDESeq2 analysis configuration."""

    alpha: float = Field(default=0.05, gt=0, lt=1, description="Significance level")
    cooks_filter: bool = Field(default=True, description="Apply Cook's distance filter")
    independent_filter: bool = Field(
        default=True, description="Apply independent filtering"
    )


class DEGConfig(BaseModel):
    """DEG classification configuration."""

    padj_threshold: float = Field(
        default=0.05, gt=0, lt=1, description="Adjusted p-value threshold"
    )
    log2fc_threshold: float = Field(
        default=1.0, ge=0, description="Log2 fold change threshold"
    )


class GOConfig(BaseModel):
    """GO enrichment configuration."""

    ontology: list[GOntology] = Field(
        default=[GOntology.BP, GOntology.MF, GOntology.CC],
        description="GO ontologies to analyze",
    )
    method: str = Field(default="offline_ora", description="Enrichment method")
    padj_method: str = Field(default="fdr_bh", description="P-value adjustment method")
    min_genes_per_term: int = Field(
        default=5, ge=1, description="Minimum genes per GO term"
    )
    max_genes_per_term: int = Field(
        default=500, ge=1, description="Maximum genes per GO term"
    )

    @field_validator("method")
    @classmethod
    def validate_offline_method(cls, value: str) -> str:
        if value != "offline_ora":
            raise ValueError("v1 only supports fully offline GO ORA: method must be offline_ora")
        return value


class PipelineConfig(BaseModel):
    """Main pipeline configuration model."""

    project_name: str = Field(..., description="Project name for output directories")
    species: Species = Field(..., description="Species for annotation")
    gene_id_type: GeneIdType = Field(default=GeneIdType.ENSEMBL, description="Gene ID type")
    condition_column: str = Field(..., description="Column name for condition variable")
    reference_level: str = Field(..., description="Reference level for condition")
    results_dir: str = Field(..., description="Base results directory")

    design: DesignConfig = Field(..., description="Design formula configuration")
    contrasts: list[ContrastConfig] = Field(
        ..., min_length=1, description="List of contrasts to perform"
    )
    filtering: FilteringConfig = Field(
        default_factory=FilteringConfig, description="Filtering configuration"
    )
    pydeseq2: PyDESeq2Config = Field(
        default_factory=PyDESeq2Config, description="PyDESeq2 configuration"
    )
    deg: DEGConfig = Field(default_factory=DEGConfig, description="DEG classification")
    go: GOConfig = Field(default_factory=GOConfig, description="GO enrichment")

    @model_validator(mode="after")
    def validate_v1_scope(self) -> "PipelineConfig":
        """Validate locked v1 scope: Ensembl IDs and simple two-group contrasts."""
        if self.gene_id_type != GeneIdType.ENSEMBL:
            raise ValueError("v1 only supports Ensembl gene IDs")

        if self.design.variables != [self.condition_column]:
            raise ValueError(
                "v1 only supports a single design variable matching condition_column"
            )

        normalized_formula = self.design.formula.replace(" ", "")
        if normalized_formula != f"~{self.condition_column}":
            raise ValueError(
                "v1 only supports simple two-group formulas like '~ condition'"
            )

        for contrast in self.contrasts:
            if contrast.variable != self.condition_column:
                raise ValueError(
                    f"Contrast {contrast.name} uses {contrast.variable}; "
                    f"v1 only supports {self.condition_column}"
                )
            if contrast.denominator != self.reference_level:
                raise ValueError(
                    f"Contrast {contrast.name} denominator must match reference_level "
                    f"{self.reference_level}"
                )
            if contrast.numerator == contrast.denominator:
                raise ValueError(f"Contrast {contrast.name} has identical levels")
        return self

    def get_results_path(self) -> Path:
        """Get the full results directory path."""
        return Path(self.results_dir)

    def get_contrast(self, name: str) -> ContrastConfig | None:
        """Get a contrast by name."""
        for contrast in self.contrasts:
            if contrast.name == name:
                return contrast
        return None


def load_and_validate_config(config_path: Path) -> PipelineConfig:
    """Load and validate a configuration file.

    Args:
        config_path: Path to the YAML configuration file.

    Returns:
        Validated PipelineConfig object.

    Raises:
        FileNotFoundError: If config file doesn't exist.
        pydantic.ValidationError: If config is invalid.
    """
    import yaml

    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path) as f:
        raw_config = yaml.safe_load(f)

    return PipelineConfig(**raw_config)
