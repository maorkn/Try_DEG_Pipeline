# Reproducible RNA-seq DEG Pipeline Implementation Plan

Status: draft for review  
Repository: `Try_DEG_Pipeline`  
Primary validation dataset: Bioconductor `airway` count dataset  
Primary statistical engine: Python-native `PyDESeq2`
Reference validation engine: R/Bioconductor `DESeq2` on the public `airway` dataset

## 1. Goal

Build a reproducible RNA-seq analysis pipeline that starts from gene-level raw count tables and sample metadata, performs input validation and QC, runs PyDESeq2 differential expression analysis, validates behavior against a known public dataset and R/Bioconductor DESeq2 reference behavior, produces DEG tables, runs offline GO enrichment, and writes a compact report bundle.

The implementation should be easy for coding agents to split into independent tasks. Each pipeline unit should be a standalone Python script with:

- A command-line interface.
- A small pure-Python core that can be unit tested.
- One unit test file.
- One smoke test that executes the script on a tiny fixture or the reduced `airway` validation fixture.
- Deterministic output paths and machine-readable logs.

## 2. Design Principles

1. Prefer explicit files over hidden state.
2. Treat the pipeline as file-in/file-out modules, not one large notebook.
3. Keep each script executable on its own.
4. Validate early, fail with actionable messages, and write validation reports.
5. Use `PyDESeq2` as the production statistical engine while validating against R/Bioconductor `DESeq2` reference outputs.
6. Use Python for orchestration, validation, QC, table handling, reporting, and enrichment.
7. Pin dependencies and record session information for every run.
8. Use `airway` as a public validation workflow, not as an assumption baked into production analysis.
9. Make tests small and fast by using toy fixtures; reserve full `airway` runs for smoke/integration validation.
10. Scope version 1 to simple two-group contrasts for human and mouse Ensembl gene IDs.

## 3. Locked Version 1 Decisions

These decisions were confirmed during planning:

1. Primary statistical engine: `PyDESeq2`.
2. Validation comparison: R/Bioconductor `DESeq2` reference output on `airway`.
3. Species scope: human and mouse.
4. Gene ID standard: Ensembl gene IDs.
5. GO enrichment: fully offline and pinned.
6. Differential expression design scope: simple two-group contrasts in v1.
7. Multi-factor, paired, and batch-aware designs are deferred.

## 3.1 Remaining Questions For Review

Please review these before coding starts:

1. Should plots be static PNG/PDF only, or also interactive HTML?
2. Should the report output be Markdown only, HTML only, or both?
3. Should we add a workflow runner such as Snakemake later, or keep the first version as standalone scripts plus a Python orchestrator?
4. For mouse and human annotation, should we vendor pinned mapping files in `data/external/`, or provide a script that downloads pinned releases by checksum?

## 4. Proposed Repository Layout

```text
Try_DEG_Pipeline/
├── README.md
├── pyproject.toml
├── environment.yml
├── configs/
│   ├── airway.yaml
│   └── example_project.yaml
├── data/
│   ├── raw/
│   │   └── .gitkeep
│   ├── external/
│   │   └── .gitkeep
│   └── fixtures/
│       ├── toy_counts.tsv
│       ├── toy_metadata.tsv
│       ├── toy_annotation.tsv
│       └── toy_go_map.tsv
├── results/
│   └── .gitkeep
├── scripts/
│   ├── 00_fetch_airway.py
│   ├── 01_validate_inputs.py
│   ├── 02_qc_counts.py
│   ├── 03_filter_normalize.py
│   ├── 04_run_pydeseq2.py
│   ├── 05_validate_airway_pydeseq2.py
│   ├── 06_make_deg_outputs.py
│   ├── 07_run_go_enrichment.py
│   ├── 08_make_report.py
│   └── run_pipeline.py
├── src/
│   └── deg_pipeline/
│       ├── __init__.py
│       ├── io.py
│       ├── validation.py
│       ├── qc.py
│       ├── pydeseq2_adapter.py
│       ├── deg.py
│       ├── enrichment.py
│       ├── plotting.py
│       └── logging_utils.py
└── tests/
    ├── conftest.py
    ├── unit/
    │   ├── test_validate_inputs.py
    │   ├── test_qc_counts.py
    │   ├── test_filter_normalize.py
    │   ├── test_pydeseq2_adapter.py
    │   ├── test_deg_outputs.py
    │   ├── test_go_enrichment.py
    │   └── test_report.py
    └── smoke/
        ├── test_00_fetch_airway_smoke.py
        ├── test_01_validate_inputs_smoke.py
        ├── test_02_qc_counts_smoke.py
        ├── test_03_filter_normalize_smoke.py
        ├── test_04_run_pydeseq2_smoke.py
        ├── test_05_validate_airway_pydeseq2_smoke.py
        ├── test_06_make_deg_outputs_smoke.py
        ├── test_07_run_go_enrichment_smoke.py
        └── test_08_make_report_smoke.py
```

The standalone scripts should import shared utilities from `src/deg_pipeline/`. This avoids copy-paste while preserving one executable script per pipeline step.

## 5. Data Contracts

### 5.1 Count Matrix

Input format: tab-separated text file.

Required shape:

- Rows are genes/features.
- Columns are samples.
- First column is a stable gene identifier.
- Count values are non-negative integers.

Example:

```text
gene_id	sample_1	sample_2	sample_3	sample_4
ENSG000001	100	120	12	10
ENSG000002	0	1	0	0
ENSG000003	500	480	950	1000
```

Validation rules:

- Gene IDs must be unique and non-empty.
- Sample columns must be unique and non-empty.
- Values must be numeric, finite, integer-like, and non-negative.
- No duplicated sample names after trimming whitespace.
- No accidental normalized values. The script should flag suspicious non-integer values.

### 5.2 Metadata Table

Input format: tab-separated text file.

Required columns:

- `sample_id`
- Condition column, default `condition`

Optional columns:

- Descriptive sample covariates for plotting or future versions.
- Batch and paired-subject covariates may be present but are not modeled in v1.

Example:

```text
sample_id	condition	batch
sample_1	control	A
sample_2	control	B
sample_3	treated	A
sample_4	treated	B
```

Validation rules:

- Every count column must have exactly one metadata row.
- Every metadata sample must exist in counts unless `--allow-extra-metadata` is passed.
- The configured condition variable must exist and must contain exactly two levels for v1.
- Contrast levels must exist.
- Sample order must be made explicit in generated output.

### 5.3 Annotation Table

Input format: tab-separated text file.

Required columns for DEG output:

- `gene_id`

Optional but recommended:

- `gene_symbol`
- `entrez_id`
- `gene_name`
- `biotype`

Required for offline GO enrichment:

- A mappable identifier from Ensembl gene ID to GO-compatible IDs.
- Species must be either `human` or `mouse` in v1.

### 5.4 Analysis Config

Config format: YAML.

Example:

```yaml
project_name: airway_demo
species: human
gene_id_type: ensembl
condition_column: dex
reference_level: untrt
results_dir: results/airway_demo

design:
  formula: "~ dex"
  variables:
    - dex

contrasts:
  - name: dex_trt_vs_untrt
    variable: dex
    numerator: trt
    denominator: untrt

filtering:
  min_count: 10
  min_samples: 3

pydeseq2:
  alpha: 0.05
  cooks_filter: true
  independent_filter: true

deg:
  padj_threshold: 0.05
  log2fc_threshold: 1.0

go:
  ontology:
    - BP
    - MF
    - CC
  method: offline_ora
  padj_method: fdr_bh
  min_genes_per_term: 5
  max_genes_per_term: 500
```

## 6. Pipeline Stages

### Stage 0: Fetch And Prepare Airway Dataset

Script: `scripts/00_fetch_airway.py`

Purpose:

- Download or extract the public Bioconductor `airway` dataset.
- Export raw gene-level counts to TSV.
- Export sample metadata to TSV.
- Export a small reduced fixture for smoke tests.
- Generate or download a pinned R/Bioconductor DESeq2 reference result file for PyDESeq2 validation.

Expected outputs:

```text
data/external/airway/counts.tsv
data/external/airway/metadata.tsv
data/external/airway/annotation.tsv
data/external/airway/reference_deseq2_dex_trt_vs_untrt.tsv
data/fixtures/airway_small_counts.tsv
data/fixtures/airway_small_metadata.tsv
results/airway_fetch/session_info.txt
results/airway_fetch/fetch_manifest.json
```

Implementation note:

- Because `airway` is a Bioconductor dataset, the most practical route is to let this setup script call `Rscript` only for dataset export and reference generation.
- The production pipeline remains Python-native. R is used for validation data preparation, not for routine differential expression.

### Stage 1: Validate Inputs

Script: `scripts/01_validate_inputs.py`

Purpose:

- Validate counts, metadata, annotation, and config before expensive analysis.
- Produce a JSON validation report and a human-readable text report.

Expected outputs:

```text
results/<project>/validation/input_validation.json
results/<project>/validation/input_validation.txt
results/<project>/intermediate/counts_ordered.tsv
results/<project>/intermediate/metadata_ordered.tsv
```

Key checks:

- File readability.
- Required columns.
- Count matrix integrity.
- Count/metadata sample matching.
- Condition variable.
- Contrast levels.
- Annotation gene ID overlap.
- Basic library size sanity.

### Stage 2: Count QC

Script: `scripts/02_qc_counts.py`

Purpose:

- Generate sample-level and gene-level QC from raw counts.
- Identify outliers before filtering and PyDESeq2.

Expected outputs:

```text
results/<project>/qc/library_sizes.tsv
results/<project>/qc/gene_detection.tsv
results/<project>/qc/sample_correlation.tsv
results/<project>/qc/pca_coordinates.tsv
results/<project>/qc/plots/library_sizes.png
results/<project>/qc/plots/detected_genes.png
results/<project>/qc/plots/sample_correlation_heatmap.png
results/<project>/qc/plots/pca_raw_or_vst.png
results/<project>/qc/qc_summary.json
```

QC metrics:

- Total reads per sample.
- Number of detected genes per sample.
- Fraction of zero counts per sample.
- Count distribution summaries.
- Sample correlation.
- PCA coordinates, preferably after variance-stabilizing transform if PyDESeq2 transform output is available; otherwise use log2(count + 1) for pre-PyDESeq2 QC.

### Stage 3: Filter And Normalize Prep

Script: `scripts/03_filter_normalize.py`

Purpose:

- Apply reproducible low-count gene filtering before PyDESeq2.
- Preserve raw counts for PyDESeq2; do not replace PyDESeq2 size factor estimation with custom normalized counts.
- Write a filtering manifest.

Expected outputs:

```text
results/<project>/intermediate/counts_filtered.tsv
results/<project>/filtering/gene_filtering_summary.tsv
results/<project>/filtering/gene_filtering_manifest.json
```

Default filter:

- Keep genes with at least `min_count` reads in at least `min_samples` samples.
- Defaults: `min_count = 10`, `min_samples = smallest condition group size`.

### Stage 4: Run PyDESeq2

Script: `scripts/04_run_pydeseq2.py`

Purpose:

- Execute Python-native DESeq2-like analysis with `PyDESeq2`.
- Support a simple two-group condition contrast in v1.
- Export normalized counts, size factors, dispersion estimates, transform matrix when available, and result tables.
- Save Python package versions and PyDESeq2 settings for reproducibility.

Expected outputs:

```text
results/<project>/pydeseq2/size_factors.tsv
results/<project>/pydeseq2/normalized_counts.tsv
results/<project>/pydeseq2/vst_or_log_counts.tsv
results/<project>/pydeseq2/dispersion_estimates.tsv
results/<project>/pydeseq2/results_<contrast>.tsv
results/<project>/pydeseq2/pydeseq2_run_manifest.json
```

PyDESeq2 defaults:

- Use raw integer counts.
- Model exactly one condition variable in v1.
- Require exactly two condition levels in v1.
- Use the configured denominator as the reference level.
- Run `DeseqDataSet.deseq2()`.
- Run `DeseqStats.summary()` for the configured contrast.
- Export result schema compatible with downstream DEG code: `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`.
- Do not claim R `lfcShrink` behavior in v1. Add shrinkage only if PyDESeq2 support is explicit and tested.

### Stage 5: Validate Airway PyDESeq2 Behavior

Script: `scripts/05_validate_airway_pydeseq2.py`

Purpose:

- Confirm that PyDESeq2 pipeline behavior on `airway` is scientifically consistent with a pinned R/Bioconductor DESeq2 reference.
- Compare generated outputs to stable reference invariants instead of requiring exact bitwise equality, because PyDESeq2 and R DESeq2 may differ slightly.

Expected outputs:

```text
results/airway_demo/validation/airway_pydeseq2_validation.json
results/airway_demo/validation/airway_pydeseq2_validation.txt
```

Validation checks:

- Expected sample count equals 8 for full airway.
- Expected primary condition is dexamethasone treatment status.
- PyDESeq2 result table contains required columns: `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`.
- Known top genes appear in the top-ranked result set within tolerance.
- Direction of selected known treatment-responsive genes is stable.
- Correlation between current and reference log2 fold changes exceeds a threshold.
- Correlation between current and reference test statistics exceeds a threshold.
- Number of significant genes under `padj < 0.05` is within a configured tolerance.
- No unexpected NA explosion in `padj`.

### Stage 6: Make DEG Outputs

Script: `scripts/06_make_deg_outputs.py`

Purpose:

- Convert raw PyDESeq2 results into clear analysis deliverables.
- Join annotation.
- Classify genes as up/down/not significant.
- Export ranked lists for downstream enrichment.

Expected outputs:

```text
results/<project>/deg/<contrast>_all_genes.tsv
results/<project>/deg/<contrast>_significant.tsv
results/<project>/deg/<contrast>_up.tsv
results/<project>/deg/<contrast>_down.tsv
results/<project>/deg/<contrast>_ranked_for_gsea.rnk
results/<project>/deg/<contrast>_summary.json
results/<project>/deg/plots/<contrast>_volcano.png
results/<project>/deg/plots/<contrast>_ma_plot.png
results/<project>/deg/plots/<contrast>_top_genes_heatmap.png
```

Default DEG classification:

- Significant if `padj <= 0.05` and `abs(log2FoldChange) >= 1.0`.
- Up if significant and `log2FoldChange > 0`.
- Down if significant and `log2FoldChange < 0`.

### Stage 7: GO Enrichment

Script: `scripts/07_run_go_enrichment.py`

Purpose:

- Run Gene Ontology enrichment for upregulated and downregulated DEG sets.
- Use an explicit universe/background, preferably all genes that passed filtering and could be mapped to GO.

Expected outputs:

```text
results/<project>/go/<contrast>_up_go_bp.tsv
results/<project>/go/<contrast>_down_go_bp.tsv
results/<project>/go/<contrast>_up_go_mf.tsv
results/<project>/go/<contrast>_down_go_mf.tsv
results/<project>/go/<contrast>_up_go_cc.tsv
results/<project>/go/<contrast>_down_go_cc.tsv
results/<project>/go/<contrast>_go_summary.json
results/<project>/go/plots/<contrast>_up_go_bp_dotplot.png
results/<project>/go/plots/<contrast>_down_go_bp_dotplot.png
```

Recommended first implementation:

- Offline over-representation analysis using `goatools`.
- Pin/download `go-basic.obo`.
- Pin/download species gene-to-GO mapping.
- Use Benjamini-Hochberg FDR correction.

Non-goal for v1:

- Online enrichment services such as Enrichr or g:Profiler. They are useful for exploration, but not part of the reproducible default.

### Stage 8: Report Bundle

Script: `scripts/08_make_report.py`

Purpose:

- Generate one compact HTML or Markdown report that links all outputs.
- Summarize input data, QC, PyDESeq2 model, contrasts, DEG counts, GO terms, warnings, and package/session info.

Expected outputs:

```text
results/<project>/report/report.md
results/<project>/report/report.html
results/<project>/report/assets/
results/<project>/run_manifest.json
```

Report sections:

- Project and run metadata.
- Input files and checksums.
- Analysis config.
- Sample metadata overview.
- QC summary.
- Filtering summary.
- PyDESeq2 condition model and contrast.
- DEG summary per contrast.
- Top genes per contrast.
- GO enrichment summary.
- Reproducibility appendix.

### Optional Stage 9: Orchestrator

Script: `scripts/run_pipeline.py`

Purpose:

- Run all stages in order from a config file.
- Stop on failure.
- Write a master manifest.

This should not hide the standalone scripts. It should simply call them with explicit arguments.

## 7. Dependency Plan

### Python

Recommended Python version:

- Python 3.11 or 3.12.

Core packages:

- `pandas`
- `numpy`
- `scipy`
- `statsmodels`
- `pydeseq2`
- `pyyaml`
- `pydantic` or `jsonschema`
- `matplotlib`
- `seaborn`
- `scikit-learn`
- `jinja2`
- `pytest`
- `goatools`

Optional packages:

- `plotly` for interactive plots.

### R/Bioconductor

R is not required for routine production analysis in v1. It is required only for the `airway` dataset export and R DESeq2 reference-generation validation path.

Validation-only R packages:

- `DESeq2`
- `airway`
- `SummarizedExperiment`
- `BiocManager`

The pipeline should capture:

- R version.
- Bioconductor version.
- R DESeq2 version.
- Full `sessionInfo()`.

Implementation note:

- Pin `pydeseq2` in `environment.yml` and write the exact version to every PyDESeq2 manifest.
- The current PyDESeq2 API supports `DeseqDataSet(..., design="~condition", ref_level=[condition, reference])`; implementation should confirm this in tests and avoid relying on deprecated parameters.

## 8. Reproducibility Requirements

Every run should write:

- Config snapshot.
- Input file checksums.
- Script name and version.
- Git commit hash when available.
- Python package versions.
- PyDESeq2 version and settings.
- R package versions only for validation runs that generate R DESeq2 reference output.
- CLI arguments.
- Start/end timestamps.
- Output file manifest.

No script should silently overwrite results unless `--force` is passed.

## 9. Testing Strategy

### Unit Tests

Each script should expose testable helper functions through `src/deg_pipeline/`.

Unit tests should avoid requiring R. R-dependent tests are limited to airway reference export and validation fixtures. Examples:

- Count matrix validation accepts valid counts.
- Count matrix validation rejects negative values.
- Metadata validation rejects missing contrast levels.
- Filtering keeps expected toy genes.
- DEG classification handles `padj = NA`.
- GO enrichment uses the configured universe.
- Report renderer includes expected sections.

### Smoke Tests

Smoke tests should execute each script as a CLI on tiny inputs.

Smoke test rules:

- Use `tmp_path` outputs.
- Keep runtime short.
- Assert output files exist.
- Assert key columns and row counts.
- Mark tests requiring R/reference DESeq2 as `pytest.mark.requires_r`.
- Full `airway` smoke validation can be slower and optional in default local testing.

Suggested commands:

```bash
pytest tests/unit
pytest tests/smoke
pytest -m requires_r
pytest -m airway
```

### Validation Tests

Validation tests differ from unit/smoke tests. They prove scientific behavior against `airway`.

Recommended validation checks:

- Pipeline completes on full `airway`.
- PyDESeq2 result schema is correct.
- Log2 fold changes correlate strongly with reference.
- Significant gene count is within tolerance.
- Selected known genes have expected direction.
- GO enrichment output is non-empty for at least one DEG set when thresholds produce enough genes.

## 10. Pseudo-code Skeletons

The following skeletons are intentionally precise enough for coding agents to implement independently.

### 10.1 Shared Script Pattern

```python
#!/usr/bin/env python3
"""One-sentence description of this pipeline stage."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from deg_pipeline.logging_utils import configure_logging, write_manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--force", action="store_true")
    return parser


def run(args: argparse.Namespace) -> int:
    configure_logging(args.outdir)
    # 1. Load config.
    # 2. Validate required inputs for this stage.
    # 3. Execute pure helper functions from src/deg_pipeline.
    # 4. Write deterministic outputs.
    # 5. Write manifest.
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
```

### 10.2 `00_fetch_airway.py`

```python
def build_parser():
    add --outdir
    add --reduced-fixture-outdir
    add --n-top-variable-genes default 200
    add --force


def render_airway_export_r_script(outdir, make_reference):
    return string containing R code that:
        library(airway)
        library(SummarizedExperiment)
        load airway SummarizedExperiment
        extract assay counts
        extract colData metadata
        write counts.tsv
        write metadata.tsv
        if make_reference:
            library(DESeq2)
            build DESeqDataSet
            run DESeq
            write reference DESeq2 results
        write sessionInfo


def run(args):
    ensure output directories exist
    check Rscript is available
    write temporary R script
    execute Rscript with subprocess.run(check=True)
    read exported counts and metadata
    create reduced fixture:
        choose all samples
        choose top variable genes plus a few low-count genes
    write fixture counts and metadata
    write fetch_manifest.json with checksums
```

Unit test:

```python
def test_render_airway_export_contains_required_packages():
    script = render_airway_export_r_script(...)
    assert "library(airway)" in script
    assert "write.table" in script
```

Smoke test:

```python
@pytest.mark.requires_r
def test_fetch_airway_cli_creates_counts_and_metadata(tmp_path):
    run script with tmp_path
    assert counts.tsv exists
    assert metadata.tsv exists
```

### 10.3 `01_validate_inputs.py`

```python
def build_parser():
    add --counts
    add --metadata
    add --config
    add --annotation optional
    add --outdir
    add --allow-extra-metadata flag


def validate_count_matrix(counts_df):
    errors = []
    warnings = []
    check gene_id column
    check duplicate gene IDs
    check sample columns
    check numeric finite values
    check integer-like values
    check non-negative values
    return ValidationResult(errors, warnings, metrics)


def validate_metadata(metadata_df, sample_names, config):
    check sample_id column
    check unique sample IDs
    check every count sample has metadata
    check condition variable exists
    check condition variable has exactly two levels
    check contrast levels exist
    return ValidationResult(...)


def validate_annotation(annotation_df, gene_ids):
    check gene_id column
    compute overlap rate
    warn if low overlap
    return ValidationResult(...)


def run(args):
    load counts, metadata, config, optional annotation
    run validators
    if errors:
        write reports
        return nonzero
    reorder metadata to match count columns
    write ordered counts and metadata
    write JSON and text report
```

Unit test:

```python
def test_validate_count_matrix_rejects_negative_counts():
    df = make_counts_with_negative_value()
    result = validate_count_matrix(df)
    assert result.has_errors
```

Smoke test:

```python
def test_validate_inputs_cli_writes_ordered_outputs(tmp_path):
    run script on toy_counts.tsv and toy_metadata.tsv
    assert input_validation.json exists
    assert metadata_ordered sample order equals counts columns
```

### 10.4 `02_qc_counts.py`

```python
def build_parser():
    add --counts
    add --metadata
    add --config
    add --outdir


def compute_library_metrics(counts_df):
    for each sample:
        total_counts
        detected_genes
        zero_fraction
        upper_quartile_count
    return sample_metrics_df


def compute_gene_metrics(counts_df):
    for each gene:
        total_counts
        mean_counts
        detected_samples
    return gene_metrics_df


def compute_sample_correlation(counts_df):
    transformed = log2(counts + 1)
    return transformed.corr(method="spearman")


def compute_pca(counts_df, metadata_df):
    transformed = log2(counts + 1)
    center and scale genes
    run sklearn PCA
    return coordinates and variance explained


def run(args):
    load ordered counts and metadata
    compute metrics
    write tables
    create library size plot
    create detected genes plot
    create sample correlation heatmap
    create PCA plot colored by condition
    write qc_summary.json
```

Unit test:

```python
def test_compute_library_metrics_counts_detected_genes():
    metrics = compute_library_metrics(toy_counts)
    assert metrics.loc["sample_1", "detected_genes"] == expected
```

Smoke test:

```python
def test_qc_counts_cli_writes_metric_tables_and_plots(tmp_path):
    run script on toy ordered counts
    assert library_sizes.tsv exists
    assert pca_coordinates.tsv exists
```

### 10.5 `03_filter_normalize.py`

```python
def build_parser():
    add --counts
    add --metadata
    add --config
    add --outdir
    add --min-count optional
    add --min-samples optional


def determine_min_samples(metadata_df, condition_column, configured_min_samples):
    if configured_min_samples is set:
        return configured_min_samples
    return size of smallest condition group


def filter_low_count_genes(counts_df, min_count, min_samples):
    keep_mask = number of samples with count >= min_count >= min_samples
    return filtered_counts_df, filtering_summary_df


def run(args):
    load counts, metadata, config
    determine thresholds
    filter genes
    fail if too few genes remain
    write counts_filtered.tsv
    write gene_filtering_summary.tsv
    write manifest
```

Unit test:

```python
def test_filter_low_count_genes_keeps_expected_gene():
    filtered, summary = filter_low_count_genes(toy_counts, min_count=10, min_samples=2)
    assert "gene_high" in filtered.index
    assert "gene_low" not in filtered.index
```

Smoke test:

```python
def test_filter_normalize_cli_writes_filtered_counts(tmp_path):
    run script on toy counts
    assert counts_filtered.tsv exists
```

### 10.6 `04_run_pydeseq2.py`

```python
def build_parser():
    add --counts-filtered
    add --metadata
    add --config
    add --outdir
    add --contrast-name optional


def prepare_pydeseq2_inputs(counts_df, metadata_df, config):
    ensure counts matrix is samples x genes for PyDESeq2
    ensure metadata index is sample_id
    ensure exactly two condition levels
    relevel condition so denominator/reference is first where PyDESeq2 requires it
    return counts_for_pydeseq2, metadata_for_pydeseq2


def run_pydeseq2_model(counts_df, metadata_df, config):
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata_df,
        design=f"~{config.condition_column}",
        ref_level=[config.condition_column, config.reference_level],
        refit_cooks=config.pydeseq2.cooks_filter,
    )
    dds.deseq2()
    return dds


def run_pydeseq2_stats(dds, config, contrast):
    stat_res = DeseqStats(
        dds,
        contrast=[contrast.variable, contrast.numerator, contrast.denominator],
        alpha=config.pydeseq2.alpha,
        cooks_filter=config.pydeseq2.cooks_filter,
        independent_filter=config.pydeseq2.independent_filter,
    )
    stat_res.summary()
    return stat_res.results_df


def export_pydeseq2_outputs(dds, results_df, outdir, contrast_name):
    write size factors if exposed by PyDESeq2 object
    write normalized counts if exposed by PyDESeq2 object
    write dispersion estimates if exposed by PyDESeq2 object
    write result table with gene_id as first column
    write package version and model settings manifest


def validate_pydeseq2_outputs(outdir, contrast_name):
    assert required files exist
    assert result files have required columns
    return output_manifest


def run(args):
    load config
    load filtered counts and ordered metadata
    prepare PyDESeq2 inputs
    run PyDESeq2 model
    run PyDESeq2 stats for configured simple contrast
    export outputs
    validate outputs
    write manifest
```

Unit test:

```python
def test_prepare_pydeseq2_inputs_rejects_three_condition_levels():
    with pytest.raises(ValueError, match="exactly two"):
        prepare_pydeseq2_inputs(counts, metadata_with_three_levels, config)
```

Smoke test:

```python
def test_run_pydeseq2_cli_writes_result_table(tmp_path):
    run script on small fixture with two groups
    assert results_test_contrast.tsv exists
    assert required PyDESeq2-compatible columns are present
```

### 10.7 `05_validate_airway_pydeseq2.py`

```python
def build_parser():
    add --pipeline-results
    add --reference-results
    add --metadata
    add --outdir
    add --contrast-name
    add --lfc-correlation-threshold default 0.95
    add --stat-correlation-threshold default 0.95
    add --significant-count-tolerance default 0.15


def load_and_align_results(current_path, reference_path):
    read both tables
    align by gene_id
    keep genes with finite compared values
    return aligned_df


def compute_airway_validation_metrics(aligned_df):
    lfc_correlation = corr(current_log2fc, reference_log2fc)
    stat_correlation = corr(current_stat, reference_stat)
    significant_counts = compare padj < 0.05 counts
    top_gene_overlap = overlap of top N by padj
    return metrics


def evaluate_validation(metrics, thresholds):
    create pass/fail checks
    return validation_report


def run(args):
    aligned = load_and_align_results(...)
    metrics = compute_airway_validation_metrics(aligned)
    report = evaluate_validation(metrics, thresholds)
    write json and text report
    return nonzero if critical checks fail
```

Unit test:

```python
def test_compute_airway_validation_metrics_detects_high_correlation():
    metrics = compute_airway_validation_metrics(fake_aligned_results)
    assert metrics["lfc_correlation"] > 0.99
```

Smoke test:

```python
def test_validate_airway_pydeseq2_cli_passes_on_tiny_reference(tmp_path):
    create tiny current and reference result files
    run script
    assert validation json says pass
```

### 10.8 `06_make_deg_outputs.py`

```python
def build_parser():
    add --pydeseq2-results-dir
    add --annotation optional
    add --config
    add --outdir


def classify_degs(results_df, padj_threshold, log2fc_threshold):
    for each gene:
        if padj is NA:
            class = "not_tested"
        elif padj <= threshold and log2fc >= threshold:
            class = "up"
        elif padj <= threshold and log2fc <= -threshold:
            class = "down"
        elif padj <= threshold:
            class = "significant_low_lfc"
        else:
            class = "not_significant"
    return annotated_results_df


def join_annotation(results_df, annotation_df):
    left join by gene_id
    preserve all PyDESeq2 rows
    return merged_df


def make_ranked_file(results_df):
    ranking_metric = sign(log2FoldChange) * -log10(pvalue)
    handle pvalue zero by clipping to smallest positive float
    write two-column RNK: gene_id_or_symbol, ranking_metric


def make_deg_plots(results_df, normalized_counts, metadata):
    volcano
    MA plot
    top genes heatmap


def run(args):
    for each contrast result file:
        load result
        join annotation
        classify DEG
        write all/significant/up/down tables
        write ranked file
        write summary json
        write plots
```

Unit test:

```python
def test_classify_degs_handles_up_down_and_na():
    classified = classify_degs(fake_results, 0.05, 1.0)
    assert classified.loc["gene_up", "deg_class"] == "up"
    assert classified.loc["gene_na", "deg_class"] == "not_tested"
```

Smoke test:

```python
def test_make_deg_outputs_cli_writes_deg_tables(tmp_path):
    run script on fake PyDESeq2 result table
    assert all_genes.tsv exists
    assert up.tsv exists
    assert volcano.png exists
```

### 10.9 `07_run_go_enrichment.py`

```python
def build_parser():
    add --deg-dir
    add --annotation
    add --go-obo
    add --gene2go
    add --config
    add --outdir


def build_gene_universe(filtered_gene_ids, annotation_df):
    map pipeline gene IDs to GO-compatible IDs
    keep genes with valid GO mapping
    return universe_ids


def load_deg_gene_sets(deg_dir, contrast):
    read up and down DEG tables
    map to GO-compatible IDs
    return {"up": set(...), "down": set(...)}


def run_go_ora(gene_set, universe, ontology, go_obo, gene2go):
    initialize goatools study object
    run enrichment
    adjust p-values
    return results_df


def make_go_dotplot(results_df, outpath):
    choose top terms by adjusted p-value
    plot term name, gene ratio, adjusted p-value


def run(args):
    load config and annotations
    universe = build_gene_universe(...)
    for each contrast and direction:
        for ontology in BP/MF/CC:
            results = run_go_ora(...)
            write table
            write plot if non-empty
    write summary json
```

Unit test:

```python
def test_build_gene_universe_excludes_unmapped_genes():
    universe = build_gene_universe(["g1", "g2"], annotation_with_one_missing)
    assert "mapped_entrez" in universe
    assert "missing" not in universe
```

Smoke test:

```python
def test_run_go_enrichment_cli_writes_results_on_toy_go(tmp_path):
    run script with tiny ontology and gene2go fixture
    assert go result table exists
```

### 10.10 `08_make_report.py`

```python
def build_parser():
    add --project-dir
    add --config
    add --outdir


def collect_run_artifacts(project_dir):
    find validation reports
    find qc summaries
    find filtering manifest
    find pydeseq2 manifest
    find deg summaries
    find go summaries
    return artifact_index


def render_markdown_report(artifact_index, config):
    use jinja2 template
    include summary tables
    link plots and output files
    include warnings
    include reproducibility appendix
    return markdown


def convert_markdown_to_html(markdown_path):
    use optional markdown package or pandoc if configured


def run(args):
    collect artifacts
    render markdown
    render html
    copy or link plot assets
    write report manifest
```

Unit test:

```python
def test_render_markdown_report_includes_deg_summary():
    markdown = render_markdown_report(fake_artifacts, fake_config)
    assert "DEG Summary" in markdown
```

Smoke test:

```python
def test_make_report_cli_writes_report(tmp_path):
    create minimal fake project artifacts
    run report script
    assert report.md exists
```

### 10.11 `run_pipeline.py`

```python
def build_parser():
    add --config
    add --counts
    add --metadata
    add --annotation optional
    add --outdir
    add --start-at optional
    add --stop-after optional
    add --force


def build_stage_commands(args, config):
    return ordered list of command arrays:
        python scripts/01_validate_inputs.py ...
        python scripts/02_qc_counts.py ...
        python scripts/03_filter_normalize.py ...
        python scripts/04_run_pydeseq2.py ...
        python scripts/06_make_deg_outputs.py ...
        python scripts/07_run_go_enrichment.py ...
        python scripts/08_make_report.py ...


def run_command(command, log_path):
    subprocess.run(command, check=True)
    capture stdout/stderr


def run(args):
    load config
    commands = build_stage_commands(...)
    apply start/stop slicing
    for command in commands:
        run command
        stop immediately on failure
    write master run manifest
```

Unit test:

```python
def test_build_stage_commands_preserves_stage_order():
    commands = build_stage_commands(args, config)
    assert "01_validate_inputs.py" in commands[0]
    assert "08_make_report.py" in commands[-1]
```

Smoke test:

```python
def test_run_pipeline_dry_run_lists_commands(tmp_path):
    run orchestrator with --dry-run
    assert command list contains all expected stages
```

## 11. Suggested Implementation Milestones

### Milestone 1: Project Scaffolding

Deliverables:

- `pyproject.toml`
- `environment.yml`
- `src/deg_pipeline/` package shell
- `tests/` layout
- Toy fixtures
- Basic logging and manifest utilities

Acceptance:

- `pytest tests/unit` runs.
- CLI scripts respond to `--help`.

### Milestone 2: Input Validation And QC

Deliverables:

- `01_validate_inputs.py`
- `02_qc_counts.py`
- Unit and smoke tests.

Acceptance:

- Toy input validation passes.
- Known malformed fixtures fail with clear errors.
- QC tables and plots are generated.

### Milestone 3: Filtering And PyDESeq2 Adapter

Deliverables:

- `03_filter_normalize.py`
- `04_run_pydeseq2.py`
- PyDESeq2 version and run settings capture.
- Unit and smoke tests.

Acceptance:

- PyDESeq2 runs on a tiny fixture.
- Output schemas are stable.

### Milestone 4: Airway Validation

Deliverables:

- `00_fetch_airway.py`
- `configs/airway.yaml`
- `05_validate_airway_pydeseq2.py`

Acceptance:

- Full `airway` run completes.
- Validation report passes configured thresholds.

### Milestone 5: DEG Deliverables

Deliverables:

- `06_make_deg_outputs.py`
- DEG summary tables.
- Volcano, MA, and heatmap plots.

Acceptance:

- DEG outputs are produced for airway contrast.
- DEG classifications are tested.

### Milestone 6: GO Enrichment

Deliverables:

- `07_run_go_enrichment.py`
- Offline GO fixture for tests.
- Human and mouse GO configs for real runs.

Acceptance:

- Toy GO smoke test passes.
- Airway DEG set produces enrichment tables or clear "not enough genes" warnings.

### Milestone 7: Report And Orchestration

Deliverables:

- `08_make_report.py`
- `run_pipeline.py`
- End-to-end smoke command.

Acceptance:

- One command can run the full pipeline.
- Report summarizes all stages.

## 12. Acceptance Criteria For Version 1

Version 1 is complete when:

- A user can provide counts, metadata, annotation, and YAML config.
- Inputs are validated before PyDESeq2 runs.
- QC outputs are generated.
- Low-count filtering is deterministic.
- PyDESeq2 results are generated through a standalone Python script.
- Airway validation against pinned R/Bioconductor DESeq2 reference behavior passes.
- DEG tables and plots are generated.
- Offline GO enrichment is available.
- A report bundle is generated.
- Every script has a unit test and smoke test.
- The pipeline records enough environment information to reproduce the run.

## 13. Remaining User Decisions Needed

Before implementation, please decide or comment on:

1. Preferred final report format: Markdown, HTML, or both.
2. Static plots only, or static plus interactive HTML plots.
3. Whether to add Snakemake/Nextflow after the standalone scripts are working.
4. Whether human/mouse annotation resources should be committed as pinned files or downloaded by a pinned-resource script.
