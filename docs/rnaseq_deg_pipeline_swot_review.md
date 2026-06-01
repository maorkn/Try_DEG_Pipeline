# SWOT Review: Reproducible RNA-seq DEG Pipeline Plan

Review target: `docs/rnaseq_deg_pipeline_implementation_plan.md`  
Review date: 2026-06-01  
Reviewer stance: independent technical review focused on reproducibility, statistical validity, implementation risk, testability, and handoff clarity.

## Scope Assumptions

This review evaluates the plan against the current v1 decisions:

- Primary DEG engine: `PyDESeq2`.
- Species: human and mouse.
- Gene IDs: Ensembl.
- GO enrichment: fully offline.
- Statistical design: simple two-group contrasts only.

## Strengths

- The plan has a strong modular architecture. The file-in/file-out stage design, standalone CLI scripts, and shared `src/deg_pipeline/` helpers make the work easy to split across coding agents while keeping testable business logic out of the scripts.
- The v1 scope is appropriately constrained. Simple two-group contrasts, Ensembl IDs, human/mouse support, and offline over-representation analysis are realistic for a first reproducible implementation.
- The locked decisions are now explicit near the top of the plan, including `PyDESeq2` as the production engine and R/Bioconductor `DESeq2` as validation-only context.
- The data contracts for counts, metadata, annotation, and config cover the main failure modes that commonly invalidate RNA-seq DEG runs: mismatched sample IDs, non-integer counts, duplicate gene IDs, missing contrast levels, and weak annotation overlap.
- The plan correctly preserves raw integer counts for `PyDESeq2` and treats filtering as deterministic preprocessing rather than custom normalization.
- The test strategy is better than typical pipeline plans. It separates unit tests, smoke tests, and scientific validation tests, and it avoids making R a default dependency for normal unit tests.
- The pseudocode skeletons give coding agents a practical handoff point. They identify CLI arguments, core helper functions, outputs, and expected tests for each stage.
- Reproducibility requirements are broad and relevant: config snapshots, checksums, Git hash, package versions, CLI args, manifests, and overwrite protection.
- The GO plan correctly names an explicit universe/background as a requirement, which is essential for statistically meaningful enrichment.
- The plan avoids claiming R `lfcShrink` equivalence in v1, which is important because shrinkage behavior and availability differ across engines.

## Weaknesses

- The offline GO strategy is not fully specified. The plan says to "pin/download" `go-basic.obo` and species gene-to-GO mapping, but fully offline enrichment requires resources to already exist locally at runtime, with exact versions, checksums, schemas, and provenance recorded. Runtime download should be outside the analysis path.
- Human and mouse Ensembl-to-GO mapping is underdefined. `goatools` commonly works with Entrez-style gene2go mappings, while the pipeline standard is Ensembl. The plan needs a concrete mapping contract for Ensembl IDs, including version suffix handling such as `ENSG... .15`, one-to-many mappings, unmapped genes, and whether the GO universe is defined before or after mapping.
- The validation approach may overstate comparability between `PyDESeq2` and R `DESeq2`. Correlation thresholds such as 0.95 for log2 fold change and statistics may be brittle across package versions, default changes, size factor behavior, independent filtering, Cook's filtering, and factor reference handling.
- R is still required for the full `airway` export/reference path. That is reasonable for validation, but the plan should distinguish "production pipeline works without R" from "full validation bundle requires R/Bioconductor" in install, CI, and smoke-test expectations.
- There is no mouse validation dataset equivalent to `airway`. Human support is validated with `airway`, but mouse support currently appears to rely on annotation/config behavior and toy GO fixtures. That is acceptable for v1 only if stated explicitly.
- Statistical validity guardrails need more detail. For simple two-group contrasts, the plan should reject or loudly warn on too few replicates per group, perfect confounding with batch-like metadata, all-zero groups for many genes, and designs where the requested contrast direction is ambiguous.
- The config schema is sketched but not yet strict enough for independent implementation. Agents need exact required fields, defaults, accepted species names, accepted ID formats, path resolution rules, and validation behavior for multiple contrasts.
- Output schemas are described by examples, but not formalized. For handoff clarity, each stage should define required columns, column types, index conventions, and whether gene IDs are stored as first column or row index.
- `03_filter_normalize.py` is slightly misleadingly named because it should filter only and preserve raw counts. The current name may encourage agents to add normalization where it does not belong.
- The QC plan mentions VST "if PyDESeq2 transform output is available" before the model stage. Pre-model QC and post-model transformed QC should be separated to avoid hidden dependencies between stages.
- The report and plot deliverables may be too broad for the same v1 milestone as statistical correctness. They are useful, but they increase implementation surface area and testing burden before the core model/enrichment path is proven.
- Smoke tests for `PyDESeq2` can be slow or fragile on very tiny fixtures. The plan should define the minimum fixture size and expected behavior when `PyDESeq2` cannot estimate dispersions reliably on toy data.
- The plan does not yet define a reproducible dependency lock strategy. `environment.yml` is listed, but there is no explicit policy for conda lockfiles, pip hashes, exact PyDESeq2 version pins, or Bioconductor snapshot/version pinning.
- The GO enrichment test plan says output should be non-empty for at least one DEG set when thresholds produce enough genes. That can become brittle. A better test should use a controlled toy ontology where enrichment is guaranteed.

## Opportunities

- Add a resource registry for offline assets, for example `configs/resources.yaml`, listing species, Ensembl release, GO release, local file paths, checksums, source URLs, date acquired, and required columns. This would turn "offline and pinned" into an enforceable contract.
- Define one canonical annotation table format for both human and mouse with columns such as `gene_id`, `gene_id_stripped`, `gene_symbol`, `entrez_id`, `go_id`, `ontology`, `evidence_code`, `resource_release`, and `species`.
- Split resource preparation from analysis. A script such as `scripts/00_prepare_reference_resources.py` could verify local pinned resources and optionally be run manually to create them, while normal pipeline execution only reads verified local files.
- Add explicit species validation: human Ensembl IDs should mostly match `ENSG`, mouse Ensembl IDs should mostly match `ENSMUSG`, with version suffix stripping controlled and logged.
- Add a small mouse fixture with Ensembl-style IDs and toy GO terms. This would validate species branching, annotation parsing, and GO universe behavior without needing a full mouse public RNA-seq dataset.
- Convert the YAML config into a typed model early, using `pydantic` or `jsonschema`. This will reduce implementation drift between agents and make unit tests more direct.
- Add a formal stage manifest schema shared by all scripts. Each manifest should record inputs, outputs, checksums, row counts, package versions, parameters, warnings, and stage status.
- Make R `DESeq2` validation a separate optional validation profile. This avoids mixing routine pipeline acceptance with a heavier external dependency.
- Use deterministic toy statistical fixtures for unit/smoke tests and reserve `airway` for integration validation. For example, create synthetic count matrices with known large treatment effects and enough replicates for stable `PyDESeq2` behavior.
- Add a `--dry-run` option to every CLI or at least the orchestrator, returning planned inputs and outputs without running expensive stages. The pseudocode mentions dry-run testing for the orchestrator but does not include it in the parser.
- Define failure semantics for low-information analyses: no genes after filtering, zero DEGs, no GO terms, too few mapped genes, or all adjusted p-values missing. These should produce valid reports with warnings, not partial crashes.
- Add cross-stage schema tests that read the output of one stage as the input contract of the next stage. This will catch handoff errors between agents early.

## Threats

- `PyDESeq2` API and statistical behavior may drift. Constructor arguments, reference-level handling, exposed size factors, normalized counts, transformed counts, and result columns should be pinned and wrapped behind an adapter with version-specific tests.
- R/Bioconductor validation can create reproducibility friction. Bioconductor versions are tied to R versions, and package installation can be slow or platform-sensitive. Without a pinned container or lockfile, `airway` reference regeneration may not be stable.
- Offline annotation assets can become the largest reproducibility risk. If GO, Ensembl, Entrez mappings, evidence-code filters, and species releases are not pinned together, enrichment results can change while all pipeline code remains unchanged.
- Ensembl mapping can silently bias GO enrichment. Dropping unmapped genes, collapsing many-to-many mappings incorrectly, or using a background universe that differs from the tested genes can inflate or deflate enrichment.
- Simple two-group design is statistically limited. Users may try to include batch, paired samples, or confounded metadata because the metadata contract allows optional covariates. The pipeline must reject unsupported modeling rather than ignoring covariates silently.
- Small-sample RNA-seq tests can be unstable. If v1 allows two samples per group or fewer, p-values, dispersion estimates, and outlier handling can be unreliable. The plan should set minimum replicate guidance and warnings.
- Plot/report generation can mask core failures if not isolated. A plotting dependency or HTML rendering issue should not obscure whether DEG and GO outputs were produced correctly.
- Multi-agent implementation may diverge unless schemas and fixtures are centralized. The plan is detailed, but agents could still make incompatible choices about ID columns, path layout, JSON keys, or warning names.
- Runtime downloads in any stage would violate the fully offline requirement and make analyses non-reproducible in locked environments.

## Prioritized Recommendations

1. Define the offline resource contract before coding GO enrichment. Specify local files, checksums, release versions, Ensembl version stripping, human/mouse ID patterns, mapping behavior, and the exact universe definition.
2. Create a strict config schema and output schema appendix. This should include required YAML fields, defaults, accepted species names, contrast constraints, and required columns for every stage output.
3. Narrow v1 validation gates for `PyDESeq2`. Keep R `DESeq2` airway comparison as optional integration validation, but use robust invariants and pinned versions rather than brittle exact or near-exact agreement assumptions.
4. Add mouse-specific test coverage. A toy mouse Ensembl + GO fixture is enough for v1, but the plan should explicitly state that full biological validation is human-only via `airway`.
5. Add statistical guardrails for simple contrasts. Enforce exactly two levels, minimum replicate warnings or failures, explicit contrast direction, no unsupported batch modeling, and clear messages when metadata contains ignored covariates.
6. Rename or clarify `03_filter_normalize.py` so it cannot be misread as custom normalization before `PyDESeq2`.
7. Centralize manifests and schemas in shared helper modules before implementing stage scripts. This will reduce handoff failures between coding agents.
8. Keep report/HTML features secondary until the core path passes validation: validation, filtering, PyDESeq2, DEG tables, offline GO, and machine-readable manifests should come first.

## Top Review Concerns

1. The offline GO and annotation resource plan is not yet reproducible enough for human + mouse Ensembl workflows.
2. The `PyDESeq2` versus R `DESeq2` airway validation thresholds and assumptions need to be made less brittle and more explicitly version-pinned.
3. Handoff clarity will depend on formal schemas for config, manifests, and stage outputs; examples alone are not enough for multiple coding agents.
