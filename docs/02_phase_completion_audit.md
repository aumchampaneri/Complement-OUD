# Phase 2 QC completion and consolidation audit

**Status: diagnostic implementation complete; Phase 2 scientific decision freeze NOT YET complete.**
Reviewed against `docs/02_human_qc_plan.md`, `docs/02a_source_provenance_audit.md`
and `docs/02b_qc_preflight.md`. This is a source- and user-output-based
audit, not a claim that the final refactored scripts have passed local execution.
All inputs, temporary extraction, annotation snapshots, generated outputs and
source checksums remain outside version control.

## Implemented and supported by observed outputs

| Goal | Evidence | Status |
| --- | --- | --- |
| Acquire/inspect all four human datasets, preserve source provenance | `01a`, `_inspect_human_data.py`, SHA-256 inventory | COMPLETE for deposited files; no FASTQ-level QC |
| Reconcile samples and donor units without pseudoreplication | `01b`, `02a`, `02d`, `02e` | COMPLETE for available metadata; midbrain donors unresolved |
| Audit bulk matrix integrity and outlier structure | `02a`, `02b`; 80 + 41 samples | COMPLETE diagnostic pass |
| Audit bulk IDs and provisional mitochondrial expression | `02c`; 59,320/60,676 Ensembl-style IDs; 37,660/39,376 Entrez-style IDs; 37 MT loci each | COMPLETE exploratory audit; reference compatibility OPEN |
| Audit inherited striatal single-nucleus quality | `02d`; 98,848 nuclei, 12 donors, 22 donor–region samples; 190 donor-region-celltype groups | COMPLETE structural pass; flag semantics OPEN |
| Audit 27 pooled midbrain 10x libraries | `02e`; 269,307 Cell Ranger filtered barcodes | COMPLETE technical pass; donor-level analysis BLOCKED |
| Preserve source counts; no automatic exclusions, pseudobulk, DE or pathway-activation claims | Phase 2 readiness reports and scripts | COMPLETE for QC stage |
| Unit tests for assay-specific checks | User reports 02a–02e tests passing before shared-helper refactor | PREVIOUSLY PASSED; rerun full suite after consolidation |

## Still required before QC freeze and Phase 3

1. **Inherited striatal flags:** inspect refreshed `striatum_qc_readiness.json`
   key `qc_flag_original_value_counts`; determine encoding and source-backed
   meaning of `dropletQC.keep`, `miQC.keep`, `scds.keep`. In prior uploaded
   output, all 98,848 values per flag appeared as `UNKNOWN` under our boolean
   parser. Do not equate UNKNOWN with failed or passed QC.
2. **Bulk technical anomalies:** record GSE194368 high mitochondrial fractions,
   especially C.18, and assess the role of library-size denominators and
   composition in multivariate profiles. No sample exclusion approved.
3. **Gene identity/source compatibility:** confirm whether GENCODE v47 and the
   captured NCBI `gene_info` snapshot correctly cover count-generation
   identifiers; resolve unmapped features and duplicated gene symbols before
   symbol-level aggregation or cross-study harmonization. Preserve the SHA-256s
   recorded in `bulk_annotation_readiness.json`.
4. **Striatal donor-celltype inclusion policy:** document counts by diagnosis,
   independent donor support, low-representation groups, anticipated donor loss
   under any proposed pseudobulk support rule, and region-specific model
   feasibility. Do not select cutoffs based on disease results.
5. **Midbrain inference gate:** recover verified 10x barcode-to-donor-to-diagnosis
   crosswalk and authors' filtering/demultiplexing provenance. Until then,
   neither 27 pools nor 269,307 barcodes are independent OUD replicates.
6. **Study designs and cross-study dependence:** source-verify GSE174409 donor
   pairs and 12 provisional cross-study GSE174409/GSE225158 donor matches;
   consider region effects, missing GSE194368 covariates, differing OUD/overdose
   clinical definitions, sample sizes, and candidate covariate model ranks.
7. **Versioned final deliverables:** the original plan called for a reviewed
   `qc_review.md`, `qc_decisions.csv` (or equivalent), hypothetical
   donor/feature attrition and sensitivity results, and an
   `analysis_design_audit.json`. These have **not** yet been produced or
   approved. Avoid labeling them complete.
8. **QC visualizations and method correspondence:** plots for study/region/
   donor-level metrics and sensitivity summaries remain to be formalized.
   Original plan mentioned Spearman correlations, whereas `02b` currently
   computes **Pearson** correlations on log1p-CPM selected genes. Either
   justify/document Pearson explicitly or add predeclared Spearman diagnostics;
   do not silently claim they are interchangeable.

## Refactor decisions

The scientific scripts remain five distinct entrypoints: `02a`–`02e`.
The small private helper `scripts/_qc_common.py` now owns repeat CSV I/O
and SHA-256 implementation. Method-specific streaming, AnnData metadata
decoding, 10x CSC sparse access, and feature provenance checks are retained in
their respective files; combining them would obscure separate acceptance
criteria. `01a` and `01b` remain standalone acquisition/validation.
`_inspect_human_data.py` **must not yet be deleted** because source H5AD
integer-value verification, initial archive inspection and cross-study
candidate matching have not all been replaced by maintained entrypoints.

The `02c` PCA-only mitochondrial sensitivity check still re-reads the
input count matrix twice; optimizing this would be a *substantive* algorithmic
change, not mere dead-code removal. Leave until parity tests and design
decisions exist.

## Final local regression and repository checks

From the repository root:

```bash
git pull origin architecture/refactor-v1
pixi run python -m unittest discover -s tests -p 'test_*.py' -v
pixi run python -m compileall -q scripts tests
git status
```

If the tests pass, inspect output parity before rerunning large jobs:
`02a/02b/02c` write the same schema/diagnostics after shared helper changes;
`02d/02e` produce the same schema. The I/O refactor is intended to be
behavior-preserving. If a local test fails, report the traceback before changing
any original source data.

**Gate:** Phase 3 may begin only when the relevant cohort's donor-aware
analysis design and QC decision register are approved, with unavailable
assays explicitly blocked or deferred.
