# Phase 2 implementation preflight — version 1

**Status:** Read-only diagnostic implementation started. No threshold-based exclusions, normalization for inference, pseudobulk, pathway testing or DE.

## Prespecified input datasets

| Cohort | Expected matrix | Count origin | Exact expected shape | Sample identifiers | Biological replication |
|---|---|---|---|---|---|
| GSE174409 | `data/raw/human/GSE174409/GSE174409_raw_counts_02102020.csv.gz` | Author-supplied gene counts | 60,676 features × 80 tissue samples | GEO title (HN/HBA) | ~40 donors across two regions; exact pairing under verification |
| GSE194368 | `data/raw/human/GSE194368/GSE194368_raw_counts_GRCh38.p13_NCBI.tsv.gz` | NCBI-generated count matrix | 39,376 features × 41 samples | GSM accession | 41 postmortem tissue samples; missing GSM5834154 must remain documented |
| GSE225158 | Author-processed H5AD `.raw.X` | Raw-style integer counts in author-selected nuclei | 98,848 nuclei × 31,611 genes | `Case` donor; `ID` donor–region | 12 donors and 22 donor–region samples |
| GSE240457 | 27 Cell Ranger filtered-feature H5 matrices | Quantified UMI after initial barcode calling | 269,307 barcodes total | GSM pool/library; not donor | Donor identity unavailable in the acquired matrices |

**Inputs for first QC script:** original deposited matrices + `outputs/inspection/GSE174409_inspection.json` and `GSE194368_inspection.json`. The script validates column headers, gene counts, integer counts, matrix dimensions and exact agreement with the preceding inspection, and records **SHA-256** of each original matrix.

## Phase 2a: bulk QC without exclusions

Run `pixi run python scripts/02a_qc_bulk.py` after `pixi run python -m unittest discover -s tests -p 'test_qc_bulk.py' -v`.

The script streams each compressed count matrix separately without building a dense matrix and writes **only**:

- `outputs/qc/human/bulk_sample_qc.csv`: per sample total **assigned** gene counts (not sequencing read depth), number of detected genes, zero-feature fraction, fraction of counts in top ten counted features, and exploratory robust outlier flags.
- `outputs/qc/human/bulk_qc_readiness.json`: exact input source and checksum, cohort sizes, limitations and explicit `filtering_performed=false`.

**Outlier flags:** median/MAD-based robust z scores of log1p assigned counts and detected-feature counts, within region for GSE174409 and within study for GSE194368. Threshold `|z| >= 3.5` is *only a flag for manual review*, **not** a statistical exclusion criterion or disease contrast. Fewer than five samples or zero MAD means no interpretable flag. A sample cannot be removed on the basis of any of these flags automatically.

**Why no single QC threshold across datasets:** sequencing platform and quantification are different; count values in bulk are assigned reads, single-nucleus measurements are UMIs, and Cell Ranger filtered barcodes are not author-retained nuclei. Common principles apply (integrity, metadata concordance, outlier review, donor/group attrition, provenance); exact measurements and any eventual decisions differ by assay.

**Deliberately deferred until we inspect these baseline results and source methods:** PCA, correlation matrices, variance stabilization, gene identifier mapping, sample exclusions, removal of lowly expressed complement genes, joint bulk cross-study analysis, and all phenotype-driven inference. No arbitrary removal or normalization is used in the first pass.

## Gating and audit trail

1. **First run:** deliver `bulk_sample_qc.csv` and `bulk_qc_readiness.json`; check agreement with preflight dimensions.
2. **Source review:** assess original author QC, model covariate completeness, paired donor mapping, NCBI-vs-author provenance, missing O28, complement gene ID mapping.
3. **Next descriptive run:** implement nucleus QC using *existing* GSE225158 flags and the author-retained H5AD population, and pooled midbrain library QC without pretending libraries are individual donors.
4. **Before any filtering or inferential DE:** freeze proposed thresholds, sample and donor attrition, sensitivity comparisons, and condition/covariate balance in a versioned decision record.

Detailed rationale: [Source provenance audit](02a_source_provenance_audit.md) and [Phase 2 QC plan](02_human_qc_plan.md).
