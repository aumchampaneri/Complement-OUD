# Phase 2: Human transcriptomics QC design and decision register

**Status:** Planning only — no additional QC thresholds, exclusions, transformations, integration, or differential expression have been applied by this phase.  
**Scope:** Four acquired human GEO series; mouse experiments excluded from this phase.  
**Inputs:** Raw acquisition files under `data/raw/human/`, GEO metadata, `outputs/inspection/`, and `outputs/validation/human/`.  
**Output policy:** Preserve immutable raw inputs; exploratory plots and metrics go into ignored `outputs/qc/human/`; version-control methods, decision records, and code.  
**Updated:** 2026-11-19.

**Required companion before implementation:** [Source provenance and QC strategy audit](02a_source_provenance_audit.md). It distinguishes sequencer reads, count-generation pipelines, original upstream QC, and author-retained objects. No universal numeric cutoffs should be applied until those provenance issues are triaged.

## 1. Scientific purpose and inferential boundaries

This phase asks whether each human expression dataset supports reliable **donor-level** comparisons of complement-associated transcription in OUD or opioid overdose, and what quality, sampling, or technical effects must be accounted for. It does **not** establish complement pathway activation at the protein, cleavage, deposition, or functional level.

The decision hierarchy is:
1. Validate **data provenance and biological sample identities** before calculating inferential statistics.
2. Describe technical metrics **without filtering**.
3. Evaluate sample/donor/region/cell-type representation and potential confounding.
4. Record candidate thresholds, their rationale, and sensitivity-analysis alternatives.
5. Freeze a versioned, donor-aware inclusion and modeling specification *before* DE testing.

The **unit of biological replication is an individual donor**. A brain region from the same donor is repeated tissue, not a new independent person. A cell/nucleus is a nested observation, not an independent replicate. For multi-donor 10x pools, a sequencing library is not equivalent to a donor.

## 2. Data inventory, origins, and what each representation means

| GEO dataset | Tissue and biological design | Acquired representation / source | Verified structural facts | Interpretation and caveats |
|---|---|---|---|---|
| **GSE174409** | Human postmortem nucleus accumbens (NAc) and dorsolateral prefrontal cortex (DLPFC); GEO describes 20 OUD and 20 controls with both regions | **Submitter-supplied bulk RNA-seq gene-count CSV:** `GSE174409_raw_counts_02102020.csv.gz`; GEO series matrix and family SOFT provide sample characteristics | 60,676 Ensembl-style gene IDs × 80 sample columns; 80/80 matched to GEO; nonnegative integer-like gene counts in inspection | A sample is a tissue-region measurement; pairing across regions must be represented by donor ID; gene identifiers require annotation/version checks; bulk expression conflates cell abundance and within-cell expression |
| **GSE225158** | Human postmortem caudate and putamen; 12 donors (6 OUD/6 controls) with 24 GEO donor–region accessions, 22 included in processed object | **Author-derived, filtered/annotated AnnData H5AD** converted from a Seurat object: `GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5ad.gz` (plus GEO series metadata). **`.raw.X`** is count-like, whereas **`.X`** is transformed | 98,848 nuclei × 31,611 raw features, 365,407,704 stored sparse entries fully checked (no negative, nonfinite or fractional values); all 1,000 nuclei tested matched `nCount_RNA` and `nFeature_RNA` exactly; `raw.var["features"]` stores gene symbols; 22 donor–region IDs | `raw.X` is technically appropriate for count aggregation, subject to upstream provenance; it is **not unfiltered droplet data**. `C-612` and `C-13291` were excluded after reported wetting/low-quality failures; these are documented exclusions, not missing unexplained cases. `Case` identifies donor, `ID` identifies donor–region. `raw.var["_index"]` is a numeric-like index, **not** gene symbols |
| **GSE194368** | Human postmortem central amygdala, GEO 21 OUD and 21 unaffected controls | **NCBI-generated GRCh38.p13 bulk gene counts**: `GSE194368_raw_counts_GRCh38.p13_NCBI.tsv.gz`, not the authors' original count-generation pipeline; GEO metadata | 39,376 features × 41 GSM-labeled sample columns; all 41 match GEO; `GSM5834154` / `O28` omitted from NCBI matrix; identifiers resemble numeric Entrez Gene IDs | Actual matrix is 20 OUD + 21 unaffected pending QC, not the nominal 21 + 21. Gene IDs must be mapped against a pinned NCBI annotation with ambiguous/missing mappings reported. Don't mix submitter counts with NCBI counts one column at a time |
| **GSE240457** | Human postmortem ventral midbrain; publication reports 95 donors (45 opioid overdose/abuse, 50 opioid-free control) | **27 Cell Ranger filtered_feature_bc_matrix.h5 files** bundled in `GSE240457_RAW.tar`, plus GEO library metadata; published demultiplexed/annotated Seurat RDS was not acquired | 27/27 H5 libraries structurally inspected, ~36,601 10x feature rows/library, sum of 269,307 filtered barcodes | Libraries contain pools of ~3–4 donors; barcodes are **not yet demultiplexed**, and do not carry donor/clinical-group assignments in our local data. Publication's ~212,713 retained nuclei are a different, processed population. Donor-aware OUD DE **blocked** until a verified barcode-to-donor/diagnosis mapping is obtained |

**Sources and recorded provenance:**
- The four local `GSE*_inspection.json` reports, `human_reconciliation.json`, `GSE225158_raw_verification.json`, `dataset_readiness.json`, and acquisition manifest are primary *local inspection evidence*. Regenerate reports if source data change.
- GSE174409: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174409
- GSE225158: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225158 ; original paper/methods: https://www.nature.com/articles/s41467-024-45165-7 ; source analysis: https://github.com/pfenninglab/Logan_Striatum_snRNA-seq
- Documented caudate exclusions: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7040228 and https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7040229
- GSE194368: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194368
- GSE240457: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE240457 ; original methods: https://www.nature.com/articles/s41467-023-41455-8 ; author's analysis code: https://github.com/piquelab/sc_brains
- Exact gene annotation versions, origin of the converted H5AD, complete original processing chain, and individual donor-to-annotation links require explicit source verification; they are **not proven** merely by integer count checks.

## 3. Clarify levels, provenance, and repeated measurements

| Field or entity | Meaning | Required action before modeling |
|---|---|---|
| GEO GSM accession | Submitted sequencing/tissue sample or pooled library record (depends on study) | Never assume GSM = donor across all studies |
| GSE174409 HN/HBA/... column | Bulk tissue sample | Match GEO accession; derive/verify a stable donor ID and region pairing; don't infer cross-study identity from label syntax alone |
| GSE225158 `Case` | Biological donor (numeric-looking ID) | Canonicalize safely as identifier, **not** quantitative covariate |
| GSE225158 `ID` | Donor + region (e.g., `C-1366`, `P-1366`) | Use to form donor × region × cell type aggregation and QC strata |
| GSE225158 `orig.ident` | Original library/object identifier, sometimes with suffix | Cross-check against `ID`; do not mistake suffixes for new donors |
| GSE240457 10x library | Multiplexed pool of individual donors | Require barcode-level demuxlet or other validated assignments; preserve library and donor in metadata |
| All studies: OUD / control labels | Different source-specific case definitions | Document phenotype contrasts separately: diagnosed OUD vs unaffected is not identical to opioid overdose vs opioid-free |
| GSE174409 × GSE225158 matching | 12 strong candidate shared donors, 24 bulk sample matches | Keep `CANDIDATE_NOT_VERIFIED` until source confirmation; never count these cohorts as independent replication |
| GSE194368 `GSM5834154` | One GEO OUD sample absent from NCBI count matrix | Track as missing; do not fill with zeros or silently generate a mixed-origin count column |

## 4. The pre-filter QC workflow

**Stage A — Source and identity audit (no statistical filtering).**
- Snapshot accession, file names, SHA-256 checksums and size, assay type, count origin, reference/annotation where documented, matrix dimensions, and count validation results.
- Assert exact one-to-one bulk sample-column ↔ GEO mapping and check duplicate gene IDs. Check design labels and region pairing, conflicting diagnoses, missing metadata, and missing samples.
- Verify raw counts (integer/nonnegative/finite), sample order, gene identity and uniqueness, total/library-size consistency. For sparse matrices distinguish stored entries from *nonzero* features.
- H5AD: retain original `obs`, `raw.var["features"]`, `raw.X`, transformed `X`, and preexisting QC fields as distinct representations, without overwriting them. Check unique names and gene symbol duplicates before pseudobulk.
- 10x H5: record genome build/reference, `features/id`, `features/name`, `features/feature_type`, and each library's barcode space; only `Gene Expression` features should enter RNA expression calculations unless explicitly justified.

**Stage B — Descriptive bulk QC (GSE174409 and GSE194368 separately).**
- Per sample: total assigned counts, number of detected genes, count distribution, zero fraction, library-size and complexity outliers; flag—not immediately discard—unexpectedly small or large profiles.
- Sample–sample Spearman correlation and ordination using an exploratory transformation appropriate to count data. Normalize for *visualization* with a documented method; raw counts must remain unchanged for any later negative-binomial DE.
- Plot QC and ordination against diagnosis, region, age, sex, PMI, RIN, brain pH, storage time, and available batch variables. **Missing variables are reported; never fabricated.**
- GSE174409: investigate 40-person pairing across NAc/DLPFC and show region-separated as well as jointly labeled plots. Region may dominate PCA; that is expected, not an automatic batch effect.
- GSE194368: sample composition is asymmetric (20 OUD/21 control) because of missing O28; assess count origin, numeric gene ID mapping and availability of covariates before choosing a regression design.
- Determine whether **any** sample is anomalous in several independent metrics, not just one z-score. Retain a proposed outlier list with rationale and before/after sensitivity plan.

**Stage C — Descriptive single-nucleus QC (GSE225158).**
- Use `raw.X` counts for direct count metrics; show distributions of UMI counts, detected genes, `percent.mt`, mitochondrial and ribosomal signal, doublet metrics (`scds.*`), droplet QC (`dropletQC.*`), mixture QC (`miQC.*`), and original annotation confidence where present.
- Important: nuclear RNA often exhibits **lower mitochondrial RNA**, different gene complexity and substantial intronic counts versus whole-cell scRNA-seq. Broad default scRNA thresholds are not prespecified as universal cutoffs.
- **Inherited selection bias:** the supplied object already excludes two failed caudate samples and may contain prior cell filtering, doublet removal and integration; we can audit included cells, but cannot reconstruct discarded droplets from this H5AD.
- Stratify all nucleus QC metrics by `ID`, `Case`, region, diagnosis and `celltype1/2/3`. Report QC flag concordance and disagreement, not one naive global threshold.
- Create donor × region × cell-type tables with number of nuclei, raw summed UMIs, per-nucleus distributions, number of expressing nuclei for predefined complement genes, and within-donor coverage; inspect class balance before selecting pseudobulk strata.
- Assess cell-type composition variability; distinguish true donor-level abundance shifts from technical depletion. Oligodendrocyte-heavy cell composition may dominate any pooled expression average.
- Save a *proposed* nucleus exclusion mask only after plots; compare with original author-filtering flags. Do not quietly refilter the entire object or assume that integrated/normalized values are counts.

**Stage D — Limited descriptive 10x QC (GSE240457), without individual-level inference.**
- Per library, compute filtered barcode number, total UMI counts, detected genes and mitochondrial proportions, and distributions across barcodes. Confirm gene-expression feature selection.
- Check for anomalous library characteristics and unequal pooling; preserve the mapping `GSM -> library -> barcode`.
- Treat this as **library-level technical** QC only. No diagnosis stratification, donor-level composition analysis or pseudobulk until author-provided demultiplexed barcode metadata are available.
- Do not compare 269,307 filtered Cell Ranger barcodes with 212,713 post-demultiplexing nuclei as though they are identical QC populations. Document how each selection step produces the denominator.

**Stage E — Summarize, specify alternatives and freeze decisions.**
- Produce QC tables, plots, threshold sensitivity summaries, ambiguous-gene mapping tables, proposed exclusions and a manually reviewed `qc_decisions.yaml` or equivalent decision table.
- Dataset, region and cell-type specific decision thresholds are chosen **only after observing the distributions and anticipated donor loss**; include retained donor counts and condition balance after each hypothetical exclusion.
- No differential-expression testing or pathway selection is allowed to drive QC choices (avoid outcome-dependent threshold tuning).

## 5. Choices with important downstream consequences

| Decision | Benefits | Risks and downstream implications | Current policy |
|---|---|---|---|
| Gene filtering by total abundance | Lowers multiplicity and unstable estimates | Complement genes can be sparse, cell-type restricted; overly stringent global cutoffs erase biologically relevant targets | Assess filter sensitivity and donor-level detection; prespecify gene-set reporting, do not exempt genes from valid statistical model requirements post hoc |
| Re-filtering nuclei | Can eliminate low-quality/doublets | Overfiltering may remove affected cell types/donors or OUD-related states; original input is already filtered | Audit inherited QC first; use sensitivity analyses if extra filters justified |
| Mitochondrial thresholds | May identify damaged nuclei | Nuclear assay has distinct mitochondrial distribution; excessive filtering is systematic cell-type bias | Plot per donor/region/cell type; do not transplant arbitrary whole-cell thresholds |
| Doublet detection and filtering | Avoids artificial coexpression | Original flags may already embody filtering; additive filters can disproportionately remove high-RNA cells | Compare `scds` and published selection; retain provenance |
| Pseudobulk aggregation | Restores valid donor-level biological replication | Rare cell types produce unstable donor profiles; aggregating over region loses anatomy | Use raw counts, group by donor × region × cell type; require prespecified minimum donor and nucleus support |
| Rare cell types | Include biologically important immune/perivascular populations | Very few nuclei/donors may prevent meaningful independent inference | Declare primary vs exploratory strata by donor coverage, not cell count alone |
| Normalization for visualization vs DE | Controls library-size effects | Log-normalized cells are not raw counts; cross-platform differences complicate direct pooling | Store separate transformations; do not use integrated `.X` in count-based DE |
| Batch adjustment or integration | May reduce technical variation | Can erase genuine disease/region effects if batch/diagnosis confounded; integrated values invalidate raw-count DE | Compare metadata balance and PCA first; do not auto-integrate assays/cohorts |
| Sample/donor covariates | Reduces confounding | 12-donor striatum has limited degrees of freedom; excessive adjustment and collinearity destabilize models | Check confounding, missingness and design-matrix rank before choosing parsimonious covariates |
| Repeated regions within donor | Enables within-person comparisons | Treating each region as independent inflates n and significance | Region-specific donor-aware analyses first; paired/mixed designs where justified |
| Cross-study harmonization | Allows concordance tests | Ensembl, gene symbols and numeric NCBI IDs may map many-to-one; bulk vs nuclei and region are not interchangeable | Use pinned annotation mapping, track one-to-many/ambiguous cases, compare study-specific effect estimates |
| Combining GSE174409 and GSE225158 | Links bulk tissue and cell types | Likely shared donors; double-counting as independent validation overstates replication | Treat as complementary modalities, not independent replication |
| OUD vs overdose phenotypes | Widens brain-region coverage | Clinical definition, cause of death, comorbid substances and tissue ascertainment differ | Report within-study contrasts and heterogeneity rather than indiscriminate meta-pooling |
| Complement pathway inference | Focuses scientific hypothesis | mRNA, cell composition, complement protein cleavage and activation are different measurements | Phrase as complement-related transcription/signatures, not activated complement without orthogonal evidence |

## 6. Proposed QC deliverables and acceptance gates

Files will be created in ignored `outputs/qc/human/`; the following filenames are **planned**, not yet produced.

1. `source_data_inventory.csv`: accession, assay, matrix layer, checksum/provenance, reference/feature ID type, sample/donor granularity.
2. `bulk_sample_qc.csv` and bulk diagnostic plots: counts, detected features, zero fraction, per-cohort PCA/correlation/outlier flags.
3. `striatum_nucleus_qc_summary.csv`: per-donor/region/cell-type QC and inherited flags; a donor × cell-type coverage matrix and report of retained versus absent regions.
4. `midbrain_library_qc.csv`: per-library technical distributions; no inferred diagnosis.
5. `feature_identifier_audit.csv`: raw Ensembl IDs vs symbols vs Entrez-style IDs; unmapped, duplicated, and ambiguous identifiers; predefined complement target coverage.
6. `qc_review.md` and `qc_decisions.csv`: alternatives, rationale, hypothetical retained donors and cells, and sensitivity tests.
7. `analysis_design_audit.json`: condition/region balance, missingness, candidate model rank, independent donors and covariate feasibility, plus explicit blockers.

**Acceptance gates before any inferential DE:**
- **GSE174409:** identify region-specific sample donor mapping, check condition/covariate balance and any outlier exclusions; record handling of paired regions and shared donors.
- **GSE225158:** use `raw.X` and `raw.var["features"]` only after unique feature mapping; prove donor/region/diagnosis labels consistent; evaluate cell-type donor coverage; document the two upstream failed libraries.
- **GSE194368:** retain the missing GSM as a documented omission; validate feature mapping and test what covariates and count provenance can actually support.
- **GSE240457:** do **not** pass donor-level DE gate without validated `barcode -> donor -> diagnosis` plus library/pool mappings, even if descriptive technical QC looks excellent.
- **Across all datasets:** preserve raw-source checksums, use donor-level uncertainty, avoid expression-driven QC decisions, and limit claims to measured RNA.

## 7. Open questions and evidence needed (not QC code blockers)

- Is there a published source donor crosswalk establishing identity across GSE174409 and GSE225158? All 12 current matches are strong but provisional.
- Can original submitter counts for GSE194368 be recovered, and is GSM5834154 omitted by design or NCBI processing? Which demographic and postmortem covariates are available?
- Can the GSE240457 authors' metadata (or a processed object) provide barcode/donor/clinical status mapping without downloading 56.7 GB? Do sample/technical batches confound the published phenotype?
- What processing and gene annotation versions generated each matrix? Note that `raw.X` checks establish technical count consistency, **not a full processing audit**.
- Which complement genes have sufficient expression in donor-level microglia, astrocytes, oligodendrocytes, neurons, and rare vascular cells to support a prespecified inferential panel?

## 8. Implementation sequence after plan approval

1. Implement lightweight **bulk QC** first: no heavy dependencies, summary distributions and transformed exploratory PCA/correlation; keep separate cohorts and stable samples.
2. Add **GSE225158 descriptive QC** using a streaming/sparse approach to avoid expensive densification; defer reprocessing/integration.
3. Implement **GSE240457 H5 streaming QC** from the tarred per-library matrices, clearly limited to technical/library-level evidence.
4. Review diagnostic outputs together and only then decide exclusion thresholds, pseudobulk minimum support, regression covariates and any batch strategy.
5. Begin human-only complement-focused inference when the appropriate dataset-specific acceptance gates pass.

**Deliberate non-goals in Phase 2:** automated sample removal, cell-type reannotation by default, uncritical normalization of `raw.X`, auto-integration across assays, scVI training, cell-level case/control p-values, inferred donor assignment from pooled library IDs, and early mouse analysis.
