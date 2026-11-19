# Phase 2a — Source provenance and QC strategy audit

**Status:** Pre-QC methodological review, 2026-11-19. **No additional filtering approved.**  
**Companion:** [Phase 2 QC decision register](02_human_qc_plan.md).  
**Evidence convention:** **VERIFIED IN FILE** (local structural/quantitative checks), **SOURCE-DOCUMENTED** (GEO/paper), **INFERRED** (reasonable but needs validation), **OPEN** (unverified). Source reports and papers must not be conflated.

## The crucial distinction: “raw data” has three meanings

1. **Raw sequencer reads** (BCL/FASTQ) permit read-level FastQC/MultiQC, adapter, base quality, mapping, rRNA/intron/exon balance, contamination, duplication/saturation and recounting. **We have not acquired FASTQs for the four cohorts.**
2. **Unnormalized gene counts** (bulk assigned read counts or 10x UMI matrices) are *analysis-ready counts* but **already processed** by alignment, quantification and usually some library/barcode selection. Nonnegative integers do **not** mean untouched reads or that all upstream QC is recoverable.
3. **Author-filtered annotated single-nucleus object** contains a selected set of nuclei, raw-count layer, transformed layers and cell labels. Its raw-count layer is not a raw droplet universe and cannot be used to reconstruct previously excluded nuclei.

**Terminology for manuscript:** “submitter-supplied gene-count matrix,” “NCBI-reprocessed gene-count matrix,” “author-filtered single-nucleus object with preserved RNA counts,” and “Cell Ranger filtered-barcode UMI matrices.” Do not call any of these “raw sequencing data.”

## Source-derived status by cohort

### GSE174409 — NAc and DLPFC bulk

- **VERIFIED IN FILE:** `GSE174409_raw_counts_02102020.csv.gz` is a 60,676-feature × 80-sample integer-count matrix, all GEO sample titles reconciled, Ensembl-like IDs, with 40 samples per region. These values are **not** FASTQs.
- **SOURCE-DOCUMENTED:** The [GEO series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174409) includes both a raw-count table and a **separate filtered log2-CPM table**. Sample records document FastQC, Trimmomatic trimming, HISAT2 alignment, GRCh38/Ensembl reference and counting (one GEO sample describes HiSeq). The [original paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC8463497/) analyzes paired DLPFC and NAc tissue from 20 OUD and 20 comparison donors.
- **INTERPRETATION:** The downloaded file is the *raw-count output of an already processed bulk RNA-seq pipeline*, apparently upstream of the separately supplied filtered log2-CPM matrix. It is not known from integer counts alone whether genes or failed libraries were excluded before deposition; library totals count **assigned gene reads**, not total sequenced reads.
- **QC achievable now:** matrix integrity, count depth, detected genes, distribution, sample similarities, ordinations, region/diagnosis contrasts, RIN/PMI/sex/pH confounding, donor pairing.
- **Not achievable without BAM/FASTQ/run metrics:** read base quality, alignment %, strandedness, insert size, duplication, rRNA fraction, 5′/3′ coverage, gene-body bias, mapping/assignment fraction.
- **Before filtering:** confirm whether raw-count and filtered-log2CPM files correspond by exact sample IDs; verify actual count-generation software/feature type and gene identifier release; check the authors' sample-level exclusions and quality screens. Do **not** impose the authors' filtered-gene list on our unfiltered-count matrix by accident. Sample exclusions should require objective concordant evidence and report retained **people and region pairs**.

### GSE225158 — caudate and putamen snRNA

- **VERIFIED IN FILE:** The acquired `GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5ad.gz` is an **author-derived, annotated, already filtered object**, not unprocessed 10x or FASTQ data. It contains 98,848 retained nuclei from 12 donors (6 OUD, 6 controls), 22 donor × region units; `.X` contains transformed noninteger expression, and `.raw.X` stores 31,611 gene features, with all 365,407,704 sparse values finite, nonnegative and integer-valued. RNA UMI totals and nonzero gene counts agreed exactly for 1,000 sampled nuclei with `nCount_RNA` and `nFeature_RNA`. Gene symbols reside in `raw.var["features"]`, **not** `raw.var["_index"]`. This establishes raw-count *representation*, not original FASTQ provenance.
- **SOURCE-DOCUMENTED:** The [original paper](https://www.nature.com/articles/s41467-024-45165-7) uses STARsolo alignment/counting (human GRCh38.p13), and multiple quality-control approaches: DropletQC nuclear fraction derived from intronic versus exonic UMIs, probabilistic miQC for mitochondrial/detected-gene quality, and scds for doublets. The paper reports removing clusters with >10% of cells flagged as low quality in the integrated embedding. The original authors analyze pseudobulk with donor block (`duplicateCorrelation`), sample quality weights (`voomWithQualityWeights`) and an unmeasured variation strategy. The [GEO sample records for C-612](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7040228) and [C-13291](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7040229) document failed caudate preparations; they are absent from the final object, while the same donors remain in putamen.
- **INTERPRETATION:** **Do not repeat miQC/scds/DropletQC as if this were a fresh unfiltered dataset.** First audit the existing `dropletQC.keep`, `miQC.keep`, `scds.keep`, `percent.mt`, `nCount_RNA` and `nFeature_RNA` against authors' methods and final inclusion logic. Check whether retained nuclei still have failed flags (a flag can be an advisory diagnostic rather than the final selection criterion), and whether `raw.var["features"]` has duplicates or ambiguous gene symbols.
- **QC achievable now:** descriptive **post-author-QC** nucleus distributions, library/region and donor quality heterogeneity, preexisting flag concordance, cell-type proportions, per-donor/region/cell-type support and pseudobulk coverage.
- **Not reconstructable from this file:** prefilter droplet universe, author-excluded nuclei and complete upstream filtering attrition, read-level alignment/intron metrics if absent.
- **Before extra filtering:** quantify how an additional threshold would affect each donor, diagnosis, region and rare complement-expressing cell type. Default is **reuse original author-retained nuclei** unless a reproducible defect is found. `raw.X` is the count aggregation source; neither transformed `.X` nor integrated reductions may be substituted.

### GSE194368 — central amygdala bulk

- **VERIFIED IN FILE:** `GSE194368_raw_counts_GRCh38.p13_NCBI.tsv.gz` contains **41** sample columns and 39,376 numeric Entrez-like gene identifiers, despite 42 GEO records. Omitted `GSM5834154` = O28 = OUD, giving 20 OUD and 21 unaffected samples in the acquired matrix. This is **an NCBI-generated count product**, not the original authors' counts.
- **SOURCE-DOCUMENTED:** The [original paper, Carmack et al., *Corticosteroid Sensitization Drives Opioid Addiction*](https://pmc.ncbi.nlm.nih.gov/articles/PMC10406162/) used rRNA-depleted total RNA sequencing on Illumina NovaSeq. It states FastQC/fastp QC, HISAT2 v2.1.0 alignment to hg38 and featureCounts v2.0.1 with `-s 2 -t transcript -g gene_name` against Ensembl GRCh38.103; hence the original counting units/annotation **differ from our NCBI Entrez-identified matrix**. See [NCBI-generated RNA-seq counts documentation](https://www.ncbi.nlm.nih.gov/geo/info/rnaseqcounts.html) for the separate reprocessing provenance and eligibility rules.
- **INTERPRETATION:** The acquired NCBI matrix is a **valid kind of raw count matrix after independent processing**, but is not interchangeable gene-for-gene or sample-for-sample with the original paper's counts. NCBI-derived assigned count totals are not a measure of original read depth; absence of O28 may reflect run-level preprocessing or eligibility and needs SRA/run-level investigation.
- **QC achievable now:** within-NCBI-matrix integrity, assigned count depth, sample similarity, complexity, data-driven outlier flags, available diagnosis group comparison. Numeric gene IDs require pinned reference mapping.
- **Not achievable without reads/metrics:** FASTQ QC and alignment quality; missing sample count reconstruction; full biological covariate assessment if unavailable in GEO.
- **Before filtering:** find whether a submitter-origin complete count matrix exists (GEO samples refer to “raw gene counts” but some display “supplementary data files not provided”); check O28 SRA/run availability, status and NCBI quantification logs; resolve feature IDs and source annotation; seek source supplemental donor covariates. **Do not mix one original-author column into 41 NCBI-reprocessed columns.** Decide *a priori* whether to use a consistent author matrix, reprocess all 42 from FASTQs, or transparently analyze the consistent 41-sample NCBI matrix.

### GSE240457 — ventral midbrain snRNA

- **VERIFIED IN FILE:** `GSE240457_RAW.tar` stores 27 `filtered_feature_bc_matrix.h5` matrices from Cell Ranger, totaling 269,307 **Cell Ranger-called barcodes** and 36,601 feature rows per library. Integer UMI values are **not raw FASTQ reads**. The archive name “RAW” means GEO supplementary archive of matrices, **not unprocessed sequencing data**.
- **SOURCE-DOCUMENTED:** The [GEO series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE240457) describes 95 donors (45 opioid history/overdose, 50 controls), pooled 3–4 per 10x library. The [original paper](https://www.nature.com/articles/s41467-023-41455-8) describes FANS/DAPI sorted nuclei, genotype-based demuxlet assignment, removal of ambiguous/doublet/unassigned barcodes, an upper mitochondrial filter of 20%, further Seurat processing and a final **212,713** nuclei across 95 donors. It reports very low observed mitochondrial fractions for most nuclei. The authors' [analysis repository](https://github.com/piquelab/sc_brains) and the processed annotated Seurat object are priority sources for missing barcode-level metadata.
- **INTERPRETATION:** The 27 filtered matrices represent **intermediate** count data after Cell Ranger cell calling, but **before** reproduced donor demultiplexing and full original QC. Our 269,307 filtered barcodes and the paper's 212,713 retained nuclei are **different denominators**; their numerical difference is not an independently established count of QC failures because the full stepwise selection overlap has not been proven.
- **QC achievable now:** technical **per-library** UMI/feature/mitochondrial distributions, feature types, unusual barcode profiles and library-to-GEO matching.
- **Not possible responsibly yet:** donor-specific OUD/control comparisons, donor-balanced QC, biological pseudobulk or donor-level cell proportions, or claiming to reproduce the authors' filtered population.
- **Before inference:** obtain/reconcile processed author Seurat object (preferred if usable) or verified `library + barcode -> donor + diagnosis + QC pass + cell type` mapping; record how doublets/ambiguous barcodes were removed and harmonize barcode suffix/prefix conventions. Only then assess whether additional QC is warranted.

## Uniform QC, individualized implementation — the policy

**Yes, uniform conceptual standards. No, not identical filtering algorithms or cutoffs.** A peer-reviewable reanalysis should apply one documented assessment framework to all four studies while accounting for their distinct assay types, count pipelines, denominators and available provenance.

| Common standard (applied to all) | Assay-specific implementation and reason |
|---|---|
| Verify identifiers, dimensions, nonnegative counts, duplicates, mapping and accession provenance | Bulk gene-count matrices vs H5AD `raw.X` vs 10x `features/id/name/type`; different sample units and gene ID namespaces |
| Report sequencing depth/complexity and outliers | Gene-assigned read counts in bulk; UMI + detected genes per nucleus in snRNA; do **not** compare these values as common depth units |
| Check metadata completeness and group balance | GSE174409 multi-region donor pairs; GSE194368 sparse covariates/omitted OUD case; GSE225158 12 donors with 2 failed caudate samples; GSE240457 pooled libraries pre-demux |
| Quantify technical artifacts and flagged observations | FASTQ/BAM quality **not accessible** with present bulk matrices; inherited `miQC`/`scds`/`DropletQC` in striatum; mitochondrial/gene filters and genetic doublets in midbrain original pipeline |
| Prespecify how proposed exclusions affect analysis | Sample donor loss for bulk; donor × region × cell-type pseudobulk support for nuclei; no disease-based filtering without donor labels |
| Preserve original assay and normalized visualization separately | No mixing raw counts, logCPM, SCT/integrated expression or source-inconsistent author-vs-NCBI counts |
| Report sensitivity and independent biological `n` | Within-donor repeated regions and cross-study donor overlap; nuclei/barcodes are not independent replicates |
| Avoid data-dependent DE QC and protect target interpretability | Flag complement-gene low detection and gene-ID loss; do not selectively lower statistical filters only for desired genes |

**Decision boundary:** A shared figure/report format and common acceptance criteria are appropriate; **different diagnostic distributions and appropriate thresholds are necessary**. This is methodological consistency, not analytical inconsistency.

## Two-stage QC model: prevent double-filtering

**Step 1: Data lineage and existing-selection audit (mandatory).** Characterize which QC steps and samples/cells the authors *already* retained or excluded; identify which metrics are intrinsically missing from secondary count products. This step may conclude that additional filtering is unnecessary.

**Step 2: Post-acquisition *descriptive* QC (mandatory), additional exclusions only when justified.** Generate sample/nucleus/library-specific diagnostics and proposed flags. Require corroborating signals or documented technical failure, show donor/group attrition, analyze whether excluded biological material is concentrated in OUD vs controls, preserve an unfiltered-vs-additionally-filtered sensitivity comparison.

## Decision pathways for the four studies

| Dataset | Default baseline for analysis | Trigger for reprocessing or additional exclusions | Gate before inferential testing |
|---|---|---|---|
| GSE174409 | Authors' unnormalized gene counts, authors' deposited samples | Reliable evidence of technical library/sample failure, gene misannotation, or invalid pairing; FASTQ reprocessing only if unresolvable | Verify donor IDs/paired-region structure and plausible covariates |
| GSE225158 | Authors' final 98,848 selected nuclei with verified `raw.X` | Retained failed QC flags with supporting defects; investigate sensitivity rather than indiscriminate refiltering | Confirm gene identities, donor/region/cell-type coverage and original-selection interpretation |
| GSE194368 | One consistent 41-sample NCBI raw-count pipeline, unless 42-sample consistent author matrix becomes available | Clear run QC failure, meaningful pipeline incompatibility, or recovered full matrix with consistent counting | Decide/count origin, explain missing O28, verify Entrez mapping and feasible covariates |
| GSE240457 | **No inferential baseline yet**; 27 Cell Ranger matrices only for technical QC | Author-processed donor annotated Seurat object or complete demux/QC metadata, which may change sample population | Verify barcode-to-donor-to-phenotype and keep/failed flags; only then select analysis baseline |

## Prespecification before running Phase 2 code

We will not hard-code numeric filtering thresholds in the first QC code. Before inspecting OUD-associated DE, write a **versioned decision register** that states:

1. Exact count files / layer choices / checksum and feature reference; original processing steps, source links and date accessed.
2. For every proposed exclusion: metric, why failure is technical (rather than possible biology), direction of effect, comparator/alternative threshold and manual-review status.
3. Both **total retained nuclei** and **number of distinct donors per diagnosis per region/cell type** before and after each candidate rule.
4. Sample/group-level effects, loss of rare complement-expressing cell types and prespecified independently supported minimum for pseudobulk comparisons.
5. Technical/biological confounder identifiability (batch, age, sex, PMI, pH, RIN, geography where available), collinearity and missingness; never invent covariate values.
6. Any deviations from original paper: reproducibility and validity rationale; prohibit seeking desired complement significance by tuning QC.
7. Replication labeling: GSE174409 ↔ GSE225158 overlapping donors cannot be presented as independent cohorts; OUD versus opioid-overdose phenotypes analyzed as different source definitions.

## Outstanding **provenance** questions before production QC implementation

These are not all hard blockers for plotting, but they **must** be triaged and explicitly marked resolved/blocked/deferred in a manuscript-facing record:

- **GSE174409:** verify exact alignment/count software and whether any submitter preprocessing excluded genes/samples from the “raw counts”; source reference GTF/version; author original QC exclusions and pair identifiers. Sample GEO text may not sufficiently describe alignment defaults.
- **GSE225158:** inspect the authors' source code/supplement for the **actual inclusion logic** (including combinations of miQC, scds, DropletQC and cluster removal). Verify donor/cell-type metadata, raw gene feature uniqueness and original sample failure evidence. Check whether 22 observed region units versus paper's “24 biological samples” is explicitly acknowledged.
- **GSE194368:** pursue the original author gene-name featureCount matrix and metadata where possible; identify reason GSM5834154 lacks an NCBI-generated column; confirm how NCBI annotation and authors' Ensembl gene-name counts differ.
- **GSE240457:** prioritize donor annotations from the published Seurat object or repository; determine whether acquiring 56+ GB is unavoidable and whether metadata can be exported separately; clarify actual post-demux/post-QC retained barcode subset.
- **Cross-study:** confirm 12 donor identities from original sample provenance; flag genomic/clinical phenotypes as potentially correlated, and document gene identifier mapping plus comparability limits.

**The go/no-go rule:** Implement **read-only, assay-appropriate metrics** only after the provenance interpretation has been reviewed. Implement **filtering or DE** only after matching the original upstream selection logic and locking a documented, sensitivity-tested analysis design. Raw FASTQs are **not required by default** for a published count-matrix reanalysis; their acquisition should be a justified exception when missing upstream quality evidence materially threatens validity.
