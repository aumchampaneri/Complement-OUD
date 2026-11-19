# %% Human sample and donor-design validation
"""Reconcile GEO inspection reports without opening large expression matrices.

Run: pixi run python scripts/01b_validate_human.py
Input: outputs/inspection/GSE*_inspection.json and GSE225158_raw_verification.json
Output: outputs/validation/human/ (gitignored)
The source datasets and inspection JSONs are never modified.
"""

# %% PATHS AND IMPORTS
from pathlib import Path
import csv
import json

ROOT = Path(__file__).resolve().parents[1]
IN = ROOT / "outputs" / "inspection"
OUT = ROOT / "outputs" / "validation" / "human"
ACCESSIONS = ("GSE174409", "GSE225158", "GSE194368", "GSE240457")


# %% FILE HELPERS
def load_json(path):
    if not path.is_file():
        raise FileNotFoundError(f"Missing {path}; run _inspect_human_data.py first")
    return json.loads(path.read_text(encoding="utf-8"))


def write_csv(name, rows, fields):
    with (OUT / name).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def issue(dataset, code, observation, next_step, priority="review"):
    return dict(dataset=dataset, code=code, priority=priority,
                observation=observation, next_step=next_step,
                status="unresolved")


# %% BULK RNA-SEQ SAMPLE AUDIT
def bulk_samples(accession, report, issues):
    geo = report["geo"]["samples"]
    matrix = report["expression"][0]
    cols = matrix["sample_names"]
    by_title = {x["title"]: x for x in geo}
    by_gsm = {x["geo_accession"]: x for x in geo}
    use_gsm = len(set(cols) & set(by_gsm)) > len(set(cols) & set(by_title))
    lookup = by_gsm if use_gsm else by_title
    expression_qc = {r["sample"]: r for r in matrix["per_sample"]}
    output = []
    for col in cols:
        s = lookup.get(col, {})
        qc = expression_qc.get(col, {})
        output.append(dict(dataset=accession, matrix_column=col,
                           gsm=s.get("geo_accession", ""),
                           geo_title=s.get("title", ""),
                           donor_id="", donor_status="not_verified",
                           region=s.get("source_name_ch1", s.get("brain_region", "")),
                           condition=s.get("diagnosis", s.get("disease_state", "")),
                           age=s.get("age", ""), sex=s.get("sex", ""), pmi=s.get("pmi", ""),
                           library_size=qc.get("library_size", ""),
                           genes_detected=qc.get("genes_detected", ""),
                           matched_to_geo=bool(s)))
    absent = sorted(set(lookup) - set(cols))
    extra = sorted(set(cols) - set(lookup))
    if absent:
        issues.append(issue(accession, "geo_missing_from_matrix",
                            f"Count matrix omits {len(absent)} GEO samples: {absent}",
                            "Investigate original author matrix or NCBI processing/exclusions"))
    if extra:
        issues.append(issue(accession, "unmatched_matrix_samples", str(extra),
                            "Resolve GEO/matrix column mapping before DE", "block"))
    for metric in ("negative_entries", "missing_entries", "malformed_rows",
                   "noninteger_or_nonnumeric_entries"):
        if matrix.get(metric):
            issues.append(issue(accession, metric, str(matrix[metric]),
                                "Investigate count matrix integrity", "block"))
    return output, dict(matrix_samples=len(cols), geo_samples=len(geo),
                        missing_geo_keys=absent,
                        missing_gsms=[s.get("geo_accession") for s in geo
                                      if s.get("geo_accession" if use_gsm else "title") in absent],
                        genes=matrix.get("genes"), source=matrix.get("file"))


# %% STRIATAL DONOR–REGION AUDIT
def striatum_samples(report, issues):
    geo = report["geo"]["samples"]
    expr = report["expression"]
    observed = dict(expr["metadata_distributions"]["ID"]["top_categories"])
    output = []
    for s in geo:
        nuclei = observed.get(s["title"], 0)
        output.append(dict(gsm=s.get("geo_accession"), sample_id=s["title"],
                           subject_id=s.get("subject"), region=s.get("tissue"),
                           condition=s.get("dx_oud"), age=s.get("age"),
                           sex=s.get("sex"), pmi=s.get("pmi"), nuclei=nuclei,
                           included_in_h5ad=bool(nuclei)))
    missing = sorted({s["title"] for s in geo} - set(observed))
    unrecognized = sorted(set(observed) - {s["title"] for s in geo})
    documented = {"C-13291", "C-612"}
    undocumented = sorted(set(missing) - documented)
    if undocumented:
        issues.append(issue("GSE225158", "unexpected_donor_region_exclusions",
                            f"Absent from processed object without documented reason: {', '.join(undocumented)}",
                            "Investigate unexpected sample exclusions", "block"))
    documented_present = sorted(set(missing) & documented)
    if unrecognized:
        issues.append(issue("GSE225158", "unknown_object_samples", str(unrecognized),
                            "Resolve sample ID mismatch", "block"))
    per_donor = {}
    for s in geo:
        donor = s.get("subject", "")
        per_donor.setdefault(donor, set()).add(s.get("dx_oud", ""))
    conflicts = sorted(k for k, groups in per_donor.items() if len(groups) > 1)
    if conflicts:
        issues.append(issue("GSE225158", "subject_diagnosis_conflict", str(conflicts),
                            "Reconcile subject diagnosis labels", "block"))
    return output, dict(geo_samples=len(geo), object_samples=len(observed),
                        donors=len(per_donor), missing_sample_ids=missing,
                        documented_exclusions=documented_present,
                        documented_exclusion_reason=("GEO sample records: wetting failures; "
                            "low mean UMI/cell and anomalously high estimated cell count. "
                            "See GSM7040228 and GSM7040229."),
                        nuclei=sum(observed.values()),
                        raw_count_dimensions=expr.get("raw/X", {}).get("shape"))


# %% MIDBRAIN LIBRARY INVENTORY
def midbrain_libraries(report, issues):
    geo = {s["geo_accession"]: s for s in report["geo"]["samples"]}
    output = []
    for lib in report["expression"]["libraries"]:
        filename = lib.get("filename", "")
        gsm = filename.split("_", 1)[0]
        shape = lib.get("matrix_shape", [None, None])
        output.append(dict(gsm=gsm, library=geo.get(gsm, {}).get("title", ""),
                           h5_filename=filename, genes=shape[0],
                           filtered_barcodes=lib.get("barcode_count", ""),
                           donor_id="", diagnosis="", donor_mapping="NOT_AVAILABLE"))
    issues.append(issue("GSE240457", "donor_demultiplexing_required",
                        "27 pooled 10x libraries do not identify biological donors",
                        "Recover author barcode-to-donor and clinical metadata; "
                        "do not use library as a donor", "block"))
    return output, dict(libraries=len(output),
                        filtered_barcodes=sum(int(r["filtered_barcodes"]) for r in output),
                        donors_with_verified_diagnosis=0)


# %% RAW-COUNT CONFIRMATION
def raw_count_verification(issues):
    path = IN / "GSE225158_raw_verification.json"
    if not path.exists():
        issues.append(issue("GSE225158", "raw_verification_missing",
                            "Corrected validation report not found",
                            "Run _inspect_human_data.py --reconcile-only --verify-h5ad-raw",
                            "block"))
        return {"status": "missing"}
    report = load_json(path)
    metrics = report.get("qc_field_comparison", {})
    good_values = report.get("all_raw_values_integer_nonnegative_finite") is True
    qc_match = all(metrics.get(k, {}).get("n_disagree") == 0
                   and metrics.get(k, {}).get("n_tested", 0) > 0
                   for k in ("nCount_RNA", "nFeature_RNA"))
    if not (good_values and qc_match):
        issues.append(issue("GSE225158", "raw_count_validation_failed",
                            "Count integrity or QC concordance not verified",
                            "Review raw/X and RNA QC metadata", "block"))
    return dict(integer_nonnegative_finite=good_values,
                count_and_feature_qc_match=qc_match,
                stored_values_scanned=report.get("raw_data_entries_checked"),
                upstream_count_provenance="not_independently_verified")


# %% MAIN
def main():
    reports = {acc: load_json(IN / f"{acc}_inspection.json")
               for acc in ACCESSIONS}
    OUT.mkdir(parents=True, exist_ok=True)
    issues = []
    bulk1, b1 = bulk_samples("GSE174409", reports["GSE174409"], issues)
    bulk2, b2 = bulk_samples("GSE194368", reports["GSE194368"], issues)
    stri, st = striatum_samples(reports["GSE225158"], issues)
    mids, mid = midbrain_libraries(reports["GSE240457"], issues)
    st["raw_verification"] = raw_count_verification(issues)
    issues.append(issue("GSE174409/GSE225158", "potential_shared_donors",
                        "12 candidate shared donors with concordant clinical metadata",
                        "Confirm cross-study identity using author/source metadata; "
                        "not independent biological replication"))
    issues.append(issue("GSE194368", "limited_covariate_metadata",
                        "GEO lacks key age/sex/PMI covariates",
                        "Search study supplementary tables and record unavailable covariates"))
    write_csv("bulk_samples.csv", bulk1 + bulk2,
              ["dataset", "matrix_column", "gsm", "geo_title", "donor_id",
               "donor_status", "region", "condition", "age", "sex", "pmi",
               "library_size", "genes_detected", "matched_to_geo"])
    write_csv("striatum_donor_region_samples.csv", stri,
              ["gsm", "sample_id", "subject_id", "region", "condition",
               "age", "sex", "pmi", "nuclei", "included_in_h5ad"])
    write_csv("midbrain_libraries.csv", mids,
              ["gsm", "library", "h5_filename", "genes", "filtered_barcodes",
               "donor_id", "diagnosis", "donor_mapping"])
    write_csv("unresolved_issues.csv", issues,
              ["dataset", "code", "priority", "observation", "next_step", "status"])
    readiness = {"GSE174409": b1, "GSE225158": st,
                 "GSE194368": b2, "GSE240457": mid,
                 "outstanding_issues": len(issues),
                 "interpretation": ("Sample design audit only. "
                                    "No cell filtering, DE, or integration performed."),
                 "midbrain_donor_level_de_ready": False}
    (OUT / "dataset_readiness.json").write_text(
        json.dumps(readiness, indent=2) + "\n", encoding="utf-8")
    print("Generated outputs/validation/human/:")
    for path in sorted(OUT.iterdir()):
        print(" ", path.name)
    for acc in ACCESSIONS:
        print(acc, readiness[acc])
    print("Open issues:", len(issues))


if __name__ == "__main__":
    main()
