# %% Phase 2b: read-only exploratory bulk QC
"""Unsupervised gene-count diagnostics. No filtering, exclusions, or DE.

Run: pixi run python scripts/02b_qc_bulk_multivariate.py
Requires Phase 2a outputs. Optional --annotation-csv has fields
dataset,feature_id,gene_symbol,chromosome (source-verified only).
"""
# %% Imports and configuration
from __future__ import annotations

# Shared standard-library I/O; support direct scripts and test imports.
try:
    from scripts._qc_common import read_csv, write_csv, sha256
except ModuleNotFoundError as exc:
    if exc.name != "scripts":
        raise
    from _qc_common import read_csv, write_csv, sha256
checksum = sha256

import argparse
import csv
import gzip
import hashlib
import heapq
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
QC = ROOT / "outputs" / "qc" / "human"
STUDIES = ("GSE174409", "GSE194368")
TOP_GENES = 2000





def load_annotation(path, dataset):
    """Only exact, externally verified feature mappings. No gene-ID guesses."""
    if path is None:
        return {}
    rows = read_csv(path)
    if not rows or not {"dataset", "feature_id", "gene_symbol", "chromosome"}.issubset(rows[0]):
        raise ValueError("Annotation CSV requires dataset,feature_id,gene_symbol,chromosome")
    relevant = [r for r in rows if r["dataset"] == dataset]
    mapping = {}
    for row in relevant:
        key = row["feature_id"].strip()
        if not key or key in mapping:
            raise ValueError(f"Empty or duplicate annotation ID: {key}")
        mapping[key] = row
    return mapping


# %% Stream one study at a time
def audit_cohort(path, samples, dataset, *, top_n=TOP_GENES, annotation=None,
                 exclude_pca_features=None):
    annotation = annotation or {}
    exclude_pca_features = set(exclude_pca_features or ())
    sample_ids = [r["matrix_sample"] for r in samples]
    if len(set(sample_ids)) != len(sample_ids):
        raise ValueError("Duplicate sample IDs in earlier QC")
    totals = np.array([int(r["assigned_gene_counts"]) for r in samples], dtype=float)
    if np.any(totals <= 0):
        raise ValueError("Zero gene-assigned count total; cannot compute logCPM")
    top = []
    mitochondrial_counts = np.zeros(len(samples), dtype=np.int64)
    found_mito = set()
    matched_annotations = 0
    counts_recomputed = np.zeros(len(samples), dtype=np.int64)
    seen = set()
    rows_read = 0
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t" if dataset == "GSE194368" else ",")
        header = next(reader)
        if header[1:] != sample_ids:
            raise ValueError(f"{dataset}: source header differs from Phase 2a QC")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if len(row) != len(sample_ids) + 1:
                raise ValueError(f"{dataset}: malformed row {rows_read + 1}")
            gene = row[0]
            if gene in seen:
                raise ValueError(f"{dataset}: duplicate feature {gene}")
            seen.add(gene)
            rows_read += 1
            values = np.array([float(v) for v in row[1:]], dtype=float)
            if not np.all(np.isfinite(values) & (values >= 0) & (values == np.floor(values))):
                raise ValueError(f"{dataset}: invalid values at feature {gene}")
            integer_values = values.astype(np.int64)
            counts_recomputed += integer_values
            logcpm = np.log1p(values / totals * 1e6)
            variance = float(logcpm.var())
            item = (variance, gene, logcpm)
            if gene not in exclude_pca_features:
                if len(top) < top_n:
                    heapq.heappush(top, item)
                elif (variance, gene) > (top[0][0], top[0][1]):
                    heapq.heapreplace(top, item)
            if gene in annotation:
                matched_annotations += 1
                chromosome = annotation[gene]["chromosome"].strip().upper().removeprefix("CHR")
                if chromosome in ("M", "MT"):
                    mitochondrial_counts += integer_values
                    found_mito.add(gene)
    if not np.array_equal(counts_recomputed, totals.astype(np.int64)):
        raise ValueError(f"{dataset}: totals mismatch Phase 2a; refusing diagnostics")
    if not top:
        raise ValueError("Empty feature matrix")
    top.sort(key=lambda r: (-r[0], r[1]))
    # Selected solely on unsupervised logCPM variance; this is not a DE filter.
    log_profiles = np.stack([item[2] for item in top], axis=1)
    x = log_profiles - log_profiles.mean(axis=0, keepdims=True)
    u, singular, _ = np.linalg.svd(x, full_matrices=False)
    scores = u * singular
    explained = (singular ** 2 / (singular ** 2).sum()).tolist() if (singular ** 2).sum() else [0.0] * len(singular)
    for axis in range(min(3, scores.shape[1])):
        pivot = int(np.argmax(np.abs(scores[:, axis])))
        if scores[pivot, axis] < 0:
            scores[:, axis] *= -1
    # Pearson correlation across selected genes: center EACH sample across genes.
    centered_profiles = log_profiles - log_profiles.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(centered_profiles, axis=1)
    denominator = np.outer(norms, norms)
    corr = np.divide(centered_profiles @ centered_profiles.T, denominator,
                     out=np.full_like(denominator, np.nan), where=denominator > 0)
    pca = []
    for i, sample in enumerate(samples):
        result = {key: sample.get(key, "") for key in
                  ("dataset", "matrix_sample", "gsm", "title", "diagnosis", "region", "age", "sex", "pmi", "rin")}
        result.update({f"PC{axis+1}": float(scores[i, axis]) if axis < scores.shape[1] else ""
                       for axis in range(3)})
        result["mt_count_fraction"] = (float(mitochondrial_counts[i] / totals[i]) if found_mito else "")
        pca.append(result)
    correlations = [
        {"dataset": dataset, "sample_a": sample_ids[i], "sample_b": sample_ids[j],
         "profile_correlation": float(corr[i, j])}
        for i in range(len(sample_ids)) for j in range(i + 1, len(sample_ids))
    ]
    features = [{"dataset": dataset, "feature_id": gene, "variance_log1p_cpm": variance,
                 "selected_for": "exploratory_PCA_and_correlation_only"}
                for variance, gene, _ in top]
    meta = {
        "dataset": dataset, "n_samples": len(samples), "n_features_read": rows_read,
        "pca_only_excluded_features": len(exclude_pca_features),
        "n_features_for_unsupervised_pca": len(top),
        "transformation": "log1p(gene_count / total_assigned_gene_counts * 1000000)",
        "pca_variance_explained_first_3": explained[:3],
        "annotation_source": "user_supplied_verified" if annotation else "not_provided",
        "annotation_matches": matched_annotations,
        "mitochondrial_mapped_features": len(found_mito),
        "mitochondrial_fraction_status": (
            "PROVISIONAL_VERIFY_ALL_37_MT_GENES_AND_ANNOTATION_BEFORE_INTERPRETATION"
            if found_mito else "NOT_CALCULATED_NO_VERIFIED_GENE_MAPPING"),
        "exclusions": 0, "differential_expression": False
    }
    return pca, correlations, features, meta


# %% Entry point
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotation-csv", type=Path, default=None,
                        help="Optional manually verified dataset,feature_id,gene_symbol,chromosome CSV")
    parser.add_argument("--top-genes", type=int, default=TOP_GENES,
                        help="Number of high-variance genes in exploratory PCA only")
    args = parser.parse_args()
    if args.top_genes < 2:
        parser.error("--top-genes must be >= 2")
    phase2a = json.loads((QC / "bulk_qc_readiness.json").read_text(encoding="utf-8"))
    by_dataset = {row["accession"]: row for row in phase2a["datasets"]}
    sample_qc = read_csv(QC / "bulk_sample_qc.csv")
    all_pca, all_corr, all_features, metadata = [], [], [], []
    for dataset in STUDIES:
        source = by_dataset[dataset]
        path = ROOT / source["source_file"]
        if checksum(path) != source["source_sha256"]:
            raise ValueError(f"{dataset}: source checksum differs from Phase 2a; STOP")
        samples = [row for row in sample_qc if row["dataset"] == dataset]
        if len(samples) != source["n_samples"]:
            raise ValueError(f"{dataset}: Phase 2a sample-count mismatch")
        annotation = load_annotation(args.annotation_csv, dataset)
        pca, corr, features, meta = audit_cohort(path, samples, dataset,
                                                  top_n=args.top_genes, annotation=annotation)
        if meta["n_features_read"] != source["n_features"]:
            raise ValueError(f"{dataset}: Phase 2a feature-count mismatch")
        meta["source_sha256"] = source["source_sha256"]
        all_pca.extend(pca)
        all_corr.extend(corr)
        all_features.extend(features)
        metadata.append(meta)
        print(f"{dataset}: {len(pca)} samples; {len(features)} exploratory PCA genes; no exclusions")
    write_csv(QC / "bulk_pca_scores.csv", all_pca,
              ["dataset", "matrix_sample", "gsm", "title", "diagnosis", "region", "age", "sex",
               "pmi", "rin", "PC1", "PC2", "PC3", "mt_count_fraction"])
    write_csv(QC / "bulk_sample_correlations.csv", all_corr,
              ["dataset", "sample_a", "sample_b", "profile_correlation"])
    write_csv(QC / "bulk_pca_feature_selection.csv", all_features,
              ["dataset", "feature_id", "variance_log1p_cpm", "selected_for"])
    summary = {
        "phase": "2b_read_only_unsupervised_diagnostics", "cohorts": metadata,
        "filtering_performed": False, "differential_expression_performed": False,
        "source_count_matrices_modified": False,
        "limitations": [
            "PCA/correlations are unsupervised exploratory diagnostics, not inferential evidence",
            "Gene selection for visualization is not a downstream DE gene filter",
            "Mitochondrial fractions require audited Ensembl/Entrez feature-ID mapping",
            "PC axes may reflect biology or technical covariates and need metadata review"
        ]
    }
    (QC / "bulk_multivariate_readiness.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print("Finished read-only Phase 2b; no sample or feature exclusions.")


if __name__ == "__main__":
    main()
