# %% Phase 2c: reference-verified annotation and mitochondrial audit
"""Read-only bulk gene identifier audit, MT fractions, and PCA sensitivity.

References must be explicit local files; no inferred mapping or auto-download:
  pixi run python scripts/02c_qc_bulk_annotation.py \
      --gencode-gtf data/reference/gencode.annotation.gtf.gz \
      --ncbi-gene-info data/reference/Homo_sapiens.gene_info.gz

The GENCODE release should be chosen/documented against original GSE174409
alignment annotation; NCBI gene_info is an ID/MT lookup, NOT proof of the
original GSE194368 counting reference. No samples/genes are removed for DE.
"""
# %% Imports
from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
QC = ROOT / "outputs" / "qc" / "human"
STUDIES = ("GSE174409", "GSE194368")
GTF_ID = re.compile(r'(?:^|;\s*)gene_id\s+"([^"]+)"')
GTF_NAME = re.compile(r'(?:^|;\s*)gene_name\s+"([^"]+)"')


def read_text(path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open("rt", encoding="utf-8")


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for buf in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(buf)
    return h.hexdigest()


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


# %% Annotation parsing: preserve ambiguity and assembly/source provenance
def gtf_records(path):
    """Parse only gene lines; do not infer chromosome from gene symbol."""
    records = {}
    collisions = set()
    with read_text(path) as source:
        for line in source:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            match_id = GTF_ID.search(fields[8])
            match_name = GTF_NAME.search(fields[8])
            if not match_id:
                continue
            raw_id = match_id.group(1)
            key = raw_id.split(".")[0]  # versionless match; collision audited
            record = {
                "gene_symbol": match_name.group(1) if match_name else "",
                "chromosome": fields[0], "reference_id": raw_id,
            }
            if key in records and records[key] != record:
                collisions.add(key)
            else:
                records[key] = record
    for key in collisions:
        records.pop(key, None)
    if not records:
        raise ValueError("No GTF gene records parsed; verify the annotation file")
    return records, collisions


def ncbi_records(path):
    """Parse NCBI gene_info, human only; exact Entrez IDs."""
    records = {}
    collisions = set()
    with read_text(path) as source:
        reader = csv.DictReader(source, delimiter="\t")
        needed = {"#tax_id", "GeneID", "Symbol", "chromosome"}
        if not reader.fieldnames or not needed.issubset(reader.fieldnames):
            raise ValueError("Expected NCBI gene_info columns: " + ",".join(sorted(needed)))
        for row in reader:
            if row["#tax_id"] != "9606":
                continue
            key = row["GeneID"]
            record = {"gene_symbol": row["Symbol"], "chromosome": row["chromosome"],
                      "reference_id": key}
            if key in records and records[key] != record:
                collisions.add(key)
            else:
                records[key] = record
    for key in collisions:
        records.pop(key, None)
    if not records:
        raise ValueError("No human gene_info records parsed")
    return records, collisions


def canonical_mito(chromosome):
    """Chromosome is source annotated, never symbol-prefix derived."""
    chrom = chromosome.strip().upper()
    if chrom.startswith("CHR"):
        chrom = chrom[3:]
    return chrom in ("M", "MT")


def map_features(matrix_ids, reference, ambiguous, dataset):
    """Return conservative feature mappings and audit states."""
    mapping = {}
    audit = []
    seen = set()
    for raw in matrix_ids:
        if raw in seen:
            raise ValueError(f"{dataset}: duplicate matrix feature {raw}")
        seen.add(raw)
        key = raw.split(".")[0] if dataset == "GSE174409" else raw
        if key in ambiguous:
            status, record = "AMBIGUOUS_REFERENCE_ID", None
        elif key not in reference:
            status, record = "UNMAPPED", None
        else:
            status, record = "MAPPED", reference[key]
            mapping[raw] = record
        audit.append({
            "dataset": dataset, "feature_id": raw, "reference_key": key,
            "gene_symbol": record["gene_symbol"] if record else "",
            "chromosome": record["chromosome"] if record else "",
            "mapping_status": status,
            "mitochondrial_by_reference": (
                "yes" if record and canonical_mito(record["chromosome"]) else
                "no" if record else "unknown"),
        })
    return mapping, audit


def feature_ids(path, dataset):
    """One pass: obtain exact matrix gene-ID universe, preserving source order."""
    values = []
    seen = set()
    with read_text(path) as inp:
        reader = csv.reader(inp, delimiter="\t" if dataset == "GSE194368" else ",")
        next(reader)
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            key = row[0]
            if key in seen:
                raise ValueError(f"{dataset}: duplicate matrix feature {key}")
            seen.add(key)
            values.append(key)
    return values


# %% Phase 2b reuse for a matched no-MT PCA sensitivity, no count changes
def load_multivariate():
    path = ROOT / "scripts" / "02b_qc_bulk_multivariate.py"
    spec = importlib.util.spec_from_file_location("bulk_qc_02b_for_sensitivity", path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gencode-gtf", type=Path, required=True)
    p.add_argument("--ncbi-gene-info", type=Path, required=True)
    args = p.parse_args()
    for path in (args.gencode_gtf, args.ncbi_gene_info):
        if not path.is_file():
            p.error(f"Reference missing: {path}. Download and record reference release first.")
    source_qc = json.loads((QC / "bulk_qc_readiness.json").read_text(encoding="utf-8"))
    prev = {r["accession"]: r for r in source_qc["datasets"]}
    sample_rows = read_csv(QC / "bulk_sample_qc.csv")
    gencode, gencode_ambiguity = gtf_records(args.gencode_gtf)
    ncbi, ncbi_ambiguity = ncbi_records(args.ncbi_gene_info)
    multiv = load_multivariate()

    annotations = []
    fractions = []
    sensitivity = []
    summaries = []
    for dataset in STUDIES:
        source = prev[dataset]
        counts_path = ROOT / source["source_file"]
        if sha256(counts_path) != source["source_sha256"]:
            raise ValueError(f"{dataset}: count file checksum mismatch; stop")
        samples = [r for r in sample_rows if r["dataset"] == dataset]
        reference, ambiguous, ref_path = (
            (gencode, gencode_ambiguity, args.gencode_gtf)
            if dataset == "GSE174409" else
            (ncbi, ncbi_ambiguity, args.ncbi_gene_info)
        )
        ids = feature_ids(counts_path, dataset)
        if len(ids) != source["n_features"]:
            raise ValueError(f"{dataset}: feature count mismatch")
        mapping, audit = map_features(ids, reference, ambiguous, dataset)
        if len(mapping) / len(ids) < 0.5:
            raise ValueError(
                f"{dataset}: fewer than 50% of count features map to the supplied reference; "
                "check assembly, ID namespace and release before proceeding")
        annotations.extend(audit)
        mt_ids = {key for key, value in mapping.items() if canonical_mito(value["chromosome"])}
        # No exclusions: excluded features affect only this optional PCA sensitivity.
        pca_all, _, all_features, all_meta = multiv.audit_cohort(
            counts_path, samples, dataset, annotation=mapping)
        pca_no_mt, _, no_mt_features, sensitivity_meta = multiv.audit_cohort(
            counts_path, samples, dataset, annotation=mapping,
            exclude_pca_features=mt_ids)
        for a, b in zip(pca_all, pca_no_mt):
            if a["matrix_sample"] != b["matrix_sample"]:
                raise ValueError("PCA sample order unexpectedly changed")
            fractions.append({
                "dataset": dataset, "matrix_sample": a["matrix_sample"],
                "mt_count_fraction": a["mt_count_fraction"],
                "mt_mapping_status": "REVIEW_FEATURE_COMPLETENESS",
                "mt_features_found_in_matrix": len(mt_ids),
            })
            sensitivity.append({
                "dataset": dataset, "matrix_sample": a["matrix_sample"],
                "PC1_all": a["PC1"], "PC2_all": a["PC2"],
                "PC1_without_MT": b["PC1"], "PC2_without_MT": b["PC2"],
                "PCs_not_directly_comparable": True,
            })
        # Report set overlap rather than comparing the directions of PCs.
        orig = {r["feature_id"] for r in all_features}
        without = {r["feature_id"] for r in no_mt_features}
        symbols = Counter(row["gene_symbol"] for row in audit
                          if row["mapping_status"] == "MAPPED" and row["gene_symbol"])
        summary = {
            "dataset": dataset, "reference_path": str(ref_path),
            "reference_sha256": sha256(ref_path), "count_source_sha256": source["source_sha256"],
            "total_features": len(ids),
            "mapped_features": len(mapping),
            "unmapped_features": len(ids) - len(mapping),
            "ambiguous_reference_ids": len(ambiguous),
            "duplicated_gene_symbols_in_mapped_features": sum(v > 1 for v in symbols.values()),
            "reference_mt_gene_records": sum(canonical_mito(r["chromosome"]) for r in reference.values()),
            "mt_genes_present_in_count_matrix": len(mt_ids),
            "mt_fraction_interpretation": "PROVISIONAL_PENDING_REFERENCE_AND_FEATURE_COVERAGE_REVIEW",
            "pca_gene_set_overlap_all_vs_no_mt": len(orig & without),
            "pca_all_pc1_variance_fraction": all_meta["pca_variance_explained_first_3"][0],
            "pca_no_mt_pc1_variance_fraction": sensitivity_meta["pca_variance_explained_first_3"][0],
            "pca_only_removed_mt_features": sensitivity_meta["pca_only_excluded_features"],
            "counts_modified": False, "DE_feature_exclusions": 0,
        }
        summaries.append(summary)
        print(f"{dataset}: {len(mapping)}/{len(ids)} IDs matched; "
              f"{len(mt_ids)} mitochondrial features in matrix; no sample exclusions")
    write_csv(QC / "bulk_feature_annotation_audit.csv", annotations,
              ["dataset", "feature_id", "reference_key", "gene_symbol", "chromosome",
               "mapping_status", "mitochondrial_by_reference"])
    write_csv(QC / "bulk_mitochondrial_qc.csv", fractions,
              ["dataset", "matrix_sample", "mt_count_fraction",
               "mt_mapping_status", "mt_features_found_in_matrix"])
    write_csv(QC / "bulk_pca_mt_sensitivity.csv", sensitivity,
              ["dataset", "matrix_sample", "PC1_all", "PC2_all", "PC1_without_MT",
               "PC2_without_MT", "PCs_not_directly_comparable"])
    (QC / "bulk_annotation_readiness.json").write_text(json.dumps({
        "phase": "2c_annotation_and_mt_sensitivity",
        "cohorts": summaries, "exclusions": 0,
        "status": "EXPLORATORY_REVIEW_ONLY",
        "limits": [
            "GENCODE / NCBI reference snapshots may differ from original count generation.",
            "MT count fractions depend on annotation and feature coverage; not a standalone degradation marker.",
            "PCA signs/axes can rotate; never subtract raw PC scores across analyses.",
            "Mitochondrial features are excluded from sensitivity PCA ONLY, not from data or DE.",
            "Final gene identifier harmonization requires reference-release review."
        ]
    }, indent=2) + "\n", encoding="utf-8")
    print("Phase 2c finished. Inputs preserved; no filtering or DE.")


if __name__ == "__main__":
    main()
