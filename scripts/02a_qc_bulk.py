# %% Phase 2: read-only human bulk count-matrix QC
"""Descriptive QC of deposited gene counts; no sample or gene filtering.

Run:
    pixi run python scripts/02a_qc_bulk.py

Requires existing local raw files and inspection JSONs. No extra dependencies.
Outputs are written under ignored outputs/qc/human/.
"""

# %% Configuration
from __future__ import annotations

from collections import Counter
from pathlib import Path
import csv
import gzip
import hashlib
import json
import math
import statistics

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "human"
INSPECT = ROOT / "outputs" / "inspection"
OUT = ROOT / "outputs" / "qc" / "human"
COHORTS = ("GSE174409", "GSE194368")
MAX_CONTRIBUTORS = 10
MIN_STRATUM_FOR_FLAGS = 5
FLAG_MAD = 3.5  # descriptive flag only, not an exclusion threshold


# %% Helpers
def read_json(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"Missing {path}; run inspection first")
    return json.loads(path.read_text(encoding="utf-8"))


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(4 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def median_mad(values: list[float]) -> tuple[float, float]:
    middle = statistics.median(values)
    deviation = statistics.median([abs(x - middle) for x in values])
    return middle, deviation


def describe_flags(rows: list[dict], attribute: str, grouped: bool) -> None:
    strata = {}
    for row in rows:
        key = row["region"] if grouped else "all"
        strata.setdefault(key, []).append(row)

    for key, subset in strata.items():
        values = [math.log1p(float(row[attribute])) for row in subset]
        # Three-sigma-like robust z-score is an aid to *review*, never exclusion.
        if len(values) < MIN_STRATUM_FOR_FLAGS:
            for row in subset:
                row[f"{attribute}_robust_z"] = ""
                row[f"{attribute}_flag"] = "INSUFFICIENT_STRATUM"
            continue
        med, mad = median_mad(values)
        for row, val in zip(subset, values):
            if mad == 0:
                row[f"{attribute}_robust_z"] = ""
                row[f"{attribute}_flag"] = "UNDEFINED_MAD"
                continue
            z = 0.67448975 * (val - med) / mad
            row[f"{attribute}_robust_z"] = round(z, 3)
            row[f"{attribute}_flag"] = "REVIEW" if abs(z) >= FLAG_MAD else ""


# %% Single-pass row streaming
def audit_matrix(path: Path, accession: str, inspection: dict) -> tuple[list[dict], dict]:
    matrix = inspection["expression"][0]
    expected = matrix["sample_names"]
    is_tab = accession == "GSE194368"
    opener = gzip.open if path.suffix == ".gz" else open
    geo = inspection["geo"]["samples"]
    meta_key = "geo_accession" if accession == "GSE194368" else "title"
    by_meta = {item[meta_key]: item for item in geo}

    if len(set(expected)) != len(expected):
        raise ValueError("Duplicate sample columns in inspection JSON")
    if set(expected) - set(by_meta):
        raise ValueError("Unmapped count matrix columns in GEO metadata")

    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t" if is_tab else ",")
        header = next(reader)
        observed = header[1:]
        if observed != expected:
            raise ValueError(f"{accession}: raw header differs from inspection, refusing QC")
        n = len(expected)
        counts = [0] * n
        detected = [0] * n
        # Track the top ten count contributions in each sample.
        top = [[] for _ in range(n)]
        genes = 0
        duplicates = set()
        seen = set()
        for line in reader:
            if not line or line[0].startswith("#"):
                continue
            if len(line) != n + 1:
                raise ValueError(f"{accession}: malformed gene row {genes + 1}")
            genes += 1
            gene = line[0]
            if gene in seen:
                duplicates.add(gene)
            seen.add(gene)
            for i, token in enumerate(line[1:]):
                if not token:
                    raise ValueError(f"{accession}: empty count at gene {gene}")
                # Accept integer-encoded floating strings only if exact and nonnegative.
                value = float(token)
                if not math.isfinite(value) or value < 0 or not value.is_integer():
                    raise ValueError(f"{accession}: invalid count at gene {gene}")
                v = int(value)
                counts[i] += v
                if v:
                    detected[i] += 1
                    vals = top[i]
                    if len(vals) < MAX_CONTRIBUTORS:
                        vals.append(v)
                        vals.sort()
                    elif v > vals[0]:
                        vals[0] = v
                        vals.sort()

    if genes != matrix["genes"]:
        raise ValueError(f"{accession}: gene-row count changed")
    if duplicates:
        raise ValueError(f"{accession}: duplicated feature IDs, examples: {sorted(duplicates)[:3]}")
    previous = {row["sample"]: row for row in matrix["per_sample"]}
    results = []
    for i, sample in enumerate(expected):
        prior = previous[sample]
        if counts[i] != prior["library_size"] or detected[i] != prior["genes_detected"]:
            raise ValueError(f"{accession}: independent QC disagrees with inspection for {sample}")
        info = by_meta[sample]
        region = info.get("source_name_ch1", "")
        results.append({
            "dataset": accession, "matrix_sample": sample,
            "gsm": info.get("geo_accession", ""), "title": info.get("title", ""),
            "region": region, "diagnosis": info.get("dx_oud", info.get("disease_state", info.get("diagnosis", ""))),
            "age": info.get("age", ""), "sex": info.get("sex", ""),
            "pmi": info.get("pmi", ""), "rin": info.get("rin", ""),
            "assigned_gene_counts": counts[i],
            "detected_features": detected[i],
            "zero_feature_fraction": round(1 - detected[i] / genes, 6),
            "top10_count_fraction": round(sum(top[i]) / counts[i], 6) if counts[i] else None,
            "technical_exclusion": "NONE_APPLIED"
        })
    for field in ("assigned_gene_counts", "detected_features"):
        describe_flags(results, field, grouped=(accession == "GSE174409"))
    metadata = {
        "accession": accession,
        "data_kind": "submitter_count_matrix" if accession == "GSE174409" else "NCBI_reprocessed_count_matrix",
        "source_file": str(path.relative_to(ROOT)),
        "source_sha256": sha256(path),
        "n_samples": n, "n_features": genes,
        "sample_id_kind": "GEO_title" if accession == "GSE174409" else "GSM",
        "independent_donor_id_verified": False if accession == "GSE174409" else "not_documented_in_manifest",
        "metrics": ["assigned_gene_counts", "detected_features", "zero_feature_fraction",
                    "top10_count_fraction", "within_stratum_log1p_robust_z"],
        "flags_are_exclusions": False,
        "read_level_qc_available": False,
        "count_total_is": "gene_assigned_counts_not_total_sequenced_reads",
        "limitations": [
            "This stage does not measure mapping %, base quality, PCR duplication, or gene body bias.",
            "Gene identifiers have not been independently mapped to symbols or complement targets.",
            "OUD-related expression was not tested or used for any QC threshold."
        ]
    }
    return results, metadata


# %% Entry point
def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    all_rows = []
    summaries = []
    for accession in COHORTS:
        report = read_json(INSPECT / f"{accession}_inspection.json")
        name = report["expression"][0]["file"]
        candidate = RAW / accession / name
        if not candidate.is_file():
            raise FileNotFoundError(f"Expected original matrix missing: {candidate}")
        rows, summary = audit_matrix(candidate, accession, report)
        all_rows.extend(rows)
        summaries.append(summary)
        print(f"{accession}: {summary['n_samples']} samples, "
              f"{summary['n_features']} gene rows; checks match inspection")

    write_csv(OUT / "bulk_sample_qc.csv", all_rows, [
        "dataset", "matrix_sample", "gsm", "title", "region", "diagnosis", "age",
        "sex", "pmi", "rin", "assigned_gene_counts", "detected_features",
        "zero_feature_fraction", "top10_count_fraction",
        "assigned_gene_counts_robust_z", "assigned_gene_counts_flag",
        "detected_features_robust_z", "detected_features_flag",
        "technical_exclusion"])
    summary = {
        "phase": "2_read_only_descriptive_bulk_QC",
        "datasets": summaries, "n_rows": len(all_rows),
        "filtering_performed": False,
        "normalization_performed": False,
        "differential_expression_performed": False,
        "next_gate": "Review flagged samples, metadata confounding, count provenance "
                     "and source-author exclusions before any filtering or inference."
    }
    (OUT / "bulk_qc_readiness.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUT.relative_to(ROOT)}")
    print("Read-only QC complete. No samples removed.")


if __name__ == "__main__":
    main()
