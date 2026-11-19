# %% Phase 2e: read-only 10x pooled midbrain QC
"""Describe GSE240457 Cell Ranger-filtered 10x libraries, not individual donors.

Run: pixi run python scripts/02e_qc_midbrain.py

Requires outputs/validation/human/midbrain_libraries.csv and
outputs/inspection/GSE240457_inspection.json. Nothing is removed or assigned
to a donor. One tar member is materialized at a time to a temporary directory.
"""
# %% Imports and setup
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import shutil
import tarfile
import tempfile

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/raw/human/GSE240457/GSE240457_RAW.tar"
INSPECTION = ROOT / "outputs/inspection/GSE240457_inspection.json"
VALIDATED = ROOT / "outputs/validation/human/midbrain_libraries.csv"
OUT = ROOT / "outputs/qc/human"
CHUNK_COLUMNS = 2048
COPY_CHUNK = 8 * 1024 * 1024
MIN_RESERVE = 2_000_000_000
MAX_LIBRARY_BYTES = 20_000_000_000


def as_text(values):
    return [x.decode("utf-8", "replace") if isinstance(x, bytes) else str(x) for x in values]


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(COPY_CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def quantiles(values):
    if len(values) == 0:
        return {"min": None, "p05": None, "median": None,
                "p95": None, "max": None}
    return {label: round(float(np.percentile(values, percentile)), 5)
            for label, percentile in (
                ("min", 0), ("p05", 5), ("median", 50),
                ("p95", 95), ("max", 100))}


# %% Safe CSC column-block summary; no dense cell × feature matrix
def summarize_10x_h5(path, *, block=CHUNK_COLUMNS):
    if block < 1:
        raise ValueError("Block size must be positive")
    with h5py.File(path, "r") as file:
        base = file["matrix"]
        required = {"data", "indices", "indptr", "shape", "barcodes", "features"}
        if not required.issubset(base.keys()):
            raise ValueError(f"Invalid 10x H5: missing {sorted(required - set(base))}")
        n_features, n_barcodes = (int(x) for x in base["shape"][:])
        features = base["features"]
        for field in ("id", "name", "feature_type"):
            if field not in features:
                raise ValueError(f"Missing 10x feature field {field}")
        if len(features["id"]) != n_features or len(features["name"]) != n_features:
            raise ValueError("10x feature metadata does not match matrix shape")
        if len(features["feature_type"]) != n_features:
            raise ValueError("10x feature types do not match matrix shape")
        if len(base["barcodes"]) != n_barcodes or len(base["indptr"]) != n_barcodes + 1:
            raise ValueError("10x barcode/indptr dimensions inconsistent")
        names = as_text(features["name"][:])
        ids = as_text(features["id"][:])
        feature_types = as_text(features["feature_type"][:])
        kinds = Counter(feature_types)
        gene_mask = np.array([x == "Gene Expression" for x in feature_types], dtype=bool)
        if not gene_mask.any():
            raise ValueError("10x matrix lacks Gene Expression features")
        if len(set(ids[i] for i in np.flatnonzero(gene_mask))) != gene_mask.sum():
            raise ValueError("Duplicate Gene Expression feature IDs; QC cannot interpret")
        # MT- symbols are a provisional naming diagnostic, NOT an annotation-verified MT fraction.
        mt_mask = np.array([bool(gene_mask[i] and name.startswith("MT-"))
                            for i, name in enumerate(names)], dtype=bool)
        genome = as_text(features["genome"][:]) if "genome" in features else []
        totals = np.zeros(n_barcodes, dtype=np.int64)
        detected = np.zeros(n_barcodes, dtype=np.int64)
        mt_counts = np.zeros(n_barcodes, dtype=np.int64)
        ptr = base["indptr"]
        values = base["data"]
        indices = base["indices"]
        n_stored = len(values)
        if len(indices) != n_stored:
            raise ValueError("10x indices/data lengths do not match")
        if int(ptr[0]) != 0 or int(ptr[-1]) != n_stored:
            raise ValueError("10x indptr boundaries invalid")
        for a in range(0, n_barcodes, block):
            b = min(a + block, n_barcodes)
            offsets = ptr[a:b+1].astype(np.int64)
            if np.any(np.diff(offsets) < 0):
                raise ValueError("10x indptr not sorted")
            start, stop = int(offsets[0]), int(offsets[-1])
            cols = np.repeat(np.arange(b-a), np.diff(offsets))
            if start == stop:
                continue
            idx = indices[start:stop]
            data = values[start:stop]
            if np.any(idx < 0) or np.any(idx >= n_features):
                raise ValueError("10x feature index outside shape")
            if not np.all(np.isfinite(data)) or np.any(data < 0) or np.any(data != np.floor(data)):
                raise ValueError("10x matrix contains invalid UMI counts")
            active = gene_mask[idx]
            gene_cols = cols[active]
            gene_data = data[active].astype(np.int64)
            totals[a:b] = np.bincount(gene_cols, weights=gene_data,
                                      minlength=b-a).astype(np.int64)
            detected[a:b] = np.bincount(gene_cols,
                                        weights=(gene_data > 0).astype(np.int64),
                                        minlength=b-a).astype(np.int64)
            mt_active = active & mt_mask[idx]
            mt_counts[a:b] = np.bincount(cols[mt_active],
                                         weights=data[mt_active],
                                         minlength=b-a).astype(np.int64)
        if (totals < 0).any() or (detected > int(gene_mask.sum())).any():
            raise ValueError("Invalid 10x barcode metric aggregation")
        fractions = np.divide(mt_counts, totals, out=np.full(n_barcodes, np.nan),
                              where=totals > 0)
        return {
            "matrix_features_all_types": n_features,
            "gene_expression_features": int(gene_mask.sum()),
            "filtered_barcodes": n_barcodes,
            "feature_type_counts": dict(kinds),
            "genomes": sorted(set(genome)),
            "stored_entries": n_stored,
            "gene_count_sum": int(totals.sum()),
            "zero_gene_umi_barcodes": int((totals == 0).sum()),
            "umi_quantiles": quantiles(totals),
            "detected_gene_quantiles": quantiles(detected),
            "putative_mt_gene_features": int(mt_mask.sum()),
            "mt_fraction_quantiles_provisional": quantiles(fractions[np.isfinite(fractions)]),
            "mt_fraction_status": "NAME_PREFIX_ONLY_NOT_SOURCE_VERIFIED",
            "donor_assignment": "NOT_AVAILABLE",
            "clinical_diagnosis": "NOT_AVAILABLE",
        }


# %% Archive inventory and validation
def audit_archive(archive, validated, inspection, *, tmp_root, max_libraries=None):
    expected = {row["h5_filename"]: row for row in validated}
    if len(expected) != len(validated):
        raise ValueError("Duplicate validated 10x filenames")
    inspected = {x["filename"]: x for x in inspection["expression"]["libraries"]}
    if set(expected) != set(inspected):
        raise ValueError("Validation inventory and inspection report disagree")
    summaries = []
    with tarfile.open(archive, "r") as tar:
        members = [m for m in tar.getmembers() if m.isfile() and m.name.endswith(".h5")]
        found = {m.name: m for m in members}
        if len(found) != len(members) or set(found) != set(expected):
            raise ValueError("Archive H5 members do not match validated library inventory")
        names = sorted(found)
        if max_libraries is not None:
            names = names[:max_libraries]
        for i, name in enumerate(names, 1):
            member = found[name]
            if not (0 < member.size <= MAX_LIBRARY_BYTES):
                raise ValueError(f"Library member has invalid size: {name}")
            free = shutil.disk_usage(tmp_root).free
            if free < member.size + MIN_RESERVE:
                raise OSError(f"Insufficient free disk space for {name}")
            # Always write to a fixed name; never trust a tar member path on disk.
            with tempfile.TemporaryDirectory(prefix="midbrain-10x-", dir=tmp_root) as temp:
                target = Path(temp) / "library.h5"
                stream = tar.extractfile(member)
                if stream is None:
                    raise ValueError(f"Cannot read tar member {name}")
                with stream, target.open("wb") as dest:
                    shutil.copyfileobj(stream, dest, length=COPY_CHUNK)
                if target.stat().st_size != member.size:
                    raise ValueError(f"Extracted byte count mismatch: {name}")
                stats = summarize_10x_h5(target)
            prior = inspected[name]
            if stats["filtered_barcodes"] != int(prior["barcode_count"]):
                raise ValueError(f"10x barcode count changed for {name}")
            if stats["matrix_features_all_types"] != int(prior["matrix_shape"][0]):
                raise ValueError(f"10x feature count changed for {name}")
            row = expected[name]
            if stats["filtered_barcodes"] != int(row["filtered_barcodes"]):
                raise ValueError(f"Validated barcode mismatch: {name}")
            summaries.append({"gsm": row["gsm"], "library": row["library"],
                              "filename": name, **stats})
            print(f"[{i}/{len(names)}] {row['gsm']}: "
                  f"{stats['filtered_barcodes']} filtered barcodes; no exclusions")
    return summaries


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive", type=Path, default=SOURCE)
    p.add_argument("--temp-dir", type=Path, default=None,
                   help="Directory for extraction of one H5 at a time")
    p.add_argument("--max-libraries", type=int, default=None,
                   help="Diagnostic subset only; partial output explicitly flagged")
    args = p.parse_args()
    if args.max_libraries is not None and args.max_libraries < 1:
        p.error("--max-libraries must be >= 1")
    if not args.archive.is_file():
        p.error(f"Missing tar archive: {args.archive}")
    tmp_root = args.temp_dir or args.archive.parent
    if not tmp_root.is_dir():
        p.error(f"Temporary directory does not exist: {tmp_root}")
    validation = read_csv(VALIDATED)
    inspection = json.loads(INSPECTION.read_text(encoding="utf-8"))
    result = audit_archive(args.archive, validation, inspection, tmp_root=tmp_root,
                           max_libraries=args.max_libraries)
    full = len(result) == len(validation)
    if full and sum(r["filtered_barcodes"] for r in result) != 269307:
        raise ValueError("Filtered barcodes inconsistent with earlier GSE240457 audit")
    OUT.mkdir(parents=True, exist_ok=True)
    write_csv(OUT / "midbrain_library_qc.csv", [{
        **row, "umi_median": row["umi_quantiles"]["median"],
        "umi_p05": row["umi_quantiles"]["p05"], "umi_p95": row["umi_quantiles"]["p95"],
        "genes_median": row["detected_gene_quantiles"]["median"],
        "genes_p05": row["detected_gene_quantiles"]["p05"],
        "genes_p95": row["detected_gene_quantiles"]["p95"],
        "mt_fraction_median_provisional": row["mt_fraction_quantiles_provisional"]["median"],
    } for row in result], [
        "gsm", "library", "filename", "filtered_barcodes", "gene_expression_features",
        "matrix_features_all_types", "stored_entries", "gene_count_sum",
        "zero_gene_umi_barcodes", "umi_median", "umi_p05", "umi_p95",
        "genes_median", "genes_p05", "genes_p95",
        "putative_mt_gene_features", "mt_fraction_median_provisional",
        "mt_fraction_status", "donor_assignment", "clinical_diagnosis"])
    summary = {
        "phase": "2e_read_only_midbrain_library_qc",
        "dataset": "GSE240457",
        "source_archive": str(args.archive),
        "archive_sha256": sha256(args.archive),
        "libraries_analyzed": len(result),
        "libraries_expected": len(validation),
        "complete": full,
        "all_filtered_barcodes": sum(r["filtered_barcodes"] for r in result),
        "libraries": result,
        "donor_level_de_ready": False,
        "sample_exclusions": 0,
        "barcode_exclusions": 0,
        "limitations": [
            "These are Cell Ranger-filtered barcodes, not authors' final annotated nuclei",
            "Library identifiers are multiplexed pools, never donor IDs",
            "MT- gene-name prefix gives provisional descriptive QC only",
            "No demuxlet mapping, doublet removal, diagnosis, or cell types inferred",
            "No sample or barcode exclusion, normalization, integration or DE"
        ],
    }
    (OUT / "midbrain_qc_readiness.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print("Midbrain technical QC complete. Donor-level inference remains blocked.")


if __name__ == "__main__":
    main()
