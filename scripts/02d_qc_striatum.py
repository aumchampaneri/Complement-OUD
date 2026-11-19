# %% Phase 2d — GSE225158 author-retained single-nucleus QC
"""Read-only audit of inherited QC and donor × region × cell type coverage.

Run: pixi run python scripts/02d_qc_striatum.py
Requires GSE225158 inspection JSON, 01b validation table, and source H5AD.gz.
No additional QC filtering, count re-normalization, or DE is performed.
"""
# %% Imports and configuration
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import shutil
import tempfile

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE_DIR = ROOT / "data" / "raw" / "human" / "GSE225158"
INSPECTION = ROOT / "outputs" / "inspection" / "GSE225158_inspection.json"
GEO_TABLE = ROOT / "outputs" / "validation" / "human" / "striatum_donor_region_samples.csv"
OUT = ROOT / "outputs" / "qc" / "human"
REQUIRED = ("Case", "ID", "Dx_OUD", "Region", "celltype1", "nCount_RNA", "nFeature_RNA")
OPTIONAL = ("celltype2", "celltype3", "orig.ident", "percent.mt",
            "dropletQC.keep", "miQC.keep", "scds.keep",
            "dropletQC", "miQC", "scds")
MIN_FREE_BYTES = 13_000_000_000
CHUNK_BYTES = 8 * 1024 * 1024


# %% Safe I/O
def text_value(value):
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", "replace")
    if isinstance(value, np.generic):
        return value.item()
    return value


def obs_vector(obs, key, size):
    """Read AnnData's categorical or 1-D HDF5 obs encoding, never expression."""
    node = obs[key]
    if isinstance(node, h5py.Group) and {"categories", "codes"} <= set(node.keys()):
        cats = node["categories"][:]
        codes = node["codes"][:]
        if len(codes) != size:
            raise ValueError(f"{key}: categorical length disagrees with nuclei")
        return np.array([text_value(cats[int(i)]) if i >= 0 else None for i in codes],
                        dtype=object)
    if isinstance(node, h5py.Dataset) and node.ndim == 1 and len(node) == size:
        values = node[:]
        if values.dtype.kind in ("S", "O", "U"):
            return np.array([text_value(v) for v in values], dtype=object)
        return values
    raise ValueError(f"{key}: unsupported AnnData obs representation")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as inp:
        for chunk in iter(lambda: inp.read(CHUNK_BYTES), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_csv(path):
    with path.open(newline="", encoding="utf-8") as file:
        return list(csv.DictReader(file))


def fmt_bool_flag(value):
    """Treat unavailable/ambiguous flags as unknown; do not guess their meaning."""
    if value is None:
        return "UNKNOWN"
    if isinstance(value, (bool, np.bool_)):
        return "TRUE" if value else "FALSE"
    token = str(text_value(value)).strip().lower()
    if token in {"true", "t", "1"}:
        return "TRUE"
    if token in {"false", "f", "0"}:
        return "FALSE"
    return "UNKNOWN"


def assert_identifiers(fields, geo):
    """Never infer one donor per library; verify donor–region IDs against GEO."""
    known = {r["sample_id"]: r for r in geo}
    if len(known) != len(geo):
        raise ValueError("Duplicate GEO sample IDs")
    by_donor = defaultdict(set)
    sample_to_donor = {}
    for sample, donor, region, condition in zip(
            fields["ID"], fields["Case"], fields["Region"], fields["Dx_OUD"]):
        sid, did = str(sample), str(donor)
        if sid not in known:
            raise ValueError(f"Object contains unknown donor–region sample: {sid}")
        expected = known[sid]
        expected_donor = str(expected["subject_id"])
        if did != expected_donor and did != sid.split("-", 1)[-1]:
            raise ValueError(f"{sid}: donor ID {did} differs from GEO {expected_donor}")
        previous = sample_to_donor.setdefault(sid, did)
        if previous != did:
            raise ValueError(f"{sid}: multiple donor identities")
        by_donor[did].add(str(condition))
        if len(by_donor[did]) > 1:
            raise ValueError(f"{did}: diagnosis disagreement across nuclei")
    missing = sorted(set(known) - set(sample_to_donor))
    undocumented = sorted(set(missing) - {"C-612", "C-13291"})
    if undocumented:
        raise ValueError(f"Unexpected absent donor–region IDs: {undocumented}")
    return {
        "included_sample_ids": len(sample_to_donor),
        "donors": len(by_donor),
        "missing_geo_sample_ids": missing,
        "documented_failed_caudate_samples": sorted(set(missing) & {"C-612", "C-13291"}),
        "sample_ids": sorted(sample_to_donor),
    }


def numeric_summary(values):
    a = np.asarray(values, dtype=float)
    finite = a[np.isfinite(a)]
    if not len(finite):
        return {"n_valid": 0, "median": None, "p05": None, "p95": None}
    return {"n_valid": int(len(finite)), "median": round(float(np.median(finite)), 4),
            "p05": round(float(np.percentile(finite, 5)), 4),
            "p95": round(float(np.percentile(finite, 95)), 4)}


# %% Audit only existing author-retained obs metadata
def analyze_obs(fields, geo):
    length = len(fields["ID"])
    if any(len(v) != length for v in fields.values()):
        raise ValueError("obs field lengths differ")
    design = assert_identifiers(fields, geo)
    categories = defaultdict(list)
    for i in range(length):
        key = (str(fields["Case"][i]), str(fields["ID"][i]),
               str(fields["Region"][i]), str(fields["Dx_OUD"][i]),
               str(fields["celltype1"][i]))
        categories[key].append(i)
    flag_cols = [name for name in fields if name.endswith(".keep")]
    metric_cols = [name for name in ("nCount_RNA", "nFeature_RNA", "percent.mt")
                   if name in fields]
    coverage = []
    flag_rows = []
    for key, indices in sorted(categories.items()):
        donor, sample, region, diagnosis, celltype = key
        entry = {"donor": donor, "sample_id": sample, "region": region,
                 "diagnosis": diagnosis, "celltype1": celltype, "n_nuclei": len(indices)}
        for metric in metric_cols:
            stats = numeric_summary(fields[metric][indices])
            stem = {"nCount_RNA": "umi", "nFeature_RNA": "detected_genes",
                    "percent.mt": "pct_mt"}[metric]
            entry[f"{stem}_median"] = stats["median"]
            entry[f"{stem}_p05"] = stats["p05"]
            entry[f"{stem}_p95"] = stats["p95"]
        coverage.append(entry)
        for flag in flag_cols:
            vals = Counter(fmt_bool_flag(fields[flag][j]) for j in indices)
            flag_rows.append({"donor": donor, "sample_id": sample,
                              "celltype1": celltype, "flag": flag,
                              "true": vals["TRUE"], "false": vals["FALSE"],
                              "unknown": vals["UNKNOWN"], "n_nuclei": len(indices)})
    if sum(row["n_nuclei"] for row in coverage) != length:
        raise AssertionError("Nucleus coverage reconciliation failed")
    per_sample = Counter(str(x) for x in fields["ID"])
    for sample, count in per_sample.items():
        expected_count = int(next(r["nuclei"] for r in geo if r["sample_id"] == sample))
        if count != expected_count:
            raise ValueError(f"{sample}: nuclei count differs from inspection-based validation")
    condition_donors = defaultdict(set)
    for donor, _, region, diagnosis, celltype in categories:
        condition_donors[(region, celltype, diagnosis)].add(donor)
    balance = [
        {"region": reg, "celltype1": cell, "diagnosis": dx,
         "independent_donors": len(donors),
         "total_nuclei": sum(r["n_nuclei"] for r in coverage
                            if (r["region"], r["celltype1"], r["diagnosis"]) == (reg, cell, dx))}
        for (reg, cell, dx), donors in sorted(condition_donors.items())
    ]
    design.update({
        "n_nuclei": length, "n_donor_region_celltype_groups": len(coverage),
        "qc_flag_columns": flag_cols,
        "metrics_present": metric_cols,
        "flag_semantics": "Descriptive only; .keep interpretation requires author source validation.",
        "source_population": "AUTHOR_RETAINED_NUCLEI_ONLY",
        "automatic_nucleus_exclusions": 0,
        "pseudobulk_or_de_performed": False,
    })
    return coverage, flag_rows, balance, design


def run(h5_path, inspection, geo):
    with h5py.File(h5_path, "r") as h5:
        if "obs" not in h5 or "raw/X" not in h5:
            raise ValueError("Expected processed AnnData object with obs and raw/X")
        raw = h5["raw/X"]
        shape = tuple(int(x) for x in raw.attrs["shape"])
        expected = tuple(int(x) for x in inspection["expression"]["raw/X"]["shape"])
        if shape != expected:
            raise ValueError(f"Raw count layer dimensions changed: {shape} vs {expected}")
        obs = h5["obs"]
        missing = sorted(set(REQUIRED) - set(obs.keys()))
        if missing:
            raise ValueError(f"Required obs fields missing: {missing}")
        size = shape[0]
        fields = {col: obs_vector(obs, col, size)
                  for col in REQUIRED + OPTIONAL if col in obs}
        coverage, flag_rows, balance, design = analyze_obs(fields, geo)
        design["raw_count_shape"] = list(shape)
        design["obs_fields_present"] = sorted(obs.keys())
        design["transformed_X_not_used"] = True
        return coverage, flag_rows, balance, design


# %% Disk-safe access for gzip-wrapped HDF5
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", type=Path, default=None,
                        help="Optional already-decompressed H5AD (avoids 8.8 GB extraction)")
    args = parser.parse_args()
    inspection = json.loads(INSPECTION.read_text(encoding="utf-8"))
    geo = load_csv(GEO_TABLE)
    if args.h5ad is not None:
        if not args.h5ad.is_file():
            parser.error(f"Missing uncompressed H5AD: {args.h5ad}")
        coverage, flags, balance, summary = run(args.h5ad, inspection, geo)
        source = args.h5ad
    else:
        archives = sorted(ARCHIVE_DIR.glob("*.h5ad.gz"))
        if len(archives) != 1:
            parser.error("Expected exactly one GSE225158 H5AD.gz archive")
        source = archives[0]
        if shutil.disk_usage(ARCHIVE_DIR).free < MIN_FREE_BYTES:
            parser.error("Need at least 13 GB free for temporary H5AD extraction. "
                         "Use --h5ad with an already-decompressed file if available.")
        with tempfile.TemporaryDirectory(prefix="striatum-qc-", dir=ARCHIVE_DIR) as tmp:
            destination = Path(tmp) / "qc_tmp.h5ad"
            with gzip.open(source, "rb") as src, destination.open("wb") as dst:
                while block := src.read(CHUNK_BYTES):
                    if shutil.disk_usage(ARCHIVE_DIR).free < 2_000_000_000:
                        raise OSError("Disk space critically low during temporary H5AD extraction")
                    dst.write(block)
            coverage, flags, balance, summary = run(destination, inspection, geo)
    OUT.mkdir(parents=True, exist_ok=True)
    write_csv(OUT / "striatum_donor_celltype_qc.csv", coverage,
              ["donor", "sample_id", "region", "diagnosis", "celltype1", "n_nuclei",
               "umi_median", "umi_p05", "umi_p95",
               "detected_genes_median", "detected_genes_p05", "detected_genes_p95",
               "pct_mt_median", "pct_mt_p05", "pct_mt_p95"])
    write_csv(OUT / "striatum_inherited_flags.csv", flags,
              ["donor", "sample_id", "celltype1", "flag", "true", "false", "unknown", "n_nuclei"])
    write_csv(OUT / "striatum_donor_coverage.csv", balance,
              ["region", "celltype1", "diagnosis", "independent_donors", "total_nuclei"])
    summary.update({"dataset": "GSE225158",
                    "original_input": str(source),
                    "source_sha256": sha256(source),
                    "gene_symbol_mapping_audited": False,
                    "raw_count_entries_rescanned": False,
                    "raw_count_verification_reused": "outputs/inspection/GSE225158_raw_verification.json",
                    "interpretation": "Inherited author-selected QC; further filtering requires methods review."})
    (OUT / "striatum_qc_readiness.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"GSE225158: {summary['n_nuclei']} retained nuclei, "
          f"{summary['donors']} donors, {summary['included_sample_ids']} donor-region samples.")
    print("Wrote read-only striatal QC outputs. No exclusions, pseudobulk or DE.")


if __name__ == "__main__":
    main()
