# %% Temporary human transcriptomics inspection
# > Read-only data audit for all four acquired human GEO studies.
# > Run: pixi run python scripts/_inspect_human_data.py
# > Optional: --datasets GSE225158 --all-h5 --sample-values 100000
# > Source files are never modified; reports are written to ignored outputs/.
#
# %% PATH SETUP
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_DIR / "data" / "raw" / "human"
REPORT_DIR = PROJECT_DIR / "outputs" / "inspection"

# %% IMPORTS
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
import logging
import math
import shutil
import tarfile
import tempfile

import h5py

LOG = logging.getLogger("human-inspection")
ACCESSIONS = ("GSE174409", "GSE225158", "GSE194368", "GSE240457")
KEY_FIELDS = ("id", "case", "donor", "subject", "orig.ident", "pair", "dx_oud",
              "region", "sex", "celltype1", "celltype2", "celltype3", "sample")
MAX_CATEGORIES = 60


# %% REPORT HELPERS
def display(value):
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if hasattr(value, "item"):
        return str(value.item())
    return str(value)


def save_json(accession, report):
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    path = REPORT_DIR / f"{accession}_inspection.json"
    with path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False, default=display)
    LOG.info("Saved %s", path.relative_to(PROJECT_DIR))


def preview(values, limit=8):
    return [display(x) for x in values[:limit]]


def read_gzip_metadata(folder):
    path = folder / (folder.name + "_series_matrix.txt.gz")
    result = {"path": str(path.relative_to(PROJECT_DIR)), "samples": [], "series": {}}
    if not path.exists():
        result["error"] = "Series matrix not found"
        return result
    sample_arrays = defaultdict(list)
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as source:
        for line in source:
            if not line.startswith("!"):
                continue
            row = next(csv.reader([line.rstrip("\n\r")], delimiter="\t"))
            if len(row) < 2:
                continue
            key = row[0]
            if key.startswith("!Sample_"):
                sample_arrays[key].append(row[1:])
            elif key in ("!Series_title", "!Series_summary", "!Series_overall_design"):
                result["series"][key] = row[1][:1500]
    widths = sorted({len(v) for blocks in sample_arrays.values() for v in blocks})
    result["sample_widths"] = widths
    n = max(widths, default=0)
    samples = [dict() for _ in range(n)]
    for key, blocks in sample_arrays.items():
        label = key.removeprefix("!Sample_")
        for index, values in enumerate(blocks):
            for sample_index, value in enumerate(values):
                if sample_index >= n:
                    continue
                if label == "characteristics_ch1":
                    if ":" in value:
                        field, detail = value.split(":", 1)
                        samples[sample_index][field.strip().lower().replace(" ", "_")] = detail.strip()
                    else:
                        samples[sample_index][f"characteristic_{index + 1}"] = value
                else:
                    samples[sample_index][label] = value
    result["samples"] = samples
    result["n_geo_samples"] = n
    result["metadata_fields"] = sorted(set().union(*(s.keys() for s in samples))) if samples else []
    return result


def write_sample_table(accession, samples):
    if not samples:
        return
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    fields = sorted(set().union(*(s.keys() for s in samples)))
    path = REPORT_DIR / f"{accession}_geo_samples.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(samples)


# %% BULK COUNT-MATRIX AUDIT
def bulk_matrix(path):
    result = {"file": path.name, "bytes": path.stat().st_size}
    opener = gzip.open if path.suffix == ".gz" else open
    sep = "\t" if ".tsv" in path.name or ".txt" in path.name else ","
    samples = []
    totals = []
    detected = []
    rows = 0
    format_errors = 0
    missing = 0
    noninteger = 0
    negative = 0
    gene_examples = []
    genes_seen = set()
    repeated_genes = 0

    with opener(path, "rt", encoding="utf-8-sig", errors="replace", newline="") as source:
        reader = csv.reader(source, delimiter=sep)
        header = next(reader, [])
        if len(header) < 2:
            result["error"] = "Count matrix header has fewer than two columns"
            return result
        samples = header[1:]
        totals = [0] * len(samples)
        detected = [0] * len(samples)
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            rows += 1
            if len(row) != len(header):
                format_errors += 1
                continue
            gene = row[0]
            if len(gene_examples) < 12:
                gene_examples.append(gene)
            if gene in genes_seen:
                repeated_genes += 1
            genes_seen.add(gene)
            for i, value in enumerate(row[1:]):
                if not value:
                    missing += 1
                    continue
                try:
                    number = float(value)
                    if not math.isfinite(number) or number != int(number):
                        noninteger += 1
                    if number < 0:
                        negative += 1
                    if number >= 0 and number == int(number):
                        totals[i] += int(number)
                        detected[i] += int(number > 0)
                except ValueError:
                    noninteger += 1
    result.update({"genes": rows, "samples": len(samples), "sample_names": samples,
                   "gene_examples": gene_examples, "duplicate_gene_ids": repeated_genes,
                   "malformed_rows": format_errors, "missing_entries": missing,
                   "noninteger_or_nonnumeric_entries": noninteger, "negative_entries": negative,
                   "per_sample": [{"sample": x, "library_size": totals[i],
                                   "genes_detected": detected[i]} for i, x in enumerate(samples)]})
    return result


# %% HDF5 METADATA AND SPARSE-MATRIX AUDIT
def h5_vector(node, length=None):
    """Small metadata columns only, not expression matrices."""
    if isinstance(node, h5py.Group) and "categories" in node and "codes" in node:
        categories = node["categories"].asstr()[:] if h5py.check_string_dtype(node["categories"].dtype) else node["categories"][:]
        codes = node["codes"][:]
        return [display(categories[int(i)]) if i >= 0 else "<missing>" for i in codes]
    if isinstance(node, h5py.Dataset) and len(node.shape) == 1:
        if length is not None and len(node) > length:
            return []
        return [display(v) for v in node[:]]
    return []


def matrix_descriptor(node, sample_values):
    result = {"encoding": display(node.attrs.get("encoding-type", "unknown"))}
    if isinstance(node, h5py.Dataset):
        result.update({"storage": "dense", "shape": list(node.shape), "dtype": str(node.dtype)})
        # Small first block only.
        vals = node[:min(20, node.shape[0]), :min(20, node.shape[1])].ravel() if len(node.shape) == 2 else node[:min(sample_values, len(node))]
    elif isinstance(node, h5py.Group) and "data" in node:
        data = node["data"]
        shape = node.attrs.get("shape", None)
        result.update({"storage": "sparse", "shape": list(shape) if shape is not None else None,
                       "dtype": str(data.dtype), "stored_values": len(data)})
        vals = data[:min(sample_values, len(data))]
    else:
        return dict(result, warning="Unrecognized matrix representation")
    import numpy as np
    arr = np.asarray(vals)
    if arr.size:
        valid = np.isfinite(arr)
        result["sampled_values"] = int(arr.size)
        result["sample_finite_fraction"] = float(valid.mean())
        result["sample_min"] = float(np.nanmin(arr))
        result["sample_max"] = float(np.nanmax(arr))
        result["sample_nonnegative_integer_fraction"] = float(
            (valid & (arr >= 0) & (arr == np.floor(arr))).mean())
        result["raw_integer_counts_confirmed"] = False
        result["interpretation"] = ("sample is integer-like; full count provenance still required"
                                    if result["sample_nonnegative_integer_fraction"] == 1
                                    else "sample contains transformed or non-count-like values")
    return result


def inspect_h5ad(archive, sample_values):
    folder = archive.parent
    result = {"file": archive.name, "bytes": archive.stat().st_size}
    # Compressed HDF5 cannot be randomly accessed by h5py.
    if shutil.disk_usage(folder).free < 14_000_000_000:
        result["error"] = "Less than 14 GB free; temporary extraction skipped"
        return result
    with tempfile.TemporaryDirectory(prefix="human-h5ad-", dir=folder) as tmp:
        target = Path(tmp) / "inspection.h5ad"
        LOG.info("Temporarily decompressing %s (approximately 8.8 GB)", archive.name)
        with gzip.open(archive, "rb") as source, target.open("wb") as output:
            while chunk := source.read(8 * 1024 * 1024):
                if shutil.disk_usage(folder).free < 2_000_000_000:
                    raise OSError("Low free disk space during H5AD extraction")
                output.write(chunk)
        result["uncompressed_bytes"] = target.stat().st_size
        with h5py.File(target, "r") as h5:
            result["root_keys"] = list(h5.keys())
            for key in ("X", "raw/X"):
                if key in h5:
                    result[key] = matrix_descriptor(h5[key], sample_values)
            result["layers"] = list(h5["layers"].keys()) if "layers" in h5 else []
            for key in ("obs", "var"):
                if key in h5:
                    result[f"{key}_columns"] = list(h5[key].keys())
            if "var" in h5:
                idx = display(h5["var"].attrs.get("_index", "_index"))
                if idx in h5["var"]:
                    result["gene_examples"] = preview(h5_vector(h5["var"][idx]), 20)
            if "raw/var" in h5:
                idx = display(h5["raw/var"].attrs.get("_index", "_index"))
                if idx in h5["raw/var"]:
                    result["raw_gene_examples"] = preview(h5_vector(h5["raw/var"][idx]), 20)

            if "obs" in h5:
                obs = h5["obs"]
                columns = {}
                for key in obs.keys():
                    if key.lower() in KEY_FIELDS:
                        values = h5_vector(obs[key])
                        if values:
                            columns[key] = values
                result["metadata_columns_inspected"] = list(columns)
                for key, values in columns.items():
                    counts = Counter(values)
                    result.setdefault("metadata_distributions", {})[key] = {
                        "n_unique": len(counts),
                        "top_categories": counts.most_common(MAX_CATEGORIES)}
                if columns:
                    n = len(next(iter(columns.values())))
                    out = REPORT_DIR / "GSE225158_nuclei_metadata.csv"
                    REPORT_DIR.mkdir(parents=True, exist_ok=True)
                    with out.open("w", newline="", encoding="utf-8") as handle:
                        writer = csv.writer(handle)
                        writer.writerow(list(columns))
                        for i in range(n):
                            writer.writerow([columns[key][i] for key in columns])
                    result["nuclei_metadata_table"] = str(out.relative_to(PROJECT_DIR))
                donor = next((k for k in ("ID", "Case", "orig.ident") if k in columns), None)
                if donor:
                    result["candidate_donor_column"] = donor
                    for field in ("Dx_OUD", "Region", "celltype1", "Sex"):
                        if field in columns:
                            summary = Counter(zip(columns[donor], columns[field]))
                            result[f"per_{donor}_{field}_counts"] = [
                                {"donor": d, "category": c, "n_nuclei": count}
                                for (d, c), count in sorted(summary.items())]
                result["donor_id_is_verified"] = False
                result["caveat"] = ("Candidate donor IDs and metadata are exploratory; "
                                    "reconcile with GEO and source publication before inference.")
    return result


# %% TAR ARCHIVE: PER-LIBRARY 10X MATRICES
def inspect_10x_h5(path):
    with h5py.File(path, "r") as h5:
        base = h5["matrix"] if "matrix" in h5 else h5
        result = {"keys": list(base.keys())}
        if "shape" in base:
            result["matrix_shape"] = [int(x) for x in base["shape"][:]]
        if "barcodes" in base:
            result["barcode_count"] = len(base["barcodes"])
            result["barcode_examples"] = preview(base["barcodes"][:min(4, len(base["barcodes"]))])
        if "data" in base:
            result["stored_nonzero_entries"] = len(base["data"])
            result["count_sample"] = matrix_descriptor(base["data"], 10000)
        if "features" in base:
            features = base["features"]
            result["feature_fields"] = list(features.keys())
            for key in ("id", "name", "feature_type", "genome"):
                if key in features:
                    result[f"{key}_examples"] = preview(features[key][:min(8, len(features[key]))])
        return result


def inspect_tar(archive, max_h5, sample_values):
    result = {"file": archive.name, "bytes": archive.stat().st_size,
              "libraries": [], "n_libraries": 0}
    with tarfile.open(archive, "r") as source:
        members = [m for m in source.getmembers() if m.isfile() and m.name.endswith(".h5")]
        result["n_libraries"] = len(members)
        result["library_filenames"] = [m.name for m in members]
        for member in members[:max_h5]:
            if member.size > 20_000_000_000:
                result.setdefault("warnings", []).append(f"Skipped oversized member: {member.name}")
                continue
            LOG.info("Inspecting library %s", member.name)
            with tempfile.TemporaryDirectory(prefix="human-10x-") as tmp:
                target = Path(tmp) / "library.h5"
                source_file = source.extractfile(member)
                if source_file is None:
                    continue
                with source_file, target.open("wb") as output:
                    shutil.copyfileobj(source_file, output, length=8 * 1024 * 1024)
                details = inspect_10x_h5(target)
            details.update({"filename": member.name, "bytes": member.size})
            result["libraries"].append(details)
    result["libraries_inspected"] = len(result["libraries"])
    result["donor_mapping_verified"] = False
    result["caveat"] = ("GSM/H5 filenames represent sequencing libraries, not necessarily "
                       "independent donors. Donor demultiplexing and diagnosis metadata "
                       "must be independently verified.")
    return result



# %% CROSS-STUDY PROVENANCE RECONCILIATION
# > Operates on existing JSON reports; it does not load count matrices.
# > All matching across studies is PROVISIONAL and requires a source crosswalk.
def candidate_matches(bulk_samples, nucleus_samples):
    """Rank informative clinical matches; never assert donor equivalence."""
    from itertools import product

    def norm(s):
        return " ".join(str(s or "").strip().casefold().split())

    # Age/sex/diagnosis are necessary but insufficient to identify a donor.
    # Exact medication/toxicology and tissue-quality concordance add evidence.
    fields = (("age", "age", 1), ("sex", "sex", 1),
              ("diagnosis", "dx_oud", 2),
              ("bloodtexicology", "blood.toxicology", 5),
              ("tissuestoragetime_month", "tissue.storage.time.mo", 2))
    matches = []
    for old, new in product(bulk_samples, nucleus_samples):
        agreements, conflicts, score = [], [], 0
        for old_key, new_key, weight in fields:
            x, y = norm(old.get(old_key)), norm(new.get(new_key))
            if old_key == "sex":
                x = {"female": "f", "male": "m"}.get(x, x)
                y = {"female": "f", "male": "m"}.get(y, y)
            if old_key == "diagnosis":
                x = "control" if x in ("cont", "none", "unaffected") else x
                y = "control" if y in ("none", "cont", "unaffected") else y
            if not x or not y:
                continue
            if x == y:
                agreements.append(old_key)
                score += weight
            else:
                conflicts.append(old_key)
        for old_key, new_key, weight, tolerance in (
            ("pmi", "pmi", 3, 0.12),
            ("rpf_rin", "rin", 2, 0.11),
            ("avgph", "ph", 2, 0.11),
        ):
            try:
                x, y = float(old[old_key]), float(new[new_key])
            except (KeyError, ValueError, TypeError):
                continue
            if abs(x-y) <= tolerance:
                agreements.append(old_key)
                score += weight
            else:
                conflicts.append(old_key)
        # A source-verified crosswalk is required regardless of numerical score.
        if (score >= 9 and "diagnosis" in agreements
                and "age" in agreements and "sex" in agreements
                and not conflicts):
            matches.append({
                "bulk_title": old.get("title", ""),
                "bulk_geo": old.get("geo_accession", ""),
                "bulk_region": old.get("source_name_ch1", ""),
                "nucleus_subject": new.get("subject", ""),
                "nucleus_title": new.get("title", ""),
                "nucleus_geo": new.get("geo_accession", ""),
                "nucleus_region": new.get("tissue", ""),
                "score": score,
                "agreements": ";".join(agreements),
                "status": "CANDIDATE_NOT_VERIFIED",
            })
    # Keep one row per bulk sample and nucleus subject, not two regions
    # of the same nucleus subject.
    output = {}
    for item in sorted(matches, key=lambda x: -x["score"]):
        k = (item["bulk_title"], item["nucleus_subject"])
        if k not in output:
            output[k] = item
    return list(output.values())


def reconcile_existing_reports():
    """Compare GEO metadata and acquired matrices without rereading large data."""
    reports = {}
    for accession in ACCESSIONS:
        path = REPORT_DIR / f"{accession}_inspection.json"
        if path.exists():
            try:
                reports[accession] = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, ValueError) as exc:
                LOG.warning("Cannot read %s: %s", path, exc)

    summary = {"available_reports": sorted(reports),
               "missing_reports": sorted(set(ACCESSIONS) - set(reports)),
               "checks": {}, "limitations": []}
    for accession, report in reports.items():
        geo = report.get("geo", {}).get("samples", [])
        expression = report.get("expression", {})
        check = {"n_geo_samples": len(geo)}
        if accession in ("GSE174409", "GSE194368"):
            matrix = (expression[0] if isinstance(expression, list) and expression else {})
            names = set(matrix.get("sample_names", []))
            sample_title_to_id = {s.get("title"): s.get("geo_accession") for s in geo}
            all_ids = {s.get("geo_accession") for s in geo}
            all_titles = set(sample_title_to_id)
            compare_ids = len(names & all_ids) >= len(names & all_titles)
            expected = all_ids if compare_ids else all_titles
            check.update({"matrix_sample_count": len(names),
                          "matching_key": "geo_accession" if compare_ids else "title",
                          "missing_from_matrix": sorted(expected - names),
                          "matrix_columns_without_geo_match": sorted(names - expected),
                          "matrix_sample_names_unique": len(names) == len(matrix.get("sample_names", [])),
                          "gene_examples": matrix.get("gene_examples", [])[:5]})
            missing = [s for s in geo if s.get("geo_accession" if compare_ids else "title") in expected - names]
            check["missing_sample_details"] = [
                {k: row.get(k, "") for k in ("title", "geo_accession", "disease_state", "diagnosis", "source_name_ch1")}
                for row in missing]
        elif accession == "GSE225158":
            ids = expression.get("metadata_distributions", {}).get("ID", {}).get("top_categories", [])
            observed = {i[0] for i in ids}
            expected = {s.get("title") for s in geo}
            check.update({"n_observed_donor_region_ids": len(observed),
                          "n_unique_geo_subjects": len({s.get("subject") for s in geo}),
                          "missing_donor_region_ids": sorted(expected - observed),
                          "unregistered_donor_region_ids": sorted(observed - expected),
                          "raw_x": expression.get("raw/X", {}),
                          "raw_gene_examples": expression.get("raw_gene_examples", []),
                          "count_provenance_verified": False})
        elif accession == "GSE240457":
            libs = expression.get("libraries", [])
            check.update({"n_geo_libraries": len(geo),
                          "n_h5_libraries": expression.get("n_libraries"),
                          "n_libraries_inspected": len(libs),
                          "n_filtered_barcodes": sum(int(i.get("barcode_count", 0)) for i in libs),
                          "donor_demultiplexing_available": False,
                          "note": "Library IDs and filtered barcodes are not individual donor assignments"})
            path = REPORT_DIR / "GSE240457_library_inventory.csv"
            fields = ["geo_accession", "library", "filename", "barcodes", "genes",
                      "donor_mapping_status"]
            with path.open("w", newline="", encoding="utf-8") as fh:
                writer = csv.DictWriter(fh, fieldnames=fields)
                writer.writeheader()
                for lib in libs:
                    filename = lib.get("filename", "")
                    geo_id = filename.split("_", 1)[0]
                    writer.writerow({"geo_accession": geo_id,
                                     "library": filename.split("_", 2)[1] if "_" in filename else "",
                                     "filename": filename, "barcodes": lib.get("barcode_count"),
                                     "genes": (lib.get("matrix_shape") or [None])[0],
                                     "donor_mapping_status": "UNRESOLVED"})
        summary["checks"][accession] = check

    if "GSE174409" in reports and "GSE225158" in reports:
        bulk = reports["GSE174409"]["geo"]["samples"]
        nuclei = reports["GSE225158"]["geo"]["samples"]
        matches = candidate_matches(bulk, nuclei)
        dest = REPORT_DIR / "cross_study_donor_candidates.csv"
        fields = ["bulk_title", "bulk_geo", "bulk_region", "nucleus_subject",
                  "nucleus_title", "nucleus_geo", "nucleus_region",
                  "score", "agreements", "status"]
        with dest.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields)
            writer.writeheader()
            writer.writerows(matches)
        summary["candidate_cross_study_matches"] = len(matches)
        summary["candidate_cross_study_donors"] = len({r["nucleus_subject"] for r in matches})
        summary["cross_study_status"] = "PROVISIONAL; NO AUTHORITATIVE CROSSWALK"

    summary["limitations"] = [
        "All cross-study donor matches are provisional, not identity proof.",
        "H5AD raw count provenance requires full matrix and RNA QC checks.",
        "GSE240457 donor-level testing requires verified barcode-to-donor annotations.",
        "GSE194368 NCBI counts may omit samples relative to submitter metadata.",
    ]
    path = REPORT_DIR / "human_reconciliation.json"
    path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    LOG.info("Saved %s", path.relative_to(PROJECT_DIR))
    for accession, check in summary["checks"].items():
        LOG.info("%s reconciliation: %s", accession, {
            k: v for k, v in check.items()
            if k.startswith("missing_") or k.startswith("n_")
        })
    return summary


# %% OPTIONAL DEEP H5AD COUNT-AND-ANNOTATION VERIFICATION
def verify_h5ad_raw(archive):
    """Full sparse-data integrity check; sample individual nuclei for QC agreement."""
    import numpy as np

    result = {"file": archive.name, "note": "Full sparse data scan; per-cell QC tested on distributed rows"}
    folder = archive.parent
    if shutil.disk_usage(folder).free < 14_000_000_000:
        return dict(result, error="Less than 14 GB free")
    with tempfile.TemporaryDirectory(prefix="human-verify-", dir=folder) as tmp:
        target = Path(tmp) / "counts.h5ad"
        with gzip.open(archive, "rb") as incoming, target.open("wb") as outgoing:
            shutil.copyfileobj(incoming, outgoing, length=8 * 1024 * 1024)
        with h5py.File(target, "r") as h5:
            node = h5["raw/X"]
            if not isinstance(node, h5py.Group) or "data" not in node:
                return dict(result, error="raw/X is not a recognized sparse matrix")
            data = node["data"]
            stats = {"nonfinite": 0, "negative": 0, "fractional": 0}
            for offset in range(0, len(data), 1_000_000):
                vals = np.asarray(data[offset:offset+1_000_000])
                stats["nonfinite"] += int((~np.isfinite(vals)).sum())
                stats["negative"] += int((vals < 0).sum())
                stats["fractional"] += int((np.isfinite(vals) & (vals != np.floor(vals))).sum())
            result["raw_data_entries_checked"] = len(data)
            result["raw_data_violations"] = stats
            result["all_raw_values_integer_nonnegative_finite"] = not any(stats.values())
            result["raw_var_columns"] = list(h5["raw/var"].keys())
            for field in ("_index", "features"):
                if field in h5["raw/var"]:
                    ds = h5["raw/var"][field]
                    result[f"raw_{field}_examples"] = preview(ds[:min(15, len(ds))])
                    result[f"raw_{field}_length"] = len(ds)
            indptr = node["indptr"]
            n_cells = len(indptr)-1
            indices = np.unique(np.linspace(0, n_cells - 1, min(1000, n_cells), dtype=np.int64))
            checks = {}
            for field, mode in (("nCount_RNA", "sum"), ("nFeature_RNA", "genes")):
                if field not in h5["obs"]:
                    continue
                obs = h5_vector(h5["obs"][field])
                errors = []
                for i in indices:
                    a, b = int(indptr[i]), int(indptr[i+1])
                    value = float(np.asarray(data[a:b]).sum()) if mode == "sum" else int(b-a)
                    try:
                        errors.append(abs(value-float(obs[int(i)])))
                    except (ValueError, IndexError):
                        pass
                checks[field] = {"n_tested": len(errors),
                                 "n_disagree": sum(v > 0.01 for v in errors),
                                 "maximum_absolute_error": max(errors, default=None)}
            result["qc_field_comparison"] = checks
            result["provenance_note"] = ("Count-like representation and QC concordance do not "
                                         "independently prove original preprocessed-count provenance.")
    dest = REPORT_DIR / "GSE225158_raw_verification.json"
    dest.write_text(json.dumps(result, indent=2), encoding="utf-8")
    LOG.info("Saved %s", dest.relative_to(PROJECT_DIR))
    return result


# %% ENTRY POINT
def main():
    parser = argparse.ArgumentParser(description="Read-only human transcriptomics audit")
    parser.add_argument("--datasets", nargs="+", choices=ACCESSIONS,
                        default=list(ACCESSIONS))
    parser.add_argument("--all-h5", action="store_true",
                        help="Inspect all H5 libraries in GSE240457 TAR")
    parser.add_argument("--max-h5", type=int, default=100000,
                        help="Maximum TAR libraries to inspect (default: all)")
    parser.add_argument("--sample-values", type=int, default=10000,
                        help="Expression values sampled per HDF5 matrix")
    parser.add_argument("--reconcile-only", action="store_true",
                        help="Use previously generated JSON reports; skip all matrix reads")
    parser.add_argument("--verify-h5ad-raw", action="store_true",
                        help="Deep scan GSE225158 raw/X and test QC field agreement (extra 8.8 GB extraction)")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    if args.reconcile_only:
        reconcile_existing_reports()
        if args.verify_h5ad_raw:
            archives = sorted((DATA_DIR / "GSE225158").glob("*.h5ad.gz"))
            if archives:
                verify_h5ad_raw(archives[0])
            else:
                LOG.warning("Cannot verify GSE225158 raw counts: archive not found")
                return 1
        return 0

    failures = 0
    for accession in args.datasets:
        LOG.info("=== %s ===", accession)
        folder = DATA_DIR / accession
        report = {"accession": accession}
        if not folder.exists():
            report["error"] = "Missing downloaded dataset directory"
            failures += 1
            save_json(accession, report)
            continue
        try:
            report["geo"] = read_gzip_metadata(folder)
            write_sample_table(accession, report["geo"].get("samples", []))
            if accession in ("GSE174409", "GSE194368"):
                candidates = list(folder.glob("*counts*.gz")) + list(folder.glob("*count*.csv"))
                report["expression"] = [bulk_matrix(path) for path in sorted(set(candidates))]
                if not candidates:
                    report["warning"] = "No count matrix found"
            elif accession == "GSE225158":
                archives = list(folder.glob("*.h5ad.gz"))
                report["expression"] = inspect_h5ad(archives[0], args.sample_values) if archives else {"error": "H5AD archive not found"}
            elif accession == "GSE240457":
                archives = list(folder.glob("*.tar"))
                report["expression"] = inspect_tar(
                    archives[0], 999999 if args.all_h5 else max(0, args.max_h5),
                    args.sample_values) if archives else {"error": "TAR archive not found"}
        except (OSError, ValueError, KeyError, RuntimeError) as exc:
            LOG.exception("Inspection failed for %s", accession)
            report["error"] = f"{type(exc).__name__}: {exc}"
            failures += 1
        save_json(accession, report)
        expression = report.get("expression", {})
        LOG.info("Status %s: %s", accession, "ERROR" if "error" in report else "report saved")
    reconcile_existing_reports()
    if args.verify_h5ad_raw:
        archive_dir = DATA_DIR / "GSE225158"
        archives = sorted(archive_dir.glob("*.h5ad.gz"))
        if archives:
            verify_h5ad_raw(archives[0])
        else:
            LOG.warning("Cannot verify GSE225158 raw counts: archive not found")
    print("Reports:", REPORT_DIR.relative_to(PROJECT_DIR))
    print("Original downloaded files were not modified.")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
