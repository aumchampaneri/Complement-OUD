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
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

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
    print("Reports:", REPORT_DIR.relative_to(PROJECT_DIR))
    print("Original downloaded files were not modified.")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
