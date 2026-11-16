# %% Temporary human dataset inspection
# > Read-only exploratory audit; delete after reviewing the results.
# > Run: pixi run python scripts/_inspect_human_data.py
#
# > Requires h5py (pixi add h5py); does not load expression matrices into RAM.

# %% PATH SETUP
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent
HUMAN_DIR = PROJECT_DIR / "data" / "raw" / "human"

# %% IMPORTS
from collections import Counter
import csv
import gzip
import json
import shutil
import tempfile


# %% BULK MATRIX
def inspect_bulk():
    folder = HUMAN_DIR / "GSE174409"
    path = folder / "GSE174409_raw_counts_02102020.csv.gz"
    print("\n=== GSE174409: bulk RNA-seq ===")
    if not path.exists():
        print("MISSING:", path)
        return

    # The whole matrix is streamed. No normalization or filtering is performed.
    values = []
    rows = 0
    noninteger = 0
    width_errors = 0
    missing = 0
    with gzip.open(path, "rt", newline="", encoding="utf-8-sig") as source:
        reader = csv.reader(source)
        header = next(reader)
        print("Header columns:", len(header))
        print("First header entries:", header[:12])
        print("Likely sample columns:", max(0, len(header) - 1))
        for row in reader:
            rows += 1
            if len(row) != len(header):
                width_errors += 1
                continue
            if len(values) < 5:
                values.append({"feature": row[0], "example_counts": row[1:6]})
            for item in row[1:]:
                if item == "":
                    missing += 1
                elif not item.isdigit():  # flags negative/fractional/non-numeric data
                    noninteger += 1
    print("Feature rows:", rows)
    print("Malformed-width rows:", width_errors)
    print("Non-negative integer-format exceptions:", noninteger)
    print("Missing expression entries:", missing)
    print("Example rows:", json.dumps(values, ensure_ascii=False))

    # GEO series metadata may reveal study/sample identifiers and covariates.
    meta = folder / "GSE174409_series_matrix.txt.gz"
    if meta.exists():
        show_series_fields(meta)


# %% GEO METADATA
def show_series_fields(path):
    print("Series metadata fields (truncated):")
    labels = (
        "!Series_title", "!Series_summary",
        "!Sample_geo_accession", "!Sample_title",
        "!Sample_characteristics_ch1", "!Sample_source_name_ch1",
    )
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as source:
        shown = Counter()
        for line in source:
            if any(line.startswith(label + "\t") for label in labels):
                label = line.partition("\t")[0]
                if shown[label] < 3:
                    print(" ", line.rstrip()[:420])
                shown[label] += 1
        print("  Metadata row counts:", dict(shown))


# %% SINGLE-NUCLEUS H5AD
def text(value):
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def describe_columns(group, limit=90):
    print(f"  Columns ({len(group)}):", list(group.keys())[:limit])
    for name in group.keys():
        column = group[name]
        lname = name.lower()
        important = any(term in lname for term in (
            "donor", "subject", "patient", "sample", "diagnos",
            "disease", "condition", "region", "brain", "cell",
            "cluster", "sex", "batch", "age", "individual", "oud",
        ))
        if not important:
            continue
        if hasattr(column, "keys") and "categories" in column:
            categories = column["categories"][:]
            print(f"  {name}: categorical levels ({len(categories)}):",
                  [text(c) for c in categories[:20]])
            if "codes" in column:
                counts = Counter()
                codes = column["codes"]
                for offset in range(0, len(codes), 100000):
                    counts.update(int(c) for c in codes[offset:offset + 100000])
                pairs = sorted(
                    ((text(categories[i]) if i >= 0 else "<missing>", n)
                     for i, n in counts.items()), key=lambda x: -x[1]
                )
                print("   Counts:", pairs[:20])
        elif hasattr(column, "shape"):
            print(f"  {name}: shape={column.shape}, dtype={column.dtype}")
            if len(column.shape) == 1:
                print("   Example values:", [text(v) for v in column[:8]])
        else:
            print(f"  {name}: group keys={list(column.keys())[:10]}")


def inspect_h5ad():
    import h5py

    folder = HUMAN_DIR / "GSE225158"
    archives = sorted(folder.glob("*.h5ad.gz"))
    print("\n=== GSE225158: snRNA-seq ===")
    if not archives:
        print("MISSING: *.h5ad.gz in", folder)
        return

    archive = archives[0]
    free = shutil.disk_usage(folder).free
    print("Compressed archive:", archive.name)
    print("Archive bytes:", archive.stat().st_size)
    print("Available disk space (GB):", round(free / 1e9, 2))
    if free < 6_000_000_000:
        print("Not inspecting h5ad: <6 GB free. Free storage before extracting.")
        return

    # h5py cannot read a gzip-compressed HDF5 file directly. Extract into an
    # automatically removed temporary directory, without changing source data.
    with tempfile.TemporaryDirectory(prefix="human-h5ad-", dir=folder) as tmp:
        target = Path(tmp) / "inspection.h5ad"
        print("Temporarily extracting archive (may need multiple GB)...")
        with gzip.open(archive, "rb") as source, target.open("wb") as out:
            while chunk := source.read(8 * 1024 * 1024):
                if shutil.disk_usage(folder).free < 2_000_000_000:
                    raise OSError("Stopping extraction: <2 GB disk space remaining")
                out.write(chunk)
        print("Decompressed bytes:", target.stat().st_size)

        with h5py.File(target, "r") as h5:
            print("H5AD root keys:", list(h5.keys()))
            print("H5AD version:", text(h5.attrs.get("encoding-version", "not specified")))
            for key in ("X", "raw", "layers", "obsm", "uns"):
                if key not in h5:
                    continue
                node = h5[key]
                if isinstance(node, h5py.Dataset):
                    print(f"{key}: shape={node.shape}, dtype={node.dtype}")
                else:
                    print(f"{key}: keys={list(node.keys())[:30]}")
                    if key == "X":
                        print(" X sparse shape:", node.attrs.get("shape", "not set"))
                        if "data" in node:
                            data = node["data"]
                            sample = data[:min(10000, len(data))]
                            print(" X sample nonnegative integers:",
                                  bool(all(x >= 0 and x == int(x) for x in sample)))
                    if key == "layers":
                        for layer in list(node.keys())[:15]:
                            item = node[layer]
                            print("  layer:", layer, "shape:",
                                  getattr(item, "shape", item.attrs.get("shape", None)))
            if "obs" in h5:
                print("obs (cell metadata):")
                describe_columns(h5["obs"])
            if "var" in h5:
                var = h5["var"]
                print("var (gene metadata) columns:", list(var.keys())[:60])
                index_name = text(var.attrs.get("_index", "_index"))
                if index_name in var:
                    print("Gene IDs examples:", [text(x) for x in var[index_name][:12]])
            if "raw" in h5 and "var" in h5["raw"]:
                print("raw.var columns:", list(h5["raw"]["var"].keys())[:30])

    meta = folder / "GSE225158_series_matrix.txt.gz"
    if meta.exists():
        show_series_fields(meta)


# %% EXECUTION
def main():
    print("READ-ONLY HUMAN DATA INSPECTION")
    inspect_bulk()
    inspect_h5ad()
    print("\nDone. No source data were modified.")


if __name__ == "__main__":
    main()
