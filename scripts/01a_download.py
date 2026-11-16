# %% Human transcriptomics data acquisition
# > Discover and download GEO supplementary files and series metadata.
# > Original source files remain unmodified; checksums document provenance.
# > Run from any directory: pixi run python scripts/01a_download.py --list

# %% PATH SETUP
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_DIR / "data" / "raw" / "human"
MANIFEST_DIR = PROJECT_DIR / "data" / "metadata" / "acquisition"

# %% IMPORTS
import argparse
import csv
from datetime import datetime, timezone
import hashlib
from html.parser import HTMLParser
import logging
import os
import re
import shutil
import sys
from urllib.error import HTTPError, URLError
from urllib.parse import unquote, urljoin, urlparse
from urllib.request import Request, urlopen

# %% DATASET REGISTRY
# Add human GEO series only after checking their assay and acquisition format.
DATASETS = {
    "GSE174409": {
        "assay": "bulk RNA-seq",
        "description": "Human DLPFC and nucleus accumbens; OUD/control",
        "file_pattern": r"raw.count|raw_count|counts.*\.csv",
    },
    "GSE225158": {
        "assay": "snRNA-seq",
        "description": "Human caudate and putamen; OUD/control",
        "file_pattern": r"\.h5ad(?:\.gz)?$",
    },
}

BASE_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series"
USER_AGENT = "Complement-OUD-acquisition/0.1 (public GEO data)"
CHUNK_SIZE = 1024 * 1024
LOGGER = logging.getLogger("acquisition")


# %% HELPERS
class LinkParser(HTMLParser):
    """Extract links from GEO's supplementary directory listing."""

    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            href = dict(attrs).get("href")
            if href:
                self.links.append(href)


def geo_paths(accession):
    if not re.fullmatch(r"GSE\d{3,}", accession):
        raise ValueError(f"Invalid GEO accession: {accession}")
    stem = accession[:-3] + "nnn"
    root = f"{BASE_URL}/{stem}/{accession}"
    return {
        "supplementary": f"{root}/suppl/",
        "soft": f"{root}/soft/{accession}_family.soft.gz",
        "series_matrix": f"{root}/matrix/{accession}_series_matrix.txt.gz",
    }


def request_url(url):
    return Request(url, headers={"User-Agent": USER_AGENT})


def discover_supplementary(url):
    """List source files without filtering for complement-related results."""
    with urlopen(request_url(url), timeout=60) as response:
        parser = LinkParser()
        parser.feed(response.read().decode("utf-8", errors="replace"))
    files = {}
    for href in parser.links:
        address = urljoin(url, href)
        parsed = urlparse(address)
        name = unquote(Path(parsed.path).name)
        if not name or name in {".", ".."} or parsed.scheme != "https":
            continue
        if parsed.netloc != "ftp.ncbi.nlm.nih.gov":
            continue
        if not parsed.path.startswith(urlparse(url).path) or parsed.path.endswith("/"):
            continue
        files[name] = address
    return dict(sorted(files.items()))


def sha256_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_file(url, destination, *, overwrite=False):
    """Stream to a temporary file and replace atomically after hashing."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and not overwrite:
        LOGGER.info("Existing file: %s (hashing)", destination.name)
        return {"status": "existing", "bytes": destination.stat().st_size,
                "sha256": sha256_file(destination)}

    temporary = destination.with_name(destination.name + ".part")
    temporary.unlink(missing_ok=True)
    LOGGER.info("Downloading: %s", destination.name)
    try:
        with urlopen(request_url(url), timeout=120) as response, temporary.open("wb") as output:
            shutil.copyfileobj(response, output, length=CHUNK_SIZE)
        if temporary.stat().st_size == 0:
            raise IOError(f"Empty download: {url}")
        checksum = sha256_file(temporary)
        length = temporary.stat().st_size
        os.replace(temporary, destination)
        return {"status": "downloaded", "bytes": length, "sha256": checksum}
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def save_manifest(rows):
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    path = MANIFEST_DIR / "human_acquisition_manifest.csv"
    fields = ["accession", "assay", "file", "source_url", "local_path",
              "bytes", "sha256", "status", "downloaded_at_utc"]
    existing = {}
    if path.exists():
        with path.open(newline="", encoding="utf-8") as handle:
            existing = {(r["accession"], r["file"]): r for r in csv.DictReader(handle)}
    for row in rows:
        existing[(row["accession"], row["file"])] = row
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(existing[key] for key in sorted(existing))
    LOGGER.info("Manifest: %s", path.relative_to(PROJECT_DIR))


# %% CLI AND ACQUISITION
def main(argv=None):
    parser = argparse.ArgumentParser(description="Acquire public human GEO expression data")
    parser.add_argument("--datasets", nargs="+", choices=sorted(DATASETS),
                        default=list(DATASETS), help="GEO accessions (default: all)")
    parser.add_argument("--list", action="store_true",
                        help="List remote files and the files selected for download")
    parser.add_argument("--all-supplementary", action="store_true",
                        help="Download every supplementary file (potentially very large)")
    parser.add_argument("--overwrite", action="store_true", help="Redownload existing files")
    parser.add_argument("--no-metadata", action="store_true", help="Skip SOFT and series matrix")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    rows = []
    failures = 0
    for accession in args.datasets:
        info = DATASETS[accession]
        paths = geo_paths(accession)
        LOGGER.info("%s: %s", accession, info["description"])
        try:
            available = discover_supplementary(paths["supplementary"])
        except (HTTPError, URLError, OSError) as exc:
            LOGGER.error("Cannot list %s files: %s", accession, exc)
            failures += 1
            continue

        selected = {name: url for name, url in available.items()
                    if args.all_supplementary or re.search(info["file_pattern"], name, re.I)}
        if not selected:
            LOGGER.warning("No supplementary files matched %s: inspect --list", accession)
            failures += 1
        for name in available:
            LOGGER.info("  %s %s", "SELECT" if name in selected else "skip  ", name)
        targets = list(selected.items())
        if not args.no_metadata:
            targets.extend((Path(urlparse(url).path).name, url)
                           for key, url in paths.items() if key != "supplementary")
        if args.list:
            continue
        for name, url in targets:
            target = DATA_DIR / accession / name
            try:
                result = download_file(url, target, overwrite=args.overwrite)
            except (HTTPError, URLError, OSError) as exc:
                LOGGER.error("Failed %s / %s: %s", accession, name, exc)
                failures += 1
                continue
            rows.append({"accession": accession, "assay": info["assay"],
                         "file": name, "source_url": url,
                         "local_path": str(target.relative_to(PROJECT_DIR)),
                         "bytes": result["bytes"], "sha256": result["sha256"],
                         "status": result["status"],
                         "downloaded_at_utc": datetime.now(timezone.utc).isoformat()})
    if rows:
        save_manifest(rows)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
