# Complement-OUD

A human-first investigation of complement-associated transcriptional changes in opioid use disorder (OUD), using public brain transcriptomic data.

We will first acquire, validate, and explore eligible **human** datasets. Mouse analyses are a possible downstream extension, guided by the strength and relevance of the human findings.

Complement-related RNA expression is not direct evidence of biochemical complement activation.

## Reproducibility

[Pixi](https://pixi.sh/) manages the computational environment. Dependencies and tasks will be added only as they become necessary.

```bash
pixi install
```

Work in progress: `architecture/refactor-v1`. Historical analyses remain available on `main`.

## Human data acquisition

The single acquisition entrypoint is `scripts/01a_download.py` (Python standard library only).

```bash
pixi run python scripts/01a_download.py --list
pixi run python scripts/01a_download.py
pixi run python scripts/01a_download.py --datasets GSE174409
```

The script discovers GEO supplementary files, selects count data for GSE174409 and an AnnData file for GSE225158 when available, retrieves GEO series metadata, and writes SHA-256 checksums to `data/metadata/acquisition/human_acquisition_manifest.csv`. Large downloads and manifests stay in ignored `data/`. Use `--all-supplementary` only after reviewing `--list`; it may download large archives. Verify each selected file and its count layer before analysis. This downloader has not yet been tested against live GEO endpoints.
