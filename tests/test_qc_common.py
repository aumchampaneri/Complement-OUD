"""Regression tests for shared read-only QC I/O helpers."""
from __future__ import annotations

import csv
import hashlib
from pathlib import Path
import tempfile
import unittest

from scripts._qc_common import read_csv, sha256, write_csv


class QCCommonTests(unittest.TestCase):
    def test_csv_round_trip_and_column_order(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "sample.csv"
            rows = [{"sample": "C1", "count": 2, "extra": "not_exported"},
                    {"sample": "O1", "count": 7, "extra": "not_exported"}]
            write_csv(path, rows, ["sample", "count"])
            with path.open(newline="", encoding="utf-8") as handle:
                self.assertEqual(next(csv.reader(handle)), ["sample", "count"])
            self.assertEqual(read_csv(path),
                             [{"sample": "C1", "count": "2"},
                              {"sample": "O1", "count": "7"}])

    def test_sha256_matches_standard_library_and_does_not_mutate(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "counts.bin"
            payload = b"test\x00counts" * 100000
            path.write_bytes(payload)
            before = path.read_bytes()
            self.assertEqual(sha256(path), hashlib.sha256(payload).hexdigest())
            self.assertEqual(path.read_bytes(), before)


if __name__ == "__main__":
    unittest.main()
