# %% Phase 2 bulk QC regression tests (standard library only)
"""Run: pixi run python -m unittest discover -s tests -p 'test_qc_bulk.py' -v"""

from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "scripts" / "02a_qc_bulk.py"
spec = importlib.util.spec_from_file_location("human_bulk_qc", MODULE)
qc = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(qc)


class BulkQcTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.file = Path(self.tmp.name) / "counts.csv"
        self.samples = [f"S{i}" for i in range(6)]
        self.genes = [
            [10, 11, 12, 13, 14, 15],
            [0, 1, 0, 2, 0, 3],
            [5, 0, 7, 0, 9, 0],
        ]
        self._write(self.genes)
        expected = []
        for i, sample in enumerate(self.samples):
            vals = [row[i] for row in self.genes]
            expected.append({"sample": sample, "library_size": sum(vals),
                             "genes_detected": sum(v > 0 for v in vals)})
        self.report = {
            "geo": {"samples": [
                {"title": sample, "geo_accession": f"GSM{i}",
                 "source_name_ch1": "NAC",
                 "diagnosis": "OUD" if i % 2 else "CONT"}
                for i, sample in enumerate(self.samples)]},
            "expression": [{
                "file": self.file.name,
                "genes": len(self.genes),
                "sample_names": self.samples,
                "per_sample": expected,
            }]
        }

    def _write(self, matrix):
        with self.file.open("w", newline="", encoding="utf-8") as out:
            writer = csv.writer(out)
            writer.writerow(["gene"] + self.samples)
            for index, row in enumerate(matrix):
                writer.writerow([f"G{index}"] + row)

    def test_qc_is_read_only_and_matches_inspector(self):
        before = self.file.read_bytes()
        rows, metadata = qc.audit_matrix(self.file, "GSE174409", self.report)
        self.assertEqual(self.file.read_bytes(), before)
        self.assertEqual(metadata["n_samples"], 6)
        self.assertEqual(metadata["n_features"], 3)
        self.assertEqual(rows[0]["assigned_gene_counts"], 15)
        self.assertEqual(rows[0]["detected_features"], 2)
        self.assertAlmostEqual(rows[0]["top10_count_fraction"], 1.0)
        self.assertEqual(rows[0]["technical_exclusion"], "NONE_APPLIED")
        self.assertFalse(metadata["flags_are_exclusions"])
        self.assertEqual(rows[1]["diagnosis"], "OUD")

    def test_rejects_mismatched_sample_order(self):
        with self.file.open("r", newline="", encoding="utf-8") as inp:
            rows = list(csv.reader(inp))
        rows[0][1], rows[0][2] = rows[0][2], rows[0][1]
        with self.file.open("w", newline="", encoding="utf-8") as out:
            csv.writer(out).writerows(rows)
        with self.assertRaisesRegex(ValueError, "raw header differs"):
            qc.audit_matrix(self.file, "GSE174409", self.report)

    def test_rejects_fractional_counts(self):
        bad = [list(row) for row in self.genes]
        bad[0][0] = "10.5"
        self._write(bad)
        with self.assertRaisesRegex(ValueError, "invalid count"):
            qc.audit_matrix(self.file, "GSE174409", self.report)

    def test_rejects_disagreement_with_existing_inspection(self):
        self.report["expression"][0]["per_sample"][0]["library_size"] += 1
        with self.assertRaisesRegex(ValueError, "disagrees with inspection"):
            qc.audit_matrix(self.file, "GSE174409", self.report)

    def test_flags_not_exclusions(self):
        rows = [{"region": "NAc", "assigned_gene_counts": 100} for _ in range(6)]
        qc.describe_flags(rows, "assigned_gene_counts", grouped=True)
        self.assertTrue(all(r["assigned_gene_counts_flag"] == "UNDEFINED_MAD" for r in rows))


if __name__ == "__main__":
    unittest.main()
