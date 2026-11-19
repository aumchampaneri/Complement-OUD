# %% Phase 2b read-only multivariate bulk QC tests
"""Run: pixi run python -m unittest discover -s tests -p 'test_qc_bulk*.py' -v"""
from pathlib import Path
import csv
import importlib.util
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("bulk_multi", ROOT / "scripts" / "02b_qc_bulk_multivariate.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class MultivariateTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name) / "counts.csv"
        self.rows = [
            {"dataset": "GSE174409", "matrix_sample": f"S{i + 1}",
             "assigned_gene_counts": n}
            for i, n in enumerate((114, 33, 58, 80))
        ]
        with self.path.open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["gene", "S1", "S2", "S3", "S4"])
            w.writerows([
                ["g1", 10, 20, 30, 10],
                ["g2", 100, 2, 8, 40],
                ["g3", 1, 10, 10, 10],
                ["g4", 3, 1, 10, 20],
            ])

    def test_read_only_pca_and_correlation(self):
        previous = self.path.read_bytes()
        pca, correlation, selected, report = mod.audit_cohort(
            self.path, self.rows, "GSE174409", top_n=3)
        self.assertEqual(self.path.read_bytes(), previous)
        self.assertEqual(len(pca), 4)
        self.assertEqual(len(correlation), 6)
        self.assertEqual(len(selected), 3)
        self.assertEqual(report["exclusions"], 0)
        self.assertEqual(report["mitochondrial_fraction_status"],
                         "NOT_CALCULATED_NO_VERIFIED_GENE_MAPPING")
        self.assertTrue(all(row["mt_count_fraction"] == "" for row in pca))
        self.assertTrue(all(-1.000001 <= row["profile_correlation"] <= 1.000001
                            for row in correlation))

    def test_source_order_guard(self):
        with self.path.open(newline="", encoding="utf-8") as f:
            rows = list(csv.reader(f))
        rows[0][1], rows[0][2] = rows[0][2], rows[0][1]
        with self.path.open("w", newline="", encoding="utf-8") as f:
            csv.writer(f).writerows(rows)
        with self.assertRaisesRegex(ValueError, "source header differs"):
            mod.audit_cohort(self.path, self.rows, "GSE174409")

    def test_annotation_duplicates_rejected(self):
        annotation = Path(self.tmp.name) / "map.csv"
        annotation.write_text(
            "dataset,feature_id,gene_symbol,chromosome\n"
            "GSE174409,g1,A,MT\nGSE174409,g1,B,MT\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "duplicate"):
            mod.load_annotation(annotation, "GSE174409")

    def test_partial_mito_annotation_explicitly_provisional(self):
        pca, _, _, report = mod.audit_cohort(
            self.path, self.rows, "GSE174409",
            annotation={"g1": {"chromosome": "MT"}})
        self.assertEqual(report["mitochondrial_mapped_features"], 1)
        self.assertIn("PROVISIONAL", report["mitochondrial_fraction_status"])
        self.assertGreater(pca[0]["mt_count_fraction"], 0)


if __name__ == "__main__":
    unittest.main()
