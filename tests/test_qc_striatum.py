# %% Tests for striatal inherited metadata QC (small synthetic H5AD fixture)
"""Run: pixi run python -m unittest discover -s tests -p 'test_qc_striatum.py' -v"""
import importlib.util
from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("striatum_qc", ROOT / "scripts" / "02d_qc_striatum.py")
qc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(qc)


class StriatumQCTests(unittest.TestCase):
    def setUp(self):
        d = tempfile.TemporaryDirectory()
        self.addCleanup(d.cleanup)
        self.path = Path(d.name) / "synthetic.h5ad"
        self.geo = [
            {"sample_id": "C-1", "subject_id": "1", "nuclei": "3", "condition": "OUD"},
            {"sample_id": "P-1", "subject_id": "1", "nuclei": "1", "condition": "OUD"},
        ]
        with h5py.File(self.path, "w") as h5:
            obs = h5.create_group("obs")
            for key, values in {
                "Case": [b"1"] * 4, "ID": [b"C-1"] * 3 + [b"P-1"],
                "Region": [b"Caudate"] * 3 + [b"Putamen"],
                "Dx_OUD": [b"OUD"] * 4,
                "celltype1": [b"Astrocyte", b"Astrocyte", b"Neuron", b"Neuron"],
                "nCount_RNA": [100, 200, 300, 400],
                "nFeature_RNA": [50, 80, 90, 120],
                "percent.mt": [0.5, 2.0, 0.8, 1.2],
                "miQC.keep": [True, False, True, True],
            }.items():
                obs.create_dataset(key, data=values)
            raw = h5.create_group("raw").create_group("X")
            raw.attrs["shape"] = [4, 5]
        self.report = {"expression": {"raw/X": {"shape": [4, 5]}}}

    def test_read_only_inherited_qc(self):
        before = self.path.read_bytes()
        coverage, flags, balance, design = qc.run(self.path, self.report, self.geo)
        self.assertEqual(self.path.read_bytes(), before)
        self.assertEqual(design["n_nuclei"], 4)
        self.assertEqual(design["donors"], 1)
        self.assertEqual(design["included_sample_ids"], 2)
        self.assertEqual(sum(r["n_nuclei"] for r in coverage), 4)
        self.assertEqual(sum(r["false"] for r in flags), 1)
        values = {(row["flag"], row["original_value"]): row["n_nuclei"]
                  for row in design["qc_flag_original_value_counts"]}
        self.assertEqual(values[("miQC.keep", "True")], 3)
        self.assertEqual(values[("miQC.keep", "False")], 1)
        self.assertEqual(design["automatic_nucleus_exclusions"], 0)
        self.assertEqual(len(balance), 3)

    def test_numeric_case_1034_point_zero_matches_geo_donor(self):
        self.assertEqual(qc.canonical_donor_id("1034.0"), "1034")
        self.assertEqual(qc.canonical_donor_id(1034.0), "1034")
        self.assertEqual(qc.canonical_donor_id("1034"), "1034")
        fields = {
            "Case": np.array(["1034.0", "1034.0"], dtype=object),
            "ID": np.array(["C-1034", "P-1034"], dtype=object),
            "Region": np.array(["Caudate", "Putamen"], dtype=object),
            "Dx_OUD": np.array(["OUD", "OUD"], dtype=object),
        }
        geo = [
            {"sample_id": "C-1034", "subject_id": "1034",
             "condition": "OUD", "nuclei": "1"},
            {"sample_id": "P-1034", "subject_id": "1034",
             "condition": "OUD", "nuclei": "1"},
        ]
        design = qc.assert_identifiers(fields, geo)
        self.assertEqual(design["donors"], 1)
        self.assertEqual(design["included_sample_ids"], 2)

    def test_reject_noninteger_donor_id(self):
        with self.assertRaisesRegex(ValueError, "Non-integer donor ID"):
            qc.canonical_donor_id("1034.5")
        with self.assertRaisesRegex(ValueError, "Non-integer donor ID"):
            qc.canonical_donor_id(float("nan"))

    def test_unknown_sample_id_rejected(self):
        with h5py.File(self.path, "r+") as h5:
            h5["obs/ID"][0] = b"C-9"
        with self.assertRaisesRegex(ValueError, "unknown donor"):
            qc.run(self.path, self.report, self.geo)

    def test_categorical_metadata(self):
        with h5py.File(self.path, "r+") as h5:
            del h5["obs/celltype1"]
            group = h5["obs"].create_group("celltype1")
            group.create_dataset("categories", data=[b"Neuron", b"Astrocyte"])
            group.create_dataset("codes", data=np.array([1, 1, 0, 0], dtype=np.int8))
        coverage, _, _, _ = qc.run(self.path, self.report, self.geo)
        self.assertEqual(len(coverage), 3)


if __name__ == "__main__":
    unittest.main()
