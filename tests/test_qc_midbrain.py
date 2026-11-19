# %% Phase 2e tests: synthetic 10x CSC HDF5
"""Run: pixi run python -m unittest discover -s tests -p 'test_qc_midbrain.py' -v"""
import importlib.util
from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("midbrain_qc", ROOT / "scripts" / "02e_qc_midbrain.py")
qc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(qc)


class MidbrainQCTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.path = Path(temp.name) / "library.h5"
        with h5py.File(self.path, "w") as h5:
            m = h5.create_group("matrix")
            # 3 features x 2 barcodes: [MT-A:5, G1:7] ; [G1:2, Antibody:100]
            m.create_dataset("shape", data=[3, 2])
            m.create_dataset("barcodes", data=[b"AA-1", b"BB-1"])
            m.create_dataset("data", data=np.array([5, 7, 2, 100], dtype=np.int32))
            m.create_dataset("indices", data=np.array([0, 1, 1, 2], dtype=np.int32))
            m.create_dataset("indptr", data=np.array([0, 2, 4], dtype=np.int32))
            ft = m.create_group("features")
            ft.create_dataset("id", data=[b"ENSGMT", b"ENSG1", b"AB1"])
            ft.create_dataset("name", data=[b"MT-A", b"G1", b"Tag"])
            ft.create_dataset("feature_type", data=[b"Gene Expression", b"Gene Expression", b"Antibody Capture"])
            ft.create_dataset("genome", data=[b"GRCh38", b"GRCh38", b""])

    def test_counts_are_gene_expression_only_and_read_only(self):
        before = self.path.read_bytes()
        r = qc.summarize_10x_h5(self.path, block=1)
        self.assertEqual(self.path.read_bytes(), before)
        self.assertEqual(r["filtered_barcodes"], 2)
        self.assertEqual(r["gene_expression_features"], 2)
        self.assertEqual(r["gene_count_sum"], 14)
        self.assertEqual(r["stored_entries"], 4)
        self.assertEqual(r["umi_quantiles"]["median"], 7.0)
        self.assertEqual(r["detected_gene_quantiles"]["median"], 1.5)
        self.assertEqual(r["putative_mt_gene_features"], 1)
        self.assertEqual(r["mt_fraction_status"], "NAME_PREFIX_ONLY_NOT_SOURCE_VERIFIED")
        self.assertEqual(r["donor_assignment"], "NOT_AVAILABLE")

    def test_refuse_invalid_umi_values(self):
        with h5py.File(self.path, "r+") as h5:
            h5["matrix/data"][0] = -1
        with self.assertRaisesRegex(ValueError, "invalid UMI"):
            qc.summarize_10x_h5(self.path)

    def test_refuse_incorrect_barcode_boundaries(self):
        with h5py.File(self.path, "r+") as h5:
            h5["matrix/indptr"][-1] = 3
        with self.assertRaisesRegex(ValueError, "boundaries invalid"):
            qc.summarize_10x_h5(self.path)


if __name__ == "__main__":
    unittest.main()
