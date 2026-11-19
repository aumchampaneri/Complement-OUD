# %% Source-verified annotation audit tests
"""Run: pixi run python -m unittest discover -s tests -p 'test_qc_bulk*.py' -v."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("qc_02c", ROOT / "scripts" / "02c_qc_bulk_annotation.py")
qc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(qc)


class AnnotationAuditTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.folder = Path(temp.name)

    def test_gtf_versionless_mapping_and_mito_reference(self):
        path = self.folder / "example.gtf"
        path.write_text(
            'chrM\tsrc\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG000001.7"; gene_name "MT-A";\n'
            'chr1\tsrc\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG000002.2"; gene_name "B";\n',
            encoding="utf-8")
        reference, ambiguity = qc.gtf_records(path)
        self.assertEqual(ambiguity, set())
        mapping, audit = qc.map_features(["ENSG000001.7", "ENSG000002.2", "ENSG999999"],
                                         reference, ambiguity, "GSE174409")
        self.assertTrue(qc.canonical_mito(mapping["ENSG000001.7"]["chromosome"]))
        self.assertEqual(audit[-1]["mapping_status"], "UNMAPPED")

    def test_reference_collision_conservatively_disallowed(self):
        path = self.folder / "dupe.gtf"
        path.write_text(
            'chrM\tsrc\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG001.1"; gene_name "A";\n'
            'chr1\tsrc\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG001.2"; gene_name "B";\n',
            encoding="utf-8")
        reference, ambiguity = qc.gtf_records(path)
        self.assertIn("ENSG001", ambiguity)
        self.assertNotIn("ENSG001", reference)

    def test_ncbi_taxonomy_and_exact_id(self):
        path = self.folder / "human.gene_info"
        path.write_text(
            "#tax_id\tGeneID\tSymbol\tchromosome\n"
            "9606\t4512\tMT-CO1\tMT\n"
            "10090\t4512\tWrong\t1\n"
            "9606\t10\tGENE\t1\n",
            encoding="utf-8")
        records, ambiguous = qc.ncbi_records(path)
        self.assertFalse(ambiguous)
        self.assertEqual(records["4512"]["gene_symbol"], "MT-CO1")
        self.assertEqual(len(records), 2)


if __name__ == "__main__":
    unittest.main()
