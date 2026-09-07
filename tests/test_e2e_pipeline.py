"""
Automated end-to-end test verifying complete data flow from sequence validation
to BLAST identification, NCBI/UniProt annotation, literature categorization,
SQLite caching, and report compilation.
"""
import os
import unittest
from modules.validation import validate_and_analyze
from modules.blast_search import execute_blast
from modules.ncbi_annotation import fetch_ncbi_annotation
from modules.uniprot_annotation import fetch_uniprot_annotation
from modules.pubmed_retrieval import fetch_pubmed_literature
from modules.cache_manager import set_cached_result, get_cached_result
from modules.report_generator import render_html_report

class TestEndToEndPipeline(unittest.TestCase):
    def test_full_pipeline_flow(self):
        # 1. Validation
        seq = "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCC"
        val = validate_and_analyze(seq)
        self.assertTrue(val["valid"])
        self.assertEqual(val["type"], "DNA")
        self.assertGreater(val["stats"]["gc_percent"], 50)

        # 2. BLAST
        blast = execute_blast(val["sequence"], val["type"], use_mock_fallback=True)
        self.assertTrue(blast["success"])
        self.assertGreater(len(blast["hits"]), 0)
        top_acc = blast["hits"][0]["accession"]
        self.assertIsNotNone(top_acc)

        # 3. NCBI Annotation
        ncbi = fetch_ncbi_annotation(top_acc, val["type"])
        self.assertIn("organism", ncbi)

        # 4. UniProt Annotation
        uniprot = fetch_uniprot_annotation(ncbi.get("gene") or "HBB")
        self.assertTrue(uniprot.get("found", False))

        # 5. PubMed Literature
        lit = fetch_pubmed_literature(ncbi.get("gene") or "HBB", ncbi.get("organism") or "Homo sapiens")
        self.assertGreater(len(lit), 0)
        self.assertIn("category", lit[0])

        # 6. Cache
        full_res = {"validation": val, "blast": blast, "ncbi": ncbi, "uniprot": uniprot, "literature": lit}
        set_cached_result(val["sequence"], val["type"], full_res)
        cached = get_cached_result(val["sequence"])
        self.assertIsNotNone(cached)
        self.assertEqual(cached["validation"]["type"], "DNA")

        # 7. Report rendering
        tpl_path = os.path.join(os.path.dirname(__file__), "..", "templates", "report_template.html")
        with open(tpl_path, "r") as f:
            tpl_str = f.read()
        html = render_html_report(full_res, tpl_str)
        self.assertIn("SeqLit Bioinformatics Report", html)

if __name__ == "__main__":
    unittest.main()
