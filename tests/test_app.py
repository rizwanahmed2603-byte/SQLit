"""
End-to-End Flask application and API route tests for SeqLit.
Verifies submission, validation, benchmark samples, hit annotation,
and report generation endpoints.
"""

import unittest
import json
from app import app

class TestSeqLitApp(unittest.TestCase):
    def setUp(self):
        self.app = app
        self.client = self.app.test_client()

    def test_index_route(self):
        response = self.client.get("/")
        self.assertEqual(response.status_code, 200)
        self.assertIn(b"SeqLit", response.data)

    def test_sample_route(self):
        response = self.client.get("/api/sample")
        self.assertEqual(response.status_code, 200)
        data = response.get_json()
        self.assertIn("sequence", data)
        self.assertIn("HBB", data.get("name", ""))

    def test_samples_list_route(self):
        response = self.client.get("/api/samples")
        self.assertEqual(response.status_code, 200)
        data = response.get_json()
        self.assertIsInstance(data, list)
        self.assertGreaterEqual(len(data), 3)

    def test_analyze_empty_sequence(self):
        response = self.client.post("/api/analyze", json={"sequence": ""})
        self.assertEqual(response.status_code, 400)
        data = response.get_json()
        self.assertFalse(data.get("success", True))

    def test_analyze_invalid_sequence(self):
        response = self.client.post("/api/analyze", json={"sequence": "ATGC$$$ZZZ"})
        self.assertEqual(response.status_code, 400)
        data = response.get_json()
        self.assertFalse(data.get("success", True))

    def test_analyze_hbb_benchmark(self):
        hbb_seq = "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAA"
        response = self.client.post("/api/analyze", json={"sequence": hbb_seq})
        self.assertEqual(response.status_code, 200)
        res = response.get_json()
        self.assertTrue(res.get("success"))
        data = res.get("data", {})

        # Validation stats
        self.assertEqual(data["validation"]["type"], "DNA")
        self.assertEqual(data["validation"]["stats"]["length"], 444)

        # BLAST Identification
        hits = data["blast"].get("hits", [])
        self.assertGreaterEqual(len(hits), 1)
        top_hit = hits[0]
        self.assertIn("HBB", top_hit["title"])

        # NCBI Annotation
        self.assertEqual(data["ncbi"]["gene"], "HBB")
        self.assertEqual(data["ncbi"]["organism"], "Homo sapiens")

        # UniProt Annotation
        self.assertEqual(data["uniprot"]["gene"], "HBB")
        self.assertIn("Hemoglobin subunit beta", data["uniprot"]["protein_name"])

        # PubMed Literature
        self.assertGreaterEqual(len(data["literature"]), 1)

    def test_annotate_hit_route(self):
        payload = {
            "accession": "NM_000518",
            "seq_type": "DNA",
            "title": "Homo sapiens hemoglobin subunit beta (HBB), mRNA",
            "gene_symbol": "HBB",
            "organism": "Homo sapiens",
            "uniprot_id": "P68871"
        }
        response = self.client.post("/api/annotate_hit", json=payload)
        self.assertEqual(response.status_code, 200)
        res = response.get_json()
        self.assertTrue(res.get("success"))
        hit_data = res.get("data", {})
        self.assertEqual(hit_data["ncbi"]["gene"], "HBB")
        self.assertEqual(hit_data["uniprot"]["accession"], "P68871")
        self.assertGreaterEqual(len(hit_data["literature"]), 1)

    def test_export_report_route(self):
        sample_analysis = {
            "validation": {
                "header": "NM_000518 HBB",
                "type": "DNA",
                "stats": {
                    "length": 444,
                    "gc_percent": 56.3,
                    "molecular_weight": 129566.12,
                    "frequencies": {"A": 20, "T": 24, "G": 30, "C": 26}
                }
            },
            "blast": {
                "hits": [
                    {
                        "rank": 1,
                        "accession": "NM_000518",
                        "title": "Homo sapiens hemoglobin subunit beta (HBB), mRNA",
                        "identity_percent": 100.0,
                        "e_value": 0.0,
                        "query_coverage": 100.0
                    }
                ]
            },
            "ncbi": {
                "organism": "Homo sapiens",
                "gene": "HBB",
                "taxonomy": "Eukaryota; Primates; Homo"
            },
            "uniprot": {
                "protein_name": "Hemoglobin subunit beta",
                "function": "Involved in oxygen transport."
            },
            "literature": [
                {
                    "pmid": "31774950",
                    "title": "CRISPR-Cas9 Gene Editing for Sickle Cell Disease and beta-Thalassemia",
                    "journal": "N Engl J Med",
                    "year": "2021",
                    "category": "Disease / Pathology",
                    "pubmed_url": "https://pubmed.ncbi.nlm.nih.gov/31774950/"
                }
            ]
        }
        response = self.client.post("/api/export/report", json={"data": sample_analysis})
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.mimetype, "text/html")
        self.assertIn(b"SeqLit Bioinformatics Report", response.data)

if __name__ == "__main__":
    unittest.main()
