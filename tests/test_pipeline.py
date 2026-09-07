import unittest
from modules.cache_manager import get_cached_result, set_cached_result, compute_hash
from modules.pubmed_retrieval import categorize_article
from modules.blast_search import execute_blast

class TestPipeline(unittest.TestCase):
    def test_cache_storage_and_retrieval(self):
        seq = "ATGCGTACGTTAGC"
        sample_data = {"test_key": "test_val", "score": 99.5}
        set_cached_result(seq, "DNA", sample_data)
        
        retrieved = get_cached_result(seq)
        self.assertIsNotNone(retrieved)
        self.assertEqual(retrieved["test_key"], "test_val")
        self.assertEqual(retrieved["score"], 99.5)

    def test_pubmed_categorization(self):
        title_struct = "Crystal structure of the catalytic domain bound to substrate"
        self.assertEqual(categorize_article(title_struct, ""), "Structure / Crystallography")

        title_drug = "Discovery of potent small-molecule inhibitors for therapeutic application"
        self.assertEqual(categorize_article(title_drug, ""), "Drug / Pharmacology")

        title_disease = "Novel missense mutations associated with hereditary cardiovascular disease"
        self.assertEqual(categorize_article(title_disease, ""), "Disease / Pathology")

        title_complex = "Direct binding interaction of the regulatory subunit in yeast"
        self.assertEqual(categorize_article(title_complex, ""), "Interaction / Complex")

    def test_blast_mock_fallback(self):
        dna_seq = "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCC"
        res = execute_blast(dna_seq, "DNA", use_mock_fallback=True)
        self.assertTrue(res["success"])
        self.assertGreater(len(res["hits"]), 0)
        self.assertGreater(res["hits"][0]["identity_percent"], 0)

if __name__ == "__main__":
    unittest.main()
