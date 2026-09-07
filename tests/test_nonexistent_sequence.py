"""
Unit and integration tests verifying that non-existent or synthetic sequences
accurately return 0 hits and never fabricate fake biological entries.
"""
import unittest
from modules.validation import validate_and_analyze
from modules.blast_search import parse_blast_xml
from modules.ncbi_annotation import fetch_ncbi_annotation
from modules.uniprot_annotation import fetch_uniprot_annotation
from modules.pubmed_retrieval import fetch_pubmed_literature

class TestAuthenticResults(unittest.TestCase):
    def test_nonexistent_accession_ncbi(self):
        # Empty or non-existent accession should return found=False
        res = fetch_ncbi_annotation("", "DNA")
        self.assertFalse(res.get("found", True))
        self.assertEqual(res.get("gene"), "N/A")

    def test_nonexistent_gene_uniprot(self):
        # Unrecognized random gene symbol should not return fake Insulin
        res = fetch_uniprot_annotation("XYZNONEXISTENTGENE999")
        self.assertFalse(res.get("found", True))
        self.assertNotEqual(res.get("protein_name"), "Insulin preproprotein")

    def test_empty_literature_for_unknown_sequence(self):
        # Should return empty list, not mock articles
        articles = fetch_pubmed_literature("XYZNONEXISTENTGENE999", "Unknown")
        self.assertEqual(len(articles), 0)

if __name__ == "__main__":
    unittest.main()
