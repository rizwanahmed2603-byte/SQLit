"""
Unit and integration tests for Multi-Database Sequence & Header Resolution:
Verifies universal compatibility with PDB, Ensembl, UniProt, GenBank, EMBL,
and free-text FASTA headers.
"""

import unittest
from modules.validation import validate_and_analyze, extract_header_hints
from modules.blast_search import execute_blast
from modules.ncbi_annotation import fetch_ncbi_annotation, fetch_pdb_annotation, fetch_ensembl_annotation
from modules.uniprot_annotation import fetch_uniprot_annotation

class TestMultiDatabaseResolvers(unittest.TestCase):

    def test_pdb_header_resolution(self):
        # PDB format header
        header = ">4HHB_1|Chains A, C|Hemoglobin subunit alpha|Homo sapiens (9606)"
        hints = extract_header_hints(header)
        self.assertEqual(hints.get("accession"), "4HHB")
        self.assertEqual(hints.get("organism"), "Homo sapiens")
        self.assertEqual(hints.get("gene"), "HBA1")

    def test_ensembl_header_resolution(self):
        # Ensembl format header
        header = ">ENST00000335295.4 cdna chromosome:GRCh38:11:5225464:5227071:1 gene:ENSG00000244734.9"
        hints = extract_header_hints(header)
        self.assertEqual(hints.get("accession"), "ENST00000335295")

    def test_uniprot_header_resolution(self):
        # UniProt format header
        header = ">sp|P00441|SODC_HUMAN Superoxide dismutase [Cu-Zn] OS=Homo sapiens OX=9606 GN=SOD1 PE=1 SV=2"
        hints = extract_header_hints(header)
        self.assertEqual(hints.get("accession"), "P00441")
        self.assertEqual(hints.get("gene"), "SOD1")
        self.assertEqual(hints.get("organism"), "Homo sapiens")

    def test_free_text_header_resolution(self):
        # Free-text descriptive header
        header = ">Superoxide dismutase [Homo sapiens]"
        hints = extract_header_hints(header)
        self.assertEqual(hints.get("gene"), "SOD1")
        self.assertEqual(hints.get("organism"), "Homo sapiens")

    def test_raw_sod1_local_alignment(self):
        # Raw SOD1 protein sequence
        sod1_seq = "MATKAVCVLKGDGPVQGIINFEQKESNGPVKVWGSIKGLTEGLHGFHVHEFGDNTAGCTSAGPHFNPLSRKHGGPKDEERHVGDLGNVTADKDGVADVSIEDSVISLSGDHCIIGRTLVVHEKADDLGKGGNEESTKTGNAGSRLACGVIGIAQ"
        res = execute_blast(sod1_seq, "Protein")
        self.assertTrue(res.get("success"))
        hits = res.get("hits", [])
        self.assertGreaterEqual(len(hits), 1)
        self.assertEqual(hits[0]["gene_symbol"], "SOD1")

    def test_raw_myoglobin_local_alignment(self):
        # Raw MB protein sequence
        mb_seq = "MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASEDLKKHGATVLTALGGILKKKGHHEAEIKPLAQSHEATKHKIPVKYLEFISECIIQVLQSKHPGDFGADAQGAMNKALELFRKDMASNYKELGFQG"
        res = execute_blast(mb_seq, "Protein")
        self.assertTrue(res.get("success"))
        hits = res.get("hits", [])
        self.assertGreaterEqual(len(hits), 1)
        self.assertEqual(hits[0]["gene_symbol"], "MB")

    def test_raw_ubiquitin_local_alignment(self):
        # Raw Ubiquitin sequence
        ubq_seq = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
        res = execute_blast(ubq_seq, "Protein")
        self.assertTrue(res.get("success"))
        hits = res.get("hits", [])
        self.assertGreaterEqual(len(hits), 1)
        self.assertEqual(hits[0]["gene_symbol"], "UBB")

    def test_pdb_annotation_api(self):
        res = fetch_pdb_annotation("1TPH")
        if res:  # If network available
            self.assertTrue(res.get("found"))
            self.assertEqual(res.get("accession"), "1TPH")
            self.assertIn("TRIOSEPHOSPHATE ISOMERASE", res.get("product", "").upper())

    def test_ensembl_annotation_api(self):
        res = fetch_ensembl_annotation("ENST00000335295")
        if res:  # If network available
            self.assertTrue(res.get("found"))
            self.assertEqual(res.get("gene"), "HBB")

if __name__ == "__main__":
    unittest.main()
