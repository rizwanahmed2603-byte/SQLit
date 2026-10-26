"""
Unit and integration tests for multi-species orthologs, spliced genomic sequences,
homopolymers, and synthetic variants.
"""

import unittest
from modules.blast_search import execute_blast
from modules.validation import detect_sequence_type

class TestOrthologsAndVariants(unittest.TestCase):
    def test_canine_insulin_ortholog(self):
        # Canine preproinsulin (110 aa)
        canine_ins = "MALWMRLLAVLALLALWRPDTAQAFVNRHLCGSHLVEALYLVCGDNGFFYTPKERRELENLQVGQAEPGMGLEAGGLQPLAQELALQKRGIVDQCCNSICSLYQLQNYCN"
        res = execute_blast(canine_ins, "Protein")
        self.assertTrue(res.get("success"))
        hits = res.get("hits", [])
        self.assertGreater(len(hits), 0)
        top = hits[0]
        self.assertIn("INS", top.get("gene_symbol", "") or top.get("title", ""))
        self.assertGreaterEqual(top.get("identity_percent", 0), 65.0)

    def test_single_chain_insulin_variant(self):
        # Single-chain insulin DNA (153 bp)
        sci_dna = "CTCGTCAGCCGGCGCCTGTGCGGCTCTAACTTAGTGGATGCGCTGTTCTTGCTATGTGAACGGAATGGCTTCTACAGACCCAACGATGGTATCGTGGACCAATGCTGCATCAATACCTGCACAACGAACCAGCTGCTGAACTACTGCAATTAG"
        res = execute_blast(sci_dna, "DNA")
        self.assertTrue(res.get("success"))
        hits = res.get("hits", [])
        self.assertGreater(len(hits), 0)
        top = hits[0]
        self.assertEqual(top.get("identity_percent"), 100.0)

    def test_ebna1_cterm_domain(self):
        # EBNA-1 C-terminal DNA-binding domain (178 aa)
        ebna1_cterm = "PADDPGEGPSTGPRGQGDGGRRKKGGWFGKHRGQGGSNPKFENIAEGLRALLARSHVERTTDEGTWVAGVFVYGGSKTSLYNLRRGTALAIPQCRLTPLSRLPFGMAPGPGPQPGPLRESIVCYFMVFLQTHIFAEVLKDAIKDLVMTKPAPTCNIRVTVCSFDDGVDLPPWFPPMVE"
        res = execute_blast(ebna1_cterm, "Protein")
        self.assertTrue(res.get("success"))
        hits = res.get("hits", [])
        self.assertGreater(len(hits), 0)
        top = hits[0]
        self.assertIn("EBNA", top.get("title", ""))
        self.assertEqual(top.get("identity_percent"), 100.0)

    def test_degus_insulin_raw_dna(self):
        # Raw degus insulin DNA without FASTA header (430 bp)
        degus_dna = "GCATTCTGAGGCATTCTCTAACAGGTTCTCGACCCTCCGCCATGGCCCCGTGGATGCATCTCCTCACCGTGCTGGCCCTGCTGGCCCTCTGGGGACCCAACTCTGTTCAGGCCTATTCCAGCCAGCACCTGTGCGGCTCCAACCTAGTGGAGGCACTGTACATGACATGTGGACGGAGTGGCTTCTATAGACCCCACGACCGCCGAGAGCTGGAGGACCTCCAGGTGGAGCAGGCAGAACTGGGTCTGGAGGCAGGCGGCCTGCAGCCTTCGGCCCTGGAGATGATTCTGCAGAAGCGCGGCATTGTGGATCAGTGCTGTAATAACATTTGCACATTTAACCAGCTGCAGAACTACTGCAATGTCCCTTAGACACCTGCCTTGGGCCTGGCCTGCTGCTCTGCCCTGGCAACCAATAAACCCCTTGAATG"
        res = execute_blast(degus_dna, "DNA")
        self.assertTrue(res.get("success"))
        hits = res.get("hits", [])
        self.assertGreater(len(hits), 0)
        top = hits[0]
        self.assertIn("INS", top.get("gene_symbol", "") or top.get("title", ""))

    def test_guinea_pig_genomic_gene(self):
        # Spliced preproinsulin genomic sequence with introns (1472 bp)
        genomic_dna = "CTGCAGACCCAGCACCAGGGAAATGATCCAGAAATTGCAACCTCAGCCCCCTGGCCATCTGCTGATGCCACCACCCCCAGGTCCCTAATGGGCCTGGTGGCAGAGTTTGGGAAGATGGGCTCAGGGCTATATAAAGTCCACAAGGACCTAAGAGCCCCCAGTGCTGCTGGGCCAGCTGTATTCTGAGGTGGTCAGCACACAGGTCTGTGTCCTCCGTGCTAGATTGGGGCTGAGAGGCTGGGGGCTCTGGGTTGGCTGGGACAGGACATGGGATTCTTCCTTGTATTGGGGGTTTTGGCTGTTACTCTGTCTCTCCATCAGGTCATCATCCTTTCATCATGGCTCTGTGGATGCATCTCCTCACCGTGCTGGCCCTGCTGGCCCTCTGGGGGCCCAACACTAATCAGGCCTTTGTCAGCCGGCATCTGTGCGGCTCCAACTTAGTGGAGACATTGTATTCAGTGTGTCAGGATGATGGCTTCTTCTATATACCCAAGGACCGTCGGGAGCTAGAGGACCCACAGGGTGAGCCCCTACCTGCCATCCCTGCTGTTTCCGTGCCAGTACCCCAGCTGGCAGGGCATAAGTAAGCAGGAAGCTAATTCCAAGGAGAGTCGATGGGTTTGTTGAAAAGGGAGGCGGCTCTCTTGGTCATTTCGTAAAGTGGTGGTGGCTTCCTATAGCTGCTTTTAAGGGTAAAGGGTAACAGCTGCACCCTTCAGCTGTGGCTTCTGAGCACAACTGGACTCTTCCCTCCACTTGCCTTCGAATGACTGCCCTGGCCTCATGGCAACAGTAGCTCCCTGGTACCAATTTTATTATGCAGATTGCATCTTGGTGTTGATAGCCTTAGGGTAGCCTGGGGGCCATTCATGGGGCGCCCCATCCCTCCTTCCTCCCTGCCTCTGGACAAATGCTCCATGGAGCTCCAAGCTCTGCCACGTGGGAGGTGTGGGTCTCCAGCGCTCTGTGTGCCCAGCATGGCAGCCTCTGTCACCTGGACCAGCTCCCTGGGAGATGCAGTGAGAGGGTGGTAGTGTGGGGCCAGTGCGCAGGCATTCTGCTGCTCCTGACAGCATCTGCCCCTGTCTCTCTCCCCACTGCTGCTGCTCCTGTATTCTGGCACCTCACCCTGCAGTGGAGCAGACAGAACTGGGCATGGGCCTGGGGGCAGGTGGACTACAGCCCTTGGCACTGGAGATGGCACTACAGAAGCGTGGCATTGTGGATCAGTGCTGTACTGGCACCTGCACACGCCACCAGCTGCAGAGCTACTGCAACTAG"
        res = execute_blast(genomic_dna, "DNA")
        self.assertTrue(res.get("success"))
        hits = res.get("hits", [])
        self.assertGreater(len(hits), 0)
        top = hits[0]
        self.assertIn("INS", top.get("gene_symbol", "") or top.get("title", ""))

    def test_poly_methionine_detection(self):
        poly_m = "M" * 50
        stype = detect_sequence_type(poly_m)
        self.assertEqual(stype, "Protein")
        res = execute_blast(poly_m, stype)
        self.assertTrue(res.get("success"))
        hits = res.get("hits", [])
        self.assertGreater(len(hits), 0)
        self.assertIn("Poly-L-methionine", hits[0].get("title", ""))

    def test_degenerate_dna_vs_homopolymer(self):
        self.assertEqual(detect_sequence_type("ATGCRYSW"), "DNA")
        self.assertEqual(detect_sequence_type("R" * 30), "Protein")
        self.assertEqual(detect_sequence_type("K" * 30), "Protein")

if __name__ == "__main__":
    unittest.main()
