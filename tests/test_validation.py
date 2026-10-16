import unittest
from modules.validation import clean_sequence, detect_sequence_type, calculate_stats, validate_and_analyze

class TestValidation(unittest.TestCase):
    def test_clean_sequence_fasta(self):
        fasta_input = """>Test_seq_header
        ATGC
        ATGC 123
        """
        seq, header = clean_sequence(fasta_input)
        self.assertEqual(seq, "ATGCATGC")
        self.assertEqual(header, "Test_seq_header")

    def test_detect_sequence_type_dna(self):
        dna = "ATGCGATCGATCGATC"
        self.assertEqual(detect_sequence_type(dna), "DNA")

    def test_detect_sequence_type_rna(self):
        rna = "AUGCGAUCGAUCGAUC"
        self.assertEqual(detect_sequence_type(rna), "RNA")

    def test_detect_sequence_type_protein(self):
        protein = "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVL"
        self.assertEqual(detect_sequence_type(protein), "Protein")

    def test_detect_sequence_type_invalid(self):
        invalid = "ATGCTZZ123$"
        cleaned, _ = clean_sequence(invalid)
        self.assertEqual(detect_sequence_type(cleaned), "Invalid")

    def test_calculate_stats_dna(self):
        dna = "GGCC"
        stats = calculate_stats(dna, "DNA")
        self.assertEqual(stats["length"], 4)
        self.assertEqual(stats["gc_percent"], 100.0)
        self.assertEqual(stats["composition"]["G"], 2)
        self.assertEqual(stats["composition"]["C"], 2)

    def test_validate_and_analyze_empty(self):
        result = validate_and_analyze("")
        self.assertFalse(result["valid"])
        self.assertIn("empty", result["error"].lower())

    def test_validate_and_analyze_valid(self):
        dna = "ATGC"
        result = validate_and_analyze(dna)
        self.assertTrue(result["valid"])
        self.assertEqual(result["type"], "DNA")
        self.assertEqual(result["stats"]["length"], 4)

    def test_degenerate_dna(self):
        degen = "ATGCYRWSKMDVHBXN"
        result = validate_and_analyze(degen)
        self.assertTrue(result["valid"])
        self.assertEqual(result["type"], "DNA")

    def test_degenerate_rna(self):
        rna_degen = "AUGCRYWSKMDVHBXN"
        result = validate_and_analyze(rna_degen)
        self.assertTrue(result["valid"])
        self.assertEqual(result["type"], "RNA")

    def test_protein_trailing_stop_codon(self):
        prot_stop = "MVHLTPEEKSAVTALW*"
        cleaned, _ = clean_sequence(prot_stop)
        self.assertEqual(cleaned, "MVHLTPEEKSAVTALW")
        result = validate_and_analyze(prot_stop)
        self.assertTrue(result["valid"])
        self.assertEqual(result["type"], "Protein")

    def test_fasta_header_hints(self):
        fasta = ">sp|P03211|EBNA1_EBVB9 Epstein-Barr nuclear antigen 1 OS=Epstein-Barr virus OX=10377 GN=EBNA1 PE=1 SV=1\nMSDEGPGTGPGN"
        result = validate_and_analyze(fasta)
        self.assertTrue(result["valid"])
        hints = result.get("hints", {})
        self.assertEqual(hints.get("accession"), "P03211")
        self.assertEqual(hints.get("gene"), "EBNA1")
        self.assertIn("Epstein-Barr virus", hints.get("organism", ""))

if __name__ == "__main__":
    unittest.main()
