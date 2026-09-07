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

if __name__ == "__main__":
    unittest.main()
