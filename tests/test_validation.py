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

    def test_genbank_flat_file_parsing(self):
        gb_text = """LOCUS       NM_000518               628 bp    mRNA    linear   PRI 01-FEB-2024
DEFINITION  Homo sapiens hemoglobin subunit beta (HBB), mRNA.
ACCESSION   NM_000518
ORGANISM    Homo sapiens
ORIGIN      
        1 acatttgctt ctgacacaac tgtgttcact agcaacctca aacagacacc atggtgcacc
       61 tgactcctga ggagaagtct gccgttactg ccctgtgggg caaggtgaac gtggatgaag
//"""
        res = validate_and_analyze(gb_text)
        self.assertTrue(res["valid"])
        self.assertEqual(res["type"], "DNA")
        self.assertEqual(len(res["sequence"]), 120)
        hints = res["hints"]
        self.assertEqual(hints["accession"], "NM_000518")
        self.assertEqual(hints["gene"], "HBB")
        self.assertEqual(hints["organism"], "Homo sapiens")

    def test_uniprot_flat_file_parsing(self):
        up_text = """ID   HBB_HUMAN               Reviewed;         147 AA.
AC   P68871;
DE   RecName: Full=Hemoglobin subunit beta;
GN   Name=HBB;
OS   Homo sapiens (Human).
SQ   SEQUENCE   147 AA;  15998 MW;  E479DE988CBEE2DA CRC64;
     MVHLTPEEKS AVTALWGKVN VDEVGGEALG RLLVVYPWTQ RFFESFGDLS TPDAVMGNPK
     VKAHGKKVLG AFSDGLAHLD NLKGTFATLS ELHCDKLHVD PENFRLLGNV LVCVLAHHFG
//"""
        res = validate_and_analyze(up_text)
        self.assertTrue(res["valid"])
        self.assertEqual(res["type"], "Protein")
        self.assertEqual(len(res["sequence"]), 120)
        hints = res["hints"]
        self.assertEqual(hints["accession"], "P68871")
        self.assertEqual(hints["gene"], "HBB")

    def test_three_letter_amino_acid_parsing(self):
        peptide_3aa = "Met-Val-His-Leu-Thr-Pro-Glu-Glu-Lys-Ser-Ala-Val-Thr-Ala-Leu"
        res = validate_and_analyze(peptide_3aa)
        self.assertTrue(res["valid"])
        self.assertEqual(res["type"], "Protein")
        self.assertEqual(res["sequence"], "MVHLTPEEKSAVTAL")

    def test_gapped_alignment_cleaning(self):
        gapped = ">alignment\nMVHLTPEEKS---AVTALWGKVNV.DEVGGEALGRLLVV~YPWTQRFFES"
        res = validate_and_analyze(gapped)
        self.assertTrue(res["valid"])
        self.assertEqual(res["type"], "Protein")
        self.assertNotIn("-", res["sequence"])
        self.assertNotIn(".", res["sequence"])
        self.assertNotIn("~", res["sequence"])

    def test_pdb_header_hints(self):
        pdb_fasta = ">4HHB_1|Chains A, C|Hemoglobin subunit alpha|Homo sapiens (9606)\nVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
        res = validate_and_analyze(pdb_fasta)
        self.assertTrue(res["valid"])
        hints = res["hints"]
        self.assertEqual(hints["accession"], "4HHB")
        self.assertEqual(hints["gene"], "HBA1")
        self.assertEqual(hints["organism"], "Homo sapiens")
        self.assertNotEqual(hints["gene"], "9606")

    def test_ensembl_header_hints(self):
        ens_fasta = ">ENST00000335295.4 cdna chromosome:GRCh38:11:5225464:5227071:1 gene:ENSG00000244734.4 gene_symbol:HBB\nATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGG"
        res = validate_and_analyze(ens_fasta)
        self.assertTrue(res["valid"])
        hints = res["hints"]
        self.assertEqual(hints["accession"], "ENST00000335295")
        self.assertEqual(hints["gene"], "HBB")

if __name__ == "__main__":
    unittest.main()
