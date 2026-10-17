import unittest
from modules.cache_manager import get_cached_result, set_cached_result, get_cached_hit, set_cached_hit, compute_hash
from modules.pubmed_retrieval import categorize_article, fetch_pubmed_literature
from modules.blast_search import extract_hit_metadata, score_curation_quality, local_benchmark_alignment
from modules.uniprot_annotation import fetch_uniprot_annotation

class TestPipeline(unittest.TestCase):
    def test_cache_storage_and_retrieval(self):
        seq = "ATGCGTACGTTAGC"
        sample_data = {"test_key": "test_val", "score": 99.5}
        set_cached_result(seq, "DNA", sample_data)
        
        retrieved = get_cached_result(seq)
        self.assertIsNotNone(retrieved)
        self.assertEqual(retrieved["test_key"], "test_val")
        self.assertEqual(retrieved["score"], 99.5)

    def test_hit_cache_storage_and_retrieval(self):
        acc = "NM_TEST_12345"
        sample_hit_data = {"ncbi": {"gene": "TESTG"}, "uniprot": {"protein_name": "Test Protein"}}
        set_cached_hit(acc, sample_hit_data)

        retrieved = get_cached_hit(acc)
        self.assertIsNotNone(retrieved)
        self.assertEqual(retrieved["ncbi"]["gene"], "TESTG")
        self.assertEqual(retrieved["uniprot"]["protein_name"], "Test Protein")

    def test_pubmed_categorization(self):
        title_struct = "Crystal structure of the catalytic domain bound to substrate"
        self.assertEqual(categorize_article(title_struct, ""), "Structure / Crystallography")

        title_drug = "Discovery of potent small-molecule inhibitors for therapeutic application"
        self.assertEqual(categorize_article(title_drug, ""), "Drug / Pharmacology")

        title_disease = "Novel missense mutations associated with hereditary cardiovascular disease"
        self.assertEqual(categorize_article(title_disease, ""), "Disease / Pathology")

        title_complex = "Direct binding interaction of the regulatory subunit in yeast"
        self.assertEqual(categorize_article(title_complex, ""), "Interaction / Complex")

    def test_extract_hit_metadata(self):
        raw_title = "PREDICTED: Pan paniscus hemoglobin subunit beta (LOC100976465), mRNA >sp|P68871.1| [Pan paniscus]"
        meta = extract_hit_metadata(raw_title, "XM_003819029", "PREDICTED: Pan paniscus hemoglobin subunit beta (LOC100976465), mRNA")
        self.assertEqual(meta["accession"], "XM_003819029")
        self.assertEqual(meta["organism"], "Pan paniscus")
        self.assertEqual(meta["gene_symbol"], "LOC100976465")
        self.assertEqual(meta["uniprot_id"], "P68871")

    def test_score_curation_quality(self):
        refseq_curated = {"accession": "NM_000518", "title": "Homo sapiens hemoglobin subunit beta (HBB), mRNA", "gene_symbol": "HBB"}
        automated_pred = {"accession": "XM_003819029", "title": "PREDICTED: Pan paniscus hemoglobin subunit beta (LOC100976465), mRNA", "gene_symbol": "LOC100976465"}
        synthetic_vector = {"accession": "AY894014", "title": "Synthetic construct Homo sapiens clone FLH130860.01L hemoglobin beta (HBB)", "gene_symbol": "HBB"}

        score_curated = score_curation_quality(refseq_curated)
        score_pred = score_curation_quality(automated_pred)
        score_synth = score_curation_quality(synthetic_vector)

        self.assertGreater(score_curated, score_pred, "Curated RefSeq should rank higher than automated predicted LOC contig")
        self.assertGreater(score_curated, score_synth, "Curated RefSeq should rank higher than synthetic cloning vector")

    def test_local_benchmark_alignment_hbb(self):
        hbb_seq = "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAA"
        hits = local_benchmark_alignment(hbb_seq, "DNA")
        self.assertGreaterEqual(len(hits), 1)
        self.assertEqual(hits[0]["accession"], "NM_000518")
        self.assertEqual(hits[0]["identity_percent"], 100.0)
        self.assertEqual(hits[0]["gene_symbol"], "HBB")

    def test_uniprot_direct_accession_resolution(self):
        # Testing EBNA-1 direct Swiss-Prot accession lookup returns correct protein (not EBNA-LP)
        res = fetch_uniprot_annotation("EBNA1", uniprot_id="P03211")
        self.assertTrue(res.get("found", False))
        self.assertIn("Epstein-Barr nuclear antigen 1", res.get("protein_name", ""))
        self.assertEqual(res.get("accession"), "P03211")

    def test_local_benchmark_alignment_protein_gfp(self):
        gfp_seq = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
        hits = local_benchmark_alignment(gfp_seq, "Protein")
        self.assertGreaterEqual(len(hits), 1)
        self.assertEqual(hits[0]["accession"], "P42212")
        self.assertEqual(hits[0]["gene_symbol"], "GFP")

    def test_local_benchmark_alignment_rna_16s(self):
        # 16S rRNA snippet containing uracil
        rna_sub = "AAAUUGAGAGUUUGAUCCUGGCUCAGGAUGAACGCUGGCGGCGUGCCUAAUACAUGCAAGUCGAACGGUAACAGGAAGAAGCUUGCUUCUUUGCUGACGAGUGGCGGACGGGUGAGUAAUGUCUGGGAAACUGCCUGAUGGAGGGGGAUAACUACUGGAAACGGUGGCUAAUACCGCAUAACGUCGCAAGACCAAAGAGGGGGACCUUCGGGCCUCUUGCCAUCGGAUGUGCCCAGAUGGGAUUAGCUAGUAGGUGGGGUAACGGCUCACCUAGGCGACGAUCCCUAGCUGGUCUGAGAGGAUGACCAGCCACACUGGAACUGAGACACGGUCCAGACUCCUACGGGAGGCAGCAGUGGGGAAUAUUGCACAAUGGGGGAAACCCUGAUGCAGCCAUGCCGCGUGUAUGAAGAAGGCCUUCGGGUUGUAAAGUACUUUCAGCGGGGAGGAAGGGAGUAAAGUUAAUACCUUUGCUCAUUGACGUUACCCGCAGAAGAAGCACCGGCUAACUCCGUGCCAGCAGCCGCGGUAAUACGGAGGGUGCAAGCGUUAAUCGGAAUUACUGGGCGUAAAGCGCACGCAGGCGGUUUGUUAAGUCAGAUGUGAAAUCCCCGGGCUCAACCUGGGAACUGCAUCUGAUACUGGCAAGCUUGAGUCUCGUAGAGGGGGGUAGAAUUCCAGGUGUAGCGGUGAAAUGCGUAGAGAUCUGGAGGAAUACCGGUGGCGAAGGCGGCCCCCUGGACGAAGACUGACGCUCAGGUGCGAAAGCGUGGGGAGCAAACAGGAUUAGAUACCCUGGUAGUCCACGCCGUAAACGAUGUCGACUUGGAGGUUGUGCCCUUGAGGCGUGGCUUCCGGAGCUAACGCGUUAAGUCGACCGCCUGGGGAGUACGGCCGCAAGGUUAAAACUCAAAUGAAUUGACGGGGGCCCGCACAAGCGGUGGAGCAUGUGGUUUAAUUCGAUGCAACGCGAAGAACCUUACCUGGUCUUGACAUCCACGGAAGUUUUCAGAGAUGAGAAUGUGCCUUCGGGAACCGUGAGACAGGUGCUGCAUGGCUGUCGUCAGCUCGUGUUGUGAAAUGUUGGGUUAAGUCCCGCAACGAGCGCAACCCUUAUCCUUUGUUGCCAGCGGUCCGGCCGGGAACUCAAAGGAGACUGCCAGUGAUAAACUGGAGGAAGGUGGGGAUGACGUCAAGUCAUCAUGGCCUUUACGACCAGGGCUACACACGUGCUACAAUGGCGCAUACAAAGAGAAGCGACCUCGCGAGAGCAAGCGGACCUCAUAAAGUGCGUCGUAGUCCGGAUUGGAGUCUGCAACUCGACUCCAUGAAGUCGGAAUCGCUAGUAAUCGUGGAUCAGAAUGCCACGGUGAAUACGUUCCCGGGCCUUGUACACACCGCCCGUCACACCAUGGGAGUGGGUUGCAAAAGAAGUAGGUAGCUUAACCUUCGGGAGGGCGCUUACCACUUUGUGAUUCAUGACUGGGGUGAAGUCGUAACAAGGUAACCGUAGGGGAACCUGCGGUUGGAUCACCUCCUUA"
        hits = local_benchmark_alignment(rna_sub, "RNA")
        self.assertGreaterEqual(len(hits), 1)
        self.assertEqual(hits[0]["accession"], "NR_024570")

    def test_local_benchmark_alignment_fragment(self):
        # 50 bp sub-sequence fragment of HBB
        hbb_full = "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAA"
        frag = hbb_full[40:90]
        hits = local_benchmark_alignment(frag, "DNA")
        self.assertGreaterEqual(len(hits), 1)
        self.assertEqual(hits[0]["gene_symbol"], "HBB")
        self.assertEqual(hits[0]["identity_percent"], 100.0)
        self.assertEqual(hits[0]["query_coverage"], 100.0)

    def test_local_benchmark_alignment_reverse_complement(self):
        # Reverse complement fragment of HBB
        hbb_full = "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAA"
        frag = hbb_full[50:100]
        comp_tab = str.maketrans("ATCG", "TAGC")
        rc_frag = frag.translate(comp_tab)[::-1]
        hits = local_benchmark_alignment(rc_frag, "DNA")
        self.assertGreaterEqual(len(hits), 1)
        self.assertEqual(hits[0]["gene_symbol"], "HBB")
        self.assertEqual(hits[0]["identity_percent"], 100.0)
        self.assertTrue(hits[0].get("is_rc", False))

    def test_local_benchmark_alignment_gapdh(self):
        gapdh_seq = "ATGGGGAAGGTGAAGGTCGGAGTCAACGGATTTGGTCGTATTGGGCGCCTGGTCACCAGGGCTGCTTTTAACTCTGGTAAAGTGGATATTGTTGCCATCAATGACCCCTTCATTGACCTCAACTACATGGTTTACATGTTCCAATATGATTCCACCCATGGCAAATTCCATGGCACCGTCAAGGCTGAGAACGGGAAGCTTGTCATCAATGGAAATCCCATCACCATCTTCCAGGAGCGAGATCCCTCCAAAATCAAGTGGGGCGATGCTGGCGCTGAGTACGTCGTGGAGTCCACTGGCGTCTTCACCACCATGGAGAAGGCTGGGGCTCATTTGCAGGGGGGAGCCAAAAGGGTCATCATCTCTGCCCCCTCTGCTGATGCCCCCATGTTCGTCATGGGTGTGAACCATGAGAAGTATGACAACAGCCTCAAGATCATCAGCAATGCCTCCTGCACCACCAACTGCTTAGCACCCCTGGCCAAGGTCATCCATGACAACTTTGGTATCGTGGAAGGACTCATGACCACAGTCCATGCCATCACTGCCACCCAGAAGACTGTGGATGGCCCCTCCGGGAAACTGTGGCGTGATGGCCGCGGGGCTCTCCAGAACATCATCCCTGCCTCTACTGGCGCTGCCAAGGCTGTGGGCAAGGTCATCCCTGAGCTGAACGGGAAGCTCACTGGCATGGCCTTCCGTGTCCCCACTGCCAACGTGTCAGTGGTGGACCTGACCTGCCGTCTAGAAAAACCTGCCAAATATGATGACATCAAGAAGGTGGTGAAGCAGGCGTCGGAGGGCCCCCTCAAGGGCATCCTGGGCTACACTGAGCACCAGGTGGTCTCCTCTGACTTCAACAGCGACACCCACTCCTCCACCTTTGACGCTGGGGCTGGCATTGCCCTCAACGACCACTTTGTCAAGCTCATTTCCTGGTATGACAACGAATTTGGCTACAGCAACAGGGTGGTGGACCTCATGGCCCACATGGCCTCCAAGGAGTAA"
        hits = local_benchmark_alignment(gapdh_seq, "DNA")
        self.assertGreaterEqual(len(hits), 1)
        self.assertEqual(hits[0]["gene_symbol"], "GAPDH")
        self.assertEqual(hits[0]["accession"], "NM_002046")

    def test_local_benchmark_alignment_actb(self):
        actb_seq = "MDDDIAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTEAPLNPKANREKMTQIMFETFNTPAMYVAIQAVLSLYASGRTTGIVMDSGDGVTHTVPIYEGYALPHAILRLDLAGRDLTDYLMKILTERGYSFTTTAEREIVRDIKEKLCYVALDFEQEMATAASSSSLEKSYELPDGQVITIGNERFRCPEALFQPSFLGMESCGIHETTFNSIMKCDVDIRKDLYANTVLSGGTTMYPGIADRMQKEITALAPSTMKIKIIAPPERKYSVWIGGSILASLSTFQQMWISKQEYDESGPSIVHRKCF"
        hits = local_benchmark_alignment(actb_seq, "Protein")
        self.assertGreaterEqual(len(hits), 1)
        self.assertEqual(hits[0]["gene_symbol"], "ACTB")
        self.assertEqual(hits[0]["accession"], "P60709")

if __name__ == "__main__":
    unittest.main()
