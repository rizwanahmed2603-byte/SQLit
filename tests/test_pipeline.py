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

    def test_local_benchmark_alignment_kras(self):
        kras_frag = "ATGACTGAATATAAACTTGTGGTAGTTGGAGCTGGTGGCGTAGGCAAGAGTGCCTTGACGATACAGCTAATTCAGAATCATTTTGTGGACGAATATGATCCAACAATAGAGGATTCCTACAGGAAG"
        hits = local_benchmark_alignment(kras_frag, "DNA")
        self.assertGreaterEqual(len(hits), 1)
        self.assertEqual(hits[0]["gene_symbol"], "KRAS")
        self.assertEqual(hits[0]["identity_percent"], 100.0)

    def test_local_benchmark_alignment_pten(self):
        pten_frag = "MTAIIKEIVSRNKRRYQEDGFDLDLTYIYPNIIAMGFPAERLEGVYRNNIDDVVRFLDSKHKNHYKIYNLCAERHYDTAKFNCRVAQYPFEDHNPPQLELIKPFCEDLDQWLSEDDNHVAAIHCK"
        hits = local_benchmark_alignment(pten_frag, "Protein")
        self.assertGreaterEqual(len(hits), 1)
        self.assertEqual(hits[0]["gene_symbol"], "PTEN")
        self.assertEqual(hits[0]["identity_percent"], 100.0)

    def test_local_benchmark_alignment_sars2(self):
        spike_rbd = "CPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNF"
        hits = local_benchmark_alignment(spike_rbd, "Protein")
        self.assertGreaterEqual(len(hits), 1)
        self.assertEqual(hits[0]["gene_symbol"], "S")
        self.assertEqual(hits[0]["identity_percent"], 100.0)

    def test_local_benchmark_alignment_hbd(self):
        hbd_seq = "ATGGTGCATCTGACTCCTGAGGAGAAGACTGCCGTCAATGCCCTGTGGGGCAAAGTGAACGTGGATGCAGTTGGTGGTGAGGCCCTGGGCAGACTGCTCGTTGTCTACCCCTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCTCTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTTCTCAGCTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCTCGCAACTTTGGCAAGGAATTCACCCCACAGATGCAGGCTGCCTATCAGAAGGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTACCACTAA"
        hits = local_benchmark_alignment(hbd_seq, "DNA")
        self.assertGreaterEqual(len(hits), 1)
        self.assertEqual(hits[0]["accession"], "NM_000519")
        self.assertEqual(hits[0]["gene_symbol"], "HBD")
        self.assertEqual(hits[0]["identity_percent"], 100.0)

    def test_local_benchmark_alignment_cd4(self):
        cd4_seq = "MNRGVPFRHLLLVLQLALLPAATQGKKVVLGKKGDTVELTCTASQKKSIQFHWKNSNQIKILGNQGSFLTKGPSKLNDRADSRRSLWDQGNFPLIIKNLKIEDSDTYICEVEDQKEEVQLLVFGLTANSDTHLLQGQSLTLTLESPPGSSPSVQCRSPRGKNIQGGKTLSVSQLELQDSGTWTCTVLQNQKKVEFKIDIVVLAFQKASSIVYKKEGEQVEFSFPLAFTVEKLTGSGELWWQAERASSSKSWITFDLKNKEVSVKRVTQDPKLQMGKKLPLHLTLPQALPQYAGSGNLTLALEAKTGKLHQEVNLVVMRATQLQKNLTCEVWGPTSPKLMLSLKLENKEAKVSKREKAVWVLNPEAGMWQCLLSDSGQVLLESNIKVLPTWSTPVQPMALIVLGGVAGLLLFIGLGIFFCVRCRHRRRQAERMSQIKRLLSEKKTCQCPHRFQKTCSPI"
        hits = local_benchmark_alignment(cd4_seq, "Protein")
        self.assertGreaterEqual(len(hits), 1)
        self.assertEqual(hits[0]["accession"], "P01730")
        self.assertEqual(hits[0]["gene_symbol"], "CD4")

    def test_blastx_translation_matching(self):
        # Coding DNA for SpCas9 N-terminal domain (MDKKYSIGLDIGTNSVGWAVITDEYKVPSKKFKVLG...)
        cas9_dna = "ATGGATAAGAAATACTCAATAGGCTTAGATATCGGCACAAATAGCGTCGGATGGGCGGTGATCACTGATGAATATAAGGTTCCGTCTAAAAAGTTCAAGGTTCTGGGAAATACAGACCGCCACAGTATCAAAAAAAATCTTATAGGGGCTCTTTTATTTGACAGTGGAGAGACAGCGGAAGCGACTCGTCTCAAACGGACAGCTCGTAGAAGGTATACACGTCGGAAGAATCGTATTTGTTATCTACAGGAGATTTTTTCAAATGAGATGGCGAAAGTAGATGATAGTTTCTTTCATCGACTTGAAGAGTCTTTTTTGGTGGAAGAAGACAAGAAGCATGAACGTCATCCTATTTTTGGAAATATAGTAGATGAAGTTGCTTATCATGAGAAATATCCAACTATCTATCATCTGCGAAAAAAATTGGTAGATTCTACTGATAAAGCGGATTTGCGCTTAATCTATTTGGCCTTAGCGCATATGATTAAGTTTCGTGGTCATTTTTTGATTGAGGGAGATTTAAATCCTGATAATAGTGATGTGGACAAACTATTTATCCAGTTGGTACAAACCTACAATCAATTATTTGAAGAAAACCCTATTAACGCAAGTGGAGTAGATGCTAAAGCGATTCTTTCTGCACGATTGAGTAAATCAAGACGATTAGAAAATCTCATTGCTCAGCTCCCCGGTGAGAAGAAAAATGGCTTATTTGGGAATCTCATTGCTTTGTCATTGGGTTTGACCCCTAATTTTAAATCAAATTTTGATTTGGCAGAAGATGCTAAATTACAGCTTTCAAAAGATACTTACGATGATGATTTAGATAATTTATTGGCGCAAATTGGAGATCAATATGCTGATTTGTTTTTGGCAGCTAAGAATTTATCAGATGCTATTTTACTTTCAGATATCCTAAGAGTAAATACTGAAATAACTAAGGCTCCCCTATCAGCTTCAATGATTAAACGCTACGATGAACATCATCAAGACTTGACTCTTTTAAAAGCTTTAGTTCGACAACAACTTCCAGAAAAGTATAAAGAAATCTTTTTTGATCAATCAAAAAACGGATATGCAGGTTATATT"
        hits = local_benchmark_alignment(cas9_dna, "DNA")
        self.assertGreaterEqual(len(hits), 1)
        self.assertEqual(hits[0]["accession"], "Q99ZW2")
        self.assertEqual(hits[0]["gene_symbol"], "CAS9")
        self.assertIn("BLASTX", hits[0]["title"])

    def test_hbd_pubmed_literature(self):
        articles = fetch_pubmed_literature("HBD", "Homo sapiens")
        self.assertGreaterEqual(len(articles), 1)
        self.assertIn("HBD", articles[0]["title"])

    def test_cd4_uniprot_annotation(self):
        res = fetch_uniprot_annotation("CD4", uniprot_id="P01730")
        self.assertTrue(res.get("found", False))
        self.assertEqual(res.get("accession"), "P01730")
        self.assertIn("CD4", res.get("protein_name", ""))

    def test_end_to_end_database_formats(self):
        from modules.validation import validate_and_analyze
        from modules.blast_search import execute_blast

        # 1. GenBank Flat file
        gb = """LOCUS       NM_000518               628 bp    mRNA    linear   PRI 01-FEB-2024
DEFINITION  Homo sapiens hemoglobin subunit beta (HBB), mRNA.
ACCESSION   NM_000518
ORIGIN      
        1 acatttgctt ctgacacaac tgtgttcact agcaacctca aacagacacc atggtgcacc
       61 tgactcctga ggagaagtct gccgttactg ccctgtgggg caaggtgaac gtggatgaag
//"""
        v_gb = validate_and_analyze(gb)
        b_gb = execute_blast(v_gb["sequence"], v_gb["type"], hints=v_gb["hints"])
        self.assertTrue(b_gb["success"])
        self.assertGreaterEqual(len(b_gb["hits"]), 1)
        self.assertIn("HBB", b_gb["hits"][0]["title"])

        # 2. UniProt Flat file
        up = """ID   HBB_HUMAN               Reviewed;         147 AA.
AC   P68871;
DE   RecName: Full=Hemoglobin subunit beta;
GN   Name=HBB;
OS   Homo sapiens (Human).
SQ   SEQUENCE   147 AA;  15998 MW;  E479DE988CBEE2DA CRC64;
     MVHLTPEEKS AVTALWGKVN VDEVGGEALG RLLVVYPWTQ RFFESFGDLS TPDAVMGNPK
     VKAHGKKVLG AFSDGLAHLD NLKGTFATLS ELHCDKLHVD PENFRLLGNV LVCVLAHHFG
//"""
        v_up = validate_and_analyze(up)
        b_up = execute_blast(v_up["sequence"], v_up["type"], hints=v_up["hints"])
        self.assertTrue(b_up["success"])
        self.assertEqual(b_up["hits"][0]["accession"], "P68871")

        # 3. 3-Letter peptide
        pep = "Met-Val-His-Leu-Thr-Pro-Glu-Glu-Lys-Ser-Ala-Val-Thr-Ala-Leu"
        v_pep = validate_and_analyze(pep)
        b_pep = execute_blast(v_pep["sequence"], v_pep["type"], hints=v_pep["hints"])
        self.assertTrue(b_pep["success"])
        self.assertGreaterEqual(len(b_pep["hits"]), 1)
        self.assertEqual(b_pep["hits"][0]["gene_symbol"], "HBB")

        # 4. PDB FASTA with taxonomy ID
        pdb = ">4HHB_1|Chains A, C|Hemoglobin subunit alpha|Homo sapiens (9606)\nVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
        v_pdb = validate_and_analyze(pdb)
        b_pdb = execute_blast(v_pdb["sequence"], v_pdb["type"], hints=v_pdb["hints"])
        self.assertTrue(b_pdb["success"])
        self.assertEqual(b_pdb["hits"][0]["gene_symbol"], "HBA1")

        # 5. Raw HBD CDS
        hbd_cds = "ATGGTGCATCTGACTCCTGAGGAGAAGACTGCCGTCAATGCCCTGTGGGGCAAAGTGAACGTGGATGCAGTTGGTGGTGAGGCCCTGGGCAGACTGCTCGTTGTCTACCCCTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCTCTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTTCTCAGCTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCTCGCAACTTTGGCAAGGAATTCACCCCACAGATGCAGGCTGCCTATCAGAAGGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTACCACTAA"
        v_hbd = validate_and_analyze(hbd_cds)
        b_hbd = execute_blast(v_hbd["sequence"], v_hbd["type"], hints=v_hbd["hints"])
        self.assertTrue(b_hbd["success"])
        self.assertEqual(b_hbd["hits"][0]["gene_symbol"], "HBD")

        # 6. Swiss-Model model.pdb (Triosephosphate Isomerase)
        tpi_seq = "RKFFVGGNWKMNGDKKSLGELIHTLNGAKLSADTEVVCGAPSIYLDFARQKLDAKIGVAAQNCYKVPKGAFTGEISPAMIKDIGAAWVILGHSERRHVFGESDELIGQKVAHALAEGLGVIACIGEKLDEREAGITEKVVFEQTKAIADNVKDWSKVVLAYEPVWAIGTGKTATPQQAQEVHEKLRGWLKSHVSDAVAQSTRIIYGGSVTGGNCKELASQHDVDGFLVGGASLKPEFVDIINAKH"
        v_tpi = validate_and_analyze(tpi_seq)
        b_tpi = execute_blast(v_tpi["sequence"], v_tpi["type"], hints=v_tpi["hints"])
        self.assertTrue(b_tpi["success"])
        self.assertEqual(b_tpi["hits"][0]["accession"], "1TPH")
        self.assertEqual(b_tpi["hits"][0]["gene_symbol"], "TPI1")
        self.assertEqual(b_tpi["hits"][0]["identity_percent"], 100.0)

if __name__ == "__main__":
    unittest.main()
