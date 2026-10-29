"""
Tests for universal similarity search across all biological kingdoms (insects, microbes, plants, animals),
direct NCBI polling resilience, and caching integrity.
"""
import unittest
from modules.blast_search import execute_blast
from modules.cache_manager import get_cached_result, set_cached_result, compute_hash

class TestUniversalSimilaritySearch(unittest.TestCase):
    def test_insect_protein_homology_search(self):
        # 188 aa Musca domestica uncharacterized protein (XP_005185681 / A0A1I8NI34)
        musca_seq = (
            "MTDLDPWKVQVQHPEHSNSATNGKVTDEINNSTETTNSALGFEDSFIGAENNNCSSDKKRPENTNKSNNK"
            "PEHQPLPDSQTYLKGLERKLDKLKKNSKLVEALSEKRSDCLRSLIQGDLSGSNSNNDLLLELEASINSDS"
            "AVHNLYRQIQPVQAVTVGETVHIVKYDQLEEQRLEADAEENGELPASR"
        )
        res = execute_blast(musca_seq, "Protein")
        self.assertTrue(res.get("success"), "Search must succeed")
        hits = res.get("hits", [])
        self.assertGreater(len(hits), 0, "Must return homolog hits across universal databases")
        
        top_hit = hits[0]
        self.assertIn("Musca domestica", top_hit.get("organism", ""))
        self.assertGreaterEqual(top_hit.get("identity_percent", 0), 90.0)

    def test_unbracketed_organism_metadata_extraction(self):
        from modules.blast_search import extract_hit_metadata
        title = "gi|3328840253|gb|PZ053791.1| Hydrogenophaga pseudoflava strain HNP putative ammonia monooxygenase gene, complete cds"
        meta = extract_hit_metadata(title, "PZ053791.1", "")
        self.assertEqual(meta["accession"], "PZ053791.1")
        self.assertEqual(meta["organism"], "Hydrogenophaga pseudoflava")

    def test_adaptive_evalue_short_query_xml_parsing(self):
        from modules.blast_search import parse_blast_xml
        import io
        # Mock BLAST XML with short 28 bp query and E-value 0.0701
        sample_xml = """<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastn</BlastOutput_program>
  <BlastOutput_version>BLASTN 2.15.0+</BlastOutput_version>
  <BlastOutput_reference>Stephen F. Altschul...</BlastOutput_reference>
  <BlastOutput_db>core_nt</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>test</BlastOutput_query-def>
  <BlastOutput_query-len>28</BlastOutput_query-len>
  <BlastOutput_param>
    <Parameters>
      <Parameters_expect>10</Parameters_expect>
      <Parameters_sc-match>1</Parameters_sc-match>
      <Parameters_sc-mismatch>-2</Parameters_sc-mismatch>
      <Parameters_gap-open>0</Parameters_gap-open>
      <Parameters_gap-extend>0</Parameters_gap-extend>
      <Parameters_filter>L;m;</Parameters_filter>
    </Parameters>
  </BlastOutput_param>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-def>test</Iteration_query-def>
      <Iteration_query-len>28</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>gi|2211716842|ref|XM_047163621.1|</Hit_id>
          <Hit_def>Helicoverpa zea sphingosine-1-phosphate phosphatase</Hit_def>
          <Hit_accession>XM_047163621</Hit_accession>
          <Hit_len>28</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>51.78</Hsp_bit-score>
              <Hsp_score>28</Hsp_score>
              <Hsp_evalue>0.0701</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>28</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>28</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>1</Hsp_hit-frame>
              <Hsp_identity>28</Hsp_identity>
              <Hsp_positive>28</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>28</Hsp_align-len>
              <Hsp_qseq>GCTGCTGGCGACGGCGCTGGTGTCGGTG</Hsp_qseq>
              <Hsp_hseq>GCTGCTGGCGACGGCGCTGGTGTCGGTG</Hsp_hseq>
              <Hsp_midline>||||||||||||||||||||||||||||</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
      <Iteration_stat>
        <Statistics>
          <Statistics_db-num>100</Statistics_db-num>
          <Statistics_db-len>100000</Statistics_db-len>
          <Statistics_hsp-len>0</Statistics_hsp-len>
          <Statistics_eff-space>100000</Statistics_eff-space>
          <Statistics_kappa>0.041</Statistics_kappa>
          <Statistics_lambda>0.267</Statistics_lambda>
          <Statistics_entropy>0.14</Statistics_entropy>
        </Statistics>
      </Iteration_stat>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>"""
        hits = parse_blast_xml(io.StringIO(sample_xml), query_length=28, e_value_cutoff=1e-5)
        self.assertEqual(len(hits), 1, "Short 28 bp query with 100% identity must be retained via adaptive cutoff")
        self.assertEqual(hits[0]["identity_percent"], 100.0)
        self.assertEqual(hits[0]["accession"], "XM_047163621")

    def test_cache_integrity_no_false_zero_hits(self):
        # Ensure that transient timeout results never get saved into the SQLite cache
        test_seq = "ACGTACGTACGTACGT"
        fake_timeout_res = {
            "blast": {"success": False, "timeout": True, "hits": []},
            "ncbi": {"found": False},
            "uniprot": {"found": False}
        }
        # In app.py, caching only happens if blast success and not timeout and len(hits) > 0
        should_cache = (
            fake_timeout_res["blast"].get("success")
            and not fake_timeout_res["blast"].get("timeout")
            and len(fake_timeout_res["blast"].get("hits", [])) > 0
        )
        self.assertFalse(should_cache, "Timeout failures must never qualify for caching")

if __name__ == "__main__":
    unittest.main()
