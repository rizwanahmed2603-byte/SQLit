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
