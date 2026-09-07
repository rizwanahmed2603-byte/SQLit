"""
BLAST similarity search and hit parsing module for SeqLit.
Executes real remote similarity search via NCBIWWW.qblast and parses
the resulting XML output.
If a sequence has no significant biological matches in NCBI nr/nt databases,
it accurately returns zero hits rather than fabricating fake homologs.
"""

import io
import logging
from typing import List, Dict, Any, Optional

try:
    from Bio.Blast import NCBIWWW, NCBIXML
    from Bio import Entrez
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

try:
    from tenacity import retry, stop_after_attempt, wait_exponential
    TENACITY_AVAILABLE = True
except ImportError:
    TENACITY_AVAILABLE = False

logger = logging.getLogger(__name__)

def parse_blast_xml(xml_handle, query_length: int, e_value_cutoff: float = 1e-5) -> List[Dict[str, Any]]:
    """
    Parses BLAST XML results and returns ranked, formatted hit dictionaries.
    """
    hits = []
    try:
        blast_record = NCBIXML.read(xml_handle)
    except Exception as e:
        logger.error(f"Error reading BLAST XML: {e}")
        return hits

    rank = 1
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect <= e_value_cutoff:
                identity_percent = round((hsp.identities / hsp.align_length) * 100, 2)
                q_len = query_length if query_length > 0 else (hsp.query_end - hsp.query_start + 1)
                coverage = round(((hsp.query_end - hsp.query_start + 1) / q_len) * 100, 2)
                coverage = min(100.0, coverage)

                raw_title = alignment.title
                accession = alignment.accession
                if not accession:
                    parts = raw_title.split("|")
                    if len(parts) >= 4:
                        accession = parts[3]
                    else:
                        accession = raw_title.split()[0]

                hit = {
                    "rank": rank,
                    "accession": accession,
                    "title": alignment.hit_def or raw_title,
                    "definition": alignment.hit_def or raw_title,
                    "e_value": hsp.expect,
                    "identity_percent": identity_percent,
                    "bit_score": round(hsp.bits, 1),
                    "alignment_length": hsp.align_length,
                    "query_coverage": coverage,
                    "query_start": hsp.query_start,
                    "query_end": hsp.query_end,
                    "sbjct_start": hsp.sbjct_start,
                    "sbjct_end": hsp.sbjct_end
                }
                hits.append(hit)
                rank += 1
                break
        if rank > 25:
            break

    return hits

def execute_blast(sequence: str, seq_type: str, e_value_cutoff: float = 1e-5) -> Dict[str, Any]:
    """
    Executes real similarity search via NCBI BLAST.
    Returns real hits or empty list if no matches exist.
    Never returns fake/mock results.
    """
    program = "blastn" if seq_type in ("DNA", "RNA") else "blastp"
    database = "nt" if program == "blastn" else "nr"
    q_len = len(sequence)

    if not BIOPYTHON_AVAILABLE:
        return {
            "success": False,
            "program": program,
            "database": database,
            "hits": [],
            "error": "Biopython is not installed on the server."
        }

    try:
        logger.info(f"Submitting real {program} query ({q_len} bp/aa) to NCBI...")
        result_handle = NCBIWWW.qblast(
            program=program,
            database=database,
            sequence=sequence,
            expect=e_value_cutoff,
            hitlist_size=20,
            format_type="XML"
        )
        hits = parse_blast_xml(result_handle, query_length=q_len, e_value_cutoff=e_value_cutoff)
        
        if not hits:
            return {
                "success": True,
                "program": program,
                "database": database,
                "hits": [],
                "message": "No significant similarity hits found in NCBI database for this sequence."
            }

        return {
            "success": True,
            "program": program,
            "database": database,
            "hits": hits,
            "message": f"Successfully retrieved {len(hits)} significant BLAST hits from NCBI."
        }

    except Exception as e:
        logger.error(f"NCBI BLAST request failed: {e}")
        return {
            "success": False,
            "program": program,
            "database": database,
            "hits": [],
            "error": f"NCBI BLAST service error: {str(e)}"
        }
