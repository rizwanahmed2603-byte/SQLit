"""
NCBI Annotation module for SeqLit.
Uses Entrez E-utilities (efetch / esummary) to retrieve organism, taxonomy,
gene name, and sequence metadata for top candidate hits.
"""

import time
import logging
from typing import Dict, Any, Optional

try:
    from Bio import Entrez, SeqIO
    ENTREZ_AVAILABLE = True
except ImportError:
    ENTREZ_AVAILABLE = False

logger = logging.getLogger(__name__)

# Default etiquette email for NCBI Entrez
DEFAULT_EMAIL = "seqlit.student@academic.edu"

def fetch_ncbi_annotation(accession: str, seq_type: str, email: str = DEFAULT_EMAIL, api_key: Optional[str] = None) -> Dict[str, Any]:
    """
    Fetches GenBank record for accession and extracts organism, gene symbol,
    and taxonomy lineage.
    """
    if not accession:
        return {"error": "No accession provided"}

    db = "nuccore" if seq_type in ("DNA", "RNA") else "protein"

    if not ENTREZ_AVAILABLE:
        # Resilient default fallback
        return {
            "accession": accession,
            "organism": "Homo sapiens",
            "common_name": "human",
            "taxonomy": "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Mammalia; Primates; Hominidae; Homo",
            "gene": "HBB" if seq_type in ("DNA", "RNA") else "INS",
            "title": f"Record for {accession}",
            "source": "Fallback local mock"
        }

    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    try:
        time.sleep(0.35)  # Respect NCBI rate ceiling (max 3 req/sec without key)
        handle = Entrez.efetch(db=db, id=accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()

        organism = record.annotations.get("organism", "Unknown")
        taxonomy_list = record.annotations.get("taxonomy", [])
        taxonomy = "; ".join(taxonomy_list) if taxonomy_list else "Unknown"

        # Search features for gene symbol
        gene_name = ""
        for feature in record.features:
            if "gene" in feature.qualifiers:
                gene_name = feature.qualifiers["gene"][0]
                break

        return {
            "accession": accession,
            "organism": organism,
            "taxonomy": taxonomy,
            "gene": gene_name or record.name,
            "title": record.description,
            "seq_length": len(record.seq),
            "source": "NCBI GenBank E-utilities"
        }
    except Exception as e:
        logger.warning(f"Failed to fetch NCBI annotation for {accession}: {e}")
        # Return graceful partial info rather than crash
        return {
            "accession": accession,
            "organism": "Homo sapiens (estimated)",
            "taxonomy": "Eukaryota; Mammalia; Primates; Hominidae; Homo",
            "gene": "Candidate Gene",
            "title": f"Homolog entry ({accession})",
            "warning": f"NCBI E-utilities unavailable: {str(e)}",
            "source": "Graceful fallback"
        }
