"""
NCBI Annotation module for SeqLit.
Uses Entrez E-utilities (efetch / esummary) to retrieve real organism, taxonomy,
gene name, and sequence metadata for confirmed hits.
If no hit or invalid accession is provided, returns empty/not-found status.
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

DEFAULT_EMAIL = "seqlit.student@academic.edu"

def fetch_ncbi_annotation(accession: str, seq_type: str, email: str = DEFAULT_EMAIL, api_key: Optional[str] = None) -> Dict[str, Any]:
    """
    Fetches real GenBank record for accession and extracts organism, gene symbol,
    and taxonomy lineage. Returns found=False if accession is empty or not in NCBI.
    """
    if not accession or accession.strip() == "":
        return {
            "found": False,
            "accession": "",
            "organism": "N/A",
            "gene": "N/A",
            "taxonomy": "N/A",
            "title": "No homolog accession identified.",
            "source": "None"
        }

    if not ENTREZ_AVAILABLE:
        return {
            "found": False,
            "accession": accession,
            "error": "Bio.Entrez is not installed."
        }

    db = "nuccore" if seq_type in ("DNA", "RNA") else "protein"
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    try:
        time.sleep(0.35)
        handle = Entrez.efetch(db=db, id=accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()

        organism = record.annotations.get("organism", "Unknown")
        taxonomy_list = record.annotations.get("taxonomy", [])
        taxonomy = "; ".join(taxonomy_list) if taxonomy_list else "Unknown"

        gene_name = ""
        for feature in record.features:
            if "gene" in feature.qualifiers:
                gene_name = feature.qualifiers["gene"][0]
                break

        return {
            "found": True,
            "accession": accession,
            "organism": organism,
            "taxonomy": taxonomy,
            "gene": gene_name or record.name,
            "title": record.description,
            "seq_length": len(record.seq),
            "source": "NCBI GenBank E-utilities"
        }
    except Exception as e:
        logger.warning(f"Failed to fetch real NCBI annotation for '{accession}': {e}")
        return {
            "found": False,
            "accession": accession,
            "organism": "Not found",
            "taxonomy": "Not found",
            "gene": "N/A",
            "title": f"Could not retrieve GenBank record for {accession}",
            "error": str(e),
            "source": "NCBI GenBank E-utilities"
        }
