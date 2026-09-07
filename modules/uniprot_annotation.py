"""
UniProt Annotation module for SeqLit.
Queries the official UniProt REST API (https://rest.uniprot.org) to retrieve
real protein function descriptions, Gene Ontology (GO) terms, domains, and cross-references.
Never produces mock data if no real UniProt entry exists.
"""

import logging
from typing import Dict, Any, List

try:
    import requests
    REQUESTS_AVAILABLE = True
except ImportError:
    REQUESTS_AVAILABLE = False

logger = logging.getLogger(__name__)

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"

def fetch_uniprot_annotation(query_term: str) -> Dict[str, Any]:
    """
    Queries UniProtKB search for gene name or accession.
    Returns real structured annotations, or found=False if not found.
    """
    if not query_term or query_term.strip() in ("", "N/A", "Unknown", "Candidate Gene"):
        return {
            "found": False,
            "protein_name": "N/A",
            "gene": "N/A",
            "organism": "N/A",
            "function": "No candidate gene identified for UniProt search.",
            "go_terms": {"biological_process": [], "molecular_function": [], "cellular_component": []},
            "features": [],
            "entry_url": ""
        }

    if not REQUESTS_AVAILABLE:
        return {"found": False, "error": "requests library not installed"}

    try:
        params = {
            "query": f"{query_term} AND reviewed:true",
            "format": "json",
            "size": 1
        }
        headers = {"Accept": "application/json"}
        response = requests.get(UNIPROT_SEARCH_URL, params=params, headers=headers, timeout=10)

        if response.status_code == 200:
            data = response.json()
            if not data.get("results"):
                # Try unreviewed
                params["query"] = query_term
                response = requests.get(UNIPROT_SEARCH_URL, params=params, headers=headers, timeout=10)
                data = response.json()

            results = data.get("results", [])
            if not results:
                return {
                    "found": False,
                    "protein_name": "Not found in UniProt",
                    "gene": query_term,
                    "function": f"No UniProt entries matched '{query_term}'.",
                    "go_terms": {"biological_process": [], "molecular_function": [], "cellular_component": []},
                    "features": [],
                    "entry_url": ""
                }

            entry = results[0]
            return _parse_uniprot_entry(entry)
        else:
            logger.warning(f"UniProt REST API error: status {response.status_code}")
            return {
                "found": False,
                "protein_name": "API Error",
                "function": f"UniProt responded with status code {response.status_code}",
                "go_terms": {"biological_process": [], "molecular_function": [], "cellular_component": []},
                "features": [],
                "entry_url": ""
            }
    except Exception as e:
        logger.warning(f"UniProt request exception: {e}")
        return {
            "found": False,
            "protein_name": "Connection Error",
            "function": f"Could not connect to UniProt API: {str(e)}",
            "go_terms": {"biological_process": [], "molecular_function": [], "cellular_component": []},
            "features": [],
            "entry_url": ""
        }

def _parse_uniprot_entry(entry: Dict[str, Any]) -> Dict[str, Any]:
    accession = entry.get("primaryAccession", "")
    protein_desc = entry.get("proteinDescription", {})
    recommended_name = (
        protein_desc.get("recommendedName", {})
        .get("fullName", {})
        .get("value", "")
    )
    if not recommended_name and protein_desc.get("submissionNames"):
        recommended_name = protein_desc["submissionNames"][0].get("fullName", {}).get("value", "")

    genes = []
    for g in entry.get("genes", []):
        if "geneName" in g:
            genes.append(g["geneName"].get("value", ""))
    gene_symbol = ", ".join(genes) if genes else "N/A"

    organism = entry.get("organism", {}).get("scientificName", "Unknown")

    function_texts = []
    for comment in entry.get("comments", []):
        if comment.get("commentType") == "FUNCTION":
            for text_obj in comment.get("texts", []):
                function_texts.append(text_obj.get("value", ""))
    function_summary = " ".join(function_texts) or "No functional description recorded."

    go_terms = {"biological_process": [], "molecular_function": [], "cellular_component": []}
    for db_ref in entry.get("uniProtKBCrossReferences", []):
        if db_ref.get("database") == "GO":
            go_id = db_ref.get("id", "")
            properties = {p.get("key"): p.get("value") for p in db_ref.get("properties", [])}
            go_term_str = properties.get("GoTerm", "")
            if go_term_str.startswith("P:"):
                go_terms["biological_process"].append({"id": go_id, "name": go_term_str[2:]})
            elif go_term_str.startswith("F:"):
                go_terms["molecular_function"].append({"id": go_id, "name": go_term_str[2:]})
            elif go_term_str.startswith("C:"):
                go_terms["cellular_component"].append({"id": go_id, "name": go_term_str[2:]})

    features = []
    for feat in entry.get("features", []):
        f_type = feat.get("type", "")
        f_desc = feat.get("description", "")
        loc = feat.get("location", {})
        start = loc.get("start", {}).get("value", "")
        end = loc.get("end", {}).get("value", "")
        if f_type in ("Domain", "Region", "Active site", "Binding site", "Chain"):
            features.append({
                "type": f_type,
                "description": f_desc or f_type,
                "start": start,
                "end": end
            })

    return {
        "found": True,
        "accession": accession,
        "protein_name": recommended_name or "Uncharacterized protein",
        "gene": gene_symbol,
        "organism": organism,
        "function": function_summary,
        "go_terms": go_terms,
        "features": features[:15],
        "entry_url": f"https://www.uniprot.org/uniprotkb/{accession}"
    }
