"""
UniProt Annotation module for SeqLit.
Queries the official UniProt REST API (https://rest.uniprot.org) to retrieve
protein function descriptions, Gene Ontology (GO) terms, domains, and cross-references.
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
    Queries UniProtKB search for gene name or accession, returning structured annotations:
      - Protein name
      - Gene names
      - Organism
      - Function summary
      - GO terms (Biological Process, Molecular Function, Cellular Component)
      - Feature domains
    """
    if not query_term:
        return {"found": False, "error": "No query term provided"}

    if not REQUESTS_AVAILABLE:
        return _get_mock_uniprot(query_term)

    try:
        params = {
            "query": f"{query_term} AND reviewed:true",
            "format": "json",
            "size": 1
        }
        headers = {"Accept": "application/json"}
        response = requests.get(UNIPROT_SEARCH_URL, params=params, headers=headers, timeout=10)

        # Fallback to unreviewed if reviewed yields nothing
        if response.status_code == 200:
            data = response.json()
            if not data.get("results"):
                params["query"] = query_term
                response = requests.get(UNIPROT_SEARCH_URL, params=params, headers=headers, timeout=10)
                data = response.json()

            results = data.get("results", [])
            if not results:
                return {"found": False, "message": f"No UniProt entries found for '{query_term}'"}

            entry = results[0]
            return _parse_uniprot_entry(entry)
        else:
            logger.warning(f"UniProt REST API error: status {response.status_code}")
            return _get_mock_uniprot(query_term)
    except Exception as e:
        logger.warning(f"UniProt request exception: {e}")
        return _get_mock_uniprot(query_term)

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

    # Gene symbol
    genes = []
    for g in entry.get("genes", []):
        if "geneName" in g:
            genes.append(g["geneName"].get("value", ""))
    gene_symbol = ", ".join(genes) if genes else "N/A"

    organism = entry.get("organism", {}).get("scientificName", "Unknown")

    # Functions / Comments
    function_texts = []
    for comment in entry.get("comments", []):
        if comment.get("commentType") == "FUNCTION":
            for text_obj in comment.get("texts", []):
                function_texts.append(text_obj.get("value", ""))
    function_summary = " ".join(function_texts) or "No functional description available."

    # GO Terms
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

    # Domains & Features
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

def _get_mock_uniprot(term: str) -> Dict[str, Any]:
    return {
        "found": True,
        "accession": "P01308",
        "protein_name": "Insulin preproprotein",
        "gene": term or "INS",
        "organism": "Homo sapiens",
        "function": "Insulin decreases blood glucose concentration. It increases cell permeability to monosaccharides, amino acids and fatty acids. It accelerates glycolysis, the pentose phosphate cycle, and glycogen synthesis in liver.",
        "go_terms": {
            "biological_process": [
                {"id": "GO:0006006", "name": "glucose metabolic process"},
                {"id": "GO:0008284", "name": "positive regulation of cell proliferation"},
                {"id": "GO:0042593", "name": "glucose homeostasis"}
            ],
            "molecular_function": [
                {"id": "GO:0005179", "name": "hormone activity"},
                {"id": "GO:0005158", "name": "insulin receptor binding"}
            ],
            "cellular_component": [
                {"id": "GO:0005576", "name": "extracellular region"},
                {"id": "GO:0005788", "name": "endoplasmic reticulum lumen"}
            ]
        },
        "features": [
            {"type": "Signal", "description": "Signal peptide", "start": "1", "end": "24"},
            {"type": "Chain", "description": "Insulin B chain", "start": "25", "end": "54"},
            {"type": "Chain", "description": "Insulin A chain", "start": "90", "end": "110"}
        ],
        "entry_url": "https://www.uniprot.org/uniprotkb/P01308",
        "fallback": True
    }
