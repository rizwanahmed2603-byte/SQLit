"""
PubMed literature retrieval & rule-based topic categorization module for SeqLit.
Uses Entrez esearch & esummary/efetch to find real articles by gene and organism,
de-duplicates by PMID, and assigns papers to one of five key research categories.
Never returns fake or fabricated articles if no real publications exist.
"""

import time
import re
import logging
from typing import List, Dict, Any, Optional

try:
    from Bio import Entrez
    ENTREZ_AVAILABLE = True
except ImportError:
    ENTREZ_AVAILABLE = False

logger = logging.getLogger(__name__)
DEFAULT_EMAIL = "seqlit.student@academic.edu"

CATEGORY_RULES = {
    "Structure / Crystallography": re.compile(
        r"\b(crystal|structure|structural|pdb|x-ray|cryo-em|nmr|conformation|folding|domain|atomic)\b",
        re.IGNORECASE
    ),
    "Drug / Pharmacology": re.compile(
        r"\b(drug|inhibitor|inhibition|therapy|therapeutic|pharmacology|compound|ligand|agonist|antagonist|target)\b",
        re.IGNORECASE
    ),
    "Disease / Pathology": re.compile(
        r"\b(disease|cancer|tumor|mutation|pathology|syndrome|patient|clinical|disorder|risk|variant|carcinoma)\b",
        re.IGNORECASE
    ),
    "Function / Mechanism": re.compile(
        r"\b(pathway|mechanism|function|activity|catalytic|signaling|regulation|kinase|expression|metabolism|synthesis)\b",
        re.IGNORECASE
    ),
    "Interaction / Complex": re.compile(
        r"\b(binding|bound|complex|interact|interaction|subunit|dimer|oligomer|partner|cross-link)\b",
        re.IGNORECASE
    )
}

def categorize_article(title: str, abstract: str = "") -> str:
    """
    Applies rule-based regex keyword matching against title and abstract text.
    """
    text = f"{title} {abstract}"
    scores = {}

    for cat_name, pattern in CATEGORY_RULES.items():
        matches = len(pattern.findall(text))
        if matches > 0:
            scores[cat_name] = matches

    if not scores:
        return "General / Uncategorized"

    return max(scores, key=scores.get)

def fetch_pubmed_literature(gene_name: str, organism: str = "", max_results: int = 15, email: str = DEFAULT_EMAIL) -> List[Dict[str, Any]]:
    """
    Retrieves real PubMed articles for given gene and organism.
    Returns empty list if gene name is empty, invalid, or no papers are found.
    """
    if not gene_name or gene_name.strip() in ("", "N/A", "Unknown", "Candidate Gene"):
        return []

    if not ENTREZ_AVAILABLE:
        return []

    Entrez.email = email
    articles = []
    seen_pmids = set()

    query_parts = [f"{gene_name}[Title/Abstract]"]
    if organism and organism not in ("Unknown", "Not found", "Homo sapiens (estimated)"):
        query_parts.append(f"{organism}[Organism]")
    term = " AND ".join(query_parts)

    try:
        time.sleep(0.35)
        search_handle = Entrez.esearch(db="pubmed", term=term, retmax=max_results, sort="relevance")
        search_results = Entrez.read(search_handle)
        search_handle.close()

        id_list = search_results.get("IdList", [])
        if not id_list:
            # Fallback to search without organism restriction if too narrow
            search_handle = Entrez.esearch(db="pubmed", term=f"{gene_name}[Title/Abstract]", retmax=max_results, sort="pub_date")
            search_results = Entrez.read(search_handle)
            search_handle.close()
            id_list = search_results.get("IdList", [])

        if not id_list:
            return []

        time.sleep(0.35)
        summary_handle = Entrez.esummary(db="pubmed", id=",".join(id_list))
        summaries = Entrez.read(summary_handle)
        summary_handle.close()

        for record in summaries:
            pmid = str(record.get("Id", ""))
            if pmid in seen_pmids:
                continue
            seen_pmids.add(pmid)

            title = record.get("Title", "No title available")
            authors_list = record.get("AuthorList", [])
            authors_str = ", ".join(authors_list[:3])
            if len(authors_list) > 3:
                authors_str += " et al."

            journal = record.get("Source", "Unknown Journal")
            pub_date = record.get("PubDate", "")
            year = pub_date[:4] if pub_date else "N/A"
            category = categorize_article(title, "")

            articles.append({
                "pmid": pmid,
                "title": title.rstrip("."),
                "authors": authors_str or "Unknown authors",
                "journal": journal,
                "year": year,
                "category": category,
                "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            })

        return articles

    except Exception as e:
        logger.warning(f"PubMed retrieval error: {e}")
        return []
