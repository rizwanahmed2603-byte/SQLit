"""
Sequence validation and statistics module for SeqLit.
Handles raw sequence and FASTA cleaning, sequence type detection (DNA, RNA, Protein),
degenerate IUPAC alphabet support, and computation of basic metrics
(GC%, molecular weight, composition frequency).
"""

import re
from typing import Dict, Any, Tuple

# Standard IUPAC Monomer Weights (approximate average g/mol)
DNA_WEIGHTS = {'A': 313.21, 'T': 304.2, 'C': 289.18, 'G': 329.21, 'N': 308.95}
RNA_WEIGHTS = {'A': 329.21, 'U': 306.17, 'C': 305.18, 'G': 345.21, 'N': 321.44}
PROTEIN_WEIGHTS = {
    'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.16,
    'E': 147.13, 'Q': 146.15, 'G': 75.07, 'H': 155.16, 'I': 131.18,
    'L': 131.18, 'K': 146.19, 'M': 149.21, 'F': 165.19, 'P': 115.13,
    'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15,
    'X': 110.0, 'B': 132.6, 'Z': 146.6, 'U': 168.05, 'O': 255.31
}

# Standard & degenerate IUPAC character sets
DNA_CHARS = set("ATCGN")
DNA_DEGENERATE = set("ATCGNRYKMSWBDHVX")
RNA_CHARS = set("AUCGN")
RNA_DEGENERATE = set("AUCGNRYKMSWBDHVX")
PROTEIN_CHARS = set("ACDEFGHIKLMNPQRSTVWY")
PROTEIN_EXTENDED = set("ACDEFGHIKLMNPQRSTVWYXBZJUO")

def clean_sequence(raw_input: str) -> Tuple[str, str]:
    """
    Strips FASTA header if present, removes whitespace, line numbers,
    and strips trailing stop codons (*).
    Returns (cleaned_sequence_upper, header_or_empty).
    """
    if not raw_input or not raw_input.strip():
        return "", ""

    lines = raw_input.strip().splitlines()
    header = ""
    seq_lines = []

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            if not header:
                header = stripped[1:].strip()
        else:
            # remove line numbers or spaces often found in genbank files
            cleaned_line = re.sub(r"[\s\d]", "", stripped)
            seq_lines.append(cleaned_line)

    sequence = "".join(seq_lines).upper()
    # Strip terminal stop codon if present e.g. from translation output
    sequence = sequence.rstrip("*")
    return sequence, header

def detect_sequence_type(seq: str) -> str:
    """
    Detects whether the sequence is DNA, RNA, Protein, or Invalid.
    Supports standard and degenerate IUPAC codes.
    """
    if not seq:
        return "Unknown"

    unique_chars = set(seq)

    # 1. RNA detection: 'U' present and no 'T'
    if 'U' in unique_chars and 'T' not in unique_chars:
        if unique_chars.issubset(RNA_DEGENERATE):
            return "RNA"

    # 2. DNA detection: ATCGN or degenerate DNA
    if unique_chars.issubset(DNA_CHARS) or unique_chars.issubset(DNA_DEGENERATE):
        return "DNA"

    # 3. If both T and U are present, check if it fits protein
    if 'T' in unique_chars and 'U' in unique_chars:
        if unique_chars.issubset(PROTEIN_EXTENDED):
            return "Protein"
        return "Invalid"

    # 4. Standard and extended Protein
    if unique_chars.issubset(PROTEIN_CHARS) or unique_chars.issubset(PROTEIN_EXTENDED):
        return "Protein"

    return "Invalid"

def calculate_stats(seq: str, seq_type: str) -> Dict[str, Any]:
    """
    Calculates sequence statistics including length, GC%, molecular weight,
    and base/residue composition.
    """
    length = len(seq)
    if length == 0:
        return {
            "length": 0,
            "gc_percent": 0.0,
            "molecular_weight": 0.0,
            "composition": {},
            "frequencies": {}
        }

    composition = {}
    for char in seq:
        composition[char] = composition.get(char, 0) + 1

    frequencies = {k: round((v / length) * 100, 2) for k, v in composition.items()}

    gc_percent = None
    molecular_weight = 0.0

    if seq_type == "DNA":
        g_count = composition.get('G', 0)
        c_count = composition.get('C', 0)
        gc_percent = round(((g_count + c_count) / length) * 100, 2)
        weight = sum(DNA_WEIGHTS.get(base, 308.0) * count for base, count in composition.items())
        if length > 1:
            weight -= (length - 1) * 18.015
        molecular_weight = round(weight, 2)

    elif seq_type == "RNA":
        g_count = composition.get('G', 0)
        c_count = composition.get('C', 0)
        gc_percent = round(((g_count + c_count) / length) * 100, 2)
        weight = sum(RNA_WEIGHTS.get(base, 321.0) * count for base, count in composition.items())
        if length > 1:
            weight -= (length - 1) * 18.015
        molecular_weight = round(weight, 2)

    elif seq_type == "Protein":
        weight = sum(PROTEIN_WEIGHTS.get(aa, 110.0) * count for aa, count in composition.items())
        if length > 1:
            weight -= (length - 1) * 18.015
        molecular_weight = round(weight, 2)

    return {
        "length": length,
        "gc_percent": gc_percent,
        "molecular_weight": molecular_weight,
        "composition": composition,
        "frequencies": frequencies
    }

def extract_header_hints(header: str) -> Dict[str, str]:
    """
    Extracts biological hints (accession, gene symbol, organism) from FASTA header.
    """
    hints = {"accession": "", "gene": "", "organism": ""}
    if not header:
        return hints

    # Match Swiss-Prot / UniProt accession: sp|P01308| or standalone UniProt ID
    sp_m = re.search(r"sp\|([A-Z0-9]{6,10})", header)
    if sp_m:
        hints["accession"] = sp_m.group(1)
    else:
        uniprot_stand = re.search(r"\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})\b", header)
        if uniprot_stand:
            hints["accession"] = uniprot_stand.group(1)

    # Match RefSeq / GenBank / Ensembl accession
    if not hints["accession"]:
        ref_m = re.search(r"\b([A-Z]{1,2}_\d{5,9}|[A-Z]{1,3}\d{5,8}|ENS[A-Z]{0,4}[GTP]\d{9,11})\b", header)
        if ref_m:
            hints["accession"] = ref_m.group(1)

    # Match Organism: OS=... or [Homo sapiens] or canonical names
    org_os = re.search(r"OS=([^=]+?)(?:\s+[A-Z]{2}=|$)", header)
    if org_os:
        hints["organism"] = org_os.group(1).strip()
    else:
        org_bracket = re.search(r"\[([A-Za-z0-9_.\s-]+)\]", header)
        if org_bracket:
            hints["organism"] = org_bracket.group(1).strip()
        else:
            org_common = re.search(r"\b(Homo sapiens|Mus musculus|Escherichia coli|Aequorea victoria|Gallus gallus|Physeter catodon|Epstein-Barr virus|SARS-CoV-2)\b", header, re.IGNORECASE)
            if org_common:
                hints["organism"] = org_common.group(1).strip()

    # Match Gene: GN=HBB, gene=HBB, (HBB), or leading gene symbol
    gn_m = re.search(r"(?:GN=|gene=|\()([A-Za-z0-9_-]{2,12})(?:\)|\s|$)", header, re.IGNORECASE)
    if gn_m:
        hints["gene"] = gn_m.group(1).strip().upper()
    else:
        # Check known common symbols in header
        known_genes = re.search(r"\b(HBB|HBA1|HBA|INS|TP53|GAPDH|ACTB|ALB|EGFR|BRCA1|MYC|CYCS|GFP|LYZ|MB|EBNA1|APOE|TNF|IL6|ACE2)\b", header, re.IGNORECASE)
        if known_genes:
            hints["gene"] = known_genes.group(1).strip().upper()

    return hints

def validate_and_analyze(raw_input: str) -> Dict[str, Any]:
    """
    Top-level validation and statistics function for SeqLit.
    """
    cleaned_seq, header = clean_sequence(raw_input)
    if not cleaned_seq:
        return {
            "valid": False,
            "error": "Sequence input is empty or contains only headers.",
            "sequence": "",
            "header": header
        }

    seq_type = detect_sequence_type(cleaned_seq)
    if seq_type == "Invalid":
        all_allowed = DNA_DEGENERATE | RNA_DEGENERATE | PROTEIN_EXTENDED
        invalid_chars = sorted(list(set(cleaned_seq) - all_allowed))
        return {
            "valid": False,
            "error": f"Sequence contains unrecognized or invalid characters: {', '.join(invalid_chars)}",
            "sequence": cleaned_seq,
            "header": header,
            "type": "Invalid"
        }

    stats = calculate_stats(cleaned_seq, seq_type)
    header_hints = extract_header_hints(header)

    return {
        "valid": True,
        "header": header or "User Sequence",
        "sequence": cleaned_seq,
        "type": seq_type,
        "stats": stats,
        "hints": header_hints
    }
