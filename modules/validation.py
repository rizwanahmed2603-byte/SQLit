"""
Sequence validation and statistics module for SeqLit.
Handles raw sequence and FASTA cleaning, sequence type detection (DNA, RNA, Protein),
and computation of basic metrics (GC%, molecular weight, composition frequency).
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
    'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15
}

DNA_CHARS = set("ATCGN")
RNA_CHARS = set("AUCGN")
PROTEIN_CHARS = set("ACDEFGHIKLMNPQRSTVWY")

def clean_sequence(raw_input: str) -> Tuple[str, str]:
    """
    Strips FASTA header if present and removes all whitespace/numbers.
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
    return sequence, header

def detect_sequence_type(seq: str) -> str:
    """
    Detects whether the sequence is DNA, RNA, Protein, or Invalid.
    """
    if not seq:
        return "Unknown"

    unique_chars = set(seq)

    # Check for RNA strictly before DNA if 'U' is present without 'T'
    if 'U' in unique_chars and 'T' not in unique_chars:
        if unique_chars.issubset(RNA_CHARS):
            return "RNA"

    # Check for DNA: only ATCGN
    if unique_chars.issubset(DNA_CHARS):
        return "DNA"

    # If both T and U are present, it's ambiguous/invalid
    if 'T' in unique_chars and 'U' in unique_chars:
        # Check if it fits protein
        if unique_chars.issubset(PROTEIN_CHARS):
            return "Protein"
        return "Invalid"

    # Check for Protein: standard 20 amino acids
    if unique_chars.issubset(PROTEIN_CHARS):
        # Even if made of A, C, G, T it fell into DNA above.
        # So reaching here means it contains non-DNA amino acids.
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
        # Approximate MW: sum of dNTPs minus water loss in phosphodiester bonds
        weight = sum(DNA_WEIGHTS.get(base, 308.0) * count for base, count in composition.items())
        # Subtract water (18.015) for each internucleotide linkage + add end groups
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
        # Sum residue weights minus water for each peptide bond
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

def validate_and_analyze(raw_input: str) -> Dict[str, Any]:
    """
    Top-level function for Week 1 validation & statistics.
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
        # Identify invalid characters
        invalid_chars = sorted(list(set(cleaned_seq) - (DNA_CHARS | RNA_CHARS | PROTEIN_CHARS)))
        return {
            "valid": False,
            "error": f"Sequence contains unrecognized or invalid characters: {', '.join(invalid_chars)}",
            "sequence": cleaned_seq,
            "header": header,
            "type": "Invalid"
        }

    stats = calculate_stats(cleaned_seq, seq_type)

    return {
        "valid": True,
        "header": header or "User Sequence",
        "sequence": cleaned_seq,
        "type": seq_type,
        "stats": stats
    }
