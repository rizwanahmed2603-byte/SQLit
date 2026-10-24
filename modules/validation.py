"""
Sequence validation and statistics module for SeqLit.
Handles raw sequence, FASTA, multi-database flat files (GenBank, UniProt, EMBL, DDBJ),
3-letter amino acid peptide notations, alignment outputs with gaps,
sequence type detection (DNA, RNA, Protein), degenerate IUPAC alphabet support,
and computation of basic metrics (GC%, molecular weight, composition frequency).
"""

import re
from typing import Dict, Any, Tuple, Optional

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

# 3-Letter to 1-Letter Amino Acid Dictionary
THREE_TO_ONE_AA = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'SEC': 'U', 'PYL': 'O', 'ASX': 'B', 'GLX': 'Z', 'XAA': 'X'
}

def parse_three_letter_peptide(raw_text: str) -> Optional[str]:
    """
    Detects and converts 3-letter amino acid sequence notation
    (e.g., Met-Val-His-Leu-Thr... or H-Ala-Cys-Asp-OH or MET VAL HIS...)
    into standard 1-letter protein sequence.
    """
    # Strip terminal modifications if present
    stripped = re.sub(r'^(?:H-|Ac-)|(?:-OH|-NH2)$', '', raw_text.strip(), flags=re.IGNORECASE)
    # Match 3-letter word tokens
    tokens = re.findall(r'\b[A-Za-z]{3}\b', stripped)
    if len(tokens) >= 3:
        matched = [THREE_TO_ONE_AA.get(t.upper()) for t in tokens if t.upper() in THREE_TO_ONE_AA]
        # If at least 75% of 3-letter tokens are valid amino acids, convert
        if len(matched) / len(tokens) >= 0.75:
            return "".join(THREE_TO_ONE_AA.get(t.upper(), "X") for t in tokens)
    return None

def extract_flat_file_sequence(raw_text: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Extracts sequence and header metadata from database flat file formats:
    - NCBI / DDBJ GenBank flat file (ORIGIN ... //)
    - UniProtKB / Swiss-Prot flat file (SQ   SEQUENCE ... //)
    - EMBL / ENA flat file (SQ   Sequence ... //)
    """
    # 1. GenBank / DDBJ format: sequence begins after ORIGIN line and ends at //
    origin_match = re.search(r'\bORIGIN\b([\s\S]*?)(?://|$)', raw_text, re.IGNORECASE)
    if origin_match:
        seq_chunk = origin_match.group(1)
        cleaned = re.sub(r'[\s\d/]', '', seq_chunk).upper()
        # Extract metadata for header
        acc_m = re.search(r'\bACCESSION\s+([A-Za-z0-9_]+)', raw_text)
        def_m = re.search(r'\bDEFINITION\s+([\s\S]+?)(?:ACCESSION|VERSION|KEYWORDS|SOURCE)', raw_text)
        org_m = re.search(r'\bORGANISM\s+([^\n]+)', raw_text)
        acc = acc_m.group(1) if acc_m else ""
        desc = re.sub(r'\s+', ' ', def_m.group(1).strip()) if def_m else ""
        org = org_m.group(1).strip() if org_m else ""
        header = f"{acc} {desc}".strip()
        gene_m = re.search(r'/(?:gene|locus_tag)="([^"]+)"', raw_text)
        if gene_m:
            header += f" gene={gene_m.group(1).strip()}"
        if org and org not in header:
            header += f" [{org}]"
        return cleaned, header or "GenBank Record"

    # 2. UniProt / EMBL format: sequence begins after SQ line and ends at //
    sq_match = re.search(r'\bSQ\s+(?:SEQUENCE|Sequence)[\s\S]*?\n([\s\S]*?)(?://|$)', raw_text, re.IGNORECASE)
    if sq_match:
        seq_chunk = sq_match.group(1)
        cleaned = re.sub(r'[\s\d/]', '', seq_chunk).upper()
        # Extract metadata
        ac_m = re.search(r'\bAC\s+([A-Z0-9]+);?', raw_text)
        de_m = re.search(r'\bDE\s+RecName:\s*Full=([^;\n]+)', raw_text)
        gn_m = re.search(r'\bGN\s+Name=([^;\s]+)', raw_text)
        os_m = re.search(r'\bOS\s+([^.\n]+)', raw_text)
        acc = ac_m.group(1) if ac_m else ""
        de = de_m.group(1).strip() if de_m else ""
        gn = gn_m.group(1).strip() if gn_m else ""
        os = os_m.group(1).strip() if os_m else ""
        header_parts = [acc]
        if gn:
            header_parts.append(f"GN={gn}")
        if de:
            header_parts.append(de)
        header = " ".join(header_parts)
        if os:
            header += f" [{os}]"
        return cleaned, header or "UniProt Record"

    return None, None

def clean_sequence(raw_input: str) -> Tuple[str, str]:
    """
    Robust multi-database sequence cleaner:
    - Parses GenBank, UniProt, EMBL, and DDBJ flat file records directly
    - Converts 3-letter amino acid notations (Met-Val-His... -> MVH...)
    - Strips FASTA headers (> or ;)
    - Removes line numbers, alignment gaps (-, ., ~), whitespace, and trailing stop codons (*).
    Returns (cleaned_sequence_upper, header_or_empty).
    """
    if not raw_input or not raw_input.strip():
        return "", ""

    text = raw_input.strip()

    # Check 1: Database flat file formats (GenBank ORIGIN, UniProt SQ)
    flat_seq, flat_header = extract_flat_file_sequence(text)
    if flat_seq:
        flat_seq = flat_seq.rstrip("*")
        return flat_seq, flat_header or ""

    # Check 2: 3-letter amino acid peptide notation (e.g. Met-Val-His-Leu-Thr...)
    three_aa_seq = parse_three_letter_peptide(text)
    if three_aa_seq:
        return three_aa_seq, "Synthetic / Natural Peptide (3-Letter Converted)"

    # Check 3: Standard FASTA or raw sequence lines
    lines = text.splitlines()
    header = ""
    seq_lines = []

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">") or stripped.startswith(";"):
            if not header:
                header = stripped.lstrip(">;").strip()
        else:
            # Remove whitespace, line numbers, alignment gaps (-, ., ~), and slashes
            cleaned_line = re.sub(r"[\s\d\-\.\~\/]", "", stripped)
            seq_lines.append(cleaned_line)

    sequence = "".join(seq_lines).upper()
    # Strip terminal stop codon if present e.g. from translation output
    sequence = sequence.rstrip("*")
    return sequence, header

def detect_sequence_type(seq: str) -> str:
    """
    Detects whether the sequence is DNA, RNA, Protein, or Invalid.
    Supports standard and degenerate IUPAC codes with intelligent
    amino acid vs. nucleotide disambiguation.
    """
    if not seq:
        return "Unknown"

    unique_chars = set(seq)

    # 1. Any characters exclusive to proteins (E, F, I, L, P, Q, Z, J, O)
    protein_only_chars = set("EFILPQZJO")
    if unique_chars & protein_only_chars:
        if unique_chars.issubset(PROTEIN_EXTENDED):
            return "Protein"
        return "Invalid"

    # 2. RNA detection: 'U' present and no 'T'
    if 'U' in unique_chars and 'T' not in unique_chars:
        if unique_chars.issubset(RNA_CHARS) or unique_chars.issubset(RNA_DEGENERATE):
            return "RNA"
        if unique_chars.issubset(PROTEIN_EXTENDED):
            return "Protein"
        return "Invalid"

    # 3. DNA detection: ATCGN or degenerate DNA
    if unique_chars.issubset(DNA_CHARS) or unique_chars.issubset(DNA_DEGENERATE):
        return "DNA"

    # 4. If both T and U are present (e.g. Threonine + Selenocysteine in protein)
    if 'T' in unique_chars and 'U' in unique_chars:
        if unique_chars.issubset(PROTEIN_EXTENDED):
            return "Protein"
        return "Invalid"

    # 5. Standard and extended Protein
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
    Extracts biological hints (accession, gene symbol, organism) from sequence headers
    across NCBI, UniProt, PDB, Ensembl, EMBL, and DDBJ formats.
    """
    hints = {"accession": "", "gene": "", "organism": ""}
    if not header:
        return hints

    # 1. Match Swiss-Prot / UniProt accession: sp|P01308| or tr|...| or standalone UniProt ID
    sp_m = re.search(r"(?:sp|tr)\|([A-Z0-9]{6,10})", header)
    if sp_m:
        hints["accession"] = sp_m.group(1)
    else:
        uniprot_stand = re.search(r"\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})\b", header)
        if uniprot_stand:
            hints["accession"] = uniprot_stand.group(1)

    # 2. Match Ensembl IDs: ENST..., ENSG..., ENSP...
    if not hints["accession"]:
        ens_m = re.search(r"\b(ENS[A-Z]{0,4}[GTP]\d{9,11})\b", header)
        if ens_m:
            hints["accession"] = ens_m.group(1)

    # 3. Match RefSeq / GenBank / PDB accession
    if not hints["accession"]:
        ref_m = re.search(r"\b([A-Z]{1,2}_\d{5,9}|[A-Z]{1,3}\d{5,8}|[A-Z]{4}\d{8})\b", header)
        if ref_m:
            hints["accession"] = ref_m.group(1)

    # 4. Match PDB Structure ID: e.g. >4HHB_1, >1A3N:A, >101M, >pdb|4HHB|A
    if not hints["accession"]:
        pdb_m = re.search(r"\b([0-9][A-Za-z0-9]{3})(?:_[A-Za-z0-9]+|:[A-Za-z0-9]+)?\b", header)
        if pdb_m:
            hints["accession"] = pdb_m.group(1).upper()

    # 5. Match Organism: OS=... or [Homo sapiens] or canonical names
    org_os = re.search(r"OS=([^=]+?)(?:\s+[A-Z]{2}=|$)", header)
    if org_os:
        hints["organism"] = org_os.group(1).strip()
    else:
        org_bracket = re.search(r"\[([A-Za-z0-9_.\s-]+)\]", header)
        if org_bracket:
            cand_org = org_bracket.group(1).strip()
            # If bracket contains Source/Acc info, ignore
            if not cand_org.startswith("Source:"):
                hints["organism"] = cand_org
        if not hints["organism"]:
            org_common = re.search(
                r"\b(Homo sapiens|Mus musculus|Rattus norvegicus|Escherichia coli|Aequorea victoria|Gallus gallus|Physeter catodon|Epstein-Barr virus|Human gammaherpesvirus 4|SARS-CoV-2|Severe acute respiratory syndrome coronavirus 2|Bos taurus|Streptococcus pyogenes)\b",
                header,
                re.IGNORECASE
            )
            if org_common:
                hints["organism"] = org_common.group(1).strip()

    # 6. Match Gene:
    # Explicit gene flags: GN=HBB, gene=HBB, gene_symbol:HBB
    gn_explicit = re.search(r"(?:GN=|gene=|gene_symbol:)\s*([A-Za-z0-9_-]{2,12})\b", header, re.IGNORECASE)
    if gn_explicit:
        hints["gene"] = gn_explicit.group(1).strip().upper()
    else:
        # Check known canonical genes first
        known_genes = re.search(
            r"\b(HBB|HBA1|HBA|HBD|HBG1|HBG2|INS|TP53|GAPDH|ACTB|ALB|EGFR|BRCA1|BRCA2|KRAS|BRAF|PTEN|MYC|CYCS|GFP|EGFP|LYZ|MB|EBNA1|APOE|TNF|IL6|ACE2|B2M|TUBB|CFTR|LDLR|CD4|CD8A|PDCD1|CD274|CTLA4|CAS9|SARS2)\b",
            header,
            re.IGNORECASE
        )
        if known_genes:
            hints["gene"] = known_genes.group(1).strip().upper()
        else:
            # Parenthetical symbols like (HBB) - but ignore taxonomy IDs like (9606) or common words like (Human)
            gn_paren = re.search(r"\(([A-Za-z][A-Za-z0-9_-]{1,11})\)", header)
            if gn_paren:
                cand = gn_paren.group(1).strip().upper()
                if not cand.isdigit() and not cand.startswith("LOC") and cand not in ("HUMAN", "MOUSE", "RAT", "CHICKEN", "YEAST", "VIRUS"):
                    hints["gene"] = cand

    # If gene still not identified, search for canonical gene names or PDB product descriptions
    if not hints["gene"]:
        header_upper = header.upper()
        if "BETA CHAIN" in header_upper or "HEMOGLOBIN SUBUNIT BETA" in header_upper:
            hints["gene"] = "HBB"
        elif "ALPHA CHAIN" in header_upper or "HEMOGLOBIN SUBUNIT ALPHA" in header_upper:
            hints["gene"] = "HBA1"
        elif "DELTA CHAIN" in header_upper or "HEMOGLOBIN SUBUNIT DELTA" in header_upper:
            hints["gene"] = "HBD"
        elif "GAMMA CHAIN" in header_upper or "HEMOGLOBIN SUBUNIT GAMMA" in header_upper:
            hints["gene"] = "HBG1"
        else:
            known_genes = re.search(
                r"\b(HBB|HBA1|HBA|HBD|HBG1|HBG2|INS|TP53|GAPDH|ACTB|ALB|EGFR|BRCA1|BRCA2|KRAS|BRAF|PTEN|MYC|CYCS|GFP|EGFP|LYZ|MB|EBNA1|APOE|TNF|IL6|ACE2|B2M|TUBB|CFTR|LDLR|CD4|CD8A|PDCD1|CD274|CTLA4|CAS9|SARS2)\b",
                header,
                re.IGNORECASE
            )
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
