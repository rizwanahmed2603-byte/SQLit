"""
SeqLit Flask Web Application Entrypoint.
Orchestrates sequence validation, BLAST search, NCBI & UniProt annotation,
PubMed retrieval with rule-based topic categorization, SQLite caching,
and report generation.
"""

import os
import json
import logging
from flask import Flask, render_template, request, jsonify, Response, send_file

from modules.validation import validate_and_analyze
from modules.blast_search import execute_blast
from modules.ncbi_annotation import fetch_ncbi_annotation
from modules.uniprot_annotation import fetch_uniprot_annotation
from modules.pubmed_retrieval import fetch_pubmed_literature
from modules.cache_manager import get_cached_result, set_cached_result
from modules.report_generator import render_html_report, export_pdf_report

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("seqlit")

app = Flask(__name__)
app.config["SECRET_KEY"] = os.environ.get("SECRET_KEY", "seqlit-secret-dev-key-2026")

# Curated sample sequence (Human HBB)
SAMPLE_DATA = {
    "name": "Human Hemoglobin Subunit Beta (HBB)",
    "type": "DNA",
    "fasta": ">NM_000518.5 Homo sapiens hemoglobin subunit beta (HBB), mRNA\nATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAA",
    "sequence": "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAA"
}

@app.route("/")
def index():
    """Renders the main single-page submission and dashboard view."""
    return render_template("index.html")

@app.route("/api/sample", methods=["GET"])
def get_sample():
    """Provides a known benchmark sample sequence."""
    return jsonify(SAMPLE_DATA)

@app.route("/api/analyze", methods=["POST"])
def analyze_sequence():
    """
    Main pipeline API:
    1. Validation & Stats (Week 1)
    2. Cache check (Week 8 optimization)
    3. BLAST Similarity search (Week 2)
    4. NCBI Entrez annotation (Week 3)
    5. UniProt annotation (Week 3)
    6. PubMed retrieval & topic categorization (Week 4)
    """
    try:
        data = request.get_json() or {}
        raw_sequence = data.get("sequence", "").strip()

        if not raw_sequence:
            return jsonify({"success": False, "error": "No sequence content provided."}), 400

        # Step 1: Validation
        val_result = validate_and_analyze(raw_sequence)
        if not val_result.get("valid"):
            return jsonify({"success": False, "error": val_result.get("error", "Invalid sequence input.")}), 400

        clean_seq = val_result["sequence"]
        seq_type = val_result["type"]

        # Step 2: Cache check
        cached = get_cached_result(clean_seq)
        if cached:
            logger.info("Found existing analysis in SQLite cache. Returning cached result.")
            return jsonify({"success": True, "cached": True, "data": cached})

        # Step 3: BLAST Similarity search
        blast_res = execute_blast(clean_seq, seq_type)
        top_hit = blast_res["hits"][0] if blast_res.get("hits") else {}
        top_accession = top_hit.get("accession", "")

        # Step 4: NCBI Annotation
        ncbi_res = fetch_ncbi_annotation(top_accession, seq_type)

        # Step 5: UniProt Functional Annotation
        gene_term = ncbi_res.get("gene") or top_accession
        uniprot_res = fetch_uniprot_annotation(gene_term)

        # Step 6: PubMed Literature
        organism = ncbi_res.get("organism") or "Homo sapiens"
        literature = fetch_pubmed_literature(gene_term, organism)

        full_analysis = {
            "validation": val_result,
            "blast": blast_res,
            "ncbi": ncbi_res,
            "uniprot": uniprot_res,
            "literature": literature
        }

        # Cache the completed result
        set_cached_result(clean_seq, seq_type, full_analysis)

        return jsonify({"success": True, "cached": False, "data": full_analysis})

    except Exception as e:
        logger.exception("Pipeline execution error")
        return jsonify({"success": False, "error": f"Server pipeline error: {str(e)}"}), 500

@app.route("/api/export/report", methods=["POST"])
def export_report():
    """Generates downloadable standalone HTML report."""
    try:
        payload = request.get_json() or {}
        analysis_data = payload.get("data")
        if not analysis_data:
            return jsonify({"error": "No data provided for report generation"}), 400

        template_path = os.path.join(app.root_path, "templates", "report_template.html")
        with open(template_path, "r", encoding="utf-8") as f:
            template_str = f.read()

        html_out = render_html_report(analysis_data, template_str)
        return Response(
            html_out,
            mimetype="text/html",
            headers={"Content-Disposition": "attachment;filename=seqlit_report.html"}
        )
    except Exception as e:
        logger.exception("Report export failed")
        return jsonify({"error": str(e)}), 500

if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))
    print(f"Starting SeqLit bioinformatics server on http://127.0.0.1:{port}")
    app.run(host="0.0.0.0", port=port, debug=True)
