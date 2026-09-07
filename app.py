"""
SeqLit Flask Web Application Entrypoint.
Orchestrates real sequence validation, NCBI BLAST search, NCBI & UniProt annotation,
PubMed retrieval with rule-based topic categorization, SQLite caching,
and report generation.
Authentic results only: if a sequence is synthetic or does not exist,
it displays zero hits and indicates that no matching biological homolog exists.
"""

import os
import ssl
import json
import logging
from flask import Flask, render_template, request, jsonify, Response

# 1. Configure robust SSL certificate validation with certifi for macOS Python
try:
    import certifi
    ssl._create_default_https_context = lambda: ssl.create_default_context(cafile=certifi.where())
except ImportError:
    pass

# 2. Ensure matplotlib writes to local writable directory
os.environ["MPLCONFIGDIR"] = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "matplotlib")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)

from modules.validation import validate_and_analyze
from modules.blast_search import execute_blast
from modules.ncbi_annotation import fetch_ncbi_annotation
from modules.uniprot_annotation import fetch_uniprot_annotation
from modules.pubmed_retrieval import fetch_pubmed_literature
from modules.cache_manager import get_cached_result, set_cached_result
from modules.report_generator import render_html_report

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("seqlit")

app = Flask(__name__)
app.config["SECRET_KEY"] = os.environ.get("SECRET_KEY", "seqlit-secret-dev-key-2026")

# Real benchmark sequence from NCBI GenBank: Human Hemoglobin Subunit Beta (HBB) mRNA
SAMPLE_DATA = {
    "name": "Human Hemoglobin Subunit Beta (HBB)",
    "type": "DNA",
    "fasta": ">NM_000518.5 Homo sapiens hemoglobin subunit beta (HBB), mRNA\nATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAA",
    "sequence": "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAA"
}

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/api/sample", methods=["GET"])
def get_sample():
    return jsonify(SAMPLE_DATA)

@app.route("/api/analyze", methods=["POST"])
def analyze_sequence():
    try:
        data = request.get_json() or {}
        raw_sequence = data.get("sequence", "").strip()

        if not raw_sequence:
            return jsonify({"success": False, "error": "No sequence content provided."}), 400

        # Step 1: Validation & Type detection
        val_result = validate_and_analyze(raw_sequence)
        if not val_result.get("valid"):
            return jsonify({"success": False, "error": val_result.get("error", "Invalid sequence input.")}), 400

        clean_seq = val_result["sequence"]
        seq_type = val_result["type"]

        # Step 2: Cache check
        cached = get_cached_result(clean_seq)
        if cached:
            logger.info("Found existing analysis in SQLite cache.")
            return jsonify({"success": True, "cached": True, "data": cached})

        # Step 3: Real BLAST Similarity search
        blast_res = execute_blast(clean_seq, seq_type)
        hits = blast_res.get("hits", [])
        
        # Only query downstream annotations if real hits were found
        if hits:
            top_hit = hits[0]
            top_accession = top_hit.get("accession", "")
            ncbi_res = fetch_ncbi_annotation(top_accession, seq_type)
            gene_term = ncbi_res.get("gene") or top_accession
            uniprot_res = fetch_uniprot_annotation(gene_term)
            organism = ncbi_res.get("organism") or ""
            literature = fetch_pubmed_literature(gene_term, organism)
        else:
            # Sequence has no homologs in biological databases
            ncbi_res = {
                "found": False,
                "accession": "None",
                "organism": "None detected",
                "gene": "None",
                "taxonomy": "No matching biological record",
                "title": "No sequence homologs found in NCBI databases."
            }
            uniprot_res = {
                "found": False,
                "protein_name": "No matching protein",
                "function": "No candidate homolog identified in databases for functional profiling.",
                "go_terms": {"biological_process": [], "molecular_function": [], "cellular_component": []},
                "features": [],
                "entry_url": ""
            }
            literature = []

        full_analysis = {
            "validation": val_result,
            "blast": blast_res,
            "ncbi": ncbi_res,
            "uniprot": uniprot_res,
            "literature": literature
        }

        # Cache only if search ran successfully
        if blast_res.get("success"):
            set_cached_result(clean_seq, seq_type, full_analysis)

        return jsonify({"success": True, "cached": False, "data": full_analysis})

    except Exception as e:
        logger.exception("Pipeline execution error")
        return jsonify({"success": False, "error": f"Server pipeline error: {str(e)}"}), 500

@app.route("/api/export/report", methods=["POST"])
def export_report():
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
    port = int(os.environ.get("PORT", 5050))
    print(f"Starting SeqLit bioinformatics server on http://127.0.0.1:{port}")
    app.run(host="0.0.0.0", port=port, debug=True)
