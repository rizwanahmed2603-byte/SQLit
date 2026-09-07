# SeqLit — Sequence Analysis & Literature Mining Dashboard

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> **SeqLit** is an end-to-end bioinformatics sequence analysis web application designed for analyzing raw DNA, RNA, and protein sequences. It performs sequence validation and statistics calculation, candidate identification via BLAST similarity search, taxonomic & gene annotation through NCBI Entrez, functional profiling with UniProt KB, and automated literature mining with topic categorization from PubMed.

---

## 🔬 Key Features

1. **Input & Validation (Objective 1 & 2):**
   - Direct raw character paste or multi-line FASTA format parsing.
   - Sequence type detection (DNA, RNA, Protein, or Invalid).
   - Real-time sequence metrics: sequence length, GC% calculation, estimated molecular weight (Da), and base/residue frequency distributions.

2. **Homology & Candidate Identification (Objective 3 & 4):**
   - NCBI BLAST similarity searching (`blastn` for nucleotide, `blastp` for protein).
   - Candidate ranking by E-value, sequence identity percentage, and alignment coverage.
   - Built-in resilience layer with fallback benchmarks for rapid testing and network rate limits.

3. **Multi-Source Functional Annotation (Objective 5 & 6):**
   - **NCBI GenBank:** Organism name, official gene symbol, and full taxonomic lineage.
   - **UniProt KB:** Recommended protein nomenclature, detailed functional summary descriptions, curated Gene Ontology (GO) terms (Biological Process, Molecular Function, Cellular Component), and structural domain annotations.

4. **Literature Mining & Categorization (Objective 7 & 8):**
   - Automated querying of PubMed via NCBI E-utilities for candidate homologs.
   - De-duplication by PMID.
   - Rule-based topic categorization:
     - 🔬 **Structure / Crystallography**
     - 💊 **Drug / Pharmacology**
     - 🧬 **Disease / Pathology**
     - ⚙️ **Function / Mechanism**
     - 🤝 **Interaction / Complex**

5. **Interactive Dashboard & Reporting (Objective 9 & 10):**
   - Responsive UI powered by Tabler UI kit and Lucide icons.
   - Interactive visualizations with Chart.js:
     - Sequence composition bar chart
     - Literature research topic distribution doughnut chart
   - Standalone styled HTML / PDF summary report generation.
   - Persistent SQLite caching layer (`seqlit_cache.db`) to avoid redundant API requests.

---

## 📁 Repository Structure

```
SQLit/
├── app.py                     # Flask entry point and pipeline API routing
├── requirements.txt           # Python dependencies
├── README.md                  # Project documentation
├── modules/
│   ├── validation.py          # Sequence cleaning, type detection, stats calculation
│   ├── blast_search.py        # BLAST submission & hit parsing with resilience fallback
│   ├── ncbi_annotation.py     # Entrez efetch for organism, gene, taxonomy
│   ├── uniprot_annotation.py  # UniProt REST API functional annotation & GO terms
│   ├── pubmed_retrieval.py    # PubMed search & rule-based topic categorization
│   ├── cache_manager.py       # SQLite caching layer for sequence queries
│   └── report_generator.py    # HTML / PDF report generation
├── templates/
│   ├── base.html              # Base layout with Tabler UI & Lucide icons
│   ├── index.html             # Sequence submission form & interactive dashboard
│   └── report_template.html   # Clean printable report layout
├── static/
│   ├── css/
│   │   └── custom.css         # Styling enhancements for dashboard & charts
│   └── js/
│       └── app.js             # Form handling, AJAX analysis trigger, Chart.js rendering
└── tests/
    ├── test_validation.py     # Unit tests for sequence validation & stats
    └── test_pipeline.py       # Unit tests for caching, categorization, & BLAST
```

---

## 🚀 Quickstart Guide

### 1. Prerequisites
- macOS, Linux, or Windows
- Python 3.11+ installed

### 2. Installation
```bash
# Clone the repository
git clone https://github.com/rizwanahmed2603-byte/SQLit.git
cd SQLit

# Create and activate virtual environment
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### 3. Run the Application
```bash
python3 app.py
```
Open your browser and navigate to:
```
http://127.0.0.1:5000
```

### 4. Running the Tests
```bash
python3 -m unittest discover tests/
```

---

## 🧪 Testing with Known Sample Sequences

Click the **"Load Example"** button on the dashboard to test the Human Hemoglobin Subunit Beta (*HBB*) benchmark sequence, or test with:

**Human Insulin (*INS*) Protein:**
```fasta
>sp|P01308|INS_HUMAN Insulin OS=Homo sapiens OX=9606 GN=INS PE=1 SV=1
MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED
LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN
```

---

## 🛡️ Resilience & Fallback Architecture

1. **SQLite Cache First:** Queries are hashed (MD5) and checked locally to prevent repeat API calls and avoid NCBI rate limits.
2. **Exponential Backoff:** Tenacity retry decorators wrap network requests to handle transient connection drops.
3. **Graceful Fallback:** If remote NCBI BLAST is rate-limited or offline, the app displays curated reference hits so that the UI, downstream annotations, and literature retrieval continue to function without interruption.
