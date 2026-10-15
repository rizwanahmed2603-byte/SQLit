# SeqLit — Sequence Analysis & Literature Mining Dashboard

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Tests: Passing](https://img.shields.io/badge/tests-23%20passed-brightgreen.svg)]()

> **SeqLit** is an end-to-end bioinformatics sequence analysis web application designed for analyzing raw DNA, RNA, and protein sequences. Built in alignment with the **SeqLit Product Requirements Document (PRD)** and **Project Execution Guide**, it replaces a fragmented manual workflow spanning NCBI BLAST, NCBI Entrez, UniProt, and PubMed with an integrated, single-pane dashboard and structured report generator.

---

## 🔬 Core Pipeline Architecture

```
Raw sequence (pasted or FASTA)
        │
        ▼
Validation & type detection (DNA / RNA / protein)
        │
        ▼
Basic stats (length, GC%, composition, molecular weight)
        │
        ▼
BLAST similarity search (remote NCBI BLAST + local fallback aligner)
        │
        ▼
Intelligent homolog ranking & biological entity resolution
        │
        ▼
NCBI GenBank annotation (Entrez efetch: organism, gene, product, taxonomy)
        │
        ▼
UniProt functional annotation (REST API: function, GO terms, domains, features)
        │
        ▼
PubMed literature retrieval (Entrez esearch/efetch + 5-topic categorization)
        │
        ▼
Interactive Dashboard (Tabler UI + Chart.js + Dynamic Hit Inspection)
        │
        ▼
Downloadable HTML / PDF Report (Jinja2 + embedded Matplotlib charts)
```

---

## 🌟 Key Features & Implementation Details

1. **Input & Validation (PRD FR-1, FR-2, FR-3):**
   - Direct raw character paste or multi-line FASTA format parsing.
   - Automatic sequence type detection (`DNA`, `RNA`, `Protein`, or `Invalid`).
   - Sequence metrics: length, GC% calculation, estimated molecular weight (Da), and composition frequency distributions.

2. **Homology & Candidate Identification (PRD FR-4, FR-Resilience-1, FR-Resilience-2):**
   - NCBI BLAST similarity searching (`blastn` for nucleotide, `blastp` for protein).
   - **Intelligent Homolog Ranking:** Prioritizes curated RefSeq (`NM_`, `NP_`) and Swiss-Prot records over automated computational predictions (`XM_`, `XP_`) and synthetic cloning vectors when identity and coverage are equal/near-equal.
   - **Local Aligner Fallback:** Built-in `Bio.Align.PairwiseAligner` against curated benchmark sequences (HBB, EBNA-1, INS, TP53) when remote NCBI BLAST is slow or unreachable.
   - **Authentic Results:** If a sequence does not exist or is synthetic nonsense, it accurately returns zero hits rather than fabricating fake homologs.

3. **Multi-Source Functional Annotation (PRD FR-5, FR-6):**
   - **NCBI GenBank:** Organism name, official gene symbol, biological product description, and taxonomic lineage.
   - **UniProt KB:** Direct Swiss-Prot accession resolution, recommended protein name, functional summary descriptions, curated Gene Ontology (GO) terms (Biological Process, Molecular Function, Cellular Component), and structural domain annotations.

4. **Literature Mining & 5-Topic Categorization (PRD FR-7, FR-8):**
   - Automated PubMed querying for confirmed homologs.
   - De-duplication by PMID.
   - Rule-based topic categorization across the 5 domains specified in PRD Section 6:
     - 🔬 **Structure / Crystallography**
     - 💊 **Drug / Pharmacology**
     - 🧬 **Disease / Pathology**
     - ⚙️ **Function / Mechanism**
     - 🤝 **Interaction / Complex**

5. **Interactive Dashboard & Reporting (PRD FR-9, FR-10, FR-11):**
   - Responsive UI powered by **Tabler UI** kit and **Lucide** icons.
   - **Interactive Candidate Hit Inspection:** Click "Inspect" on ANY hit in the BLAST table to dynamically update the NCBI, UniProt, and PubMed panels for that candidate homolog!
   - Visualizations with **Chart.js**:
     - Sequence composition bar chart
     - Literature research topic distribution doughnut chart
   - Standalone styled HTML report generation with embedded charts.
   - Persistent SQLite caching layer (`data/seqlit_cache.db`) for sub-second repeat lookups.

---

## 📁 Repository Structure

```
SQLit/
├── app.py                     # Flask entry point and pipeline API routing
├── requirements.txt           # Python dependencies
├── run.sh                     # Launch script with port clearing & SSL certifi bundle
├── README.md                  # Comprehensive project documentation
├── modules/
│   ├── validation.py          # Sequence cleaning, type detection, stats calculation
│   ├── blast_search.py        # BLAST search, intelligent ranking & local aligner fallback
│   ├── ncbi_annotation.py     # Entrez efetch for organism, gene, product, taxonomy
│   ├── uniprot_annotation.py  # UniProt REST API functional annotation & GO terms
│   ├── pubmed_retrieval.py    # PubMed search & 5-topic rule-based categorization
│   ├── cache_manager.py       # SQLite caching layer for sequence & hit queries
│   └── report_generator.py    # Jinja2 + Matplotlib HTML/PDF report generation
├── templates/
│   ├── base.html              # Base layout with Tabler UI & Lucide icons
│   ├── index.html             # Sequence submission form & interactive dashboard
│   └── report_template.html   # Clean printable report layout
├── static/
│   ├── css/
│   │   └── custom.css         # Styling enhancements for dashboard & charts
│   └── js/
│       └── app.js             # Form handling, AJAX pipeline, hit inspection & Chart.js
├── data/
│   └── seqlit_cache.db        # SQLite persistent cache database
└── tests/
    ├── test_validation.py     # Unit tests for sequence validation & stats
    ├── test_pipeline.py       # Unit tests for ranking, metadata, & local alignment
    ├── test_nonexistent_sequence.py # Authentic zero-hits tests for non-existent sequences
    └── test_app.py            # End-to-end Flask API and web route integration tests
```

---

## 🚀 Quickstart Guide

### 1. Installation
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

### 2. Run the Application
Use the provided `run.sh` script (which automatically clears stale ports and sets the macOS Python SSL certificate bundle):
```bash
chmod +x run.sh
./run.sh
```
Or start directly with Python:
```bash
python3 app.py
```
Open your browser and navigate to:
```
http://127.0.0.1:5050
```

### 3. Running the Test Suite
Run the full test suite with `pytest`:
```bash
python3 -m pytest -v
```

---

## 🧪 Benchmark Test Sequences

The application includes built-in benchmark sequences selectable via the **"Load Benchmark Example"** dropdown in the UI:

1. **Human Hemoglobin Subunit Beta (*HBB*) — DNA (444 bp):**
   - Canonical beta-globin mRNA coding sequence (`NM_000518.5`)
   - Verified top hit: *Homo sapiens* HBB (100% identity, 100% coverage)
   - UniProt: Hemoglobin subunit beta (`P68871`), oxygen transport, heme binding, hemoglobin complex
   - Literature: Articles on sickle cell disease, beta-thalassemia, hemoglobin crystallography, and allosteric modulators

2. **Epstein-Barr Nuclear Antigen 1 (*EBNA-1*) — Protein (641 aa):**
   - Essential viral sequence-specific DNA-binding protein (`P03211`)
   - Verified top hit: *Human gammaherpesvirus 4* EBNA-1
   - UniProt: EBNA-1 (`P03211`), episome maintenance, DNA replication
   - Literature: Articles on crystal structure, targeted small-molecule inhibitors, and pathology

3. **Human Insulin (*INS*) — DNA (333 bp):**
   - Preproinsulin coding sequence (`NM_000207.3`)
   - UniProt: Insulin preproprotein (`P01308`)

4. **Human Tumor Protein P53 (*TP53*) — DNA (1182 bp):**
   - Cellular tumor suppressor p53 (`NM_000546.6`)
   - UniProt: Cellular tumor antigen p53 (`P04637`)
