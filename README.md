# SeqLit — Sequence Analysis & Literature Mining Dashboard

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Tests: Passing](https://img.shields.io/badge/tests-60%20passed-brightgreen.svg)]()

> **SeqLit** is an end-to-end bioinformatics sequence analysis web application designed for analyzing raw DNA, RNA, and protein sequences. Built in alignment with the **SeqLit Product Requirements Document (PRD)** and **Project Execution Guide**, it replaces a fragmented manual workflow spanning NCBI BLAST, NCBI Entrez, UniProt, and PubMed with an integrated, single-pane dashboard and structured report generator.

---

## 🔬 Core Pipeline Architecture

```
Raw sequence (pasted, FASTA, multi-database flat files, 3-letter codes, or gapped alignments)
        │
        ▼
Multi-Format Parser & Cleaner (GenBank ORIGIN, UniProt SQ, 3-letter to 1-letter, gap stripping)
        │
        ▼
Validation & Intelligent Type Detection (DNA / RNA / Protein / Degenerate IUPAC)
        │
        ▼
Basic Metrics (length, GC%, composition frequency, molecular weight)
        │
        ▼
Fast Homolog Search:
 ├─ Sub-millisecond k-mer pre-screening
 ├─ Local Smith-Waterman aligner (forward & reverse-complement)
 ├─ 6-Frame Translation (BLASTX mode for nucleotide queries against protein references)
 ├─ FASTA/PDB/Ensembl Header Reference Resolver
 └─ Remote NCBI BLAST with non-blocking timeout handling
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

1. **Universal Multi-Database Input & Validation (PRD FR-1, FR-2, FR-3):**
   - **NCBI GenBank & DDBJ Flat Files:** Directly paste full records; automatically parses sequences between `ORIGIN` and `//`, stripping line numbers and extracting metadata.
   - **UniProtKB & EMBL Flat Files:** Parses `SQ   SEQUENCE ... //` blocks, extracting accession, gene, product, and organism.
   - **3-Letter Amino Acid Peptides:** Automatically parses and translates 3-letter code strings (`Met-Val-His-Leu-Thr...` or `H-Ala-Cys-Asp-OH` or `MET VAL HIS...`) into 1-letter standard amino acids.
   - **Alignment & Gapped Sequences:** Automatically strips alignment gap characters (`-`, `.`, `~`), spaces, and terminal stop codons (`*`).
   - **PDB & Ensembl Header Resolver:** Extracts PDB IDs (`4HHB`, `1A3N`) and chain names while ignoring numeric taxonomy tags like `(9606)`, and parses Ensembl transcript/gene IDs (`ENST...`, `ENSG...`).
   - **Intelligent Type Detection:** Automatically differentiates standard and degenerate `DNA`, `RNA`, and `Protein` with frequency-based amino acid disambiguation.
   - **Physicochemical Metrics:** Calculates sequence length, GC% ratio, molecular weight (Da), and base/residue composition.

2. **Sub-Millisecond Search & 6-Frame Translation (PRD FR-4, FR-Resilience-1, FR-Resilience-2):**
   - **k-mer Pre-Screening:** Sub-millisecond candidate lookup ($<0.01\text{ s}$) before Smith-Waterman alignment.
   - **6-Frame Translation (BLASTX Mode):** Translates nucleotide queries in all 6 frames ($+1, +2, +3, -1, -2, -3$) to identify protein homolog counterparts.
   - **Reverse Complement Strand Detection:** Automatically tests forward and reverse-complement strands for all nucleotide sequences.
   - **Comprehensive Reference Library:** Over 38 canonical reference targets covering:
     - *Hemoglobins:* `HBB`, `HBA1`, `HBD`, `HBG1`, `HBG2`
     - *Immune Checkpoints & Receptors:* `CD4`, `CD8A`, `PDCD1` (PD-1), `CD274` (PD-L1), `CTLA4`, `B2M`
     - *Clinical & Cancer Markers:* `TP53`, `KRAS`, `BRAF`, `PTEN`, `EGFR`, `BRCA1`, `BRCA2`, `MYC`, `CFTR`, `LDLR`, `APOE`, `TNF`, `IL6`, `ACE2`
     - *Model Organisms & Laboratory Enzymes:* `Cas9` (SpCas9), `GFP`, `EGFP`, `LYZ` (Hen HEWL & Human), `MB` (Myoglobin), `ALB` (Bovine BSA & Human), `CYCS`, `GAPDH`, `ACTB`, `EBNA1`, `SARS2_S`, `SARS2_N`, `16S rRNA`
   - **Remote NCBI BLAST:** Falls back to remote NCBI `blastn` and `blastp` with non-blocking timeout handling.
   - **Authentic Zero-Hit Guarantee:** Non-existent or synthetic random sequences authentically return zero hits with zero fabrication.

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
   - **Interactive Candidate Hit Inspection:** Click "Inspect" on ANY hit in the BLAST table to dynamically update the NCBI, UniProt, and PubMed panels for that candidate homolog.
   - Visualizations with **Chart.js**:
     - Sequence composition frequency chart
     - Literature research topic distribution doughnut chart
   - Standalone styled HTML report generation with embedded charts.
   - Persistent SQLite caching layer (`data/seqlit_cache.db`) for sub-second repeat lookups.

---

## 📚 Project Documentation & References

- 📋 **[Product Requirements Document (PRD)](./SeqLit_PRD.md):** Complete specifications, functional requirements (FR-1 to FR-11), acceptance criteria, and system architecture.
- 🛠️ **[Project Execution Guide](./SeqLit_Project_Execution_Guide.md):** Week-by-week build plan, tech stack analysis, macOS setup instructions, and fallback layers.
- 🧪 **[Test Log & Ground Truth Verification](./TEST_LOG.md):** Comprehensive test report covering benchmark test sequences (HBB, EBNA-1, INS, TP53), edge cases, and bug fix logs.

---

## 📁 Repository Structure

```
SQLit/
├── app.py                     # Flask entry point and pipeline API routing
├── requirements.txt           # Python dependencies
├── run.sh                     # Launch script with port clearing & SSL certifi bundle
├── SeqLit_PRD.md              # Product Requirements Document
├── SeqLit_Project_Execution_Guide.md # Execution Guide & architecture roadmap
├── TEST_LOG.md                # Comprehensive testing & validation log
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

---

## 🚀 Free Web Deployment (Live Website with Free Domain)

SeqLit is pre-configured with a [`Procfile`](Procfile), [`render.yaml`](render.yaml), and a multi-platform [`Dockerfile`](Dockerfile) for instant 1-click cloud deployment with a free public domain and automatic HTTPS:

| Platform | Free Domain Example | Auto-Deploy on Push | Free Tier Specs | Setup Time |
|---|---|---|---|---|
| **[Render](https://render.com)** *(Recommended)* | `https://seqlit.onrender.com` | ✅ Yes (GitHub) | Free Web Service | 2 minutes |
| **[Hugging Face Spaces](https://huggingface.co/spaces)** | `https://<user>-seqlit.hf.space` | ✅ Yes (Git / Docker) | 2 vCPU, 16 GB RAM | 2 minutes |
| **[Koyeb](https://koyeb.com)** | `https://seqlit.koyeb.app` | ✅ Yes (GitHub) | Free Eco nano instance | 2 minutes |

👉 See the complete step-by-step walkthrough in [**`DEPLOYMENT.md`**](DEPLOYMENT.md).

---

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
