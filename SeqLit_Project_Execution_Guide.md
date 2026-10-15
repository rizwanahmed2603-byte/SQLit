# SeqLit — Project Execution Guide
### A free, macOS-ready roadmap for your MSc Bioinformatics minor project

This turns your synopsis's proposed pipeline and 8-week timeline into an actual build plan: what to install, in what order, and what to build each week. Everything listed below — every tool, library, and data source — is free. Nothing here needs a subscription, license key, or credit card.

---

## 1. The Complete Tech Stack

### Core development environment
| Tool | Role | Cost |
|---|---|---|
| Homebrew | macOS package manager — installs most things below with one command | Free |
| Python 3.11+ | Language for the whole backend | Free |
| Visual Studio Code | Code editor (Python, HTML/CSS/JS, Markdown support) | Free |
| Git + GitHub | Version control and backup/portfolio hosting | Free |
| Terminal.app (built-in) | Command line — nothing extra to install | Free |

### Backend & bioinformatics (Python)
| Tool | Role | Cost |
|---|---|---|
| Flask | Web framework — serves the app and API routes | Free |
| Biopython | Sequence parsing, BLAST calls, Entrez (NCBI) access, PairwiseAligner | Free |
| requests | HTTP calls to the UniProt REST API | Free |
| Pandas | Tabular handling of hits, annotations, literature | Free |
| NCBI BLAST+ (optional, local) | Local blastn/blastp instead of the web service | Free |

### Frontend
| Tool | Role | Cost |
|---|---|---|
| HTML5 / CSS3 / vanilla JavaScript | Submission form + dashboard | Free (browser-native) |
| [Tabler](https://github.com/tabler/tabler) | Free, MIT-licensed HTML/CSS dashboard UI kit built on Bootstrap — real sidebar nav, cards, tables and form components out of the box, so the interface reads as a built product instead of a bare Bootstrap starter page | Free |
| Chart.js (via CDN) | Interactive charts on the dashboard | Free |
| [Lucide](https://github.com/lucide-icons/lucide) | Clean, consistent SVG icon set (MIT), 1,600+ icons, for panel headers, buttons and status indicators | Free |

> Tabler's own `dist/libs` folder ships ApexCharts, which is dual-licensed (free only under $2M annual org revenue). You don't need it — pull in just Tabler's CSS/HTML components and use Chart.js (already free, MIT) for charts, as planned above.

### External data sources (public, no login required to use)
| Resource | What you get | Cost |
|---|---|---|
| NCBI BLAST | Sequence similarity search | Free |
| NCBI Entrez / E-utilities | Sequence, gene, taxonomy records + PubMed search | Free (email required; free API key optional — see Section 2) |
| UniProt REST API | Protein function, domains, GO terms, cross-references | Free, no key needed |
| PubMed | Literature metadata and abstracts | Free, via Entrez |

### Reporting & visualization
| Tool | Role | Cost |
|---|---|---|
| Jinja2 (ships with Flask) | HTML report templating | Free |
| Matplotlib | Static charts baked into report | Free |

### Storage, testing, docs
| Tool | Role | Cost |
|---|---|---|
| SQLite (built into Python) | Lightweight local cache/results store — no server to install | Free |
| pytest | Unit tests for validation/parsing logic | Free |
| draw.io / diagrams.net | Architecture diagram for your report or viva | Free |

### Resilience & fallback
| Tool | Role | Cost |
|---|---|---|
| [tenacity](https://github.com/jd/tenacity) | Retry decorator with exponential backoff — wraps every external API call | Free |
| Biopython PairwiseAligner | Built-in local sequence aligner — fallback identification path against benchmark library when remote BLAST is slow or unreachable | Free |

See Section 7 for how these fit together as an actual fallback chain.

Nothing above has a paid tier that matters here — NCBI, UniProt and PubMed are publicly funded resources with no commercial access tier, and every library listed is open-source.

---

## 2. macOS Setup, Step by Step

Open Terminal (Applications → Utilities → Terminal) and work through these in order.

**1. Install Xcode Command Line Tools** (gives you `git`, `make`, a C compiler — needed by packages below)
```bash
xcode-select --install
```

**2. Install Homebrew**
```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

**3. Install Python**
```bash
brew install python
python3 --version
```

**4. Create the project folder and a virtual environment**
```bash
mkdir seqlit && cd seqlit
python3 -m venv venv
source venv/bin/activate
```
You'll run `source venv/bin/activate` every time you reopen Terminal to work on the project — it keeps this project's packages separate from the rest of your Mac.

**5. Install the core Python packages**
```bash
pip install -r requirements.txt
```

**6. Start the server**
```bash
chmod +x run.sh
./run.sh
```

---

## 3. The Pipeline, End to End

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
BLAST search (NCBI BLASTN/BLASTP + Local PairwiseAligner fallback) ──▶ ranked hits
        │
        ▼
NCBI annotation (Entrez efetch: organism, gene, product, taxonomy)
        │
        ▼
UniProt annotation (REST API: function, GO terms, domains, cross-refs)
        │
        ▼
PubMed search (Entrez esearch/efetch) ──▶ de-duplicated, 5-topic categorized publications
        │
        ▼
Dashboard (Tabler UI + tables + Chart.js visuals + interactive hit inspection)
        │
        ▼
Report generator (Jinja2 → HTML report with embedded charts)
```

---

## 4. Week-by-Week Build Plan

### Week 1 — Environment, UI skeleton, sequence validation
- Complete the setup in Section 2.
- Sketch the dashboard layout first: one input panel, four result panels (Identification / Annotation / Literature / Report).
- Build a minimal Flask app: one route (`/`) rendering the submission form, one route (`/api/analyze`) accepting a POST with the sequence.
- **Validation module:** strip whitespace/headers, confirm only valid characters appear (ATCG(N) for DNA, AUCG for RNA, the 20 standard amino-acid letters for protein), and auto-detect the sequence type from character composition.
- *Done when:* pasting a sequence into the form returns its detected type and basic stats on screen.

### Week 2 — BLAST integration
- Remote BLAST via `Bio.Blast.NCBIWWW.qblast()`.
- Parse the returned XML with `Bio.Blast.NCBIXML` and extract accession, percent identity, E-value, and query coverage per hit.
- Rank hits intelligently, prioritizing curated RefSeq over automated predictions.
- *Done when:* submitting a known test sequence (a well-characterized gene) returns a correct, ranked top hit.

### Week 3 — NCBI annotation + UniProt integration
- For the top BLAST accession(s), call `Bio.Entrez.efetch` against the nucleotide/protein database to pull organism, gene name, product, and taxonomy from the GenBank record.
- For protein annotation, query UniProt with direct Swiss-Prot accession resolution and multi-tier query fallback. Pull protein name, function description, GO terms, and cross-references.
- *Done when:* the Annotation panel shows organism, gene, function, and GO terms for a submitted sequence.

### Week 4 — PubMed literature retrieval
- Build search terms from gene symbol, organism, and biological product.
- `Bio.Entrez.esearch(db="pubmed", term=...)` to get PMIDs, then `esummary`/`efetch` to pull title, authors, journal, year, and abstract.
- De-duplicate by PMID.
- Categorize with rule-based keyword matching across the 5 PRD categories (Structure, Drug, Disease, Function, Interaction).
- *Done when:* the Literature panel lists real, correctly categorized, de-duplicated papers linking out to PubMed.

### Week 5 — Dashboard, tables, filters, visualizations
- Wire all four result panels to real backend data.
- Add sortable/filterable tables.
- Add Chart.js visuals: composition bar chart, literature topic doughnut chart.
- Add interactive candidate hit inspection: clicking "Inspect" on any BLAST row updates all panels for that candidate homolog!
- *Done when:* a full run — pasted sequence to populated dashboard — works end to end in the browser.

### Week 6 — Report generation
- Build an HTML report template with Jinja2, reusing the same data the dashboard displays.
- Generate static chart images with Matplotlib for report inclusion.
- *Done when:* clicking "Download report" produces a clean, complete report.

### Week 7 — Testing
- Test with a handful of well-characterized sequences (HBB, EBNA-1, INS, TP53).
- Test edge cases deliberately: empty input, garbage characters, short sequences, ambiguous bases, synthetic sequences with no hits.
- Write unit tests for validation, ranking, caching, and API routes.
- Document in `TEST_LOG.md`.
- *Done when:* you have a written test log with all tests passing.

### Week 8 — Optimization, documentation, presentation
- Add SQLite cache (`data/seqlit_cache.db`) keyed by sequence MD5 hash.
- Wrap external API calls in tenacity retries with exponential backoff.
- Write the README and Project Execution Guide.
- *Done when:* the repo is documented well enough that someone else could clone it and run it from the README alone.

---

## 5. Suggested Project Structure

```
SQLit/
├── app.py                     # Flask entry point and API routes
├── requirements.txt           # Python dependencies
├── run.sh                     # Launch script
├── SeqLit_PRD.md              # Product Requirements Document
├── SeqLit_Project_Execution_Guide.md # Execution Guide & architecture roadmap
├── TEST_LOG.md                # Comprehensive testing & validation log
├── README.md                  # Main documentation
├── modules/
│   ├── validation.py          # sequence cleaning, type detection, stats
│   ├── blast_search.py        # BLAST submission, ranking & local aligner fallback
│   ├── ncbi_annotation.py     # Entrez efetch for organism, gene, product, taxonomy
│   ├── uniprot_annotation.py  # UniProt REST calls, Swiss-Prot lookups, GO terms
│   ├── pubmed_retrieval.py    # Entrez esearch/esummary & 5-topic categorization
│   ├── cache_manager.py       # SQLite caching layer for sequence & hit queries
│   └── report_generator.py    # Jinja2 + Matplotlib HTML report generation
├── templates/
│   ├── base.html              # Base layout with Tabler UI & Lucide icons
│   ├── index.html             # submission form, dashboard, hit inspector
│   └── report_template.html   # printable report layout
├── static/
│   ├── css/custom.css         # Styling enhancements
│   └── js/app.js              # Dashboard client logic, Chart.js & hit inspection
├── data/
│   └── seqlit_cache.db        # SQLite database
└── tests/
    ├── test_validation.py     # Unit tests for sequence validation & stats
    ├── test_pipeline.py       # Unit tests for ranking, metadata, & local alignment
    ├── test_nonexistent_sequence.py # Tests verifying authentic zero-hits
    └── test_app.py            # End-to-end Flask API and web route tests
```

---

## 6. Practical Notes Worth Knowing Upfront

- **NCBI etiquette:** always set `Entrez.email`; keep calls under 3/second (`time.sleep(0.34)` is implemented).
- **Remote vs. local BLAST:** remote (`qblast`) is primary. Biopython's built-in `PairwiseAligner` serves as a lightweight local fallback aligner for benchmark sequences when remote servers time out.
- **Caching pays off fast:** SQLite caching prevents repeat BLAST and Entrez calls, making demoing and re-analysis instantaneous (< 10ms).
- **Zero recurring cost:** Entire stack remains 100% free and open-source.

---

## 7. Fallback & Resilience Layer

Every external dependency (BLAST, Entrez, UniProt) has a complete fallback chain: **Cache → Retry → Secondary source → Graceful degradation.**

### BLAST (identification)
1. **Cache first** — check SQLite cache for the sequence MD5 hash.
2. **Retry** — tenacity `@retry` decorator with exponential backoff for network timeouts.
3. **Secondary path** — `Bio.Align.PairwiseAligner` against curated benchmark library if remote NCBI is unreachable.
4. **Degrade gracefully** — if no matches exist (e.g. synthetic or non-existent sequence), display authentic zero hits and explanation.

### NCBI Entrez (annotation) and PubMed
1. **Cache first**, same as above.
2. **Retry** with tenacity, respecting NCBI rate limits.
3. **Secondary path** — resolve canonical gene and product from BLAST headers if GenBank locus tags are incomplete.
4. **Degrade gracefully** — mark failed fields clearly without breaking the dashboard.

### UniProt
1. **Cache first.**
2. **Retry** with tenacity.
3. **Secondary path** — direct Swiss-Prot accession lookup or biological product name query.
4. **Degrade gracefully** — display clear unavailable status if no reviewed entry exists.

---

## 8. Checklist Against Your Synopsis's Objectives

- [x] Objective 1 — HTML/CSS/JS submission interface
- [x] Objective 2 — FASTA/raw input validation + sequence-type detection
- [x] Objective 3 — NCBI BLAST similarity search
- [x] Objective 4 — Confidence-ranked candidate homologs (E-value/coverage/identity)
- [x] Objective 5 — NCBI sequence/gene/taxonomy retrieval
- [x] Objective 6 — UniProt functional annotation (domains, GO terms, cross-refs)
- [x] Objective 7 — PubMed literature retrieval
- [x] Objective 8 — Literature organized by topic (function/disease/structure/drug/interaction)
- [x] Objective 9 — Single dashboard presenting all results
- [x] Objective 10 — Structured HTML/PDF report
