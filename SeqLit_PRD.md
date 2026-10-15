# SeqLit — Product Requirements Document

| | |
|---|---|
| **Status** | Approved v1.0 |
| **Owner** | Rizwan Ahmed |
| **Type** | MSc Bioinformatics Minor Project |
| **Related docs** | `SeqLit_Project_Execution_Guide.md` (setup, week-by-week build plan), `TEST_LOG.md` |

---

## 1. Overview

SeqLit is a web application that takes a raw DNA, RNA, or protein sequence and returns, in one place, what it is, what it does, and what's been published about it. It replaces a manual workflow that currently spans NCBI BLAST, NCBI Entrez, UniProt, and PubMed as separate steps with a single submission → dashboard → report pipeline.

## 2. Problem Statement

A researcher or student with an unidentified sequence has to move manually between at least four separate tools — BLAST for identification, Entrez for sequence/gene/taxonomy data, UniProt for functional annotation, and PubMed for literature — reconciling results by hand each time. This is slow, repetitive, and error-prone, particularly for students doing preliminary characterization work who don't yet have a practiced routine across all four resources.

## 3. Goals & Success Metrics

| Goal | Success metric |
|---|---|
| Reduce time-to-first-insight for an unknown sequence | A full run (submission → populated dashboard) completes without manual intervention |
| Correct, confidence-ranked identification | Top hit matches ground truth on a benchmark set of known test sequences (Week 7 test log) |
| Usable, single-pane result view | All four result categories (Identification, Annotation, Literature, Report) visible from one dashboard, no tab-switching to external sites required to get the summary |
| Zero recurring cost | Entire stack — software and data sources — runs at $0, indefinitely |
| Resilient to third-party outages | A BLAST/Entrez/UniProt outage degrades the report (missing section, clearly marked) rather than crashing the app (see Section 10) |

## 4. Target Users

- **Primary:** the student/researcher submitting sequences for preliminary characterization (also the developer, in this case).
- **Secondary:** anyone evaluating the project — instructor, viva panel — who needs to see a working, coherent tool with a sensible architecture, not just a script.

## 5. Scope

### In scope
- Single-sequence submission (paste or FASTA), DNA/RNA/protein, auto-typed.
- BLAST-based similarity search and ranked hit list.
- NCBI-sourced sequence/gene/taxonomy annotation.
- UniProt-sourced functional annotation (function, GO terms, domains, cross-references).
- PubMed literature search, de-duplication, rule-based topic categorization.
- Single interactive dashboard; downloadable HTML/PDF report.
- Local caching of prior results.
- Fallback behavior when any external data source is unavailable (Section 10).

### Out of scope (this iteration)
- Batch/multi-sequence upload.
- User accounts, saved history across sessions, multi-user support.
- Structural prediction or 3D visualization.
- Integration with InterPro, PDB, STRING, or pathway databases (flagged as future extensions in the synopsis).
- NLP/ML-based literature summarization (categorization stays rule-based/keyword-driven).
- Deployment to a public server — this is a local, single-user application for the project's purposes.

## 6. Functional Requirements

| ID | Requirement | Priority | Acceptance criteria | Status |
|---|---|---|---|---|
| FR-1 | Accept raw or FASTA-formatted sequence input via a web form | P0 | Submitting valid FASTA or plain-text sequence returns no validation error | Completed |
| FR-2 | Validate input and auto-detect sequence type (DNA/RNA/protein) | P0 | Invalid characters are rejected with a specific error; valid sequences are correctly typed on a 10-sequence spot-check | Completed |
| FR-3 | Compute basic sequence statistics (length, GC%, molecular weight as applicable) | P0 | Stats display alongside the detected type | Completed |
| FR-4 | Run a BLAST similarity search and return ranked hits (accession, %identity, E-value, coverage) | P0 | A known test sequence returns its correct top hit within the confidence cutoff | Completed |
| FR-5 | Retrieve organism, gene, and taxonomy data for top hit(s) from NCBI | P0 | Annotation panel is populated for any sequence with a successful BLAST hit | Completed |
| FR-6 | Retrieve function, GO terms, domains, and cross-references from UniProt | P0 | Annotation panel shows UniProt fields for protein-coding hits | Completed |
| FR-7 | Search PubMed using identified gene/protein/organism terms, de-duplicate by PMID | P1 | No duplicate PMIDs in the literature list across repeated runs | Completed |
| FR-8 | Categorize literature by topic (function/disease/structure/drug/interaction) via keyword rules | P1 | Each listed publication carries at least one category label | Completed |
| FR-9 | Present all results (identification, annotation, literature, stats) on one dashboard | P0 | All four panels render from a single submission, no page navigation required | Completed |
| FR-10 | Generate a downloadable structured HTML/PDF report | P0 | Report contains sequence info, top hits, annotation summary, and literature list | Completed |
| FR-11 | Cache results per unique input sequence | P1 | Re-submitting an identical sequence returns instantly from cache, with no external calls made | Completed |

## 7. Non-Functional Requirements

- **Cost:** the entire stack — software, libraries, and data sources — must remain free at every layer; no dependency introduces a paid tier at this project's scale (see the Execution Guide's tech stack tables for the full free/open-source inventory).
- **Performance:** a cached lookup returns in under 1 second; an uncached full pipeline run completes within the time a remote BLAST search takes (typically well under a few minutes), with a visible loading state throughout.
- **Reliability:** no single external dependency failure should take down the whole app — see Section 10, Fallback & Resilience.
- **Portability:** runs locally on macOS with the setup in the Execution Guide; no cloud account or paid infrastructure required.
- **Usability:** results are legible without needing to open any external site — annotation and literature summaries are self-contained on the dashboard, with outbound links offered as a convenience, not a requirement.
- **Maintainability:** each pipeline stage (validation, BLAST, NCBI annotation, UniProt annotation, PubMed, report) is an independently testable module (see `modules/` structure in the Execution Guide).

## 8. System Architecture & Tech Stack

Full install instructions and the complete free/open-source tool list live in `SeqLit_Project_Execution_Guide.md`. Summary:

- **Backend:** Python, Flask, Biopython, `requests`, Pandas
- **Frontend:** HTML5/CSS3/vanilla JS, [Tabler](https://github.com/tabler/tabler) (UI kit), [Lucide](https://github.com/lucide-icons/lucide) (icons), Chart.js (charts)
- **Data sources:** NCBI BLAST, NCBI Entrez/E-utilities, UniProt REST API, PubMed
- **Resilience layer:** [tenacity](https://github.com/jd/tenacity) (retry/backoff), Biopython PairwiseAligner (curated local reference aligner fallback)
- **Reporting:** Jinja2, Matplotlib
- **Storage/testing:** SQLite, pytest

## 9. UI/UX Requirements

- Layout: sidebar navigation + four content panels (Identification, Annotation, Literature, Report), following Tabler's standard dashboard patterns rather than a single unstyled page — the goal is an interface that reads as a considered product, not a default framework scaffold.
- Use Tabler's card, table, and form components directly rather than rebuilding them from scratch in custom CSS — this also keeps styling consistent panel-to-panel with minimal extra work.
- Use Lucide icons (not emoji or mismatched icon sets) for panel headers and status indicators, for visual consistency.
- Every panel must have a defined loading state and a defined empty/failed state (see Section 10) — no panel should ever render blank with no explanation.
- Charts (Chart.js) supplement tables, never replace them — every chart's underlying data must also be visible in tabular form.
- Responsive down to tablet width at minimum; phone support is a nice-to-have, not a requirement, given this is a working tool rather than a public product.

## 10. Fallback & Resilience Requirements

Full technical detail (per-dependency fallback chains, retry parameters, logging) is specified in Section 7, "Fallback & Resilience Layer," of `SeqLit_Project_Execution_Guide.md`. Summary of the requirement:

- **FR-Resilience-1:** every external call (BLAST, Entrez, UniProt) goes through a cache check first, then a `tenacity`-managed retry with exponential backoff, before being treated as failed.
- **FR-Resilience-2:** BLAST has a secondary local path (curated reference aligner) for when the remote service is slow or unreachable.
- **FR-Resilience-3:** NCBI/UniProt annotation has biological product and direct accession fallbacks when gene tags are incomplete.
- **FR-Resilience-4:** any pipeline stage that ultimately fails degrades its own panel only — with a clear, specific message — and never blocks the other panels from rendering.
- **FR-Resilience-5:** every fallback trigger is logged (stage, fallback used, reason) for later review.

## 11. Milestones

Maps to the 8-week plan in the Execution Guide:

| Week | Milestone | Status |
|---|---|---|
| 1 | Environment set up; validation module working end to end | Completed |
| 2 | BLAST integration returns correct ranked hits | Completed |
| 3 | NCBI + UniProt annotation populated | Completed |
| 4 | PubMed retrieval, de-duplication, categorization working | Completed |
| 5 | Full dashboard wired to live data | Completed |
| 6 | HTML/PDF report generation working | Completed |
| 7 | Test log complete; known/edge-case sequences verified | Completed |
| 8 | Caching, fallback layer, docs, and presentation complete | Completed |

## 12. Risks & Assumptions

| Risk | Mitigation |
|---|---|
| NCBI/UniProt/PubMed downtime or rate-limiting during development or demo | Cache + retry + fallback layer (Section 10) |
| Remote BLAST latency makes live demos slow | Cache warm results ahead of a demo; local benchmark aligner fallback |
| Scrapling/API schema shifts | Adaptive entity extraction & multi-tier query fallback |
| Scope creep toward future extensions (InterPro, PDB, STRING, pathway DBs) | Explicitly out of scope for this iteration — see Section 5 |

## 13. Appendix — Reference Repositories

| Repository | Purpose |
|---|---|
| [biopython/biopython](https://github.com/biopython/biopython) | Sequence parsing, BLAST calls, Entrez access, PairwiseAligner |
| [pallets/flask](https://github.com/pallets/flask) | Backend web framework |
| [tabler/tabler](https://github.com/tabler/tabler) | Dashboard UI kit |
| [lucide-icons/lucide](https://github.com/lucide-icons/lucide) | Icon set |
| [chartjs/Chart.js](https://github.com/chartjs/Chart.js) | Dashboard charts |
| [jd/tenacity](https://github.com/jd/tenacity) | Retry/backoff for external calls |
| [pandas-dev/pandas](https://github.com/pandas-dev/pandas) | Tabular data handling |
