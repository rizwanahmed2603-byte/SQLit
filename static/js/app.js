// Client-side interactions, AJAX execution, Chart.js updates, candidate hit inspection, and report export for SeqLit

let compositionChartInstance = null;
let literatureChartInstance = null;
let lastAnalysisResult = null;
let currentActiveHitRank = 1;

document.addEventListener("DOMContentLoaded", function() {
  const form = document.getElementById("sequenceForm");
  const clearBtn = document.getElementById("clearBtn");
  const exportPdfBtn = document.getElementById("exportPdfBtn");
  const categoryFilter = document.getElementById("categoryFilter");

  if (form) {
    form.addEventListener("submit", handleSequenceSubmit);
  }

  if (clearBtn) {
    clearBtn.addEventListener("click", () => {
      document.getElementById("sequenceInput").value = "";
      hideError();
      hideNoHitsAlert();
      document.getElementById("resultsDashboard").classList.add("d-none");
      exportPdfBtn.classList.add("d-none");
      document.getElementById("activeSampleBadge").textContent = "Ready";
    });
  }

  // Setup benchmark sample options
  document.querySelectorAll(".sample-option").forEach(item => {
    item.addEventListener("click", function(e) {
      e.preventDefault();
      const sampleId = this.dataset.id || "hbb";
      loadSampleSequence(sampleId);
    });
  });

  if (exportPdfBtn) {
    exportPdfBtn.addEventListener("click", downloadReport);
  }

  if (categoryFilter) {
    categoryFilter.addEventListener("change", (e) => {
      filterLiteratureTable(e.target.value);
    });
  }
});

async function loadSampleSequence(sampleId = "hbb") {
  try {
    const res = await fetch(`/api/sample?id=${sampleId}`);
    const data = await res.json();
    if (data && (data.fasta || data.sequence)) {
      document.getElementById("sequenceInput").value = data.fasta || data.sequence;
      const badge = document.getElementById("activeSampleBadge");
      if (badge) {
        badge.textContent = data.name || sampleId.toUpperCase();
      }
    }
  } catch (err) {
    console.error("Could not load sample:", err);
  }
}

async function handleSequenceSubmit(e) {
  e.preventDefault();
  hideError();
  hideNoHitsAlert();

  const seqText = document.getElementById("sequenceInput").value.trim();
  if (!seqText) {
    showError("Please enter or paste a valid DNA, RNA, or Protein sequence.");
    return;
  }

  showProgress("Querying NCBI BLAST similarity database (or local aligner fallback)...", 30);

  try {
    const response = await fetch("/api/analyze", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ sequence: seqText })
    });

    const result = await response.json();

    if (!response.ok || !result.success) {
      hideProgress();
      showError(result.error || "Sequence analysis could not be completed.");
      return;
    }

    lastAnalysisResult = result.data;
    currentActiveHitRank = 1;
    renderDashboard(result.data);
    hideProgress();

    // Show report export button
    document.getElementById("exportPdfBtn").classList.remove("d-none");

  } catch (err) {
    hideProgress();
    showError("Network or server connection failed: " + err.message);
  }
}

let progressTimer = null;
let progressStartTime = null;

function showProgress(text, pct) {
  const pSection = document.getElementById("progressSection");
  const pText = document.getElementById("pipelineStatusText");
  const pPct = document.getElementById("progressPercent");
  const submitBtn = document.getElementById("submitBtn");

  pSection.classList.remove("d-none");
  pPct.textContent = pct + "%";
  submitBtn.disabled = true;

  if (progressTimer) clearInterval(progressTimer);
  progressStartTime = Date.now();

  const updateStatus = () => {
    const elapsed = Math.floor((Date.now() - progressStartTime) / 1000);
    pText.textContent = `${text} (elapsed: ${elapsed}s)`;
  };
  updateStatus();
  progressTimer = setInterval(updateStatus, 1000);
}

function hideProgress() {
  if (progressTimer) {
    clearInterval(progressTimer);
    progressTimer = null;
  }
  document.getElementById("progressSection").classList.add("d-none");
  document.getElementById("submitBtn").disabled = false;
}

function showError(msg) {
  const alertEl = document.getElementById("errorAlert");
  document.getElementById("errorMessage").textContent = msg;
  alertEl.classList.remove("d-none");
}

function hideError() {
  document.getElementById("errorAlert").classList.add("d-none");
}

function showNoHitsAlert(msg) {
  const alertEl = document.getElementById("noHitsAlert");
  if (msg) {
    document.getElementById("noHitsMessage").textContent = msg;
  }
  alertEl.classList.remove("d-none");
}

function hideNoHitsAlert() {
  document.getElementById("noHitsAlert").classList.add("d-none");
}

function renderDashboard(data) {
  const dashboard = document.getElementById("resultsDashboard");
  dashboard.classList.remove("d-none");

  // Panel 1: Validation & Stats
  const val = data.validation || {};
  const stats = val.stats || {};
  document.getElementById("statSeqType").textContent = val.type || "Unknown";
  document.getElementById("statSeqLength").innerHTML = `${stats.length || 0} <span class="fs-5 fw-normal text-muted">${val.type === 'Protein' ? 'aa' : 'bp'}</span>`;
  document.getElementById("statGcPercent").textContent = (stats.gc_percent !== null && stats.gc_percent !== undefined) ? `${stats.gc_percent}%` : "N/A";
  document.getElementById("statMolWeight").innerHTML = `${(stats.molecular_weight || 0).toLocaleString()} <span class="fs-5 fw-normal text-muted">Da</span>`;

  // Render Composition Chart
  renderCompositionChart(stats.frequencies || {});

  // Panel 2: BLAST Identification
  const blast = data.blast || {};
  const hits = blast.hits || [];
  document.getElementById("blastProgramBadge").textContent = blast.program || "blastn";
  document.getElementById("blastHitCountBadge").textContent = `${hits.length} hit${hits.length === 1 ? '' : 's'}`;
  renderBlastTable(hits, currentActiveHitRank);

  if (hits.length === 0) {
    showNoHitsAlert("No significant biological homologs were found in NCBI nr/nt databases for this sequence. Real sequence stats have been computed, but no annotations or literature exist for unrecorded sequences.");
  } else {
    hideNoHitsAlert();
  }

  // Panel 3 & 4: Annotations & Literature for current candidate
  renderAnnotationsAndLiterature(data.ncbi || {}, data.uniprot || {}, data.literature || [], currentActiveHitRank);

  if (window.lucide) {
    lucide.createIcons();
  }

  dashboard.scrollIntoView({ behavior: 'smooth' });
}

function renderAnnotationsAndLiterature(ncbi, uniprot, literature, rank) {
  // NCBI Card
  document.getElementById("ncbiCandidateBadge").textContent = `Candidate #${rank}`;
  document.getElementById("ncbiAccession").textContent = ncbi.accession || "N/A";
  document.getElementById("ncbiOrganism").textContent = ncbi.organism || "N/A";
  document.getElementById("ncbiGene").textContent = ncbi.gene || "N/A";
  document.getElementById("ncbiProduct").textContent = ncbi.product || "N/A";
  document.getElementById("ncbiTaxonomy").textContent = ncbi.taxonomy || "N/A";
  document.getElementById("ncbiTitle").textContent = ncbi.title || "N/A";

  // UniProt Card
  document.getElementById("uniprotProteinName").textContent = uniprot.protein_name || "Uncharacterized";
  document.getElementById("uniprotFunction").textContent = uniprot.function || "No functional annotation recorded.";
  const uniprotLink = document.getElementById("uniprotLink");
  if (uniprot.entry_url) {
    uniprotLink.href = uniprot.entry_url;
    uniprotLink.classList.remove("d-none");
  } else {
    uniprotLink.classList.add("d-none");
  }

  // GO Terms
  const goContainer = document.getElementById("uniprotGoTerms");
  goContainer.innerHTML = "";
  const goGroups = uniprot.go_terms || {};
  const allGo = [
    ...(goGroups.biological_process || []),
    ...(goGroups.molecular_function || []),
    ...(goGroups.cellular_component || [])
  ];

  if (allGo.length > 0) {
    allGo.slice(0, 10).forEach(go => {
      const span = document.createElement("span");
      span.className = "badge-go";
      span.title = go.id;
      span.textContent = `${go.name} (${go.id})`;
      goContainer.appendChild(span);
    });
  } else {
    goContainer.innerHTML = '<span class="text-muted small">No direct GO terms mapped.</span>';
  }

  // PubMed Literature
  document.getElementById("literatureCountBadge").textContent = `${literature.length} paper${literature.length === 1 ? '' : 's'}`;
  renderLiteratureTable(literature);
  renderLiteratureChart(literature);
}

function renderCompositionChart(frequencies) {
  const ctx = document.getElementById("compositionChart");
  if (!ctx) return;

  if (compositionChartInstance) {
    compositionChartInstance.destroy();
  }

  const labels = Object.keys(frequencies);
  const data = Object.values(frequencies);

  compositionChartInstance = new Chart(ctx, {
    type: "bar",
    data: {
      labels: labels,
      datasets: [{
        label: "Composition (%)",
        data: data,
        backgroundColor: "#206bc4",
        borderRadius: 4
      }]
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      plugins: {
        legend: { display: false }
      },
      scales: {
        y: {
          beginAtZero: true,
          max: 100,
          ticks: { callback: v => v + "%" }
        }
      }
    }
  });
}

function renderBlastTable(hits, activeRank = 1) {
  const tbody = document.getElementById("blastHitsTableBody");
  tbody.innerHTML = "";

  if (!hits || hits.length === 0) {
    tbody.innerHTML = `<tr><td colspan="7" class="text-center text-muted py-4">No significant biological hits found in NCBI database.</td></tr>`;
    return;
  }

  hits.forEach(hit => {
    const tr = document.createElement("tr");
    const isActive = (hit.rank === activeRank);
    if (isActive) {
      tr.className = "table-primary fw-medium";
    }

    tr.innerHTML = `
      <td><span class="badge ${isActive ? 'bg-primary text-white' : 'bg-light text-dark'}">${hit.rank}</span></td>
      <td>
        <span class="badge bg-blue-lt fw-bold font-monospace">${hit.accession}</span>
        ${hit.is_primary ? '<span class="badge bg-green text-white ms-1" style="font-size: 10px;">Primary</span>' : ''}
      </td>
      <td class="text-truncate" style="max-width: 300px;" title="${hit.title}">${hit.title}</td>
      <td><span class="fw-bold ${hit.identity_percent >= 98 ? 'text-success' : 'text-primary'}">${hit.identity_percent}%</span></td>
      <td class="font-monospace small">${hit.e_value}</td>
      <td>
        <div class="d-flex align-items-center gap-2">
          <span>${hit.query_coverage}%</span>
          <div class="progress progress-xs" style="width: 50px;">
            <div class="progress-bar bg-success" style="width: ${hit.query_coverage}%"></div>
          </div>
        </div>
      </td>
      <td>
        <button class="btn btn-sm ${isActive ? 'btn-primary' : 'btn-outline-primary'} py-0 px-2 inspect-btn" data-rank="${hit.rank}">
          ${isActive ? 'Active' : 'Inspect'}
        </button>
      </td>
    `;

    const inspectBtn = tr.querySelector(".inspect-btn");
    inspectBtn.addEventListener("click", () => inspectHit(hit));

    tbody.appendChild(tr);
  });
}

async function inspectHit(hit) {
  if (!hit || !lastAnalysisResult) return;

  currentActiveHitRank = hit.rank;
  lastAnalysisResult.active_hit_rank = hit.rank;

  // Re-render table to update active row styling
  renderBlastTable(lastAnalysisResult.blast.hits, currentActiveHitRank);

  // If this hit's annotation is already in lastAnalysisResult and rank == 1
  if (hit.rank === 1 && lastAnalysisResult.ncbi && lastAnalysisResult.ncbi.accession === hit.accession) {
    renderAnnotationsAndLiterature(
      lastAnalysisResult.ncbi,
      lastAnalysisResult.uniprot,
      lastAnalysisResult.literature,
      hit.rank
    );
    return;
  }

  // Fetch annotation for selected hit via API
  try {
    const res = await fetch("/api/annotate_hit", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        accession: hit.accession,
        seq_type: lastAnalysisResult.validation.type || "DNA",
        title: hit.title,
        gene_symbol: hit.gene_symbol || "",
        organism: hit.organism || "",
        uniprot_id: hit.uniprot_id || ""
      })
    });

    const resData = await res.json();
    if (resData.success && resData.data) {
      const d = resData.data;
      lastAnalysisResult.ncbi = d.ncbi;
      lastAnalysisResult.uniprot = d.uniprot;
      lastAnalysisResult.literature = d.literature;

      renderAnnotationsAndLiterature(d.ncbi, d.uniprot, d.literature, hit.rank);
      if (window.lucide) lucide.createIcons();
    }
  } catch (err) {
    console.error("Failed to inspect hit:", err);
  }
}

function renderLiteratureTable(articles) {
  const tbody = document.getElementById("literatureTableBody");
  tbody.innerHTML = "";

  if (!articles || articles.length === 0) {
    tbody.innerHTML = `<tr><td colspan="4" class="text-center text-muted py-4">No PubMed articles retrieved for this sequence.</td></tr>`;
    return;
  }

  articles.forEach(art => {
    const tr = document.createElement("tr");
    tr.dataset.category = art.category;
    tr.innerHTML = `
      <td><a href="${art.pubmed_url}" target="_blank" class="fw-bold font-monospace text-primary">${art.pmid}</a></td>
      <td>
        <div class="fw-medium">${art.title}</div>
        <div class="text-muted small">${art.authors}</div>
      </td>
      <td class="text-muted small">${art.journal} (${art.year})</td>
      <td><span class="badge bg-azure-lt">${art.category}</span></td>
    `;
    tbody.appendChild(tr);
  });
}

function renderLiteratureChart(articles) {
  const ctx = document.getElementById("literatureChart");
  if (!ctx) return;

  if (literatureChartInstance) {
    literatureChartInstance.destroy();
  }

  const catCounts = {};
  articles.forEach(a => {
    const c = a.category || "General / Biology";
    catCounts[c] = (catCounts[c] || 0) + 1;
  });

  const labels = Object.keys(catCounts);
  const data = Object.values(catCounts);

  const colors = ["#206bc4", "#4299e1", "#2fb344", "#f76707", "#d63939", "#ae3ec9"];

  literatureChartInstance = new Chart(ctx, {
    type: "doughnut",
    data: {
      labels: labels,
      datasets: [{
        data: data,
        backgroundColor: colors.slice(0, labels.length)
      }]
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      plugins: {
        legend: { position: "bottom", labels: { boxWidth: 12, font: { size: 11 } } }
      }
    }
  });
}

function filterLiteratureTable(category) {
  const rows = document.querySelectorAll("#literatureTableBody tr");
  rows.forEach(row => {
    if (category === "ALL" || row.dataset.category === category) {
      row.style.display = "";
    } else {
      row.style.display = "none";
    }
  });
}

async function downloadReport() {
  if (!lastAnalysisResult) return;

  try {
    const res = await fetch("/api/export/report", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ data: lastAnalysisResult })
    });

    if (res.ok) {
      const blob = await res.blob();
      const url = window.URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = `seqlit_report_${Date.now()}.html`;
      document.body.appendChild(a);
      a.click();
      a.remove();
    }
  } catch (err) {
    alert("Report download failed: " + err.message);
  }
}
