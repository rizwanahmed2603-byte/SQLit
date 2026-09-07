// Client-side interactions, AJAX execution, Chart.js updates and filter logic for SeqLit

let compositionChartInstance = null;
let literatureChartInstance = null;
let lastAnalysisResult = null;

document.addEventListener("DOMContentLoaded", function() {
  const form = document.getElementById("sequenceForm");
  const clearBtn = document.getElementById("clearBtn");
  const loadSampleBtn = document.getElementById("loadSampleBtn");
  const exportPdfBtn = document.getElementById("exportPdfBtn");
  const categoryFilter = document.getElementById("categoryFilter");

  if (form) {
    form.addEventListener("submit", handleSequenceSubmit);
  }

  if (clearBtn) {
    clearBtn.addEventListener("click", () => {
      document.getElementById("sequenceInput").value = "";
      hideError();
      document.getElementById("resultsDashboard").classList.add("d-none");
      exportPdfBtn.classList.add("d-none");
    });
  }

  if (loadSampleBtn) {
    loadSampleBtn.addEventListener("click", loadSampleSequence);
  }

  if (exportPdfBtn) {
    exportPdfBtn.addEventListener("click", downloadReport);
  }

  if (categoryFilter) {
    categoryFilter.addEventListener("change", (e) => {
      filterLiteratureTable(e.target.value);
    });
  }
});

async function loadSampleSequence() {
  try {
    const res = await fetch("/api/sample");
    const data = await res.json();
    if (data && data.sequence) {
      document.getElementById("sequenceInput").value = data.fasta || data.sequence;
    }
  } catch (err) {
    // Default fallback sample: Human HBB fragment
    document.getElementById("sequenceInput").value =
`>Human_HBB_hemoglobin_beta_segment
ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAA`;
  }
}

async function handleSequenceSubmit(e) {
  e.preventDefault();
  hideError();

  const seqText = document.getElementById("sequenceInput").value.trim();
  if (!seqText) {
    showError("Please enter or paste a valid DNA, RNA, or Protein sequence.");
    return;
  }

  showProgress("Validating sequence and inspecting type...", 15);

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
    renderDashboard(result.data);
    hideProgress();

    // Show report export button
    document.getElementById("exportPdfBtn").classList.remove("d-none");

  } catch (err) {
    hideProgress();
    showError("Network or server connection failed: " + err.message);
  }
}

function showProgress(text, pct) {
  const pSection = document.getElementById("progressSection");
  const pText = document.getElementById("pipelineStatusText");
  const pPct = document.getElementById("progressPercent");
  const submitBtn = document.getElementById("submitBtn");

  pSection.classList.remove("d-none");
  pText.textContent = text;
  pPct.textContent = pct + "%";
  submitBtn.disabled = true;
}

function hideProgress() {
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

function renderDashboard(data) {
  const dashboard = document.getElementById("resultsDashboard");
  dashboard.classList.remove("d-none");

  // Panel 1: Validation & Stats
  const val = data.validation || {};
  const stats = val.stats || {};
  document.getElementById("statSeqType").textContent = val.type || "Unknown";
  document.getElementById("statSeqLength").innerHTML = `${stats.length || 0} <span class="fs-5 fw-normal text-muted">${val.type === 'Protein' ? 'aa' : 'bp'}</span>`;
  document.getElementById("statGcPercent").textContent = stats.gc_percent !== null && stats.gc_percent !== undefined ? `${stats.gc_percent}%` : "N/A";
  document.getElementById("statMolWeight").innerHTML = `${(stats.molecular_weight || 0).toLocaleString()} <span class="fs-5 fw-normal text-muted">Da</span>`;

  // Render Composition Chart
  renderCompositionChart(stats.frequencies || {});

  // Panel 2: BLAST Identification
  const blast = data.blast || {};
  document.getElementById("blastProgramBadge").textContent = blast.program || "blastn";
  renderBlastTable(blast.hits || []);

  // Panel 3: Annotations
  const ncbi = data.ncbi || {};
  document.getElementById("ncbiAccession").textContent = ncbi.accession || "N/A";
  document.getElementById("ncbiOrganism").textContent = ncbi.organism || "N/A";
  document.getElementById("ncbiGene").textContent = ncbi.gene || "N/A";
  document.getElementById("ncbiTaxonomy").textContent = ncbi.taxonomy || "N/A";
  document.getElementById("ncbiTitle").textContent = ncbi.title || "N/A";

  const uniprot = data.uniprot || {};
  document.getElementById("uniprotProteinName").textContent = uniprot.protein_name || "Uncharacterized";
  document.getElementById("uniprotFunction").textContent = uniprot.function || "No functional annotation.";
  if (uniprot.entry_url) {
    document.getElementById("uniprotLink").href = uniprot.entry_url;
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
    allGo.slice(0, 8).forEach(go => {
      const span = document.createElement("span");
      span.className = "badge-go";
      span.title = go.id;
      span.textContent = go.name;
      goContainer.appendChild(span);
    });
  } else {
    goContainer.innerHTML = '<span class="text-muted small">No direct GO terms mapped.</span>';
  }

  // Panel 4: Literature
  const literature = data.literature || [];
  renderLiteratureTable(literature);
  renderLiteratureChart(literature);

  // Refresh icons for dynamically added contents
  if (window.lucide) {
    lucide.createIcons();
  }

  // Smooth scroll down to dashboard
  dashboard.scrollIntoView({ behavior: 'smooth' });
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

function renderBlastTable(hits) {
  const tbody = document.getElementById("blastHitsTableBody");
  tbody.innerHTML = "";

  if (!hits || hits.length === 0) {
    tbody.innerHTML = `<tr><td colspan="6" class="text-center text-muted py-3">No significant hits found.</td></tr>`;
    return;
  }

  hits.forEach(hit => {
    const tr = document.createElement("tr");
    tr.innerHTML = `
      <td><span class="badge bg-light text-dark">${hit.rank}</span></td>
      <td><span class="badge bg-blue-lt fw-bold">${hit.accession}</span></td>
      <td class="text-truncate" style="max-width: 320px;" title="${hit.title}">${hit.title}</td>
      <td><span class="fw-bold text-success">${hit.identity_percent}%</span></td>
      <td class="font-monospace small">${hit.e_value}</td>
      <td>
        <div class="d-flex align-items-center gap-2">
          <span>${hit.query_coverage}%</span>
          <div class="progress progress-xs w-100">
            <div class="progress-bar bg-success" style="width: ${hit.query_coverage}%"></div>
          </div>
        </div>
      </td>
    `;
    tbody.appendChild(tr);
  });
}

function renderLiteratureTable(articles) {
  const tbody = document.getElementById("literatureTableBody");
  tbody.innerHTML = "";

  if (!articles || articles.length === 0) {
    tbody.innerHTML = `<tr><td colspan="4" class="text-center text-muted py-3">No PubMed articles retrieved.</td></tr>`;
    return;
  }

  articles.forEach(art => {
    const tr = document.createElement("tr");
    tr.dataset.category = art.category;
    tr.innerHTML = `
      <td><a href="${art.pubmed_url}" target="_blank" class="fw-bold">${art.pmid}</a></td>
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
    const c = a.category || "General";
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
