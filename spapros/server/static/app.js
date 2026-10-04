// Browser GUI for the spapros server: start runs, follow them, and read the feasibility report.
// Plain JavaScript with no build step; all requests use relative URLs so the app also works behind a path prefix.
"use strict";

const view = document.getElementById("view");
const STAGES = ["loading", "selecting", "selecting_reference_sets", "evaluating", "panel_size_curve", "reporting"];
const STAGE_LABEL = {
  loading: "Loading data", selecting: "Selecting genes", selecting_reference_sets: "Baseline gene sets",
  evaluating: "Evaluating", panel_size_curve: "Panel size curve", reporting: "Writing report",
};
const VERDICT_STYLE = {
  "FEASIBLE": "go", "FEASIBLE WITH CAVEATS": "caveat", "NOT FEASIBLE": "fail", "INCONCLUSIVE": "unknown",
};
const STATUS_STYLE = { queued: "unknown", running: "caveat", succeeded: "go", failed: "fail", cancelled: "unknown" };
const N_CELLS_MIN = 40;
// Feasibility thresholds shown in the advanced settings, with the server defaults.
const THRESHOLDS = [
  ["acc_resolved", "Recall for a resolved cell type", 0.8],
  ["acc_marginal", "Recall for a marginal cell type", 0.6],
  ["frac_resolved_go", "Share of types resolved for a go", 0.9],
  ["frac_resolved_caveat", "Share of types resolved below which it fails", 0.7],
  ["mean_acc_go", "Mean recall for a go", 0.85],
  ["mean_acc_caveat", "Mean recall below which it fails", 0.75],
  ["beats_random_margin", "Margin over random gene sets", 0.05],
  ["baseline_tolerance", "Allowed gap to the best baseline", 0.02],
  ["headroom", "Share of the budget to keep free", 0.2],
  ["secondary_gap_vs_pca", "Allowed gap to PCA in clustering / kNN", 0.1],
  ["gene_corr_min", "Minimum share of non-redundant genes", 0.9],
  ["marker_corr_min", "Minimum marker correlation", 0.5],
  ["max_excluded_types", "Max share of cell types too small to assess", 0.2],
  ["max_excluded_cells", "Max share of cells in such types", 0.05],
];

let timer = null;

function esc(s) {
  return String(s ?? "").replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[c]);
}
function badge(text, style) { return `<span class="badge ${style}">${esc(text)}</span>`; }
function when(iso) { return iso ? new Date(iso).toLocaleString() : ""; }
function duration(a, b) {
  if (!a) return "";
  const s = Math.max(0, ((b ? new Date(b) : new Date()) - new Date(a)) / 1000);
  return s < 90 ? `${Math.round(s)} s` : s < 5400 ? `${Math.round(s / 60)} min` : `${(s / 3600).toFixed(1)} h`;
}

async function api(path, opts = {}) {
  const r = await fetch(path, opts);
  if (!r.ok) {
    let msg = `${r.status} ${r.statusText}`;
    try {
      const body = await r.json();
      msg = typeof body.detail === "string" ? body.detail
        : Array.isArray(body.detail) ? body.detail.map((d) => `${(d.loc || []).join(".")}: ${d.msg}`).join("; ")
        : JSON.stringify(body.detail);
    } catch (e) { /* not JSON */ }
    throw new Error(msg);
  }
  return r.status === 204 ? null : r.json();
}

// Runs list --------------------------------------------------------------------------------------------------------

async function showRuns() {
  const jobs = (await api("jobs")).reverse();
  if (!jobs.length) {
    view.innerHTML = `<h1>Runs</h1><div class="card"><p>No runs yet.</p>
      <p class="muted">A run selects a gene panel with spapros, evaluates it against simple baselines and ends with a
      feasibility report.</p><a class="button primary" href="#/new">Start a run</a></div>`;
    return;
  }
  const rows = jobs.map((j) => {
    const v = j.verdict ? badge(j.verdict.verdict, VERDICT_STYLE[j.verdict.verdict] || "unknown") : "";
    const stage = j.status === "running" && j.stage ? ` <span class="muted">${esc(STAGE_LABEL[j.stage] || j.stage)}</span>` : "";
    return `<tr><td><a href="#/jobs/${esc(j.id)}">${esc(j.input_filename)}</a></td>
      <td>${esc(j.options.celltype_key)}</td><td class="num">${esc(j.options.n)}</td>
      <td>${badge(j.status, STATUS_STYLE[j.status])}${stage}</td><td>${v}</td>
      <td class="muted">${esc(when(j.created_at))}</td></tr>`;
  }).join("");
  view.innerHTML = `<h1>Runs</h1><div class="scroll"><table><tr><th>Dataset</th><th>Cell type column</th>
    <th class="num">Genes</th><th>Status</th><th>Verdict</th><th>Started</th></tr>${rows}</table></div>`;
  if (jobs.some((j) => j.status === "queued" || j.status === "running")) timer = setTimeout(route, 5000);
}

// New run ----------------------------------------------------------------------------------------------------------

async function showNew() {
  const uploads = await api("uploads");
  const reuse = uploads.length ? `<label for="reuse">Or reuse an earlier upload</label>
    <select id="reuse"><option value="">–</option>${uploads.slice().reverse().map((u) =>
      `<option value="${esc(u.id)}">${esc(u.filename)} (${u.n_cells} cells, ${esc(when(u.created_at))})</option>`).join("")}
    </select>` : "";
  view.innerHTML = `<h1>New run</h1>
    <div class="card">
      <label for="file">Single-cell reference (.h5ad)
        <span class="hint">Raw counts in X and a cell type column in obs.</span></label>
      <input type="file" id="file" accept=".h5ad">
      <div class="progress" hidden><div></div></div><div id="upload-msg" class="muted"></div>
      ${reuse}
    </div>
    <form id="form" hidden></form>`;
  document.getElementById("file").addEventListener("change", (e) => e.target.files[0] && upload(e.target.files[0]));
  const sel = document.getElementById("reuse");
  if (sel) sel.addEventListener("change", async () => sel.value && renderForm(await api(`uploads/${sel.value}`)));
}

function upload(file) {
  const bar = document.querySelector(".progress");
  const msg = document.getElementById("upload-msg");
  bar.hidden = false;
  msg.textContent = "Uploading…";
  msg.className = "muted";
  const xhr = new XMLHttpRequest();
  xhr.open("POST", "uploads");
  xhr.upload.onprogress = (e) => {
    if (e.lengthComputable) {
      bar.firstElementChild.style.width = `${(100 * e.loaded) / e.total}%`;
      if (e.loaded === e.total) msg.textContent = "Reading the dataset…";
    }
  };
  xhr.onload = () => {
    let body = {};
    try { body = JSON.parse(xhr.responseText); } catch (e) { /* ignore */ }
    if (xhr.status === 201) { msg.textContent = ""; bar.hidden = true; renderForm(body); }
    else { msg.textContent = `Upload failed: ${body.detail || xhr.statusText}`; msg.className = "error"; }
  };
  xhr.onerror = () => { msg.textContent = "Upload failed: network error"; msg.className = "error"; };
  const fd = new FormData();
  fd.append("file", file);
  xhr.send(fd);
}

function guessCelltypeKey(cols) {
  const names = cols.map((c) => c.name);
  for (const pat of [/^cell_?type$/i, /cell_?type/i, /annot/i, /cluster|leiden|louvain/i]) {
    const hit = names.find((n) => pat.test(n));
    if (hit) return hit;
  }
  return names[0];
}

function renderForm(ds) {
  const form = document.getElementById("form");
  form.hidden = false;
  const cols = ds.obs_columns.filter((c) => c.values);
  if (!cols.length) {
    form.innerHTML = `<p class="error">No column in obs looks like a cell type annotation (2 to 500 distinct values).</p>`;
    return;
  }
  const boolCols = ds.var_bool_columns;
  const defaultGenesKey = boolCols.find((c) => c.name === "highly_variable") ? "highly_variable" : "";
  form.innerHTML = `
    <div class="card"><strong>${esc(ds.filename)}</strong>
      <span class="muted">${ds.n_cells} cells × ${ds.n_genes} genes</span></div>
    <div class="card">
      <label for="celltype_key">Cell type column</label>
      <select id="celltype_key">${cols.map((c) => `<option>${esc(c.name)}</option>`).join("")}</select>
      <label>Critical cell types <span class="hint">Types the panel must tell apart. Any unresolved one makes the
        run not feasible. Leave all unticked if none is critical.</span></label>
      <div class="checks" id="critical"></div>
    </div>
    <div class="card">
      <div class="row">
        <div><label for="panel_capacity">Panel capacity <span class="hint">genes the platform holds</span></label>
          <input type="number" id="panel_capacity" min="1" value="300"></div>
        <div><label for="reserved_slots">Reserved slots <span class="hint">controls, genes added by hand</span></label>
          <input type="number" id="reserved_slots" min="0" value="0"></div>
      </div>
      <div class="budget" id="budget"></div>
      <label for="genes_key">Select genes from</label>
      <select id="genes_key"><option value="">all genes (${ds.n_genes})</option>${boolCols.map((c) =>
        `<option value="${esc(c.name)}">${esc(c.name)} (${c.n_true} genes)</option>`).join("")}</select>
      <label><input type="checkbox" id="normalize" ${ds.looks_like_counts === false ? "" : "checked"}>
        Normalise and log-transform first <span class="hint">${ds.looks_like_counts === false
          ? "X does not look like raw counts, so this is off" : "X looks like raw counts"}</span></label>
    </div>
    <div class="card">
      <label for="marker_list">Marker list CSV <span class="hint">optional; one column per cell type, genes as rows
        </span></label><input type="file" id="marker_list" accept=".csv">
      <label for="preselected">Genes that must be on the panel <span class="hint">optional; comma or newline
        separated</span></label><textarea id="preselected"></textarea>
      <details><summary>Advanced settings</summary>
        <div class="row">
          <div><label for="evaluation_scheme">Evaluation</label><select id="evaluation_scheme">
            <option value="quick">quick</option><option value="full">full (adds clustering similarity)</option></select></div>
          <div><label for="seed">Seed</label><input type="number" id="seed" value="0"></div>
          <div><label for="n_jobs">CPUs <span class="hint">-1 uses all</span></label>
            <input type="number" id="n_jobs" value="-1"></div>
        </div>
        <label><input type="checkbox" id="evaluate_reference_sets" checked>Compare with PCA, DE, HVG and random gene
          sets</label>
        <label><input type="checkbox" id="panel_size_curve" checked>Estimate how many genes are needed (panel size
          curve)</label>
        <h2>Feasibility thresholds</h2>
        <div class="row">${THRESHOLDS.map(([k, label, d]) => `<div><label for="th_${k}">${esc(label)}</label>
          <input type="number" id="th_${k}" step="0.01" min="0" max="1" value="${d}"></div>`).join("")}</div>
      </details>
    </div>
    <div class="actions"><button class="button primary" type="submit" id="submit">Start run</button></div>
    <div id="form-msg"></div>`;

  const ctSel = document.getElementById("celltype_key");
  ctSel.value = guessCelltypeKey(cols);
  const renderCritical = () => {
    const col = cols.find((c) => c.name === ctSel.value);
    document.getElementById("critical").innerHTML = col.values.map((v) => {
      const small = v.count < N_CELLS_MIN;
      return `<label class="${small ? "disabled" : ""}" title="${small ? "Too few cells to be assessed" : ""}">
        <input type="checkbox" value="${esc(v.name)}" ${small ? "disabled" : ""}>${esc(v.name)}
        <span class="hint">&nbsp;${v.count} cells${small ? ", too few to assess" : ""}</span></label>`;
    }).join("");
  };
  ctSel.addEventListener("change", renderCritical);
  renderCritical();

  const gk = document.getElementById("genes_key");
  gk.value = defaultGenesKey;
  const renderBudget = () => {
    const c = +document.getElementById("panel_capacity").value;
    const r = +document.getElementById("reserved_slots").value;
    const pool = gk.value ? boolCols.find((x) => x.name === gk.value).n_true : ds.n_genes;
    const b = c - r;
    const el = document.getElementById("budget");
    el.innerHTML = b < 1 ? `<span class="error">Reserved slots leave no genes for spapros.</span>`
      : b > pool ? `<span class="error">spapros may use ${b} genes but only ${pool} candidate genes are available.</span>`
      : `spapros will select <strong>${b}</strong> genes from ${pool} candidates.`;
  };
  ["panel_capacity", "reserved_slots", "genes_key"].forEach((id) =>
    document.getElementById(id).addEventListener("input", renderBudget));
  renderBudget();

  form.onsubmit = async (e) => {
    e.preventDefault();
    const val = (id) => document.getElementById(id).value;
    const feasibility = {};
    THRESHOLDS.forEach(([k]) => { feasibility[k] = parseFloat(val(`th_${k}`)); });
    const options = {
      celltype_key: ctSel.value,
      panel_capacity: parseInt(val("panel_capacity"), 10),
      reserved_slots: parseInt(val("reserved_slots"), 10),
      critical_celltypes: [...document.querySelectorAll("#critical input:checked")].map((x) => x.value),
      genes_key: gk.value || null,
      normalize: document.getElementById("normalize").checked,
      preselected_genes: val("preselected").split(/[\s,;]+/).filter(Boolean),
      evaluation_scheme: val("evaluation_scheme"),
      evaluate_reference_sets: document.getElementById("evaluate_reference_sets").checked,
      panel_size_curve: document.getElementById("panel_size_curve").checked,
      seed: parseInt(val("seed"), 10),
      n_jobs: parseInt(val("n_jobs"), 10),
      feasibility,
    };
    const fd = new FormData();
    fd.append("upload_id", ds.id);
    fd.append("options", JSON.stringify(options));
    const markers = document.getElementById("marker_list").files[0];
    if (markers) fd.append("marker_list", markers);
    const btn = document.getElementById("submit");
    const msg = document.getElementById("form-msg");
    btn.disabled = true;
    msg.textContent = "";
    try {
      const job = await api("jobs", { method: "POST", body: fd });
      location.hash = `#/jobs/${job.id}`;
    } catch (err) {
      msg.innerHTML = `<p class="error">${esc(err.message)}</p>`;
      btn.disabled = false;
    }
  };
}

// One run ----------------------------------------------------------------------------------------------------------

async function showJob(id) {
  const job = await api(`jobs/${id}`);
  const o = job.options;
  const active = job.status === "queued" || job.status === "running";
  const stageIdx = STAGES.indexOf(job.stage);
  const steps = STAGES.filter((s) => (s !== "selecting_reference_sets" || o.evaluate_reference_sets)
      && (s !== "panel_size_curve" || o.panel_size_curve))
    .map((s) => {
      const i = STAGES.indexOf(s);
      const cls = job.status === "succeeded" || (stageIdx > i) ? "done" : (job.status === "running" && i === stageIdx ? "current" : "");
      return `<li class="${cls}">${esc(STAGE_LABEL[s])}</li>`;
    }).join("");
  const crit = o.critical_celltypes.length ? o.critical_celltypes.map(esc).join(", ") : "none";
  let html = `<p><a href="#/">← Runs</a></p><h1>${esc(job.input_filename)}</h1>
    <div>${badge(job.status, STATUS_STYLE[job.status])}
      <span class="muted">started ${esc(when(job.created_at))}${job.started_at ? ` · ran ${duration(job.started_at, job.finished_at)}` : ""}</span></div>
    <ul class="steps">${steps}</ul>
    <div class="card muted">Cell type column <strong>${esc(o.celltype_key)}</strong> · ${o.n} genes of a
      ${o.panel_capacity}-gene panel (${o.reserved_slots} reserved) · critical cell types: ${crit}</div>
    <div class="actions">
      ${active ? `<button class="button danger" data-act="cancel">Cancel</button>` : ""}
      ${job.status === "failed" || job.status === "cancelled" ? `<button class="button" data-act="retry">Retry (resumes)</button>` : ""}
      ${!active ? `<button class="button danger" data-act="delete">Delete</button>` : ""}
    </div>`;
  if (job.status === "failed") html += `<p class="error">${esc(job.error)}</p>`;
  if (job.status === "succeeded") {
    const v = job.verdict;
    if (v) {
      html += `<div class="verdict ${VERDICT_STYLE[v.verdict] || "unknown"}">
        <div class="label">${esc(v.verdict)}</div><div>${esc(v.reason)}</div></div>`;
    }
    html += `<div class="actions"><a class="button primary" href="jobs/${esc(id)}/report" target="_blank">Open report</a>
      <a class="button" href="jobs/${esc(id)}/report?download=1">Download report</a>
      <a class="button" href="jobs/${esc(id)}/files/verdict.json" download>verdict.json</a>
      <a class="button" href="jobs/${esc(id)}/files/probeset.csv" download>probeset.csv</a>
      <a class="button" href="jobs/${esc(id)}/files/evaluation_summary.csv" download>evaluation_summary.csv</a></div>
      <iframe class="report" src="jobs/${esc(id)}/report" title="Feasibility report"></iframe>`;
  }
  if (job.status !== "succeeded") {
    const log = await fetch(`jobs/${id}/log?lines=80`).then((r) => r.text());
    html += `<h2>Log</h2><pre class="log">${esc(log) || "(empty)"}</pre>`;
  }
  view.innerHTML = html;
  view.querySelectorAll("[data-act]").forEach((b) => b.addEventListener("click", async () => {
    const act = b.dataset.act;
    if (act === "delete" && !confirm("Delete this run and its results?")) return;
    try {
      if (act === "delete") { await api(`jobs/${id}`, { method: "DELETE" }); location.hash = "#/"; return; }
      await api(`jobs/${id}/${act}`, { method: "POST" });
      route();
    } catch (err) { alert(err.message); }
  }));
  const log = view.querySelector("pre.log");
  if (log) log.scrollTop = log.scrollHeight;
  if (active) timer = setTimeout(route, 3000);
}

// Router -----------------------------------------------------------------------------------------------------------

async function route() {
  clearTimeout(timer);
  const hash = location.hash.replace(/^#/, "") || "/";
  try {
    const m = hash.match(/^\/jobs\/([A-Za-z0-9]+)$/);
    if (m) await showJob(m[1]);
    else if (hash === "/new") await showNew();
    else await showRuns();
  } catch (err) {
    view.innerHTML = `<p class="error">${esc(err.message)}</p><p><a href="#/">Back to runs</a></p>`;
  }
}

window.addEventListener("hashchange", route);
route();
