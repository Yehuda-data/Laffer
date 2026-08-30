"use strict";

const numericFields = [
  "tau_c", "tau_n", "tau_k", "eta", "phi", "theta", "delta", "kappa",
  "n_target", "k_y", "x_y", "R", "psi", "gamma", "debt_y", "g_y",
  "s_y", "external_balance_y", "other_waste_y", "grid_min", "grid_max", "grid_step"
];
const selectFields = ["closure", "calibration", "external_balance_convention", "kappa_mode"];
const colors = {blue: "#2166ac", red: "#b2182b", green: "#26734d", purple: "#7b3294", black: "#17212b"};
const plotConfig = {responsive: true, displaylogo: false, modeBarButtonsToRemove: ["lasso2d", "select2d"]};
const layoutBase = {
  margin: {l: 62, r: 24, t: 52, b: 55}, paper_bgcolor: "#fff", plot_bgcolor: "#fff",
  hovermode: "x unified", font: {family: "Segoe UI, Arial, sans-serif", color: "#243442"},
  xaxis: {gridcolor: "#e5eaee", zeroline: false}, yaxis: {gridcolor: "#e5eaee", zeroline: false},
  legend: {orientation: "h", y: -0.2}
};

const state = {
  currentName: "Custom", baseline: null, curve: null, sensitivity: null,
  comparison: null, scenarioA: null, scenarioB: null, requestId: 0, timer: null
};

const $ = id => document.getElementById(id);
const finite = value => typeof value === "number" && Number.isFinite(value);
const fmt = (value, digits = 4) => finite(value) ? value.toLocaleString(undefined, {maximumFractionDigits: digits}) : "N/A";
const pct = value => finite(value) ? `${(100 * value).toFixed(2)}%` : "N/A";

async function api(path, options = {}) {
  if (window.LafferStaticApi) return window.LafferStaticApi.request(path, options);
  const response = await fetch(path, {headers: {"Content-Type": "application/json"}, ...options});
  let data;
  try { data = await response.json(); } catch { data = {detail: response.statusText}; }
  if (!response.ok) throw new Error(typeof data.detail === "string" ? data.detail : JSON.stringify(data.detail));
  return data;
}

function readSpecification() {
  const specification = {name: state.currentName};
  selectFields.forEach(id => specification[id] = $(id).value);
  numericFields.forEach(id => {
    const value = $(id).value.trim();
    specification[id] = value === "" ? null : Number(value);
  });
  return specification;
}

function populate(specification) {
  state.currentName = specification.name || "Custom";
  selectFields.forEach(id => { if (specification[id] != null) $(id).value = specification[id]; });
  numericFields.forEach(id => { $(id).value = specification[id] ?? ""; });
  updateVisibility();
}

function updateVisibility() {
  const closure = $("closure").value;
  const implied = $("calibration").value === "model_implied";
  const directKappa = $("kappa_mode").value === "kappa";
  $("g_y").disabled = closure !== "s_laffer";
  $("s_y").disabled = closure !== "g_laffer";
  $("gy-wrap").style.opacity = closure === "s_laffer" ? "1" : ".48";
  $("sy-wrap").style.opacity = closure === "g_laffer" ? "1" : ".48";
  $("closure-help").textContent = closure === "s_laffer"
    ? "Government spending g is fixed in levels; transfers s are endogenous."
    : "Transfers s are fixed in levels; government spending g is endogenous.";
  ["theta", "delta"].forEach(id => { $(id).disabled = implied; $(id + "-wrap").style.opacity = implied ? ".48" : "1"; });
  ["k_y", "x_y"].forEach(id => { $(id).disabled = !implied; $(id.replace("_", "") + "-wrap"); });
  $("ky-wrap").style.opacity = implied ? "1" : ".48";
  $("xy-wrap").style.opacity = implied ? "1" : ".48";
  $("kappa").disabled = !directKappa;
  $("n_target").disabled = directKappa;
  $("kappa-wrap").style.opacity = directKappa ? "1" : ".48";
  $("n-target-wrap").style.opacity = directKappa ? ".48" : "1";
  const isImports = $("external_balance_convention").value === "net_imports";
  $("balance-value-label").firstChild.textContent = isImports ? "Net imports m/y " : "Trade balance tb/y ";
  const value = Number($("external_balance_y").value || 0);
  $("balance-conversion").textContent = isImports
    ? `Backend conversion: tb/y = −m/y = ${fmt(-value, 5)}`
    : `Visible conversion: m/y = −tb/y = ${fmt(-value, 5)}`;
  updateSensitivityOptions();
  if (closure === "g_laffer" && $("tax-chart-kind").value === "capital") $("tax-chart-kind").value = "labor";
}

function updateSensitivityOptions() {
  const implied = $("calibration").value === "model_implied";
  const closure = $("closure").value;
  const options = ["R", "psi", "eta", "phi", "kappa", "m_y"];
  if (implied) options.unshift("k_y", "x_y"); else options.unshift("theta", "delta");
  options.push(closure === "s_laffer" ? "g_y" : "s_y");
  const selected = $("sensitivity-parameter").value;
  $("sensitivity-parameter").innerHTML = options.map(item => `<option value="${item}">${item}</option>`).join("");
  if (options.includes(selected)) $("sensitivity-parameter").value = selected;
  suggestSensitivityRange();
}

function suggestSensitivityRange() {
  const parameter = $("sensitivity-parameter").value;
  let field = parameter;
  if (parameter === "m_y") field = "external_balance_y";
  const current = Number($(field)?.value || 0.1);
  const spread = Math.max(Math.abs(current) * 0.25, parameter.includes("_y") ? 0.02 : 0.05);
  $("sensitivity-min").value = (current - spread).toPrecision(5);
  $("sensitivity-max").value = (current + spread).toPrecision(5);
}

function showError(message = "") {
  $("error-banner").textContent = message;
  $("error-banner").classList.toggle("hidden", !message);
}

function statusBadge(status) {
  const css = String(status).toLowerCase().replaceAll(" ", "-");
  return `<span class="badge ${css}">${status}</span>`;
}

function renderBaseline(data) {
  state.baseline = data;
  const eq = data.equilibrium || {};
  const valid = data.validity?.valid;
  $("baseline-subtitle").textContent = `${state.currentName} · ${$("closure").value.replace("_", "-")} · ${valid ? "valid equilibrium" : "invalid equilibrium"}`;
  const cards = [
    ["Labor n", eq.n], ["Output y", eq.y], ["Capital/output", eq.k_y],
    ["Consumption/output", eq.c_y], ["Tax revenue/output", eq.T_total_y]
  ];
  $("summary-cards").innerHTML = cards.map(([label, value]) =>
    `<div class="card ${valid ? "" : "invalid"}"><div class="label">${label}</div><div class="value">${fmt(value)}</div></div>`
  ).join("");

  const parameterRows = Object.entries(data.parameters || {}).map(([name, item]) => [name, item.value, item.status]);
  const macroKeys = ["n", "y", "k_y", "x_y", "c_y", "g_y", "s_y", "m_y", "w", "d",
    "labor_tax_base_y", "capital_tax_base_y", "consumption_tax_base_y", "T_n_y", "T_k_y", "T_c_y", "T_total_y"];
  const outputRows = macroKeys.map(name => [name, eq[name], "EQUILIBRIUM OUTPUT"]);
  $("baseline-table").innerHTML = `<thead><tr><th>Variable</th><th>Value</th><th>Status</th></tr></thead><tbody>${
    [...parameterRows, ...outputRows].map(([name, value, status]) => `<tr><td>${name}</td><td>${fmt(value, 7)}</td><td>${statusBadge(status)}</td></tr>`).join("")
  }</tbody>`;
  renderDiagnostics(data.diagnostics || []);
}

function renderDiagnostics(items) {
  $("diagnostics").innerHTML = items.map(item =>
    `<div class="diagnostic ${item.level.toLowerCase().replaceAll(" ", "-")}"><strong>${item.level}</strong>${item.message}</div>`
  ).join("");
}

function chartValue(point, key) { return point.valid && finite(point[key]) ? point[key] : null; }
function trace(points, xKey, yKey, name, color, dash = "solid") {
  return {x: points.map(p => p[xKey]), y: points.map(p => chartValue(p, yKey)), name, mode: "lines", line: {color, width: 2.2, dash}, connectgaps: false};
}

function renderCurve(data, kind = "labor") {
  state.curve = data;
  const points = data.curve || [];
  const xKey = kind === "labor" ? "tau_n" : "tau_k";
  const titleTax = kind === "labor" ? "Labor" : "Capital";
  const summary = data.summary || {};
  const base = summary.baseline_tax;
  const peak = summary.peak_tax;
  const peakRevenue = summary.peak_revenue;
  const annotations = [];
  if (finite(peak) && finite(peakRevenue)) annotations.push({x: peak, y: peakRevenue, text: `Peak ${fmt(peak, 3)}<br>${fmt(peakRevenue, 2)}`, showarrow: true, arrowhead: 2});
  const xTitle = kind === "labor" ? "τₙ" : "τₖ";
  Plotly.react("laffer-chart", [trace(points, xKey, "T_total_index", "Total tax revenue", colors.blue)], {
    ...layoutBase, title: {text: `${titleTax}-tax Laffer curve`}, xaxis: {...layoutBase.xaxis, title: {text: xTitle}},
    yaxis: {...layoutBase.yaxis, title: {text: "Baseline total revenue = 100"}}, annotations,
    shapes: finite(base) ? [{type: "line", x0: base, x1: base, y0: 0, y1: 1, yref: "paper", line: {color: "#555", dash: "dot"}}] : []
  }, plotConfig);
  Plotly.react("revenue-chart", [
    trace(points, xKey, "T_n_index", "Tₙ", colors.blue), trace(points, xKey, "T_k_index", "Tₖ", colors.red, "dash"),
    trace(points, xKey, "T_c_index", "T꜀", colors.green, "dot"), trace(points, xKey, "T_total_index", "Total", colors.black)
  ], {...layoutBase, title: {text: "Revenue decomposition"}, xaxis: {...layoutBase.xaxis, title: {text: xTitle}}, yaxis: {...layoutBase.yaxis, title: {text: "Total baseline revenue = 100"}}}, plotConfig);
  Plotly.react("macro-chart", [
    trace(points, xKey, "n_index", "n", colors.blue), trace(points, xKey, "y_index", "y", colors.red, "dash"),
    trace(points, xKey, "k_index", "k", colors.purple, "dot"), trace(points, xKey, "c_index", "c", colors.green, "dashdot")
  ], {...layoutBase, title: {text: "Macro response"}, xaxis: {...layoutBase.xaxis, title: {text: xTitle}}, yaxis: {...layoutBase.yaxis, title: {text: "Own baseline = 100"}}}, plotConfig);
  Plotly.react("fiscal-chart", [
    trace(points, xKey, "g_y", "g/y", colors.blue), trace(points, xKey, "s_y", "s/y", colors.red, "dash"), trace(points, xKey, "T_total_y", "T/y", colors.green, "dot")
  ], {...layoutBase, title: {text: "Fiscal response"}, xaxis: {...layoutBase.xaxis, title: {text: xTitle}}, yaxis: {...layoutBase.yaxis, title: {text: "Ratio to output"}}}, plotConfig);
  const laborDecompositionChart = $("labor-decomposition-chart");
  laborDecompositionChart.classList.toggle("hidden", kind !== "labor");
  if (kind === "labor") {
    Plotly.react("labor-decomposition-chart", [
      trace(points, "tau_n", "n_index", "n — hours worked", colors.blue),
      trace(points, "tau_n", "w_index", "w — wage", colors.purple, "dash"),
      {...trace(points, "tau_n", "T_n_own_index", "Tₙ — total labor-income tax revenue", colors.red), line: {color: colors.red, width: 3.2}}
    ], {...layoutBase, title: {text: "Labor-income tax revenue decomposition: Tₙ = τₙ · w · n"},
      xaxis: {...layoutBase.xaxis, title: {text: "Labor-income tax rate τₙ"}},
      yaxis: {...layoutBase.yaxis, title: {text: "Each series at baseline = 100"}}}, plotConfig);
  }

  const diagnostics = [
    ["Valid τ range", summary.valid_tau_min == null ? "N/A" : `${fmt(summary.valid_tau_min, 3)}–${fmt(summary.valid_tau_max, 3)}`],
    ["First g ≤ 0", fmt(summary.first_g_nonpositive, 3)], ["First n ≥ 1", fmt(summary.first_n_at_least_one, 3)],
    ["First c ≤ 0", fmt(summary.first_c_nonpositive, 3)], ["Solver failures", summary.solver_failure_rates?.length || 0],
    ["Possible branch switches", summary.possible_branch_switches?.length || 0]
  ];
  $("curve-diagnostics").innerHTML = `<table><tbody>${diagnostics.map(([a, b]) => `<tr><td>${a}</td><td>${b}</td></tr>`).join("")}</tbody></table>`;
}

function renderEquations(data) {
  $("equations").innerHTML = (data.equations || []).map(item =>
    `<div class="equation"><strong>${item.name}</strong><code>${item.latex}</code><small>${item.source}</small></div>`
  ).join("");
}

async function recompute() {
  const id = ++state.requestId;
  showError();
  $("recompute").disabled = true;
  $("api-status").textContent = "Computing…";
  try {
    const spec = readSpecification();
    const [base, curve, equationData] = await Promise.all([
      api("/api/baseline", {method: "POST", body: JSON.stringify(spec)}),
      api("/api/laffer/labor", {method: "POST", body: JSON.stringify(spec)}),
      api(`/api/equations/${spec.closure}`)
    ]);
    if (id !== state.requestId) return;
    renderBaseline(base);
    renderCurve(curve, "labor");
    renderEquations(equationData);
    $("tax-chart-kind").value = "labor";
    $("unsupported-capital").classList.add("hidden");
    $("api-status").textContent = "Model ready";
    $("api-status").className = "status-chip online";
  } catch (error) {
    showError(error.message);
    $("api-status").textContent = "Request failed";
    $("api-status").className = "status-chip offline";
  } finally { $("recompute").disabled = false; }
}

function scheduleRecompute() {
  clearTimeout(state.timer);
  state.timer = setTimeout(recompute, 650);
}

async function switchTaxChart() {
  const kind = $("tax-chart-kind").value;
  if (kind === "labor") { if (state.curve) renderCurve(state.curve, "labor"); return; }
  if ($("closure").value === "g_laffer") {
    $("unsupported-capital").classList.remove("hidden");
    $("tax-chart-kind").value = "labor";
    return;
  }
  $("unsupported-capital").classList.add("hidden");
  try {
    const data = await api("/api/laffer/capital", {method: "POST", body: JSON.stringify(readSpecification())});
    renderCurve(data, "capital");
  } catch (error) { showError(error.message); }
}

async function runSensitivity() {
  showError();
  const request = {
    specification: readSpecification(), parameter: $("sensitivity-parameter").value,
    minimum: Number($("sensitivity-min").value), maximum: Number($("sensitivity-max").value),
    scenarios: Number($("sensitivity-count").value)
  };
  try {
    const data = await api("/api/sensitivity", {method: "POST", body: JSON.stringify(request)});
    state.sensitivity = data;
    const traces = data.scenarios.map((scenario, index) => {
      const points = scenario.result.curve || [];
      return {...trace(points, "tau_n", "T_total_index", scenario.label, [colors.blue, colors.red, colors.green, colors.purple, "#e08214"][index % 5]), connectgaps: false};
    });
    Plotly.react("sensitivity-chart", traces, {...layoutBase, title: {text: `Sensitivity: ${data.parameter}`}, xaxis: {...layoutBase.xaxis, title: {text: "τₙ"}}, yaxis: {...layoutBase.yaxis, title: {text: "Scenario baseline revenue = 100"}}}, plotConfig);
  } catch (error) { showError(error.message); }
}

function saveScenario(which) {
  if (which === "A") state.scenarioA = readSpecification(); else state.scenarioB = readSpecification();
  $("scenario-status").textContent = `Scenario A: ${state.scenarioA ? "saved" : "not saved"}; Scenario B: ${state.scenarioB ? "saved" : "not saved"}.`;
}

async function runComparison() {
  if (!state.scenarioA || !state.scenarioB) { showError("Save both Scenario A and Scenario B before comparing."); return; }
  try {
    const data = await api("/api/compare", {method: "POST", body: JSON.stringify({scenario_a: state.scenarioA, scenario_b: state.scenarioB})});
    state.comparison = data;
    const a = data.scenario_a.curve.curve || [], b = data.scenario_b.curve.curve || [];
    Plotly.react("comparison-chart", [trace(a, "tau_n", "T_total_index", "Scenario A", colors.blue), trace(b, "tau_n", "T_total_index", "Scenario B", colors.red, "dash")],
      {...layoutBase, title: {text: "Scenario A vs B"}, xaxis: {...layoutBase.xaxis, title: {text: "τₙ"}}, yaxis: {...layoutBase.yaxis, title: {text: "Own baseline revenue = 100"}}}, plotConfig);
    const rowsA = data.scenario_a.key_rates, rowsB = data.scenario_b.key_rates;
    $("comparison-table").innerHTML = `<thead><tr><th>τₙ</th><th>A: n</th><th>B: n</th><th>A: T index</th><th>B: T index</th></tr></thead><tbody>${rowsA.map((row, index) => `<tr><td>${row.tau_n.toFixed(2)}</td><td>${fmt(row.n)}</td><td>${fmt(rowsB[index].n)}</td><td>${fmt(row.T_total_index)}</td><td>${fmt(rowsB[index].T_total_index)}</td></tr>`).join("")}</tbody>`;
  } catch (error) { showError(error.message); }
}

function flatten(object, prefix = "", target = {}) {
  Object.entries(object || {}).forEach(([key, value]) => {
    const name = prefix ? `${prefix}.${key}` : key;
    if (value && typeof value === "object" && !Array.isArray(value)) flatten(value, name, target);
    else if (!Array.isArray(value)) target[name] = value;
  });
  return target;
}

function downloadCsv(filename, rows) {
  if (!rows?.length) { showError("No results are available for this export."); return; }
  const keys = [...new Set(rows.flatMap(row => Object.keys(row)))];
  const escape = value => `"${String(value ?? "").replaceAll('"', '""')}"`;
  const csv = [keys.map(escape).join(","), ...rows.map(row => keys.map(key => escape(row[key])).join(","))].join("\r\n");
  const link = document.createElement("a");
  link.href = URL.createObjectURL(new Blob(["\ufeff", csv], {type: "text/csv;charset=utf-8"}));
  link.download = filename;
  link.click();
  URL.revokeObjectURL(link.href);
}

function exportData(kind) {
  const slug = `${state.currentName}-${$("closure").value}-${$("calibration").value}`.replaceAll(/[^a-zA-Z0-9_-]/g, "_");
  if (kind === "baseline") downloadCsv(`${slug}-baseline.csv`, [flatten(state.baseline)]);
  if (kind === "curve") downloadCsv(`${slug}-curve.csv`, state.curve?.curve || []);
  if (kind === "sensitivity") {
    const rows = (state.sensitivity?.scenarios || []).flatMap(s => (s.result.curve || []).map(row => ({scenario: s.label, scenario_value: s.value, ...row})));
    downloadCsv(`${slug}-sensitivity.csv`, rows);
  }
  if (kind === "comparison") {
    const rows = ["scenario_a", "scenario_b"].flatMap(key => (state.comparison?.[key]?.curve?.curve || []).map(row => ({scenario: key, ...row})));
    downloadCsv(`${slug}-comparison.csv`, rows);
  }
}

async function loadPresets() {
  try {
    const data = await api("/api/presets");
    $("preset-buttons").innerHTML = data.presets.map(item => `<button class="secondary" data-preset="${item.id}">${item.label}</button>`).join("");
    data.presets.forEach(item => document.querySelector(`[data-preset="${item.id}"]`).addEventListener("click", () => { populate(item.specification); recompute(); }));
    populate(data.presets[0].specification);
    await recompute();
  } catch (error) { showError(error.message); $("api-status").textContent = "Backend unavailable"; $("api-status").className = "status-chip offline"; }
}

function bindEvents() {
  ["closure", "calibration", "external_balance_convention", "kappa_mode"].forEach(id => $(id).addEventListener("change", () => { updateVisibility(); scheduleRecompute(); }));
  numericFields.forEach(id => $(id).addEventListener("change", () => { if (id === "external_balance_y") updateVisibility(); scheduleRecompute(); }));
  $("recompute").addEventListener("click", recompute);
  $("tax-chart-kind").addEventListener("change", switchTaxChart);
  $("sensitivity-parameter").addEventListener("change", suggestSensitivityRange);
  $("run-sensitivity").addEventListener("click", runSensitivity);
  $("save-a").addEventListener("click", () => saveScenario("A"));
  $("save-b").addEventListener("click", () => saveScenario("B"));
  $("compare").addEventListener("click", runComparison);
  $("custom-preset").addEventListener("click", () => { state.currentName = "Custom"; $("baseline-subtitle").textContent = "Custom editable scenario"; });
  document.querySelectorAll("[data-export]").forEach(button => button.addEventListener("click", () => exportData(button.dataset.export)));
  $("download-chart").addEventListener("click", () => Plotly.downloadImage("laffer-chart", {format: "png", filename: "laffer-curve", width: 1200, height: 700}));
}

bindEvents();
loadPresets();
