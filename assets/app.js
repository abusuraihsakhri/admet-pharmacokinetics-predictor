"use strict";

const PYODIDE_INDEX = "https://cdn.jsdelivr.net/pyodide/v0.29.5/full/";
let pyodide = null;
let evaluatePython = null;
let simulatePython = null;

const runtimeStatus = document.getElementById("runtime-status");
const screenForm = document.getElementById("screen-form");
const pkForm = document.getElementById("pk-form");
const screenSubmit = screenForm.querySelector('button[type="submit"]');
const pkSubmit = pkForm.querySelector('button[type="submit"]');

function setRuntimeStatus(text, state = "") {
  runtimeStatus.textContent = text;
  runtimeStatus.className = "status" + (state ? " " + state : "");
}

function finiteNumber(formData, name, optional = false) {
  const raw = formData.get(name);
  if ((raw === null || raw === "") && optional) return null;
  const value = Number(raw);
  if (!Number.isFinite(value)) throw new Error(`${name} must be a finite number.`);
  return value;
}

async function initPython() {
  try {
    pyodide = await loadPyodide({ indexURL: PYODIDE_INDEX });
    const response = await fetch("./admet_predictor/__init__.py", { cache: "no-store" });
    if (!response.ok) throw new Error(`Could not load Python engine (HTTP ${response.status}).`);
    const source = await response.text();

    try { pyodide.FS.mkdir("/home/pyodide/admet_predictor"); } catch (_) {}
    pyodide.FS.writeFile("/home/pyodide/admet_predictor/__init__.py", source);

    await pyodide.runPythonAsync(`
import json
import sys
from dataclasses import asdict

if "/home/pyodide" not in sys.path:
    sys.path.insert(0, "/home/pyodide")

from admet_predictor import MoleculeProperties, ADMETPredictor, PharmacokineticSimulator

def web_evaluate(payload_json):
    payload = json.loads(payload_json)
    molecule = MoleculeProperties(**payload)
    return ADMETPredictor.evaluate_candidate(molecule).to_json()

def web_simulate(payload_json):
    payload = json.loads(payload_json)
    route = payload.pop("route")
    if route == "oral":
        result = PharmacokineticSimulator.simulate_oral_single(**payload)
    elif route == "iv":
        result = PharmacokineticSimulator.simulate_iv_bolus(**payload)
    elif route == "multi":
        result = PharmacokineticSimulator.simulate_oral_multiple(**payload)
    else:
        raise ValueError("Unsupported route")
    return json.dumps(asdict(result))
`);

    evaluatePython = pyodide.globals.get("web_evaluate");
    simulatePython = pyodide.globals.get("web_simulate");
    screenSubmit.disabled = false;
    pkSubmit.disabled = false;
    setRuntimeStatus("Python ready", "ready");
  } catch (error) {
    console.error(error);
    setRuntimeStatus("Runtime failed", "error");
    showError(document.getElementById("screen-results"), error);
    showError(document.getElementById("pk-results"), error);
  }
}

function setTheme(theme) {
  document.documentElement.dataset.theme = theme;
  localStorage.setItem("admet-theme", theme);
  document.getElementById("theme-toggle").setAttribute("aria-label", `Switch to ${theme === "dark" ? "light" : "dark"} theme`);
  redrawVisibleChart();
}

function initTheme() {
  const stored = localStorage.getItem("admet-theme");
  const preferred = window.matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light";
  setTheme(stored || preferred);
}

document.getElementById("theme-toggle").addEventListener("click", () => {
  setTheme(document.documentElement.dataset.theme === "dark" ? "light" : "dark");
});

const tabs = Array.from(document.querySelectorAll(".tab"));
tabs.forEach((tab, index) => {
  tab.addEventListener("click", () => activateTab(tab));
  tab.addEventListener("keydown", (event) => {
    if (!["ArrowLeft", "ArrowRight"].includes(event.key)) return;
    event.preventDefault();
    const delta = event.key === "ArrowRight" ? 1 : -1;
    const next = tabs[(index + delta + tabs.length) % tabs.length];
    next.focus();
    activateTab(next);
  });
});

function activateTab(tab) {
  tabs.forEach(item => {
    const selected = item === tab;
    item.classList.toggle("active", selected);
    item.setAttribute("aria-selected", String(selected));
    const panel = document.getElementById(item.dataset.tab);
    panel.hidden = !selected;
  });
  redrawVisibleChart();
}

document.getElementById("load-example").addEventListener("click", () => {
  const values = {
    name: "Aspirin",
    mw: 180.16,
    logp: 1.19,
    hbd: 1,
    hba: 3,
    tpsa: 63.6,
    rotatable_bonds: 3,
    aromatic_rings: 1,
    heavy_atoms: 13,
    molar_refractivity: 43.8,
    fsp3: 0.11,
    pka_base: "",
    pka_acid: 3.5,
  };
  Object.entries(values).forEach(([key, value]) => {
    const input = screenForm.elements.namedItem(key);
    if (input) input.value = value;
  });
});

screenForm.addEventListener("submit", async (event) => {
  event.preventDefault();
  if (!evaluatePython) return;

  const resultCard = screenForm.nextElementSibling;
  resultCard.setAttribute("aria-busy", "true");
  screenSubmit.disabled = true;
  screenSubmit.textContent = "Analyzing…";

  try {
    const data = new FormData(screenForm);
    const payload = {
      name: String(data.get("name") || "Candidate").trim() || "Candidate",
      mw: finiteNumber(data, "mw"),
      logp: finiteNumber(data, "logp"),
      hbd: Math.trunc(finiteNumber(data, "hbd")),
      hba: Math.trunc(finiteNumber(data, "hba")),
      tpsa: finiteNumber(data, "tpsa"),
      rotatable_bonds: Math.trunc(finiteNumber(data, "rotatable_bonds")),
      aromatic_rings: Math.trunc(finiteNumber(data, "aromatic_rings")),
      heavy_atoms: Math.trunc(finiteNumber(data, "heavy_atoms")),
      molar_refractivity: finiteNumber(data, "molar_refractivity"),
      fsp3: finiteNumber(data, "fsp3"),
      pka_base: finiteNumber(data, "pka_base", true),
      pka_acid: finiteNumber(data, "pka_acid", true),
    };
    const report = JSON.parse(evaluatePython(JSON.stringify(payload)));
    renderScreening(report);
  } catch (error) {
    showError(document.getElementById("screen-results"), error);
  } finally {
    resultCard.setAttribute("aria-busy", "false");
    screenSubmit.disabled = false;
    screenSubmit.textContent = "Analyze molecule";
  }
});

function el(tag, className, text) {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text !== undefined) node.textContent = text;
  return node;
}

function metric(label, value) {
  const node = el("div", "metric");
  node.append(el("span", "", label), el("strong", "", value));
  return node;
}

function renderScreening(report) {
  document.getElementById("screen-empty").classList.add("hidden");
  const root = document.getElementById("screen-results");
  root.classList.remove("hidden");
  root.replaceChildren();

  const heading = el("div", "result-title");
  const titleWrap = el("div");
  titleWrap.append(el("h2", "", report.molecule.name), el("p", "", report.overall_assessment));
  const score = el("div", "score-badge");
  score.append(el("strong", "", report.overall_druglikeness_score.toFixed(1)), el("span", "", "Rule score / 100"));
  heading.append(titleWrap, score);
  root.append(heading);

  const metrics = el("div", "metric-grid");
  metrics.append(
    metric("QED-like", report.qed.qed_score.toFixed(3)),
    metric("CNS MPO", report.cns_mpo.score.toFixed(2) + " / 6"),
    metric("HIA heuristic", report.admet_prediction.hia_pct.toFixed(1) + "%"),
    metric("PPB heuristic", report.admet_prediction.ppb_pct.toFixed(1) + "%"),
    metric("Vd heuristic", report.admet_prediction.vd_ss_l_kg.toFixed(2) + " L/kg"),
    metric("Half-life heuristic", report.admet_prediction.elimination_half_life_hr.toFixed(2) + " h"),
    metric("hERG flag", report.admet_prediction.herg_risk),
    metric("DILI flag", report.admet_prediction.dili_risk)
  );
  root.append(metrics);

  root.append(el("h3", "section-title", "Rule-based filters"));
  const table = el("table", "rule-table");
  const thead = el("thead");
  const hr = el("tr");
  ["Filter", "Result", "Violations"].forEach(label => hr.append(el("th", "", label)));
  thead.append(hr);
  table.append(thead);
  const tbody = el("tbody");
  [
    ["Lipinski", report.lipinski],
    ["Veber", report.veber],
    ["Egan", report.egan],
    ["Ghose", report.ghose],
    ["Muegge", report.muegge],
    ["Lead-likeness", report.lead_likeness],
  ].forEach(([name, result]) => {
    const row = el("tr");
    row.append(
      el("td", "", name),
      el("td", result.passes ? "pass" : "fail", result.passes ? "PASS" : "FAIL"),
      el("td", "", String(result.violations))
    );
    tbody.append(row);
  });
  table.append(tbody);
  root.append(table);

  if (report.recommendations.length) {
    root.append(el("h3", "section-title", "Screening notes"));
    const list = el("ul", "recommendations");
    report.recommendations.forEach(item => list.append(el("li", "", item)));
    root.append(list);
  }
}

function updateRouteFields() {
  const route = pkForm.elements.namedItem("route").value;
  document.querySelectorAll(".oral-field").forEach(node => node.classList.toggle("hidden", route === "iv"));
  document.querySelectorAll(".single-field").forEach(node => node.classList.toggle("hidden", route === "multi"));
  document.querySelectorAll(".multi-field").forEach(node => node.classList.toggle("hidden", route !== "multi"));
}

pkForm.querySelectorAll('input[name="route"]').forEach(input => input.addEventListener("change", updateRouteFields));
updateRouteFields();

pkForm.addEventListener("submit", async (event) => {
  event.preventDefault();
  if (!simulatePython) return;

  const resultCard = pkForm.nextElementSibling;
  resultCard.setAttribute("aria-busy", "true");
  pkSubmit.disabled = true;
  pkSubmit.textContent = "Simulating…";

  try {
    const data = new FormData(pkForm);
    const route = String(data.get("route"));
    let payload;
    if (route === "iv") {
      payload = {
        route,
        dose_mg: finiteNumber(data, "dose_mg"),
        ke_hr: finiteNumber(data, "ke_hr"),
        vd_l: finiteNumber(data, "vd_l"),
        duration_hr: finiteNumber(data, "duration_hr"),
      };
    } else if (route === "multi") {
      payload = {
        route,
        dose_mg: finiteNumber(data, "dose_mg"),
        bioavailability_f: finiteNumber(data, "bioavailability_f"),
        ka_hr: finiteNumber(data, "ka_hr"),
        ke_hr: finiteNumber(data, "ke_hr"),
        vd_l: finiteNumber(data, "vd_l"),
        dosing_interval_tau_hr: finiteNumber(data, "tau"),
        num_doses: Math.trunc(finiteNumber(data, "num_doses")),
      };
    } else {
      payload = {
        route,
        dose_mg: finiteNumber(data, "dose_mg"),
        bioavailability_f: finiteNumber(data, "bioavailability_f"),
        ka_hr: finiteNumber(data, "ka_hr"),
        ke_hr: finiteNumber(data, "ke_hr"),
        vd_l: finiteNumber(data, "vd_l"),
        duration_hr: finiteNumber(data, "duration_hr"),
      };
    }
    const result = JSON.parse(simulatePython(JSON.stringify(payload)));
    renderPk(result);
  } catch (error) {
    showError(document.getElementById("pk-results"), error);
  } finally {
    resultCard.setAttribute("aria-busy", "false");
    pkSubmit.disabled = false;
    pkSubmit.textContent = "Run simulation";
  }
});

let latestCurve = null;

function renderPk(result) {
  document.getElementById("pk-empty").classList.add("hidden");
  const root = document.getElementById("pk-results");
  root.classList.remove("hidden");
  root.replaceChildren();

  const heading = el("div", "result-title");
  const titleWrap = el("div");
  titleWrap.append(el("h2", "", result.dosing_route.replaceAll("_", " ")), el("p", "", `Dose ${result.dose_mg.toFixed(1)} mg · Vd ${result.volume_distribution_l.toFixed(2)} L`));
  heading.append(titleWrap);
  root.append(heading);

  const metrics = el("div", "metric-grid");
  metrics.append(
    metric("Cmax", result.cmax_mg_l.toFixed(4) + " mg/L"),
    metric("Tmax", result.tmax_hr.toFixed(2) + " h"),
    metric("Half-life", result.half_life_hr.toFixed(2) + " h"),
    metric("AUC 0–∞", result.auc_0_inf_mg_hr_l.toFixed(2) + " mg·h/L")
  );
  if (result.c_ss_avg_mg_l !== null) {
    metrics.append(
      metric("Css avg", result.c_ss_avg_mg_l.toFixed(4) + " mg/L"),
      metric("Css max", result.c_ss_max_mg_l.toFixed(4) + " mg/L"),
      metric("Css min", result.c_ss_min_mg_l.toFixed(4) + " mg/L"),
      metric("Peak accumulation", result.accumulation_ratio.toFixed(2) + "×")
    );
  }
  root.append(metrics);

  const chart = el("div", "chart-wrap");
  const canvas = document.createElement("canvas");
  canvas.id = "pk-chart";
  canvas.setAttribute("aria-label", "Plasma concentration versus time");
  canvas.setAttribute("role", "img");
  chart.append(canvas);
  root.append(chart);

  latestCurve = result.concentration_curve;
  drawCurve(canvas, latestCurve);
}

function showError(root, error) {
  const emptyId = root.id === "screen-results" ? "screen-empty" : "pk-empty";
  document.getElementById(emptyId).classList.add("hidden");
  root.classList.remove("hidden");
  root.replaceChildren();
  const box = el("div", "inline-note");
  box.textContent = error instanceof Error ? error.message : String(error);
  root.append(box);
}

function drawCurve(canvas, points) {
  if (!canvas || !points || !points.length) return;
  const rect = canvas.getBoundingClientRect();
  const dpr = Math.max(1, window.devicePixelRatio || 1);
  canvas.width = Math.max(320, Math.floor(rect.width * dpr));
  canvas.height = Math.floor(210 * dpr);
  const ctx = canvas.getContext("2d");
  ctx.scale(dpr, dpr);

  const width = canvas.width / dpr;
  const height = canvas.height / dpr;
  const pad = { left: 48, right: 14, top: 14, bottom: 32 };
  const plotW = width - pad.left - pad.right;
  const plotH = height - pad.top - pad.bottom;
  const maxX = Math.max(...points.map(p => p.time_hr), 1);
  const maxY = Math.max(...points.map(p => p.plasma_conc_mg_l), 1e-6);
  const styles = getComputedStyle(document.documentElement);
  const border = styles.getPropertyValue("--border").trim();
  const muted = styles.getPropertyValue("--muted").trim();
  const accent = styles.getPropertyValue("--accent").trim();

  ctx.clearRect(0, 0, width, height);
  ctx.strokeStyle = border;
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(pad.left, pad.top);
  ctx.lineTo(pad.left, pad.top + plotH);
  ctx.lineTo(pad.left + plotW, pad.top + plotH);
  ctx.stroke();

  ctx.fillStyle = muted;
  ctx.font = "11px system-ui";
  ctx.textAlign = "center";
  ctx.fillText("Time (h)", pad.left + plotW / 2, height - 7);
  ctx.save();
  ctx.translate(12, pad.top + plotH / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText("Concentration (mg/L)", 0, 0);
  ctx.restore();

  ctx.strokeStyle = accent;
  ctx.lineWidth = 2;
  ctx.beginPath();
  points.forEach((point, index) => {
    const x = pad.left + (point.time_hr / maxX) * plotW;
    const y = pad.top + plotH - (point.plasma_conc_mg_l / maxY) * plotH;
    if (index === 0) ctx.moveTo(x, y);
    else ctx.lineTo(x, y);
  });
  ctx.stroke();

  ctx.fillStyle = muted;
  ctx.textAlign = "left";
  ctx.fillText(maxY.toFixed(maxY < 1 ? 3 : 2), 4, pad.top + 4);
  ctx.textAlign = "right";
  ctx.fillText(maxX.toFixed(1), pad.left + plotW, pad.top + plotH + 17);
}

function redrawVisibleChart() {
  if (!latestCurve) return;
  const canvas = document.getElementById("pk-chart");
  if (canvas && !document.getElementById("pk-panel").hidden) drawCurve(canvas, latestCurve);
}

let resizeTimer = null;
window.addEventListener("resize", () => {
  window.clearTimeout(resizeTimer);
  resizeTimer = window.setTimeout(redrawVisibleChart, 100);
});

initTheme();
initPython();
