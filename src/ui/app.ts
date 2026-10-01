// The page: problem form, Compute (in a Web Worker), result panel, plots, CSV download.

import type { ProblemInput, SolveResult } from "../solver/api.ts";
import { solutionJson } from "../solver/api.ts";
import { problemLoad, problemToml } from "../solver/problemfile.ts";
import { PORT_VERSION, SNAPSHOT_COMMIT, SNAPSHOT_REPO, SNAPSHOT_TAG, SNAPSHOT_VERSION } from "../solver/snapshot.ts";
import { refusalDetail } from "./domain.ts";
import { renderPlots } from "./plots.ts";
import { PRESETS } from "./presets.ts";

const VARS = ["ρ", "vx", "vy", "vz", "Bx", "By", "Bz", "p"] as const;
const HARD_TIMEOUT_MS = 60_000;
/**
 * Shipped homotopy budget. Julia's default is 5 s; the port is 5–10× slower on hard
 * homotopy cases (reference/BASELINE.md), so the page allows more (plan §4.2).
 */
const DEFAULT_TIME_LIMIT = 20;

/** Plain-language explanation of a retcode/reason (plan §7). */
const REASONS: Record<string, string> = {
  none: "",
  trivial: "Left and right states are equal.",
  quasi_euler: "",
  nonfinite: "An input value is not a finite number.",
  gamma: "γ must be larger than 1.",
  density: "Densities must be positive.",
  pressure: "Pressures must be positive.",
  bn_jump: "The normal field Bx must be the same on both sides.",
  extreme_ratio: "Density or pressure ratio beyond 1:10⁶, outside the supported range.",
  switch_on_off: "Transverse field too small (Bt/√p below 10⁻⁴ on the larger and 10⁻⁶ on the smaller side): switch-on/off waves are not supported in v0.1.",
  vacuum: "The solution would need (near-)vacuum generation, a limit of the regular waves.",
  fast_switch_off: "The solution needs a fast switch-off, a limit of the regular waves.",
  slow_limit: "The slow shock reaches its regular limit: a compound or intermediate wave would be needed (not in v0.1).",
  bt_small: "The transverse field of a middle state vanishes (Bt/√p < 10⁻⁸).",
  no_convergence: "The nonlinear solver did not converge from any start, nor by continuation.",
  assembly_failed: "The waves could not be assembled after convergence.",
  check: "The independent checks failed. This would be a bug; please report the problem.",
};

const SUP: Record<string, string> = { "-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };
/** 1e-10 → "10⁻¹⁰", 2.5e-8 → "2.5·10⁻⁸" (the notation of the domain panel). */
const pow10 = (v: number) => {
  if (v === 0) return "0";
  const [m, e] = v.toExponential().split("e") as [string, string];
  const exp = [...String(Number(e))].map((c) => SUP[c] ?? c).join("");
  return `${m === "1" ? "" : `${m}·`}10${exp}`;
};

const el = <T extends HTMLElement = HTMLElement>(id: string) => document.getElementById(id) as T;
const fmt = (v: number, d = 6) => (Number.isFinite(v) ? v.toFixed(d) : String(v));
const sci = (v: number) => (Number.isFinite(v) ? v.toExponential(2) : String(v));

export function startApp(): void {
  buildForm();
  el<HTMLSelectElement>("preset").addEventListener("change", () => loadPreset(el<HTMLSelectElement>("preset").value));
  // editing the problem after choosing an example: it is no longer that example
  el("problem").addEventListener("input", (e) => {
    const t = e.target as HTMLElement;
    if (!(t instanceof HTMLInputElement) || t.type === "file" || t.closest(".advanced")) return;
    clearSelection();
  });
  el("problem").addEventListener("submit", (e) => { e.preventDefault(); compute(); });
  el<HTMLInputElement>("file").addEventListener("change", loadFile);
  el<HTMLInputElement>("bn_euler").value = "1e-10";
  el<HTMLInputElement>("time_limit").value = String(DEFAULT_TIME_LIMIT);
  el<HTMLInputElement>("homotopy").checked = true;
  // the domain panel shows the quasi-Euler threshold currently set in the advanced options
  const showBnEuler = () => {
    const s = el<HTMLInputElement>("bn_euler").value.trim();
    const v = s === "" ? NaN : Number(s);
    el("domain-bn-euler").textContent = Number.isFinite(v) && v >= 0 ? pow10(v) : "bn_euler";
  };
  el("bn_euler").addEventListener("input", showBnEuler);
  showBnEuler();
  const link = document.createElement("a");
  link.href = `${SNAPSHOT_REPO}/releases/tag/${SNAPSHOT_TAG}`;
  link.textContent = `ExactMHDRiemannSolver.jl ${SNAPSHOT_VERSION} (${SNAPSHOT_COMMIT.slice(0, 7)})`;
  el("snapshot").replaceChildren(`TypeScript port ${PORT_VERSION} of `, link, ". Computed in your browser.");
  showIdle();
}

/** Start state: no result, a hint what to do. */
function showIdle(): void {
  const status = el("status");
  status.className = "status idle";
  status.textContent = "Choose an example, load a problem file or enter the two states, then press Compute.";
  for (const id of ["facts", "messages"]) el(id).replaceChildren();
  for (const id of ["waves", "downloads", "plots-section"]) el(id).hidden = true;
}

/** The fields no longer show a selected example (or loaded file). */
function clearSelection(): void {
  el<HTMLSelectElement>("preset").value = "";
  el("preset-note").textContent = "";
  const msg = el("file-msg");
  msg.className = "muted";
  msg.textContent = "";
}

function buildForm(): void {
  const sel = el<HTMLSelectElement>("preset");
  sel.add(new Option("Select an example…", "", true, true));
  for (const p of PRESETS) sel.add(new Option(p.label, p.id));
  const grid = el("states");
  for (const side of ["L", "R"] as const) {
    const fs = document.createElement("fieldset");
    fs.className = "state";
    const lg = document.createElement("legend");
    lg.textContent = side === "L" ? "Left state" : "Right state";
    fs.append(lg);
    VARS.forEach((v, i) => {
      const lab = document.createElement("label");
      const span = document.createElement("span");
      span.textContent = v;
      const inp = document.createElement("input");
      inp.id = `${side}${i}`;
      inp.inputMode = "decimal";
      inp.autocomplete = "off";
      inp.spellcheck = false;
      lab.append(span, inp);
      fs.append(lab);
    });
    grid.append(fs);
  }
}

function loadPreset(id: string): void {
  const p = PRESETS.find((q) => q.id === id);
  if (p === undefined) { clearSelection(); return; }
  el<HTMLSelectElement>("preset").value = p.id;
  el("preset-note").textContent = p.note;
  el("file-msg").textContent = "";
  fillForm(p.problem);
}

async function loadFile(): Promise<void> {
  const inp = el<HTMLInputElement>("file");
  const f = inp.files?.[0];
  const msg = el("file-msg");
  if (!f) return;
  try {
    const p = problemLoad(await f.text(), f.name);
    fillForm(p);
    clearSelection();
    msg.textContent = `Loaded ${f.name}. Press Compute to solve it.`;
  } catch (e) {
    msg.className = "muted bad";
    msg.textContent = e instanceof Error ? e.message : String(e);
  } finally {
    inp.value = "";
  }
}

function fillForm(pr: ProblemInput): void {
  for (const inp of el("problem").querySelectorAll<HTMLInputElement>("input[aria-invalid]")) inp.removeAttribute("aria-invalid");
  VARS.forEach((_, i) => {
    el<HTMLInputElement>(`L${i}`).value = String(pr.L[i]);
    el<HTMLInputElement>(`R${i}`).value = String(pr.R[i]);
  });
  el<HTMLInputElement>("gamma").value = String(pr.gamma);
  el<HTMLInputElement>("t").value = String(pr.t);
  el<HTMLInputElement>("xmin").value = String(pr.x[0]);
  el<HTMLInputElement>("xmax").value = String(pr.x[1]);
  el<HTMLInputElement>("n").value = String(pr.n);
}

function readNumber(id: string, ok: (v: number) => boolean = Number.isFinite): number | null {
  const inp = el<HTMLInputElement>(id);
  const s = inp.value.trim();
  const v = s === "" ? NaN : Number(s);
  const good = ok(v);
  inp.toggleAttribute("aria-invalid", !good);
  return good ? v : null;
}

function readProblem(): ProblemInput | null {
  const L = VARS.map((_, i) => readNumber(`L${i}`));
  const R = VARS.map((_, i) => readNumber(`R${i}`));
  const gamma = readNumber("gamma");
  const t = readNumber("t", (v) => v > 0 && Number.isFinite(v));
  const xmin = readNumber("xmin"), xmax = readNumber("xmax");
  const n = readNumber("n", (v) => Number.isInteger(v) && v >= 2 && v <= 100_000);
  const bnEuler = readNumber("bn_euler", (v) => v >= 0 && Number.isFinite(v));
  const timeLimit = readNumber("time_limit", (v) => v > 0);
  const gs = el<HTMLInputElement>("guess").value.trim();
  const guess = gs === "" ? undefined : gs.split(/[\s,;]+/).map(Number);
  const guessOk = guess === undefined || (guess.length === 5 && guess.every(Number.isFinite));
  el<HTMLInputElement>("guess").toggleAttribute("aria-invalid", !guessOk);
  const all = [...L, ...R, gamma, t, xmin, xmax, n, bnEuler, timeLimit];
  if (all.some((v) => v === null) || !guessOk) return null;
  if (xmin! >= xmax!) {
    el<HTMLInputElement>("xmax").setAttribute("aria-invalid", "");
    return null;
  }
  const opts = { bn_euler: bnEuler!, time_limit: timeLimit!, homotopy: el<HTMLInputElement>("homotopy").checked };
  return { L: L as number[], R: R as number[], gamma: gamma!, t: t!, x: [xmin!, xmax!], n: n!, opts,
    ...(guess ? { guess } : {}) };
}

// ---------------------------------------------------------------- worker

let worker: Worker | null = null;
let nextId = 1;

function newWorker(): Worker {
  return new Worker(new URL("../worker/solver.worker.ts", import.meta.url), { type: "module" });
}

function runInWorker(problem: ProblemInput): Promise<SolveResult> {
  worker ??= newWorker();
  const w = worker;
  const id = nextId++;
  return new Promise((resolve, reject) => {
    const timer = setTimeout(() => {
      w.terminate();
      worker = null;                                         // a fresh worker for the next solve
      reject(new Error(`no result after ${HARD_TIMEOUT_MS / 1000} s`));
    }, HARD_TIMEOUT_MS);
    w.onmessage = (ev: MessageEvent) => {
      if (ev.data.id !== id) return;
      clearTimeout(timer);
      if (ev.data.ok) resolve(ev.data.result as SolveResult);
      else reject(new Error(ev.data.error));
    };
    w.onerror = (ev) => { clearTimeout(timer); reject(new Error(ev.message)); };
    w.postMessage({ id, problem });
  });
}

let busy = false;

async function compute(): Promise<void> {
  if (busy) return;
  const problem = readProblem();
  const status = el("status");
  if (problem === null) {
    status.className = "status bad";
    status.textContent = "Please correct the highlighted fields. (Or try reloading the page.)";
    return;
  }
  busy = true;
  const btn = el<HTMLButtonElement>("compute");
  btn.disabled = true;
  status.className = "status";
  status.textContent = "Computing…";
  try {
    const r = await runInWorker(problem);
    showResult(problem, r);
  } catch (e) {
    status.className = "status bad";
    status.textContent = `The computation failed: ${e instanceof Error ? e.message : String(e)}`;
  } finally {
    busy = false;
    btn.disabled = false;
  }
}

// ---------------------------------------------------------------- result

let lastCsv = "";
let lastJson = "";
let lastToml = "";

function showResult(problem: ProblemInput, r: SolveResult): void {
  const status = el("status");
  const good = r.retcode === "Success";
  status.className = `status ${good ? "good" : r.retcode === "RegularLimit" || r.retcode === "Unsupported" ? "warn" : "bad"}`;
  const head = document.createElement("strong");
  head.textContent = r.retcode;
  const detail = [r.reason !== "none" ? r.reason : "", r.method !== "none" ? `method ${r.method}` : ""].filter(Boolean).join(" · ");
  status.replaceChildren(head, detail ? `  ${detail}` : "");
  const why = REASONS[r.reason] ?? "";
  const measured = refusalDetail(r.reason, problem.L, problem.R, problem.opts?.bn_euler);
  for (const text of [why, measured]) {
    if (!text) continue;
    const p = document.createElement("p");
    p.textContent = text;
    status.append(p);
  }

  const facts: [string, string][] = [
    ["residual", sci(r.residual)],
    ["worst check error", r.check ? sci(r.check.maxerr) : "–"],
    ["time", `${r.ms.toFixed(0)} ms`],
  ];
  el("facts").replaceChildren(...facts.map(([k, v]) => {
    const d = document.createElement("div");
    const dt = document.createElement("dt"); dt.textContent = k;
    const dd = document.createElement("dd"); dd.textContent = v;
    d.append(dt, dd);
    return d;
  }));
  const msgs = el("messages");
  msgs.replaceChildren(...(r.check?.messages ?? []).map((m) => { const li = document.createElement("li"); li.textContent = m; return li; }));

  // wave table
  const tbody = el("waves-body");
  const rows = [
    ["left state", "", "", ...problem.L.map((v) => fmt(v))],
    ...r.table.map((w) => [w.kind.replace("_", " "), fmt(w.s_left), fmt(w.s_right),
      fmt(w.rho), fmt(w.vx), fmt(w.vy), fmt(w.vz), fmt(w.Bx), fmt(w.By), fmt(w.Bz), fmt(w.p)]),
  ];
  tbody.replaceChildren(...rows.map((cells) => {
    const tr = document.createElement("tr");
    cells.forEach((c, j) => { const td = document.createElement(j === 0 ? "th" : "td"); td.textContent = c; tr.append(td); });
    return tr;
  }));
  el("waves").hidden = r.table.length === 0;

  lastCsv = r.csv;
  lastJson = solutionJson(problem, r);
  lastToml = problemToml(problem);
  el("downloads").hidden = false;
  el("dl-csv").hidden = r.csv === "";
  el("plots-section").hidden = r.x.length === 0;
  if (r.x.length > 0) {
    el("plots-caption").textContent = `Primitive variables at t = ${problem.t}, ${problem.n} points.`;
    renderPlots(el("plots"), r.x, r.W, [...VARS]);
  }
}

function save(text: string, name: string, type: string): void {
  if (text === "") return;
  const url = URL.createObjectURL(new Blob([text], { type }));
  const a = document.createElement("a");
  a.href = url;
  a.download = name;
  a.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

export function bindDownloads(): void {
  el("dl-csv").addEventListener("click", () => save(lastCsv, "solution.csv", "text/csv"));
  el("dl-json").addEventListener("click", () => save(lastJson, "solution.json", "application/json"));
  el("dl-toml").addEventListener("click", () => save(lastToml, "problem.toml", "application/toml"));
}
