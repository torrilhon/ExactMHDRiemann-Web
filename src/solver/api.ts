// What the page asks the solver for and gets back (plain data, structured-clonable,
// so it crosses the Worker boundary).

import { toCSV, sampleGrid, wavetable, type WaveRow } from "./output.ts";
import { juliaRange } from "./format.ts";
import { solve } from "./solve.ts";
import type { CheckReport, RetCode, SolverOptions } from "./types.ts";
import { riemannProblem } from "./types.ts";

export interface ProblemInput {
  L: number[];
  R: number[];
  gamma: number;
  /** output time and grid (problem file `[output]`) */
  t: number;
  x: [number, number];
  n: number;
  opts?: Partial<SolverOptions>;
  /** start vector Ψ (canonical frame), tried first */
  guess?: number[];
}

export interface SolveResult {
  retcode: RetCode;
  reason: string;
  method: string;
  residual: number;
  psi: number[];
  check: CheckReport | null;
  table: WaveRow[];
  /** sampled profile: x and the 8 primitive variables (empty if there are no waves) */
  x: number[];
  W: number[][];
  csv: string;
  ms: number;
}

export function runProblem(p: ProblemInput): SolveResult {
  const t0 = performance.now();
  const sol = solve(riemannProblem(p.L, p.R, p.gamma), { opts: p.opts ?? {}, ...(p.guess ? { guess: p.guess } : {}) });
  const table = wavetable(sol);
  const hasWaves = sol.waves.length > 0 && sol.frame !== null;
  const x = hasWaves ? juliaRange(p.x[0], p.x[1], p.n) : [];
  const W = hasWaves ? sampleGrid(sol, x, p.t) : [];
  const csv = hasWaves ? toCSV(sol, p.t, p.x[0], p.x[1], p.n) : "";
  return { retcode: sol.retcode, reason: sol.reason, method: sol.method, residual: sol.residual, psi: [...sol.psi],
    check: sol.check, table, x, W, csv, ms: performance.now() - t0 };
}

import { PORT_VERSION, SNAPSHOT_COMMIT, SNAPSHOT_REPO, SNAPSHOT_VERSION } from "./snapshot.ts";
import { defaultOptions } from "./types.ts";

/** JSON has no Inf/NaN: written as the strings "Inf", "-Inf", "NaN" (plan §4.7). */
function nonfinite(_k: string, v: unknown): unknown {
  if (typeof v === "number" && !Number.isFinite(v)) return Number.isNaN(v) ? "NaN" : v > 0 ? "Inf" : "-Inf";
  return v;
}

/** solution.json: input, options, result and provenance. */
export function solutionJson(p: ProblemInput, r: SolveResult): string {
  const doc = {
    format: "ExactMHDRiemann-Web solution, v1",
    port_version: PORT_VERSION,
    snapshot: { repository: SNAPSHOT_REPO, commit: SNAPSHOT_COMMIT, version: SNAPSHOT_VERSION },
    input: { gamma: p.gamma, left: p.L, right: p.R, output: { t: p.t, x: p.x, n: p.n } },
    options: defaultOptions(p.opts ?? {}),
    retcode: r.retcode,
    reason: r.reason,
    method: r.method,
    residual: r.residual,
    psi: r.psi,
    waves: r.table,
    check: r.check,
  };
  return JSON.stringify(doc, nonfinite, 2) + "\n";
}
