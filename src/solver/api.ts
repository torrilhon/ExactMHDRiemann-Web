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
  const sol = solve(riemannProblem(p.L, p.R, p.gamma), { opts: p.opts ?? {} });
  const table = wavetable(sol);
  const hasWaves = sol.waves.length > 0 && sol.frame !== null;
  const x = hasWaves ? juliaRange(p.x[0], p.x[1], p.n) : [];
  const W = hasWaves ? sampleGrid(sol, x, p.t) : [];
  const csv = hasWaves ? toCSV(sol, p.t, p.x[0], p.x[1], p.n) : "";
  return { retcode: sol.retcode, reason: sol.reason, method: sol.method, residual: sol.residual, psi: [...sol.psi],
    check: sol.check, table, x, W, csv, ms: performance.now() - t0 };
}
