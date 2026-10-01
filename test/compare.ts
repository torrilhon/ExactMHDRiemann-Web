// Comparison of TS solves with the golden Julia solves (plan §5.3).
//
// - retcode, reason, method must be identical;
// - wave tables: same kinds; speeds and states relerr ≤ max(1e-9, 1e4 × residual);
// - coplanar problems (canonical twist |α| = π): the mirror image z → -z of the canonical
//   solution is an equally valid solution, picked by round-off (reference/BASELINE.md);
//   they are compared on mirror-invariant quantities (ρ, vx, p, components of Bt and
//   vt - vt_R along the common field direction, magnitudes of the normal components);
// - Ψ only as a diagnostic.

import { canonicalStates, makeFrame, toHState } from "../src/solver/canonical.ts";
import { wavetable } from "../src/solver/output.ts";
import { solve } from "../src/solver/solve.ts";
import type { RiemannSolution, SolverOptions } from "../src/solver/types.ts";
import { riemannProblem } from "../src/solver/types.ts";

export interface GoldenCase {
  name?: string;
  seed?: number;
  index?: number;
  problem: { L: number[]; R: number[]; gamma: number };
  opts?: Record<string, unknown>;
  retcode: string;
  reason: string;
  method: string;
  psi: number[];
  residual: number;
  check: { ok: boolean; maxerr: number; messages: string[] } | null;
  table: (string | number)[][];
  time_ms: number;
}

export interface CaseResult {
  label: string;
  sol: RiemannSolution;
  ms: number;
  codeMatch: boolean;
  /** worst relerr of the wave table (NaN if the tables have different shape) */
  tableErr: number;
  tableTol: number;
  coplanar: boolean;
  psiErr: number;
  checkOk: boolean;
  messages: string[];
}

const OPT_KEYS = ["bn_euler", "bt_min_larger", "bt_min_smaller", "bt_floor", "ratio_max", "delta", "s_vac", "sat",
  "rho_floor", "abstol", "maxiters", "homotopy", "homotopy_maxsteps", "ode_tol", "check_tol"] as const;

function optsOf(c: GoldenCase): Partial<SolverOptions> {
  const o: Record<string, unknown> = { time_limit: Infinity };
  if (c.opts) {
    for (const k of OPT_KEYS) if (k in c.opts) o[k] = c.opts[k];
    if ("δ" in c.opts) o.delta = c.opts["δ"];
  }
  return o as Partial<SolverOptions>;
}

const rel = (a: number, b: number) => Math.abs(a - b) / Math.max(Math.abs(a), Math.abs(b), 1);

/** Twist angle of the canonical right state (coplanar if |α| = π). */
function canonicalAlpha(L: number[], R: number[]): number {
  const F = makeFrame(L, R);
  return toHState(canonicalStates(L, R, F)[1]).phi;
}

/** Mirror-invariant row: ρ, vx, p, parallel and |normal| parts of Bt and vt - vtR along e. */
function invariantRow(row: readonly number[], e: [number, number], vtR: [number, number]): number[] {
  const [, sl, sr, rho, vx, vy, vz, , By, Bz, p] = row as number[];
  const bpar = By! * e[0] + Bz! * e[1], bnor = Math.abs(-By! * e[1] + Bz! * e[0]);
  const wy = vy! - vtR[0], wz = vz! - vtR[1];
  const vpar = wy * e[0] + wz * e[1], vnor = Math.abs(-wy * e[1] + wz * e[0]);
  return [sl!, sr!, rho!, vx!, p!, bpar, bnor, vpar, vnor];
}

export function runCase(c: GoldenCase, label: string): CaseResult {
  const prob = riemannProblem(c.problem.L, c.problem.R, c.problem.gamma);
  const t0 = performance.now();
  const sol = solve(prob, { opts: optsOf(c) });
  const ms = performance.now() - t0;
  const codeMatch = sol.retcode === c.retcode && sol.reason === c.reason && sol.method === c.method;
  const rows = wavetable(sol).map((r) => [r.kind, r.s_left, r.s_right, r.rho, r.vx, r.vy, r.vz, r.Bx, r.By, r.Bz, r.p]);
  const tableTol = Math.max(1e-9, 1e4 * Math.max(Number.isFinite(c.residual) ? c.residual : 0,
    Number.isFinite(sol.residual) ? sol.residual : 0));
  let coplanar = false;
  let tableErr = 0;
  // a vanishing wave (ψ ~ 1e-17) is a shock or a fan by the sign of round-off: same
  // family and zero width count as the same wave
  const family = (k: unknown) => String(k).replace(/_(shock|fan)$/, "");
  const vanishing = (r: readonly unknown[]) => Math.abs((r[2] as number) - (r[1] as number)) <= 1e-9 * Math.max(1, Math.abs(r[1] as number));
  const sameKind = (a: readonly unknown[], b: readonly unknown[]) =>
    a[0] === b[0] || (family(a[0]) === family(b[0]) && vanishing(a) && vanishing(b));
  if (rows.length !== c.table.length || rows.some((r, i) => !sameKind(r, c.table[i]!))) {
    tableErr = NaN;
  } else if (rows.length > 0) {
    const direct = Math.max(0, ...rows.flatMap((r, i) => r.slice(1).map((v, j) => rel(v as number, c.table[i]![j + 1] as number))));
    tableErr = direct;
    if (direct > tableTol && sol.prob.L.every(Number.isFinite)) {
      const alpha = canonicalAlpha(c.problem.L, c.problem.R);
      if (Math.abs(alpha) > Math.PI - 1e-6) {
        coplanar = true;
        const L = c.problem.L, R = c.problem.R;
        const bl = Math.hypot(L[5]!, L[6]!);
        const e: [number, number] = bl > 0 ? [L[5]! / bl, L[6]! / bl] : [1, 0];
        const vtR: [number, number] = [R[2]!, R[3]!];
        tableErr = Math.max(0, ...rows.flatMap((r, i) => {
          const a = invariantRow(r as number[], e, vtR), b = invariantRow(c.table[i] as number[], e, vtR);
          return a.map((v, j) => rel(v, b[j]!));
        }));
      }
    }
  }
  const psiErr = Math.max(0, ...sol.psi.map((v, i) => {
    const g = c.psi[i]!;
    if (i === 2) {                                             // twist angle: modulo 2π
      const d = Math.abs(v - g) % (2 * Math.PI);
      return Math.min(d, 2 * Math.PI - d) / Math.max(1, Math.abs(g));
    }
    return rel(v, g);
  }));
  return { label, sol, ms, codeMatch, tableErr, tableTol, coplanar, psiErr,
    checkOk: sol.check?.ok ?? true, messages: sol.check?.messages ?? [] };
}

export function caseLabel(c: GoldenCase, i: number): string {
  return c.name ?? (c.seed !== undefined ? `seed ${c.seed}` : c.index !== undefined ? `#${c.index}` : `#${i}`);
}
