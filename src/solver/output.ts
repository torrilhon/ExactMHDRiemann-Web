// Sampling of the self-similar solution and tabular output, all in the user frame.
// Mirrors output.jl.

import { findZeroBrent } from "./brent.ts";
import { speedToUser } from "./canonical.ts";
import { userState } from "./check.ts";
import { fanSpeed, fanState } from "./fans.ts";
import { juliaRange, juliaRepr } from "./format.ts";
import { ArgumentError, sqrtD } from "./math.ts";
import type { RiemannSolution, WaveKind } from "./types.ts";

/**
 * Primitive state (ρ, vx, vy, vz, Bx, By, Bz, p) at similarity coordinate ξ = x/t.
 * Fan interiors are evaluated exactly from the dense ODE output.
 */
export function sampleXi(sol: RiemannSolution, xi: number): number[] {
  if (!Number.isFinite(xi)) throw new ArgumentError("sample needs a finite ξ = x/t");
  if (sol.waves.length === 0 || sol.frame === null || sol.ctx === null) {
    throw new ArgumentError(`solution has no waves (retcode ${sol.retcode})`);
  }
  const F = sol.frame, ctx = sol.ctx;
  const xic = (F.mirror ? -xi : xi) / sqrtD(F.p0 / F.rho0);    // canonical speed
  for (const w of sol.waves) {
    if (xic < w.s_left) return userState(w.left, ctx, F);
    if (w.fan !== null && xic <= w.s_right) {
      const f = w.fan;
      const tau = findZeroBrent((t) => fanSpeed(f, t, ctx) - xic, 0.0, 1.0, { xatol: 1e-14 });
      return userState(fanState(f, tau, ctx), ctx, F);
    }
  }
  return userState(sol.waves[sol.waves.length - 1]!.right, ctx, F);
}

/** State at position x and time t > 0. */
export function sampleXT(sol: RiemannSolution, x: number, t: number): number[] {
  if (!(t > 0)) throw new ArgumentError("sample needs t > 0");
  return sampleXi(sol, x / t);
}

/** States at the positions x at time t > 0 (rows of 8). */
export function sampleGrid(sol: RiemannSolution, x: readonly number[], t: number): number[][] {
  if (!(t > 0)) throw new ArgumentError("sample needs t > 0");
  return x.map((xi) => sampleXi(sol, xi / t));
}

export interface WaveRow {
  kind: WaveKind;
  s_left: number;
  s_right: number;
  rho: number; vx: number; vy: number; vz: number; Bx: number; By: number; Bz: number; p: number;
}

/** One row per wave, left to right: kind, speeds (user frame) and the state right of the wave. */
export function wavetable(sol: RiemannSolution): WaveRow[] {
  const F = sol.frame, ctx = sol.ctx;
  if (F === null || ctx === null) return [];
  // in a mirrored frame the canonical waves run right to left in the user frame
  const ws = F.mirror ? [...sol.waves].reverse() : sol.waves;
  return ws.map((w) => {
    const W = userState(F.mirror ? w.left : w.right, ctx, F);
    const a = speedToUser(w.s_left, F), b = speedToUser(w.s_right, F);
    return { kind: w.kind, s_left: Math.min(a, b), s_right: Math.max(a, b),
      rho: W[0]!, vx: W[1]!, vy: W[2]!, vz: W[3]!, Bx: W[4]!, By: W[5]!, Bz: W[6]!, p: W[7]! };
  });
}

export const CSV_HEADER = "x,rho,vx,vy,vz,Bx,By,Bz,p";

/**
 * The text of Julia's `write_csv(path, sol; t, x = range(xmin, xmax; length = n))`:
 * header and one line per point, numbers as Julia's `repr`.
 */
export function toCSV(sol: RiemannSolution, t = 1.0, xmin = -1.0, xmax = 1.0, n = 2001): string {
  const x = juliaRange(xmin, xmax, n);
  const M = sampleGrid(sol, x, t);
  const lines = [CSV_HEADER];
  x.forEach((xi, i) => lines.push([xi, ...M[i]!].map(juliaRepr).join(",")));
  return lines.join("\n") + "\n";
}
