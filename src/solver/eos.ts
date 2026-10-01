// Characteristic speeds, conserved variables and fluxes of 1D ideal MHD, γ-law gas.
// Mirrors eos.jl.

import { hypot, sqrtD } from "./math.ts";
import type { Ctx, HState, Prim } from "./types.ts";

/**
 * Julia `speeds(ρ, p, bt, Bn, γ) -> (cf, cA, cs)`. cs is evaluated as a·cA/cf, which
 * avoids the cancellation of the textbook formula.
 */
export function speeds(rho: number, p: number, bt: number, Bn: number, gamma: number): [number, number, number] {
  const a2 = (gamma * p) / rho;
  const ca2 = (Bn * Bn) / rho;
  const b2 = ca2 + (bt * bt) / rho;
  const ab = a2 + b2;
  const d = sqrtD(Math.max(ab * ab - 4 * a2 * ca2, 0));
  const cf2 = (a2 + b2 + d) / 2;
  const cf = sqrtD(cf2);
  const cs = sqrtD((a2 * ca2) / cf2);
  return [cf, sqrtD(ca2), cs];
}

/** Julia `speeds(h::HState, ctx)`. */
export function speedsH(h: HState, ctx: Ctx): [number, number, number] {
  return speeds(h.rho, h.p, h.bt, ctx.Bn, ctx.gamma);
}

/** Julia `speeds(W::AbstractVector, γ)`: speeds of a primitive state. */
export function speedsW(W: Prim, gamma: number): [number, number, number] {
  return speeds(W[0]!, W[7]!, hypot(W[5]!, W[6]!), W[4]!, gamma);
}

/** Conserved variables (ρ, ρvx, ρvy, ρvz, By, Bz, E) of a primitive state. */
export function conserved(W: Prim, gamma: number): number[] {
  const [rho, u, v, w, bx, by, bz, p] = W as [number, number, number, number, number, number, number, number];
  const E = p / (gamma - 1) + (rho * (u * u + v * v + w * w)) / 2 + (bx * bx + by * by + bz * bz) / 2;
  return [rho, rho * u, rho * v, rho * w, by, bz, E];
}

/** Physical flux in x of a primitive state (Bx constant, its equation omitted). */
export function flux(W: Prim, gamma: number): number[] {
  const [rho, u, v, w, bx, by, bz, p] = W as [number, number, number, number, number, number, number, number];
  const B2 = bx * bx + by * by + bz * bz;
  const pt = p + B2 / 2;
  const E = p / (gamma - 1) + (rho * (u * u + v * v + w * w)) / 2 + B2 / 2;
  const vB = u * bx + v * by + w * bz;
  return [rho * u, rho * (u * u) + pt - bx * bx, rho * u * v - bx * by, rho * u * w - bx * bz,
    by * u - bx * v, bz * u - bx * w, (E + pt) * u - bx * vB];
}

/** Quasi-linear matrix of the primitive system q = (ρ, u, v, w, By, Bz, p), 7×7 row-major. */
export function primitiveJacobian(W: Prim, gamma: number): number[][] {
  const [rho, u, , , bx, by, bz, p] = W as [number, number, number, number, number, number, number, number];
  const A = Array.from({ length: 7 }, () => new Array<number>(7).fill(0));
  const set = (i: number, j: number, v: number) => { A[i - 1]![j - 1] = v; };
  set(1, 1, u); set(1, 2, rho);
  set(2, 2, u); set(2, 5, by / rho); set(2, 6, bz / rho); set(2, 7, 1 / rho);
  set(3, 3, u); set(3, 5, -bx / rho);
  set(4, 4, u); set(4, 6, -bx / rho);
  set(5, 2, by); set(5, 3, -bx); set(5, 5, u);
  set(6, 2, bz); set(6, 4, -bx); set(6, 6, u);
  set(7, 2, gamma * p); set(7, 7, u);
  return A;
}
