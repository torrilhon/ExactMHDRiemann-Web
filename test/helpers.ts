// Comparison helpers for the golden tests (plan §5.3): Julia's relerr, relative to
// max(|a|, |b|, 1).

import type { HState } from "../src/solver/types.ts";

export function relerr(a: readonly number[], b: readonly number[]): number {
  let m = 0;
  for (let i = 0; i < a.length; i++) {
    const x = a[i]!, y = b[i]!;
    const e = Math.abs(x - y) / Math.max(Math.abs(x), Math.abs(y), 1);
    if (Number.isNaN(e)) return NaN;
    m = Math.max(m, e);
  }
  return m;
}

export interface HStateJson { rho: number; u: number; p: number; bt: number; phi: number; vt: number[] }

export const hvec = (h: HState | HStateJson) => [h.rho, h.u, h.p, h.bt, h.phi, h.vt[0]!, h.vt[1]!];

/** Difference of two angles, modulo 2π (φ = π and φ = -π are the same direction). */
export function angleDiff(a: number, b: number): number {
  const d = Math.abs(a - b) % (2 * Math.PI);
  return Math.min(d, 2 * Math.PI - d);
}

/** relerr of two HStates (plus extra values), with φ compared modulo 2π. */
export function hrelerr(a: HState | HStateJson, b: HState | HStateJson, xa: number[] = [], xb: number[] = []): number {
  const va = hvec(a), vb = hvec(b);
  va[4] = 0; vb[4] = angleDiff(a.phi, b.phi);
  return relerr([...va, ...xa], [...vb, ...xb]);
}

export function fromJson(h: HStateJson): HState {
  return { rho: h.rho, u: h.u, p: h.p, bt: h.bt, phi: h.phi, vt: [h.vt[0]!, h.vt[1]!] };
}

/** Collects the worst error and the failing cases of a comparison loop. */
export class Worst {
  max = 0;
  bad: string[] = [];
  constructor(readonly tol: number) {}
  add(label: string, err: number): void {
    if (!(err <= this.max)) this.max = Number.isNaN(err) ? Infinity : err;
    if (!(err <= this.tol)) this.bad.push(`${label}: ${err.toExponential(3)}`);
  }
}
