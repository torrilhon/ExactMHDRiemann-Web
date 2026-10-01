// Finite-difference Jacobian of the 5×5 residual (plan §4.3).
//
// - central differences, h_j = ε^{1/3} max(1, |x_j|);
// - branch-aware: if the central stencil of component j would cross one of its branch
//   points s (where build_side switches between the shock and the fan parametrization,
//   and the residual has a kink), a one-sided difference on the base point's side is
//   used instead (x > s: forward, x ≤ s: backward), with h = √ε max(1, |x_j|);
// - domain edges: a perturbed point whose residual is invalid (the BIG sentinel) is
//   replaced by the one-sided difference on the other side; if both sides are invalid
//   the column is zero and marked unreliable;
// - `before(kind)` lets the caller switch a step tape to record (base point, already
//   evaluated by the caller) or replay (all perturbed points).

import type { Mat, Vec } from "./linalg.ts";
import { zeros } from "./linalg.ts";
import { sqrtD } from "./math.ts";

const H_CENTRAL = Math.cbrt(Number.EPSILON);
const H_ONESIDED = sqrtD(Number.EPSILON);

export interface FdOptions {
  /** branch points per component (e.g. ψf: [0], ψs: [SLOW_EPS]) */
  switches?: readonly (readonly number[])[];
  /** true if a residual value is invalid (domain error sentinel) */
  invalid?: (r: readonly number[]) => boolean;
  /** called before the perturbed evaluations (to replay a step tape) */
  beforePerturbed?: () => void;
}

export interface FdResult {
  J: Mat;
  /** per column: false if no valid difference could be formed */
  ok: boolean[];
  nf: number;
}

export function fdJacobian(f: (x: readonly number[]) => Vec, x: readonly number[], fx: readonly number[],
  opts: FdOptions = {}): FdResult {
  const n = x.length, m = fx.length;
  const J = zeros(m, n);
  const ok = new Array<boolean>(n).fill(true);
  const invalid = opts.invalid ?? ((r: readonly number[]) => !r.every(Number.isFinite));
  let nf = 0;
  const at = (j: number, xj: number): Vec | null => {
    const y = x.slice();
    y[j] = xj;
    opts.beforePerturbed?.();
    nf++;
    const r = f(y);
    return invalid(r) ? null : r;
  };
  for (let j = 0; j < n; j++) {
    const xj = x[j]!;
    const scale = Math.max(1, Math.abs(xj));
    const hc = H_CENTRAL * scale, h1 = H_ONESIDED * scale;
    const sw = opts.switches?.[j] ?? [];
    const crosses = sw.some((s) => xj - hc <= s && s < xj + hc);
    const above = sw.some((s) => xj - hc <= s && s < xj + hc && xj > s);
    let col: Vec | null = null;
    if (!crosses) {
      const fp = at(j, xj + hc), fm = at(j, xj - hc);
      if (fp !== null && fm !== null) col = fp.map((v, i) => (v - fm[i]!) / (2 * hc));
      else if (fp !== null) col = fp.map((v, i) => (v - fx[i]!) / hc);
      else if (fm !== null) col = fm.map((v, i) => (fx[i]! - v) / hc);
    } else {
      // stay on the branch of the base point
      const first = above ? at(j, xj + h1) : at(j, xj - h1);
      if (first !== null) {
        col = above ? first.map((v, i) => (v - fx[i]!) / h1) : first.map((v, i) => (fx[i]! - v) / h1);
      }
    }
    if (col === null) { ok[j] = false; continue; }
    for (let i = 0; i < m; i++) J[i]![j] = col[i]!;
  }
  return { J, ok, nf };
}
