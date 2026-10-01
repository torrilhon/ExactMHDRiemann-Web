// Small dense linear algebra: LU with partial pivoting (the 5×5 Newton/trust-region
// systems), matrix-vector products and ∞-norms (plan §4.5). Matrices are row-major
// arrays of rows.

export type Vec = number[];
export type Mat = number[][];

export function zeros(n: number, m: number): Mat {
  return Array.from({ length: n }, () => new Array<number>(m).fill(0));
}

export function matvec(A: Mat, x: readonly number[]): Vec {
  return A.map((row) => {
    let s = 0;
    for (let j = 0; j < row.length; j++) s += (row[j] ?? 0) * (x[j] ?? 0);
    return s;
  });
}

/** Aᵀ x */
export function matTvec(A: Mat, x: readonly number[]): Vec {
  const m = A[0]?.length ?? 0;
  const out = new Array<number>(m).fill(0);
  A.forEach((row, i) => {
    const xi = x[i] ?? 0;
    for (let j = 0; j < m; j++) out[j] = (out[j] ?? 0) + (row[j] ?? 0) * xi;
  });
  return out;
}

/** ‖x‖∞ (NaN if any component is NaN, like Julia's maximum(abs, x)). */
export function normInf(x: readonly number[]): number {
  let m = 0;
  for (const v of x) {
    const a = Math.abs(v);
    if (Number.isNaN(a)) return NaN;
    if (a > m) m = a;
  }
  return m;
}

export function norm2(x: readonly number[]): number {
  return Math.hypot(...x);
}

export function dot(x: readonly number[], y: readonly number[]): number {
  let s = 0;
  for (let i = 0; i < x.length; i++) s += (x[i] ?? 0) * (y[i] ?? 0);
  return s;
}

/** Julia `opnorm(A, Inf)`: maximum absolute row sum. */
export function opnormInf(A: Mat): number {
  let m = 0;
  for (const row of A) {
    let s = 0;
    for (const v of row) s += Math.abs(v);
    if (Number.isNaN(s)) return NaN;
    if (s > m) m = s;
  }
  return m;
}

export interface LU {
  lu: Mat;
  perm: number[];
  /** true if a zero (or non-finite) pivot was met */
  singular: boolean;
}

/** LU factorization with partial pivoting (LAPACK getrf order: pivot = largest |a| in the column). */
export function luFactor(A: Mat): LU {
  const n = A.length;
  const lu = A.map((r) => r.slice());
  const perm = Array.from({ length: n }, (_, i) => i);
  let singular = false;
  for (let k = 0; k < n; k++) {
    let p = k;
    let pmax = Math.abs(lu[k]![k]!);
    for (let i = k + 1; i < n; i++) {
      const v = Math.abs(lu[i]![k]!);
      if (v > pmax) { pmax = v; p = i; }
    }
    if (p !== k) {
      [lu[p], lu[k]] = [lu[k]!, lu[p]!];
      [perm[p], perm[k]] = [perm[k]!, perm[p]!];
    }
    const rk = lu[k]!;
    const piv = rk[k]!;
    if (piv === 0 || !Number.isFinite(piv)) { singular = true; continue; }
    for (let i = k + 1; i < n; i++) {
      const ri = lu[i]!;
      const l = ri[k]! / piv;
      ri[k] = l;
      if (l === 0) continue;
      for (let j = k + 1; j < n; j++) ri[j] = ri[j]! - l * rk[j]!;
    }
  }
  return { lu, perm, singular };
}

/** Solve A x = b from a factorization (no check for singularity: see `singular`). */
export function luSolve(f: LU, b: readonly number[]): Vec {
  const { lu, perm } = f;
  const n = lu.length;
  const x = perm.map((p) => b[p] ?? 0);
  for (let i = 1; i < n; i++) {
    const r = lu[i]!;
    let s = x[i]!;
    for (let j = 0; j < i; j++) s -= r[j]! * x[j]!;
    x[i] = s;
  }
  for (let i = n - 1; i >= 0; i--) {
    const r = lu[i]!;
    let s = x[i]!;
    for (let j = i + 1; j < n; j++) s -= r[j]! * x[j]!;
    x[i] = s / r[i]!;
  }
  return x;
}

/** Solve A x = b; returns null if A is singular to working precision. */
export function solve(A: Mat, b: readonly number[]): Vec | null {
  const f = luFactor(A);
  if (f.singular) return null;
  const x = luSolve(f, b);
  return x.every(Number.isFinite) ? x : null;
}
