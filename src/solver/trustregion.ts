// Trust-region solver for small dense square systems (plan §4.2): a port of what
// `NonlinearSolve.TrustRegion()` does at the pinned version (NonlinearSolve 4.32.0,
// NonlinearSolveFirstOrder 2.10.0, NonlinearSolveBase 2.54.1):
//
// - descent `MoreTrustRegionDescent` (scaling None, D = I) in its normal-equation form,
//   which NonlinearSolve uses for static arrays (the solver's SVector{5} problems): the
//   subproblem min ‖J p + F‖, ‖p‖ ≤ Δ is solved nearly exactly by Moré's safeguarded
//   iteration on the damping λ: Gauss-Newton step from JᵀJ p = -JᵀF if it lies inside
//   the region, otherwise at most 10 λ updates of (JᵀJ + λI) p = -JᵀF until
//   |‖p‖ - Δ| ≤ θΔ with θ = 1e-4; λ ≥ eps·max diag(JᵀJ), warm-started from the previous
//   solve; predicted reduction ½‖Jp‖² + λ‖p‖² (MINPACK form). The linear systems are
//   solved by LU with partial pivoting, which — as StaticArrays' `lu` — throws
//   SingularException on an exactly zero pivot (e.g. a zero column of J); `nl_solve`
//   treats that as a failed start;
// - radius update `RadiusUpdateSchemes.More` (MINPACK `lmder`): accept if ρ > 1e-3;
//   ρ < 1/4 → Δ = ¼ min(Δ, 10‖p‖); ρ ≥ 3/4 or λ = 0 → Δ = 2‖p‖; initial Δ = ‖u₀‖₂ if
//   > 1e-4, else 1; a non-finite trial residual is rejected (ρ = -1);
// - Jacobian recomputed after accepted steps only;
// - termination `AbsNormSafeBestTerminationMode(‖·‖∞; max_stalled_steps = 32)`:
//   success at ‖F‖∞ ≤ abstol, stop on a non-finite residual, after `maxiters` steps
//   (accepted or not) or more than 32 consecutive shrinks; the best iterate (smallest
//   ‖F‖∞) is returned.
//
// Bit-for-bit agreement with NonlinearSolve is not required (plan §4.2); the Jacobian
// comes from finite differences, not ForwardDiff.

import type { Mat, Vec } from "./linalg.ts";
import { dot, luFactor, luSolve, matvec, matTvec, norm2, normInf } from "./linalg.ts";
import { sqrtD } from "./math.ts";

export type NlRetcode = "Success" | "MaxIters" | "Unstable" | "Stalled" | "ShrinkThresholdExceeded"
  | "InternalLinearSolveFailed";

export interface TrustRegionOptions {
  abstol: number;
  maxiters: number;
  maxShrinkTimes?: number;
}

export interface TrustRegionResult {
  u: Vec;
  /** ‖F(u)‖∞ of the returned (best) iterate */
  resnorm: number;
  retcode: NlRetcode;
  nsteps: number;
  nf: number;
  njac: number;
}

export type Residual = (u: readonly number[]) => Vec;
export type JacobianFn = (u: readonly number[], fu: readonly number[]) => Mat;

const EPS = Number.EPSILON;

/** Julia `LinearAlgebra.SingularException`: exactly zero pivot in an LU factorization. */
export class SingularException extends Error {
  readonly info: number;
  constructor(info: number) {
    super(`SingularException(${info})`);
    this.name = "SingularException";
    this.info = info;
  }
}

/** Solve A x = b by LU with partial pivoting; throws SingularException on a zero pivot. */
function luSolveChecked(A: Mat, b: readonly number[]): { x: Vec; f: ReturnType<typeof luFactor> } {
  const f = luFactor(A);
  for (let k = 0; k < A.length; k++) {
    if (f.lu[k]![k] === 0) throw new SingularException(k + 1);
  }
  return { x: luSolve(f, b), f };
}
const THETA = 1e-4;
const LAMBDA_ITERS = 10;

const finite = (v: readonly number[]) => v.every(Number.isFinite);

// ---------------------------------------------------------------- descent (lmpar)

interface DescentState {
  lambda: number;
  boundLambda: number;
  boundNorm: number;
  boundValid: boolean;
  gn: Vec | null;
  gnNorm: number;
  gnValid: boolean;
  JtJ?: Mat;
}

interface DescentResult {
  ok: boolean;
  p: Vec;
  lambda: number;
  predicted: number;
  stepNorm: number;
}

function extras(J: Mat, p: Vec, lambda: number): Omit<DescentResult, "ok"> {
  const Jp = matvec(J, p);
  return { p, lambda, predicted: dot(Jp, Jp) / 2 + lambda * dot(p, p), stepNorm: norm2(p) };
}

function descent(st: DescentState, J: Mat, fu: readonly number[], delta: number, newJacobian: boolean): DescentResult {
  const n = J[0]!.length;
  if (newJacobian) {
    st.gnValid = false; st.boundValid = false;
    st.JtJ = Array.from({ length: n }, (_, a) => Array.from({ length: n }, (_, b) => {
      let x = 0;
      for (let k = 0; k < J.length; k++) x += J[k]![a]! * J[k]![b]!;
      return x;
    }));
  }
  const JtJ = st.JtJ!;
  const Jtf = matTvec(J, fu);
  if (norm2(Jtf) === 0) {
    // stationary point of the model: p = 0 for every Δ
    return { ok: true, ...extras(J, new Array<number>(n).fill(0), 0) };
  }
  if (!st.gnValid) {
    // Gauss-Newton: JᵀJ p = -JᵀF (throws SingularException on a zero pivot, as Julia)
    const { x } = luSolveChecked(JtJ, Jtf);
    const p = x.map((v) => -v);
    st.gn = finite(p) ? p : null;
    st.gnNorm = st.gn === null ? Infinity : norm2(st.gn);
    st.gnValid = true;
  }
  if (st.gn !== null && st.gnNorm <= delta) return { ok: true, ...extras(J, st.gn, 0) };

  // Moré's safeguarded Newton iteration on λ (normal form, `_more_generic_λloop!`)
  const ubound = norm2(Jtf) / delta;
  let lambda = st.lambda === 0 ? 1e-3 * ubound : Math.min(st.lambda, ubound);
  let l = 0, u = ubound;
  if (st.boundValid) {
    if (st.boundNorm >= delta) {
      l = Math.max(l, st.boundLambda);
      lambda = Math.min(st.boundLambda * (st.boundNorm / delta), ubound);
    } else {
      u = Math.min(u, st.boundLambda);
      lambda = Math.min(lambda, u);
    }
  }
  // below ~eps·maxdiag(JᵀJ) the normal-equations factorization cannot succeed
  lambda = Math.max(lambda, EPS * Math.max(...JtJ.map((r, k) => r[k]!)));
  u = Math.max(u, lambda);
  st.lambda = lambda;
  let got: Vec | null = null;
  let posBound = l > 0;
  let phiPrev = 0, lambdaPrev = lambda;
  for (let i = 1; i <= LAMBDA_ITERS; i++) {
    const A = JtJ.map((r, a) => r.map((v, b) => (a === b ? v + lambda : v)));
    const { x, f } = luSolveChecked(A, Jtf);
    if (!finite(x)) {
      l = Math.max(l, lambda); lambda *= 10; u = Math.max(u, lambda);
      continue;
    }
    const p = x.map((v) => -v);
    got = p;
    st.lambda = lambda;
    const pnorm = norm2(p);
    const phi = pnorm - delta;
    st.boundLambda = lambda; st.boundNorm = pnorm; st.boundValid = true;
    if (Math.abs(phi) <= THETA * delta || i === LAMBDA_ITERS) break;
    // MINPACK's `parl == 0` exit (rank-deficient hard case)
    if (!posBound && phiPrev < 0 && phi <= phiPrev + THETA * delta && lambda <= lambdaPrev) break;
    if (phi > 0) posBound = true;
    phiPrev = phi; lambdaPrev = lambda;
    // q = (JᵀJ + λI)⁻¹ p
    const q = luSolve(f, p);
    if (!finite(q)) {
      l = Math.max(l, lambda); lambda *= 10; u = Math.max(u, lambda);
      continue;
    }
    // _more_update_λ
    if (phi < 0) u = lambda; else l = lambda;
    lambda += (phi / delta) * dot(p, p) / dot(p, q);
    if (!(l <= lambda && lambda <= u)) lambda = Math.max(l + 0.01 * (u - l), sqrtD(l * u));
  }
  if (got === null) return { ok: false, p: new Array<number>(n).fill(0), lambda: NaN, predicted: NaN, stepNorm: NaN };
  return { ok: true, ...extras(J, got, st.lambda) };
}

// ---------------------------------------------------------------- driver

/** Solve F(u) = 0 from u0. `jac(u, F(u))` returns the Jacobian at u. */
export function trustRegion(F: Residual, jac: JacobianFn, u0: readonly number[], opts: TrustRegionOptions): TrustRegionResult {
  const maxShrink = opts.maxShrinkTimes ?? 32;
  let nf = 0, njac = 0;
  const f = (u: readonly number[]) => { nf++; return F(u); };
  let u = u0.slice();
  let fu = f(u);
  let best = u.slice(), bestObj = normInf(fu);
  const done = (retcode: NlRetcode, nsteps: number): TrustRegionResult => {
    // safe-best: roll back to the best iterate
    return { u: best, resnorm: bestObj, retcode, nsteps, nf, njac };
  };
  if (!Number.isFinite(bestObj)) return done("Unstable", 0);
  if (bestObj <= opts.abstol) return done("Success", 0);

  const u0norm = norm2(u);
  let delta = u0norm > 1e-4 ? u0norm : 1;
  const st: DescentState = {
    lambda: 0, boundLambda: 0, boundNorm: 0, boundValid: false, gn: null, gnNorm: Infinity, gnValid: false,
  };
  let J: Mat = [];
  let newJacobian = true;
  let shrink = 0;
  const stepNorms: number[] = [];
  for (let nsteps = 1; nsteps <= opts.maxiters; nsteps++) {
    if (newJacobian) { J = jac(u, fu); njac++; }
    const d = descent(st, J, fu, delta, newJacobian);
    if (!d.ok) {
      if (newJacobian) return done("InternalLinearSolveFailed", nsteps);
      newJacobian = true;
      nsteps--;                                   // the retry is part of this step
      continue;
    }
    const uNew = u.map((v, i) => v + d.p[i]!);
    const fNew = f(uNew);
    const num = (dot(fNew, fNew) - dot(fu, fu)) / 2;
    const denom = -d.predicted;
    const finiteTrial = Number.isFinite(num) && finite(fNew);
    const rho = denom < 0 && finiteTrial ? num / denom : -1;
    const accepted = rho > 1e-3;
    if (rho < 0.25) {
      delta = 0.25 * Math.min(delta, 10 * d.stepNorm);
      shrink++;
    } else {
      shrink = 0;
      if (rho >= 0.75 || d.lambda === 0) delta = 2 * d.stepNorm;
    }
    const uprev = u;
    if (accepted) { u = uNew; fu = fNew; newJacobian = true; } else { newJacobian = false; }
    // the shrink limit stops the solve, but the termination check below still runs and,
    // if it triggers, its retcode wins (NonlinearSolveFirstOrder step! order)
    const shrinkStop = shrink > maxShrink;
    // termination (AbsNormSafeBest, ∞-norm)
    const obj = normInf(fu);
    if (!Number.isFinite(obj)) return done("Unstable", nsteps);
    if (obj < bestObj) { bestObj = obj; best = u.slice(); }
    if (obj <= opts.abstol) return done("Success", nsteps);
    stepNorms.push(norm2(u.map((v, i) => v - uprev[i]!)));
    if (stepNorms.length > 32) stepNorms.shift();
    if (nsteps > 32 && Math.max(...stepNorms) <= opts.abstol) return done("Stalled", nsteps);
    if (shrinkStop) return done("ShrinkThresholdExceeded", nsteps);
  }
  return done("MaxIters", opts.maxiters);
}
