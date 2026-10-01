// Brent's method with the semantics of Roots.jl 3.0.9 `find_zero(f, (a, b), Roots.Brent(); kw...)`
// (src/Bracketing/brent.jl, bracketing.jl, convergence.jl, find_zero.jl), plan §4.4.
//
// What the solver relies on:
// - default tolerances xatol = eps^3, xrtol = eps, atol = rtol = 0, at most 60 steps;
// - with atol = rtol = 0 the stopping test is on the bracket width only ("fexact");
// - bisection steps use the midpoint of the bit patterns (`__middle`), so a bracket
//   like (1e-300, 1) is halved in exponent, not linearly;
// - a bracket without sign change throws ArgumentError, a failed search throws
//   ConvergenceFailed (neither is a DomainError: callers rethrow them as in Julia).

import { ArgumentError, EPS, sqrtD } from "./math.ts";

/** Roots.jl `ConvergenceFailed`. */
export class ConvergenceFailed extends Error {
  constructor(message = "Algorithm failed to converge") {
    super(message);
    this.name = "ConvergenceFailed";
  }
}

export interface BrentOptions {
  xatol?: number;
  xrtol?: number;
  /** maximum number of update steps (Roots.jl `maxevals`, default 60) */
  maxevals?: number;
}

const f64 = new Float64Array(1);
const u64 = new BigUint64Array(f64.buffer);

function bitsOf(x: number): bigint {
  f64[0] = x;
  return u64[0]!;
}
function fromBits(b: bigint): number {
  u64[0] = b;
  return f64[0]!;
}

/** Julia `nextfloat` / `prevfloat` for finite doubles. */
export function nextfloat(x: number): number {
  if (Number.isNaN(x) || x === Infinity) return x;
  if (x === 0) return Number.MIN_VALUE;
  const b = bitsOf(x);
  return fromBits(x > 0 ? b + 1n : b - 1n);
}
export function prevfloat(x: number): number {
  return -nextfloat(-x);
}

/** Roots.jl `__middle` for Float64: midpoint of the bit patterns of |x|, |y|. */
function bitMiddle(x: number, y: number): number {
  const mid = (bitsOf(Math.abs(x)) + bitsOf(Math.abs(y))) >> 1n;
  return Math.sign(x + y) * fromBits(mid);
}

/** Roots.jl `_middle`. */
function middle(x: number, y: number): number {
  const a = x === -Infinity || x === Infinity ? nextfloat(x) : x;
  const b = y === -Infinity || y === Infinity ? prevfloat(y) : y;
  if (Math.sign(a) * Math.sign(b) < 0) return 0;
  return bitMiddle(a, b);
}

interface State {
  xn1: number; xn0: number; c: number; d: number;
  fxn1: number; fxn0: number; fc: number; mflag: boolean;
}

function inverseQuadraticStep(a: number, b: number, c: number, fa: number, fb: number, fc: number): number {
  let s = 0;
  s += (a * fb * fc) / (fa - fb) / (fa - fc);
  s += (b * fa * fc) / (fb - fa) / (fb - fc);
  s += (c * fa * fb) / (fc - fa) / (fc - fb);
  return s;
}

const secantStep = (a: number, b: number, fa: number, fb: number) => a - (fa * (b - a)) / (fb - fa);

/**
 * Root of f in the bracket (a, b), as Roots.jl `find_zero(f, (a, b), Brent(); xatol, xrtol)`.
 */
export function findZeroBrent(f: (x: number) => number, a: number, b: number, opts: BrentOptions = {}): number {
  const xatol = opts.xatol ?? EPS * EPS * EPS;
  const xrtol = opts.xrtol ?? EPS;
  const maxiters = opts.maxevals ?? 60;
  // atol = rtol = 0 here; Roots.jl picks the stopping rule from which tolerances vanish
  const exact: "xexact" | "fexact" = xatol === 0 && xrtol === 0 ? "xexact" : "fexact";

  // adjust_bracket, then Brent's init_state
  if (a === b) throw new ArgumentError("Need extrema to return two distinct values");
  let x0 = Math.min(a, b), x1 = Math.max(a, b);
  if (x0 === -Infinity) x0 = nextfloat(x0);
  if (x1 === Infinity) x1 = prevfloat(x1);
  let u = x0, v = x1;
  let fu = f(u), fv = f(v);
  if (Math.abs(fu) > Math.abs(fv)) [u, v, fu, fv] = [v, u, fv, fu];
  if (!(fu === 0 || fv === 0) && !(Math.sign(fu) * Math.sign(fv) < 0)) {
    throw new ArgumentError(
      "The interval [a,b] is not a bracketing interval. You need f(a) and f(b) to have different signs (f(a) * f(b) < 0).",
    );
  }
  let st: State = { xn1: u, xn0: v, c: v, d: v, fxn1: fu, fxn0: fv, fc: fv, mflag: true };

  const assess = (s: State): [string, boolean] => {
    // AbstractBracketingMethod with exact = :xexact / :fexact
    if (s.fxn0 === 0 || s.fxn1 === 0) return ["exact_zero", true];
    if (Number.isNaN(s.fxn1) || Number.isNaN(s.fxn0)) return ["nan", true];
    if (exact === "xexact") {
      let lo = s.xn0, hi = s.xn1;
      if (hi < lo) [lo, hi] = [hi, lo];
      if (nextfloat(lo) === hi) return ["x_converged", true];
    } else {
      const uu = Math.abs(s.fxn0) < Math.abs(s.fxn1) ? s.xn0 : s.xn1;     // choose_smallest
      const dx = Math.max(xatol, 2 * (Math.abs(uu) * xrtol));
      if (Math.abs(s.xn1 - s.xn0) <= dx) return ["x_converged", true];
    }
    return ["not_converged", false];
  };

  const update = (s: State): [State, boolean] => {
    let mflag = s.mflag;
    let a = s.xn0, b = s.xn1, c = s.c, d = s.d;
    let fa = s.fxn0, fb = s.fxn1, fc = s.fc;

    let sx = inverseQuadraticStep(a, b, c, fa, fb, fc);
    if (!Number.isFinite(sx)) sx = secantStep(a, b, fa, fb);

    let lo = (3 * a + b) / 4, hi = b;
    if (lo > hi) [lo, hi] = [hi, lo];

    const tol = Math.max(xatol, Math.max(Math.abs(b), Math.abs(c), Math.abs(d)) * xrtol);
    if (
      !(lo < sx && sx < hi) ||
      (mflag && Math.abs(sx - b) >= Math.abs(b - c) / 2) ||
      (!mflag && Math.abs(sx - b) >= Math.abs(c - d) / 2) ||
      (mflag && Math.abs(b - c) <= tol) ||
      (!mflag && Math.abs(c - d) <= tol)
    ) {
      sx = middle(a, b);
      mflag = true;
    } else {
      mflag = false;
    }

    const fs = f(sx);
    if (fs === 0) return [{ ...s, xn1: sx, fxn1: fs }, true];
    if (!Number.isFinite(fs)) return [s, true];

    d = c;
    c = b; fc = fb;
    if (Math.sign(fa) * Math.sign(fs) < 0) { b = sx; fb = fs; } else { a = sx; fa = fs; }
    if (Math.abs(fa) < Math.abs(fb)) [a, b, fa, fb] = [b, a, fb, fa];
    return [{ xn1: b, fxn1: fb, xn0: a, fxn0: fa, c, d, fc, mflag }, false];
  };

  // solve!
  let ctr = 1;
  for (;;) {
    const [, stopped] = assess(st);
    if (stopped) break;
    if (ctr > maxiters) break;
    const [ns, stop2] = update(st);
    st = ns;
    ctr++;
    if (stop2) break;
  }
  const [val] = assess(st);

  // decide_convergence (AbstractBracketingMethod)
  const bb = st.xn1, aa = st.xn0, fb = st.fxn1, fa = st.fxn0;
  const xs = [bb, aa] as const;
  const fxs = [fb, fa] as const;
  // findmin(abs, (fb, fa)): first minimum; a NaN is returned if present
  const i = Number.isNaN(fb) ? 0 : Number.isNaN(fa) ? 1 : Math.abs(fb) <= Math.abs(fa) ? 0 : 1;
  const alpha = xs[i];
  const m = Math.abs(fxs[i]);
  let result: number;
  if (val === "exact_zero") {
    result = alpha;
  } else if (val === "x_converged") {
    const [uu, vv, fuu, fvv] = aa < bb ? [aa, bb, fa, fb] : [bb, aa, fb, fa];
    const up = nextfloat(uu);
    if (vv === nextfloat(up)) {
      const fup = f(up);
      const cand = [uu, up, vv], fc = [Math.abs(fuu), Math.abs(fup), Math.abs(fvv)];
      let k = 0;
      for (let j = 1; j < 3; j++) if (fc[j]! < fc[k]!) k = j;
      result = cand[k]!;
    } else {
      result = alpha;
    }
  } else if (val === "not_converged") {
    // strict = false: relaxed convergence on Δx and f(x)
    const mx = Math.max(Math.abs(xs[0]), Math.abs(xs[1]));
    const dx = Math.abs(xs[0] - xs[1]);
    const eps = 256 * Math.max(xatol, mx * sqrtD(xrtol));
    if (dx <= eps) {
      result = alpha;
    } else {
      const atol = EPS, rtol = Math.max(EPS, xrtol);
      const delta = 16 * Math.min(16 * atol, mx * rtol);
      result = m <= delta ? alpha : NaN;
    }
  } else if (Number.isNaN(fa)) result = aa;
  else if (Number.isNaN(fb)) result = bb;
  else if (fa === 0) result = aa;
  else if (fb === 0) result = bb;
  else result = alpha;

  if (Number.isNaN(result)) throw new ConvergenceFailed();
  return result;
}
