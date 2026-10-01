// Checked elementary functions (plan §4.6). Julia's sqrt, log and ^ throw DomainError
// for arguments outside their real domain, where JavaScript returns NaN. The solver's
// control flow depends on that exception (safe_residual, slow_limit, the quasi-Euler
// driver), so every sqrt/log/pow in src/solver goes through these wrappers; the
// lint:numerics check enforces it.

import { exactRational, ratToDouble } from "./exact.ts";

/** Mirrors Julia's `DomainError`: an argument outside the domain of a function. */
export class DomainError extends Error {
  readonly value: unknown;
  constructor(value: unknown, message = "") {
    super(message === "" ? `DomainError with ${String(value)}` : `DomainError with ${String(value)}: ${message}`);
    this.name = "DomainError";
    this.value = value;
  }
}

/** Julia `ArgumentError`. */
export class ArgumentError extends Error {
  constructor(message: string) {
    super(message);
    this.name = "ArgumentError";
  }
}

/** Julia `eps()` = eps(Float64). */
export const EPS = Number.EPSILON;

/** Julia `sqrt`: throws DomainError for x < 0 (-0.0 and NaN pass through). */
export function sqrtD(x: number): number {
  if (x < 0) throw new DomainError(x, "sqrt was called with a negative real argument");
  return Math.sqrt(x);
}

/** Julia `log`: throws DomainError for x < 0 (log(0) = -Inf, NaN passes through). */
export function logD(x: number): number {
  if (x < 0) throw new DomainError(x, "log was called with a negative real argument");
  return Math.log(x);
}

/**
 * Julia `x^y` for Float64: throws DomainError for x < 0 and non-integer y. Literal
 * small integer powers in the Julia source (`x^2`, `x^3`) are written as products
 * instead, matching Julia's `literal_pow`.
 */
export function powD(x: number, y: number): number {
  if (x < 0 && Number.isFinite(y) && !Number.isInteger(y)) {
    throw new DomainError(x, "exponentiation of a negative number to a non-integer power");
  }
  return Math.pow(x, y);
}

// ---------------------------------------------------------------- fused multiply-add

const SPLITTER = 134217729;            // 2^27 + 1 (Veltkamp)
const f64 = new Float64Array(1);
const u64 = new BigUint64Array(f64.buffer);

function twoSum(a: number, b: number): [number, number] {
  const s = a + b;
  const bb = s - a;
  return [s, (a - (s - bb)) + (b - bb)];
}

function twoProd(a: number, b: number): [number, number] {
  const p = a * b;
  let t = SPLITTER * a;
  const ah = t - (t - a), al = a - ah;
  t = SPLITTER * b;
  const bh = t - (t - b), bl = b - bh;
  return [p, ((ah * bh - p) + ah * bl + al * bh) + al * bl];
}

/** Round-to-odd of a + b (the exact sum is s + e). */
function addOdd(a: number, b: number): number {
  const [s, e] = twoSum(a, b);
  if (e === 0 || s === 0) return s;
  f64[0] = s;
  if ((u64[0]! & 1n) === 1n) return s;            // already odd
  // move one ulp towards the exact value: the neighbour with odd last bit
  u64[0] = (s > 0) === (e > 0) ? u64[0]! + 1n : u64[0]! - 1n;
  return f64[0]!;
}

/**
 * Correctly rounded a·b + c, as Julia's `fma` (and `muladd`, which Julia compiles to a
 * hardware FMA on x86-64 with FMA3 and on ARM64). JavaScript has no fma, so it is
 * emulated: Boldo & Melquiond, "Emulation of FMA and correctly rounded sums: proved
 * algorithms using rounding to odd", IEEE Trans. Comput. 57 (2008), with an exact
 * BigInt evaluation where the error-free transformations could under- or overflow.
 */
export function fma(a: number, b: number, c: number): number {
  const p = a * b;
  if (!Number.isFinite(a) || !Number.isFinite(b)) return p + c;   // Inf/NaN as in IEEE fma
  if (!Number.isFinite(c)) return c;                               // exact a·b is finite
  if (a === 0 || b === 0) return p + c;           // exact product: one rounding, signed zeros as fma
  const ap = Math.abs(p), aa = Math.abs(a), ab = Math.abs(b), ac = Math.abs(c);
  const SAFE_LO = 2 ** -900, SAFE_HI = 2 ** 900;
  if (ap < SAFE_LO || ap > SAFE_HI || aa > SAFE_HI || ab > SAFE_HI || aa < SAFE_LO || ab < SAFE_LO ||
      (ac !== 0 && (ac < SAFE_LO || ac > SAFE_HI))) {
    return fmaExact(a, b, c);
  }
  const [uh, ul] = twoProd(a, b);
  const [th, tl] = twoSum(c, uh);
  const v = addOdd(tl, ul);
  const r = th + v;
  // exact cancellation to zero: IEEE gives +0 in round-to-nearest
  return r === 0 ? 0 : r;
}

function fmaExact(a: number, b: number, c: number): number {
  const [an, ad] = exactRational(a);
  const [bn, bd] = exactRational(b);
  const [cn, cd] = exactRational(c);
  // a·b + c = (an·bn·cd + cn·ad·bd) / (ad·bd·cd)
  const num = an * bn * cd + cn * ad * bd;
  if (num === 0n) return 0;
  return ratToDouble(num, ad * bd * cd);
}

/**
 * Julia `evalpoly(x, (c0, c1, …))`: Horner with `muladd`, i.e. ex = fma(x, ex, c[i]),
 * coefficients in ascending order (base/math.jl `_evalpoly`).
 */
export function evalpoly(x: number, c: readonly number[]): number {
  let n = c.length - 1;
  let s = c[n] ?? 0;
  while (n > 0) {
    n--;
    s = fma(x, s, c[n] ?? 0);
  }
  return s;
}

/** Julia `atan(y, x)`. */
export const atan2 = Math.atan2;

/** Julia `hypot(x, y)`. */
export function hypot(x: number, y: number): number {
  return Math.hypot(x, y);
}

/** Julia `sign`: -1, 0 or 1 (NaN for NaN, keeps -0.0). */
export const sign = Math.sign;
