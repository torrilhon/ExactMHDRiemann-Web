// Checked elementary functions (plan §4.6). Julia's sqrt, log and ^ throw DomainError
// for arguments outside their real domain, where JavaScript returns NaN. The solver's
// control flow depends on that exception (safe_residual, slow_limit, the quasi-Euler
// driver), so every sqrt/log/pow in src/solver goes through these wrappers; the
// lint:numerics check enforces it.


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

/**
 * Julia `evalpoly(x, (c0, c1, …))`: Horner, coefficients in ascending order.
 *
 * Julia writes it with `muladd`, which LLVM may or may not fuse into an fma depending on
 * the compilation context. Inside the solver's inlined code it mostly does not: plain
 * multiply-add reproduces Julia's fast-shock roots bit for bit in 600 of 621 golden
 * cases, fused fma in 476 (reference/BASELINE.md). The rest differ in the last bits.
 */
export function evalpoly(x: number, c: readonly number[]): number {
  let n = c.length - 1;
  let s = c[n] ?? 0;
  while (n > 0) {
    n--;
    s = s * x + (c[n] ?? 0);
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
