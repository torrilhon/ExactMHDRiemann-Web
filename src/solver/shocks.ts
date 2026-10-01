// Fast and slow shocks (paper Sec. 4.2 and Appendices A/B), in local dimensionless
// variables v̂ = v1/v0, p̂ = p1/p0, B̂t = Bt/√p0, A = B̂t of the upstream state,
// B = Bn/√p0, X = γM² (M = shock speed relative to the gas, over the sound speed).
// Mirrors shocks.jl.

import { findZeroBrent } from "./brent.ts";
import { speedsH } from "./eos.ts";
import { DomainError, EPS, evalpoly, sqrtD } from "./math.ts";
import type { Ctx, HState } from "./types.ts";
import { hstate } from "./types.ts";

const CBRT_EPS = Math.cbrt(EPS);

/**
 * Julia `ift_polish(f, fv, x0)` for Float64: one Newton polish step x0 - f(x0)/f'(x0).
 * (In Julia, with dual-valued parameters, the same step makes derivatives of the root
 * exact by the implicit function theorem; the finite-difference Jacobian of the port
 * does not need that.) `dfv` is the derivative of the value-level function at x0.
 */
export function iftPolish(f: (x: number) => number, x0: number, dfv: number): number {
  return x0 - f(x0) / dfv;
}

/**
 * Derivative of f at x0 by central differences, for `ift_polish` where Julia uses
 * ForwardDiff.derivative. Falls back to one side if the other leaves the domain.
 */
function derivative(f: (x: number) => number, x0: number): number {
  const h = CBRT_EPS * Math.max(Math.abs(x0), 1e-300);
  const ev = (x: number): number | null => {
    try {
      const v = f(x);
      return Number.isFinite(v) ? v : null;
    } catch (e) {
      if (e instanceof DomainError) return null;
      throw e;
    }
  };
  const fp = ev(x0 + h), fm = ev(x0 - h);
  if (fp !== null && fm !== null) return (fp - fm) / (2 * h);
  const f0 = f(x0);
  if (fp !== null) return (fp - f0) / h;
  if (fm !== null) return (f0 - fm) / h;
  return NaN;
}

/**
 * Coefficients (c0, c1, c2, c3) of the fast-shock cubic in y = B̂t1 (Julia `fast_cubic_bt`):
 *   c3 = -B²(k-1),  c2 = A(2D - X(k-1)),
 *   c1 = -D(A²(k+1) + 2Dk - 2X(k-1) + 2(k+1)),  c0 = 2 A D² k.
 */
export function fastCubicBt(A: number, B2: number, X: number, D: number, k: number): [number, number, number, number] {
  return [
    2 * A * (D * D) * k,
    -D * ((A * A) * (k + 1) + 2 * D * k - 2 * X * (k - 1) + 2 * (k + 1)),
    A * (2 * D - X * (k - 1)),
    -B2 * (k - 1),
  ];
}

/**
 * Largest real root of a real cubic in the interval (lo, hi], by bracketing on its
 * monotone pieces (Julia `largest_root_in`). The interval matters when Bn → 0: then the
 * leading coefficient c3 ∝ B² vanishes and a spurious root of size ~1/B² appears far
 * outside the physical range.
 */
export function largestRootIn(c: readonly [number, number, number, number], lo: number, hi: number): number {
  const [, c1, c2, c3] = c;
  const f = (x: number) => evalpoly(x, c);
  const disc = c2 * c2 - 3 * c3 * c1;                        // critical points of the cubic
  const crit = disc > 0 ? [(-c2 - sqrtD(disc)) / (3 * c3), (-c2 + sqrtD(disc)) / (3 * c3)] : [];
  const pts = [lo, ...crit.filter((x) => lo < x && x < hi), hi].sort((a, b) => a - b);
  for (let i = pts.length - 2; i >= 0; i--) {
    const a = pts[i]!, b = pts[i + 1]!;
    const fa = f(a), fb = f(b);
    if (fb === 0) return b;
    if (fa * fb < 0) return findZeroBrent(f, a, b, { xatol: 0.0, xrtol: 4 * EPS });
  }
  throw new DomainError(c, "fast shock: no root of the Hugoniot cubic in the physical range");
}

/**
 * Fast shock on the state U travelling in direction σ (-1 left, +1 right), with path
 * variable ψ = √(γ(M² - ĉA²)) - √(γ(ĉf² - ĉA²)) ≥ 0 (M: shock Mach number).
 * Returns (downstream, speed).
 */
export function fastShock(U: HState, psi: number, sigma: number, ctx: Ctx): [HState, number] {
  const { gamma, kappa, Bn } = ctx;
  const sp = sqrtD(U.p);
  const A = U.bt / sp, B = Bn / sp;
  const a0 = sqrtD((gamma * U.p) / U.rho);
  const [cf, cA, cs] = speedsH(U, ctx);
  // c_f² - c_A² without cancellation: (cf² - cA²)(cA² - cs²) = cA² bt²/ρ
  const df2 = cA > sqrtD((gamma * U.p) / U.rho)
    ? ((cA * cA) * (U.bt * U.bt)) / (U.rho * (cA * cA - cs * cs))
    : cf * cf - cA * cA;
  // Path variable ψ = √D - √D_f with D = X - B² = γ(M² - ĉA²) and D_f its value at M = ĉ_f
  const Df = (gamma * df2) / (a0 * a0);
  const sD = sqrtD(Df) + psi;
  const D = sD * sD;
  const X = B * B + D;
  const M = sqrtD(X / gamma);
  const c = fastCubicBt(A, B * B, X, D, kappa);
  // physical range of the downstream field: 0 < y and p̂ ≥ 0, i.e. y² ≤ A² + 2(1 + X)
  const yhi = sqrtD(A * A + 2 * (1 + X)) * (1 + 1e-12);
  const y0 = largestRootIn(c, 1e-300, yhi);
  const dc = [c[1], 2 * c[2], 3 * c[3]];                      // derivative of the cubic
  const y = iftPolish((x) => evalpoly(x, c), y0, evalpoly(y0, dc));
  const t = (D * (y - A)) / (X * y);                          // 1 - v̂
  const ph = 1 + X * t - (y * y - A * A) / 2;
  const bt1 = y * sp;
  const Cc = (-sigma * Bn) / (U.rho * a0 * M);
  const e0 = Math.cos(U.phi), e1 = Math.sin(U.phi);
  const k = Cc * (bt1 - U.bt);
  const Dn = hstate(U.rho / (1 - t), U.u + sigma * a0 * M * t, U.p * ph, bt1, U.phi,
    [U.vt[0] + k * e0, U.vt[1] + k * e1]);
  return [Dn, U.u + sigma * a0 * M];
}

/**
 * Slow shock lowering B̂t from A to A - Δ (paper eqs. 78-83; Julia `slow_volume`),
 * rewritten for t = 1 - v̂ = Δ z:  a z² - e1 z + e0 = 0. Returns (v̂, X = γM²).
 */
export function slowVolume(Delta: number, A: number, B2: number, k: number): [number, number] {
  const D = Delta;
  const A2 = A * A, D2 = D * D, D3 = D * D * D;
  const a = A2 * D * k - (3 * A * D2 * k) / 2 + (A * D2) / 2 + A * k + A + B2 * D * k + (D3 * k) / 2 - D3 / 2 - D * k - D;
  const e1 = (2 * A2 * k - 2 * A2 - 5 * A * D * k + 3 * A * D + 2 * B2 * k - 2 * B2 + 2 * D2 * k - 2 * D2 - 2 * k - 2) / 2;
  const e0 = (-(2 * A - D) * (k - 1)) / 2;
  const disc = e1 * e1 - 4 * a * e0;
  if (!(a > 0 && e0 < 0)) throw new DomainError(Delta, "slow shock outside the regular branch");
  const z = e1 > 0 ? (e1 + sqrtD(disc)) / (2 * a) : (-2 * e0) / (sqrtD(disc) - e1);
  return [1 - D * z, B2 / (1 + (A - D) * z)];
}

/** Slow shock that lowers B̂t by Δ > 0. Returns (downstream, speed, M, v̂). */
export function slowShockState(U: HState, Delta: number, sigma: number, ctx: Ctx): [HState, number, number, number] {
  const { gamma, kappa, Bn } = ctx;
  const sp = sqrtD(U.p);
  const A = U.bt / sp;
  const Bs = Bn / sp;
  const B2 = Bs * Bs;
  const a0 = sqrtD((gamma * U.p) / U.rho);
  const Bth = A - Delta;
  const [v, X] = slowVolume(Delta, A, B2, kappa);
  const M = sqrtD(X / gamma);
  // momentum balance (Rayleigh line) instead of the Hugoniot form of p̂
  const ph = 1 - X * (v - 1) - (Bth * Bth - A * A) / 2;
  const bt1 = Bth * sp;
  const Cc = (-sigma * Bn) / (U.rho * a0 * M);
  const e0 = Math.cos(U.phi), e1 = Math.sin(U.phi);
  const k = Cc * (bt1 - U.bt);
  const D = hstate(U.rho / v, U.u + sigma * a0 * M * (1 - v), U.p * ph, bt1, U.phi,
    [U.vt[0] + k * e0, U.vt[1] + k * e1]);
  return [D, U.u + sigma * a0 * M, M, v];
}

/**
 * Lax defect of a slow shock: (downstream relative speed)² - (downstream slow speed)².
 * Negative for a regular 3→4 shock, zero at the maximal Mach number.
 */
export function slowDefect(U: HState, Delta: number, ctx: Ctx): number {
  const [D, , M, v] = slowShockState(U, Delta, 1, ctx);
  const w1 = sqrtD((ctx.gamma * U.p) / U.rho) * M * v;
  const [, , cs1] = speedsH(D, ctx);
  return (w1 * w1 - cs1 * cs1) / ((ctx.gamma * U.p) / U.rho);
}

/**
 * Largest drop Δ of B̂t for which the slow shock on U stays a regular 3→4 shock with
 * B̂t > 0: either the point of maximal Mach number (paper eqs. 86-88) or Δ = A
 * (Julia `slow_limit`).
 */
export function slowLimit(U: HState, ctx: Ctx, N = 64): number {
  const A = U.bt / sqrtD(U.p);
  const fv = (Delta: number) => slowDefect(U, Delta, ctx);
  let Dprev = 1e-4 * A;
  const fprev = fv(Dprev);
  if (fprev >= 0) return Dprev;                               // degenerate: no regular slow shock
  for (let i = 1; i <= N; i++) {
    const Di = i === N ? A * (1 - 1e-9) : (A * i) / N;
    let fi: number;
    try {
      fi = fv(Di);
    } catch (e) {
      if (!(e instanceof DomainError)) throw e;
      fi = NaN;
    }
    if (Number.isNaN(fi)) {                                   // quadratic ceased to exist: end of the branch
      return Dprev;
    } else if (fi >= 0) {
      const D0 = findZeroBrent(fv, Dprev, Di, { xatol: 0.0, xrtol: 4 * EPS });
      return iftPolish(fv, D0, derivative(fv, D0));
    }
    Dprev = Di;
  }
  return U.bt / sqrtD(U.p);                                   // regular up to B̂t = 0
}
