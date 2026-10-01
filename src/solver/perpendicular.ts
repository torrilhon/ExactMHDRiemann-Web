// Quasi-Euler solver for a vanishing normal field (Bn = 0, used for Bn/√p ≤ bn_euler).
// Mirrors perpendicular.jl.
//
// With Bn = 0 the slow and Alfvén waves merge into the contact, which becomes a
// tangential discontinuity: u and the total pressure P = p + Bt²/2 are continuous,
// while ρ, p, |Bt|, the direction of Bt and vt may jump. Across the fast waves Bt/ρ,
// the direction of Bt and vt stay constant, and the fast speed is cf² = (γp + Bt²)/ρ.
// One unknown, P*, is found from u*_L(P*) = u*_R(P*) by bracketing.

import { findZeroBrent } from "./brent.ts";
import { toHState } from "./canonical.ts";
import { checkWaves } from "./check.ts";
import { DomainError, EPS, logD, powD, sqrtD } from "./math.ts";
import type { Rhs } from "./rk/engine.ts";
import { integrate } from "./rk/engine.ts";
import { VERN9 } from "./rk/vern9.ts";
import type { Ctx, FanData, Frame, HState, RiemannProblem, RiemannSolution, SolverOptions, Wave } from "./types.ts";
import { hstate, makeCtx } from "./types.ts";

export const ptot = (rho: number, p: number, bt: number) => p + (bt * bt) / 2;

// Fast fan on U down to total pressure P < P0, parametrized by s = log(ρ0/ρ) = s* τ.
function perpFanRhs(rho0: number, p0: number, bt0: number, send: number, gamma: number, sigma: number): Rhs {
  return (_y, tau) => {
    const s = send * tau;
    const rho = rho0 * Math.exp(-s), p = p0 * Math.exp(-gamma * s), bt = bt0 * Math.exp(-s);
    const cf = sqrtD((gamma * p + bt * bt) / rho);
    return [send * (-sigma * cf)];
  };
}

export function perpFan(U: HState, P: number, sigma: number, ctx: Ctx, record = false): [HState, FanData | null] {
  const gamma = ctx.gamma;
  const g = (rho: number) => ptot(rho, U.p * powD(rho / U.rho, gamma), (U.bt * rho) / U.rho) - P;
  const rhos = findZeroBrent(g, 1e-300, U.rho, { xatol: 0.0, xrtol: 4 * EPS });
  const send = logD(U.rho / rhos);
  const sol = integrate(VERN9, perpFanRhs(U.rho, U.p, U.bt, send, gamma, sigma), [U.u], 0.0, 1.0,
    { abstol: ctx.opts.ode_tol, reltol: ctx.opts.ode_tol, dense: record });
  if (!sol.successful) throw new DomainError(P, "quasi-Euler fan ODE failed");
  const D = hstate(rhos, sol.end[0]!, U.p * powD(rhos / U.rho, gamma), (U.bt * rhos) / U.rho, U.phi, U.vt);
  return [D, record ? { family: "fast0", sigma, up: U, par: send, sol } : null];
}

// Fast shock on U up to total pressure P > P0: compression r = ρ1/ρ0 ∈ (1, κ) from the
// energy jump condition with Bt ∝ ρ, then the mass flux from the momentum condition.
export function perpShock(U: HState, P: number, sigma: number, ctx: Ctx): [HState, number] {
  const { gamma, kappa } = ctx;
  const A = U.bt / sqrtD(U.p);
  const ph = (r: number) => {
    const v = 1 / r;
    const d = A * r - A;
    return ((1 - v) * (d * d) / 4 + 1 / (gamma - 1) - (v - 1) / 2) / (v / (gamma - 1) + (v - 1) / 2);
  };
  const f = (r: number) => { const b = U.bt * r; return U.p * ph(r) + (b * b) / 2 - P; };
  const rmax = kappa * (1 - 1e-14);
  // beyond P/P0 ~ 1e14 the compression is indistinguishable from κ in Float64
  if (!(f(rmax) > 0)) throw new DomainError(P, "quasi-Euler shock too strong");
  const r = findZeroBrent(f, 1.0, rmax, { xatol: 0.0, xrtol: 4 * EPS });
  const rho1 = U.rho * r;
  const P0 = ptot(U.rho, U.p, U.bt);
  const m = sqrtD((P - P0) / (1 / U.rho - 1 / rho1));            // mass flux
  const D = hstate(rho1, U.u + (sigma * (P - P0)) / m, U.p * ph(r), U.bt * r, U.phi, U.vt);
  return [D, U.u + (sigma * m) / U.rho];
}

export function perpSide(U: HState, P: number, sigma: number, ctx: Ctx, record = false): [HState, Wave | null] {
  const P0 = ptot(U.rho, U.p, U.bt);
  if (P > P0 * (1 + 1e-14)) {
    const [D, s] = perpShock(U, P, sigma, ctx);
    const w: Wave | null = record
      ? (sigma < 0 ? { kind: "fast_shock", side: sigma, s_left: s, s_right: s, left: U, right: D, fan: null }
        : { kind: "fast_shock", side: sigma, s_left: s, s_right: s, left: D, right: U, fan: null })
      : null;
    return [D, w];
  } else if (P < P0 * (1 - 1e-14)) {
    const [D, fan] = perpFan(U, P, sigma, ctx, record);
    if (record) {
      const c0 = sqrtD((ctx.gamma * U.p + U.bt * U.bt) / U.rho), c1 = sqrtD((ctx.gamma * D.p + D.bt * D.bt) / D.rho);
      const w: Wave = sigma < 0
        ? { kind: "fast_fan", side: sigma, s_left: U.u - c0, s_right: D.u - c1, left: U, right: D, fan }
        : { kind: "fast_fan", side: sigma, s_left: D.u + c1, s_right: U.u + c0, left: D, right: U, fan };
      return [D, w];
    }
    return [D, null];
  }
  return [U, null];                                             // no fast wave on this side
}

const IDENTITY: Frame = { mirror: false, sB: 1, theta: 0.0, vtR: [0.0, 0.0], rho0: 1.0, p0: 1.0 };

function result(prob: RiemannProblem, retcode: RiemannSolution["retcode"], reason: string, ctx: Ctx,
  psi: number[], residual: number, waves: Wave[], check: RiemannSolution["check"]): RiemannSolution {
  return { retcode, reason, prob, psi, residual, waves, frame: IDENTITY, ctx, check, method: "quasi_euler" };
}

/**
 * Quasi-Euler solution for Bn ≈ 0 (fast waves and a tangential discontinuity), in the
 * user's units. Refuses only vacuum generation.
 */
export function solvePerpendicular(prob: RiemannProblem, opts: SolverOptions): RiemannSolution {
  const { L, R, gamma } = prob;
  const F = IDENTITY;
  const ctx = makeCtx(gamma, L[4]!, opts);
  const UL = toHState(L), UR = toHState(R);
  const PL = ptot(UL.rho, UL.p, UL.bt), PR = ptot(UR.rho, UR.p, UR.bt);
  const g = (P: number) => perpSide(UL, P, -1, ctx)[0].u - perpSide(UR, P, +1, ctx)[0].u;
  const lo = 1e-12 * Math.min(PL, PR);
  let hi = 2 * Math.max(PL, PR);
  const zero = [0, 0, 0, 0, 0];
  if (g(lo) < 0) return result(prob, "RegularLimit", "vacuum", ctx, zero, Infinity, [], null);
  let k = 0;
  while (g(hi) > 0) {
    hi *= 2; k++;
    if (k > 200) return result(prob, "NoConvergence", "no_convergence", ctx, zero, Infinity, [], null);
  }
  const Ps = findZeroBrent(g, lo, hi, { xatol: 0.0, xrtol: 4 * EPS });
  const [SL, wl] = perpSide(UL, Ps, -1, ctx, true);
  const [SR, wr] = perpSide(UR, Ps, +1, ctx, true);
  const us = (SL.u + SR.u) / 2;
  const waves: Wave[] = [];
  if (wl !== null) waves.push(wl);
  waves.push({ kind: "tangential", side: 0, s_left: us, s_right: us, left: SL, right: SR, fan: null });
  if (wr !== null) waves.push(wr);
  const btmax = Math.max(UL.bt, UR.bt);
  const r = Math.abs(SL.u - SR.u) /
    Math.max(sqrtD((gamma * Math.max(L[7]!, R[7]!) + btmax * btmax) / Math.min(L[0]!, R[0]!)), EPS);
  const chk = checkWaves(waves, ctx, F, gamma);
  const lim = Math.min(SL.rho, SR.rho, SL.p, SR.p) < opts.rho_floor * Math.min(L[0]!, R[0]!, L[7]!, R[7]!) ? "vacuum" : "none";
  const psi = [Ps, 0.0, 0.0, 0.0, 0.0];
  if (lim !== "none") return result(prob, "RegularLimit", lim, ctx, psi, r, waves, chk);
  if (!chk.ok) return result(prob, "CheckFailed", "check", ctx, psi, r, waves, chk);
  return result(prob, "Success", "none", ctx, psi, r, waves, chk);
}
