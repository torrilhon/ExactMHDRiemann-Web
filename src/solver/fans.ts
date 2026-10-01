// Fast and slow rarefaction fans: integral curves of the fast/slow eigenvectors.
// Along a fan ρ = ρ0 e^{-s}, p = p0 e^{-γ s} (paper eq. 53); u, Bt and the transverse
// velocity follow ODEs (paper eq. 54). Mirrors fans.jl.
//
// Fast fan: parametrized by ℓ = log(Bt/Bt0) = ψ τ, τ ∈ [0,1], ψ ≤ 0. Then
//   ds/dℓ = cA²/cf² - 1 ∈ (-1, 0),  du/dℓ = -σ cf ds/dℓ,  dw/dℓ = -σ Bt Bn/(ρ cf).
// Slow fan: parametrized by the arc length ς of the curve (s, Bt/√p0),
// ς = ς_end τ with ς_end = s_vac tanh(-ψ/s_vac); in terms of s,
//   du/ds = -σ cs,  dBt/ds = Bt cs²/(cA² - cs²),  dw/ds = -σ cs Bt cA²/(Bn (cA² - cs²)).
// w is the transverse velocity change along the (fixed) direction of Bt.

import { speeds, speedsH } from "./eos.ts";
import { DomainError, sqrtD } from "./math.ts";
import type { OdeSolution, Rhs } from "./rk/engine.ts";
import { integrate } from "./rk/engine.ts";
import { VERN9 } from "./rk/vern9.ts";
import type { Ctx, FanData, HState } from "./types.ts";
import { hstate } from "./types.ts";

export function fastFanRhs(rho0: number, p0: number, bt0: number, psi: number, gamma: number, Bn: number, sigma: number): Rhs {
  return (y, tau) => {
    const s = y[0]!;
    const rho = rho0 * Math.exp(-s);
    const p = p0 * Math.exp(-gamma * s);
    const bt = bt0 * Math.exp(psi * tau);
    const [cf, cA] = speeds(rho, p, bt, Bn, gamma);
    const dsdl = (cA * cA) / (cf * cf) - 1;
    return [psi * dsdl, psi * (-sigma * cf * dsdl), psi * ((-sigma * bt * Bn) / (rho * cf))];
  };
}

export function slowFanRhs(rho0: number, p0: number, send: number, gamma: number, Bn: number, sigma: number): Rhs {
  return (y) => {
    const s = y[0]!, bt = y[2]!;
    const rho = rho0 * Math.exp(-s);
    const p = p0 * Math.exp(-gamma * s);
    const [cf, cA, cs] = speeds(rho, p, bt, Bn, gamma);
    let Bp: number, Wp: number;
    if (cA * cA < (gamma * p) / rho) {
      // a > cA: cs → cA as Bt → 0; cA² - cs² = cA² bt²/(ρ (cf² - cA²)) avoids the cancellation
      const g = rho * (cf * cf - cA * cA);
      Bp = ((cs * cs) * g) / ((cA * cA) * bt);                // dBt/ds
      Wp = (-sigma * cs * g) / (Bn * bt);                     // dw/ds
    } else {
      const den = cA * cA - cs * cs;                           // cs < a ≤ cA: well separated
      Bp = (bt * (cs * cs)) / den;
      Wp = (-sigma * cs * bt * (cA * cA)) / (Bn * den);
    }
    const q = Bp / sqrtD(p0);
    const N = sqrtD(1 + q * q);                                // dς/ds
    return [send * (1 / N), send * ((-sigma * cs) / N), send * (Bp / N), send * (Wp / N)];
  };
}

/** Julia `_fan_solve`: Vern9 on τ ∈ [0, 1]; an unsuccessful retcode becomes DomainError. */
export function fanSolve(rhs: Rhs, y0: number[], ctx: Ctx, dense: boolean, par: number): OdeSolution {
  const opts = { abstol: ctx.opts.ode_tol, reltol: ctx.opts.ode_tol, dense, maxiters: 100_000 };
  // the residual (dense = false) integrates through the step tape of the Jacobian
  const sol = !dense && ctx.tape !== undefined
    ? ctx.tape.integrate(VERN9, rhs, y0, 0.0, 1.0, opts)
    : integrate(VERN9, rhs, y0, 0.0, 1.0, opts);
  if (!sol.successful) throw new DomainError(par, `fan ODE failed: ${sol.retcode}`);
  return sol;
}

/** Fast fan on U. Returns (downstream, FanData or null). */
export function fastFan(U: HState, psi: number, sigma: number, ctx: Ctx, record = false): [HState, FanData | null] {
  const rhs = fastFanRhs(U.rho, U.p, U.bt, psi, ctx.gamma, ctx.Bn, sigma);
  const sol = fanSolve(rhs, [0, U.u, 0], ctx, record, psi);
  const [s, u, w] = sol.end as [number, number, number];
  const e0 = Math.cos(U.phi), e1 = Math.sin(U.phi);
  const D = hstate(U.rho * Math.exp(-s), u, U.p * Math.exp(-ctx.gamma * s), U.bt * Math.exp(psi), U.phi,
    [U.vt[0] + w * e0, U.vt[1] + w * e1]);
  const fan: FanData | null = record ? { family: "fast", sigma, up: U, par: psi, sol } : null;
  return [D, fan];
}

export function slowFanSend(psi: number, ctx: Ctx): number {
  return ctx.opts.s_vac * Math.tanh(-psi / ctx.opts.s_vac);
}

/** Slow fan on U. Returns (downstream, FanData or null). */
export function slowFan(U: HState, psi: number, sigma: number, ctx: Ctx, record = false): [HState, FanData | null] {
  const send = slowFanSend(psi, ctx);
  const rhs = slowFanRhs(U.rho, U.p, send, ctx.gamma, ctx.Bn, sigma);
  const sol = fanSolve(rhs, [0, U.u, U.bt, 0], ctx, record, send);
  const [s, u, bt, w] = sol.end as [number, number, number, number];
  const e0 = Math.cos(U.phi), e1 = Math.sin(U.phi);
  const D = hstate(U.rho * Math.exp(-s), u, U.p * Math.exp(-ctx.gamma * s), bt, U.phi,
    [U.vt[0] + w * e0, U.vt[1] + w * e1]);
  const fan: FanData | null = record ? { family: "slow", sigma, up: U, par: send, sol } : null;
  return [D, fan];
}

/** State inside a fan at parameter τ ∈ [0, 1]. */
export function fanState(f: FanData, tau: number, ctx: Ctx): HState {
  const U = f.up;
  const y = f.sol.at(tau);
  const e0 = Math.cos(U.phi), e1 = Math.sin(U.phi);
  if (f.family === "fast") {
    const [s, u, w] = y as [number, number, number];
    return hstate(U.rho * Math.exp(-s), u, U.p * Math.exp(-ctx.gamma * s), U.bt * Math.exp(f.par * tau), U.phi,
      [U.vt[0] + w * e0, U.vt[1] + w * e1]);
  } else if (f.family === "fast0") {                      // quasi-Euler fan (Bn = 0): s = s* τ, Bt ∝ ρ
    const s = f.par * tau;
    return hstate(U.rho * Math.exp(-s), y[0]!, U.p * Math.exp(-ctx.gamma * s), U.bt * Math.exp(-s), U.phi, U.vt);
  }
  const [s, u, bt, w] = y as [number, number, number, number];
  return hstate(U.rho * Math.exp(-s), u, U.p * Math.exp(-ctx.gamma * s), bt, U.phi,
    [U.vt[0] + w * e0, U.vt[1] + w * e1]);
}

/** Characteristic speed of the fan family at τ. */
export function fanSpeed(f: FanData, tau: number, ctx: Ctx): number {
  const h = fanState(f, tau, ctx);
  const [cf, , cs] = speedsH(h, ctx);
  return h.u + f.sigma * (f.family === "slow" ? cs : cf);
}
