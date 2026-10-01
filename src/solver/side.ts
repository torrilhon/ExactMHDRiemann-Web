// The single wave-sequence builder. Residual, output and checks all use it, so what
// Newton solves and what is reported cannot differ. Mirrors side.jl.

import { speedsH } from "./eos.ts";
import { fastFan, slowFan } from "./fans.ts";
import { sqrtD } from "./math.ts";
import { fastShock, slowLimit, slowShockState } from "./shocks.ts";
import type { Ctx, FanData, HState, Wave, WaveKind } from "./types.ts";
import { hstate } from "./types.ts";

/** Alfvén (rotational) discontinuity turning Bt from U.φ to α (requires Bn > 0). */
export function rotation(U: HState, alpha: number, sigma: number, ctx: Ctx): [HState, number] {
  const e00 = Math.cos(U.phi), e01 = Math.sin(U.phi);
  const e10 = Math.cos(alpha), e11 = Math.sin(alpha);
  const f = -sigma * U.bt;
  const r = sqrtD(U.rho);
  const D = hstate(U.rho, U.u, U.p, U.bt, alpha, [U.vt[0] + (f * (e10 - e00)) / r, U.vt[1] + (f * (e11 - e01)) / r]);
  const [, cA] = speedsH(U, ctx);
  return [D, U.u + sigma * cA];
}

/** slow ψ below this is treated as a (vanishing) fan */
export const SLOW_EPS = 1e-9;

/** Store a wave in physical left-right order (`_push!`). */
function push(waves: Wave[], kind: WaveKind, sigma: number, up: HState, down: HState, sUp: number, sDown: number,
  fan: FanData | null = null): void {
  if (sigma < 0) waves.push({ kind, side: sigma, s_left: sUp, s_right: sDown, left: up, right: down, fan });
  else waves.push({ kind, side: sigma, s_left: sDown, s_right: sUp, left: down, right: up, fan });
}

/**
 * Apply fast wave, rotation to angle αR and slow wave to the outer state U of side σ
 * (-1 left, +1 right). With record = true the waves are returned in outer-to-inner
 * order, otherwise null.
 */
export function buildSide(U: HState, psif: number, psis: number, alphaR: number, sigma: number, ctx: Ctx,
  record = false): [HState, Wave[] | null] {
  const waves: Wave[] | null = record ? [] : null;
  // fast wave
  let F: HState;
  if (psif > 0) {
    let s: number;
    [F, s] = fastShock(U, psif, sigma, ctx);
    if (waves) push(waves, "fast_shock", sigma, U, F, s, s);
  } else {
    let fan: FanData | null;
    [F, fan] = fastFan(U, psif, sigma, ctx, record);
    if (waves) {
      const [cf0] = speedsH(U, ctx), [cf1] = speedsH(F, ctx);
      push(waves, "fast_fan", sigma, U, F, U.u + sigma * cf0, F.u + sigma * cf1, fan);
    }
  }
  // rotation
  const [Rt, sr] = rotation(F, alphaR, sigma, ctx);
  if (waves) push(waves, "rotation", sigma, F, Rt, sr, sr);
  // slow wave
  let S: HState;
  if (psis > SLOW_EPS) {
    const Dmax = (1 - ctx.opts.delta) * slowLimit(Rt, ctx);
    const Delta = Dmax * Math.tanh(psis / Dmax);
    let s: number;
    [S, s] = slowShockState(Rt, Delta, sigma, ctx);
    if (waves) push(waves, "slow_shock", sigma, Rt, S, s, s);
  } else {
    let fan: FanData | null;
    [S, fan] = slowFan(Rt, psis, sigma, ctx, record);
    if (waves) {
      const cs0 = speedsH(Rt, ctx)[2], cs1 = speedsH(S, ctx)[2];
      push(waves, "slow_fan", sigma, Rt, S, Rt.u + sigma * cs0, S.u + sigma * cs1, fan);
    }
  }
  return [S, waves];
}

/**
 * Mismatch of (p, u, |Bt|, vt) between the inner states of the two sides, in units of
 * the problem's own scales. Ψ = (ψf⁻, ψs⁻, αR, ψs⁺, ψf⁺).
 */
export function residual(Psi: readonly number[], UL: HState, UR: HState, ctx: Ctx): number[] {
  const [SL] = buildSide(UL, Psi[0]!, Psi[1]!, Psi[2]!, -1, ctx);
  const [SR] = buildSide(UR, Psi[4]!, Psi[3]!, Psi[2]!, +1, ctx);
  const [P, C, B] = ctx.scale;
  return [(SL.p - SR.p) / P, (SL.u - SR.u) / C, (SL.bt - SR.bt) / B,
    (SL.vt[0] - SR.vt[0]) / C, (SL.vt[1] - SR.vt[1]) / C];
}

/** Assemble all waves of a solution, left to right, including the contact. */
export function assemble(Psi: readonly number[], UL: HState, UR: HState, ctx: Ctx): Wave[] {
  const [SL, wl] = buildSide(UL, Psi[0]!, Psi[1]!, Psi[2]!, -1, ctx, true);
  const [SR, wr] = buildSide(UR, Psi[4]!, Psi[3]!, Psi[2]!, +1, ctx, true);
  const uc = (SL.u + SR.u) / 2;
  return [...wl!, { kind: "contact", side: 0, s_left: uc, s_right: uc, left: SL, right: SR, fan: null }, ...wr!.reverse()];
}
