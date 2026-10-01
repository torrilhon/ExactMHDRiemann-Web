// The residual as the nonlinear solver sees it (solve.jl: BIG, safe_residual,
// problem_scale, resnorm) and its finite-difference Jacobian (plan §4.3).

import { speeds } from "./eos.ts";
import { fdJacobian } from "./jacobian.ts";
import type { Mat } from "./linalg.ts";
import { DomainError, sqrtD } from "./math.ts";
import { StepTape } from "./rk/tape.ts";
import { residual, SLOW_EPS } from "./side.ts";
import type { Ctx, HState } from "./types.ts";

export const BIG = 1e6;
const BIGV = [BIG, BIG, BIG, BIG, BIG];

/**
 * Residual that never throws for numerical reasons: domain errors and non-finite values
 * map to the BIG vector, so the trust region simply rejects the step. Other exceptions
 * (programming errors, ConvergenceFailed, ArgumentError) propagate, as in Julia.
 */
export function safeResidual(Psi: readonly number[], UL: HState, UR: HState, ctx: Ctx): number[] {
  try {
    const r = residual(Psi, UL, UR, ctx);
    if (r.every(Number.isFinite)) return r;
  } catch (e) {
    if (!(e instanceof DomainError)) throw e;
  }
  return BIGV.slice();
}

export const isBig = (r: readonly number[]) => r.every((v) => v === BIG);

/** Natural scales (total pressure, fast speed, field strength) of the two input states. */
export function problemScale(UL: HState, UR: HState, Bn: number, gamma: number): [number, number, number] {
  const P = Math.max(UL.p + (Bn * Bn + UL.bt * UL.bt) / 2, UR.p + (Bn * Bn + UR.bt * UR.bt) / 2);
  const C = Math.max(speeds(UL.rho, UL.p, UL.bt, Bn, gamma)[0], speeds(UR.rho, UR.p, UR.bt, Bn, gamma)[0]);
  const bmax = Math.max(UL.bt, UR.bt);
  const B = sqrtD(Bn * Bn + bmax * bmax);
  return [P, C, B];
}

export function resnorm(Psi: readonly number[], UL: HState, UR: HState, ctx: Ctx): number {
  let m = 0;
  for (const v of safeResidual(Psi, UL, UR, ctx)) {
    const a = Math.abs(v);
    if (Number.isNaN(a)) return NaN;
    m = Math.max(m, a);
  }
  return m;
}

/** branch points of build_side per component of Ψ = (ψf⁻, ψs⁻, αR, ψs⁺, ψf⁺) */
const SWITCHES: readonly (readonly number[])[] = [[0], [SLOW_EPS], [], [SLOW_EPS], [0]];

/**
 * Finite-difference Jacobian of the safe residual at Ψ: branch-aware stencils, BIG
 * points treated as invalid, and frozen fan-ODE steps (the base point is evaluated once
 * more to record its step sequences; every perturbed point replays them).
 */
export function residualJacobian(Psi: readonly number[], UL: HState, UR: HState, ctx: Ctx): Mat {
  const tape = new StepTape();
  const c: Ctx = { ...ctx, tape };
  tape.record();
  const r0 = safeResidual(Psi, UL, UR, c);
  const { J } = fdJacobian((x) => safeResidual(x, UL, UR, c), Psi, r0, {
    switches: SWITCHES,
    invalid: isBig,
    beforePerturbed: () => tape.replay(),
  });
  return J;
}
