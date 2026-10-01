// Canonical frame: larger Bt/√p on the left, Bn > 0, left Bt along +y, right vt = 0,
// units ρ_L = p_L = 1. Every step is an exact symmetry of ideal MHD, so the map is
// invertible. Mirrors canonical.jl.
//
// Side swap: reflection x → -x (swap L and R, vx → -vx, Bx → -Bx) combined with field
// reversal B → -B gives W'(ξ) = T(W(-ξ)) with T(ρ, vx, vy, vz, Bx, By, Bz, p) =
// (ρ, -vx, vy, vz, Bx, -By, -Bz, p): Bn keeps its sign and the sides trade places.

import { atan2, hypot, sqrtD } from "./math.ts";
import type { Ctx, Frame, HState, Prim, RetCode, RiemannProblem, SolverOptions } from "./types.ts";
import { hstate } from "./types.ts";

/** rot2(θ) * v for the 2×2 rotation [cos θ  -sin θ; sin θ  cos θ]. */
function rot(theta: number, v0: number, v1: number): [number, number] {
  const c = Math.cos(theta), s = Math.sin(theta);
  return [c * v0 + -s * v1, s * v0 + c * v1];
}

export function mirrorT(W: Prim): number[] {
  return [W[0]!, -W[1]!, W[2]!, W[3]!, W[4]!, -W[5]!, -W[6]!, W[7]!];
}

export function btp(W: Prim): number {
  return hypot(W[5]!, W[6]!) / sqrtD(W[7]!);
}

export function makeFrame(L0: Prim, R0: Prim): Frame {
  const mirror = btp(L0) < btp(R0);
  const [L, R] = mirror ? [mirrorT(R0), mirrorT(L0)] : [L0, R0];
  const sB = L[4]! < 0 ? -1 : 1;
  const theta = atan2(sB * L[6]!, sB * L[5]!);
  const vtR = rot(-theta, R[2]!, R[3]!);
  return { mirror, sB, theta, vtR, rho0: L[0]!, p0: L[7]! };
}

/** Canonical 8-vectors (left, right) of a problem in the frame F. */
export function canonicalStates(L: Prim, R: Prim, F: Frame): [number[], number[]] {
  if (F.mirror) [L, R] = [mirrorT(R), mirrorT(L)];
  return [toCanonical(L, F), toCanonical(R, F)];
}

/** The same frame without the side swap (used by the independent checks). */
export function unmirrored(F: Frame): Frame {
  return { ...F, mirror: false };
}

/** Primitive 8-vector of the (already side-swapped) problem → canonical 8-vector (`_to_canonical`). */
export function toCanonical(W: Prim, F: Frame): number[] {
  const r = rot(-F.theta, W[2]!, W[3]!);
  const vt = [r[0] - F.vtR[0], r[1] - F.vtR[1]];
  const bt = rot(-F.theta, F.sB * W[5]!, F.sB * W[6]!);
  const c0 = sqrtD(F.p0 / F.rho0), b0 = sqrtD(F.p0);
  return [W[0]! / F.rho0, W[1]! / c0, vt[0]! / c0, vt[1]! / c0, (F.sB * W[4]!) / b0, bt[0] / b0, bt[1] / b0, W[7]! / F.p0];
}

/** Canonical 8-vector → user-frame primitive 8-vector. */
export function fromCanonical(W: Prim, F: Frame): number[] {
  const c0 = sqrtD(F.p0 / F.rho0), b0 = sqrtD(F.p0);
  const vt = rot(F.theta, W[2]! * c0 + F.vtR[0], W[3]! * c0 + F.vtR[1]);
  const b = rot(F.theta, W[5]! * b0, W[6]! * b0);
  const bt = [F.sB * b[0], F.sB * b[1]];
  const U = [W[0]! * F.rho0, W[1]! * c0, vt[0], vt[1], F.sB * W[4]! * b0, bt[0]!, bt[1]!, W[7]! * F.p0];
  return F.mirror ? mirrorT(U) : U;
}

export function speedToUser(s: number, F: Frame): number {
  return (F.mirror ? -s : s) * sqrtD(F.p0 / F.rho0);
}

/** HState (canonical) → canonical primitive 8-vector. */
export function toPrim(h: HState, ctx: Ctx): number[] {
  return [h.rho, h.u, h.vt[0], h.vt[1], ctx.Bn, h.bt * Math.cos(h.phi), h.bt * Math.sin(h.phi), h.p];
}

/** Canonical primitive 8-vector → HState (angle in (-π, π]). */
export function toHState(W: Prim): HState {
  return hstate(W[0]!, W[1]!, W[7]!, hypot(W[5]!, W[6]!), atan2(W[6]!, W[5]!), [W[2]!, W[3]!]);
}

/** Checks of the v1 input domain, applied before any solve (`validate_input`). */
export function validateInput(prob: RiemannProblem, opts: SolverOptions): [RetCode, string] {
  const { L, R, gamma } = prob;
  if (!(L.every(Number.isFinite) && R.every(Number.isFinite) && Number.isFinite(gamma))) return ["InvalidInput", "nonfinite"];
  if (!(gamma > 1)) return ["InvalidInput", "gamma"];
  if (!(L[0]! > 0 && R[0]! > 0)) return ["InvalidInput", "density"];
  if (!(L[7]! > 0 && R[7]! > 0)) return ["InvalidInput", "pressure"];
  if (!(Math.abs(L[4]! - R[4]!) <= 1e-12 * Math.max(Math.abs(L[4]!), 1.0))) return ["InvalidInput", "bn_jump"];
  const rr = R[0]! / L[0]!, pr = R[7]! / L[7]!;
  if (!(1 / opts.ratio_max <= rr && rr <= opts.ratio_max && 1 / opts.ratio_max <= pr && pr <= opts.ratio_max)) {
    return ["Unsupported", "extreme_ratio"];
  }
  // vanishing normal field: quasi-Euler solver, any Bt (including 0) is fine
  if (Math.max(Math.abs(L[4]!) / sqrtD(L[7]!), Math.abs(R[4]!) / sqrtD(R[7]!)) <= opts.bn_euler) return ["Success", "quasi_euler"];
  // transverse field: one side may be very small as long as the other is not
  const btL = btp(L), btR = btp(R);
  const tolr = 1 - 1e-12;       // a value set exactly at a threshold must not be refused by round-off
  if (!(Math.min(btL, btR) >= tolr * opts.bt_min_smaller && Math.max(btL, btR) >= tolr * opts.bt_min_larger)) {
    return ["Unsupported", "switch_on_off"];
  }
  return ["Success", "none"];
}
