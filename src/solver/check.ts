// Independent checks of an assembled solution, in the user frame, using the full
// conservative flux and a separate eigenvector integration for fans. They share no
// code with the wave kernels except the characteristic speeds. Mirrors check.jl,
// including the noise floor of the pointwise fan defect (Julia commit 0ca567a).

import { EigenvalueDecomposition, Matrix } from "ml-matrix";
import { fromCanonical, speedToUser, toPrim } from "./canonical.ts";
import { conserved, flux, primitiveJacobian, speedsW } from "./eos.ts";
import { fanSpeed, fanState } from "./fans.ts";
import { juliaRepr } from "./format.ts";
import { hypot } from "./math.ts";
import { opnormInf } from "./linalg.ts";
import { integrate } from "./rk/engine.ts";
import { VERN9 } from "./rk/vern9.ts";
import type { CheckReport, Ctx, FanData, Frame, HState, Prim, Wave, WaveKind } from "./types.ts";

export function userState(h: HState, ctx: Ctx, F: Frame): number[] {
  return fromCanonical(toPrim(h, ctx), F);
}

/** Region index of a relative flow speed w: 1 super-fast … 4 sub-slow. */
export function laxRegion(w: number, cf: number, cA: number, cs: number, tol = 1e-9): number {
  if (w > cf * (1 + tol)) return 1;
  if (w > cA * (1 + tol)) return 2;
  if (w > cs * (1 + tol)) return 3;
  return 4;
}

/** Eigenvector of family `fam` (1-based index into the sorted real parts), normalized by r[1]. */
export function fanEigenRhs(q: readonly number[], bx: number, gamma: number, fam: number): number[] {
  const W = [q[0]!, q[1]!, q[2]!, q[3]!, bx, q[4]!, q[5]!, q[6]!];
  const E = new EigenvalueDecomposition(new Matrix(primitiveJacobian(W, gamma)));
  const vals = E.realEigenvalues;
  // sortperm(real.(E.values))[fam]: stable sort, numeric comparator
  const idx = vals.map((v, i) => [v, i] as const).sort((a, b) => a[0] - b[0])[fam - 1]![1];
  const V = E.eigenvectorMatrix;
  const r = Array.from({ length: 7 }, (_, i) => V.get(i, idx));
  return r.map((v) => v / r[0]!);
}

/**
 * Integrate the fan family `fam` (1, 3, 5 or 7) from state W0 to density ρ1. Julia uses
 * Vern7 at tol 1e-11 without checking the retcode; the port uses the Vern9 engine at the
 * same tolerance (plan §4.1). A step whose stages leave the domain of the eigenvector
 * field (NaN) is retried with a smaller step, where Julia's Vern7 happens to take steps
 * that stay inside (perpendicular scan problem 280).
 */
export function fanEigenIntegrate(W0: Prim, rho1: number, fam: number, gamma: number): number[] {
  const q0 = [W0[0]!, W0[1]!, W0[2]!, W0[3]!, W0[5]!, W0[6]!, W0[7]!];
  const s = integrate(VERN9, (q) => fanEigenRhs(q, W0[4]!, gamma, fam), q0, W0[0]!, rho1,
    { abstol: 1e-11, reltol: 1e-11, retryNaN: true });
  const q = s.end;
  return [q[0]!, q[1]!, q[2]!, q[3]!, W0[4]!, q[4]!, q[5]!, q[6]!];
}

/**
 * Distance of a fan's family speed to the nearest other characteristic speed, relative
 * to the fast speed. For the slow family both the Alfvén speed and the contact count.
 */
export function speedGap(W: Prim, gamma: number, kind: WaveKind): number {
  const [cf, cA, cs] = speedsW(W, gamma);
  return kind === "fast_fan" ? (cf - cA) / cf : Math.min(cA - cs, cs) / cf;
}

const prim7 = (W: Prim) => [W[0]!, W[1]!, W[2]!, W[3]!, W[5]!, W[6]!, W[7]!];
const maxAbs = (v: readonly number[]) => v.reduce((m, x) => Math.max(m, Math.abs(x)), 0);

/** Largest relative defect of the eigen relation (A - λI) dq/dτ = 0 inside a fan. */
export function fanPointwiseDefect(f: FanData, ctx: Ctx, F: Frame, gamma: number): number {
  let d = 0.0;
  const h = 1e-3;
  const q = (t: number) => prim7(userState(fanState(f, t, ctx), ctx, F));
  for (const tau of [0.02, 0.1, 0.25, 0.4, 0.55, 0.7, 0.85, 0.98]) {
    const W = userState(fanState(f, tau, ctx), ctx, F);
    const q2p = q(tau + 2 * h), q1p = q(tau + h), q1m = q(tau - h), q2m = q(tau - 2 * h);
    const dq = q2p.map((_, i) => (-q2p[i]! + 8 * q1p[i]! - 8 * q1m[i]! + q2m[i]!) / (12 * h));   // 4th-order difference
    const lambda = speedToUser(fanSpeed(f, tau, ctx), F);
    const A = primitiveJacobian(W, gamma);
    const nq = maxAbs(dq);
    if (!(nq > 0)) continue;
    // dq of a vanishingly weak fan is as small as the noise of the dense output divided
    // by h; measuring against a floor of 1e-6 of the state size keeps that noise from
    // being reported as a mismatch (the defect of such a fan is bounded by its strength)
    const nfloor = 1e-6 * Math.max(maxAbs(prim7(W)), 1.0);
    const Ad = A.map((row, i) => row.reduce((s, a, j) => s + a * dq[j]!, 0) - lambda * dq[i]!);
    d = Math.max(d, maxAbs(Ad) / (opnormInf(A) * Math.max(nq, nfloor)));
  }
  return d;
}

export function relerr(a: readonly number[], b: readonly number[]): number {
  let m = -Infinity;
  for (let i = 0; i < a.length; i++) {
    m = Math.max(m, Math.abs(a[i]! - b[i]!) / Math.max(Math.abs(a[i]!), Math.abs(b[i]!), 1.0));
  }
  return m;
}

const fmt = juliaRepr;

/** Julia `check_waves(waves, ctx, F, γ)`: returns (ok, maxerr, messages). */
export function checkWaves(waves: readonly Wave[], ctx: Ctx, F: Frame, gamma: number): CheckReport {
  const tol = ctx.opts.check_tol;
  const msgs: string[] = [];
  let maxerr = 0.0;
  const note = (ok: boolean, m: () => string) => { if (!ok) msgs.push(m()); };
  waves.forEach((wv, k) => {
    const i = k + 1;                                          // Julia's 1-based wave number in messages
    const Wl = userState(wv.left, ctx, F), Wr = userState(wv.right, ctx, F);
    const sl = speedToUser(wv.s_left, F), sr = speedToUser(wv.s_right, F);
    const sigma = wv.side;
    const [Wup, Wdn] = sigma < 0 ? [Wl, Wr] : [Wr, Wl];
    if (wv.kind === "fast_shock" || wv.kind === "slow_shock" || wv.kind === "rotation" || wv.kind === "contact" ||
        wv.kind === "tangential") {
      const s = sl;
      const Ur = conserved(Wr, gamma), Ul = conserved(Wl, gamma), Fr = flux(Wr, gamma), Fl = flux(Wl, gamma);
      const rh = Ur.map((_, j) => s * (Ur[j]! - Ul[j]!) - (Fr[j]! - Fl[j]!));
      const scale = Math.max(maxAbs(Fl), maxAbs(Fr), Math.abs(s) * maxAbs(Ul), 1.0);
      const e = maxAbs(rh) / scale;
      maxerr = Math.max(maxerr, e);
      note(e <= tol, () => `wave ${i} (${wv.kind}): Rankine-Hugoniot error ${fmt(e)}`);
    }
    const weak = Math.abs(Wdn[0]! - Wup[0]!) <= 1e-7 * Wup[0]!;      // vanishing wave: only RH applies
    if ((wv.kind === "fast_shock" || wv.kind === "slow_shock") && !weak) {
      const s = sl;
      const [cfu, cAu, csu] = speedsW(Wup, gamma), [cfd, cAd, csd] = speedsW(Wdn, gamma);
      const wu = sigma * (s - Wup[1]!), wd = sigma * (s - Wdn[1]!);
      const ru = laxRegion(wu, cfu, cAu, csu), rd = laxRegion(wd, cfd, cAd, csd);
      const want = wv.kind === "fast_shock" ? [1, 2] : [3, 4];
      note(ru === want[0] && rd === want[1], () => `wave ${i} (${wv.kind}): shock type ${ru}→${rd}, expected ${want[0]}→${want[1]}`);
      note(Wdn[0]! > Wup[0]!, () => `wave ${i} (${wv.kind}): density does not increase (entropy)`);
    } else if (wv.kind === "rotation") {
      let e = relerr([Wl[0]!, Wl[1]!, Wl[7]!, hypot(Wl[5]!, Wl[6]!)], [Wr[0]!, Wr[1]!, Wr[7]!, hypot(Wr[5]!, Wr[6]!)]);
      const [, cA] = speedsW(Wl, gamma);
      e = Math.max(e, Math.abs(sl - (Wl[1]! + sigma * cA)) / Math.max(Math.abs(sl), 1.0));
      maxerr = Math.max(maxerr, e);
      note(e <= tol, () => `wave ${i} (rotation): not a rotational discontinuity (${fmt(e)})`);
    } else if (wv.kind === "tangential") {
      // Bn = 0: u and total pressure continuous, everything else may jump
      const Pl = Wl[7]! + (Wl[5]! * Wl[5]! + Wl[6]! * Wl[6]!) / 2, Pr = Wr[7]! + (Wr[5]! * Wr[5]! + Wr[6]! * Wr[6]!) / 2;
      const e = Math.max(Math.abs(Wl[1]! - Wr[1]!) / Math.max(1.0, Math.abs(Wl[1]!)), Math.abs(Pl - Pr) / Math.max(Pl, Pr));
      maxerr = Math.max(maxerr, e);
      note(e <= tol, () => `tangential discontinuity: jump in u or total pressure of ${fmt(e)}`);
      note(Math.abs(sl - Wl[1]!) <= tol * Math.max(1.0, Math.abs(sl)), () => "tangential discontinuity: speed differs from flow speed");
    } else if (wv.kind === "contact") {
      const e = relerr(Wl.slice(1, 8), Wr.slice(1, 8));
      maxerr = Math.max(maxerr, e);
      note(e <= tol, () => `contact: jump in (v, B, p) of ${fmt(e)}`);
      note(Math.abs(sl - Wl[1]!) <= tol * Math.max(1.0, Math.abs(sl)), () => "contact: speed differs from flow speed");
    } else if (wv.kind === "fast_fan" || wv.kind === "slow_fan") {
      const fam = wv.kind === "fast_fan" ? (sigma < 0 ? 1 : 7) : (sigma < 0 ? 3 : 5);
      const [cfu, , csu] = speedsW(Wup, gamma), [cfd, , csd] = speedsW(Wdn, gamma);
      const [cu, cd] = wv.kind === "fast_fan" ? [cfu, cfd] : [csu, csd];
      const [shead, stail] = sigma < 0 ? [sl, sr] : [sr, sl];
      let e = Math.max(Math.abs(shead - (Wup[1]! + sigma * cu)), Math.abs(stail - (Wdn[1]! + sigma * cd))) /
        Math.max(Math.abs(shead), Math.abs(stail), 1.0);
      note(weak || sl <= sr + tol * Math.max(1.0, Math.abs(sr)), () => `wave ${i} (${wv.kind}): fan speeds not increasing`);
      if (Math.abs(Wdn[0]! - Wup[0]!) > 1e-12 * Wup[0]! && wv.fan !== null) {
        // (a) pointwise: along the fan, (A(q) - λ I) dq/dτ = 0 with λ the family speed
        e = Math.max(e, fanPointwiseDefect(wv.fan, ctx, F, gamma));
        // (b) where the family is well separated: independent integration along
        //     eigenvectors of the numerical Jacobian, compared at the fan's end
        if (speedGap(Wup, gamma, wv.kind) > 1e-3 && speedGap(Wdn, gamma, wv.kind) > 1e-3) {
          e = Math.max(e, relerr(fanEigenIntegrate(Wup, Wdn[0]!, fam, gamma), Wdn));
        }
      }
      maxerr = Math.max(maxerr, e);
      note(e <= 100 * tol, () => `wave ${i} (${wv.kind}): fan mismatch ${fmt(e)}`);
    }
    if (k > 0) {
      const prev = waves[k - 1]!;
      const ok = speedToUser(prev.s_right, F) <= sl + tol * Math.max(1.0, Math.abs(sl));
      note(ok, () => `waves ${i - 1} and ${i} overlap`);
    }
  });
  return { ok: msgs.length === 0, maxerr, messages: msgs };
}
