// Driver: input checks, canonical frame, nonlinear solve with a small fixed set of
// starts and a homotopy fallback, limit detection and independent checks.
// Mirrors solve.jl.

import { canonicalStates, makeFrame, toHState, unmirrored, validateInput } from "./canonical.ts";
import { checkWaves } from "./check.ts";
import { ArgumentError, DomainError, logD, sqrtD } from "./math.ts";
import { solvePerpendicular } from "./perpendicular.ts";
import { problemScale, residualJacobian, resnorm, safeResidual } from "./residual.ts";
import { assemble, SLOW_EPS } from "./side.ts";
import { slowLimit } from "./shocks.ts";
import { trustRegion } from "./trustregion.ts";
import { FanCache } from "./rk/tape.ts";
import type { CheckReport, Ctx, Frame, HState, RetCode, RiemannProblem, RiemannSolution, SolverOptions, Wave } from "./types.ts";
import { defaultOptions, hstate, makeCtx } from "./types.ts";

const ZERO5 = [0, 0, 0, 0, 0];

/** Errors that indicate a programming error and must surface (Julia's rethrow list in nl_solve). */
function isProgrammingError(e: unknown): boolean {
  return e instanceof ArgumentError || e instanceof TypeError || e instanceof RangeError ||
    e instanceof ReferenceError || e instanceof SyntaxError || !(e instanceof Error);
}

/** Trust-region solve from Ψ0; returns (Ψ, ‖residual‖∞), with Inf on numerical breakdown. */
export function nlSolve(Psi0: readonly number[], UL: HState, UR: HState, ctx: Ctx): [number[], number] {
  let u: number[];
  try {
    const s = trustRegion((x) => safeResidual(x, UL, UR, ctx), (x) => residualJacobian(x, UL, UR, ctx), Psi0,
      { abstol: ctx.opts.abstol, maxiters: ctx.opts.maxiters });
    u = s.u;
  } catch (e) {
    // numerical breakdown only; programming errors must surface
    if (isProgrammingError(e)) throw e;
    return [Psi0.slice(), Infinity];
  }
  return [u, resnorm(u, UL, UR, ctx)];
}

export function accept(r: number, ctx: Ctx): boolean {
  return r <= Math.max(1e3 * ctx.opts.abstol, 1e-10);
}

/** Starting guesses in the canonical frame (α = twist angle of the right state). */
export function starts(alpha: number): number[][] {
  const g = [[0.01, 0.01, alpha / 2, 0.01, 0.01]];
  if (Math.abs(alpha) > Math.PI - 1e-6) {                     // coplanar: α/2 = ±π/2 is a poor start
    g.push([0.01, 0.01, 0.0, 0.01, 0.01], [0.01, 0.01, Math.PI, 0.01, 0.01]);
  }
  g.push([0.01, 0.01, alpha / 2 + Math.PI, 0.01, 0.01]);
  // slow fans as start: on a side with tiny Bt the slow-shock branch is capped at
  // ΔB̂t ≤ A and its path variable saturates at once, so Newton needs ψs ≤ 0 there
  g.push([0.01, -0.01, alpha / 2, -0.01, 0.01], [0.01, 0.01, alpha / 2, -0.01, 0.01]);
  return g;
}

const now = (): number => (typeof performance !== "undefined" ? performance.now() : Date.now()) / 1000;

/**
 * Homotopy from the trivial problem (right state = left state) to the real one.
 * Returns (Ψ, residual, λ, URλ): on failure λ < 1 is the last value reached.
 */
export function homotopy(UL: HState, UR: HState, ctx: Ctx): [number[], number, number, HState] {
  const lUL = logD(UL.rho), lUR = logD(UR.rho), lpL = logD(UL.p), lpR = logD(UR.p);
  const Rl = (l: number) => hstate(Math.exp((1 - l) * lUL + l * lUR), (1 - l) * UL.u + l * UR.u,
    Math.exp((1 - l) * lpL + l * lpR), (1 - l) * UL.bt + l * UR.bt, l * UR.phi,
    [(1 - l) * UL.vt[0] + l * UR.vt[0], (1 - l) * UL.vt[1] + l * UR.vt[1]]);
  let lam = 0.0, h = 0.05;
  let Psi = ZERO5.slice(), PsiOld = ZERO5.slice(), lamOld = 0.0;
  const t0 = now();
  for (let it = 0; it < ctx.opts.homotopy_maxsteps; it++) {
    if (now() - t0 > ctx.opts.time_limit) break;
    const lamN = Math.min(lam + h, 1.0);
    // linear predictor from the last two points
    const pred = lam > 0
      ? Psi.map((v, i) => v + (v - PsiOld[i]!) * ((lamN - lam) / (lam - lamOld)))
      : Psi.map((v, i) => v + [1e-4, 1e-4, 0.0, 1e-4, 1e-4][i]!);
    const [PsiN, r] = nlSolve(pred, UL, Rl(lamN), ctx);
    if (accept(r, ctx)) {
      [PsiOld, lamOld, Psi, lam] = [Psi, lam, PsiN, lamN];
      if (lam >= 1.0) return [Psi, r, 1.0, Rl(1.0)];
      h = Math.min(2 * h, 0.25);
    } else {
      h /= 2;
      if (h < 1e-4) break;
    }
  }
  return [Psi, Infinity, lam, Rl(lam)];
}

interface FailedOpts {
  frame?: Frame | null; ctx?: Ctx | null; psi?: number[]; r?: number; waves?: Wave[]; method?: string; chk?: CheckReport | null;
}

function failed(prob: RiemannProblem, rc: RetCode, reason: string, o: FailedOpts = {}): RiemannSolution {
  return { retcode: rc, reason, prob, psi: o.psi ?? ZERO5.slice(), residual: o.r ?? Infinity, waves: o.waves ?? [],
    frame: o.frame ?? null, ctx: o.ctx ?? null, check: o.chk ?? null, method: o.method ?? "none" };
}

/** Regular-limit diagnostics of a converged Ψ. */
export function limitReason(Psi: readonly number[], waves: readonly Wave[], _UL: HState, _UR: HState, ctx: Ctx): string {
  const o = ctx.opts;
  for (const w of waves) {
    for (const h of [w.left, w.right]) {
      if (h.rho < o.rho_floor || h.p < o.rho_floor) return "vacuum";
      if (h.bt / sqrtD(h.p) < o.bt_floor) {
        return waves.some((x) => x.kind === "fast_fan") ? "fast_switch_off" : "bt_small";
      }
    }
  }
  for (const [psis, sigma] of [[Psi[1]!, -1], [Psi[3]!, +1]] as const) {
    if (psis > SLOW_EPS) {
      const rot = waves.find((x) => x.kind === "rotation" && x.side === sigma)!;
      const Rt = sigma < 0 ? rot.right : rot.left;
      const Dmax = (1 - o.delta) * slowLimit(Rt, ctx);
      if (psis / Dmax > o.sat) return "slow_limit";
    } else if (-psis / o.s_vac > o.sat) {
      return "vacuum";
    }
  }
  return "none";
}

export interface SolveArgs {
  opts?: Partial<SolverOptions>;
  /** 5 path variables in the canonical frame, tried first */
  guess?: readonly number[];
}

/**
 * Solve the Riemann problem. Never throws for physical reasons; inspect `retcode`.
 */
export function solve(prob: RiemannProblem, args: SolveArgs = {}): RiemannSolution {
  const opts = defaultOptions(args.opts);
  const [rc, reason] = validateInput(prob, opts);
  if (rc !== "Success") return failed(prob, rc, reason);
  if (reason === "quasi_euler") {
    try {
      return solvePerpendicular(prob, opts);
    } catch (e) {
      if (!(e instanceof DomainError)) throw e;
      return failed(prob, "NoConvergence", "no_convergence", { method: "quasi_euler" });
    }
  }
  const F = makeFrame(prob.L, prob.R);
  const [Lc, Rc] = canonicalStates(prob.L, prob.R, F);
  const UL = toHState(Lc), UR = toHState(Rc);
  const ctx: Ctx = { ...makeCtx(prob.gamma, Lc[4]!, opts, problemScale(UL, UR, Lc[4]!, prob.gamma)), fanCache: new FanCache() };
  if (Math.max(...Lc.map((v, i) => Math.abs(v - Rc[i]!))) <= 1e-14) {
    const w: Wave[] = [{ kind: "contact", side: 0, s_left: UL.u, s_right: UL.u, left: UL, right: UR, fan: null }];
    return { retcode: "Success", reason: "trivial", prob, psi: ZERO5.slice(), residual: 0.0, waves: w, frame: F, ctx,
      check: checkWaves(w, ctx, unmirrored(F), prob.gamma), method: "trivial" };
  }
  const alpha = UR.phi;
  const cand = args.guess === undefined ? starts(alpha) : [args.guess.slice(), ...starts(alpha)];
  let Psi = ZERO5.slice(), r = Infinity, method = "none";
  for (const Psi0 of cand) {
    [Psi, r] = nlSolve(Psi0, UL, UR, ctx);
    if (accept(r, ctx)) {
      method = "direct";
      break;
    }
  }
  if (!accept(r, ctx) && opts.homotopy) {
    let lam: number, URl: HState;
    [Psi, r, lam, URl] = homotopy(UL, UR, ctx);
    method = "homotopy";
    if (!accept(r, ctx) && lam > 0) {
      // why did the continuation stall? A limit reached on the way is reported as such
      let lim: string;
      try {
        lim = limitReason(Psi, assemble(Psi, UL, URl, ctx), UL, URl, ctx);
      } catch (e) {
        if (!(e instanceof DomainError)) throw e;
        lim = "none";
      }
      if (lim !== "none") return failed(prob, "RegularLimit", lim, { frame: F, ctx, psi: Psi, r, method: "homotopy_stalled" });
    }
  }
  if (!accept(r, ctx)) return failed(prob, "NoConvergence", "no_convergence", { frame: F, ctx, psi: Psi, r, method });
  // recording the waves (dense fan output) and the checks can still break down
  // numerically in degenerate cases; that must never escape as an exception
  let waves: Wave[], chk: CheckReport, lim: string;
  try {
    waves = assemble(Psi, UL, UR, ctx);
    chk = checkWaves(waves, ctx, unmirrored(F), prob.gamma);
    lim = limitReason(Psi, waves, UL, UR, ctx);
  } catch (e) {
    if (!(e instanceof DomainError)) throw e;
    return failed(prob, "NoConvergence", "assembly_failed", { frame: F, ctx, psi: Psi, r, method });
  }
  const sol = (rc2: RetCode, reason2: string): RiemannSolution =>
    ({ retcode: rc2, reason: reason2, prob, psi: Psi, residual: r, waves, frame: F, ctx, check: chk, method });
  if (lim !== "none") return sol("RegularLimit", lim);
  if (!chk.ok) return sol("CheckFailed", "check");
  return sol("Success", "none");
}

/** Re-run the independent checks on a solution (Julia `check(sol)`). */
export function check(sol: RiemannSolution): CheckReport {
  if (sol.ctx === null || sol.frame === null) throw new ArgumentError("solution has no frame");
  return checkWaves(sol.waves, sol.ctx, unmirrored(sol.frame), sol.prob.gamma);
}
