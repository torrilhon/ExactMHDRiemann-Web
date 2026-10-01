// Core types: options, return codes, internal half-states, waves, solutions.
// Mirrors types.jl.

import type { OdeSolution } from "./rk/engine.ts";
import type { StepTape } from "./rk/tape.ts";
import { logD } from "./math.ts";

/**
 * Outcome of `solve` (Julia `RetCode`). Only Success carries a checked solution.
 * - Success        converged and every independent check passed
 * - InvalidInput   non-finite values, ρ ≤ 0, p ≤ 0, γ ≤ 1 or Bx differing between L and R
 * - Unsupported    input outside the supported domain (see `reason`)
 * - RegularLimit   converged, but the solution sits at a regular-wave limit
 * - NoConvergence  the nonlinear solve failed from every start
 * - CheckFailed    solver converged but the independent checks failed (a bug)
 */
export type RetCode = "Success" | "InvalidInput" | "Unsupported" | "RegularLimit" | "NoConvergence" | "CheckFailed";

/** Julia `SolverOptions`: all thresholds of the solver (canonical units). */
export interface SolverOptions {
  /** Bn/√p at or below which the quasi-Euler solver (Bn = 0) is used */
  bn_euler: number;
  /** min Bt/√p of the input state with the larger value */
  bt_min_larger: number;
  /** min Bt/√p of the input state with the smaller value */
  bt_min_smaller: number;
  /** min Bt/√p of all middle states (else RegularLimit) */
  bt_floor: number;
  /** max ρ_R/ρ_L, p_R/p_L (and inverse) */
  ratio_max: number;
  /** safety margin below the slow-shock limit */
  delta: number;
  /** cap of slow-fan strength s = log(ρ0/ρ) */
  s_vac: number;
  /** tanh argument above which a map counts as saturated */
  sat: number;
  /** min ρ, p of middle states (canonical) */
  rho_floor: number;
  /** residual tolerance (∞-norm, canonical) */
  abstol: number;
  /** nonlinear iterations per start */
  maxiters: number;
  homotopy: boolean;
  homotopy_maxsteps: number;
  /** wall-clock budget of the homotopy fallback [s] */
  time_limit: number;
  ode_tol: number;
  /** relative tolerance of the post-solve checks */
  check_tol: number;
  verbose: boolean;
}

export function defaultOptions(overrides: Partial<SolverOptions> = {}): SolverOptions {
  return {
    bn_euler: 1e-10,
    bt_min_larger: 1e-4,
    bt_min_smaller: 1e-6,
    bt_floor: 1e-8,
    ratio_max: 1e6,
    delta: 1e-6,
    s_vac: logD(1e8),
    sat: 5.0,
    rho_floor: 1e-8,
    abstol: 1e-12,
    maxiters: 60,
    homotopy: true,
    homotopy_maxsteps: 400,
    time_limit: 5.0,
    ode_tol: 1e-12,
    check_tol: 1e-8,
    verbose: false,
    ...overrides,
  };
}

/**
 * Internal state along one side of the solution, canonical frame (Julia `HState`):
 * ρ, u (normal velocity), p, bt = |Bt| ≥ 0, φ = direction of Bt, vt = (vy, vz).
 */
export interface HState {
  readonly rho: number;
  readonly u: number;
  readonly p: number;
  readonly bt: number;
  readonly phi: number;
  readonly vt: readonly [number, number];
}

export function hstate(rho: number, u: number, p: number, bt: number, phi: number, vt: readonly [number, number]): HState {
  return { rho, u, p, bt, phi, vt };
}

/** Julia `Ctx`. `tape` (TS only) records/replays fan ODE steps for the Jacobian. */
export interface Ctx {
  readonly gamma: number;
  /** (γ+1)/(γ-1) */
  readonly kappa: number;
  /** canonical, > 0 */
  readonly Bn: number;
  readonly opts: SolverOptions;
  /** residual units: pressure, velocity, magnetic field */
  readonly scale: readonly [number, number, number];
  tape?: StepTape;
}

export function makeCtx(gamma: number, Bn: number, opts: SolverOptions,
  scale: readonly [number, number, number] = [1.0, 1.0, 1.0]): Ctx {
  return { gamma, kappa: (gamma + 1) / (gamma - 1), Bn, opts, scale };
}

export type FanFamily = "fast" | "slow" | "fast0";

/** Data of a rarefaction fan (dense ODE output in the fan parameter τ ∈ [0, 1]). */
export interface FanData {
  readonly family: FanFamily;
  readonly sigma: number;
  /** state at τ = 0 (outer end) */
  readonly up: HState;
  /** ψ (fast: ℓ = ψτ) or s_end (slow: s = s_end τ; fast0: s = s* τ) */
  readonly par: number;
  readonly sol: OdeSolution;
}

export type WaveKind = "fast_shock" | "fast_fan" | "rotation" | "slow_shock" | "slow_fan" | "contact" | "tangential";

/**
 * One elementary wave, canonical frame, in physical left-to-right order: `left` is the
 * state on its left, `right` on its right. For a fan s_left < s_right.
 */
export interface Wave {
  readonly kind: WaveKind;
  readonly side: number;
  readonly s_left: number;
  readonly s_right: number;
  readonly left: HState;
  readonly right: HState;
  readonly fan: FanData | null;
}

/** Symmetry map between the user frame and the canonical frame. */
export interface Frame {
  /** x → -x combined with B → -B (sides swapped, Bn keeps its sign) */
  readonly mirror: boolean;
  /** sign flip of B */
  readonly sB: number;
  /** rotation angle of the transverse plane */
  readonly theta: number;
  /** transverse velocity of R after flip/rotation (unscaled) */
  readonly vtR: readonly [number, number];
  readonly rho0: number;
  readonly p0: number;
}

/** Primitive state (ρ, vx, vy, vz, Bx, By, Bz, p). */
export type Prim = readonly number[];

/** Index constants of primitive 8-vectors (Julia W[1] … W[8]). */
export const RHO = 0, VX = 1, VY = 2, VZ = 3, BX = 4, BY = 5, BZ = 6, P = 7;

export interface RiemannProblem {
  readonly L: Prim;
  readonly R: Prim;
  readonly gamma: number;
}

export function riemannProblem(L: readonly number[], R: readonly number[], gamma = 5 / 3): RiemannProblem {
  return { L: L.map(Number), R: R.map(Number), gamma: Number(gamma) };
}

export interface CheckReport {
  readonly ok: boolean;
  readonly maxerr: number;
  readonly messages: string[];
}

/**
 * Result of `solve`: retcode, reason, prob, Ψ (path variables), residual, waves
 * (canonical frame), frame, ctx, check (report of the independent checks), method.
 */
export interface RiemannSolution {
  readonly retcode: RetCode;
  readonly reason: string;
  readonly prob: RiemannProblem;
  readonly psi: readonly number[];
  readonly residual: number;
  readonly waves: readonly Wave[];
  readonly frame: Frame | null;
  readonly ctx: Ctx | null;
  readonly check: CheckReport | null;
  readonly method: string;
}

export const successful = (sol: RiemannSolution) => sol.retcode === "Success";
