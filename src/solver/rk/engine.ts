// Adaptive explicit Runge-Kutta integrator with dense output (plan §4.1).
//
// The step and its error estimate follow OrdinaryDiffEqVerner's perform_step for a
// constant cache; step-size control follows OrdinaryDiffEqCore's defaults for an
// explicit method of order p: PI controller with β₁ = 7/(10p), β₂ = 2/(5p), γ = 9/10,
// q ∈ [1/qmax, 1/qmin] with qmin = 1/5, qmax = 10 (10000 on the first step),
// qoldinit = 1e-4; the error norm is the RMS of err_i / (abstol + reltol·max(|y_i|, |ŷ_i|));
// the initial step is Hairer's heuristic (initdt.jl). Exact agreement with Julia is not
// required (plan §4.1): at tol 1e-12 any correct controller agrees far below the
// acceptance limits.
//
// Two modes:
// - adaptive: as above, recording the accepted step sizes;
// - replay: integrate with a given step sequence and no error control. Replaying the
//   steps of a base run makes the solution a smooth function of the ODE's parameters,
//   which the finite-difference Jacobian needs (frozen-step "internal numerical
//   differentiation", plan §4.3).

export interface InterpWeight {
  /** stage index (0-based) */
  stage: number;
  /** b_stage(Θ) = Θ^pow · Σ_k coeffs[k] Θ^k */
  pow: number;
  coeffs: number[];
}

export interface Tableau {
  name: string;
  order: number;
  /** number of stages of one step; stages beyond are dense-output stages */
  stages: number;
  c: number[];
  /** a[i] = sparse row of the Butcher matrix: [j, a_ij] */
  a: [number, number][][];
  /** weights of the solution */
  b: [number, number][];
  /** weights of the error estimate: err = dt Σ btilde_j k_j */
  btilde: [number, number][];
  /** dense output: y(t + Θ dt) = y + dt Σ b_j(Θ) k_j */
  interp: InterpWeight[];
}

export type Rhs = (y: readonly number[], t: number) => number[];

export type OdeRetcode = "Success" | "MaxIters" | "DtLessThanMin" | "DtNaN" | "Unstable";

export interface OdeOptions {
  abstol: number;
  reltol: number;
  /** maximum number of steps (accepted + rejected), default 100_000 as OrdinaryDiffEq */
  maxiters?: number;
  /** keep the stages of every step for dense output */
  dense?: boolean;
  /** replay these step sizes instead of adaptive stepping */
  steps?: readonly number[];
}

interface StepData {
  t: number;
  dt: number;
  y: number[];
  /** stage derivatives; dense-output stages are filled on first use */
  k: (number[] | undefined)[];
}

/** Result of an integration. `at(t)` evaluates the dense output (requires `dense`). */
export class OdeSolution {
  readonly retcode: OdeRetcode;
  /** accepted time nodes, t[0] = t0 */
  readonly t: number[];
  /** solution at the nodes */
  readonly u: number[][];
  /** accepted step sizes (signed) */
  readonly steps: number[];
  /** number of right-hand side evaluations */
  nf: number;
  private readonly data: StepData[] | null;
  private readonly tab: Tableau;
  private readonly f: Rhs;

  constructor(retcode: OdeRetcode, t: number[], u: number[][], steps: number[], data: StepData[] | null,
    tab: Tableau, f: Rhs, nf: number) {
    this.retcode = retcode;
    this.t = t;
    this.u = u;
    this.steps = steps;
    this.data = data;
    this.tab = tab;
    this.f = f;
    this.nf = nf;
  }

  get successful(): boolean {
    return this.retcode === "Success";
  }

  /** solution at the last node */
  get end(): number[] {
    return this.u[this.u.length - 1]!;
  }

  /** Dense output at time τ (within the integrated interval). */
  at(tau: number): number[] {
    const data = this.data;
    if (data === null) throw new Error("OdeSolution: no dense output (integrate with dense = true)");
    const n = data.length;
    if (n === 0) return this.u[0]!.slice();
    const t0 = this.t[0]!, t1 = this.t[n]!;
    const dir = t1 >= t0 ? 1 : -1;
    if (tau === t1) return this.u[n]!.slice();
    if (tau === t0) return this.u[0]!.slice();
    // interval i with t[i] <= τ < t[i+1] (in the direction of integration)
    let lo = 0, hi = n - 1;
    while (lo < hi) {
      const mid = (lo + hi + 1) >> 1;
      if (dir * (tau - this.t[mid]!) >= 0) lo = mid; else hi = mid - 1;
    }
    const s = data[lo]!;
    const theta = (tau - s.t) / s.dt;
    this.ensureDenseStages(s);
    const y = s.y.slice();
    for (const w of this.tab.interp) {
      let p = 0;
      for (let k = w.coeffs.length - 1; k >= 0; k--) p = p * theta + w.coeffs[k]!;
      let bw = p;
      for (let k = 0; k < w.pow; k++) bw *= theta;
      const kk = s.k[w.stage]!;
      const f = s.dt * bw;
      for (let i = 0; i < y.length; i++) y[i] = y[i]! + f * kk[i]!;
    }
    return y;
  }

  private ensureDenseStages(s: StepData): void {
    const tab = this.tab;
    for (let i = tab.stages; i < tab.c.length; i++) {
      if (s.k[i] !== undefined) continue;
      s.k[i] = this.f(stageArg(tab, i, s.y, s.dt, s.k), s.t + tab.c[i]! * s.dt);
      this.nf++;
    }
  }
}

function stageArg(tab: Tableau, i: number, y: readonly number[], dt: number, k: readonly (number[] | undefined)[]): number[] {
  const g = y.slice();
  for (const [j, aij] of tab.a[i]!) {
    const kj = k[j]!;
    const f = dt * aij;
    for (let m = 0; m < g.length; m++) g[m] = g[m]! + f * kj[m]!;
  }
  return g;
}

function rms(x: readonly number[]): number {
  let s = 0;
  for (const v of x) s += v * v;
  return Math.sqrt(s / x.length);
}

/** One step: new state, error estimate (scaled RMS), stages. */
function rkStep(tab: Tableau, f: Rhs, t: number, y: readonly number[], dt: number,
  abstol: number, reltol: number): { y: number[]; err: number; k: (number[] | undefined)[] } {
  const k: (number[] | undefined)[] = new Array(tab.c.length);
  k[0] = f(y, t);
  for (let i = 1; i < tab.stages; i++) k[i] = f(stageArg(tab, i, y, dt, k), t + tab.c[i]! * dt);
  const ynew = y.slice();
  for (const [j, bj] of tab.b) {
    const kj = k[j]!, fct = dt * bj;
    for (let m = 0; m < ynew.length; m++) ynew[m] = ynew[m]! + fct * kj[m]!;
  }
  const e = new Array<number>(y.length).fill(0);
  for (const [j, bj] of tab.btilde) {
    const kj = k[j]!, fct = dt * bj;
    for (let m = 0; m < e.length; m++) e[m] = e[m]! + fct * kj[m]!;
  }
  for (let m = 0; m < e.length; m++) {
    e[m] = e[m]! / (abstol + Math.max(Math.abs(y[m]!), Math.abs(ynew[m]!)) * reltol);
  }
  return { y: ynew, err: rms(e), k };
}

/** Julia `eps(x)` for a finite double: distance to the next larger double in magnitude. */
function epsOf(x: number): number {
  const a = Math.abs(x);
  if (a < 2.2250738585072014e-308) return 5e-324;
  return Math.pow(2, Math.floor(Math.log2(a)) - 52);
}

/** OrdinaryDiffEqCore `_ode_initdt_oop` (Hairer & Wanner, Sec. II.4). */
function initialDt(tab: Tableau, f: Rhs, t0: number, y0: readonly number[], tdir: number, dtmax: number,
  abstol: number, reltol: number, count: () => void): number {
  const dtmin = Math.max(epsOf(t0), 0) + 5e-324;
  const smalldt = Math.max(dtmin, 1e-6);
  const sk = y0.map((v) => abstol + Math.abs(v) * reltol);
  const d0 = rms(y0.map((v, i) => v / sk[i]!));
  const f0 = f(y0, t0); count();
  if (f0.some(Number.isNaN)) return tdir * dtmin;
  const d1 = rms(f0.map((v, i) => v / sk[i]!));
  if (Number.isNaN(d1)) return tdir * dtmin;
  let dt0 = d0 < 1e-5 || d1 < 1e-5 ? smalldt : d0 / d1 / 100;
  dt0 = Math.min(dt0, dtmax);
  const y1 = y0.map((v, i) => v + tdir * dt0 * f0[i]!);
  const f1 = f(y1, t0 + tdir * dt0); count();
  if (f0.every((v, i) => v === f1[i])) return tdir * Math.max(dtmin, 100 * dt0);
  const d2 = rms(f1.map((v, i) => (v - f0[i]!) / sk[i]!)) / dt0;
  const m = Math.max(d1, d2);
  const dt1 = m <= 1e-15 ? Math.max(smalldt, dt0 * 1e-3) : Math.pow(10, -(2 + Math.log10(m)) / tab.order);
  return tdir * Math.max(dtmin, Math.min(100 * dt0, dt1, dtmax));
}

/**
 * Integrate y' = f(y, t) from t0 to t1. Exceptions thrown by f (DomainError, …)
 * propagate to the caller, as in Julia.
 */
export function integrate(tab: Tableau, f: Rhs, y0: readonly number[], t0: number, t1: number,
  opts: OdeOptions): OdeSolution {
  const { abstol, reltol } = opts;
  const maxiters = opts.maxiters ?? 100_000;
  const dense = opts.dense ?? false;
  const tdir = t1 >= t0 ? 1 : -1;
  const ts = [t0], us = [y0.slice()], steps: number[] = [];
  const data: StepData[] | null = dense ? [] : null;
  let nf = 0;
  const count = () => { nf++; };
  const counted: Rhs = (y, t) => { nf++; return f(y, t); };
  const done = (rc: OdeRetcode) => new OdeSolution(rc, ts, us, steps, data, tab, f, nf);
  if (t0 === t1) return done("Success");

  let t = t0;
  let y = y0.slice();

  if (opts.steps !== undefined) {
    // replay: fixed steps, no error control
    for (const dt of opts.steps) {
      const r = rkStep(tab, counted, t, y, dt, abstol, reltol);
      data?.push({ t, dt, y, k: r.k });
      t = t + dt; y = r.y;
      if (!y.every(Number.isFinite)) return done("Unstable");
      ts.push(t); us.push(y); steps.push(dt);
    }
    return done("Success");
  }

  const order = tab.order;
  const beta1 = 7 / (10 * order), beta2 = 2 / (5 * order), gamma = 0.9;
  const qmin = 0.2, qmax = 10, qmaxFirst = 10000, qoldinit = 1e-4;
  const dtmax = Math.abs(t1 - t0);
  let dt = initialDt(tab, f, t0, y, tdir, dtmax, abstol, reltol, count);
  let errold = qoldinit;
  let iters = 0, accepted = 0;
  for (;;) {
    if (iters >= maxiters) return done("MaxIters");
    iters++;
    // land exactly on t1 (modify_dt_for_tstops!)
    const dist = Math.abs(t1 - t);
    const tol = 100 * epsOf(Math.max(Math.abs(t), Math.abs(t1)));
    let last = false;
    if (Math.abs(dt) + tol >= dist) { dt = tdir * dist; last = true; }
    const r = rkStep(tab, counted, t, y, dt, abstol, reltol);
    const err = r.err;
    if (Number.isNaN(err)) return done("DtNaN");
    let q: number, q11: number;
    if (err === 0) {
      q = 1 / (accepted === 0 ? qmaxFirst : qmax);
      q11 = 1;
    } else {
      q11 = Math.pow(err, beta1);
      q = q11 / Math.pow(errold, beta2);
      q = Math.min(Math.max(q / gamma, 1 / (accepted === 0 ? qmaxFirst : qmax)), 1 / qmin);
    }
    if (err <= 1) {
      data?.push({ t, dt, y, k: r.k });
      t = last ? t1 : t + dt;
      y = r.y;
      if (y.some(Number.isNaN)) return done("Unstable");
      ts.push(t); us.push(y); steps.push(dt);
      accepted++;
      if (last) return done("Success");
      errold = Math.max(err, qoldinit);
      dt = dt / q;
    } else {
      dt = dt / Math.min(1 / qmin, q11 / gamma);
    }
    if (Math.abs(dt) <= epsOf(t)) return done("DtLessThanMin");
  }
}
