import { describe, expect, it } from "vitest";
import { integrate } from "../../src/solver/rk/engine.ts";
import { VERN9 } from "../../src/solver/rk/vern9.ts";
import { DomainError } from "../../src/solver/math.ts";

const T = VERN9;
const dot = (w: [number, number][], f: (j: number) => number) => w.reduce((s, [j, x]) => s + x * f(j), 0);

describe("Vern9 tableau (generated from OrdinaryDiffEqVerner)", () => {
  it("is consistent", () => {
    expect(T.c.length).toBe(26);
    T.a.forEach((row, i) => {
      const s = row.reduce((acc, [, x]) => acc + x, 0);
      const amax = Math.max(1, ...row.map(([, x]) => Math.abs(x)));
      expect(Math.abs(s - T.c[i]!)).toBeLessThanOrEqual(64 * Number.EPSILON * amax);
      row.forEach(([j]) => expect(j).toBeLessThan(i));        // explicit
    });
    expect(dot(T.b, () => 1)).toBeCloseTo(1, 14);
    expect(Math.abs(dot(T.btilde, () => 1))).toBeLessThan(1e-14);
  });
  it("satisfies the order conditions up to order 4", () => {
    const c = (j: number) => T.c[j]!;
    const Ac = (i: number) => T.a[i]!.reduce((s, [j, x]) => s + x * c(j), 0);
    const Ac2 = (i: number) => T.a[i]!.reduce((s, [j, x]) => s + x * c(j) ** 2, 0);
    const AAc = (i: number) => T.a[i]!.reduce((s, [j, x]) => s + x * Ac(j), 0);
    expect(dot(T.b, c)).toBeCloseTo(1 / 2, 13);
    expect(dot(T.b, (j) => c(j) ** 2)).toBeCloseTo(1 / 3, 13);
    expect(dot(T.b, Ac)).toBeCloseTo(1 / 6, 13);
    expect(dot(T.b, (j) => c(j) ** 3)).toBeCloseTo(1 / 4, 13);
    expect(dot(T.b, (j) => c(j) * Ac(j))).toBeCloseTo(1 / 8, 13);
    expect(dot(T.b, Ac2)).toBeCloseTo(1 / 12, 13);
    expect(dot(T.b, AAc)).toBeCloseTo(1 / 24, 13);
  });
  it("has an interpolant that reproduces the step at Θ = 1", () => {
    const at1 = new Map<number, number>();
    for (const w of T.interp) at1.set(w.stage, w.coeffs.reduce((s, x) => s + x, 0));
    // to the rounding of the (large) interpolation coefficients
    for (const [j, bj] of T.b) expect(Math.abs((at1.get(j) ?? 0) - bj)).toBeLessThan(1e-11);
    for (const [j, v] of at1) if (!T.b.some(([k]) => k === j)) expect(Math.abs(v)).toBeLessThan(1e-11);
  });
});

// y' = cos(t) y, y(0) = 1: y = exp(sin t)
const f = (y: readonly number[], t: number) => [Math.cos(t) * y[0]!, -y[1]!];
const exact = (t: number) => [Math.exp(Math.sin(t)), Math.exp(-t)];

describe("adaptive integration", () => {
  it("meets the tolerance and converges as tol is lowered", () => {
    let prev = Infinity;
    for (const tol of [1e-6, 1e-8, 1e-10, 1e-12]) {
      const s = integrate(T, f, [1, 1], 0, 3, { abstol: tol, reltol: tol });
      expect(s.retcode).toBe("Success");
      expect(s.t[s.t.length - 1]).toBe(3);
      const e = Math.max(...s.end.map((v, i) => Math.abs(v - exact(3)[i]!)));
      expect(e).toBeLessThan(50 * tol);
      expect(e).toBeLessThan(prev);
      prev = e;
    }
  });
  it("has 9th-order convergence in fixed steps", () => {
    const err = (n: number) => {
      const s = integrate(T, f, [1, 1], 0, 6, { abstol: 1, reltol: 1, steps: new Array(n).fill(6 / n) });
      return Math.abs(s.end[0]! - exact(6)[0]!);
    };
    expect(Math.log2(err(6) / err(12))).toBeGreaterThan(9);
    expect(Math.log2(err(12) / err(24))).toBeGreaterThan(9);
  });
  it("integrates backwards", () => {
    const s = integrate(T, f, exact(2), 2, 0.5, { abstol: 1e-12, reltol: 1e-12 });
    expect(s.end[0]).toBeCloseTo(exact(0.5)[0]!, 11);
    expect(s.t[s.t.length - 1]).toBe(0.5);
  });
});

describe("dense output", () => {
  const s = integrate(T, f, [1, 1], 0, 3, { abstol: 1e-12, reltol: 1e-12, dense: true });
  it("is accurate to ~tol inside the steps", () => {
    let worst = 0;
    for (let i = 0; i < 1000; i++) {
      const tau = (3 * (i + 0.5)) / 1000;
      const y = s.at(tau), e = exact(tau);
      worst = Math.max(worst, Math.abs(y[0]! - e[0]!), Math.abs(y[1]! - e[1]!));
    }
    expect(worst).toBeLessThan(1e-11);
  });
  it("returns the nodes exactly", () => {
    expect(s.at(3)).toEqual(s.end);
    expect(s.at(0)).toEqual([1, 1]);
    expect(s.at(s.t[2]!)).toEqual(s.u[2]);
  });
});

describe("replay", () => {
  it("reproduces an adaptive run bit for bit from its own steps", () => {
    const a = integrate(T, f, [1, 1], 0, 3, { abstol: 1e-12, reltol: 1e-12 });
    const b = integrate(T, f, [1, 1], 0, 3, { abstol: 1e-12, reltol: 1e-12, steps: a.steps });
    expect(b.end).toEqual(a.end);
  });
  it("is smooth in a parameter", () => {
    // y' = λ y: the replayed solution is differentiable in λ to round-off
    const g = (lam: number) => (y: readonly number[]) => [lam * y[0]!];
    const base = integrate(T, g(-1.3), [1], 0, 1, { abstol: 1e-12, reltol: 1e-12 });
    const h = 1e-6;
    const yp = integrate(T, g(-1.3 + h), [1], 0, 1, { abstol: 1, reltol: 1, steps: base.steps }).end[0]!;
    const ym = integrate(T, g(-1.3 - h), [1], 0, 1, { abstol: 1, reltol: 1, steps: base.steps }).end[0]!;
    expect((yp - ym) / (2 * h)).toBeCloseTo(Math.exp(-1.3), 9);       // d/dλ e^λ = e^λ
  });
});

describe("failures", () => {
  it("passes exceptions of the right-hand side through", () => {
    const bad = () => { throw new DomainError(-1); };
    expect(() => integrate(T, bad, [1], 0, 1, { abstol: 1e-12, reltol: 1e-12 })).toThrow(DomainError);
  });
  it("reports NaN and step-count failures", () => {
    const nan = (y: readonly number[], t: number) => [t > 0.5 ? NaN : y[0]!];
    expect(integrate(T, nan, [1], 0, 1, { abstol: 1e-12, reltol: 1e-12 }).successful).toBe(false);
    const stiff = (y: readonly number[]) => [-1e6 * y[0]!];
    expect(integrate(T, stiff, [1], 0, 1, { abstol: 1e-12, reltol: 1e-12, maxiters: 50 }).retcode).toBe("MaxIters");
  });
});
