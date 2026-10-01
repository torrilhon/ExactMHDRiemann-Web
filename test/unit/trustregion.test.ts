import { describe, expect, it } from "vitest";
import { trustRegion } from "../../src/solver/trustregion.ts";
import { fdJacobian } from "../../src/solver/jacobian.ts";
import { loadGolden } from "../golden.ts";

type F = (u: readonly number[]) => number[];
const sq = (x: number) => x * x;
const PROBLEMS: Record<string, F> = {
  rosenbrock: (u) => [10 * (u[1]! - sq(u[0]!)), 1 - u[0]!],
  powell_badly_scaled: (u) => [1e4 * u[0]! * u[1]! - 1, Math.exp(-u[0]!) + Math.exp(-u[1]!) - 1.0001],
  helical_valley: (u) => {
    const th = Math.atan2(u[1]!, u[0]!) / (2 * Math.PI) + (u[0]! < 0 ? 0.5 : 0);
    return [10 * (u[2]! - 10 * th), 10 * (Math.hypot(u[0]!, u[1]!) - 1), u[2]!];
  },
  powell_singular: (u) => [u[0]! + 10 * u[1]!, Math.sqrt(5) * (u[2]! - u[3]!), sq(u[1]! - 2 * u[2]!), Math.sqrt(10) * sq(u[0]! - u[3]!)],
  trigonometric5: (u) => {
    const s = u.reduce((a, x) => a + Math.cos(x), 0);
    return u.map((x, i) => 5 - s + (i + 1) * (1 - Math.cos(x)) - Math.sin(x));
  },
  broyden_tridiagonal5: (u) => u.map((x, i) => (3 - 2 * x) * x - (u[i - 1] ?? 0) - 2 * (u[i + 1] ?? 0) + 1),
  singular_start: (u) => [sq(u[0]!) + sq(u[1]!) - 1, u[0]! - u[1]!],
  freudenstein_roth: (u) => [-13 + u[0]! + ((5 - u[1]!) * u[1]! - 2) * u[1]!, -29 + u[0]! + ((u[1]! + 1) * u[1]! - 14) * u[1]!],
  linear5: (u) => [0, 1, 2, 3, 4].map((i) => u.reduce((a, x, j) => a + (i === j ? 4 : 1 / (i + j + 2)) * x, 0) - (i + 1)),
  big_plateau: (u) => (Math.abs(u[0]!) > 2 ? [1e6, 1e6] : [u[0]! ** 3 - 0.5, u[1]! - u[0]!]),
};

// helical valley starts on the branch cut of atan(y, x): a central difference across
// it is meaningless, ForwardDiff differentiates one side (the problem, not the solver)
const EXACT_J: Record<string, (u: readonly number[]) => number[][]> = {
  helical_valley: (u) => {
    const r2 = sq(u[0]!) + sq(u[1]!), r = Math.sqrt(r2);
    const t1 = -u[1]! / (2 * Math.PI * r2), t2 = u[0]! / (2 * Math.PI * r2);
    return [[-100 * t1, -100 * t2, 10], [10 * u[0]! / r, 10 * u[1]! / r, 0], [0, 0, 1]];
  },
};

interface TrCase { name: string; u0: number[]; u: number[]; retcode: string; resnorm: number; nsteps: number }
const G = loadGolden<{ cases: TrCase[] }>("trustregion");

describe("trust region vs NonlinearSolve.TrustRegion (MGH test problems)", () => {
  const rows = G.cases.map((c) => {
    const F = PROBLEMS[c.name]!;
    const jac = (u: readonly number[], fu: readonly number[]) => (EXACT_J[c.name] ?? ((x) => fdJacobian(F, x, fu).J))(u);
    const r = trustRegion(F, jac, c.u0, { abstol: 1e-12, maxiters: 60 });
    return { c, r };
  });
  it.each(rows.map(({ c, r }) => [`${c.name} ${JSON.stringify(c.u0)}`, c, r] as const))("%s", (_n, c, r) => {
    expect(r.retcode).toBe(c.retcode);
    if (c.retcode === "Success") {
      expect(r.resnorm).toBeLessThanOrEqual(1e-12);
      r.u.forEach((v, i) => expect(Math.abs(v - c.u[i]!)).toBeLessThan(1e-6 * Math.max(1, Math.abs(c.u[i]!))));
    }
  });
  it("needs a similar number of steps (report)", () => {
    const table = rows.map(({ c, r }) => `${c.name.padEnd(22)} julia ${String(c.nsteps).padStart(2)}  ts ${String(r.nsteps).padStart(2)}`);
    console.log(table.join("\n"));
    const ratio = rows.filter(({ c }) => c.retcode === "Success").map(({ c, r }) => r.nsteps / Math.max(1, c.nsteps));
    expect(Math.max(...ratio)).toBeLessThan(2.5);
  });
});
